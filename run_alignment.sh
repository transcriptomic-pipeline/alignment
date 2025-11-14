#!/bin/bash

# Alignment Module - Main Execution Script
# Supports STAR and HISAT2 with auto-installation, auto-indexing, and "both" mode
# Uses pre-installed python3 or creates Miniconda python3 environment for HISAT2

set -euo pipefail

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CONFIG_FILE="${SCRIPT_DIR}/config/install_paths.conf"
REF_CONFIG="${SCRIPT_DIR}/config/reference_paths.conf"

log_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Default parameters
INPUT_DIR=""
OUTPUT_DIR=""
SAMPLE_FILE=""
THREADS=12
ALIGNER=""
REFERENCE_GENOME=""
REFERENCE_GTF=""
GENOME_INDEX_DIR=""
PAIRED_END=true

# Reference selection flags
USE_DEFAULT_REF="no"
USE_PREVIOUS_REF="no"

# Tool paths
STAR_BIN=""
HISAT2_BIN=""
HISAT2_BUILD_BIN=""
SAMTOOLS_BIN=""
STAR_INDEX_DIR=""
HISAT2_INDEX_DIR=""
HISAT2_INDEX_PREFIX=""

# Conda settings (loaded from config)
CONDA_DIR=""
CONDA_ENV_NAME=""

# Timeout for auto-continue prompts (in seconds)
PROMPT_TIMEOUT=30

usage() {
    cat << EOF
Usage: $0 -i <input_dir> -o <output_dir> [options]

Required:
  -i, --input <dir>           Input directory with trimmed FASTQ files (e.g., qc_results/trimmed/PE)
  -o, --output <dir>          Output directory for BAM files

Optional:
  --aligner <star|hisat2|both> Which aligner to use (default: interactive prompt)
  -t, --threads <num>         Number of threads (default: 12)
  
Reference Genome Options:
  --reference <file>          Path to custom reference FASTA file
  --gtf <file>                Path to custom GTF annotation file
  --reference-default         Use default GRCh38 reference (skip prompt)
  --reference-previous        Use previously configured reference (skip prompt)
  --genome-index <dir>        Pre-built genome index directory
  
Other Options:
  --single-end                Process single-end reads (default: paired-end)
  -h, --help                  Show this help

Examples:
  # Interactive mode (prompts for aligner and reference)
  $0 -i qc_output/trimmed/PE -o alignment_results

  # Use default reference with STAR (non-interactive)
  $0 -i qc_output/trimmed/PE -o alignment_results \\
    --aligner star --reference-default -t 20

  # Use previous reference with both aligners
  $0 -i qc_output/trimmed/PE -o alignment_results \\
    --aligner both --reference-previous -t 16

  # Use custom reference genome
  $0 -i qc_output/trimmed/PE -o alignment_results \\
    --aligner star \\
    --reference /path/to/genome.fa \\
    --gtf /path/to/genes.gtf \\
    -t 20

  # Use pre-built HISAT2 index
  $0 -i qc_output/trimmed/PE -o alignment_results \\
    --aligner hisat2 \\
    --genome-index /path/to/hisat2_index \\
    -t 12

EOF
    exit 1
}

# Load configuration
check_configuration() {
    if [ ! -f "$CONFIG_FILE" ]; then
        log_error "Configuration not found: $CONFIG_FILE"
        log_info "Please run: bash install.sh"
        exit 1
    fi
    
    source "$CONFIG_FILE"
    
    STAR_BIN="${STAR_BIN:-}"
    HISAT2_BIN="${HISAT2_BIN:-}"
    HISAT2_BUILD_BIN="${HISAT2_BUILD_BIN:-}"
    SAMTOOLS_BIN="${SAMTOOLS_BIN:-}"
    
    log_success "Configuration loaded"
}

# Prompt user to select/confirm reference genome
prompt_reference_selection() {
    echo ""
    echo "========================================"
    echo "  Reference Genome Selection"
    echo "========================================"
    echo ""
    
    # Check if there's a remembered reference
    local HAS_PREVIOUS="no"
    local PREV_TYPE=""
    local PREV_FASTA=""
    local PREV_GTF=""
    
    if [ -f "$REF_CONFIG" ]; then
        source "$REF_CONFIG"
        if [ -n "$REFERENCE_FASTA" ] && [ -f "$REFERENCE_FASTA" ]; then
            HAS_PREVIOUS="yes"
            PREV_TYPE="$USE_CUSTOM_REFERENCE"
            PREV_FASTA="$REFERENCE_FASTA"
            PREV_GTF="$REFERENCE_GTF"
        fi
    fi
    
    # Show options
    log_info "Select reference genome to use:"
    echo ""
    echo "  1) Default: Human GRCh38 (Ensembl 113)"
    
    if [ "$HAS_PREVIOUS" = "yes" ]; then
        if [ "$PREV_TYPE" = "yes" ]; then
            echo "  2) Previous custom: $(basename $PREV_FASTA)"
        else
            echo "  2) Previous default: GRCh38 (Ensembl 113)"
        fi
        echo "  3) New custom reference"
    else
        echo "  2) Custom reference"
    fi
    
    echo ""
    
    # Get user choice with timeout
    local timeout=30
    echo -ne "Enter choice [1-"
    [ "$HAS_PREVIOUS" = "yes" ] && echo -ne "3" || echo -ne "2"
    echo -ne "] (default: "
    [ "$HAS_PREVIOUS" = "yes" ] && echo -ne "2" || echo -ne "1"
    echo -ne ", auto-select in ${timeout}s): "
    
    read -t $timeout -n 1 -r REPLY || true
    echo
    
    # Default behavior
    if [ -z "$REPLY" ]; then
        if [ "$HAS_PREVIOUS" = "yes" ]; then
            REPLY="2"
            log_info "No response in ${timeout}s, using previous reference"
        else
            REPLY="1"
            log_info "No response in ${timeout}s, using default (GRCh38)"
        fi
    fi
    
    # Handle selection
    case "${REPLY}" in
        1)
            # Use default GRCh38
            USE_CUSTOM_REFERENCE="no"
            REFERENCE_GENOME="${SCRIPT_DIR}/references/GRCh38_ensembl113/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
            REFERENCE_GTF="${SCRIPT_DIR}/references/GRCh38_ensembl113/Homo_sapiens.GRCh38.113.gtf"
            STAR_INDEX_DIR="${SCRIPT_DIR}/references/GRCh38_ensembl113/STAR_index"
            HISAT2_INDEX_DIR="${SCRIPT_DIR}/references/GRCh38_ensembl113/HISAT2_index"
            HISAT2_INDEX_PREFIX="${HISAT2_INDEX_DIR}/genome"
            
            # Update config
            update_reference_config
            
            log_success "Selected: Default GRCh38 (Ensembl 113)"
            log_info "FASTA: $(basename $REFERENCE_GENOME)"
            log_info "GTF: $(basename $REFERENCE_GTF)"
            ;;
            
        2)
            if [ "$HAS_PREVIOUS" = "yes" ]; then
                # Use previous reference
                USE_CUSTOM_REFERENCE="$PREV_TYPE"
                REFERENCE_GENOME="$PREV_FASTA"
                REFERENCE_GTF="$PREV_GTF"
                
                # Load index paths from config
                source "$REF_CONFIG"
                
                log_success "Using previous reference"
                log_info "FASTA: $REFERENCE_GENOME"
                log_info "GTF: $REFERENCE_GTF"
            else
                # Prompt for custom reference
                prompt_custom_reference_interactive
            fi
            ;;
            
        3)
            if [ "$HAS_PREVIOUS" = "yes" ]; then
                # Prompt for new custom reference
                prompt_custom_reference_interactive
            else
                log_error "Invalid choice"
                exit 1
            fi
            ;;
            
        *)
            log_error "Invalid choice. Please select a valid option."
            exit 1
            ;;
    esac
    
    echo ""
}

# Prompt for custom reference paths
prompt_custom_reference_interactive() {
    echo ""
    log_info "Custom reference genome configuration"
    echo ""
    
    # Prompt for FASTA
    read -p "Enter path to genome FASTA file: " CUSTOM_FASTA
    CUSTOM_FASTA="${CUSTOM_FASTA/#\~/$HOME}"
    
    # Validate FASTA
    if [ ! -f "$CUSTOM_FASTA" ]; then
        log_error "FASTA file not found: $CUSTOM_FASTA"
        echo ""
        log_info "Please provide a valid path and re-run."
        exit 1
    fi
    
    # Prompt for GTF
    read -p "Enter path to GTF annotation file: " CUSTOM_GTF
    CUSTOM_GTF="${CUSTOM_GTF/#\~/$HOME}"
    
    # Validate GTF
    if [ ! -f "$CUSTOM_GTF" ]; then
        log_error "GTF file not found: $CUSTOM_GTF"
        echo ""
        log_info "Please provide a valid path and re-run."
        exit 1
    fi
    
    # Set variables
    USE_CUSTOM_REFERENCE="yes"
    REFERENCE_GENOME="$CUSTOM_FASTA"
    REFERENCE_GTF="$CUSTOM_GTF"
    
    # Create index paths based on genome basename
    local GENOME_BASENAME=$(basename "${REFERENCE_GENOME}" | sed 's/\.[^.]*$//')
    local REF_DIR="$(dirname "$REFERENCE_GENOME")"
    STAR_INDEX_DIR="${REF_DIR}/STAR_index_${GENOME_BASENAME}"
    HISAT2_INDEX_DIR="${REF_DIR}/HISAT2_index_${GENOME_BASENAME}"
    HISAT2_INDEX_PREFIX="${HISAT2_INDEX_DIR}/genome"
    
    # Update config
    update_reference_config
    
    log_success "Custom reference configured:"
    log_info "FASTA: $REFERENCE_GENOME"
    log_info "GTF: $REFERENCE_GTF"
    log_info "Index base: ${GENOME_BASENAME}"
}

# Update reference configuration file
update_reference_config() {
    mkdir -p "${SCRIPT_DIR}/config"
    
    cat > "$REF_CONFIG" << EOF
# Reference Genome Configuration
# Generated: $(date)
# Last modified by: run_alignment.sh

USE_CUSTOM_REFERENCE="${USE_CUSTOM_REFERENCE}"
REFERENCE_DIR="$(dirname "$REFERENCE_GENOME")"
REFERENCE_FASTA="${REFERENCE_GENOME}"
REFERENCE_GTF="${REFERENCE_GTF}"
STAR_INDEX_DIR="${STAR_INDEX_DIR}"
HISAT2_INDEX_DIR="${HISAT2_INDEX_DIR}"
HISAT2_INDEX_PREFIX="${HISAT2_INDEX_PREFIX}"
EOF
}

# Check reference genome configuration
check_reference() {
    # Priority 1: Custom reference provided via CLI
    if [ -n "$REFERENCE_GENOME" ] && [ -n "$REFERENCE_GTF" ]; then
        if [ ! -f "$REFERENCE_GENOME" ]; then
            log_error "Reference genome not found: $REFERENCE_GENOME"
            exit 1
        fi
        if [ ! -f "$REFERENCE_GTF" ]; then
            log_error "Reference GTF not found: $REFERENCE_GTF"
            exit 1
        fi
        
        # Setup index paths for custom reference
        local GENOME_BASENAME=$(basename "${REFERENCE_GENOME}" | sed 's/\.[^.]*$//')
        local REF_DIR="$(dirname "$REFERENCE_GENOME")"
        STAR_INDEX_DIR="${REF_DIR}/STAR_index_${GENOME_BASENAME}"
        HISAT2_INDEX_DIR="${REF_DIR}/HISAT2_index_${GENOME_BASENAME}"
        HISAT2_INDEX_PREFIX="${HISAT2_INDEX_DIR}/genome"
        
        # Update config
        USE_CUSTOM_REFERENCE="yes"
        update_reference_config
        
        log_success "Using CLI-specified custom reference:"
        log_info "FASTA: $REFERENCE_GENOME"
        log_info "GTF: $REFERENCE_GTF"
        return 0
    fi
    
    # Priority 2: Use default reference (--reference-default flag)
    if [ "$USE_DEFAULT_REF" = "yes" ]; then
        USE_CUSTOM_REFERENCE="no"
        REFERENCE_GENOME="${SCRIPT_DIR}/references/GRCh38_ensembl113/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        REFERENCE_GTF="${SCRIPT_DIR}/references/GRCh38_ensembl113/Homo_sapiens.GRCh38.113.gtf"
        STAR_INDEX_DIR="${SCRIPT_DIR}/references/GRCh38_ensembl113/STAR_index"
        HISAT2_INDEX_DIR="${SCRIPT_DIR}/references/GRCh38_ensembl113/HISAT2_index"
        HISAT2_INDEX_PREFIX="${HISAT2_INDEX_DIR}/genome"
        
        # Update config
        update_reference_config
        
        log_success "Using default reference: GRCh38 (Ensembl 113)"
        log_info "FASTA: $(basename $REFERENCE_GENOME)"
        log_info "GTF: $(basename $REFERENCE_GTF)"
        
        # Validate files exist
        if [ ! -f "$REFERENCE_GENOME" ]; then
            log_error "Default reference not found: $REFERENCE_GENOME"
            log_info "Please run: bash install.sh --download-reference"
            exit 1
        fi
        if [ ! -f "$REFERENCE_GTF" ]; then
            log_error "Default GTF not found: $REFERENCE_GTF"
            log_info "Please run: bash install.sh --download-reference"
            exit 1
        fi
        
        return 0
    fi
    
    # Priority 3: Use previous reference (--reference-previous flag)
    if [ "$USE_PREVIOUS_REF" = "yes" ]; then
        if [ ! -f "$REF_CONFIG" ]; then
            log_error "No previous reference configuration found"
            log_info "Please run without --reference-previous to configure reference"
            exit 1
        fi
        
        source "$REF_CONFIG"
        
        if [ ! -f "$REFERENCE_FASTA" ]; then
            log_error "Previous reference not found: $REFERENCE_FASTA"
            log_info "Please reconfigure reference"
            exit 1
        fi
        
        REFERENCE_GENOME="$REFERENCE_FASTA"
        REFERENCE_GTF="$REFERENCE_GTF"
        
        if [ "$USE_CUSTOM_REFERENCE" = "yes" ]; then
            log_success "Using previous custom reference:"
        else
            log_success "Using previous default reference:"
        fi
        log_info "FASTA: $REFERENCE_GENOME"
        log_info "GTF: $REFERENCE_GTF"
        
        return 0
    fi
    
    # Priority 4: Interactive mode - prompt for reference selection
    prompt_reference_selection
    
    # Validate selected reference
    if [ ! -f "$REFERENCE_GENOME" ]; then
        log_error "Reference genome not found: $REFERENCE_GENOME"
        log_info "Please run: bash install.sh --download-reference"
        exit 1
    fi
    
    if [ ! -f "$REFERENCE_GTF" ]; then
        log_error "Reference GTF not found: $REFERENCE_GTF"
        log_info "Please run: bash install.sh --download-reference"
        exit 1
    fi
    
    log_success "Reference genome: $REFERENCE_GENOME"
    log_success "Reference GTF: $REFERENCE_GTF"
}

# Check STAR installation
check_star() {
    if [ -z "$STAR_BIN" ] || [ ! -f "$STAR_BIN" ]; then
        log_error "STAR not found"
        log_info "Please run: bash install.sh --aligner star"
        exit 1
    fi
    log_success "STAR is available"
}

# Check HISAT2 installation
check_hisat2() {
    if [ -z "$HISAT2_BIN" ] || [ ! -f "$HISAT2_BIN" ]; then
        log_error "HISAT2 not found"
        log_info "Please run: bash install.sh --aligner hisat2"
        exit 1
    fi
    log_success "HISAT2 is available"
}

# Build STAR genome index
build_star_index() {
    log_warning "STAR index not found. Building STAR genome index..."
    log_info "This may take 30-60 minutes and requires ~32 GB RAM"
    
    echo -ne "Build STAR index now? [Y/n] (auto-continue in ${PROMPT_TIMEOUT}s): "
    read -t $PROMPT_TIMEOUT -n 1 -r || true
    echo
    
    if [[ ! $REPLY =~ ^[Nn]$ ]]; then
        log_info "Running STAR genome index build..."
        
        mkdir -p "$STAR_INDEX_DIR"
        
        "$STAR_BIN" \
            --runThreadN "$THREADS" \
            --runMode genomeGenerate \
            --genomeDir "$STAR_INDEX_DIR" \
            --genomeFastaFiles "$REFERENCE_GENOME" \
            --sjdbGTFfile "$REFERENCE_GTF" \
            --sjdbOverhang 100
        
        if [ $? -eq 0 ]; then
            log_success "STAR index built successfully: $STAR_INDEX_DIR"
        else
            log_error "STAR index build failed"
            exit 1
        fi
    else
        log_error "STAR index required for alignment"
        exit 1
    fi
}

# Check or build STAR index
check_or_build_star_index() {
    log_info "Checking STAR genome index..."
    
    # Show which reference is being used
    if [ -f "$REF_CONFIG" ]; then
        source "$REF_CONFIG"
        if [ "$USE_CUSTOM_REFERENCE" = "yes" ]; then
            log_info "Custom reference: $(basename $REFERENCE_FASTA)"
        else
            log_info "Default reference: GRCh38 (Ensembl 113)"
        fi
        log_info "Index directory: $STAR_INDEX_DIR"
    fi
    
    # Check if index directory exists AND contains the required genomeParameters.txt file
    if [ ! -d "$STAR_INDEX_DIR" ] || [ -z "$(ls -A $STAR_INDEX_DIR 2>/dev/null)" ] || [ ! -f "${STAR_INDEX_DIR}/genomeParameters.txt" ]; then
        log_warning "STAR index not found or incomplete"
        echo ""
        log_info "Index location: $STAR_INDEX_DIR"
        log_info "Reference genome: $REFERENCE_GENOME"
        log_info "Reference GTF: $REFERENCE_GTF"
        echo ""
        build_star_index
    else
        log_success "STAR index found: $STAR_INDEX_DIR"
    fi
    
    GENOME_INDEX_DIR="$STAR_INDEX_DIR"
}

# Build HISAT2 genome index
build_hisat2_index() {
    log_warning "HISAT2 index not found. Building HISAT2 genome index..."
    log_info "This may take 15-30 minutes and requires ~8 GB RAM"
    
    echo -ne "Build HISAT2 index now? [Y/n] (auto-continue in ${PROMPT_TIMEOUT}s): "
    read -t $PROMPT_TIMEOUT -n 1 -r || true
    echo
    
    if [[ ! $REPLY =~ ^[Nn]$ ]]; then
        log_info "Running HISAT2 genome index build..."
        
        mkdir -p "$HISAT2_INDEX_DIR"
        
        "$HISAT2_BUILD_BIN" \
            -p "$THREADS" \
            "$REFERENCE_GENOME" \
            "$HISAT2_INDEX_PREFIX"
        
        if [ $? -eq 0 ]; then
            log_success "HISAT2 index built successfully: $HISAT2_INDEX_DIR"
        else
            log_error "HISAT2 index build failed"
            exit 1
        fi
    else
        log_error "HISAT2 index required for alignment"
        exit 1
    fi
}

# Check or build HISAT2 index
check_or_build_hisat2_index() {
    log_info "Checking HISAT2 genome index..."
    
    # Show which reference is being used
    if [ -f "$REF_CONFIG" ]; then
        source "$REF_CONFIG"
        if [ "$USE_CUSTOM_REFERENCE" = "yes" ]; then
            log_info "Custom reference: $(basename $REFERENCE_FASTA)"
        else
            log_info "Default reference: GRCh38 (Ensembl 113)"
        fi
        log_info "Index directory: $HISAT2_INDEX_DIR"
    fi
    
    # Check if index exists
    if [ ! -d "$HISAT2_INDEX_DIR" ] || [ -z "$(ls -A $HISAT2_INDEX_DIR 2>/dev/null)" ]; then
        log_warning "HISAT2 index not found"
        echo ""
        log_info "Index location: $HISAT2_INDEX_DIR"
        log_info "Reference genome: $REFERENCE_GENOME"
        log_info "Reference GTF: $REFERENCE_GTF"
        echo ""
        build_hisat2_index
    else
        log_success "HISAT2 index found: $HISAT2_INDEX_DIR"
    fi
    
    GENOME_INDEX_DIR="$HISAT2_INDEX_PREFIX"
}

# Detect file pairing pattern
detect_pairing_pattern() {
    local INPUT_DIRECTORY="$1"
    
    log_info "Detecting file pairing pattern..."
    
    # Count files with different patterns
    local R1_UNDERSCORE=$(find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "*_R1_*.fastq" -o -name "*_R1_*.fq" -o -name "*_R1_*.fastq.gz" -o -name "*_R1_*.fq.gz" \) 2>/dev/null | wc -l)
    local R2_UNDERSCORE=$(find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "*_R2_*.fastq" -o -name "*_R2_*.fq" -o -name "*_R2_*.fastq.gz" -o -name "*_R2_*.fq.gz" \) 2>/dev/null | wc -l)
    
    local R1_DOT=$(find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "*_R1.fastq" -o -name "*_R1.fq" -o -name "*_R1.fastq.gz" -o -name "*_R1.fq.gz" \) 2>/dev/null | wc -l)
    local R2_DOT=$(find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "*_R2.fastq" -o -name "*_R2.fq" -o -name "*_R2.fastq.gz" -o -name "*_R2.fq.gz" \) 2>/dev/null | wc -l)
    
    local NUM_UNDERSCORE=$(find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "*_1.fastq" -o -name "*_1.fq" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) 2>/dev/null | wc -l)
    local NUM2_UNDERSCORE=$(find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "*_2.fastq" -o -name "*_2.fq" -o -name "*_2.fastq.gz" -o -name "*_2.fq.gz" \) 2>/dev/null | wc -l)
    
    # Determine pattern based on counts
    if [ "$R1_UNDERSCORE" -gt 0 ] && [ "$R2_UNDERSCORE" -gt 0 ] && [ "$R1_UNDERSCORE" -eq "$R2_UNDERSCORE" ]; then
        PAIRING_PATTERN="_R1_/_R2_"
        log_success "Detected pattern: _R1_/_R2_"
        echo "$R1_UNDERSCORE"
    elif [ "$R1_DOT" -gt 0 ] && [ "$R2_DOT" -gt 0 ] && [ "$R1_DOT" -eq "$R2_DOT" ]; then
        PAIRING_PATTERN="_R1/_R2"
        log_success "Detected pattern: _R1/_R2"
        echo "$R1_DOT"
    elif [ "$NUM_UNDERSCORE" -gt 0 ] && [ "$NUM2_UNDERSCORE" -gt 0 ] && [ "$NUM_UNDERSCORE" -eq "$NUM2_UNDERSCORE" ]; then
        PAIRING_PATTERN="_1/_2"
        log_success "Detected pattern: _1/_2"
        echo "$NUM_UNDERSCORE"
    else
        log_error "Could not detect consistent pairing pattern"
        log_info "Supported patterns: _R1_/_R2_, _R1/_R2, _1/_2"
        log_info "Found: R1_=$R1_UNDERSCORE, R2_=$R2_UNDERSCORE, R1=$R1_DOT, R2=$R2_DOT, 1=$NUM_UNDERSCORE, 2=$NUM2_UNDERSCORE"
        exit 1
    fi
}

# Get sample list based on pattern
get_sample_list() {
    local INPUT_DIRECTORY="$1"
    local PATTERN="$2"
    
    case "$PATTERN" in
        "_R1_/_R2_")
            find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "*_R1_*.fastq" -o -name "*_R1_*.fq" -o -name "*_R1_*.fastq.gz" -o -name "*_R1_*.fq.gz" \) | while read R1_FILE; do
                local BASENAME=$(basename "$R1_FILE")
                echo "${BASENAME%%_R1_*}"
            done | sort -u
            ;;
        "_R1/_R2")
            find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "*_R1.fastq" -o -name "*_R1.fq" -o -name "*_R1.fastq.gz" -o -name "*_R1.fq.gz" \) | while read R1_FILE; do
                local BASENAME=$(basename "$R1_FILE")
                echo "${BASENAME%_R1.*}"
            done | sort -u
            ;;
        "_1/_2")
            find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "*_1.fastq" -o -name "*_1.fq" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) | while read R1_FILE; do
                local BASENAME=$(basename "$R1_FILE")
                echo "${BASENAME%_1.*}"
            done | sort -u
            ;;
    esac
}

# Find paired files
find_paired_files() {
    local SAMPLE_ID="$1"
    local INPUT_DIRECTORY="$2"
    local PATTERN="$3"
    
    local R1_FILE=""
    local R2_FILE=""
    
    case "$PATTERN" in
        "_R1_/_R2_")
            R1_FILE=$(find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "${SAMPLE_ID}_R1_*.fastq" -o -name "${SAMPLE_ID}_R1_*.fq" -o -name "${SAMPLE_ID}_R1_*.fastq.gz" -o -name "${SAMPLE_ID}_R1_*.fq.gz" \) | head -n 1)
            R2_FILE=$(find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "${SAMPLE_ID}_R2_*.fastq" -o -name "${SAMPLE_ID}_R2_*.fq" -o -name "${SAMPLE_ID}_R2_*.fastq.gz" -o -name "${SAMPLE_ID}_R2_*.fq.gz" \) | head -n 1)
            ;;
        "_R1/_R2")
            R1_FILE=$(find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "${SAMPLE_ID}_R1.fastq" -o -name "${SAMPLE_ID}_R1.fq" -o -name "${SAMPLE_ID}_R1.fastq.gz" -o -name "${SAMPLE_ID}_R1.fq.gz" \) | head -n 1)
            R2_FILE=$(find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "${SAMPLE_ID}_R2.fastq" -o -name "${SAMPLE_ID}_R2.fq" -o -name "${SAMPLE_ID}_R2.fastq.gz" -o -name "${SAMPLE_ID}_R2.fq.gz" \) | head -n 1)
            ;;
        "_1/_2")
            R1_FILE=$(find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "${SAMPLE_ID}_1.fastq" -o -name "${SAMPLE_ID}_1.fq" -o -name "${SAMPLE_ID}_1.fastq.gz" -o -name "${SAMPLE_ID}_1.fq.gz" \) | head -n 1)
            R2_FILE=$(find "$INPUT_DIRECTORY" -maxdepth 1 -type f \( -name "${SAMPLE_ID}_2.fastq" -o -name "${SAMPLE_ID}_2.fq" -o -name "${SAMPLE_ID}_2.fastq.gz" -o -name "${SAMPLE_ID}_2.fq.gz" \) | head -n 1)
            ;;
    esac
    
    echo "$R1_FILE|$R2_FILE"
}

# Confirm detected samples
confirm_samples() {
    local SAMPLE_COUNT="$1"
    
    echo ""
    echo "========================================"
    echo "  Pairing Pattern Detection"
    echo "========================================"
    log_success "Detected pattern: $PAIRING_PATTERN"
    echo ""
    log_info "Found $SAMPLE_COUNT R1 files and $(find "$INPUT_DIR" -maxdepth 1 -type f -name "*R2*" 2>/dev/null | wc -l) R2 files"
    echo ""
    log_info "Sample preview (first 3):"
    
    local SAMPLES=($(get_sample_list "$INPUT_DIR" "$PAIRING_PATTERN"))
    local PREVIEW_COUNT=0
    for SAMPLE_ID in "${SAMPLES[@]}"; do
        [ $PREVIEW_COUNT -ge 3 ] && break
        PREVIEW_COUNT=$((PREVIEW_COUNT + 1))
        
        local FILES=$(find_paired_files "$SAMPLE_ID" "$INPUT_DIR" "$PAIRING_PATTERN")
        local R1=$(echo "$FILES" | cut -d'|' -f1)
        local R2=$(echo "$FILES" | cut -d'|' -f2)
        
        echo "  Sample $PREVIEW_COUNT: $SAMPLE_ID"
        echo "    R1: $(basename "$R1")"
        if [ -n "$R2" ]; then
            echo "    R2: $(basename "$R2") ✓"
        else
            echo "    R2: NOT FOUND ✗"
        fi
    done
    
    echo ""
    log_warning "Please verify the detected pattern is correct"
    echo ""
    
    echo -ne "Continue with detected pattern? [Y/n] (auto-continue in ${PROMPT_TIMEOUT}s): "
    read -t $PROMPT_TIMEOUT -n 1 -r || true
    echo
    
    if [[ $REPLY =~ ^[Nn]$ ]]; then
        log_error "Pattern confirmation cancelled by user"
        log_info "Please rename files to match expected patterns and re-run"
        exit 1
    fi
    
    log_info "Proceeding with detected pattern: $PAIRING_PATTERN"
}

# Run STAR alignment
run_star() {
    local SAMPLE_ID="$1"
    local R1_FILE="$2"
    local R2_FILE="$3"
    local OUT_DIR="$4"
    local LOG_FILE="$5"
    
    log_info "Running STAR alignment for: $SAMPLE_ID"
    
    mkdir -p "$OUT_DIR"
    
    "$STAR_BIN" \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX_DIR" \
        --readFilesIn "$R1_FILE" "$R2_FILE" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${OUT_DIR}/${SAMPLE_ID}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --limitBAMsortRAM 31000000000 \
        2>&1 | tee "$LOG_FILE"
    
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        # Index BAM file
        "$SAMTOOLS_BIN" index "${OUT_DIR}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam"
        log_success "STAR alignment completed: $SAMPLE_ID"
        return 0
    else
        log_error "STAR alignment failed: $SAMPLE_ID"
        return 1
    fi
}

# Run HISAT2 alignment
run_hisat2() {
    local SAMPLE_ID="$1"
    local R1_FILE="$2"
    local R2_FILE="$3"
    local OUT_DIR="$4"
    local LOG_FILE="$5"
    
    log_info "Running HISAT2 alignment for: $SAMPLE_ID"
    
    mkdir -p "$OUT_DIR"
    
    "$HISAT2_BIN" \
        -p "$THREADS" \
        -x "$HISAT2_INDEX_PREFIX" \
        -1 "$R1_FILE" \
        -2 "$R2_FILE" \
        --rna-strandness RF \
        2>"$LOG_FILE" \
        | "$SAMTOOLS_BIN" view -bS - \
        | "$SAMTOOLS_BIN" sort -@ "$THREADS" -o "${OUT_DIR}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam"
    
    if [ ${PIPESTATUS[0]} -eq 0 ] && [ ${PIPESTATUS[2]} -eq 0 ]; then
        # Index BAM file
        "$SAMTOOLS_BIN" index "${OUT_DIR}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam"
        log_success "HISAT2 alignment completed: $SAMPLE_ID"
        return 0
    else
        log_error "HISAT2 alignment failed: $SAMPLE_ID"
        return 1
    fi
}

# Process all samples
process_all_samples() {
    local SAMPLES=($(get_sample_list "$INPUT_DIR" "$PAIRING_PATTERN"))
    local TOTAL=${#SAMPLES[@]}
    local COMPLETED=0
    local FAILED=0
    
    log_info "Processing $TOTAL samples with $ALIGNER..."
    echo ""
    
    for SAMPLE_ID in "${SAMPLES[@]}"; do
        local FILES=$(find_paired_files "$SAMPLE_ID" "$INPUT_DIR" "$PAIRING_PATTERN")
        local R1=$(echo "$FILES" | cut -d'|' -f1)
        local R2=$(echo "$FILES" | cut -d'|' -f2)
        
        if [ -z "$R1" ] || [ -z "$R2" ]; then
            log_warning "Skipping $SAMPLE_ID: Missing paired files"
            FAILED=$((FAILED + 1))
            continue
        fi
        
        case "$ALIGNER" in
            star)
                local OUT_DIR="$OUTPUT_DIR/$SAMPLE_ID"
                local LOG_FILE="$OUTPUT_DIR/logs/${SAMPLE_ID}_star.log"
                mkdir -p "$OUTPUT_DIR/logs"
                
                if [ -f "${OUT_DIR}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam" ]; then
                    log_info "Skipping $SAMPLE_ID: Already aligned"
                else
                    if run_star "$SAMPLE_ID" "$R1" "$R2" "$OUT_DIR" "$LOG_FILE"; then
                        COMPLETED=$((COMPLETED + 1))
                    else
                        FAILED=$((FAILED + 1))
                    fi
                fi
                ;;
                
            hisat2)
                local OUT_DIR="$OUTPUT_DIR/$SAMPLE_ID"
                local LOG_FILE="$OUTPUT_DIR/logs/${SAMPLE_ID}_hisat2.log"
                mkdir -p "$OUTPUT_DIR/logs"
                
                if [ -f "${OUT_DIR}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam" ]; then
                    log_info "Skipping $SAMPLE_ID: Already aligned"
                else
                    if run_hisat2 "$SAMPLE_ID" "$R1" "$R2" "$OUT_DIR" "$LOG_FILE"; then
                        COMPLETED=$((COMPLETED + 1))
                    else
                        FAILED=$((FAILED + 1))
                    fi
                fi
                ;;
                
            both)
                # Run both STAR and HISAT2
                local STAR_OUT="$OUTPUT_DIR/STAR/$SAMPLE_ID"
                local HISAT2_OUT="$OUTPUT_DIR/HISAT2/$SAMPLE_ID"
                local STAR_LOG="$OUTPUT_DIR/logs/${SAMPLE_ID}_star.log"
                local HISAT2_LOG="$OUTPUT_DIR/logs/${SAMPLE_ID}_hisat2.log"
                mkdir -p "$OUTPUT_DIR/logs"
                
                local STAR_SUCCESS=0
                local HISAT2_SUCCESS=0
                
                if [ -f "${STAR_OUT}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam" ]; then
                    log_info "Skipping STAR for $SAMPLE_ID: Already aligned"
                    STAR_SUCCESS=1
                else
                    if run_star "$SAMPLE_ID" "$R1" "$R2" "$STAR_OUT" "$STAR_LOG"; then
                        STAR_SUCCESS=1
                    fi
                fi
                
                if [ -f "${HISAT2_OUT}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam" ]; then
                    log_info "Skipping HISAT2 for $SAMPLE_ID: Already aligned"
                    HISAT2_SUCCESS=1
                else
                    if run_hisat2 "$SAMPLE_ID" "$R1" "$R2" "$HISAT2_OUT" "$HISAT2_LOG"; then
                        HISAT2_SUCCESS=1
                    fi
                fi
                
                if [ $STAR_SUCCESS -eq 1 ] && [ $HISAT2_SUCCESS -eq 1 ]; then
                    COMPLETED=$((COMPLETED + 1))
                else
                    FAILED=$((FAILED + 1))
                fi
                ;;
        esac
    done
    
    echo ""
    echo "========================================"
    echo "  Alignment Summary"
    echo "========================================"
    log_info "Total samples: $TOTAL"
    log_success "Completed: $COMPLETED"
    [ $FAILED -gt 0 ] && log_error "Failed: $FAILED" || log_info "Failed: 0"
    echo ""
}

# Prompt for aligner selection
prompt_aligner_selection() {
    echo ""
    echo "========================================"
    echo "  Aligner Selection"
    echo "========================================"
    echo ""
    log_info "Both STAR and HISAT2 are installed"
    echo ""
    echo "Which aligner do you want to use?"
    echo "  1) STAR (high memory ~32GB, best for splice detection)"
    echo "  2) HISAT2 (low memory ~8GB, faster)"
    echo "  3) BOTH (run alignment with both aligners for comparison)"
    echo ""
    read -p "Enter choice [1-3]: " ALIGNER_CHOICE
    
    case "${ALIGNER_CHOICE}" in
        1) ALIGNER="star" ;;
        2) ALIGNER="hisat2" ;;
        3) ALIGNER="both" ;;
        *) log_error "Invalid choice"; exit 1 ;;
    esac
    
    log_success "Selected aligner: $ALIGNER"
}

# Main execution
main() {
    echo "========================================"
    echo "  Alignment Module"
    echo "========================================"
    echo ""
    
    log_info "Checking for installed tools..."
    check_configuration
    
    echo ""
    echo "========================================"
    echo "  Aligner Selection"
    echo "========================================"
    echo ""
    
    # Check which aligners are installed
    local STAR_AVAILABLE="no"
    local HISAT2_AVAILABLE="no"
    
    [ -n "$STAR_BIN" ] && [ -f "$STAR_BIN" ] && STAR_AVAILABLE="yes"
    [ -n "$HISAT2_BIN" ] && [ -f "$HISAT2_BIN" ] && HISAT2_AVAILABLE="yes"
    
    if [ "$STAR_AVAILABLE" = "no" ] && [ "$HISAT2_AVAILABLE" = "no" ]; then
        log_error "No aligners installed"
        log_info "Please run: bash install.sh"
        exit 1
    fi
    
    # If aligner not specified, prompt or use only available one
    if [ -z "$ALIGNER" ]; then
        if [ "$STAR_AVAILABLE" = "yes" ] && [ "$HISAT2_AVAILABLE" = "yes" ]; then
            prompt_aligner_selection
        elif [ "$STAR_AVAILABLE" = "yes" ]; then
            ALIGNER="star"
            log_info "Only STAR is installed, using STAR"
        else
            ALIGNER="hisat2"
            log_info "Only HISAT2 is installed, using HISAT2"
        fi
    fi
    
    # Check selected aligner(s)
    case "$ALIGNER" in
        star)
            check_star
            ;;
        hisat2)
            check_hisat2
            # Check Python for HISAT2
            if [ -n "$USE_SYSTEM_PYTHON" ] && [ "$USE_SYSTEM_PYTHON" = "yes" ]; then
                log_info "Using system Python for HISAT2 (no conda needed)"
            fi
            ;;
        both)
            check_star
            check_hisat2
            if [ -n "$USE_SYSTEM_PYTHON" ] && [ "$USE_SYSTEM_PYTHON" = "yes" ]; then
                log_info "Using system Python for HISAT2 (no conda needed)"
            fi
            ;;
        *)
            log_error "Invalid aligner: $ALIGNER"
            exit 1
            ;;
    esac
    
    # Check reference genome
    log_info "Checking reference genome..."
    check_reference
    
    # Check/build indexes
    case "$ALIGNER" in
        star)
            check_or_build_star_index
            ;;
        hisat2)
            check_or_build_hisat2_index
            ;;
        both)
            check_or_build_star_index
            check_or_build_hisat2_index
            ;;
    esac
    
    # Detect pairing pattern
    local SAMPLE_COUNT=$(detect_pairing_pattern "$INPUT_DIR")
    
    # Confirm samples
    confirm_samples "$SAMPLE_COUNT"
    
    # Process all samples
    process_all_samples
    
    log_success "Alignment pipeline completed!"
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --aligner)
            case "$2" in
                star|hisat2|both)
                    ALIGNER="$2"
                    ;;
                *)
                    log_error "Invalid aligner: $2. Choose star, hisat2, or both"
                    exit 1
                    ;;
            esac
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --reference)
            REFERENCE_GENOME="$2"
            REFERENCE_GENOME="${REFERENCE_GENOME/#\~/$HOME}"
            shift 2
            ;;
        --gtf)
            REFERENCE_GTF="$2"
            REFERENCE_GTF="${REFERENCE_GTF/#\~/$HOME}"
            shift 2
            ;;
        --reference-default)
            USE_DEFAULT_REF="yes"
            shift
            ;;
        --reference-previous)
            USE_PREVIOUS_REF="yes"
            shift
            ;;
        --genome-index)
            GENOME_INDEX_DIR="$2"
            shift 2
            ;;
        --single-end)
            PAIRED_END=false
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    log_error "Missing required arguments"
    usage
fi

# Run main
main
