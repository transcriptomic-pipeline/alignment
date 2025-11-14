#!/bin/bash

# Alignment Module - Main Execution Script
# Supports STAR and HISAT2 with auto-installation, auto-indexing, and "both" mode

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

# Tool paths
STAR_BIN=""
HISAT2_BIN=""
HISAT2_BUILD_BIN=""
SAMTOOLS_BIN=""
STAR_INDEX_DIR=""
HISAT2_INDEX_DIR=""
HISAT2_INDEX_PREFIX=""
PYTHON_BIN=""
USE_SYSTEM_PYTHON="no"
CONDA_DIR=""
CONDA_ENV_NAME=""

# Default reference paths (loaded from config)
DEFAULT_FASTA=""
DEFAULT_GTF=""
DEFAULT_REFERENCE_DIR=""

# Usage
show_usage() {
    cat << EOF
Usage: $(basename "$0") -i INPUT_DIR -o OUTPUT_DIR [OPTIONS]

Required Arguments:
  -i, --input DIR          Input directory with trimmed FASTQ files
  -o, --output DIR         Output directory for aligned BAM files

Optional Arguments:
  -t, --threads NUM        Number of threads (default: 12)
  -a, --aligner TYPE       Aligner: star, hisat2, or both (default: prompt)
  -s, --samples FILE       Sample list file (auto-detected if not provided)
  
  --reference FASTA        Custom reference genome FASTA
  --gtf GTF                Custom reference GTF annotation
  
  -h, --help               Show this help message

Examples:
  # Interactive mode (prompts for aligner and reference)
  $(basename "$0") -i qc_output/trimmed/PE -o alignment_output -t 12
  
  # With specific aligner
  $(basename "$0") -i qc_output/trimmed/PE -o alignment_output -a star
  
  # With custom reference
  $(basename "$0") -i qc_output/trimmed/PE -o alignment_output \\
    --reference /path/to/genome.fa --gtf /path/to/genes.gtf
    
  # Run both aligners for comparison
  $(basename "$0") -i qc_output/trimmed/PE -o alignment_output -a both

EOF
}

# Parse command line arguments
parse_args() {
    if [ $# -eq 0 ]; then
        show_usage
        exit 1
    fi
    
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
            -t|--threads)
                THREADS="$2"
                shift 2
                ;;
            -a|--aligner)
                ALIGNER="$2"
                shift 2
                ;;
            -s|--samples)
                SAMPLE_FILE="$2"
                shift 2
                ;;
            --reference)
                REFERENCE_GENOME="$2"
                shift 2
                ;;
            --gtf)
                REFERENCE_GTF="$2"
                shift 2
                ;;
            -h|--help)
                show_usage
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                show_usage
                exit 1
                ;;
        esac
    done
    
    # Validate required arguments
    if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
        log_error "Input and output directories are required"
        show_usage
        exit 1
    fi
}

# Check if tools are installed
check_installation() {
    log_info "Checking for installed tools..."
    
    if [ ! -f "$CONFIG_FILE" ]; then
        log_warning "Alignment tools not installed"
        log_info "Installation required"
        echo ""
        
        # Determine which aligner to install
        local INSTALL_ALIGNER=""
        if [ -n "$ALIGNER" ]; then
            INSTALL_ALIGNER="$ALIGNER"
        else
            echo "Which aligner would you like to install?"
            echo "  1) STAR (recommended for high-memory systems)"
            echo "  2) HISAT2 (recommended for low-memory systems)"
            echo "  3) Both (choose at runtime)"
            echo ""
            read -p "Enter choice [1-3]: " INSTALL_CHOICE
            
            case "${INSTALL_CHOICE}" in
                1) INSTALL_ALIGNER="star" ;;
                2) INSTALL_ALIGNER="hisat2" ;;
                3) INSTALL_ALIGNER="both" ;;
                *) INSTALL_ALIGNER="both" ;;
            esac
        fi
        
        read -p "Run installation now? [Y/n] " -n 1 -r
        echo
        
        if [[ ! $REPLY =~ ^[Nn]$ ]]; then
            log_info "Running installer..."
            bash "${SCRIPT_DIR}/install.sh" --aligner "$INSTALL_ALIGNER"
            
            if [ $? -eq 0 ]; then
                log_success "Installation completed"
            else
                log_error "Installation failed"
                exit 1
            fi
        else
            log_error "Installation is required to proceed"
            exit 1
        fi
    fi
    
    # Load configuration
    if [ -f "$CONFIG_FILE" ]; then
        source "$CONFIG_FILE"
        log_success "Configuration loaded"
    else
        log_error "Configuration file not found: $CONFIG_FILE"
        exit 1
    fi
}

# Prompt for aligner selection
prompt_aligner_selection() {
    echo ""
    echo "========================================"
    echo "  Aligner Selection"
    echo "========================================"
    echo ""
    
    local STAR_AVAILABLE="no"
    local HISAT2_AVAILABLE="no"
    
    [ "$STAR_INSTALLED" = "yes" ] && [ -x "$STAR_BIN" ] && STAR_AVAILABLE="yes"
    [ "$HISAT2_INSTALLED" = "yes" ] && [ -x "$HISAT2_BIN" ] && HISAT2_AVAILABLE="yes"
    
    if [ "$STAR_AVAILABLE" = "yes" ] && [ "$HISAT2_AVAILABLE" = "yes" ]; then
        log_info "Both STAR and HISAT2 are installed"
        echo ""
        echo "Which aligner do you want to use?"
        echo "  1) STAR (high memory ~32GB, best for splice detection)"
        echo "  2) HISAT2 (low memory ~8GB, faster)"
        echo "  3) BOTH (run alignment with both aligners for comparison)"
        echo ""
        read -p "Enter choice [1-3]: " -n 1 -r CHOICE
        echo
        
        case "${CHOICE}" in
            1) ALIGNER="star" ;;
            2) ALIGNER="hisat2" ;;
            3) ALIGNER="both" ;;
            *) ALIGNER="both" ;;
        esac
    elif [ "$STAR_AVAILABLE" = "yes" ]; then
        ALIGNER="star"
        log_info "Using STAR (only aligner installed)"
    elif [ "$HISAT2_AVAILABLE" = "yes" ]; then
        ALIGNER="hisat2"
        log_info "Using HISAT2 (only aligner installed)"
    else
        log_error "No aligners are installed or available"
        exit 1
    fi
    
    log_success "Selected aligner: ${ALIGNER^^}"
}

# Check aligner availability
check_aligner_availability() {
    if [[ "$ALIGNER" == "star" || "$ALIGNER" == "both" ]]; then
        if [ ! -x "$STAR_BIN" ]; then
            log_error "STAR not found: $STAR_BIN"
            exit 1
        fi
        log_success "STAR is available"
    fi
    
    if [[ "$ALIGNER" == "hisat2" || "$ALIGNER" == "both" ]]; then
        if [ ! -x "$HISAT2_BIN" ]; then
            log_error "HISAT2 not found: $HISAT2_BIN"
            exit 1
        fi
        log_success "HISAT2 is available"
        
        # Check Python for HISAT2
        if [ "$USE_SYSTEM_PYTHON" = "yes" ]; then
            log_info "Using system Python for HISAT2 (no conda needed)"
        elif [ -d "$CONDA_DIR" ]; then
            log_info "Using Miniconda Python for HISAT2"
        else
            log_warning "Python environment not properly configured for HISAT2"
        fi
    fi
}

# Prompt for custom reference (interactive)
prompt_custom_reference_interactive() {
    echo ""
    log_info "Custom reference genome configuration"
    echo ""
    log_warning "You need to provide:"
    echo "  1. Genome FASTA file (.fa or .fasta)"
    echo "  2. Gene annotation GTF file (.gtf)"
    echo ""
    
    read -p "Enter path to genome FASTA file: " CUSTOM_FASTA
    read -p "Enter path to GTF annotation file: " CUSTOM_GTF
    
    # Expand tilde
    CUSTOM_FASTA="${CUSTOM_FASTA/#\~/$HOME}"
    CUSTOM_GTF="${CUSTOM_GTF/#\~/$HOME}"
    
    # Validate files
    if [ ! -f "$CUSTOM_FASTA" ]; then
        log_error "Reference FASTA not found: $CUSTOM_FASTA"
        exit 1
    fi
    
    if [ ! -f "$CUSTOM_GTF" ]; then
        log_error "Reference GTF not found: $CUSTOM_GTF"
        exit 1
    fi
    
    REFERENCE_GENOME="$CUSTOM_FASTA"
    REFERENCE_GTF="$CUSTOM_GTF"
    
    # Setup index paths for custom reference
    local GENOME_BASENAME=$(basename "${REFERENCE_GENOME}" | sed 's/\.[^.]*$//')
    local REF_DIR="$(dirname "$REFERENCE_GENOME")"
    STAR_INDEX_DIR="${REF_DIR}/STAR_index_${GENOME_BASENAME}"
    HISAT2_INDEX_DIR="${REF_DIR}/HISAT2_index_${GENOME_BASENAME}"
    HISAT2_INDEX_PREFIX="${HISAT2_INDEX_DIR}/genome"
    
    log_success "Custom reference validated:"
    log_info "FASTA: $REFERENCE_GENOME"
    log_info "GTF: $REFERENCE_GTF"
}

# SMART reference selection - ALWAYS prompts, shows options based on what exists
prompt_reference_selection() {
    echo ""
    echo "========================================"
    echo "  Reference Genome Selection"
    echo "========================================"
    echo ""
    
    # Check what exists
    local HAS_DEFAULT="no"
    local HAS_PREVIOUS="no"
    local PREV_FASTA=""
    local PREV_GTF=""
    
    # Check for default reference
    if [ -n "$DEFAULT_FASTA" ] && [ -f "$DEFAULT_FASTA" ] && [ -n "$DEFAULT_GTF" ] && [ -f "$DEFAULT_GTF" ]; then
        HAS_DEFAULT="yes"
    fi
    
    # Check for previous custom reference
    if [ -f "$REF_CONFIG" ]; then
        source "$REF_CONFIG"
        if [ -n "$REFERENCE_FASTA" ] && [ -f "$REFERENCE_FASTA" ] && [ -n "$REFERENCE_GTF" ] && [ -f "$REFERENCE_GTF" ]; then
            # Only consider it "previous" if it's NOT the default
            if [ "$REFERENCE_FASTA" != "$DEFAULT_FASTA" ]; then
                HAS_PREVIOUS="yes"
                PREV_FASTA="$REFERENCE_FASTA"
                PREV_GTF="$REFERENCE_GTF"
            fi
        fi
    fi
    
    # Build options dynamically based on what exists
    log_info "Select reference genome to use:"
    echo ""
    
    local option_num=1
    declare -A option_map
    
    # Option: Default (if exists)
    if [ "$HAS_DEFAULT" = "yes" ]; then
        echo "  ${option_num}) Default: GRCh38 (Ensembl 113) - found at ${DEFAULT_REFERENCE_DIR}"
        option_map[$option_num]="default"
        option_num=$((option_num + 1))
    fi
    
    # Option: Previous custom (if exists)
    if [ "$HAS_PREVIOUS" = "yes" ]; then
        echo "  ${option_num}) Previous custom: $(basename $PREV_FASTA) and $(basename $PREV_GTF)"
        option_map[$option_num]="previous"
        option_num=$((option_num + 1))
    fi
    
    # Option: New custom (always available)
    echo "  ${option_num}) New custom: Enter new genome paths"
    option_map[$option_num]="custom"
    
    echo ""
    local max_option=$((option_num))
    read -p "Enter choice [1-${max_option}]: " -n 1 -r CHOICE
    echo
    
    # Validate choice
    if [ -z "$CHOICE" ] || [ "$CHOICE" -lt 1 ] || [ "$CHOICE" -gt "$max_option" ]; then
        log_error "Invalid choice"
        exit 1
    fi
    
    # Process choice
    local selected_option="${option_map[$CHOICE]}"
    
    case "$selected_option" in
        default)
            REFERENCE_GENOME="$DEFAULT_FASTA"
            REFERENCE_GTF="$DEFAULT_GTF"
            STAR_INDEX_DIR="${DEFAULT_REFERENCE_DIR}/STAR_index"
            HISAT2_INDEX_DIR="${DEFAULT_REFERENCE_DIR}/HISAT2_index"
            HISAT2_INDEX_PREFIX="${HISAT2_INDEX_DIR}/genome"
            log_success "Using default GRCh38 reference"
            ;;
            
        previous)
            REFERENCE_GENOME="$PREV_FASTA"
            REFERENCE_GTF="$PREV_GTF"
            
            # Load index paths from previous config
            if [ -f "$REF_CONFIG" ]; then
                source "$REF_CONFIG"
                STAR_INDEX_DIR="$STAR_INDEX_DIR"
                HISAT2_INDEX_DIR="$HISAT2_INDEX_DIR"
                HISAT2_INDEX_PREFIX="$HISAT2_INDEX_PREFIX"
            else
                # Fallback: construct index paths
                local GENOME_BASENAME=$(basename "${REFERENCE_GENOME}" | sed 's/\.[^.]*$//')
                local REF_DIR="$(dirname "$REFERENCE_GENOME")"
                STAR_INDEX_DIR="${REF_DIR}/STAR_index_${GENOME_BASENAME}"
                HISAT2_INDEX_DIR="${REF_DIR}/HISAT2_index_${GENOME_BASENAME}"
                HISAT2_INDEX_PREFIX="${HISAT2_INDEX_DIR}/genome"
            fi
            
            log_success "Using previous custom reference:"
            log_info "FASTA: $REFERENCE_GENOME"
            log_info "GTF: $REFERENCE_GTF"
            ;;
            
        custom)
            prompt_custom_reference_interactive
            ;;
            
        *)
            log_error "Invalid selection"
            exit 1
            ;;
    esac
    
    # Save selection to config
    update_reference_config
}

# Update reference config
update_reference_config() {
    mkdir -p "${SCRIPT_DIR}/config"
    
    cat > "$REF_CONFIG" << EOF
# Reference Genome Configuration
# Last modified: $(date)

REFERENCE_FASTA="${REFERENCE_GENOME}"
REFERENCE_GTF="${REFERENCE_GTF}"
STAR_INDEX_DIR="${STAR_INDEX_DIR}"
HISAT2_INDEX_DIR="${HISAT2_INDEX_DIR}"
HISAT2_INDEX_PREFIX="${HISAT2_INDEX_PREFIX}"
EOF
}

# Check reference genome - ALWAYS prompts unless CLI args provided
check_reference() {
    log_info "Checking reference genome..."
    
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
        
        # Setup index paths
        local GENOME_BASENAME=$(basename "${REFERENCE_GENOME}" | sed 's/\.[^.]*$//')
        local REF_DIR="$(dirname "$REFERENCE_GENOME")"
        STAR_INDEX_DIR="${REF_DIR}/STAR_index_${GENOME_BASENAME}"
        HISAT2_INDEX_DIR="${REF_DIR}/HISAT2_index_${GENOME_BASENAME}"
        HISAT2_INDEX_PREFIX="${HISAT2_INDEX_DIR}/genome"
        
        log_success "Using CLI-specified reference:"
        log_info "FASTA: $REFERENCE_GENOME"
        log_info "GTF: $REFERENCE_GTF"
        
        update_reference_config
        return 0
    fi
    
    # Priority 2: ALWAYS prompt for reference selection (interactive)
    prompt_reference_selection
    
    # Validate selected reference
    if [ ! -f "$REFERENCE_GENOME" ]; then
        log_error "Reference genome not found: $REFERENCE_GENOME"
        exit 1
    fi
    
    if [ ! -f "$REFERENCE_GTF" ]; then
        log_error "Reference GTF not found: $REFERENCE_GTF"
        exit 1
    fi
}

# Check and build STAR genome index
check_star_index() {
    log_info "Checking STAR genome index..."
    
    local INDEX_COMPLETE="${STAR_INDEX_DIR}/SA"
    
    if [ -f "$INDEX_COMPLETE" ]; then
        log_success "STAR index found"
        return 0
    fi
    
    log_warning "STAR index not found or incomplete"
    log_info "Building STAR genome index..."
    log_warning "This may take 30-60 minutes and requires ~32 GB RAM"
    
    read -p "Build STAR index now? [Y/n] (auto-continue in 30s): " -t 30 -n 1 -r
    echo
    
    if [[ ! $REPLY =~ ^[Nn]$ ]]; then
        mkdir -p "$STAR_INDEX_DIR"
        
        log_info "Running STAR genome index build..."
        "$STAR_BIN" --runThreadN "$THREADS" \
            --runMode genomeGenerate \
            --genomeDir "$STAR_INDEX_DIR" \
            --genomeFastaFiles "$REFERENCE_GENOME" \
            --sjdbGTFfile "$REFERENCE_GTF" \
            --sjdbOverhang 100 2>&1 | tee "${SCRIPT_DIR}/star_index_build.log"
        
        if [ -f "$INDEX_COMPLETE" ]; then
            log_success "STAR index built successfully"
        else
            log_error "STAR index build failed"
            exit 1
        fi
    else
        log_error "STAR index is required for alignment"
        exit 1
    fi
}

# Check and build HISAT2 genome index
check_hisat2_index() {
    log_info "Checking HISAT2 genome index..."
    
    local INDEX_BASE="${HISAT2_INDEX_PREFIX}"
    
    if [ -f "${INDEX_BASE}.1.ht2" ]; then
        log_success "HISAT2 index found"
        return 0
    fi
    
    log_warning "HISAT2 index not found"
    log_info "Building HISAT2 genome index..."
    log_warning "This may take 20-40 minutes"
    
    read -p "Build HISAT2 index now? [Y/n] (auto-continue in 30s): " -t 30 -n 1 -r
    echo
    
    if [[ ! $REPLY =~ ^[Nn]$ ]]; then
        mkdir -p "$HISAT2_INDEX_DIR"
        
        log_info "Running HISAT2 genome index build..."
        "$HISAT2_BUILD_BIN" -p "$THREADS" "$REFERENCE_GENOME" "$INDEX_BASE" 2>&1 | tee "${SCRIPT_DIR}/hisat2_index_build.log"
        
        if [ -f "${INDEX_BASE}.1.ht2" ]; then
            log_success "HISAT2 index built successfully"
        else
            log_error "HISAT2 index build failed"
            exit 1
        fi
    else
        log_error "HISAT2 index is required for alignment"
        exit 1
    fi
}

# Detect read layout
detect_read_layout() {
    log_info "Detecting read layout..."
    
    local FIRST_FILE=$(find "$INPUT_DIR" -name "*.fastq.gz" -o -name "*.fq.gz" | head -n 1)
    
    if [ -z "$FIRST_FILE" ]; then
        log_error "No FASTQ files found in $INPUT_DIR"
        exit 1
    fi
    
    # Check for paired-end pattern
    if [[ "$FIRST_FILE" =~ _R1[_.]|_1[_.]|\.1[_.]|_F[_.] ]]; then
        PAIRED_END=true
        log_success "Detected: Paired-end reads"
    else
        PAIRED_END=false
        log_success "Detected: Single-end reads"
    fi
}

# Generate sample list
generate_sample_list() {
    if [ -n "$SAMPLE_FILE" ] && [ -f "$SAMPLE_FILE" ]; then
        log_info "Using provided sample list: $SAMPLE_FILE"
        return 0
    fi
    
    log_info "Generating sample list..."
    
    SAMPLE_FILE="${OUTPUT_DIR}/sample_list.txt"
    mkdir -p "$(dirname "$SAMPLE_FILE")"
    
    if [ "$PAIRED_END" = true ]; then
        find "$INPUT_DIR" -name "*_R1*.fastq.gz" -o -name "*_1.fastq.gz" | \
            sed 's/_R1.*//;s/_1\..*//' | sort -u > "$SAMPLE_FILE"
    else
        find "$INPUT_DIR" -name "*.fastq.gz" | \
            sed 's/\.fastq\.gz$//' | sort -u > "$SAMPLE_FILE"
    fi
    
    local SAMPLE_COUNT=$(wc -l < "$SAMPLE_FILE")
    log_success "Found $SAMPLE_COUNT samples"
}

# Run STAR alignment
run_star_alignment() {
    local SAMPLE=$1
    local SAMPLE_NAME=$(basename "$SAMPLE")
    
    log_info "Running STAR alignment for: $SAMPLE_NAME"
    
    local STAR_OUTPUT="${OUTPUT_DIR}/STAR/${SAMPLE_NAME}"
    mkdir -p "$STAR_OUTPUT"
    
    local READ1="${SAMPLE}_R1.fastq.gz"
    local READ2="${SAMPLE}_R2.fastq.gz"
    
    if [ "$PAIRED_END" = true ]; then
        "$STAR_BIN" --runThreadN "$THREADS" \
            --genomeDir "$STAR_INDEX_DIR" \
            --readFilesIn "$READ1" "$READ2" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${STAR_OUTPUT}/" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM NM MD \
            --quantMode GeneCounts
    else
        "$STAR_BIN" --runThreadN "$THREADS" \
            --genomeDir "$STAR_INDEX_DIR" \
            --readFilesIn "$READ1" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${STAR_OUTPUT}/" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM NM MD \
            --quantMode GeneCounts
    fi
    
    # Index BAM
    "$SAMTOOLS_BIN" index "${STAR_OUTPUT}/Aligned.sortedByCoord.out.bam"
    
    log_success "STAR alignment complete: $SAMPLE_NAME"
}

# Run HISAT2 alignment
run_hisat2_alignment() {
    local SAMPLE=$1
    local SAMPLE_NAME=$(basename "$SAMPLE")
    
    log_info "Running HISAT2 alignment for: $SAMPLE_NAME"
    
    local HISAT2_OUTPUT="${OUTPUT_DIR}/HISAT2/${SAMPLE_NAME}"
    mkdir -p "$HISAT2_OUTPUT"
    
    local READ1="${SAMPLE}_R1.fastq.gz"
    local READ2="${SAMPLE}_R2.fastq.gz"
    local BAM_FILE="${HISAT2_OUTPUT}/aligned.bam"
    
    if [ "$PAIRED_END" = true ]; then
        "$HISAT2_BIN" -p "$THREADS" \
            -x "$HISAT2_INDEX_PREFIX" \
            -1 "$READ1" \
            -2 "$READ2" \
            --summary-file "${HISAT2_OUTPUT}/alignment_summary.txt" | \
            "$SAMTOOLS_BIN" view -bS - | \
            "$SAMTOOLS_BIN" sort -@ "$THREADS" -o "$BAM_FILE"
    else
        "$HISAT2_BIN" -p "$THREADS" \
            -x "$HISAT2_INDEX_PREFIX" \
            -U "$READ1" \
            --summary-file "${HISAT2_OUTPUT}/alignment_summary.txt" | \
            "$SAMTOOLS_BIN" view -bS - | \
            "$SAMTOOLS_BIN" sort -@ "$THREADS" -o "$BAM_FILE"
    fi
    
    # Index BAM
    "$SAMTOOLS_BIN" index "$BAM_FILE"
    
    log_success "HISAT2 alignment complete: $SAMPLE_NAME"
}

# Main execution
main() {
    echo "========================================"
    echo "  Alignment Module"
    echo "========================================"
    echo ""
    
    # Check installation
    check_installation
    
    # Prompt for aligner if not specified
    if [ -z "$ALIGNER" ]; then
        prompt_aligner_selection
    fi
    
    # Check aligner availability
    check_aligner_availability
    
    # Check reference genome
    check_reference
    
    # Check/build genome indices
    if [[ "$ALIGNER" == "star" || "$ALIGNER" == "both" ]]; then
        check_star_index
    fi
    
    if [[ "$ALIGNER" == "hisat2" || "$ALIGNER" == "both" ]]; then
        check_hisat2_index
    fi
    
    # Detect read layout
    detect_read_layout
    
    # Generate sample list
    generate_sample_list
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    
    # Run alignment for each sample
    log_info "Starting alignment..."
    
    while IFS= read -r SAMPLE; do
        if [[ "$ALIGNER" == "star" || "$ALIGNER" == "both" ]]; then
            run_star_alignment "$SAMPLE"
        fi
        
        if [[ "$ALIGNER" == "hisat2" || "$ALIGNER" == "both" ]]; then
            run_hisat2_alignment "$SAMPLE"
        fi
    done < "$SAMPLE_FILE"
    
    echo ""
    log_success "Alignment complete!"
    log_info "Output directory: $OUTPUT_DIR"
    echo ""
}

# Parse arguments and run
parse_args "$@"
main
