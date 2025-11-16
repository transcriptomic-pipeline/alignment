#!/bin/bash

# Alignment Module - Main Execution Script
# Supports STAR and HISAT2 with comprehensive paired-end and single-end alignment

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

# Default reference paths
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
  # Interactive mode
  $(basename "$0") -i qc_output/trimmed/PE -o alignment_output -t 12
  
  # With specific aligner
  $(basename "$0") -i qc_output/trimmed/PE -o alignment_output -a star
  
  # With custom reference
  $(basename "$0") -i qc_output/trimmed/PE -o alignment_output \\
    --reference /path/to/genome.fa --gtf /path/to/genes.gtf
    
  # Run both aligners
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
            -i|--input) INPUT_DIR="$2"; shift 2 ;;
            -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
            -t|--threads) THREADS="$2"; shift 2 ;;
            -a|--aligner) ALIGNER="$2"; shift 2 ;;
            -s|--samples) SAMPLE_FILE="$2"; shift 2 ;;
            --reference) REFERENCE_GENOME="$2"; shift 2 ;;
            --gtf) REFERENCE_GTF="$2"; shift 2 ;;
            -h|--help) show_usage; exit 0 ;;
            *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
        esac
    done
    
    if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
        log_error "Input and output directories are required"
        show_usage
        exit 1
    fi
}

# Check installation
check_installation() {
    log_info "Checking for installed tools..."
    
    if [ ! -f "$CONFIG_FILE" ]; then
        log_warning "Alignment tools not installed"
        log_info "Installation required"
        echo ""
        
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
    
    if [ -f "$CONFIG_FILE" ]; then
        source "$CONFIG_FILE"
        log_success "Configuration loaded"
    else
        log_error "Configuration file not found: $CONFIG_FILE"
        exit 1
    fi
}

# Prompt aligner selection
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
    
    CUSTOM_FASTA="${CUSTOM_FASTA/#\~/$HOME}"
    CUSTOM_GTF="${CUSTOM_GTF/#\~/$HOME}"
    
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
    
    local GENOME_BASENAME=$(basename "${REFERENCE_GENOME}" | sed 's/\.[^.]*$//')
    local REF_DIR="$(dirname "$REFERENCE_GENOME")"
    STAR_INDEX_DIR="${REF_DIR}/STAR_index_${GENOME_BASENAME}"
    HISAT2_INDEX_DIR="${REF_DIR}/HISAT2_index_${GENOME_BASENAME}"
    HISAT2_INDEX_PREFIX="${HISAT2_INDEX_DIR}/genome"
    
    log_success "Custom reference validated:"
    log_info "FASTA: $REFERENCE_GENOME"
    log_info "GTF: $REFERENCE_GTF"
}

# Smart reference selection - ALWAYS prompts
prompt_reference_selection() {
    echo ""
    echo "========================================"
    echo "  Reference Genome Selection"
    echo "========================================"
    echo ""
    
    local HAS_DEFAULT="no"
    local HAS_PREVIOUS="no"
    local PREV_FASTA=""
    local PREV_GTF=""
    
    if [ -n "$DEFAULT_FASTA" ] && [ -f "$DEFAULT_FASTA" ] && [ -n "$DEFAULT_GTF" ] && [ -f "$DEFAULT_GTF" ]; then
        HAS_DEFAULT="yes"
    fi
    
    if [ -f "$REF_CONFIG" ]; then
        source "$REF_CONFIG"
        if [ -n "$REFERENCE_FASTA" ] && [ -f "$REFERENCE_FASTA" ] && [ -n "$REFERENCE_GTF" ] && [ -f "$REFERENCE_GTF" ]; then
            if [ "$REFERENCE_FASTA" != "$DEFAULT_FASTA" ]; then
                HAS_PREVIOUS="yes"
                PREV_FASTA="$REFERENCE_FASTA"
                PREV_GTF="$REFERENCE_GTF"
            fi
        fi
    fi
    
    log_info "Select reference genome to use:"
    echo ""
    
    local option_num=1
    declare -A option_map
    
    if [ "$HAS_DEFAULT" = "yes" ]; then
        echo "  ${option_num}) Default: GRCh38 (Ensembl 113) - found at ${DEFAULT_REFERENCE_DIR}"
        option_map[$option_num]="default"
        option_num=$((option_num + 1))
    fi
    
    if [ "$HAS_PREVIOUS" = "yes" ]; then
        echo "  ${option_num}) Previous custom: $(basename $PREV_FASTA) and $(basename $PREV_GTF)"
        option_map[$option_num]="previous"
        option_num=$((option_num + 1))
    fi
    
    echo "  ${option_num}) New custom: Enter new genome paths"
    option_map[$option_num]="custom"
    
    echo ""
    local max_option=$((option_num))
    read -p "Enter choice [1-${max_option}]: " -n 1 -r CHOICE
    echo
    
    if [ -z "$CHOICE" ] || [ "$CHOICE" -lt 1 ] || [ "$CHOICE" -gt "$max_option" ]; then
        log_error "Invalid choice"
        exit 1
    fi
    
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
            
            if [ -f "$REF_CONFIG" ]; then
                source "$REF_CONFIG"
                STAR_INDEX_DIR="$STAR_INDEX_DIR"
                HISAT2_INDEX_DIR="$HISAT2_INDEX_DIR"
                HISAT2_INDEX_PREFIX="$HISAT2_INDEX_PREFIX"
            else
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

# Check reference genome
check_reference() {
    log_info "Checking reference genome..."
    
    if [ -n "$REFERENCE_GENOME" ] && [ -n "$REFERENCE_GTF" ]; then
        if [ ! -f "$REFERENCE_GENOME" ]; then
            log_error "Reference genome not found: $REFERENCE_GENOME"
            exit 1
        fi
        if [ ! -f "$REFERENCE_GTF" ]; then
            log_error "Reference GTF not found: $REFERENCE_GTF"
            exit 1
        fi
        
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
    
    prompt_reference_selection
    
    if [ ! -f "$REFERENCE_GENOME" ]; then
        log_error "Reference genome not found: $REFERENCE_GENOME"
        exit 1
    fi
    
    if [ ! -f "$REFERENCE_GTF" ]; then
        log_error "Reference GTF not found: $REFERENCE_GTF"
        exit 1
    fi
}

# Check and build STAR index
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
    
    read -t 30 -p "Build STAR index now? [Y/n] (auto-continue in 30s): " -n 1 -r REPLY || REPLY=""
    echo
    
    if [[ "$REPLY" =~ ^[Nn]$ ]]; then
        log_error "STAR index is required for alignment"
        exit 1
    fi
    
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
}

# Check and build HISAT2 index
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
    
    read -t 30 -p "Build HISAT2 index now? [Y/n] (auto-continue in 30s): " -n 1 -r REPLY || REPLY=""
    echo
    
    if [[ "$REPLY" =~ ^[Nn]$ ]]; then
        log_error "HISAT2 index is required for alignment"
        exit 1
    fi
    
    mkdir -p "$HISAT2_INDEX_DIR"
    
    log_info "Running HISAT2 genome index build..."
    "$HISAT2_BUILD_BIN" -p "$THREADS" "$REFERENCE_GENOME" "$INDEX_BASE" 2>&1 | tee "${SCRIPT_DIR}/hisat2_index_build.log"
    
    if [ -f "${INDEX_BASE}.1.ht2" ]; then
        log_success "HISAT2 index built successfully"
    else
        log_error "HISAT2 index build failed"
        exit 1
    fi
}

# Detect read layout
detect_read_layout() {
    log_info "Detecting read layout..."
    
    local FIRST_FILE=$(find "$INPUT_DIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) 2>/dev/null | head -n 1)
    
    if [ -z "$FIRST_FILE" ]; then
        log_error "No FASTQ files found in $INPUT_DIR"
        log_info "Expected patterns: *.fastq.gz, *.fq.gz, *.fastq, *.fq"
        exit 1
    fi
    
    if [[ "$FIRST_FILE" =~ _R1[_.]|_1[_.]|\.1[_.]|_F[_.]|_R1_|_1\.|_paired ]]; then
        PAIRED_END=true
        log_success "Detected: Paired-end reads"
    else
        PAIRED_END=false
        log_warning "Detected: Single-end reads"
        
        # Prompt for single-end continuation with 30s timeout
        echo ""
        read -t 30 -p "Continue with single-end alignment? [Y/n] (auto-continue in 30s): " -n 1 -r REPLY || REPLY=""
        echo
        
        if [[ "$REPLY" =~ ^[Nn]$ ]]; then
            log_error "Alignment cancelled by user"
            exit 1
        fi
        
        log_info "Proceeding with single-end alignment..."
    fi
}

# Helper function to find R2 file
find_r2_file() {
    local r1_file=$1
    local basename=$(basename "$r1_file")
    local dir=$(dirname "$r1_file")
    
    # Try different R2 patterns
    local r2_patterns=(
        "${basename/_R1/_R2}"
        "${basename/_r1/_r2}"
        "${basename/_1/_2}"
        "${basename/.1./.2.}"
    )
    
    for pattern in "${r2_patterns[@]}"; do
        local r2_candidate="${dir}/${pattern}"
        if [ -f "$r2_candidate" ]; then
            echo "$r2_candidate"
            return 0
        fi
    done
    
    return 1
}

# Detect pattern used
detect_pattern() {
    local file=$1
    local basename=$(basename "$file")
    
    if [[ "$basename" =~ _R1 ]]; then
        echo "_R1/_R2"
    elif [[ "$basename" =~ _r1 ]]; then
        echo "_r1/_r2"
    elif [[ "$basename" =~ _1\. ]]; then
        echo "_1./_2."
    elif [[ "$basename" =~ \.1\. ]]; then
        echo ".1./.2."
    elif [[ "$basename" =~ _paired ]]; then
        echo "_paired"
    else
        echo "unknown"
    fi
}

# Create output directories
create_output_dirs() {
    log_info "Creating output directory structure..."
    
    mkdir -p "${OUTPUT_DIR}/logs"
    mkdir -p "${OUTPUT_DIR}/temp"
    
    if [ "$ALIGNER" = "both" ]; then
        mkdir -p "${OUTPUT_DIR}/STAR"
        mkdir -p "${OUTPUT_DIR}/HISAT2"
    fi
    
    log_success "Output directories created"
}

# STAR alignment for paired-end (with your exact parameters)
align_star_pe() {
    local sample_id=$1
    local r1_file=$2
    local r2_file=$3
    local subdir=$4
    
    log_info "Aligning sample: ${sample_id} (STAR)"
    
    local base_output="${OUTPUT_DIR}"
    [ -n "$subdir" ] && base_output="${OUTPUT_DIR}/${subdir}"
    
    local output_prefix="${base_output}/${sample_id}/${sample_id}_"
    mkdir -p "${base_output}/${sample_id}"
    
    "$STAR_BIN" \
        --runThreadN "$THREADS" \
        --limitBAMsortRAM 10000000000 \
        --runMode alignReads \
        --genomeDir "$STAR_INDEX_DIR" \
        --readFilesIn "$r1_file" "$r2_file" \
        --sjdbGTFfile "$REFERENCE_GTF" \
        --sjdbOverhang 100 \
        --alignSJDBoverhangMin 1 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "$output_prefix" \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFilterMatchNminOverLread 0.66 \
        --outFilterScoreMinOverLread 0.66 \
        --outFilterType BySJout \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 12 \
        --outSAMstrandField intronMotif \
        --outSAMattributes NH HI NM MD AS XS \
        --outSAMunmapped Within \
        --outReadsUnmapped Fastx \
        --outSAMheaderHD @HD VN:1.4 \
        2>&1 | tee "${OUTPUT_DIR}/logs/${sample_id}_star.log"
    
    if [ -f "${output_prefix}Aligned.sortedByCoord.out.bam" ]; then
        log_info "Indexing BAM file..."
        samtools index "${output_prefix}Aligned.sortedByCoord.out.bam"
        log_success "STAR alignment completed: ${sample_id}"
    else
        log_error "STAR alignment failed: ${sample_id}"
        return 1
    fi
}

# STAR alignment for single-end
align_star_se() {
    local sample_id=$1
    local r1_file=$2
    local subdir=$3
    
    log_info "Aligning sample: ${sample_id} (STAR - single-end)"
    
    local base_output="${OUTPUT_DIR}"
    [ -n "$subdir" ] && base_output="${OUTPUT_DIR}/${subdir}"
    
    local output_prefix="${base_output}/${sample_id}/${sample_id}_"
    mkdir -p "${base_output}/${sample_id}"
    
    "$STAR_BIN" \
        --runThreadN "$THREADS" \
        --limitBAMsortRAM 10000000000 \
        --runMode alignReads \
        --genomeDir "$STAR_INDEX_DIR" \
        --readFilesIn "$r1_file" \
        --sjdbGTFfile "$REFERENCE_GTF" \
        --sjdbOverhang 100 \
        --alignSJDBoverhangMin 1 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "$output_prefix" \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFilterMatchNminOverLread 0.66 \
        --outFilterScoreMinOverLread 0.66 \
        --outFilterType BySJout \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 12 \
        --outSAMstrandField intronMotif \
        --outSAMattributes NH HI NM MD AS XS \
        --outSAMunmapped Within \
        --outReadsUnmapped Fastx \
        --outSAMheaderHD @HD VN:1.4 \
        2>&1 | tee "${OUTPUT_DIR}/logs/${sample_id}_star.log"
    
    if [ -f "${output_prefix}Aligned.sortedByCoord.out.bam" ]; then
        log_info "Indexing BAM file..."
        samtools index "${output_prefix}Aligned.sortedByCoord.out.bam"
        log_success "STAR alignment completed: ${sample_id}"
    else
        log_error "STAR alignment failed: ${sample_id}"
        return 1
    fi
}

# HISAT2 alignment for paired-end
align_hisat2_pe() {
    local sample_id=$1
    local r1_file=$2
    local r2_file=$3
    local subdir=$4
    
    log_info "Aligning sample: ${sample_id} (HISAT2)"
    
    local base_output="${OUTPUT_DIR}"
    [ -n "$subdir" ] && base_output="${OUTPUT_DIR}/${subdir}"
    
    mkdir -p "${base_output}/${sample_id}"
    local output_sam="${base_output}/${sample_id}/${sample_id}.sam"
    local output_bam="${base_output}/${sample_id}/${sample_id}.bam"
    local output_sorted="${base_output}/${sample_id}/${sample_id}_Aligned.sortedByCoord.out.bam"
    
    "$HISAT2_BIN" \
        -q \
        -x "$HISAT2_INDEX_PREFIX" \
        -1 "$r1_file" \
        -2 "$r2_file" \
        --phred33 \
        --summary-file "${OUTPUT_DIR}/logs/${sample_id}_hisat2_summary.txt" \
        --met-file "${OUTPUT_DIR}/logs/${sample_id}_hisat2_metrics.txt" \
        -p "$THREADS" \
        --reorder \
        -S "$output_sam" \
        2>&1 | tee "${OUTPUT_DIR}/logs/${sample_id}_hisat2.log"
    
    if [ $? -eq 0 ] && [ -f "$output_sam" ]; then
        log_info "Converting and sorting BAM..."
        samtools view -bS "$output_sam" > "$output_bam"
        samtools sort "$output_bam" -o "$output_sorted"
        samtools index "$output_sorted"
        rm -f "$output_sam" "$output_bam"
        
        log_success "HISAT2 alignment completed: ${sample_id}"
    else
        log_error "HISAT2 alignment failed: ${sample_id}"
        return 1
    fi
}

# HISAT2 alignment for single-end
align_hisat2_se() {
    local sample_id=$1
    local r1_file=$2
    local subdir=$3
    
    log_info "Aligning sample: ${sample_id} (HISAT2 - single-end)"
    
    local base_output="${OUTPUT_DIR}"
    [ -n "$subdir" ] && base_output="${OUTPUT_DIR}/${subdir}"
    
    mkdir -p "${base_output}/${sample_id}"
    local output_sam="${base_output}/${sample_id}/${sample_id}.sam"
    local output_bam="${base_output}/${sample_id}/${sample_id}.bam"
    local output_sorted="${base_output}/${sample_id}/${sample_id}_Aligned.sortedByCoord.out.bam"
    
    "$HISAT2_BIN" \
        -q \
        -x "$HISAT2_INDEX_PREFIX" \
        -U "$r1_file" \
        --phred33 \
        --summary-file "${OUTPUT_DIR}/logs/${sample_id}_hisat2_summary.txt" \
        --met-file "${OUTPUT_DIR}/logs/${sample_id}_hisat2_metrics.txt" \
        -p "$THREADS" \
        --reorder \
        -S "$output_sam" \
        2>&1 | tee "${OUTPUT_DIR}/logs/${sample_id}_hisat2.log"
    
    if [ $? -eq 0 ] && [ -f "$output_sam" ]; then
        log_info "Converting and sorting BAM..."
        samtools view -bS "$output_sam" > "$output_bam"
        samtools sort "$output_bam" -o "$output_sorted"
        samtools index "$output_sorted"
        rm -f "$output_sam" "$output_bam"
        
        log_success "HISAT2 alignment completed: ${sample_id}"
    else
        log_error "HISAT2 alignment failed: ${sample_id}"
        return 1
    fi
}

# Process all samples in directory
process_all_samples() {
    log_info "Processing all samples in: $INPUT_DIR"
    
    if [ "$PAIRED_END" = true ]; then
        local r1_files=$(find "$INPUT_DIR" -type f \( -name "*_R1*.fastq" -o -name "*_R1*.fq" -o -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*_1.fastq" -o -name "*_1.fq" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" -o -name "*_paired.fastq" -o -name "*_paired.fq" \) | sort)
        
        if [ -z "$r1_files" ]; then
            log_error "No paired-end FASTQ files found in $INPUT_DIR"
            exit 1
        fi
        
        local sample_count=0
        
        for r1_file in $r1_files; do
            sample_count=$((sample_count + 1))
            
            local basename=$(basename "$r1_file")
            local sample_id=$(echo "$basename" | sed -E 's/_(R1|r1|1|paired).*//g')
            
            local r2_file=$(find_r2_file "$r1_file")
            local pattern=$(detect_pattern "$r1_file")
            
            if [ ! -f "$r2_file" ]; then
                log_warning "R2 not found for $r1_file (pattern: ${pattern}), skipping..."
                continue
            fi
            
            echo ""
            log_info "========================================="
            log_info "Sample ${sample_count}: ${sample_id}"
            log_info "========================================="
            log_info "Pattern: ${pattern}"
            log_info "R1: $r1_file"
            log_info "R2: $r2_file"
            
            case "$ALIGNER" in
                star)
                    align_star_pe "$sample_id" "$r1_file" "$r2_file" ""
                    ;;
                hisat2)
                    align_hisat2_pe "$sample_id" "$r1_file" "$r2_file" ""
                    ;;
                both)
                    log_info "Running alignment with BOTH aligners..."
                    align_star_pe "$sample_id" "$r1_file" "$r2_file" "STAR"
                    align_hisat2_pe "$sample_id" "$r1_file" "$r2_file" "HISAT2"
                    ;;
            esac
        done
        
        log_success "Processed ${sample_count} samples"
    else
        # Single-end processing
        local se_files=$(find "$INPUT_DIR" -type f \( -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)
        
        if [ -z "$se_files" ]; then
            log_error "No FASTQ files found in $INPUT_DIR"
            exit 1
        fi
        
        local sample_count=0
        
        for se_file in $se_files; do
            sample_count=$((sample_count + 1))
            
            local basename=$(basename "$se_file")
            local sample_id=$(echo "$basename" | sed -E 's/\.(fastq|fq)(\.gz)?$//g')
            
            echo ""
            log_info "========================================="
            log_info "Sample ${sample_count}: ${sample_id}"
            log_info "========================================="
            log_info "File: $se_file"
            
            case "$ALIGNER" in
                star)
                    align_star_se "$sample_id" "$se_file" ""
                    ;;
                hisat2)
                    align_hisat2_se "$sample_id" "$se_file" ""
                    ;;
                both)
                    log_info "Running alignment with BOTH aligners..."
                    align_star_se "$sample_id" "$se_file" "STAR"
                    align_hisat2_se "$sample_id" "$se_file" "HISAT2"
                    ;;
            esac
        done
        
        log_success "Processed ${sample_count} samples"
    fi
}

# Display final summary
display_final_summary() {
    echo ""
    echo "========================================"
    echo "  Alignment Complete"
    echo "========================================"
    echo ""
    log_success "All alignments completed!"
    echo ""
    log_info "Aligner(s) used: ${ALIGNER^^}"
    log_info "Results: ${OUTPUT_DIR}"
    echo ""
    
    if [ "$ALIGNER" = "both" ]; then
        log_info "Output structure (BOTH aligners):"
        echo "  STAR results:   ${OUTPUT_DIR}/STAR/[sample]/"
        echo "  HISAT2 results: ${OUTPUT_DIR}/HISAT2/[sample]/"
    else
        log_info "Output structure:"
        echo "  BAM files: ${OUTPUT_DIR}/[sample]/*Aligned.sortedByCoord.out.bam"
        echo "  BAM index: ${OUTPUT_DIR}/[sample]/*Aligned.sortedByCoord.out.bam.bai"
    fi
    
    echo "  Logs:      ${OUTPUT_DIR}/logs/"
    echo ""
    log_info "Next steps:"
    echo "  - Visualize: IGV or samtools view"
    echo "  - Quantify: featureCounts, RSEM, or Salmon"
    echo "  - QC: samtools flagstat, Qualimap"
    echo ""
    echo "========================================"
}

# Main execution
main() {
    echo "========================================"
    echo "  Alignment Module"
    echo "========================================"
    echo ""
    
    check_installation
    
    if [ -z "$ALIGNER" ]; then
        prompt_aligner_selection
    fi
    
    check_aligner_availability
    check_reference
    
    if [[ "$ALIGNER" == "star" || "$ALIGNER" == "both" ]]; then
        check_star_index
    fi
    
    if [[ "$ALIGNER" == "hisat2" || "$ALIGNER" == "both" ]]; then
        check_hisat2_index
    fi
    
    detect_read_layout
    create_output_dirs
    
    log_info "Starting alignment..."
    process_all_samples
    
    display_final_summary
}

parse_args "$@"
main
