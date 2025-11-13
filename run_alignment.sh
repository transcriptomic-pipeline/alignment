#!/bin/bash

# Alignment Module - Main Execution Script
# Supports STAR and HISAT2 with system Python OR conda environment

set -euo pipefail

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
PAIRED_END=true

# Tool paths
STAR_BIN=""
HISAT2_BIN=""
HISAT2_BUILD_BIN=""
SAMTOOLS_BIN=""
STAR_INDEX_DIR=""
HISAT2_INDEX_DIR=""
HISAT2_INDEX_PREFIX=""

# Python/Conda settings (loaded from config)
USE_SYSTEM_PYTHON=""
PYTHON_BIN=""
CONDA_DIR=""
CONDA_ENV_NAME=""

# Timeout for auto-continue prompts
PROMPT_TIMEOUT=30

usage() {
    cat << EOF
Usage: $0 -i <input_dir> -o <output_dir> [options]

Required:
    -i, --input         Input directory with trimmed FASTQ files
    -o, --output        Output directory for aligned BAM files

Optional:
    --aligner           Aligner: star, hisat2, both (default: prompt)
    -s, --samples       Sample list file
    -t, --threads       Number of threads (default: 12)
    --reference         Reference FASTA file
    --gtf               Reference GTF file
    --single-end        Single-end reads (default: paired-end)
    -h, --help          Show this help

Examples:
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner star
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner both -t 20

EOF
    exit 1
}

command_exists() { command -v "$1" >/dev/null 2>&1; }

# Check if conda environment exists and can be activated (only if needed)
setup_python_for_hisat2() {
    # If using system Python, nothing to do
    if [ "$USE_SYSTEM_PYTHON" = "yes" ]; then
        log_info "Using system Python: $PYTHON_BIN"
        return 0
    fi
    
    # If using conda, verify environment
    if [ "$USE_SYSTEM_PYTHON" = "no" ]; then
        if [ -z "$CONDA_DIR" ] || [ ! -d "$CONDA_DIR" ]; then
            log_error "Conda directory not found: $CONDA_DIR"
            log_error "Run: ./install.sh --aligner hisat2"
            exit 1
        fi
        
        if [ ! -f "${CONDA_DIR}/etc/profile.d/conda.sh" ]; then
            log_error "Conda profile script not found"
            log_error "Run: ./install.sh --aligner hisat2"
            exit 1
        fi
        
        log_info "Using conda environment: $CONDA_ENV_NAME"
        return 0
    fi
    
    # Unknown configuration
    log_error "Invalid Python configuration in config file"
    exit 1
}

# Check installation
check_installation() {
    log_info "Checking for installed tools..."
    
    if [ ! -f "$CONFIG_FILE" ]; then
        log_error "Configuration not found"
        log_info "Please run: ./install.sh"
        exit 1
    fi
    
    source "$CONFIG_FILE"
    [ -f "$REF_CONFIG" ] && source "$REF_CONFIG"
    
    log_success "Configuration loaded"
}

# Prompt for aligner selection
prompt_aligner_selection() {
    echo ""
    echo "========================================"
    echo "  Aligner Selection"
    echo "========================================"
    echo ""
    
    local STAR_AVAILABLE=$([[ "$STAR_INSTALLED" == "yes" && -f "$STAR_BIN" ]] && echo "yes" || echo "no")
    local HISAT2_AVAILABLE=$([[ "$HISAT2_INSTALLED" == "yes" && -f "$HISAT2_BIN" ]] && echo "yes" || echo "no")
    
    if [ "$STAR_AVAILABLE" = "yes" ] && [ "$HISAT2_AVAILABLE" = "yes" ]; then
        echo "  1) STAR (high memory ~32GB)"
        echo "  2) HISAT2 (low memory ~8GB)"
        echo "  3) BOTH (for comparison)"
        read -p "Enter choice [1-3]: " ALIGNER_CHOICE
        case "$ALIGNER_CHOICE" in
            1) ALIGNER="star" ;;
            2) ALIGNER="hisat2" ;;
            3) ALIGNER="both" ;;
            *) log_error "Invalid choice"; exit 1 ;;
        esac
    elif [ "$STAR_AVAILABLE" = "yes" ]; then
        ALIGNER="star"
    elif [ "$HISAT2_AVAILABLE" = "yes" ]; then
        ALIGNER="hisat2"
    else
        log_error "No aligner installed"
        exit 1
    fi
    
    log_success "Selected aligner: ${ALIGNER^^}"
}

# Verify aligner availability
verify_aligner() {
    case "$ALIGNER" in
        star)
            [ "$STAR_INSTALLED" != "yes" ] || [ ! -f "$STAR_BIN" ] && { log_error "STAR not installed"; exit 1; }
            ;;
        hisat2)
            [ "$HISAT2_INSTALLED" != "yes" ] || [ ! -f "$HISAT2_BIN" ] && { log_error "HISAT2 not installed"; exit 1; }
            setup_python_for_hisat2
            ;;
        both)
            [ "$STAR_INSTALLED" != "yes" ] || [ ! -f "$STAR_BIN" ] && { log_error "STAR not installed"; exit 1; }
            [ "$HISAT2_INSTALLED" != "yes" ] || [ ! -f "$HISAT2_BIN" ] && { log_error "HISAT2 not installed"; exit 1; }
            setup_python_for_hisat2
            ;;
    esac
}

# Check reference genome
check_reference() {
    log_info "Checking reference genome..."
    
    [ -z "$REFERENCE_GENOME" ] && [ -n "$REFERENCE_FASTA" ] && [ -f "$REFERENCE_FASTA" ] && REFERENCE_GENOME="$REFERENCE_FASTA"
    [ -z "$REFERENCE_GTF" ] && [ -f "${REFERENCE_DIR}/Homo_sapiens.GRCh38.113.gtf" ] && REFERENCE_GTF="${REFERENCE_DIR}/Homo_sapiens.GRCh38.113.gtf"
    
    [ ! -f "$REFERENCE_GENOME" ] && { log_error "Reference genome not found"; exit 1; }
    [ ! -f "$REFERENCE_GTF" ] && { log_error "Reference GTF not found"; exit 1; }
    
    log_success "Reference genome: $REFERENCE_GENOME"
    log_success "Reference GTF: $REFERENCE_GTF"
}

# Check or build index
check_or_build_index() {
    case "$ALIGNER" in
        star) check_or_build_star_index ;;
        hisat2) check_or_build_hisat2_index ;;
        both) check_or_build_star_index; check_or_build_hisat2_index ;;
    esac
}

check_or_build_star_index() {
    log_info "Checking STAR genome index..."
    [ ! -d "$STAR_INDEX_DIR" ] || [ -z "$(ls -A $STAR_INDEX_DIR 2>/dev/null)" ] && build_star_index || log_success "STAR index found"
}

check_or_build_hisat2_index() {
    log_info "Checking HISAT2 genome index..."
    [ ! -f "${HISAT2_INDEX_PREFIX}.1.ht2" ] && build_hisat2_index || log_success "HISAT2 index found"
}

build_star_index() {
    log_info "Building STAR genome index (may take 30-60 min, needs ~32GB RAM)..."
    echo -ne "Build now? [Y/n] (auto-continue in ${PROMPT_TIMEOUT}s): "
    read -t $PROMPT_TIMEOUT -n 1 -r REPLY || REPLY="y"
    echo
    [[ $REPLY =~ ^[Nn]$ ]] && { log_error "STAR index required"; exit 1; }
    
    mkdir -p "$STAR_INDEX_DIR"
    "$STAR_BIN" --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir "$STAR_INDEX_DIR" \
        --genomeFastaFiles "$REFERENCE_GENOME" --sjdbGTFfile "$REFERENCE_GTF" --sjdbOverhang 100 \
        2>&1 | tee "${SCRIPT_DIR}/star_index_build.log"
    
    [ $? -eq 0 ] && log_success "STAR index built" || { log_error "STAR index build failed"; exit 1; }
}

build_hisat2_index() {
    log_info "Building HISAT2 genome index (may take 30-45 min, needs ~8GB RAM)..."
    echo -ne "Build now? [Y/n] (auto-continue in ${PROMPT_TIMEOUT}s): "
    read -t $PROMPT_TIMEOUT -n 1 -r REPLY || REPLY="y"
    echo
    [[ $REPLY =~ ^[Nn]$ ]] && { log_error "HISAT2 index required"; exit 1; }
    
    mkdir -p "$HISAT2_INDEX_DIR"
    
    # Use system Python OR conda environment
    if [ "$USE_SYSTEM_PYTHON" = "yes" ]; then
        log_info "Building HISAT2 index with system Python..."
        "$HISAT2_BUILD_BIN" -f "$REFERENCE_GENOME" -p "$THREADS" "$HISAT2_INDEX_PREFIX" \
            2>&1 | tee "${SCRIPT_DIR}/hisat2_index_build.log"
    else
        log_info "Building HISAT2 index with conda environment '$CONDA_ENV_NAME'..."
        set +u
        source "${CONDA_DIR}/etc/profile.d/conda.sh"
        conda activate "$CONDA_ENV_NAME" || { log_error "Failed to activate conda"; exit 1; }
        
        "$HISAT2_BUILD_BIN" -f "$REFERENCE_GENOME" -p "$THREADS" "$HISAT2_INDEX_PREFIX" \
            2>&1 | tee "${SCRIPT_DIR}/hisat2_index_build.log"
        
        local STATUS=$?
        conda deactivate
        set -u
        [ $STATUS -ne 0 ] && { log_error "HISAT2 index build failed"; exit 1; }
    fi
    
    log_success "HISAT2 index built"
}

create_output_dirs() {
    mkdir -p "${OUTPUT_DIR}/logs" "${OUTPUT_DIR}/temp"
    [ "$ALIGNER" = "both" ] && mkdir -p "${OUTPUT_DIR}/STAR" "${OUTPUT_DIR}/HISAT2"
}

detect_pattern() {
    local f=$1
    [[ "$f" =~ _R1|_r1 ]] && echo "_R1/_R2" || [[ "$f" =~ _1[\._] ]] && echo "_1/_2" || echo "unknown"
}

find_r2_file() {
    local r1=$1
    echo "$r1" | sed -E 's/_(R1|r1|1)([\._])/_\22\2/; s/_R1/_R2/; s/_r1/_r2/; s/_1\./_2./; s/_1_/_2_/'
}

align_star_pe() {
    local sid=$1 r1=$2 r2=$3 sub=$4
    local out="${OUTPUT_DIR}${sub:+/$sub}/${sid}/${sid}_"
    mkdir -p "$(dirname "$out")"
    
    log_info "Aligning $sid (STAR)..."
    "$STAR_BIN" --runThreadN "$THREADS" --runMode alignReads --genomeDir "$STAR_INDEX_DIR" \
        --readFilesIn "$r1" "$r2" --sjdbGTFfile "$REFERENCE_GTF" --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "$out" --outFilterMultimapNmax 20 --alignIntronMax 1000000 \
        2>&1 | tee "${OUTPUT_DIR}/logs/${sid}_star.log"
    
    [ -f "${out}Aligned.sortedByCoord.out.bam" ] && samtools index "${out}Aligned.sortedByCoord.out.bam" && log_success "STAR: $sid" || log_error "STAR failed: $sid"
}

align_hisat2_pe() {
    local sid=$1 r1=$2 r2=$3 sub=$4
    local base="${OUTPUT_DIR}${sub:+/$sub}/${sid}"
    mkdir -p "$base"
    
    log_info "Aligning $sid (HISAT2)..."
    "$HISAT2_BIN" -q -x "$HISAT2_INDEX_PREFIX" -1 "$r1" -2 "$r2" -p "$THREADS" \
        -S "${base}/${sid}.sam" 2>&1 | tee "${OUTPUT_DIR}/logs/${sid}_hisat2.log"
    
    samtools view -bS "${base}/${sid}.sam" | samtools sort -o "${base}/${sid}Aligned.sortedByCoord.out.bam"
    samtools index "${base}/${sid}Aligned.sortedByCoord.out.bam"
    rm -f "${base}/${sid}.sam"
    log_success "HISAT2: $sid"
}

process_all_samples() {
    local r1_files=$(find "$INPUT_DIR" -type f \( -name "*_R1*.f*q*" -o -name "*_1.f*q*" \) | sort)
    [ -z "$r1_files" ] && { log_error "No FASTQ files found"; exit 1; }
    
    local count=0
    for r1 in $r1_files; do
        count=$((count + 1))
        local sid=$(basename "$r1" | sed -E 's/_(R1|r1|1).*//g')
        local r2=$(find_r2_file "$r1")
        [ ! -f "$r2" ] && { log_warning "R2 not found for $sid"; continue; }
        
        echo ""
        log_info "Sample $count: $sid"
        
        case "$ALIGNER" in
            star) align_star_pe "$sid" "$r1" "$r2" "" ;;
            hisat2) align_hisat2_pe "$sid" "$r1" "$r2" "" ;;
            both) align_star_pe "$sid" "$r1" "$r2" "STAR"; align_hisat2_pe "$sid" "$r1" "$r2" "HISAT2" ;;
        esac
    done
}

generate_summary() {
    cat > "${OUTPUT_DIR}/alignment_summary.txt" << EOF
Alignment Module Summary
Date: $(date)
Aligner: ${ALIGNER^^}
Reference: ${REFERENCE_GENOME}
Threads: ${THREADS}
Results: ${OUTPUT_DIR}
EOF
    log_success "Summary: ${OUTPUT_DIR}/alignment_summary.txt"
}

parse_arguments() {
    [ $# -eq 0 ] && usage
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input) INPUT_DIR="$2"; shift 2 ;;
            -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
            -s|--samples) SAMPLE_FILE="$2"; shift 2 ;;
            -t|--threads) THREADS="$2"; shift 2 ;;
            --aligner) ALIGNER="$2"; shift 2 ;;
            --reference) REFERENCE_GENOME="$2"; shift 2 ;;
            --gtf) REFERENCE_GTF="$2"; shift 2 ;;
            --single-end) PAIRED_END=false; shift ;;
            -h|--help) usage ;;
            *) log_error "Unknown option: $1"; usage ;;
        esac
    done
    
    [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] && { log_error "Input and output required"; usage; }
    [ ! -d "$INPUT_DIR" ] && { log_error "Input directory not found: $INPUT_DIR"; exit 1; }
}

main() {
    echo "========================================"
    echo "  Alignment Module"
    echo "========================================"
    echo ""
    
    parse_arguments "$@"
    check_installation
    
    [ -z "$ALIGNER" ] && prompt_aligner_selection
    verify_aligner
    check_reference
    check_or_build_index
    create_output_dirs
    
    local start=$(date +%s)
    log_info "Started: $(date)"
    
    [ -n "$SAMPLE_FILE" ] && log_error "Sample file mode not yet implemented" || process_all_samples
    
    generate_summary
    
    local end=$(date +%s)
    local runtime=$((end - start))
    log_info "Completed in ${runtime}s"
    
    echo ""
    echo "========================================"
    echo "  Alignment Complete"
    echo "========================================"
    log_success "Results: ${OUTPUT_DIR}"
    [ "$ALIGNER" = "both" ] && echo "  STAR:   ${OUTPUT_DIR}/STAR/" && echo "  HISAT2: ${OUTPUT_DIR}/HISAT2/"
    echo "========================================"
}

main "$@"
