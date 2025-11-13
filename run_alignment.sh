#!/bin/bash

# Alignment Module - Main Execution Script
# Works standalone without requiring install.sh to be run first

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

log_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Default parameters
INPUT_DIR=""
OUTPUT_DIR=""
THREADS=12
ALIGNER=""
REFERENCE_GENOME=""
REFERENCE_GTF=""

# Tool paths (auto-detected)
STAR_BIN=""
HISAT2_BIN=""
HISAT2_BUILD_BIN=""

# Reference paths
STAR_INDEX_DIR=""
HISAT2_INDEX_DIR=""
HISAT2_INDEX_PREFIX=""

PROMPT_TIMEOUT=30

usage() {
    cat << EOF
Usage: $0 -i <input_dir> -o <output_dir> [options]

Required:
    -i, --input         Input directory with trimmed FASTQ files
    -o, --output        Output directory for aligned BAM files

Optional:
    --aligner           Aligner: star, hisat2, both (default: prompt)
    -t, --threads       Number of threads (default: 12)
    --reference         Reference FASTA file
    --gtf               Reference GTF file
    -h, --help          Show this help

Examples:
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner star -t 20
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner hisat2

EOF
    exit 1
}

command_exists() { command -v "$1" >/dev/null 2>&1; }

# Auto-detect installed tools
detect_tools() {
    log_info "Detecting installed alignment tools..."
    
    # Check for STAR
    if command_exists STAR; then
        STAR_BIN=$(command -v STAR)
        log_success "STAR found: $STAR_BIN"
    elif [ -f "${HOME}/softwares/bin/STAR" ]; then
        STAR_BIN="${HOME}/softwares/bin/STAR"
        log_success "STAR found: $STAR_BIN"
    fi
    
    # Check for HISAT2
    if command_exists hisat2; then
        HISAT2_BIN=$(command -v hisat2)
        HISAT2_BUILD_BIN=$(command -v hisat2-build)
        log_success "HISAT2 found: $HISAT2_BIN"
    elif [ -f "${HOME}/softwares/bin/hisat2" ]; then
        HISAT2_BIN="${HOME}/softwares/bin/hisat2"
        HISAT2_BUILD_BIN="${HOME}/softwares/bin/hisat2-build"
        log_success "HISAT2 found: $HISAT2_BIN"
    fi
    
    # Check if at least one aligner is available
    if [ -z "$STAR_BIN" ] && [ -z "$HISAT2_BIN" ]; then
        log_error "No aligner found (STAR or HISAT2)"
        log_info "Please install with: ./install.sh"
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
    
    if [ -n "$STAR_BIN" ] && [ -n "$HISAT2_BIN" ]; then
        echo "  1) STAR (high memory ~32GB)"
        echo "  2) HISAT2 (low memory ~8GB)"
        echo "  3) BOTH (for comparison)"
        echo ""
        read -p "Enter choice [1-3]: " ALIGNER_CHOICE
        
        case "$ALIGNER_CHOICE" in
            1) ALIGNER="star" ;;
            2) ALIGNER="hisat2" ;;
            3) ALIGNER="both" ;;
            *) log_error "Invalid choice"; exit 1 ;;
        esac
    elif [ -n "$STAR_BIN" ]; then
        ALIGNER="star"
        log_info "Using STAR (only installed aligner)"
    elif [ -n "$HISAT2_BIN" ]; then
        ALIGNER="hisat2"
        log_info "Using HISAT2 (only installed aligner)"
    fi
    
    log_success "Selected aligner: ${ALIGNER^^}"
}

# Auto-detect reference files
detect_reference() {
    log_info "Looking for reference genome..."
    
    local POSSIBLE_REF_DIRS=(
        "${HOME}/references/GRCh38_ensembl113"
        "${HOME}/Desktop/pipeline/references/GRCh38_ensembl113"
        "${SCRIPT_DIR}/../references/GRCh38_ensembl113"
    )
    
    for ref_dir in "${POSSIBLE_REF_DIRS[@]}"; do
        if [ -f "${ref_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
            REFERENCE_GENOME="${ref_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
            REFERENCE_GTF="${ref_dir}/Homo_sapiens.GRCh38.113.gtf"
            STAR_INDEX_DIR="${ref_dir}/STAR_index"
            HISAT2_INDEX_DIR="${ref_dir}/HISAT2_index"
            HISAT2_INDEX_PREFIX="${HISAT2_INDEX_DIR}/genome"
            
            log_success "Reference found: $REFERENCE_GENOME"
            log_success "GTF found: $REFERENCE_GTF"
            return 0
        fi
    done
    
    log_error "Reference genome not found"
    log_info "Please download with: ./install.sh --download-reference"
    exit 1
}

# Check or build index
check_or_build_index() {
    case "$ALIGNER" in
        star) check_or_build_star_index ;;
        hisat2) check_or_build_hisat2_index ;;
        both)
            check_or_build_star_index
            check_or_build_hisat2_index
            ;;
    esac
}

check_or_build_star_index() {
    log_info "Checking STAR genome index..."
    
    if [ ! -d "$STAR_INDEX_DIR" ] || [ -z "$(ls -A $STAR_INDEX_DIR 2>/dev/null)" ]; then
        log_warning "STAR index not found"
        build_star_index
    else
        log_success "STAR index found"
    fi
}

check_or_build_hisat2_index() {
    log_info "Checking HISAT2 genome index..."
    
    if [ ! -f "${HISAT2_INDEX_PREFIX}.1.ht2" ]; then
        log_warning "HISAT2 index not found"
        build_hisat2_index
    else
        log_success "HISAT2 index found"
    fi
}

build_star_index() {
    log_info "Building STAR genome index (30-60 min, ~32GB RAM)..."
    echo -ne "Build now? [Y/n] (auto in ${PROMPT_TIMEOUT}s): "
    read -t $PROMPT_TIMEOUT -n 1 -r REPLY || REPLY="y"
    echo
    
    [[ $REPLY =~ ^[Nn]$ ]] && { log_error "Index required"; exit 1; }
    
    mkdir -p "$STAR_INDEX_DIR"
    "$STAR_BIN" --runThreadN "$THREADS" --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX_DIR" --genomeFastaFiles "$REFERENCE_GENOME" \
        --sjdbGTFfile "$REFERENCE_GTF" --sjdbOverhang 100
    
    log_success "STAR index built"
}

build_hisat2_index() {
    log_info "Building HISAT2 genome index (30-45 min, ~8GB RAM)..."
    echo -ne "Build now? [Y/n] (auto in ${PROMPT_TIMEOUT}s): "
    read -t $PROMPT_TIMEOUT -n 1 -r REPLY || REPLY="y"
    echo
    
    [[ $REPLY =~ ^[Nn]$ ]] && { log_error "Index required"; exit 1; }
    
    mkdir -p "$HISAT2_INDEX_DIR"
    "$HISAT2_BUILD_BIN" -f "$REFERENCE_GENOME" -p "$THREADS" "$HISAT2_INDEX_PREFIX"
    
    log_success "HISAT2 index built"
}

find_r2_file() {
    local r1=$1
    echo "$r1" | sed -E 's/_(R1|r1|1)([\._])/_\22\2/; s/_R1/_R2/; s/_r1/_r2/; s/_1\./_2./; s/_1_/_2_/'
}

align_star_pe() {
    local sid=$1 r1=$2 r2=$3 sub=$4
    log_info "Aligning $sid (STAR)..."
    
    local base="${OUTPUT_DIR}${sub:+/$sub}"
    mkdir -p "${base}/${sid}"
    local out="${base}/${sid}/${sid}_"
    
    "$STAR_BIN" --runThreadN "$THREADS" --runMode alignReads \
        --genomeDir "$STAR_INDEX_DIR" --readFilesIn "$r1" "$r2" \
        --sjdbGTFfile "$REFERENCE_GTF" --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "$out" --outFilterMultimapNmax 20 \
        --alignIntronMax 1000000 2>&1 | tee "${OUTPUT_DIR}/logs/${sid}_star.log"
    
    [ -f "${out}Aligned.sortedByCoord.out.bam" ] && samtools index "${out}Aligned.sortedByCoord.out.bam" && log_success "STAR: $sid"
}

align_hisat2_pe() {
    local sid=$1 r1=$2 r2=$3 sub=$4
    log_info "Aligning $sid (HISAT2)..."
    
    local base="${OUTPUT_DIR}${sub:+/$sub}/${sid}"
    mkdir -p "$base"
    
    "$HISAT2_BIN" -q -x "$HISAT2_INDEX_PREFIX" -1 "$r1" -2 "$r2" -p "$THREADS" \
        -S "${base}/${sid}.sam" 2>&1 | tee "${OUTPUT_DIR}/logs/${sid}_hisat2.log"
    
    samtools view -bS "${base}/${sid}.sam" | samtools sort -o "${base}/${sid}Aligned.sortedByCoord.out.bam"
    samtools index "${base}/${sid}Aligned.sortedByCoord.out.bam"
    rm -f "${base}/${sid}.sam"
    log_success "HISAT2: $sid"
}

process_all_samples() {
    log_info "Scanning for samples in: $INPUT_DIR"
    
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
            both)
                align_star_pe "$sid" "$r1" "$r2" "STAR"
                align_hisat2_pe "$sid" "$r1" "$r2" "HISAT2"
                ;;
        esac
    done
}

parse_arguments() {
    [ $# -eq 0 ] && usage
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input) INPUT_DIR="$2"; shift 2 ;;
            -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
            -t|--threads) THREADS="$2"; shift 2 ;;
            --aligner) ALIGNER="$2"; shift 2 ;;
            --reference) REFERENCE_GENOME="$2"; shift 2 ;;
            --gtf) REFERENCE_GTF="$2"; shift 2 ;;
            -h|--help) usage ;;
            *) log_error "Unknown: $1"; usage ;;
        esac
    done
    
    [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] && { log_error "Input/output required"; usage; }
    [ ! -d "$INPUT_DIR" ] && { log_error "Input not found: $INPUT_DIR"; exit 1; }
}

main() {
    echo "========================================"
    echo "  Alignment Module"
    echo "========================================"
    echo ""
    
    parse_arguments "$@"
    detect_tools
    
    [ -z "$ALIGNER" ] && prompt_aligner_selection
    [ -z "$REFERENCE_GENOME" ] && detect_reference
    
    check_or_build_index
    mkdir -p "${OUTPUT_DIR}/logs"
    [ "$ALIGNER" = "both" ] && mkdir -p "${OUTPUT_DIR}/STAR" "${OUTPUT_DIR}/HISAT2"
    
    local start=$(date +%s)
    log_info "Started: $(date)"
    
    process_all_samples
    
    local end=$(date +%s)
    log_info "Completed in $((end - start))s"
    
    echo ""
    echo "========================================"
    echo "  Alignment Complete"
    echo "========================================"
    log_success "Results: ${OUTPUT_DIR}"
    echo "========================================"
}

main "$@"
