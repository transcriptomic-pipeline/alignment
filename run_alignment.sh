#!/bin/bash

# Alignment Module - Main Execution Script
# Supports STAR and HISAT2 without Python dependency issues

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

# Python settings from config
USE_SYSTEM_PYTHON=""
PYTHON_BIN=""

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
    $0 -i qc_results/trimmed/PE -o alignment_results
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner star -t 20
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner both

EOF
    exit 1
}

command_exists() { command -v "$1" >/dev/null 2>&1; }

check_installation() {
    log_info "Checking for installed tools..."
    
    if [ ! -f "$CONFIG_FILE" ]; then
        log_warning "Configuration not found. Running from outside alignment directory or not installed yet."
        
        # Try to find config in common locations
        if [ -f "${HOME}/softwares/alignment/config/install_paths.conf" ]; then
            CONFIG_FILE="${HOME}/softwares/alignment/config/install_paths.conf"
            REF_CONFIG="${HOME}/softwares/alignment/config/reference_paths.conf"
            log_info "Found configuration at: $CONFIG_FILE"
        else
            log_error "Configuration not found. Please run: ./install.sh"
            exit 1
        fi
    fi
    
    source "$CONFIG_FILE"
    [ -f "$REF_CONFIG" ] && source "$REF_CONFIG"
    
    log_success "Configuration loaded"
    
    # Log Python configuration (informational only - no conda enforcement)
    if [ "$USE_SYSTEM_PYTHON" = "yes" ]; then
        log_info "Python: System ($PYTHON_BIN)"
    else
        log_info "Python: Not configured (not required for alignment)"
    fi
}

prompt_aligner_selection() {
    echo ""
    echo "========================================"
    echo "  Aligner Selection"
    echo "========================================"
    echo ""
    
    local STAR_AVAILABLE=$([[ "$STAR_INSTALLED" == "yes" && -f "$STAR_BIN" ]] && echo "yes" || echo "no")
    local HISAT2_AVAILABLE=$([[ "$HISAT2_INSTALLED" == "yes" && -f "$HISAT2_BIN" ]] && echo "yes" || echo "no")
    
    if [ "$STAR_AVAILABLE" = "yes" ] && [ "$HISAT2_AVAILABLE" = "yes" ]; then
        log_info "Both STAR and HISAT2 are installed"
        echo ""
        echo "  1) STAR (high memory ~32GB)"
        echo "  2) HISAT2 (low memory ~8GB)"
        echo "  3) BOTH (run with both aligners)"
        echo ""
        read -p "Enter choice [1-3]: " ALIGNER_CHOICE
        
        case "$ALIGNER_CHOICE" in
            1) ALIGNER="star" ;;
            2) ALIGNER="hisat2" ;;
            3) ALIGNER="both" ;;
            *) log_error "Invalid choice"; exit 1 ;;
        esac
    elif [ "$STAR_AVAILABLE" = "yes" ]; then
        log_info "Using STAR (only installed aligner)"
        ALIGNER="star"
    elif [ "$HISAT2_AVAILABLE" = "yes" ]; then
        log_info "Using HISAT2 (only installed aligner)"
        ALIGNER="hisat2"
    else
        log_error "No aligner installed"
        exit 1
    fi
    
    log_success "Selected aligner: ${ALIGNER^^}"
}

verify_aligner() {
    log_info "Verifying aligner: ${ALIGNER^^}"
    
    case "$ALIGNER" in
        star)
            if [ "$STAR_INSTALLED" != "yes" ] || [ ! -f "$STAR_BIN" ]; then
                log_error "STAR not installed"
                exit 1
            fi
            log_success "STAR is available"
            ;;
        hisat2)
            if [ "$HISAT2_INSTALLED" != "yes" ] || [ ! -f "$HISAT2_BIN" ]; then
                log_error "HISAT2 not installed"
                exit 1
            fi
            log_success "HISAT2 is available"
            ;;
        both)
            if [ "$STAR_INSTALLED" != "yes" ] || [ ! -f "$STAR_BIN" ]; then
                log_error "STAR not installed"
                exit 1
            fi
            if [ "$HISAT2_INSTALLED" != "yes" ] || [ ! -f "$HISAT2_BIN" ]; then
                log_error "HISAT2 not installed"
                exit 1
            fi
            log_success "Both STAR and HISAT2 are available"
            ;;
    esac
}

check_reference() {
    log_info "Checking reference genome..."
    
    [ -z "$REFERENCE_GENOME" ] && [ -n "$REFERENCE_FASTA" ] && [ -f "$REFERENCE_FASTA" ] && REFERENCE_GENOME="$REFERENCE_FASTA"
    [ -z "$REFERENCE_GTF" ] && [ -n "${REFERENCE_DIR}" ] && [ -f "${REFERENCE_DIR}/Homo_sapiens.GRCh38.113.gtf" ] && REFERENCE_GTF="${REFERENCE_DIR}/Homo_sapiens.GRCh38.113.gtf"
    
    if [ ! -f "$REFERENCE_GENOME" ]; then
        log_error "Reference genome not found: $REFERENCE_GENOME"
        exit 1
    fi
    
    if [ ! -f "$REFERENCE_GTF" ]; then
        log_error "Reference GTF not found: $REFERENCE_GTF"
        exit 1
    fi
    
    log_success "Reference genome: $REFERENCE_GENOME"
    log_success "Reference GTF: $REFERENCE_GTF"
}

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
    log_info "Building STAR genome index..."
    log_warning "This may take 30-60 minutes and requires ~32 GB RAM"
    
    echo -ne "Build STAR index now? [Y/n] (auto-continue in ${PROMPT_TIMEOUT}s): "
    read -t $PROMPT_TIMEOUT -n 1 -r REPLY || REPLY="y"
    echo
    
    if [[ $REPLY =~ ^[Nn]$ ]]; then
        log_error "STAR index required"
        exit 1
    fi
    
    mkdir -p "$STAR_INDEX_DIR"
    
    "$STAR_BIN" \
        --runThreadN "$THREADS" \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX_DIR" \
        --genomeFastaFiles "$REFERENCE_GENOME" \
        --sjdbGTFfile "$REFERENCE_GTF" \
        --sjdbOverhang 100 \
        2>&1 | tee "${SCRIPT_DIR}/star_index_build.log"
    
    if [ $? -eq 0 ]; then
        log_success "STAR index built"
    else
        log_error "STAR index build failed"
        exit 1
    fi
}

build_hisat2_index() {
    log_info "Building HISAT2 genome index..."
    log_warning "This may take 30-45 minutes and requires ~8 GB RAM"
    
    echo -ne "Build HISAT2 index now? [Y/n] (auto-continue in ${PROMPT_TIMEOUT}s): "
    read -t $PROMPT_TIMEOUT -n 1 -r REPLY || REPLY="y"
    echo
    
    if [[ $REPLY =~ ^[Nn]$ ]]; then
        log_error "HISAT2 index required"
        exit 1
    fi
    
    mkdir -p "$HISAT2_INDEX_DIR"
    
    # HISAT2 doesn't require Python for indexing
    "$HISAT2_BUILD_BIN" -f "$REFERENCE_GENOME" -p "$THREADS" "$HISAT2_INDEX_PREFIX" \
        2>&1 | tee "${SCRIPT_DIR}/hisat2_index_build.log"
    
    if [ $? -eq 0 ]; then
        log_success "HISAT2 index built"
    else
        log_error "HISAT2 index build failed"
        exit 1
    fi
}

create_output_dirs() {
    log_info "Creating output directories..."
    mkdir -p "${OUTPUT_DIR}/logs" "${OUTPUT_DIR}/temp"
    
    if [ "$ALIGNER" = "both" ]; then
        mkdir -p "${OUTPUT_DIR}/STAR" "${OUTPUT_DIR}/HISAT2"
    fi
    
    log_success "Output directories created"
}

find_r2_file() {
    local r1=$1
    local r2=""
    
    if [[ "$r1" =~ _R1 ]]; then
        r2=$(echo "$r1" | sed 's/_R1/_R2/g')
    elif [[ "$r1" =~ _r1 ]]; then
        r2=$(echo "$r1" | sed 's/_r1/_r2/g')
    elif [[ "$r1" =~ _1\. ]]; then
        r2=$(echo "$r1" | sed 's/_1\./_2./g')
    elif [[ "$r1" =~ _1_ ]]; then
        r2=$(echo "$r1" | sed 's/_1_/_2_/g')
    fi
    
    echo "$r2"
}

align_star_pe() {
    local sample_id=$1
    local r1_file=$2
    local r2_file=$3
    local subdir=$4
    
    log_info "Aligning: ${sample_id} (STAR)"
    
    local base_output="${OUTPUT_DIR}"
    [ -n "$subdir" ] && base_output="${OUTPUT_DIR}/${subdir}"
    
    local output_prefix="${base_output}/${sample_id}/${sample_id}_"
    mkdir -p "${base_output}/${sample_id}"
    
    "$STAR_BIN" \
        --runThreadN "$THREADS" \
        --runMode alignReads \
        --genomeDir "$STAR_INDEX_DIR" \
        --readFilesIn "$r1_file" "$r2_file" \
        --sjdbGTFfile "$REFERENCE_GTF" \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "$output_prefix" \
        --outFilterMultimapNmax 20 \
        --alignIntronMax 1000000 \
        2>&1 | tee "${OUTPUT_DIR}/logs/${sample_id}_star.log"
    
    if [ -f "${output_prefix}Aligned.sortedByCoord.out.bam" ]; then
        samtools index "${output_prefix}Aligned.sortedByCoord.out.bam"
        log_success "STAR completed: ${sample_id}"
        return 0
    else
        log_error "STAR failed: ${sample_id}"
        return 1
    fi
}

align_hisat2_pe() {
    local sample_id=$1
    local r1_file=$2
    local r2_file=$3
    local subdir=$4
    
    log_info "Aligning: ${sample_id} (HISAT2)"
    
    local base_output="${OUTPUT_DIR}"
    [ -n "$subdir" ] && base_output="${OUTPUT_DIR}/${subdir}"
    
    mkdir -p "${base_output}/${sample_id}"
    local output_sam="${base_output}/${sample_id}/${sample_id}.sam"
    local output_bam="${base_output}/${sample_id}/${sample_id}.bam"
    local output_sorted="${base_output}/${sample_id}/${sample_id}Aligned.sortedByCoord.out.bam"
    
    "$HISAT2_BIN" \
        -q \
        -x "$HISAT2_INDEX_PREFIX" \
        -1 "$r1_file" \
        -2 "$r2_file" \
        -p "$THREADS" \
        -S "$output_sam" \
        2>&1 | tee "${OUTPUT_DIR}/logs/${sample_id}_hisat2.log"
    
    if [ $? -eq 0 ] && [ -f "$output_sam" ]; then
        samtools view -bS "$output_sam" > "$output_bam"
        samtools sort "$output_bam" -o "$output_sorted"
        samtools index "$output_sorted"
        rm -f "$output_sam" "$output_bam"
        
        log_success "HISAT2 completed: ${sample_id}"
        return 0
    else
        log_error "HISAT2 failed: ${sample_id}"
        return 1
    fi
}

process_all_samples() {
    if [ "$PAIRED_END" = false ]; then
        log_error "Single-end mode not yet implemented"
        exit 1
    fi
    
    log_info "Scanning for paired-end samples in: $INPUT_DIR"
    
    local r1_files=$(find "$INPUT_DIR" -type f \( -name "*_R1*.f*q*" -o -name "*_r1*.f*q*" -o -name "*_1.f*q*" \) | sort)
    
    if [ -z "$r1_files" ]; then
        log_error "No paired-end FASTQ files found in $INPUT_DIR"
        exit 1
    fi
    
    local total_count=0
    local processed_count=0
    
    for r1_file in $r1_files; do
        total_count=$((total_count + 1))
        
        local basename=$(basename "$r1_file")
        local sample_id=$(echo "$basename" | sed -E 's/_(R1|r1|1).*//g')
        
        local r2_file=$(find_r2_file "$r1_file")
        
        if [ ! -f "$r2_file" ]; then
            log_warning "R2 not found for $sample_id, skipping..."
            continue
        fi
        
        processed_count=$((processed_count + 1))
        
        echo ""
        log_info "Sample ${processed_count}: ${sample_id}"
        
        case "$ALIGNER" in
            star)
                align_star_pe "$sample_id" "$r1_file" "$r2_file" ""
                ;;
            hisat2)
                align_hisat2_pe "$sample_id" "$r1_file" "$r2_file" ""
                ;;
            both)
                align_star_pe "$sample_id" "$r1_file" "$r2_file" "STAR"
                align_hisat2_pe "$sample_id" "$r1_file" "$r2_file" "HISAT2"
                ;;
        esac
    done
    
    if [ $processed_count -eq 0 ]; then
        log_error "No valid sample pairs processed"
        exit 1
    fi
    
    log_success "Processed ${processed_count} of ${total_count} samples"
}

generate_summary() {
    local summary_file="${OUTPUT_DIR}/alignment_summary.txt"
    
    cat > "$summary_file" << EOF
========================================
Alignment Module Summary
========================================
Date: $(date)
Input: ${INPUT_DIR}
Output: ${OUTPUT_DIR}
Aligner: ${ALIGNER^^}
Reference: ${REFERENCE_GENOME}
GTF: ${REFERENCE_GTF}
Threads: ${THREADS}
========================================
EOF
    
    log_success "Summary: $summary_file"
}

display_final_summary() {
    echo ""
    echo "========================================"
    echo "  Alignment Complete"
    echo "========================================"
    echo ""
    log_success "All alignments completed!"
    echo ""
    log_info "Results: ${OUTPUT_DIR}"
    
    if [ "$ALIGNER" = "both" ]; then
        echo "  STAR:   ${OUTPUT_DIR}/STAR/"
        echo "  HISAT2: ${OUTPUT_DIR}/HISAT2/"
    fi
    
    echo "  Logs:   ${OUTPUT_DIR}/logs/"
    echo ""
    echo "========================================"
}

parse_arguments() {
    if [ $# -eq 0 ]; then
        usage
    fi
    
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
    
    if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
        log_error "Input and output directories required"
        usage
    fi
    
    if [ ! -d "$INPUT_DIR" ]; then
        log_error "Input directory not found: $INPUT_DIR"
        exit 1
    fi
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
    
    local start_time=$(date +%s)
    log_info "Started: $(date)"
    
    if [ -n "$SAMPLE_FILE" ]; then
        log_error "Sample file mode not yet implemented"
        exit 1
    else
        process_all_samples
    fi
    
    generate_summary
    
    local end_time=$(date +%s)
    local runtime=$((end_time - start_time))
    log_info "Runtime: ${runtime}s"
    
    display_final_summary
}

main "$@"
