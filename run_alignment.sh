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

# Python/Conda settings
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
    # Interactive aligner selection
    $0 -i qc_results/trimmed/PE -o alignment_results

    # Use STAR only
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner star -t 20

    # Use HISAT2 only
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner hisat2 -t 12

    # Use both aligners for comparison
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner both

EOF
    exit 1
}

command_exists() { command -v "$1" >/dev/null 2>&1; }

# Setup Python for HISAT2 (system Python OR conda)
setup_python_for_hisat2() {
    if [ "$USE_SYSTEM_PYTHON" = "yes" ]; then
        log_info "Using system Python: $PYTHON_BIN"
        return 0
    fi
    
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
    
    log_error "Invalid Python configuration"
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
        log_info "Both STAR and HISAT2 are installed"
        echo ""
        echo "  1) STAR (high memory ~32GB, best for splice detection)"
        echo "  2) HISAT2 (low memory ~8GB, faster)"
        echo "  3) BOTH (run with both aligners for comparison)"
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
        log_info "Please run: ./install.sh"
        exit 1
    fi
    
    log_success "Selected aligner: ${ALIGNER^^}"
}

# Verify aligner availability
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
            setup_python_for_hisat2
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
            setup_python_for_hisat2
            log_success "Both STAR and HISAT2 are available"
            ;;
        *)
            log_error "Invalid aligner: $ALIGNER"
            exit 1
            ;;
    esac
}

# Check reference genome
check_reference() {
    log_info "Checking reference genome..."
    
    [ -z "$REFERENCE_GENOME" ] && [ -n "$REFERENCE_FASTA" ] && [ -f "$REFERENCE_FASTA" ] && REFERENCE_GENOME="$REFERENCE_FASTA"
    [ -z "$REFERENCE_GTF" ] && [ -n "${REFERENCE_DIR}" ] && [ -f "${REFERENCE_DIR}/Homo_sapiens.GRCh38.113.gtf" ] && REFERENCE_GTF="${REFERENCE_DIR}/Homo_sapiens.GRCh38.113.gtf"
    
    if [ ! -f "$REFERENCE_GENOME" ]; then
        log_error "Reference genome not found: $REFERENCE_GENOME"
        log_info "Download with: ./install.sh --download-reference"
        exit 1
    fi
    
    if [ ! -f "$REFERENCE_GTF" ]; then
        log_error "Reference GTF not found: $REFERENCE_GTF"
        log_info "Download with: ./install.sh --download-reference"
        exit 1
    fi
    
    log_success "Reference genome: $REFERENCE_GENOME"
    log_success "Reference GTF: $REFERENCE_GTF"
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
        log_success "STAR index found: $STAR_INDEX_DIR"
    fi
}

check_or_build_hisat2_index() {
    log_info "Checking HISAT2 genome index..."
    
    if [ ! -f "${HISAT2_INDEX_PREFIX}.1.ht2" ]; then
        log_warning "HISAT2 index not found"
        build_hisat2_index
    else
        log_success "HISAT2 index found: $HISAT2_INDEX_PREFIX"
    fi
}

build_star_index() {
    log_info "Building STAR genome index..."
    log_warning "This may take 30-60 minutes and requires ~32 GB RAM"
    
    echo -ne "Build STAR index now? [Y/n] (auto-continue in ${PROMPT_TIMEOUT}s): "
    read -t $PROMPT_TIMEOUT -n 1 -r REPLY || REPLY="y"
    echo
    
    if [[ $REPLY =~ ^[Nn]$ ]]; then
        log_error "STAR index required for alignment"
        exit 1
    fi
    
    mkdir -p "$STAR_INDEX_DIR"
    log_info "Running STAR genome index build..."
    
    "$STAR_BIN" \
        --runThreadN "$THREADS" \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX_DIR" \
        --genomeFastaFiles "$REFERENCE_GENOME" \
        --sjdbGTFfile "$REFERENCE_GTF" \
        --sjdbOverhang 100 \
        2>&1 | tee "${SCRIPT_DIR}/star_index_build.log"
    
    if [ $? -eq 0 ]; then
        log_success "STAR index built successfully"
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
        log_error "HISAT2 index required for alignment"
        exit 1
    fi
    
    mkdir -p "$HISAT2_INDEX_DIR"
    
    # Use system Python OR conda environment
    if [ "$USE_SYSTEM_PYTHON" = "yes" ]; then
        log_info "Building HISAT2 index with system Python..."
        "$HISAT2_BUILD_BIN" -f "$REFERENCE_GENOME" -p "$THREADS" "$HISAT2_INDEX_PREFIX" \
            2>&1 | tee "${SCRIPT_DIR}/hisat2_index_build.log"
        local STATUS=$?
    else
        log_info "Building HISAT2 index with conda environment '$CONDA_ENV_NAME'..."
        set +u
        source "${CONDA_DIR}/etc/profile.d/conda.sh"
        
        if ! conda activate "$CONDA_ENV_NAME" 2>/dev/null; then
            log_error "Failed to activate conda environment"
            exit 1
        fi
        
        "$HISAT2_BUILD_BIN" -f "$REFERENCE_GENOME" -p "$THREADS" "$HISAT2_INDEX_PREFIX" \
            2>&1 | tee "${SCRIPT_DIR}/hisat2_index_build.log"
        local STATUS=$?
        
        conda deactivate
        set -u
    fi
    
    if [ $STATUS -eq 0 ]; then
        log_success "HISAT2 index built successfully"
    else
        log_error "HISAT2 index build failed"
        exit 1
    fi
}

create_output_dirs() {
    log_info "Creating output directory structure..."
    mkdir -p "${OUTPUT_DIR}/logs" "${OUTPUT_DIR}/temp"
    
    if [ "$ALIGNER" = "both" ]; then
        mkdir -p "${OUTPUT_DIR}/STAR" "${OUTPUT_DIR}/HISAT2"
        log_info "Created subdirectories for STAR and HISAT2"
    fi
    
    log_success "Output directories created"
}

detect_pattern() {
    local filename=$1
    if [[ "$filename" =~ _R1|_r1 ]]; then
        echo "_R1/_R2"
    elif [[ "$filename" =~ _1[\._] ]]; then
        echo "_1/_2"
    else
        echo "unknown"
    fi
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
        2>&1 | tee "${OUTPUT_DIR}/logs/${sample_id}_star.log"
    
    if [ -f "${output_prefix}Aligned.sortedByCoord.out.bam" ]; then
        log_info "Indexing BAM file..."
        samtools index "${output_prefix}Aligned.sortedByCoord.out.bam"
        log_success "STAR alignment completed: ${sample_id}"
        return 0
    else
        log_error "STAR alignment failed: ${sample_id}"
        return 1
    fi
}

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
    local output_sorted="${base_output}/${sample_id}/${sample_id}Aligned.sortedByCoord.out.bam"
    
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
        return 0
    else
        log_error "HISAT2 alignment failed: ${sample_id}"
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
        log_info "Expected patterns: *_R1*.fastq, *_r1*.fastq, *_1.fastq (and .fq, .gz variants)"
        exit 1
    fi
    
    local total_count=0
    local processed_count=0
    
    for r1_file in $r1_files; do
        total_count=$((total_count + 1))
        
        local basename=$(basename "$r1_file")
        local sample_id=$(echo "$basename" | sed -E 's/_(R1|r1|1).*//g')
        
        local r2_file=$(find_r2_file "$r1_file")
        local pattern=$(detect_pattern "$r1_file")
        
        if [ ! -f "$r2_file" ]; then
            log_warning "R2 not found for $sample_id (pattern: ${pattern}), skipping..."
            continue
        fi
        
        processed_count=$((processed_count + 1))
        
        echo ""
        log_info "========================================="
        log_info "Sample ${processed_count}: ${sample_id}"
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
    
    if [ $processed_count -eq 0 ]; then
        log_error "No valid sample pairs found and processed"
        exit 1
    fi
    
    log_success "Processed ${processed_count} of ${total_count} samples"
}

generate_summary() {
    log_info "Generating summary..."
    
    local summary_file="${OUTPUT_DIR}/alignment_summary.txt"
    
    cat > "$summary_file" << EOF
========================================
Alignment Module Summary
========================================
Date: $(date)
Input Directory: ${INPUT_DIR}
Output Directory: ${OUTPUT_DIR}

Aligner: ${ALIGNER^^}
Reference: ${REFERENCE_GENOME}
GTF: ${REFERENCE_GTF}
Threads: ${THREADS}

Output Structure:
$(if [ "$ALIGNER" = "both" ]; then
    echo "  STAR results:   ${OUTPUT_DIR}/STAR/[sample]/"
    echo "  HISAT2 results: ${OUTPUT_DIR}/HISAT2/[sample]/"
else
    echo "  Aligned BAM files: ${OUTPUT_DIR}/[sample]/"
fi)
  Logs: ${OUTPUT_DIR}/logs/

Next Steps:
  - Check BAM files in sample directories
  - Review alignment logs
  - Proceed with quantification (featureCounts, RSEM, etc.)

========================================
EOF
    
    log_success "Summary saved: $summary_file"
}

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

parse_arguments() {
    if [ $# -eq 0 ]; then
        usage
    fi
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input)
                INPUT_DIR="$2"
                shift 2 ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2 ;;
            -s|--samples)
                SAMPLE_FILE="$2"
                shift 2 ;;
            -t|--threads)
                THREADS="$2"
                shift 2 ;;
            --aligner)
                ALIGNER="$2"
                shift 2 ;;
            --reference)
                REFERENCE_GENOME="$2"
                shift 2 ;;
            --gtf)
                REFERENCE_GTF="$2"
                shift 2 ;;
            --single-end)
                PAIRED_END=false
                shift ;;
            -h|--help)
                usage ;;
            *)
                log_error "Unknown option: $1"
                usage ;;
        esac
    done
    
    if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
        log_error "Input and output directories are required"
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
    log_info "Aligner: ${ALIGNER^^}"
    log_info "Threads: $THREADS"
    
    echo ""
    if [ -n "$SAMPLE_FILE" ]; then
        log_error "Sample file mode not yet implemented"
        exit 1
    else
        process_all_samples
    fi
    
    echo ""
    generate_summary
    
    local end_time=$(date +%s)
    local runtime=$((end_time - start_time))
    local hours=$((runtime / 3600))
    local minutes=$(((runtime % 3600) / 60))
    local seconds=$((runtime % 60))
    
    log_info "Completed: $(date)"
    log_info "Runtime: ${hours}h ${minutes}m ${seconds}s"
    
    display_final_summary
}

main "$@"
