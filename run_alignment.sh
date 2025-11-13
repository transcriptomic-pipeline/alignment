#!/bin/bash

# Alignment Module - Main Execution Script
# Supports STAR and HISAT2 with auto-installation, auto-indexing, and "both" mode
# Uses system Python when available, conda only when necessary

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

# Timeout for auto-continue prompts (in seconds)
PROMPT_TIMEOUT=30

usage() {
    cat << EOF
Usage: $0 -i <input_dir> -o <output_dir> [options]

Required:
    -i, --input         Input directory with trimmed FASTQ files (e.g., qc_results/trimmed/PE)
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
    # Interactive mode (will prompt for aligner)
    $0 -i qc_results/trimmed/PE -o alignment_results

    # Use STAR with 20 threads
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner star -t 20

    # Use HISAT2
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner hisat2

    # Use both aligners for comparison
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner both

EOF
    exit 1
}

command_exists() { command -v "$1" >/dev/null 2>&1; }

# Check installation and load configuration
check_installation() {
    log_info "Checking for installed tools..."
    
    if [ ! -f "$CONFIG_FILE" ]; then
        log_error "Configuration not found"
        log_info "Please run: ./install.sh"
        exit 1
    fi
    
    # Load configuration
    source "$CONFIG_FILE"
    
    if [ -f "$REF_CONFIG" ]; then
        source "$REF_CONFIG"
    fi
    
    log_success "Configuration loaded"
}

# Check if conda environment exists and is properly set up (ONLY IF NEEDED)
check_conda_env() {
    # If using system Python, skip conda check entirely
    if [ "$USE_SYSTEM_PYTHON" = "yes" ]; then
        log_info "Using system Python for HISAT2 (no conda needed)"
        return 0
    fi
    
    log_info "Checking conda environment for HISAT2..."
    if [ -z "$CONDA_DIR" ] || [ ! -d "$CONDA_DIR" ]; then
        log_error "Conda directory not found: $CONDA_DIR"
        log_error "Please run: ./install.sh --aligner hisat2"
        exit 1
    fi
    
    if [ ! -f "$CONDA_DIR/bin/conda" ]; then
        log_error "Conda not properly installed at: $CONDA_DIR"
        log_error "Please run: ./install.sh --aligner hisat2"
        exit 1
    fi
    
    # Initialize conda in this shell session
    if [ ! -f "${CONDA_DIR}/etc/profile.d/conda.sh" ]; then
        log_error "Conda profile script not found at: ${CONDA_DIR}/etc/profile.d/conda.sh"
        log_error "Please reinstall: ./install.sh --aligner hisat2"
        exit 1
    fi
    
    # Source conda profile
    source "${CONDA_DIR}/etc/profile.d/conda.sh"
    
    # Check if environment exists
    if ! conda env list | grep -q "^${CONDA_ENV_NAME} "; then
        log_error "Conda environment '$CONDA_ENV_NAME' not found"
        log_error "Available environments:"
        conda env list
        log_error "Please run: ./install.sh --aligner hisat2"
        exit 1
    fi
    
    log_success "Conda environment '$CONDA_ENV_NAME' is available"
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
        echo "Which aligner do you want to use?"
        echo "  1) STAR (high memory ~32GB, best for splice detection)"
        echo "  2) HISAT2 (low memory ~8GB, faster)"
        echo "  3) BOTH (run alignment with both aligners for comparison)"
        echo ""
        read -p "Enter choice [1-3]: " ALIGNER_CHOICE
        
        case "$ALIGNER_CHOICE" in
            1) ALIGNER="star" ;;
            2) ALIGNER="hisat2" ;;
            3) ALIGNER="both" ;;
            *)
                log_error "Invalid choice. Please select 1, 2, or 3."
                exit 1
                ;;
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
                log_info "Please run: ./install.sh --aligner star"
                exit 1
            fi
            log_success "STAR is available"
            ;;
        hisat2)
            if [ "$HISAT2_INSTALLED" != "yes" ] || [ ! -f "$HISAT2_BIN" ]; then
                log_error "HISAT2 not installed"
                log_info "Please run: ./install.sh --aligner hisat2"
                exit 1
            fi
            # Check for conda/python environment
            check_conda_env
            log_success "HISAT2 is available"
            ;;
        both)
            if [ "$STAR_INSTALLED" != "yes" ] || [ ! -f "$STAR_BIN" ]; then
                log_error "STAR not installed"
                log_info "Please run: ./install.sh --aligner both"
                exit 1
            fi
            if [ "$HISAT2_INSTALLED" != "yes" ] || [ ! -f "$HISAT2_BIN" ]; then
                log_error "HISAT2 not installed"
                log_info "Please run: ./install.sh --aligner both"
                exit 1
            fi
            # Check for conda/python environment
            check_conda_env
            log_success "Both STAR and HISAT2 are available"
            ;;
        *)
            log_error "Invalid aligner: $ALIGNER"
            exit 1
            ;;
    esac
}

# Check reference genome files
check_reference() {
    log_info "Checking reference genome..."
    
    # Use configured reference if not provided via command line
    if [ -z "$REFERENCE_GENOME" ]; then
        if [ -n "$REFERENCE_FASTA" ] && [ -f "$REFERENCE_FASTA" ]; then
            REFERENCE_GENOME="$REFERENCE_FASTA"
        fi
    fi
    
    if [ -z "$REFERENCE_GTF" ]; then
        if [ -n "$REFERENCE_DIR" ] && [ -f "${REFERENCE_DIR}/Homo_sapiens.GRCh38.113.gtf" ]; then
            REFERENCE_GTF="${REFERENCE_DIR}/Homo_sapiens.GRCh38.113.gtf"
        fi
    fi
    
    if [ ! -f "$REFERENCE_GENOME" ]; then
        log_error "Reference genome not found: $REFERENCE_GENOME"
        log_info "Please download with: ./install.sh --download-reference"
        exit 1
    fi
    
    if [ ! -f "$REFERENCE_GTF" ]; then
        log_error "Reference GTF not found: $REFERENCE_GTF"
        log_info "Please download with: ./install.sh --download-reference"
        exit 1
    fi
    
    log_success "Reference genome: $REFERENCE_GENOME"
    log_success "Reference GTF: $REFERENCE_GTF"
}

# Check or build genome indexes
check_or_build_index() {
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
}

# Check STAR index
check_or_build_star_index() {
    log_info "Checking STAR genome index..."
    
    # Check if index directory exists and contains required files
    if [ ! -d "$STAR_INDEX_DIR" ] || [ -z "$(ls -A $STAR_INDEX_DIR 2>/dev/null)" ] || [ ! -f "${STAR_INDEX_DIR}/genomeParameters.txt" ]; then
        log_warning "STAR index not found or incomplete"
        build_star_index
    else
        log_success "STAR index found: $STAR_INDEX_DIR"
    fi
}

# Check HISAT2 index
check_or_build_hisat2_index() {
    log_info "Checking HISAT2 genome index..."
    
    if [ ! -f "${HISAT2_INDEX_PREFIX}.1.ht2" ]; then
        log_warning "HISAT2 index not found"
        build_hisat2_index
    else
        log_success "HISAT2 index found: $HISAT2_INDEX_PREFIX"
    fi
}

# Build STAR index with AUTO-CONTINUE
build_star_index() {
    log_info "Building STAR genome index..."
    log_warning "This may take 30-60 minutes and requires ~32 GB RAM"
    
    # AUTO-CONTINUE: Wait ${PROMPT_TIMEOUT} seconds, then proceed with 'y'
    echo -ne "Build STAR index now? [Y/n] (auto-continue in ${PROMPT_TIMEOUT}s): "
    read -t $PROMPT_TIMEOUT -n 1 -r REPLY || true
    echo
    
    # If no response (timeout), default to 'y'
    if [ -z "$REPLY" ]; then
        REPLY="y"
        log_warning "No response in ${PROMPT_TIMEOUT}s, automatically proceeding with 'y'"
    fi
    
    if [[ $REPLY =~ ^[Nn]$ ]]; then
        log_error "STAR index required for alignment"
        exit 1
    fi
    
    mkdir -p "$STAR_INDEX_DIR"
    log_info "Running STAR genome index build..."
    log_info "Index directory: $STAR_INDEX_DIR"
    
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
        log_info "Index location: $STAR_INDEX_DIR"
    else
        log_error "STAR index build failed"
        log_info "Check log: ${SCRIPT_DIR}/star_index_build.log"
        exit 1
    fi
}

# Build HISAT2 index with AUTO-CONTINUE and proper conda activation (ONLY IF NEEDED)
build_hisat2_index() {
    log_info "Building HISAT2 genome index..."
    log_warning "This may take 30-45 minutes and requires ~8 GB RAM"
    
    # AUTO-CONTINUE: Wait ${PROMPT_TIMEOUT} seconds, then proceed with 'y'
    echo -ne "Build HISAT2 index now? [Y/n] (auto-continue in ${PROMPT_TIMEOUT}s): "
    read -t $PROMPT_TIMEOUT -n 1 -r REPLY || true
    echo
    
    if [ -z "$REPLY" ]; then
        REPLY="y"
        log_warning "No response in ${PROMPT_TIMEOUT}s, automatically proceeding with 'y'"
    fi
    
    if [[ $REPLY =~ ^[Nn]$ ]]; then
        log_error "HISAT2 index required for alignment"
        exit 1
    fi
    
    mkdir -p "$HISAT2_INDEX_DIR"
    
    # Check if we need conda or can use system Python
    if [ "$USE_SYSTEM_PYTHON" = "yes" ]; then
        log_info "Running HISAT2 genome index build with system Python..."
        "$HISAT2_BUILD_BIN" \
            -f "$REFERENCE_GENOME" \
            -p "$THREADS" \
            "$HISAT2_INDEX_PREFIX" \
            2>&1 | tee "${SCRIPT_DIR}/hisat2_index_build.log"
        local BUILD_STATUS=$?
    else
        log_info "Activating conda environment for HISAT2 index build..."
        set +u # Temporarily disable unset variable check for conda
        
        if [ ! -f "${CONDA_DIR}/etc/profile.d/conda.sh" ]; then
            log_error "Conda profile script not found"
            exit 1
        fi
        
        source "${CONDA_DIR}/etc/profile.d/conda.sh"
        
        if ! conda activate "$CONDA_ENV_NAME" 2>/dev/null; then
            log_error "Failed to activate conda environment: $CONDA_ENV_NAME"
            log_error "Please run: ./install.sh --aligner hisat2"
            exit 1
        fi
        
        log_info "Running HISAT2 genome index build in conda environment '$CONDA_ENV_NAME'..."
        "$HISAT2_BUILD_BIN" \
            -f "$REFERENCE_GENOME" \
            -p "$THREADS" \
            "$HISAT2_INDEX_PREFIX" \
            2>&1 | tee "${SCRIPT_DIR}/hisat2_index_build.log"
        
        local BUILD_STATUS=$?
        conda deactivate
        set -u # Re-enable unset variable check
    fi
    
    if [ $BUILD_STATUS -eq 0 ]; then
        log_success "HISAT2 index built successfully"
        log_info "Index location: $HISAT2_INDEX_DIR"
    else
        log_error "HISAT2 index build failed"
        log_info "Check log: ${SCRIPT_DIR}/hisat2_index_build.log"
        exit 1
    fi
}

# Create output directory structure
create_output_dirs() {
    log_info "Creating output directory structure..."
    
    mkdir -p "${OUTPUT_DIR}/logs"
    mkdir -p "${OUTPUT_DIR}/temp"
    
    if [ "$ALIGNER" = "both" ]; then
        mkdir -p "${OUTPUT_DIR}/STAR"
        mkdir -p "${OUTPUT_DIR}/HISAT2"
        log_info "Created subdirectories for STAR and HISAT2"
    fi
    
    log_success "Output directories created"
}

# Detect Read 2 file from Read 1 filename
find_r2_file() {
    local r1_file=$1
    local r2_file=""
    
    # Try different naming patterns
    if [[ "$r1_file" =~ _R1 ]]; then
        r2_file=$(echo "$r1_file" | sed 's/_R1/_R2/g')
    elif [[ "$r1_file" =~ _r1 ]]; then
        r2_file=$(echo "$r1_file" | sed 's/_r1/_r2/g')
    elif [[ "$r1_file" =~ _1\. ]]; then
        r2_file=$(echo "$r1_file" | sed 's/_1\./_2./g')
    elif [[ "$r1_file" =~ _1_ ]]; then
        r2_file=$(echo "$r1_file" | sed 's/_1_/_2_/g')
    fi
    
    echo "$r2_file"
}

# STAR alignment (paired-end)
align_star_pe() {
    local sample_id=$1
    local r1_file=$2
    local r2_file=$3
    local subdir=$4
    
    log_info "Aligning sample: ${sample_id} (STAR)"
    
    local base_output="${OUTPUT_DIR}"
    if [ -n "$subdir" ]; then
        base_output="${OUTPUT_DIR}/${subdir}"
    fi
    
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
        --alignMatesGapMax 1000000 \
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

# HISAT2 alignment (paired-end)
align_hisat2_pe() {
    local sample_id=$1
    local r1_file=$2
    local r2_file=$3
    local subdir=$4
    
    log_info "Aligning sample: ${sample_id} (HISAT2)"
    
    local base_output="${OUTPUT_DIR}"
    if [ -n "$subdir" ]; then
        base_output="${OUTPUT_DIR}/${subdir}"
    fi
    
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

# Process all samples
process_all_samples() {
    if [ "$PAIRED_END" = false ]; then
        log_error "Single-end mode not yet implemented"
        log_info "Currently only paired-end reads are supported"
        exit 1
    fi
    
    log_info "Scanning for paired-end samples in: $INPUT_DIR"
    log_info "Processing all samples in: $INPUT_DIR"
    
    # Find all R1 files
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
        
        if [ ! -f "$r2_file" ]; then
            log_warning "R2 not found for $sample_id, skipping..."
            log_warning "Expected: $r2_file"
            continue
        fi
        
        processed_count=$((processed_count + 1))
        
        echo ""
        log_info "========================================="
        log_info "Sample ${processed_count}: ${sample_id}"
        log_info "========================================="
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
    
    echo ""
    log_success "Processed ${processed_count} of ${total_count} samples"
}

# Generate summary
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
  - Review alignment logs in logs/ directory
  - Proceed with quantification (featureCounts, RSEM, etc.)

========================================
EOF
    
    log_success "Summary saved: $summary_file"
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

# Parse command-line arguments
parse_arguments() {
    if [ $# -eq 0 ]; then
        usage
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
            -s|--samples)
                SAMPLE_FILE="$2"
                shift 2
                ;;
            -t|--threads)
                THREADS="$2"
                shift 2
                ;;
            --aligner)
                ALIGNER="$2"
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
    
    if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
        log_error "Input and output directories are required"
        usage
    fi
    
    if [ ! -d "$INPUT_DIR" ]; then
        log_error "Input directory not found: $INPUT_DIR"
        exit 1
    fi
}

# Main function
main() {
    echo "========================================"
    echo "  Alignment Module"
    echo "========================================"
    echo ""
    
    parse_arguments "$@"
    check_installation
    
    # Prompt for aligner if not specified
    if [ -z "$ALIGNER" ]; then
        prompt_aligner_selection
    fi
    
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
        log_info "Currently processes all samples in input directory"
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

# Run main function
main "$@"
