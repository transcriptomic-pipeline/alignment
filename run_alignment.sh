#!/bin/bash

# Alignment Module - Main Execution Script
# Supports STAR and HISAT2 with auto-installation, auto-indexing, and "both" mode
# Uses pre-installed Miniconda environment for HISAT2

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

# Conda settings (loaded from config)
CONDA_DIR=""
CONDA_ENV_NAME=""

# Timeout for auto-continue prompts (in seconds)
PROMPT_TIMEOUT=30

# Usage
usage() {
    cat << EOF
Usage: $0 -i <input_dir> -o <output_dir> [--aligner <star|hisat2|both>] [options]

Required:
    -i, --input         Input directory with trimmed FASTQ files
                        (e.g., qc_results/trimmed/PE/)
    -o, --output        Output directory for aligned BAM files

Optional:
    --aligner           Aligner to use: star, hisat2, both, or omit to prompt at runtime
    -s, --samples       Sample list file (one sample ID per line)
    -t, --threads       Number of threads (default: 12)
    --reference         Reference FASTA file (auto-detected if downloaded)
    --gtf               Reference GTF file (auto-detected if downloaded)
    --single-end        Process as single-end reads (default: paired-end)
    -h, --help          Show this help message

Aligner Selection:
    --aligner star      Use STAR only (high memory ~32GB, best for splice detection)
    --aligner hisat2    Use HISAT2 only (low memory ~8GB, faster)
    --aligner both      Use BOTH aligners (creates separate outputs for comparison)
    (no --aligner)      Prompt at runtime to choose interactively

Examples:
    # Ask at runtime which aligner to use
    $0 -i qc_results/trimmed/PE -o alignment_results -t 20

    # Use STAR only
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner star -t 20

    # Use HISAT2 only
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner hisat2 -t 12

    # Use BOTH aligners (for comparison)
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner both -t 20

Integration with QC Module:
    # Direct from QC output
    $0 -i qc_results/trimmed/PE -o alignment_results --aligner star

EOF
    exit 1
}

command_exists() { command -v "$1" >/dev/null 2>&1; }

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

# Check installation
check_installation() {
    log_info "Checking for installed tools..."
    
    if [ ! -f "$CONFIG_FILE" ]; then
        log_warning "Alignment tools not installed"
        log_info "Installation required"
        
        # Ask which aligner to install
        if [ "$ALIGNER" = "" ]; then
            echo ""
            echo "Which aligner would you like to install?"
            echo "  1) STAR (recommended for high-memory systems)"
            echo "  2) HISAT2 (recommended for low-memory systems)"
            echo "  3) Both (choose at runtime)"
            echo ""
            read -p "Enter choice [1-3]: " INSTALL_CHOICE
            
            case "$INSTALL_CHOICE" in
                1) INSTALL_ALIGNER="star" ;;
                2) INSTALL_ALIGNER="hisat2" ;;
                3) INSTALL_ALIGNER="both" ;;
                *) log_error "Invalid choice"; exit 1 ;;
            esac
        elif [ "$ALIGNER" = "both" ]; then
            INSTALL_ALIGNER="both"
        else
            INSTALL_ALIGNER="$ALIGNER"
        fi
        
        read -p "Run installation now? [Y/n] " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Nn]$ ]]; then
            if [ -f "${SCRIPT_DIR}/install.sh" ]; then
                log_info "Running installer..."
                bash "${SCRIPT_DIR}/install.sh" --aligner "$INSTALL_ALIGNER" --download-reference
                
                if [ $? -eq 0 ] && [ -f "$CONFIG_FILE" ]; then
                    log_success "Installation completed"
                else
                    log_error "Installation failed"
                    exit 1
                fi
            else
                log_error "Installer not found: ${SCRIPT_DIR}/install.sh"
                exit 1
            fi
        else
            log_error "Tools must be installed to proceed"
            exit 1
        fi
    fi
    
    # Load configuration
    source "$CONFIG_FILE"
    [ -f "$REF_CONFIG" ] && source "$REF_CONFIG"
    
    log_success "Configuration loaded"
}

# Prompt for aligner selection at runtime
prompt_aligner_selection() {
    echo ""
    echo "========================================"
    echo "  Aligner Selection"
    echo "========================================"
    echo ""
    
    # Check which aligners are installed
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

# Verify aligner is available
verify_aligner() {
    local need_star=false
    local need_hisat2=false
    
    case "$ALIGNER" in
        star)
            need_star=true
            ;;
        hisat2)
            need_hisat2=true
            ;;
        both)
            need_star=true
            need_hisat2=true
            ;;
        *)
            log_error "Invalid aligner: $ALIGNER"
            exit 1
            ;;
    esac
    
    # Check STAR if needed
    if [ "$need_star" = true ]; then
        if [ "$STAR_INSTALLED" != "yes" ] || [ ! -f "$STAR_BIN" ]; then
            log_error "STAR not installed"
            log_info "Run: ./install.sh --aligner star"
            exit 1
        fi
        log_success "STAR is available"
    fi
    
    # Check HISAT2 if needed
    if [ "$need_hisat2" = true ]; then
        if [ "$HISAT2_INSTALLED" != "yes" ] || [ ! -f "$HISAT2_BIN" ]; then
            log_error "HISAT2 not installed"
            log_info "Run: ./install.sh --aligner hisat2"
            exit 1
        fi
        log_success "HISAT2 is available"
        
        # Check conda environment for HISAT2
        check_conda_env
    fi
}

# Check reference genome
check_reference() {
    log_info "Checking reference genome..."
    
    # Use provided reference or auto-detect
    if [ -z "$REFERENCE_GENOME" ] && [ -n "$REFERENCE_FASTA" ] && [ -f "$REFERENCE_FASTA" ]; then
        REFERENCE_GENOME="$REFERENCE_FASTA"
    fi
    
    if [ -z "$REFERENCE_GTF" ] && [ -f "${REFERENCE_DIR}/Homo_sapiens.GRCh38.113.gtf" ]; then
        REFERENCE_GTF="${REFERENCE_DIR}/Homo_sapiens.GRCh38.113.gtf"
    fi
    
    if [ ! -f "$REFERENCE_GENOME" ]; then
        log_error "Reference genome not found"
        log_info "Download with: ./install.sh --download-reference"
        exit 1
    fi
    
    if [ ! -f "$REFERENCE_GTF" ]; then
        log_error "Reference GTF not found"
        log_info "Download with: ./install.sh --download-reference"
        exit 1
    fi
    
    log_success "Reference genome: $REFERENCE_GENOME"
    log_success "Reference GTF: $REFERENCE_GTF"
}

# Check or build genome index
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

# Check or build STAR index
check_or_build_star_index() {
    log_info "Checking STAR genome index..."
    
    # Check if index directory exists AND contains the required genomeParameters.txt file
    if [ ! -d "$STAR_INDEX_DIR" ] || [ -z "$(ls -A $STAR_INDEX_DIR 2>/dev/null)" ] || [ ! -f "${STAR_INDEX_DIR}/genomeParameters.txt" ]; then
        log_warning "STAR index not found or incomplete"
        build_star_index
    else
        log_success "STAR index found: $STAR_INDEX_DIR"
    fi
    
    GENOME_INDEX_DIR="$STAR_INDEX_DIR"
}

# Check or build HISAT2 index
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
        log_error "STAR index build failed (check RAM - requires ~32GB)"
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

# Create output directories
create_output_dirs() {
    log_info "Creating output directory structure..."
    
    mkdir -p "${OUTPUT_DIR}/logs"
    mkdir -p "${OUTPUT_DIR}/temp"
    
    # Create aligner-specific directories if using both
    if [ "$ALIGNER" = "both" ]; then
        mkdir -p "${OUTPUT_DIR}/STAR"
        mkdir -p "${OUTPUT_DIR}/HISAT2"
        log_info "Created subdirectories for STAR and HISAT2"
    fi
    
    log_success "Output directories created"
}

# Detect file naming pattern
detect_pattern() {
    local filename=$1
    
    if [[ "$filename" =~ _R1 ]] || [[ "$filename" =~ _r1 ]]; then
        echo "_R1/_R2"
    elif [[ "$filename" =~ _1\. ]] || [[ "$filename" =~ _1_ ]]; then
        echo "_1/_2"
    elif [[ "$filename" =~ _paired ]]; then
        echo "_paired"
    else
        echo "unknown"
    fi
}

# Find R2 file based on R1
find_r2_file() {
    local r1_file=$1
    local r2_file=""
    
    # Try different patterns
    if [[ "$r1_file" =~ _R1 ]]; then
        r2_file=$(echo "$r1_file" | sed 's/_R1/_R2/g')
    elif [[ "$r1_file" =~ _r1 ]]; then
        r2_file=$(echo "$r1_file" | sed 's/_r1/_r2/g')
    elif [[ "$r1_file" =~ _1\. ]]; then
        r2_file=$(echo "$r1_file" | sed 's/_1\./_2./g')
    elif [[ "$r1_file" =~ _1_ ]]; then
        r2_file=$(echo "$r1_file" | sed 's/_1_/_2_/g')
    elif [[ "$r1_file" =~ _paired ]]; then
        r2_file=$(echo "$r1_file" | sed 's/_R1_paired/_R2_paired/g')
    fi
    
    echo "$r2_file"
}

# STAR alignment for paired-end
align_star_pe() {
    local sample_id=$1
    local r1_file=$2
    local r2_file=$3
    local subdir=$4  # Optional subdirectory (e.g., "STAR" when using both)
    
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

# HISAT2 alignment for paired-end
align_hisat2_pe() {
    local sample_id=$1
    local r1_file=$2
    local r2_file=$3
    local subdir=$4  # Optional subdirectory
    
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
    else
        log_error "HISAT2 alignment failed: ${sample_id}"
        return 1
    fi
}

# Process all samples with intelligent pattern detection
process_all_samples() {
    if [ "$PAIRED_END" = false ]; then
        log_error "Single-end mode not yet implemented"
        log_info "Currently only paired-end reads are supported"
        exit 1
    fi
    
    log_info "Scanning for paired-end samples in: $INPUT_DIR"
    echo ""
    
    # Find ALL fastq files first to analyze naming patterns
    local all_fastq_files=$(find "$INPUT_DIR" -type f \( -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)
    
    if [ -z "$all_fastq_files" ]; then
        log_error "No FASTQ files found in $INPUT_DIR"
        exit 1
    fi
    
    # Detect naming pattern from first file
    local first_file=$(echo "$all_fastq_files" | head -n 1)
    local detected_pattern=""
    
    if [[ "$first_file" =~ _R1_ ]]; then
        detected_pattern="_R1_/_R2_"
    elif [[ "$first_file" =~ _r1_ ]]; then
        detected_pattern="_r1_/_r2_"
    elif [[ "$first_file" =~ _1\. ]] || [[ "$first_file" =~ _1_ ]]; then
        detected_pattern="_1/_2"
    else
        log_error "Unable to detect read pairing pattern"
        log_info "Expected patterns: *_R1_*.fastq/*_R2_*.fastq, *_r1_*.fastq/*_r2_*.fastq, or *_1.fastq/*_2.fastq"
        echo ""
        log_info "Files found:"
        echo "$all_fastq_files" | head -n 5
        echo ""
        log_error "Please rename your files to follow one of the expected patterns"
        exit 1
    fi
    
    # Display detection results
    echo "========================================"
    echo "  Pairing Pattern Detection"
    echo "========================================"
    log_success "Detected pattern: ${detected_pattern}"
    echo ""
    
    # Find R1 files based on detected pattern
    local r1_files=""
    case "$detected_pattern" in
        "_R1_/_R2_")
            r1_files=$(find "$INPUT_DIR" -type f \( -name "*_R1_*.fastq" -o -name "*_R1_*.fq" -o -name "*_R1_*.fastq.gz" -o -name "*_R1_*.fq.gz" \) | sort)
            ;;
        "_r1_/_r2_")
            r1_files=$(find "$INPUT_DIR" -type f \( -name "*_r1_*.fastq" -o -name "*_r1_*.fq" -o -name "*_r1_*.fastq.gz" -o -name "*_r1_*.fq.gz" \) | sort)
            ;;
        "_1/_2")
            r1_files=$(find "$INPUT_DIR" -type f \( -name "*_1.fastq" -o -name "*_1.fq" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) | sort)
            ;;
    esac
    
    if [ -z "$r1_files" ]; then
        log_error "No R1 files found matching pattern: ${detected_pattern}"
        exit 1
    fi
    
    # Count total samples
    local total_r1=$(echo "$r1_files" | wc -l)
    local total_r2=$(echo "$all_fastq_files" | grep -E "(_R2_|_r2_|_2\.)" | wc -l)
    
    log_info "Found ${total_r1} R1 files and ${total_r2} R2 files"
    echo ""
    
    # Show sample list preview
    log_info "Sample preview (first 3):"
    local count=0
    for r1_file in $r1_files; do
        count=$((count + 1))
        if [ $count -le 3 ]; then
            local basename=$(basename "$r1_file")
            local sample_id=$(echo "$basename" | sed -E 's/_(R1|r1|1).*//g')
            local r2_file=$(find_r2_file "$r1_file")
            
            echo "  Sample ${count}: ${sample_id}"
            echo "    R1: $(basename $r1_file)"
            if [ -f "$r2_file" ]; then
                echo "    R2: $(basename $r2_file) ✓"
            else
                echo "    R2: NOT FOUND ✗"
            fi
        fi
    done
    
    echo ""
    echo "========================================"
    echo ""
    
    # Confirmation prompt with timeout
    log_warning "Please verify the detected pattern is correct"
    echo ""
    echo "If the pattern detection is INCORRECT:"
    echo "  1. Press Ctrl+C to abort"
    echo "  2. Rename your files to follow the expected pattern"
    echo "  3. Re-run this script"
    echo ""
    echo "Expected naming:"
    echo "  - Pattern: SampleID_R1_*.fastq / SampleID_R2_*.fastq"
    echo "  - Example: SRR1045522_R1_paired.fastq / SRR1045522_R2_paired.fastq"
    echo ""
    
    # Auto-continue after timeout
    local continue_timeout=30
    echo -ne "Continue with detected pattern? [Y/n] (auto-continue in ${continue_timeout}s): "
    read -t $continue_timeout -n 1 -r REPLY || true
    echo
    
    if [ -z "$REPLY" ]; then
        REPLY="y"
        log_info "Auto-continuing after ${continue_timeout}s timeout"
    fi
    
    if [[ $REPLY =~ ^[Nn]$ ]]; then
        log_info "Alignment cancelled by user"
        echo ""
        log_info "To fix naming issues:"
        echo "  1. Ensure R1 files contain: _R1_ (uppercase) or _r1_ (lowercase)"
        echo "  2. Ensure R2 files contain: _R2_ (uppercase) or _r2_ (lowercase)"
        echo "  3. Alternative: Use _1.fastq / _2.fastq naming"
        exit 0
    fi
    
    echo ""
    log_info "Starting alignment for ${total_r1} samples..."
    echo ""
    
    # Process each sample
    local processed_count=0
    local failed_count=0
    
    for r1_file in $r1_files; do
        local basename=$(basename "$r1_file")
        local sample_id=$(echo "$basename" | sed -E 's/_(R1|r1|1).*//g')
        local r2_file=$(find_r2_file "$r1_file")
        
        if [ ! -f "$r2_file" ]; then
            log_warning "R2 not found for $sample_id, skipping..."
            log_warning "Expected: $r2_file"
            failed_count=$((failed_count + 1))
            continue
        fi
        
        processed_count=$((processed_count + 1))
        
        echo ""
        log_info "========================================="
        log_info "Sample ${processed_count}/${total_r1}: ${sample_id}"
        log_info "========================================="
        log_info "Pattern: ${detected_pattern}"
        log_info "R1: $(basename $r1_file)"
        log_info "R2: $(basename $r2_file)"
        echo ""
        
        case "$ALIGNER" in
            star)
                align_star_pe "$sample_id" "$r1_file" "$r2_file" "" || failed_count=$((failed_count + 1))
                ;;
            hisat2)
                align_hisat2_pe "$sample_id" "$r1_file" "$r2_file" "" || failed_count=$((failed_count + 1))
                ;;
            both)
                log_info "Running alignment with BOTH aligners..."
                align_star_pe "$sample_id" "$r1_file" "$r2_file" "STAR" || failed_count=$((failed_count + 1))
                align_hisat2_pe "$sample_id" "$r1_file" "$r2_file" "HISAT2" || failed_count=$((failed_count + 1))
                ;;
        esac
    done
    
    echo ""
    log_success "Processed ${processed_count} of ${total_r1} samples"
    
    if [ $failed_count -gt 0 ]; then
        log_warning "${failed_count} samples skipped or failed"
    fi
    
    if [ $processed_count -eq 0 ]; then
        log_error "No valid sample pairs found and processed"
        exit 1
    fi
}

# Process samples from file
process_samples_from_file() {
    log_info "Processing samples from: $SAMPLE_FILE"
    
    if [ ! -f "$SAMPLE_FILE" ]; then
        log_error "Sample file not found: $SAMPLE_FILE"
        exit 1
    fi
    
    local sample_count=0
    
    while IFS= read -r sample_id; do
        [[ -z "$sample_id" || "$sample_id" =~ ^#.*$ ]] && continue
        
        sample_count=$((sample_count + 1))
        
        echo ""
        log_info "========================================="
        log_info "Sample ${sample_count}: ${sample_id}"
        log_info "========================================="
        
        # Find R1 file
        local r1_file=$(find "$INPUT_DIR" -type f \( -name "${sample_id}*_R1*.fastq" -o -name "${sample_id}*_R1*.fq" -o -name "${sample_id}*_R1*.fastq.gz" -o -name "${sample_id}*_1.fastq" -o -name "${sample_id}*_1.fq" \) | head -n 1)
        
        if [ -z "$r1_file" ] || [ ! -f "$r1_file" ]; then
            log_error "R1 file not found for: ${sample_id}"
            continue
        fi
        
        # Find R2 file
        local r2_file=$(find_r2_file "$r1_file")
        local pattern=$(detect_pattern "$r1_file")
        
        if [ ! -f "$r2_file" ]; then
            log_error "R2 not found for: ${sample_id} (pattern: ${pattern})"
            continue
        fi
        
        log_info "Pattern: ${pattern}"
        log_info "R1: $r1_file"
        log_info "R2: $r2_file"
        
        # Align
        case "$ALIGNER" in
            star)
                align_star_pe "$sample_id" "$r1_file" "$r2_file" ""
                ;;
            hisat2)
                align_hisat2_pe "$sample_id" "$r1_file" "$r2_file" ""
                ;;
            both)
                log_info "Running with BOTH aligners..."
                align_star_pe "$sample_id" "$r1_file" "$r2_file" "STAR"
                align_hisat2_pe "$sample_id" "$r1_file" "$r2_file" "HISAT2"
                ;;
        esac
        
    done < "$SAMPLE_FILE"
    
    log_success "Processed ${sample_count} samples"
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
  - Review alignment logs
  - Proceed with quantification

========================================
EOF
    
    log_success "Summary: $summary_file"
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

# Parse arguments
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
    
    # Validate required arguments
    if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
        log_error "Input and output directories are required"
        usage
    fi
    
    if [ ! -d "$INPUT_DIR" ]; then
        log_error "Input directory not found: $INPUT_DIR"
        exit 1
    fi
}

# Main
main() {
    echo "========================================"
    echo "  Alignment Module"
    echo "========================================"
    echo ""
    
    # Parse arguments
    parse_arguments "$@"
    
    # Check installation
    check_installation
    
    # Prompt for aligner if not specified
    if [ -z "$ALIGNER" ]; then
        prompt_aligner_selection
    fi
    
    # Verify selected aligner is available
    verify_aligner
    
    # Check reference genome
    check_reference
    
    # Check or build index
    check_or_build_index
    
    # Create output directories
    create_output_dirs
    
    # Log start
    local start_time=$(date +%s)
    log_info "Started: $(date)"
    log_info "Aligner: ${ALIGNER^^}"
    log_info "Threads: $THREADS"
    
    # Process samples
    echo ""
    if [ -n "$SAMPLE_FILE" ]; then
        process_samples_from_file
    else
        process_all_samples
    fi
    
    # Generate summary
    echo ""
    generate_summary
    
    # Calculate runtime
    local end_time=$(date +%s)
    local runtime=$((end_time - start_time))
    local hours=$((runtime / 3600))
    local minutes=$(((runtime % 3600) / 60))
    local seconds=$((runtime % 60))
    
    log_info "Completed: $(date)"
    log_info "Runtime: ${hours}h ${minutes}m ${seconds}s"
    
    # Display summary
    display_final_summary
}

# Run main
main "$@"
