#!/bin/bash

# Alignment Module Installation Script
# Aligns with QC module structure for seamless pipeline integration

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Tool versions and URLs
STAR_VERSION="2.7.11b"
STAR_URL="https://github.com/alexdobin/STAR/releases/download/${STAR_VERSION}/STAR_${STAR_VERSION}.zip"

HISAT2_VERSION="2.2.1"
HISAT2_URL="https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download/hisat2-${HISAT2_VERSION}-Linux_x86_64.zip"

# Reference URLs
REFERENCE_FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
REFERENCE_GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz"

# Defaults (aligned with QC module)
DEFAULT_INSTALL_DIR="${HOME}/softwares"
DEFAULT_REFERENCE_DIR="${HOME}/references/GRCh38_ensembl113"

INSTALL_BASE_DIR=""
REFERENCE_DIR=""
STAR_DIR=""
HISAT2_DIR=""
BIN_DIR=""

# Installation flags
INSTALL_STAR=""
INSTALL_HISAT2=""
INSTALL_SAMTOOLS=""
DOWNLOAD_REFERENCE=""

log_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

command_exists() { command -v "$1" >/dev/null 2>&1; }

detect_distro() {
    if [ -f /etc/os-release ]; then
        . /etc/os-release
        echo "$ID"
    elif [ -f /etc/redhat-release ]; then
        echo "redhat"
    else
        echo "unknown"
    fi
}

# Check if QC module is installed (for alignment detection)
check_qc_module() {
    local QC_CONFIG="${HOME}/qc/config/install_paths.conf"
    if [ -f "$QC_CONFIG" ]; then
        log_info "QC module detected, will use compatible paths"
        # Only use QC path if user hasn't provided custom path
        if [ -z "$INSTALL_BASE_DIR" ]; then
            source "$QC_CONFIG"
            if [ -n "$INSTALL_BASE_DIR" ]; then
                log_info "Using QC module installation directory: $INSTALL_BASE_DIR"
            fi
        fi
    fi
}

prompt_install_directory() {
    echo ""
    echo "========================================"
    echo "  Installation Directory"
    echo "========================================"
    echo ""
    log_info "Choose installation directory"
    echo ""
    echo "  1) ${HOME}/softwares (recommended, matches QC module)"
    echo "  2) /opt/alignment-tools (system-wide, requires sudo)"
    echo "  3) Custom directory"
    echo ""
    read -p "Enter choice [1-3] (default: 1): " DIR_CHOICE
    
    case "${DIR_CHOICE:-1}" in
        1) INSTALL_BASE_DIR="${HOME}/softwares" ;;
        2) INSTALL_BASE_DIR="/opt/alignment-tools" ;;
        3) read -p "Enter custom directory: " CUSTOM_DIR
           INSTALL_BASE_DIR="${CUSTOM_DIR}" ;;
        *) INSTALL_BASE_DIR="${DEFAULT_INSTALL_DIR}" ;;
    esac
    
    INSTALL_BASE_DIR="${INSTALL_BASE_DIR/#\~/$HOME}"
    INSTALL_BASE_DIR="${INSTALL_BASE_DIR%/}"
    
    STAR_DIR="${INSTALL_BASE_DIR}/STAR"
    HISAT2_DIR="${INSTALL_BASE_DIR}/HISAT2"
    BIN_DIR="${INSTALL_BASE_DIR}/bin"
    
    mkdir -p "${INSTALL_BASE_DIR}" "${BIN_DIR}"
    log_success "Installation directory: ${INSTALL_BASE_DIR}"
}

prompt_aligner_choice() {
    echo ""
    echo "========================================"
    echo "  Aligner Selection"
    echo "========================================"
    echo ""
    log_info "Which aligner(s) do you want to install?"
    echo ""
    echo "  1) STAR only (recommended for high-memory systems)"
    echo "     - Memory: 32+ GB"
    echo "     - Best for: Splice-aware alignment, isoform detection"
    echo ""
    echo "  2) HISAT2 only (recommended for limited resources)"
    echo "     - Memory: 8-16 GB"
    echo "     - Best for: Fast alignment, low memory usage"
    echo ""
    echo "  3) Both STAR and HISAT2 (flexibility)"
    echo "     - Choose aligner at runtime"
    echo ""
    read -p "Enter choice [1-3] (default: 3): " ALIGNER_CHOICE
    
    case "${ALIGNER_CHOICE:-3}" in
        1) 
            INSTALL_STAR="yes"
            INSTALL_HISAT2="no"
            log_info "Selected: STAR only"
            ;;
        2) 
            INSTALL_STAR="no"
            INSTALL_HISAT2="yes"
            log_info "Selected: HISAT2 only"
            ;;
        3) 
            INSTALL_STAR="yes"
            INSTALL_HISAT2="yes"
            log_info "Selected: Both aligners"
            ;;
        *) 
            INSTALL_STAR="yes"
            INSTALL_HISAT2="yes"
            ;;
    esac
}

prompt_reference_directory() {
    echo ""
    log_info "Reference genome directory"
    echo ""
    echo "  1) ${HOME}/references/GRCh38_ensembl113 (recommended)"
    echo "  2) Custom directory"
    echo ""
    read -p "Enter choice [1-2] (default: 1): " REF_CHOICE
    
    case "${REF_CHOICE:-1}" in
        1) REFERENCE_DIR="${DEFAULT_REFERENCE_DIR}" ;;
        2) read -p "Enter custom reference directory: " CUSTOM_REF
           REFERENCE_DIR="${CUSTOM_REF}" ;;
        *) REFERENCE_DIR="${DEFAULT_REFERENCE_DIR}" ;;
    esac
    
    REFERENCE_DIR="${REFERENCE_DIR/#\~/$HOME}"
    REFERENCE_DIR="${REFERENCE_DIR%/}"
    
    mkdir -p "${REFERENCE_DIR}"
    log_success "Reference directory: ${REFERENCE_DIR}"
}

check_star() {
    if [ -f "${STAR_DIR}/STAR" ] && "${STAR_DIR}/STAR" --version >/dev/null 2>&1; then
        log_success "STAR already installed"
        return 0
    fi
    return 1
}

install_star() {
    log_info "Installing STAR ${STAR_VERSION}..."
    
    cd /tmp
    if ! wget -q --show-progress "${STAR_URL}" -O "STAR.zip"; then
        log_error "Failed to download STAR"
        return 1
    fi
    
    unzip -q "STAR.zip"
    mkdir -p "${STAR_DIR}"
    
    # Find and copy STAR binary
    local STAR_BINARY=$(find . -name "STAR" -type f -executable 2>/dev/null | grep "Linux_x86_64" | head -n 1)
    if [ -z "$STAR_BINARY" ]; then
        log_error "STAR binary not found"
        return 1
    fi
    
    cp "$STAR_BINARY" "${STAR_DIR}/"
    ln -sf "${STAR_DIR}/STAR" "${BIN_DIR}/STAR"
    rm -rf STAR* "STAR.zip"
    
    if "${STAR_DIR}/STAR" --version >/dev/null 2>&1; then
        log_success "STAR installed successfully"
        return 0
    else
        log_error "STAR verification failed"
        return 1
    fi
}

check_hisat2() {
    if [ -f "${HISAT2_DIR}/hisat2" ] && "${HISAT2_DIR}/hisat2" --version >/dev/null 2>&1; then
        log_success "HISAT2 already installed"
        return 0
    fi
    return 1
}

install_hisat2() {
    log_info "Installing HISAT2 ${HISAT2_VERSION}..."
    
    cd /tmp
    if ! wget -q --show-progress "${HISAT2_URL}" -O "hisat2.zip"; then
        log_error "Failed to download HISAT2"
        return 1
    fi
    
    unzip -q "hisat2.zip"
    
    # FIXED: Find extracted directory - matches actual structure "hisat2-2.2.1"
    local HISAT2_EXTRACTED=$(find . -maxdepth 1 -type d -iname "hisat2-${HISAT2_VERSION}*" -o -iname "hisat2-*" | head -n 1)
    if [ -z "$HISAT2_EXTRACTED" ]; then
        log_error "HISAT2 extracted directory not found"
        log_info "Available directories:"
        ls -la /tmp/hisat2* || true
        return 1
    fi
    
    log_info "Found HISAT2 directory: $HISAT2_EXTRACTED"
    mkdir -p "${HISAT2_DIR}"
    cp -r "$HISAT2_EXTRACTED/"* "${HISAT2_DIR}/"
    
    ln -sf "${HISAT2_DIR}/hisat2" "${BIN_DIR}/hisat2"
    ln -sf "${HISAT2_DIR}/hisat2-build" "${BIN_DIR}/hisat2-build"
    ln -sf "${HISAT2_DIR}/hisat2-inspect" "${BIN_DIR}/hisat2-inspect"
    
    rm -rf hisat2* "hisat2.zip"
    
    if "${HISAT2_DIR}/hisat2" --version >/dev/null 2>&1; then
        log_success "HISAT2 installed successfully"
        return 0
    else
        log_error "HISAT2 verification failed"
        return 1
    fi
}

check_samtools() {
    if command_exists samtools; then
        log_success "samtools already installed"
        return 0
    fi
    return 1
}

install_samtools() {
    log_info "Installing samtools..."
    
    DISTRO=$(detect_distro)
    case "$DISTRO" in
        ubuntu|debian|linuxmint|mint)
            sudo apt-get update && sudo apt-get install -y samtools ;;
        centos|rhel|redhat|fedora)
            sudo yum install -y samtools ;;
        *)
            log_error "Unsupported distribution"
            return 1 ;;
    esac
    
    if command_exists samtools; then
        log_success "samtools installed"
        return 0
    fi
    return 1
}

check_dependencies() {
    log_info "Checking dependencies..."
    
    local MISSING=()
    command_exists wget || MISSING+=("wget")
    command_exists unzip || MISSING+=("unzip")
    command_exists gunzip || MISSING+=("gunzip")
    
    if [ ${#MISSING[@]} -gt 0 ]; then
        log_warning "Missing: ${MISSING[*]}"
        
        DISTRO=$(detect_distro)
        case "$DISTRO" in
            ubuntu|debian|linuxmint|mint)
                sudo apt-get update && sudo apt-get install -y wget unzip gzip ;;
            centos|rhel|redhat|fedora)
                sudo yum install -y wget unzip gzip ;;
        esac
    fi
    
    log_success "Dependencies ready"
}

download_reference() {
    log_info "Downloading reference genome (GRCh38 Ensembl 113)..."
    log_warning "This will download ~1 GB of data"
    
    cd "${REFERENCE_DIR}"
    
    # Download FASTA
    if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
        log_info "Downloading reference FASTA (~900 MB)..."
        if wget -q --show-progress "${REFERENCE_FASTA_URL}" -O Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz; then
            log_info "Extracting FASTA..."
            gunzip -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
            log_success "Reference FASTA ready"
        else
            log_error "Failed to download FASTA"
            return 1
        fi
    else
        log_success "Reference FASTA already exists"
    fi
    
    # Download GTF
    if [ ! -f "Homo_sapiens.GRCh38.113.gtf" ]; then
        log_info "Downloading reference GTF (~50 MB)..."
        if wget -q --show-progress "${REFERENCE_GTF_URL}" -O Homo_sapiens.GRCh38.113.gtf.gz; then
            log_info "Extracting GTF..."
            gunzip -f Homo_sapiens.GRCh38.113.gtf.gz
            log_success "Reference GTF ready"
        else
            log_error "Failed to download GTF"
            return 1
        fi
    else
        log_success "Reference GTF already exists"
    fi
    
    log_success "Reference genome ready: ${REFERENCE_DIR}"
}

update_path() {
    local SHELL_RC="${HOME}/.bashrc"
    
    if grep -q "# Alignment Module" "$SHELL_RC" 2>/dev/null; then
        log_warning "PATH already configured"
        return 0
    fi
    
    echo "" >> "$SHELL_RC"
    echo "# Alignment Module - added by installer" >> "$SHELL_RC"
    echo "export PATH=\"${BIN_DIR}:\$PATH\"" >> "$SHELL_RC"
    export PATH="${BIN_DIR}:$PATH"
    
    log_success "PATH updated"
}

save_config() {
    mkdir -p "${SCRIPT_DIR}/config"
    
    cat > "${SCRIPT_DIR}/config/install_paths.conf" << EOF
# Alignment Module Configuration
# Generated: $(date)

INSTALL_BASE_DIR="${INSTALL_BASE_DIR}"
BIN_DIR="${BIN_DIR}"
STAR_DIR="${STAR_DIR}"
HISAT2_DIR="${HISAT2_DIR}"
STAR_INSTALLED="${INSTALL_STAR}"
HISAT2_INSTALLED="${INSTALL_HISAT2}"
STAR_BIN="${BIN_DIR}/STAR"
HISAT2_BIN="${BIN_DIR}/hisat2"
HISAT2_BUILD_BIN="${BIN_DIR}/hisat2-build"
SAMTOOLS_BIN="$(command -v samtools 2>/dev/null || echo '')"
EOF
    
    cat > "${SCRIPT_DIR}/config/reference_paths.conf" << EOF
# Reference Genome Configuration
# Generated: $(date)

REFERENCE_DIR="${REFERENCE_DIR}"
REFERENCE_FASTA="${REFERENCE_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REFERENCE_GTF="${REFERENCE_DIR}/Homo_sapiens.GRCh38.113.gtf"
STAR_INDEX_DIR="${REFERENCE_DIR}/STAR_index"
HISAT2_INDEX_DIR="${REFERENCE_DIR}/HISAT2_index"
HISAT2_INDEX_PREFIX="${REFERENCE_DIR}/HISAT2_index/genome"
EOF
    
    log_success "Configuration saved"
}

main() {
    echo "========================================"
    echo "  Alignment Module Installation"
    echo "========================================"
    echo ""
    
    DISTRO=$(detect_distro)
    log_info "Detected OS: ${DISTRO}"
    
    # Parse command-line arguments FIRST
    # (this is done before main() is called, see bottom of script)
    
    # Check for QC module only if no custom path provided
    check_qc_module
    
    # Prompt for directories if not set
    if [ -z "$INSTALL_BASE_DIR" ]; then
        prompt_install_directory
    else
        # Ensure derived directories are set
        STAR_DIR="${INSTALL_BASE_DIR}/STAR"
        HISAT2_DIR="${INSTALL_BASE_DIR}/HISAT2"
        BIN_DIR="${INSTALL_BASE_DIR}/bin"
        mkdir -p "${INSTALL_BASE_DIR}" "${BIN_DIR}"
        log_info "Using installation directory: ${INSTALL_BASE_DIR}"
    fi
    
    # Prompt for aligner choice
    if [ -z "$INSTALL_STAR" ] && [ -z "$INSTALL_HISAT2" ]; then
        prompt_aligner_choice
    fi
    
    [ -z "$INSTALL_SAMTOOLS" ] && INSTALL_SAMTOOLS="yes"
    
    # Check dependencies
    check_dependencies
    
    # Install selected aligners
    echo ""
    log_info "Installing selected tools..."
    echo ""
    
    if [ "$INSTALL_STAR" = "yes" ]; then
        check_star || install_star || log_error "STAR installation failed"
    else
        log_info "Skipping STAR"
    fi
    
    echo ""
    
    if [ "$INSTALL_HISAT2" = "yes" ]; then
        check_hisat2 || install_hisat2 || log_error "HISAT2 installation failed"
    else
        log_info "Skipping HISAT2"
    fi
    
    echo ""
    
    if [ "$INSTALL_SAMTOOLS" = "yes" ]; then
        check_samtools || install_samtools || log_warning "samtools issues"
    fi
    
    # Download reference
    if [ -z "$DOWNLOAD_REFERENCE" ]; then
        echo ""
        read -p "Download reference genome (GRCh38 Ensembl 113)? [Y/n] " -n 1 -r
        echo
        DOWNLOAD_REFERENCE=$([[ ! $REPLY =~ ^[Nn]$ ]] && echo "yes" || echo "no")
    fi
    
    if [ "$DOWNLOAD_REFERENCE" = "yes" ]; then
        [ -z "$REFERENCE_DIR" ] && prompt_reference_directory
        download_reference || log_warning "Reference incomplete"
    else
        REFERENCE_DIR="${DEFAULT_REFERENCE_DIR}"
    fi
    
    # Finalize
    echo ""
    update_path
    save_config
    
    echo ""
    echo "========================================"
    echo "  Installation Complete"
    echo "========================================"
    echo ""
    log_success "Installation successful!"
    echo ""
    log_info "Installed tools:"
    [ "$INSTALL_STAR" = "yes" ] && echo "  ✓ STAR ${STAR_VERSION}"
    [ "$INSTALL_HISAT2" = "yes" ] && echo "  ✓ HISAT2 ${HISAT2_VERSION}"
    [ "$INSTALL_SAMTOOLS" = "yes" ] && echo "  ✓ samtools"
    echo ""
    log_info "Paths (aligned with QC module):"
    log_info "  Tools: ${INSTALL_BASE_DIR}"
    [ "$DOWNLOAD_REFERENCE" = "yes" ] && log_info "  Reference: ${REFERENCE_DIR}"
    echo ""
    log_warning "Restart terminal or run: source ~/.bashrc"
    echo ""
    log_info "Next: ./run_alignment.sh -i qc_results/trimmed/PE/ -o alignment_results/ --aligner [star|hisat2]"
    echo ""
}

# Parse arguments BEFORE calling main
while [[ $# -gt 0 ]]; do
    case $1 in
        --install-dir)
            INSTALL_BASE_DIR="$2"
            INSTALL_BASE_DIR="${INSTALL_BASE_DIR/#\~/$HOME}"
            INSTALL_BASE_DIR="${INSTALL_BASE_DIR%/}"
            STAR_DIR="${INSTALL_BASE_DIR}/STAR"
            HISAT2_DIR="${INSTALL_BASE_DIR}/HISAT2"
            BIN_DIR="${INSTALL_BASE_DIR}/bin"
            mkdir -p "${INSTALL_BASE_DIR}" "${BIN_DIR}"
            shift 2 ;;
        --reference-dir)
            REFERENCE_DIR="$2"
            REFERENCE_DIR="${REFERENCE_DIR/#\~/$HOME}"
            REFERENCE_DIR="${REFERENCE_DIR%/}"
            mkdir -p "${REFERENCE_DIR}"
            shift 2 ;;
        --aligner)
            case "$2" in
                star) INSTALL_STAR="yes"; INSTALL_HISAT2="no" ;;
                hisat2) INSTALL_STAR="no"; INSTALL_HISAT2="yes" ;;
                both) INSTALL_STAR="yes"; INSTALL_HISAT2="yes" ;;
                *) log_error "Invalid aligner: $2"; exit 1 ;;
            esac
            shift 2 ;;
        --download-reference) DOWNLOAD_REFERENCE="yes"; shift ;;
        --skip-reference) DOWNLOAD_REFERENCE="no"; shift ;;
        -h|--help)
            cat << EOF
Usage: $0 [options]

Options:
  --install-dir <path>              Installation directory (default: ~/softwares)
  --reference-dir <path>            Reference genome directory
  --aligner <star|hisat2|both>      Which aligner(s) to install
  --download-reference              Download reference genome
  --skip-reference                  Skip reference download
  -h, --help                        Show help

Examples:
  # Install only STAR
  $0 --aligner star

  # Install both, same dir as QC module
  $0 --install-dir ~/softwares --aligner both

  # Quick install with reference
  $0 --aligner hisat2 --download-reference
EOF
            exit 0 ;;
        *)
            log_error "Unknown option: $1"
            exit 1 ;;
    esac
done

main
