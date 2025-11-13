#!/bin/bash

# Alignment Module Installation Script
# Robust system Python detection with fallback to Miniconda

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

# Defaults
DEFAULT_INSTALL_DIR="${HOME}/softwares"
DEFAULT_REFERENCE_DIR="${HOME}/references/GRCh38_ensembl113"

INSTALL_BASE_DIR=""
REFERENCE_DIR=""
STAR_DIR=""
HISAT2_DIR=""
BIN_DIR=""

# Python/Conda settings
CONDA_DIR=""
CONDA_ENV_NAME="bioinfo_python3"
PYTHON_VERSION="3.9"
USE_SYSTEM_PYTHON="no"
PYTHON_BIN=""

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

# Robust system Python check
check_system_python() {
    log_info "Checking for system Python 3..."
    
    # Try python3 first
    if command_exists python3; then
        local PYTHON_VERSION_OUTPUT=$(python3 --version 2>&1)
        if [[ "$PYTHON_VERSION_OUTPUT" =~ Python\ ([0-9]+)\.([0-9]+) ]]; then
            local MAJOR="${BASH_REMATCH[1]}"
            local MINOR="${BASH_REMATCH[2]}"
            
            if [ "$MAJOR" -eq 3 ] && [ "$MINOR" -ge 6 ]; then
                PYTHON_BIN=$(command -v python3)
                USE_SYSTEM_PYTHON="yes"
                log_success "System Python 3 found: $PYTHON_VERSION_OUTPUT"
                log_info "Using system Python: $PYTHON_BIN"
                return 0
            fi
        fi
    fi
    
    # Try python (might be python3)
    if command_exists python; then
        local PYTHON_VERSION_OUTPUT=$(python --version 2>&1)
        if [[ "$PYTHON_VERSION_OUTPUT" =~ Python\ ([0-9]+)\.([0-9]+) ]]; then
            local MAJOR="${BASH_REMATCH[1]}"
            local MINOR="${BASH_REMATCH[2]}"
            
            if [ "$MAJOR" -eq 3 ] && [ "$MINOR" -ge 6 ]; then
                PYTHON_BIN=$(command -v python)
                USE_SYSTEM_PYTHON="yes"
                log_success "System Python 3 found: $PYTHON_VERSION_OUTPUT"
                log_info "Using system Python: $PYTHON_BIN"
                return 0
            fi
        fi
    fi
    
    log_warning "No suitable system Python 3 found (requires >= 3.6)"
    log_info "Will install Miniconda with Python ${PYTHON_VERSION}"
    USE_SYSTEM_PYTHON="no"
    PYTHON_BIN=""
    return 1
}

check_qc_module() {
    local QC_CONFIG="${HOME}/qc/config/install_paths.conf"
    if [ -f "$QC_CONFIG" ]; then
        log_info "QC module detected"
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
    echo "  1) ${HOME}/softwares (recommended)"
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
    CONDA_DIR="${INSTALL_BASE_DIR}/miniconda"
    
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
    echo "  1) STAR only (high memory ~32GB)"
    echo "  2) HISAT2 only (low memory ~8GB)"
    echo "  3) Both (choose at runtime)"
    echo ""
    read -p "Enter choice [1-3] (default: 3): " ALIGNER_CHOICE
    
    case "${ALIGNER_CHOICE:-3}" in
        1) INSTALL_STAR="yes"; INSTALL_HISAT2="no" ;;
        2) INSTALL_STAR="no"; INSTALL_HISAT2="yes" ;;
        3) INSTALL_STAR="yes"; INSTALL_HISAT2="yes" ;;
        *) INSTALL_STAR="yes"; INSTALL_HISAT2="yes" ;;
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
    [ -f "${STAR_DIR}/STAR" ] && "${STAR_DIR}/STAR" --version >/dev/null 2>&1
}

install_star() {
    log_info "Installing STAR ${STAR_VERSION}..."
    cd /tmp
    wget -q --show-progress "${STAR_URL}" -O "STAR.zip" || return 1
    unzip -q "STAR.zip"
    mkdir -p "${STAR_DIR}"
    
    local STAR_BINARY=$(find . -name "STAR" -type f -executable 2>/dev/null | grep "Linux_x86_64" | head -n 1)
    [ -z "$STAR_BINARY" ] && return 1
    
    cp "$STAR_BINARY" "${STAR_DIR}/"
    ln -sf "${STAR_DIR}/STAR" "${BIN_DIR}/STAR"
    rm -rf STAR* "STAR.zip"
    
    "${STAR_DIR}/STAR" --version >/dev/null 2>&1
}

check_hisat2() {
    [ -f "${HISAT2_DIR}/hisat2" ] && "${HISAT2_DIR}/hisat2" --version >/dev/null 2>&1
}

install_hisat2() {
    log_info "Installing HISAT2 ${HISAT2_VERSION}..."
    cd /tmp
    wget -q --show-progress "${HISAT2_URL}" -O "hisat2.zip" || return 1
    unzip -q "hisat2.zip"
    
    local HISAT2_EXTRACTED=$(find . -maxdepth 1 -type d -name "hisat2-${HISAT2_VERSION}*" 2>/dev/null | head -n 1)
    [ -z "$HISAT2_EXTRACTED" ] && return 1
    
    mkdir -p "${HISAT2_DIR}"
    cp -r "$HISAT2_EXTRACTED/"* "${HISAT2_DIR}/"
    
    ln -sf "${HISAT2_DIR}/hisat2" "${BIN_DIR}/hisat2"
    ln -sf "${HISAT2_DIR}/hisat2-build" "${BIN_DIR}/hisat2-build"
    ln -sf "${HISAT2_DIR}/hisat2-inspect" "${BIN_DIR}/hisat2-inspect"
    
    rm -rf hisat2* "hisat2.zip"
    "${HISAT2_DIR}/hisat2" --version >/dev/null 2>&1
}

# Setup Python environment (only if system Python not suitable)
setup_python_environment() {
    # Check if system Python is suitable
    if check_system_python; then
        return 0
    fi
    
    log_info "Setting up Python environment via Miniconda..."
    
    if [ ! -d "$CONDA_DIR" ]; then
        log_info "Downloading Miniconda..."
        local CONDA_INSTALLER="${INSTALL_BASE_DIR}/miniconda_installer.sh"
        wget -q --show-progress https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O "$CONDA_INSTALLER" || return 1
        
        bash "$CONDA_INSTALLER" -b -p "$CONDA_DIR"
        rm -f "$CONDA_INSTALLER"
        log_success "Miniconda installed"
    fi
    
    # Initialize conda
    export PATH="$CONDA_DIR/bin:$PATH"
    [ ! -f "${CONDA_DIR}/etc/profile.d/conda.sh" ] && return 1
    source "${CONDA_DIR}/etc/profile.d/conda.sh"
    
    # Create environment
    if ! conda env list | grep -qw "^${CONDA_ENV_NAME}"; then
        log_info "Creating conda environment '$CONDA_ENV_NAME'..."
        conda create -y -n "$CONDA_ENV_NAME" python="$PYTHON_VERSION" 2>&1 | grep -E "(Collecting|Solving|done)" || true
    fi
    
    # Verify activation
    set +u
    if conda activate "$CONDA_ENV_NAME" 2>/dev/null; then
        python --version && log_success "Python environment ready"
        conda deactivate
    else
        log_error "Failed to activate conda environment"
        set -u
        return 1
    fi
    set -u
    return 0
}

check_samtools() {
    command_exists samtools
}

install_samtools() {
    log_info "Installing samtools..."
    DISTRO=$(detect_distro)
    case "$DISTRO" in
        ubuntu|debian|linuxmint|mint)
            sudo apt-get update && sudo apt-get install -y samtools ;;
        centos|rhel|redhat|fedora)
            sudo yum install -y samtools ;;
        *) return 1 ;;
    esac
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
}

download_reference() {
    log_info "Downloading reference genome..."
    cd "${REFERENCE_DIR}"
    
    if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
        wget -q --show-progress "${REFERENCE_FASTA_URL}" -O Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        gunzip -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    fi
    
    if [ ! -f "Homo_sapiens.GRCh38.113.gtf" ]; then
        wget -q --show-progress "${REFERENCE_GTF_URL}" -O Homo_sapiens.GRCh38.113.gtf.gz
        gunzip -f Homo_sapiens.GRCh38.113.gtf.gz
    fi
}

update_path() {
    local SHELL_RC="${HOME}/.bashrc"
    grep -q "# Alignment Module" "$SHELL_RC" 2>/dev/null && return 0
    
    echo "" >> "$SHELL_RC"
    echo "# Alignment Module" >> "$SHELL_RC"
    echo "export PATH=\"${BIN_DIR}:\$PATH\"" >> "$SHELL_RC"
    
    if [ "$USE_SYSTEM_PYTHON" = "no" ] && [ -d "$CONDA_DIR" ]; then
        echo "export PATH=\"${CONDA_DIR}/bin:\$PATH\"" >> "$SHELL_RC"
    fi
    export PATH="${BIN_DIR}:$PATH"
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

# Python environment
USE_SYSTEM_PYTHON="${USE_SYSTEM_PYTHON}"
PYTHON_BIN="${PYTHON_BIN}"
CONDA_DIR="${CONDA_DIR}"
CONDA_ENV_NAME="${CONDA_ENV_NAME}"
EOF
    
    cat > "${SCRIPT_DIR}/config/reference_paths.conf" << EOF
# Reference Genome Configuration
REFERENCE_DIR="${REFERENCE_DIR}"
REFERENCE_FASTA="${REFERENCE_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REFERENCE_GTF="${REFERENCE_DIR}/Homo_sapiens.GRCh38.113.gtf"
STAR_INDEX_DIR="${REFERENCE_DIR}/STAR_index"
HISAT2_INDEX_DIR="${REFERENCE_DIR}/HISAT2_index"
HISAT2_INDEX_PREFIX="${REFERENCE_DIR}/HISAT2_index/genome"
EOF
}

main() {
    echo "========================================"
    echo "  Alignment Module Installation"
    echo "========================================"
    echo ""
    
    check_qc_module
    
    [ -z "$INSTALL_BASE_DIR" ] && prompt_install_directory || {
        STAR_DIR="${INSTALL_BASE_DIR}/STAR"
        HISAT2_DIR="${INSTALL_BASE_DIR}/HISAT2"
        BIN_DIR="${INSTALL_BASE_DIR}/bin"
        CONDA_DIR="${INSTALL_BASE_DIR}/miniconda"
        mkdir -p "${INSTALL_BASE_DIR}" "${BIN_DIR}"
    }
    
    [ -z "$INSTALL_STAR" ] && [ -z "$INSTALL_HISAT2" ] && prompt_aligner_choice
    [ -z "$INSTALL_SAMTOOLS" ] && INSTALL_SAMTOOLS="yes"
    
    check_dependencies
    
    echo ""
    [ "$INSTALL_STAR" = "yes" ] && { check_star && log_success "STAR already installed" || install_star && log_success "STAR installed"; }
    
    if [ "$INSTALL_HISAT2" = "yes" ]; then
        check_hisat2 && log_success "HISAT2 already installed" || install_hisat2 && log_success "HISAT2 installed"
        setup_python_environment
    fi
    
    [ "$INSTALL_SAMTOOLS" = "yes" ] && { check_samtools && log_success "samtools already installed" || install_samtools; }
    
    if [ -z "$DOWNLOAD_REFERENCE" ]; then
        echo ""
        read -p "Download reference genome? [Y/n] " -n 1 -r
        echo
        DOWNLOAD_REFERENCE=$([[ ! $REPLY =~ ^[Nn]$ ]] && echo "yes" || echo "no")
    fi
    
    [ "$DOWNLOAD_REFERENCE" = "yes" ] && { [ -z "$REFERENCE_DIR" ] && prompt_reference_directory; download_reference; } || REFERENCE_DIR="${DEFAULT_REFERENCE_DIR}"
    
    update_path
    save_config
    
    echo ""
    echo "========================================"
    echo "  Installation Complete"
    echo "========================================"
    echo ""
    log_success "Installation successful!"
    echo ""
    [ "$INSTALL_STAR" = "yes" ] && echo "  ✓ STAR ${STAR_VERSION}"
    [ "$INSTALL_HISAT2" = "yes" ] && echo "  ✓ HISAT2 ${HISAT2_VERSION}"
    [ "$USE_SYSTEM_PYTHON" = "yes" ] && echo "  ✓ Python (system: $PYTHON_BIN)" || echo "  ✓ Python (conda: $CONDA_ENV_NAME)"
    echo ""
    log_warning "Restart terminal or run: source ~/.bashrc"
    echo ""
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --install-dir)
            INSTALL_BASE_DIR="$2"
            INSTALL_BASE_DIR="${INSTALL_BASE_DIR/#\~/$HOME}"
            shift 2 ;;
        --reference-dir)
            REFERENCE_DIR="$2"
            shift 2 ;;
        --aligner)
            case "$2" in
                star) INSTALL_STAR="yes"; INSTALL_HISAT2="no" ;;
                hisat2) INSTALL_STAR="no"; INSTALL_HISAT2="yes" ;;
                both) INSTALL_STAR="yes"; INSTALL_HISAT2="yes" ;;
            esac
            shift 2 ;;
        --download-reference) DOWNLOAD_REFERENCE="yes"; shift ;;
        --skip-reference) DOWNLOAD_REFERENCE="no"; shift ;;
        -h|--help) echo "Usage: $0 [--install-dir <path>] [--aligner star|hisat2|both] [--download-reference]"; exit 0 ;;
        *) log_error "Unknown option: $1"; exit 1 ;;
    esac
done

main
