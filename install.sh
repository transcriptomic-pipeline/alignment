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
CONDA_DIR=""
CONDA_ENV_NAME="hisat2_env"

# Python detection
USE_SYSTEM_PYTHON="no"
PYTHON_BIN=""

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

# Check for system Python 3 (>= 3.6)
check_system_python() {
    log_info "Checking for system Python 3..."
    
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
    log_info "Will need to install Miniconda if HISAT2 is selected"
    USE_SYSTEM_PYTHON="no"
    PYTHON_BIN=""
    return 1
}

# Create python symlink if needed (for HISAT2 compatibility)
setup_python_symlink() {
    if [ "$USE_SYSTEM_PYTHON" = "yes" ]; then
        if ! command_exists python; then
            log_warning "'python' command not found, creating symlink..."
            mkdir -p "${BIN_DIR}"
            
            if [ -n "$PYTHON_BIN" ]; then
                ln -sf "$PYTHON_BIN" "${BIN_DIR}/python"
                log_success "Created python symlink: ${BIN_DIR}/python -> $PYTHON_BIN"
            fi
        else
            log_info "'python' command already exists"
        fi
    fi
}

# Check if QC module is already installed
check_qc_module() {
    local QC_CONFIG="${HOME}/qc/config/install_paths.conf"
    
    if [ -f "$QC_CONFIG" ]; then
        log_info "QC module detected at ${HOME}/qc"
        
        # Try to use same installation directory
        if [ -z "$INSTALL_BASE_DIR" ]; then
            source "$QC_CONFIG" 2>/dev/null || true
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
        3) 
            read -p "Enter custom directory: " CUSTOM_DIR
            INSTALL_BASE_DIR="${CUSTOM_DIR}"
            ;;
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
        2) 
            read -p "Enter custom reference directory: " CUSTOM_REF
            REFERENCE_DIR="${CUSTOM_REF}"
            ;;
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
    if [ -z "$STAR_BINARY" ]; then
        log_error "STAR binary not found in downloaded package"
        return 1
    fi
    
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
    if [ -z "$HISAT2_EXTRACTED" ]; then
        log_error "HISAT2 directory not found in downloaded package"
        return 1
    fi
    
    mkdir -p "${HISAT2_DIR}"
    cp -r "$HISAT2_EXTRACTED/"* "${HISAT2_DIR}/"
    
    ln -sf "${HISAT2_DIR}/hisat2" "${BIN_DIR}/hisat2"
    ln -sf "${HISAT2_DIR}/hisat2-build" "${BIN_DIR}/hisat2-build"
    ln -sf "${HISAT2_DIR}/hisat2-inspect" "${BIN_DIR}/hisat2-inspect"
    
    rm -rf hisat2* "hisat2.zip"
    
    # Patch HISAT2 to use python3 if system Python is being used
    if [ "$USE_SYSTEM_PYTHON" = "yes" ]; then
        patch_hisat2_python
    fi
    
    "${HISAT2_DIR}/hisat2" --version >/dev/null 2>&1
}

# Patch HISAT2 scripts to use python3 instead of python
patch_hisat2_python() {
    log_info "Patching HISAT2 scripts to use python3..."
    
    local HISAT2_SCRIPTS=(
        "${HISAT2_DIR}/hisat2"
        "${HISAT2_DIR}/hisat2-build"
        "${HISAT2_DIR}/hisat2-inspect"
        "${HISAT2_DIR}/hisat2-align-s"
        "${HISAT2_DIR}/hisat2-align-l"
        "${HISAT2_DIR}/hisat2-build-s"
        "${HISAT2_DIR}/hisat2-build-l"
    )
    
    for script in "${HISAT2_SCRIPTS[@]}"; do
        if [ -f "$script" ]; then
            # Check if script uses '#!/usr/bin/env python'
            if head -n 1 "$script" | grep -q "#!/usr/bin/env python$"; then
                log_info "Patching: $(basename $script)"
                # Replace first line with python3
                sed -i '1s|#!/usr/bin/env python$|#!/usr/bin/env python3|' "$script"
            fi
        fi
    done
    
    log_success "HISAT2 scripts patched to use python3"
}

check_samtools() {
    command_exists samtools
}

install_samtools() {
    log_info "Installing samtools..."
    
    local DISTRO=$(detect_distro)
    
    case "$DISTRO" in
        ubuntu|debian|linuxmint|mint)
            sudo apt-get update && sudo apt-get install -y samtools
            ;;
        centos|rhel|redhat|fedora)
            sudo yum install -y samtools
            ;;
        *)
            log_error "Unsupported distribution: $DISTRO"
            return 1
            ;;
    esac
}

check_dependencies() {
    log_info "Checking dependencies..."
    
    local MISSING=()
    command_exists wget || MISSING+=("wget")
    command_exists unzip || MISSING+=("unzip")
    command_exists gunzip || MISSING+=("gunzip")
    
    if [ ${#MISSING[@]} -gt 0 ]; then
        log_warning "Missing dependencies: ${MISSING[*]}"
        log_info "Installing missing dependencies..."
        
        local DISTRO=$(detect_distro)
        case "$DISTRO" in
            ubuntu|debian|linuxmint|mint)
                sudo apt-get update && sudo apt-get install -y wget unzip gzip
                ;;
            centos|rhel|redhat|fedora)
                sudo yum install -y wget unzip gzip
                ;;
        esac
    fi
    
    log_success "All dependencies available"
}

download_reference() {
    log_info "Downloading reference genome..."
    
    cd "${REFERENCE_DIR}"
    
    if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
        log_info "Downloading reference FASTA..."
        wget -q --show-progress "${REFERENCE_FASTA_URL}" -O Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        gunzip -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        log_success "Reference FASTA downloaded"
    else
        log_info "Reference FASTA already exists"
    fi
    
    if [ ! -f "Homo_sapiens.GRCh38.113.gtf" ]; then
        log_info "Downloading reference GTF..."
        wget -q --show-progress "${REFERENCE_GTF_URL}" -O Homo_sapiens.GRCh38.113.gtf.gz
        gunzip -f Homo_sapiens.GRCh38.113.gtf.gz
        log_success "Reference GTF downloaded"
    else
        log_info "Reference GTF already exists"
    fi
}

update_path() {
    local SHELL_RC="${HOME}/.bashrc"
    
    if grep -q "# Alignment Module" "$SHELL_RC" 2>/dev/null; then
        log_info "PATH already configured in .bashrc"
        return 0
    fi
    
    echo "" >> "$SHELL_RC"
    echo "# Alignment Module" >> "$SHELL_RC"
    echo "export PATH=\"${BIN_DIR}:\$PATH\"" >> "$SHELL_RC"
    
    if [ "$USE_SYSTEM_PYTHON" = "no" ] && [ -d "$CONDA_DIR" ]; then
        echo "export PATH=\"${CONDA_DIR}/bin:\$PATH\"" >> "$SHELL_RC"
    fi
    
    export PATH="${BIN_DIR}:$PATH"
    
    log_success "PATH updated in .bashrc"
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

# Python configuration
USE_SYSTEM_PYTHON="${USE_SYSTEM_PYTHON}"
PYTHON_BIN="${PYTHON_BIN}"
CONDA_DIR="${CONDA_DIR}"
CONDA_ENV_NAME="${CONDA_ENV_NAME}"
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
    
    log_success "Configuration saved to: ${SCRIPT_DIR}/config/"
}

main() {
    echo "========================================"
    echo "  Alignment Module Installation"
    echo "========================================"
    echo ""
    
    # Check system Python FIRST
    check_system_python
    
    # Check for QC module
    check_qc_module
    
    # Prompt for installation directory if not set
    if [ -z "$INSTALL_BASE_DIR" ]; then
        prompt_install_directory
    else
        STAR_DIR="${INSTALL_BASE_DIR}/STAR"
        HISAT2_DIR="${INSTALL_BASE_DIR}/HISAT2"
        BIN_DIR="${INSTALL_BASE_DIR}/bin"
        CONDA_DIR="${INSTALL_BASE_DIR}/miniconda"
        mkdir -p "${INSTALL_BASE_DIR}" "${BIN_DIR}"
        log_info "Using installation directory: $INSTALL_BASE_DIR"
    fi
    
    # Setup python symlink if using system Python
    setup_python_symlink
    
    # Prompt for aligner selection if not set
    if [ -z "$INSTALL_STAR" ] && [ -z "$INSTALL_HISAT2" ]; then
        prompt_aligner_choice
    fi
    
    # Set samtools installation flag
    [ -z "$INSTALL_SAMTOOLS" ] && INSTALL_SAMTOOLS="yes"
    
    # Check and install dependencies
    check_dependencies
    
    echo ""
    log_info "Starting installation..."
    echo ""
    
    # Install STAR if requested
    if [ "$INSTALL_STAR" = "yes" ]; then
        if check_star; then
            log_success "STAR already installed"
        else
            if install_star; then
                log_success "STAR installed successfully"
            else
                log_error "STAR installation failed"
            fi
        fi
    fi
    
    # Install HISAT2 if requested
    if [ "$INSTALL_HISAT2" = "yes" ]; then
        if check_hisat2; then
            log_success "HISAT2 already installed"
        else
            if install_hisat2; then
                log_success "HISAT2 installed successfully"
            else
                log_error "HISAT2 installation failed"
            fi
        fi
    fi
    
    # Install samtools if requested
    if [ "$INSTALL_SAMTOOLS" = "yes" ]; then
        if check_samtools; then
            log_success "samtools already installed"
        else
            install_samtools
        fi
    fi
    
    # Download reference genome
    if [ -z "$DOWNLOAD_REFERENCE" ]; then
        echo ""
        read -p "Download reference genome? [Y/n] " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Nn]$ ]]; then
            DOWNLOAD_REFERENCE="yes"
        else
            DOWNLOAD_REFERENCE="no"
        fi
    fi
    
    if [ "$DOWNLOAD_REFERENCE" = "yes" ]; then
        [ -z "$REFERENCE_DIR" ] && prompt_reference_directory
        download_reference
    else
        REFERENCE_DIR="${DEFAULT_REFERENCE_DIR}"
        log_info "Skipping reference download"
    fi
    
    # Update PATH
    update_path
    
    # Save configuration
    save_config
    
    echo ""
    echo "========================================"
    echo "  Installation Complete"
    echo "========================================"
    echo ""
    log_success "Installation successful!"
    echo ""
    
    # Display installed tools
    [ "$INSTALL_STAR" = "yes" ] && echo "  ✓ STAR ${STAR_VERSION}"
    [ "$INSTALL_HISAT2" = "yes" ] && echo "  ✓ HISAT2 ${HISAT2_VERSION}"
    
    if [ "$USE_SYSTEM_PYTHON" = "yes" ]; then
        echo "  ✓ Python (system: $PYTHON_BIN)"
    else
        echo "  ✓ Python (conda: $CONDA_ENV_NAME)"
    fi
    
    echo ""
    log_info "Installation directory: ${INSTALL_BASE_DIR}"
    [ "$DOWNLOAD_REFERENCE" = "yes" ] && log_info "Reference directory: ${REFERENCE_DIR}"
    echo ""
    log_warning "Restart terminal or run: source ~/.bashrc"
    echo ""
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --install-dir)
            INSTALL_BASE_DIR="$2"
            INSTALL_BASE_DIR="${INSTALL_BASE_DIR/#\~/$HOME}"
            INSTALL_BASE_DIR="${INSTALL_BASE_DIR%/}"
            shift 2
            ;;
        --reference-dir)
            REFERENCE_DIR="$2"
            REFERENCE_DIR="${REFERENCE_DIR/#\~/$HOME}"
            REFERENCE_DIR="${REFERENCE_DIR%/}"
            shift 2
            ;;
        --aligner)
            case "$2" in
                star)
                    INSTALL_STAR="yes"
                    INSTALL_HISAT2="no"
                    ;;
                hisat2)
                    INSTALL_STAR="no"
                    INSTALL_HISAT2="yes"
                    ;;
                both)
                    INSTALL_STAR="yes"
                    INSTALL_HISAT2="yes"
                    ;;
                *)
                    log_error "Invalid aligner: $2"
                    exit 1
                    ;;
            esac
            shift 2
            ;;
        --download-reference)
            DOWNLOAD_REFERENCE="yes"
            shift
            ;;
        --skip-reference)
            DOWNLOAD_REFERENCE="no"
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --install-dir <path>        Installation directory (default: ~/softwares)"
            echo "  --reference-dir <path>      Reference genome directory"
            echo "  --aligner <star|hisat2|both> Which aligner(s) to install"
            echo "  --download-reference        Download reference genome"
            echo "  --skip-reference            Skip reference download"
            echo "  -h, --help                  Show this help"
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Run main installation
main
