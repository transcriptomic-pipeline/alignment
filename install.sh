#!/bin/bash

set -euo pipefail

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Versions and URLs
STAR_VERSION="2.7.11b"
STAR_URL="https://github.com/alexdobin/STAR/releases/download/${STAR_VERSION}/STAR_${STAR_VERSION}.zip"
HISAT2_VERSION="2.2.1"
HISAT2_URL="https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download/hisat2-${HISAT2_VERSION}-Linux_x86_64.zip"

# Defaults
DEFAULT_INSTALL_DIR="${HOME}/softwares"
INSTALL_BASE_DIR=""
STAR_DIR=""
HISAT2_DIR=""
BIN_DIR=""

INSTALL_STAR=""
INSTALL_HISAT2=""
INSTALL_SAMTOOLS="yes"

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

prompt_install_directory() {
    echo ""
    echo "========================================"
    echo "  Installation Directory"
    echo "========================================"
    echo ""
    log_info "Choose installation directory"
    echo ""
    echo "  1) ${HOME}/softwares (recommended, matches QC module)"
    echo "  2) /opt/alignment-tools (requires sudo)"
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
    echo "  1) STAR only (high memory ~32GB)"
    echo "  2) HISAT2 only (low memory ~8GB)"
    echo "  3) Both STAR and HISAT2"
    echo ""
    read -p "Enter choice [1-3] (default: 3): " ALIGNER_CHOICE

    case "${ALIGNER_CHOICE:-3}" in
        1) INSTALL_STAR="yes"; INSTALL_HISAT2="no" ;;
        2) INSTALL_STAR="no"; INSTALL_HISAT2="yes" ;;
        3) INSTALL_STAR="yes"; INSTALL_HISAT2="yes" ;;
        *) INSTALL_STAR="yes"; INSTALL_HISAT2="yes" ;;
    esac
}

check_star() {
    [ -f "${STAR_DIR}/STAR" ] && "${STAR_DIR}/STAR" --version >/dev/null 2>&1
}
install_star() {
    log_info "Installing STAR ${STAR_VERSION}..."
    cd /tmp
    wget -q --show-progress "${STAR_URL}" -O "STAR.zip"
    unzip -q "STAR.zip"
    local STAR_BINARY=$(find . -name "STAR" -type f -executable 2>/dev/null | grep "Linux_x86_64" | head -n 1)
    if [ -z "$STAR_BINARY" ]; then
        log_error "STAR binary not found"
        return 1
    fi
    mkdir -p "${STAR_DIR}"
    cp "$STAR_BINARY" "${STAR_DIR}/"
    ln -sf "${STAR_DIR}/STAR" "${BIN_DIR}/STAR"
    rm -rf STAR* "STAR.zip"
    "${STAR_DIR}/STAR" --version >/dev/null 2>&1 && log_success "STAR installed successfully" || { log_error "STAR verification failed"; return 1; }
}

check_hisat2() {
    [ -f "${HISAT2_DIR}/hisat2" ] && "${HISAT2_DIR}/hisat2" --version >/dev/null 2>&1
}
install_hisat2() {
    log_info "Installing HISAT2 ${HISAT2_VERSION}..."
    cd /tmp
    wget -q --show-progress "${HISAT2_URL}" -O "hisat2.zip"
    unzip -q "hisat2.zip"
    # Find actual extracted directory, case-insensitive
    HISAT2_EXTRACTED=$(find . -maxdepth 1 -type d -iname "hisat2-*-linux_x86_64" | head -n 1)
    if [ -z "$HISAT2_EXTRACTED" ]; then
        log_error "HISAT2 extracted directory not found (fix: check extracted folder name)"
        return 1
    fi
    mkdir -p "${HISAT2_DIR}"
    cp -r "$HISAT2_EXTRACTED/"* "${HISAT2_DIR}/"
    ln -sf "${HISAT2_DIR}/hisat2" "${BIN_DIR}/hisat2"
    ln -sf "${HISAT2_DIR}/hisat2-build" "${BIN_DIR}/hisat2-build"
    ln -sf "${HISAT2_DIR}/hisat2-inspect" "${BIN_DIR}/hisat2-inspect"
    rm -rf hisat2* "hisat2.zip"
    "${HISAT2_DIR}/hisat2" --version >/dev/null 2>&1 && log_success "HISAT2 installed successfully" || { log_error "HISAT2 verification failed"; return 1; }
}

check_samtools() {
    command_exists samtools
}
install_samtools() {
    log_info "Installing samtools..."
    DISTRO=$(detect_distro)
    case "$DISTRO" in
        ubuntu|debian|linuxmint|mint)
            sudo apt-get update
            sudo apt-get install -y samtools
            ;;
        centos|rhel|redhat|fedora)
            sudo yum install -y samtools
            ;;
        *)
            log_error "Unsupported distribution for samtools"
            return 1 ;;
    esac
    command_exists samtools && log_success "samtools installed" || { log_error "samtools installation failed"; return 1; }
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
                sudo apt-get update
                sudo apt-get install -y wget unzip gzip ;;
            centos|rhel|redhat|fedora)
                sudo yum install -y wget unzip gzip ;;
        esac
    fi
    log_success "Dependencies ready"
}

update_path() {
    local SHELL_RC="${HOME}/.bashrc"
    if ! grep -q "# Alignment Module" "$SHELL_RC" 2>/dev/null; then
        echo "" >> "$SHELL_RC"
        echo "# Alignment Module - added by installer" >> "$SHELL_RC"
        echo "export PATH=\"${BIN_DIR}:\$PATH\"" >> "$SHELL_RC"
    fi
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
    log_success "Configuration saved"
}

main() {
    echo "========================================"
    echo "  Alignment Module Installation"
    echo "========================================"
    echo ""
    DISTRO=$(detect_distro)
    log_info "Detected OS: ${DISTRO}"
    [ -z "$INSTALL_BASE_DIR" ] && prompt_install_directory
    if [ -z "$INSTALL_STAR" ] && [ -z "$INSTALL_HISAT2" ]; then prompt_aligner_choice; fi
    check_dependencies
    echo ""
    log_info "Installing selected tools..."
    echo ""
    if [ "$INSTALL_STAR" = "yes" ]; then
        check_star || install_star || log_error "STAR installation failed"
    else log_info "Skipping STAR"; fi
    echo ""
    if [ "$INSTALL_HISAT2" = "yes" ]; then
        check_hisat2 || install_hisat2 || log_error "HISAT2 installation failed"
    else log_info "Skipping HISAT2"; fi
    echo ""
    if [ "$INSTALL_SAMTOOLS" = "yes" ]; then
        check_samtools || install_samtools || log_warning "samtools issues"
    fi
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
    log_info "Installation directory: ${INSTALL_BASE_DIR}"
    log_warning "Restart terminal or run: source ~/.bashrc"
    echo ""
}

# Parse argument flags
while [[ $# -gt 0 ]]; do
    case $1 in
        --install-dir)
            INSTALL_BASE_DIR="$2"
            STAR_DIR="${INSTALL_BASE_DIR}/STAR"
            HISAT2_DIR="${INSTALL_BASE_DIR}/HISAT2"
            BIN_DIR="${INSTALL_BASE_DIR}/bin"
            mkdir -p "${INSTALL_BASE_DIR}" "${BIN_DIR}"
            shift 2 ;;
        --aligner)
            case "$2" in
                star) INSTALL_STAR="yes"; INSTALL_HISAT2="no" ;;
                hisat2) INSTALL_STAR="no"; INSTALL_HISAT2="yes" ;;
                both) INSTALL_STAR="yes"; INSTALL_HISAT2="yes" ;;
                *) log_error "Invalid aligner: $2"; exit 1 ;;
            esac
            shift 2 ;;
        --skip-samtools) INSTALL_SAMTOOLS="no"; shift ;;
        -h|--help)
            cat << EOF
Usage: $0 [options]

Options:
  --install-dir <path>           Installation directory (default: ~/softwares)
  --aligner <star|hisat2|both>   Which aligner(s) to install
  --skip-samtools                Skip samtools installation
  -h, --help                     Show this help
EOF
            exit 0 ;;
        *)
            log_error "Unknown option: $1"; exit 1 ;;
    esac
done

main
