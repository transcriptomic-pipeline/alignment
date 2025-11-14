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

# Reference type
USE_CUSTOM_REFERENCE="no"
CUSTOM_FASTA=""
CUSTOM_GTF=""

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
    
    # First, try python3 command
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
    
    # Second, try python command (must be Python 3.6+)
    if command_exists python; then
        local PYTHON_VERSION_OUTPUT=$(python --version 2>&1)
        if [[ "$PYTHON_VERSION_OUTPUT" =~ Python\ ([0-9]+)\.([0-9]+) ]]; then
            local MAJOR="${BASH_REMATCH[1]}"
            local MINOR="${BASH_REMATCH[2]}"
            
            if [ "$MAJOR" -eq 3 ] && [ "$MINOR" -ge 6 ]; then
                PYTHON_BIN=$(command -v python)
                USE_SYSTEM_PYTHON="yes"
                log_success "System Python 3 found: $PYTHON_VERSION_OUTPUT"
                log_info "Using system Python: $PYTHON_BIN (command: python)"
                log_warning "HISAT2 will use 'python' command (not patching to python3)"
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

prompt_reference_selection() {
    echo ""
    echo "========================================"
    echo "  Reference Genome Selection"
    echo "========================================"
    echo ""
    
    # Check if there's a previous reference config
    local HAS_PREVIOUS="no"
    local PREV_TYPE=""
    local PREV_FASTA=""
    local PREV_GTF=""
    
    local REF_CONFIG="${SCRIPT_DIR}/config/reference_paths.conf"
    if [ -f "$REF_CONFIG" ]; then
        source "$REF_CONFIG"
        if [ -n "$REFERENCE_FASTA" ] && [ -f "$REFERENCE_FASTA" ] && [ -n "$REFERENCE_GTF" ] && [ -f "$REFERENCE_GTF" ]; then
            HAS_PREVIOUS="yes"
            PREV_TYPE="$USE_CUSTOM_REFERENCE"
            PREV_FASTA="$REFERENCE_FASTA"
            PREV_GTF="$REFERENCE_GTF"
        fi
    fi
    
    log_info "Select reference genome type"
    echo ""
    echo "  1) Default: Human GRCh38 (Ensembl 113) [recommended]"
    echo "  2) Custom: Provide your own genome FASTA and GTF"
    
    if [ "$HAS_PREVIOUS" = "yes" ]; then
        if [ "$PREV_TYPE" = "yes" ]; then
            echo "  3) Previous custom: $(basename $PREV_FASTA)"
        else
            echo "  3) Previous: GRCh38 (Ensembl 113)"
        fi
    fi
    
    echo ""
    
    # Auto-timeout prompt
    local timeout=30
    echo -ne "Enter choice [1-"
    [ "$HAS_PREVIOUS" = "yes" ] && echo -ne "3" || echo -ne "2"
    echo -ne "] (default: 1, auto-select in ${timeout}s): "
    
    read -t $timeout -n 1 -r REPLY || true
    echo
    
    # Default to 1 if no response
    if [ -z "$REPLY" ]; then
        REPLY="1"
        log_info "No response in ${timeout}s, using default (GRCh38)"
    fi
    
    case "${REPLY}" in
        1)
            USE_CUSTOM_REFERENCE="no"
            prompt_reference_directory
            ;;
        2)
            USE_CUSTOM_REFERENCE="yes"
            prompt_custom_reference
            ;;
        3)
            if [ "$HAS_PREVIOUS" = "yes" ]; then
                # Use previous reference
                USE_CUSTOM_REFERENCE="$PREV_TYPE"
                CUSTOM_FASTA="$PREV_FASTA"
                CUSTOM_GTF="$PREV_GTF"
                REFERENCE_DIR="$(dirname "$PREV_FASTA")"
                
                log_success "Using previous reference"
                log_info "FASTA: $CUSTOM_FASTA"
                log_info "GTF: $CUSTOM_GTF"
                log_info "Reference directory: $REFERENCE_DIR"
            else
                log_error "Invalid choice"
                exit 1
            fi
            ;;
        *)
            log_error "Invalid choice. Please select a valid option."
            exit 1
            ;;
    esac
}

prompt_reference_directory() {
    echo ""
    log_info "Default reference genome directory"
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

prompt_custom_reference() {
    echo ""
    log_info "Custom reference genome configuration"
    echo ""
    log_warning "You need to provide:"
    echo "  1. Genome FASTA file (.fa or .fasta)"
    echo "  2. Gene annotation GTF file (.gtf)"
    echo ""
    
    # Prompt for FASTA
    read -p "Enter path to genome FASTA file: " CUSTOM_FASTA
    CUSTOM_FASTA="${CUSTOM_FASTA/#\~/$HOME}"
    
    # Prompt for GTF
    read -p "Enter path to GTF annotation file: " CUSTOM_GTF
    CUSTOM_GTF="${CUSTOM_GTF/#\~/$HOME}"
    
    # Validate files
    if [ ! -f "$CUSTOM_FASTA" ]; then
        log_error "FASTA file not found: $CUSTOM_FASTA"
        echo ""
        log_info "Please place your genome FASTA file at the specified path and re-run the installer."
        log_info "Example: cp /path/to/your/genome.fa $CUSTOM_FASTA"
        exit 1
    fi
    
    if [ ! -f "$CUSTOM_GTF" ]; then
        log_error "GTF file not found: $CUSTOM_GTF"
        echo ""
        log_info "Please place your GTF annotation file at the specified path and re-run the installer."
        log_info "Example: cp /path/to/your/genes.gtf $CUSTOM_GTF"
        exit 1
    fi
    
    # Set reference directory to parent directory of FASTA
    REFERENCE_DIR="$(dirname "$CUSTOM_FASTA")"
    
    log_success "Custom reference validated:"
    log_info "FASTA: $CUSTOM_FASTA"
    log_info "GTF: $CUSTOM_GTF"
    log_info "Reference directory: $REFERENCE_DIR"
}

# Validate custom reference from CLI arguments
validate_custom_reference_cli() {
    if [ "$USE_CUSTOM_REFERENCE" = "yes" ]; then
        # Validate FASTA
        if [ -z "$CUSTOM_FASTA" ]; then
            log_error "Custom reference selected but no FASTA file specified"
            log_info "Use: --custom-reference /path/to/genome.fa"
            exit 1
        fi
        
        if [ ! -f "$CUSTOM_FASTA" ]; then
            log_error "FASTA file not found: $CUSTOM_FASTA"
            echo ""
            log_info "Please place your genome FASTA file at the specified path and re-run."
            log_info "Example: cp /path/to/your/genome.fa $CUSTOM_FASTA"
            exit 1
        fi
        
        # Validate GTF
        if [ -z "$CUSTOM_GTF" ]; then
            log_error "Custom reference selected but no GTF file specified"
            log_info "Use: --custom-gtf /path/to/genes.gtf"
            exit 1
        fi
        
        if [ ! -f "$CUSTOM_GTF" ]; then
            log_error "GTF file not found: $CUSTOM_GTF"
            echo ""
            log_info "Please place your GTF annotation file at the specified path and re-run."
            log_info "Example: cp /path/to/your/genes.gtf $CUSTOM_GTF"
            exit 1
        fi
        
        # Set reference directory
        REFERENCE_DIR="$(dirname "$CUSTOM_FASTA")"
        
        log_success "Custom reference validated (CLI):"
        log_info "FASTA: $CUSTOM_FASTA"
        log_info "GTF: $CUSTOM_GTF"
    fi
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

# Patch HISAT2 scripts to use python3 instead of python (ONLY if needed)
patch_hisat2_python() {
    # Only patch if PYTHON_BIN is python3 (not python)
    if [[ "$PYTHON_BIN" == *"python3"* ]]; then
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
    else
        log_info "HISAT2 scripts will use 'python' command (no patching needed)"
    fi
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
    if [ "$USE_CUSTOM_REFERENCE" = "yes" ]; then
        log_info "Using custom reference genome (skipping download)"
        return 0
    fi
    
    log_info "Downloading default reference genome (GRCh38 Ensembl 113)..."
    
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
    
    # Determine reference files based on custom vs default
    if [ "$USE_CUSTOM_REFERENCE" = "yes" ]; then
        local REF_FASTA="$CUSTOM_FASTA"
        local REF_GTF="$CUSTOM_GTF"
        local GENOME_BASENAME=$(basename "${REF_FASTA}" | sed 's/\.[^.]*$//')
        local STAR_IDX="${REFERENCE_DIR}/STAR_index_${GENOME_BASENAME}"
        local HISAT2_IDX="${REFERENCE_DIR}/HISAT2_index_${GENOME_BASENAME}"
        local HISAT2_PREFIX="${HISAT2_IDX}/genome"
    else
        local REF_FASTA="${REFERENCE_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        local REF_GTF="${REFERENCE_DIR}/Homo_sapiens.GRCh38.113.gtf"
        local STAR_IDX="${REFERENCE_DIR}/STAR_index"
        local HISAT2_IDX="${REFERENCE_DIR}/HISAT2_index"
        local HISAT2_PREFIX="${HISAT2_IDX}/genome"
    fi
    
    cat > "${SCRIPT_DIR}/config/reference_paths.conf" << EOF
# Reference Genome Configuration
# Generated: $(date)

USE_CUSTOM_REFERENCE="${USE_CUSTOM_REFERENCE}"
REFERENCE_DIR="${REFERENCE_DIR}"
REFERENCE_FASTA="${REF_FASTA}"
REFERENCE_GTF="${REF_GTF}"
STAR_INDEX_DIR="${STAR_IDX}"
HISAT2_INDEX_DIR="${HISAT2_IDX}"
HISAT2_INDEX_PREFIX="${HISAT2_PREFIX}"
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
    
    # Validate custom reference if provided via CLI
    validate_custom_reference_cli
    
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
    
    # Download reference genome OR prompt for selection
    if [ -z "$DOWNLOAD_REFERENCE" ]; then
        echo ""
        if [ "$USE_CUSTOM_REFERENCE" = "no" ]; then
            read -p "Download default reference genome (GRCh38)? [Y/n] " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Nn]$ ]]; then
                DOWNLOAD_REFERENCE="yes"
                # Only prompt for directory if not already set
                if [ -z "$REFERENCE_DIR" ]; then
                    prompt_reference_selection
                fi
            else
                DOWNLOAD_REFERENCE="no"
            fi
        else
            DOWNLOAD_REFERENCE="no"
        fi
    else
        # CLI argument provided, prompt for reference selection if needed
        if [ "$DOWNLOAD_REFERENCE" = "yes" ] && [ -z "$REFERENCE_DIR" ]; then
            prompt_reference_selection
        fi
    fi
    
    if [ "$DOWNLOAD_REFERENCE" = "yes" ]; then
        download_reference
    else
        if [ "$USE_CUSTOM_REFERENCE" = "no" ] && [ -z "$REFERENCE_DIR" ]; then
            REFERENCE_DIR="${DEFAULT_REFERENCE_DIR}"
        fi
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
    if [ "$USE_CUSTOM_REFERENCE" = "yes" ]; then
        log_info "Custom reference: $(basename $CUSTOM_FASTA)"
    elif [ "$DOWNLOAD_REFERENCE" = "yes" ]; then
        log_info "Reference directory: ${REFERENCE_DIR}"
    fi
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
        --custom-reference)
            USE_CUSTOM_REFERENCE="yes"
            CUSTOM_FASTA="$2"
            CUSTOM_FASTA="${CUSTOM_FASTA/#\~/$HOME}"
            shift 2
            ;;
        --custom-gtf)
            CUSTOM_GTF="$2"
            CUSTOM_GTF="${CUSTOM_GTF/#\~/$HOME}"
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
            echo "  --install-dir <path>         Installation directory (default: ~/softwares)"
            echo "  --reference-dir <path>       Reference genome directory"
            echo "  --aligner <star|hisat2|both> Which aligner(s) to install"
            echo "  --download-reference         Download default reference (GRCh38)"
            echo "  --skip-reference             Skip reference download"
            echo "  --custom-reference <path>    Use custom genome FASTA"
            echo "  --custom-gtf <path>          Use custom GTF annotation"
            echo "  -h, --help                   Show this help"
            echo ""
            echo "Examples:"
            echo "  # Install with default reference"
            echo "  $0 --aligner both --download-reference"
            echo ""
            echo "  # Install with custom reference"
            echo "  $0 --aligner star \\"
            echo "    --custom-reference /path/to/genome.fa \\"
            echo "    --custom-gtf /path/to/genes.gtf"
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
