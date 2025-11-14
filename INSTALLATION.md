# Installation Guide

This document provides detailed installation instructions for the RNA-seq Alignment Module v1.0.0.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Quick Installation](#quick-installation)
- [Detailed Installation Steps](#detailed-installation-steps)
- [Installation Options](#installation-options)
- [Verifying Installation](#verifying-installation)
- [Post-Installation](#post-installation)
- [Troubleshooting](#troubleshooting)

## Prerequisites

### System Requirements

**Minimum Requirements:**
- **Operating System**: Linux (Ubuntu 18.04+, CentOS 7+, Debian 10+, Red Hat 7+)
- **RAM**: 
  - 8 GB minimum (HISAT2 only)
  - 32 GB recommended (STAR)
  - 40 GB for running both aligners
- **Storage**: 
  - 50 GB for reference genome and indexes
  - Additional space for your sequencing data and results
- **CPU**: 4 cores minimum, 12+ cores recommended

**Software Dependencies:**
- Python 3.6+ (will be auto-detected or installed)
- wget
- unzip
- gunzip
- samtools (will be installed if missing)
- sudo access (for installing system packages)

### Checking Prerequisites

Check Python version
python3 --version # Should be 3.6 or higher

Check available RAM
free -h

Check available disk space
df -h ~

Check CPU cores
nproc

## Quick Installation

### Option 1: Interactive Installation (Recommended for First-Time Users)

Clone the repository
git clone https://github.com/transcriptomic-pipeline/alignment.git
cd alignment

Make scripts executable
chmod +x install.sh run_alignment.sh

Run the alignment module directly (it will perform "a." and "b." below, automatically)
bash run_alignment.sh

OR

a. Run interactive installer
bash install.sh

b. Run interactive installer
bash run_alignment.sh

The installer will:
1. Detect system Python or install Miniconda
2. Prompt for installation directory
3. Ask which aligner(s) to install (STAR, HISAT2, or both)
4. Offer to download reference genome

### Option 2: Automated Installation

Install with STAR only
bash install.sh --aligner star --download-reference --install-dir ~/softwares

Install with HISAT2 only
bash install.sh --aligner hisat2 --download-reference --install-dir ~/softwares

Install both aligners
bash install.sh --aligner both --download-reference --install-dir ~/softwares

## Detailed Installation Steps

### Step 1: Download the Module

Via Git
git clone https://github.com/transcriptomic-pipeline/alignment.git
cd alignment

Or download ZIP
wget https://github.com/transcriptomic-pipeline/alignment/archive/alignment_v1.0.0.zip
unzip alignment_v1.0.0.zip
cd alignment

### Step 2: System Dependencies

The installer checks for required system packages. If missing, install manually:

**Ubuntu/Debian:**

sudo apt-get update
sudo apt-get install -y wget unzip gzip samtools python3

**CentOS/Red Hat:**

sudo yum update
sudo yum install -y wget unzip gzip samtools python3

### Step 3: Run Installer

**Interactive Mode:**

bash install.sh

**Command-Line Mode:**

bash install.sh
--aligner both
--install-dir ~/softwares
--reference-dir ~/references/GRCh38_ensembl113
--download-reference

### Step 4: Installation Process

The installer will:

1. **Detect Python** (3.6+)
   - Uses system Python if available
   - Creates minimal conda environment if needed

2. **Install Aligners**
   - STAR 2.7.11b (~150 MB)
   - HISAT2 2.2.1 (~10 MB)

3. **Download Reference Genome** (optional)
   - Human GRCh38 FASTA (~800 MB compressed, ~3.1 GB uncompressed)
   - Ensembl 113 GTF (~50 MB compressed, ~1.5 GB uncompressed)

4. **Setup Environment**
   - Adds tools to PATH in ~/.bashrc
   - Creates configuration files

**Expected Installation Time:**
- Aligner installation: 2-5 minutes
- Reference download: 5-15 minutes (depending on connection)
- Total: ~10-20 minutes

## Installation Options

### Command-Line Arguments

bash install.sh [OPTIONS]

Options:
--install-dir <path>   Installation directory (default: ~/softwares)

--reference-dir <path>   Reference genome directory (default: ~/references/GRCh38_ensembl113)

--aligner <star|hisat2|both>   Which aligner(s) to install

--download-reference   Download reference genome

--skip-reference   Skip reference download

-h, --help   Show help message

### Example Installations

**1. Minimal STAR Installation (no reference):**

bash install.sh --aligner star --skip-reference

**2. Complete HISAT2 Installation:**

bash install.sh --aligner hisat2 --download-reference

**3. Custom Directories:**

bash install.sh
--aligner both
--install-dir /softwares/alignment-tools
--reference-dir /data/references/GRCh38
--download-reference

**4. System-Wide Installation (requires sudo):**

sudo bash install.sh
--aligner both
--install-dir /opt/alignment-tools
--download-reference

## Verifying Installation

### Check Installed Tools

Source updated environment
source ~/.bashrc

Check STAR
STAR --version

Check HISAT2
hisat2 --version

### Verify Directory Structure

Installation directory structure
~/softwares/

├── bin/

│ ├── STAR

│ ├── hisat2

│ ├── hisat2-build

│ └── python -> /usr/bin/python3 # symlink if needed

├── STAR/

│ └── STAR

└── HISAT2/

├── hisat2

├── hisat2-build

└── ...

Reference directory structure
~/references/GRCh38_ensembl113/

├── Homo_sapiens.GRCh38.dna.primary_assembly.fa

├── Homo_sapiens.GRCh38.113.gtf

├── STAR_index/

│ └── (will be built on first run)

└── HISAT2_index/

└── (will be built on first run)

## Post-Installation

### Update PATH

If the installer didn't automatically update your PATH:

Add to ~/.bashrc
echo 'export PATH="$HOME/softwares/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

### Test Run (Optional)

Create test data directory
mkdir -p test_data

Download small test dataset (optional)
... add test data instructions ...
Run alignment test
bash run_alignment.sh -i test_data -o test_output --aligner star -t 4

## Troubleshooting

### Installation Fails to Detect Python

**Problem**: "No suitable system Python 3 found"

**Solution 1**: Install Python 3.6+

sudo apt-get install python3 python3-dev # Ubuntu/Debian
sudo yum install python3 python3-devel # CentOS/RHEL

**Solution 2**: Let installer create conda environment
- Re-run installer, it will automatically install Miniconda

### Permission Denied Errors

**Problem**: Cannot write to installation directory

**Solution 1**: Use home directory

bash install.sh --install-dir ~/softwares

**Solution 2**: Use sudo for system-wide installation

sudo bash install.sh --install-dir /opt/alignment-tools

### Download Failures

**Problem**: Reference genome download fails

**Solution**: Manual download

cd ~/references/GRCh38
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz

Decompress
gunzip *.gz

Re-run installer
cd ~/alignment-module
bash install.sh --skip-reference

### Insufficient Disk Space

**Problem**: Not enough space for reference genome

**Solution**: Use external/mounted drive

bash install.sh
--reference-dir /mnt/large_drive/references/GRCh38
--download-reference

### HISAT2 Python Compatibility Issues

**Problem**: HISAT2 scripts fail with Python errors

**Solution**: The installer automatically patches HISAT2 scripts, but if issues persist:

Check Python version
python3 --version # Must be 3.6+

Reinstall HISAT2
bash install.sh --aligner hisat2

## Uninstallation

To remove the alignment module:

Remove installation directory
rm -rf ~/softwares/STAR ~/softwares/HISAT2 ~/softwares/bin/STAR ~/softwares/bin/hisat2*

Remove references (optional)
rm -rf ~/references/GRCh38_ensembl113

Remove conda environment (if created)
conda env remove -n hisat2_env

Remove PATH entry from ~/.bashrc
Edit ~/.bashrc and remove the "# Alignment Module" section

## Getting Help

If you encounter issues not covered here:

1. Check [TROUBLESHOOTING.md](TROUBLESHOOTING.md)
2. Search [GitHub Issues](https://github.com/transcriptomic-pipeline/alignment/issues)
3. Open a new issue with:
   - Error message
   - Installation command used
   - System information (OS, RAM, Python version)
   - Installation log

## Next Steps

After successful installation:

1. Read [USAGE.md](USAGE.md) for usage examples
2. Test with sample data
3. Integrate with upstream QC module (optional)

---

**Installation Support**: Issues? Email your query to email@example.com or open a [GitHub Issue](https://github.com/transcriptomic-pipeline/alignment/issues)

