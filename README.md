# alignment
Alignment of transcriptomic reads with human reference genome using STAR and HISAT2 aligner

# RNA-seq Alignment Module v1.0.0

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/yourusername/alignment-module/releases/tag/v1.0.0)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

A production-ready RNA-seq alignment module supporting STAR and HISAT2 aligners with automatic installation, intelligent file detection, and seamless integration with upstream QC modules.

## Features

- **Dual Aligner Support**: STAR (splice-aware, high-accuracy) and HISAT2 (fast, low-memory)
- **Smart Python Handling**: Auto-detects system Python 3.6+ or creates minimal conda environment
- **Intelligent Pattern Detection**: Automatically detects `_R1_/_R2_` or `_1/_2` naming patterns
- **Auto-Index Building**: Builds genome indexes on-demand with progress tracking
- **Pipeline Integration**: Seamlessly works with QC module output
- **Production Ready**: Comprehensive error handling, logging, and user feedback

## Quick Start

Clone repository
git clone https://github.com/yourusername/alignment-module.git
cd alignment-module

Install (interactive mode)
bash install.sh

Run alignment (interactive aligner selection)
bash run_alignment.sh -i qc_output/trimmed/PE -o alignment_results -t 12

Or specify aligner
bash run_alignment.sh -i qc_output/trimmed/PE -o alignment_results --aligner star -t 20

## System Requirements

### Minimum
- **OS**: Linux (Ubuntu 18.04+, CentOS 7+, Debian 10+)
- **RAM**: 8 GB (HISAT2) or 32 GB (STAR)
- **Disk**: 50 GB for reference genome and indexes
- **CPU**: 4 cores minimum, 12+ recommended

### Dependencies
- Python 3.6+ (auto-detected or installed via conda)
- samtools
- wget, unzip, gunzip

## Installation

See [INSTALLATION.md](INSTALLATION.md) for detailed installation instructions.

Quick install with STAR
bash install.sh --aligner star --download-reference

Quick install with HISAT2
bash install.sh --aligner hisat2 --download-reference

Install both aligners
bash install.sh --aligner both --download-reference

## Usage

See [USAGE.md](USAGE.md) for complete usage examples.

### Basic Usage

Interactive mode - prompts for aligner selection
bash run_alignment.sh -i input_directory -o output_directory

Use STAR only
bash run_alignment.sh -i qc_output/trimmed/PE -o alignment_output --aligner star -t 20

Use HISAT2 only
bash run_alignment.sh -i qc_output/trimmed/PE -o alignment_output --aligner hisat2 -t 12

Use BOTH aligners for comparison
bash run_alignment.sh -i qc_output/trimmed/PE -o alignment_output --aligner both -t 20

### Input Format

Expected paired-end FASTQ naming patterns:
- `SampleID_R1_paired.fastq` / `SampleID_R2_paired.fastq` (QC module output)
- `SampleID_R1.fastq` / `SampleID_R2.fastq`
- `SampleID_1.fastq` / `SampleID_2.fastq`

### Output Structure

alignment_output/

├── SampleID1/

│ ├── SampleID1_Aligned.sortedByCoord.out.bam

│ ├── SampleID1_Aligned.sortedByCoord.out.bam.bai

│ └── SampleID1_Log.final.out

├── SampleID2/

│ └── ...

├── logs/

│ ├── SampleID1_star.log

│ └── ...

└── alignment_summary.txt

## Reference Genome

- **Species**: Homo sapiens
- **Assembly**: GRCh38 (Ensembl release 113)
- **FASTA**: `Homo_sapiens.GRCh38.dna.primary_assembly.fa`
- **GTF**: `Homo_sapiens.GRCh38.113.gtf`

## Integration with QC Module

This module is designed to work seamlessly with the [QC Module](https://github.com/yourusername/qc-module):

Complete workflow
bash qc/run_qc.sh -i raw_fastq/ -o qc_output
bash alignment/run_alignment.sh -i qc_output/trimmed/PE -o alignment_output --aligner star

## Performance Benchmarks

| Aligner | Memory Usage | Time (6 samples) | Splice Detection |
|---------|--------------|------------------|------------------|
| STAR    | ~32 GB       | ~45 min          | Excellent        |
| HISAT2  | ~8 GB        | ~30 min          | Good             |

*Tested on Intel Xeon 12-core, 64GB RAM*

## Troubleshooting

### Common Issues

**1. "Configuration not found"**

Solution: Run installer first
bash install.sh

**2. "STAR index not found"**

Solution: Script will auto-build on first run (30-60 min)
Or manually build:
bash install.sh --download-reference

**3. "Python not found" (HISAT2)**

Solution: Install Python 3.6+
sudo apt-get install python3 # Ubuntu/Debian

See [TROUBLESHOOTING.md](TROUBLESHOOTING.md) for more solutions.

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file.

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Support

- **Issues**: [GitHub Issues](https://github.com/transcriptomic-pipeline/alignment/issues)
- **Discussions**: [GitHub Discussions](https://github.com/transcriptomic-pipeline/alignment/discussions)
- **Email**: babul@ibioinformatics.com

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history.

## Acknowledgments

- **STAR**: Dobin et al., Bioinformatics 2013
- **HISAT2**: Kim et al., Nature Methods 2019
- **Ensembl**: Human reference genome GRCh38

---

**Developed by** [Babul Pradhan and Dr. Jyoti Sharma, Institute of Bioinformatics, Bangalore] | **Version** 1.0.0 | **Release Date** November 14, 2025

## Citation

If you use this module in your research, please cite:

@software{alignment_module_2025,
author = {Babul Pradhan, Dr. Jyoti Sharma},
title = {RNA-seq Alignment Module},
year = {2025},
version = {1.0.0},
url = {https://github.com/transcriptomic-pipeline/alignment}
}
