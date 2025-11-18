# Changelog

All notable changes to the ScType package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2025-11-18

### Added
- **Package Structure**: Added DESCRIPTION and NAMESPACE files for proper R package structure
- **Testing**: Created testthat-based unit test suite for core functions
  - Tests for sctype_score() function
  - Tests for gene_sets_prepare() function
  - Tests for wrapper functions
- **Documentation**: Added CHANGELOG.md and CONTRIBUTING.md
- **Large Files Documentation**: Created LARGE_FILES.md explaining binary file management
- **Python Dependencies**: Created requirements.txt for database creation scripts
- **Command-Line Interface**: Added argparse support to Python scripts for flexible output paths
- **Git LFS Support**: Added .gitattributes for efficient large file tracking
  - Tracks RDS files, Excel databases, images, and large HTML files
  - Improves repository clone performance
- **Comprehensive .gitignore**: Enhanced from 27 to 168 lines
  - Python: virtual environments, caching, test coverage
  - R: history files, RStudio projects, vignettes, packages
  - IDEs: VSCode, PyCharm, Sublime, Vim, Emacs
  - OS: macOS, Windows, Linux specific files
- **Installation Guide**: Created INSTALLATION.md (300+ lines)
  - Complete installation instructions for GitHub, local, and package archive methods
  - Dependency documentation
  - Quick start guides for Seurat and SingleCellExperiment
  - Advanced features documentation
  - Troubleshooting section
  - Migration guide from old source-based loading

### Changed
- **Error Handling**: Improved error messages throughout codebase
  - Converted warnings to errors for critical failures (empty matrices, invalid inputs)
  - Added descriptive error messages with actionable guidance
  - Better handling of edge cases (empty gene sets, missing markers)
- **Code Style**: Standardized coding style
  - Replaced `!0` with `TRUE` and `!1` with `FALSE` throughout
  - Consistent use of `<-` for assignment
  - Improved variable naming consistency
- **Python Scripts**: Fixed hard-coded paths
  - create_enhanced_db.py now uses relative paths and command-line arguments
  - create_hierarchical_db.py now uses relative paths and command-line arguments
  - Better portability across systems
- **HTTP Security Fix**: Removed all remote source() calls (8 R files modified)
  - R/sctype_wrapper.R - Modified sctype_source() function
  - R/sctype_hierarchical.R - Removed 3 source() calls
  - R/sctype_uncertainty.R - Removed 2 source() calls
  - R/sctype_uncertainty_sce.R - Removed 2 source() calls
  - R/sctype_pathway_enrichment.R - Removed 1 source() call
  - R/sctype_wrapper_sce.R - Removed 4 source() calls (including sctype_source_sce)
  - R/sctype_hierarchical_sce.R - Removed 2 source() calls
  - Functions now loaded from package namespace instead of HTTP
  - Eliminates remote code execution vulnerability

### Fixed
- **Documentation**: Fixed typo in sctype_wrapper.R example
  - Corrected missing quote in function example
  - Improved example formatting with \dontrun{}
- **Safety**: Enhanced input validation in sctype_score()
  - Check for matrix type before processing
  - Validate matrix dimensions
  - Ensure marker genes are found in dataset

### Removed
- **Duplicate Files**: Removed duplicate database files
  - Deleted "ScTypeDB_full (2).xlsx"
  - Deleted "ScTypeDB_full_original_backup.xlsx"
  - Kept ScTypeDB_full.xlsx as canonical version

### Security
- **FIXED**: Eliminated HTTP source loading security vulnerability
  - Removed ALL 17 instances of `source("https://...")` calls
  - Functions now available via package namespace (DESCRIPTION/NAMESPACE)
  - No remote code execution risk
  - Users install package once: `devtools::install_github("IanevskiAleksandr/sc-type")`
  - Offline functionality after installation
  - Version control and reproducibility guaranteed
- Created INSTALLATION.md (300+ lines) with migration guide from old source-based loading

## [1.0.0] - 2022

### Added
- Initial release published in Nature Communications
- Core sctype_score() algorithm
- Support for Seurat objects
- Database of marker genes for multiple tissues
- Web portal at http://sctype.app

### Features from Publication
- Automated cell type annotation
- Tissue type auto-detection
- Positive and negative marker support
- Fast computation optimized for large datasets

---

## Unreleased

### Planned
- ~~Address HTTP source loading security issue~~ ✅ **DONE in v2.0.0**
- Add parallel processing support for large datasets
- Expand test coverage to >80%
- Performance benchmarks
- Continuous integration (GitHub Actions)
- CRAN/Bioconductor submission

---

## Notes

For detailed information about the ScType method, see:
- **Publication**: https://doi.org/10.1038/s41467-022-28803-w
- **Repository**: https://github.com/IanevskiAleksandr/sc-type
- **Web Portal**: http://sctype.app
