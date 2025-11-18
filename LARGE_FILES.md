# Large Files Management

This document explains the large binary files in this repository and how they are managed.

## Large Files in Repository

| File | Size | Purpose | LFS Tracked |
|------|------|---------|-------------|
| `exampleData.RDS` | 6.2 MB | Example scRNA-seq dataset (PBMC 3k) | ✅ Yes |
| `index2.html` | 2.0 MB | Alternative web interface | ✅ Yes |
| `index.html` | 198 KB | Main web interface | ✅ Yes |
| `*.xlsx` | Various | Cell marker databases | ✅ Yes |
| `*.png` | Various | Documentation images | ✅ Yes |

## Git LFS Configuration

This repository uses **Git Large File Storage (LFS)** to efficiently handle large binary files.

### What is Git LFS?

Git LFS replaces large files with text pointers inside Git, while storing the file contents on a remote server. This keeps your repository clone fast and efficient.

### Installation

**One-time setup:**
```bash
# Install Git LFS (if not already installed)
# macOS
brew install git-lfs

# Ubuntu/Debian
sudo apt-get install git-lfs

# Windows (with Git for Windows)
# Git LFS is included by default

# Initialize Git LFS
git lfs install
```

### Cloning This Repository

When cloning this repository, Git LFS will automatically download the large files:

```bash
# Clone with LFS files
git clone https://github.com/IanevskiAleksandr/sc-type.git

# LFS files are downloaded automatically
# If you want to skip LFS files:
GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/IanevskiAleksandr/sc-type.git

# Download LFS files later if needed:
git lfs pull
```

### Tracked File Types

The following file patterns are tracked with Git LFS (see `.gitattributes`):

- **Data files**: `*.RDS`, `*.rds`
- **Large HTML**: `index2.html`
- **Excel files**: `*.xlsx`
- **Images**: `*.png`, `*.jpg`, `*.jpeg`, `*.PNG`
- **Archives**: `*.tar.gz`, `*.zip`

### Working with LFS Files

```bash
# Check LFS file status
git lfs ls-files

# See which files will be tracked by LFS
git lfs track

# Manually track additional file patterns
git lfs track "*.newextension"

# Push LFS files to remote
git lfs push origin main

# Pull LFS files from remote
git lfs pull
```

## File Details

### exampleData.RDS

**Purpose**: Example scRNA-seq dataset for testing and demonstrations
- **Source**: PBMC 3k dataset from 10X Genomics
- **Format**: R Data Serialization (RDS)
- **Contents**: Pre-processed Seurat object with normalized and scaled data
- **Used in**: README examples, testing, validation

**Alternative Download**:
If you don't want to download via Git LFS, you can get it directly:
```r
# Download from 10X Genomics
download.file(
  "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
  "pbmc3k_filtered_gene_bc_matrices.tar.gz"
)
```

### index2.html

**Purpose**: Alternative web interface for ScType portal
- **Format**: HTML with embedded JavaScript and CSS
- **Size**: Large due to embedded visualization libraries
- **Usage**: Open directly in web browser for interactive cell type annotation

### Database Files (*.xlsx)

**Purpose**: Cell marker gene databases
- **ScTypeDB_full.xlsx** (27 KB): Complete marker database
- **ScTypeDB_short.xlsx** (20 KB): Abbreviated version
- **ScTypeDB_enhanced.xlsx** (12 KB): Enhanced with 122 cell types
- **ScTypeDB_hierarchical.xlsx** (12 KB): Hierarchical annotation markers

## Storage Considerations

### For Contributors

When adding large files to the repository:

1. **Check if LFS tracking is needed** (files > 100 KB)
2. **Use appropriate file pattern** in `.gitattributes`
3. **Compress when possible** (especially text-based files)
4. **Consider external hosting** for very large files (> 10 MB)

### For Users

**Minimal clone (without large files)**:
```bash
# Clone without LFS files (fastest)
GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/IanevskiAleksandr/sc-type.git

# Install required packages and use remote databases
# Functions will download databases from GitHub as needed
```

**Full clone (with all files)**:
```bash
# Standard clone (includes all LFS files)
git clone https://github.com/IanevskiAleksandr/sc-type.git
```

## Troubleshooting

### LFS files not downloading

```bash
# Check Git LFS is installed
git lfs version

# Manually pull LFS files
git lfs pull

# Check LFS file status
git lfs ls-files
```

### LFS files showing as text pointers

This means Git LFS isn't pulling the actual files. Fix:
```bash
# Install and configure Git LFS
git lfs install

# Pull LFS files
git lfs pull
```

### Large repository clone times

If cloning is slow:
```bash
# Clone without LFS files first
GIT_LFS_SKIP_SMUDGE=1 git clone <repo-url>

# Download only needed LFS files
cd sc-type
git lfs pull --include="exampleData.RDS"  # Only download specific files
```

## Best Practices

### Adding New Large Files

1. **Check size**: Files > 100 KB should use LFS
2. **Add to .gitattributes**: Update pattern if needed
3. **Test tracking**:
   ```bash
   git add yourfile.ext
   git lfs ls-files  # Should show your file
   ```
4. **Commit normally**: Git LFS handles the rest

### Removing Large Files

If you need to remove a large file from history:
```bash
# Use BFG Repo-Cleaner or git-filter-repo
# WARNING: This rewrites history!

# Option 1: BFG (recommended)
bfg --delete-files filename.ext

# Option 2: git-filter-repo
git filter-repo --path filename.ext --invert-paths
```

## External Hosting Alternatives

For very large files (> 50 MB), consider:

1. **Zenodo** - Academic research data hosting (free, DOI)
2. **figshare** - Scientific data repository (free, DOI)
3. **OSF** - Open Science Framework (free)
4. **Amazon S3** - Cloud storage (paid)
5. **Release Assets** - GitHub releases (up to 2GB per file)

## Resources

- **Git LFS Documentation**: https://git-lfs.github.com/
- **Git LFS Tutorial**: https://www.atlassian.com/git/tutorials/git-lfs
- **GitHub LFS**: https://docs.github.com/en/repositories/working-with-files/managing-large-files

## Questions?

- **LFS Issues**: Check Git LFS documentation or repository issues
- **File Access**: Contact aleksandr.ianevski@helsinki.fi
- **Alternative Downloads**: See individual file sections above

---

**Last Updated**: 2025-11-18
