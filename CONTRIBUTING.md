# Contributing to ScType

Thank you for your interest in contributing to ScType! This document provides guidelines for contributing to the project.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [How to Contribute](#how-to-contribute)
- [Development Workflow](#development-workflow)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Documentation](#documentation)
- [Submitting Changes](#submitting-changes)

## Code of Conduct

This project adheres to a code of conduct. By participating, you are expected to uphold this code:
- Be respectful and inclusive
- Welcome newcomers
- Focus on what is best for the community
- Show empathy towards other community members

## Getting Started

### Prerequisites

**R Development:**
- R >= 4.0.0
- RStudio (recommended)
- Required packages: `dplyr`, `Seurat`, `HGNChelper`, `openxlsx`, `testthat`

**Python Development:**
- Python >= 3.8
- Required packages: `pandas`, `openpyxl` (see requirements.txt)

### Installation for Development

```bash
# Clone the repository
git clone https://github.com/IanevskiAleksandr/sc-type.git
cd sc-type

# Install R dependencies
R -e "install.packages(c('dplyr', 'Seurat', 'HGNChelper', 'openxlsx', 'testthat'))"

# Install Python dependencies
pip install -r requirements.txt
```

## How to Contribute

### Reporting Bugs

Before creating bug reports, please check existing issues. When creating a bug report, include:

- **Clear title and description**
- **Steps to reproduce**
- **Expected vs actual behavior**
- **R/Python version and package versions**
- **Example code or data** (if possible)

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. When suggesting an enhancement:

- **Use a clear, descriptive title**
- **Provide detailed description** of the proposed feature
- **Explain why this enhancement would be useful**
- **List any alternative solutions** you've considered

### Pull Requests

1. **Fork the repository** and create your branch from `master`
2. **Make your changes** following our coding standards
3. **Add tests** for new functionality
4. **Update documentation** as needed
5. **Ensure all tests pass**
6. **Submit a pull request**

## Development Workflow

### Branching Strategy

- `master` - Stable release branch
- `develop` - Development branch (if used)
- `feature/your-feature-name` - Feature branches
- `bugfix/issue-number` - Bug fix branches

### Making Changes

```bash
# Create a new branch
git checkout -b feature/my-new-feature

# Make your changes
# ... edit files ...

# Run tests
R -e "testthat::test_dir('tests/testthat')"

# Commit your changes
git commit -m "Add feature: brief description"

# Push to your fork
git push origin feature/my-new-feature
```

## Coding Standards

### R Code Style

Follow the [tidyverse style guide](https://style.tidyverse.org/):

```r
# Good
calculate_score <- function(data, threshold = 0.5) {
  result <- data %>%
    filter(value > threshold) %>%
    summarize(mean_value = mean(value))
  return(result)
}

# Bad
calculate.score=function(data,threshold=0.5){
result<-data%>%filter(value>threshold)%>%summarize(mean_value=mean(value))
return(result)}
```

**Key Points:**
- Use `<-` for assignment, not `=`
- Use `TRUE`/`FALSE`, not `T`/`F` or `!0`/`!1`
- Consistent indentation (2 spaces)
- Meaningful variable names
- Add roxygen2 documentation for all exported functions

### Python Code Style

Follow [PEP 8](https://pep8.org/):

```python
# Good
def calculate_score(data, threshold=0.5):
    """Calculate score from data above threshold.

    Args:
        data: Pandas DataFrame with numerical values
        threshold: Minimum value for inclusion (default: 0.5)

    Returns:
        float: Mean score
    """
    filtered = data[data['value'] > threshold]
    return filtered['value'].mean()

# Bad
def CalculateScore(data,threshold=.5):
    filtered=data[data['value']>threshold]
    return filtered['value'].mean()
```

**Key Points:**
- Use 4 spaces for indentation
- Maximum line length of 100 characters
- Add docstrings to all functions
- Use type hints where appropriate

### Naming Conventions

**R Functions:**
- Use snake_case: `sctype_score()`, `gene_sets_prepare()`
- Descriptive names: `auto_detect_tissue_type()` not `auto_detect()`

**Python Functions:**
- Use snake_case: `create_enhanced_database()`
- Descriptive names: `create_hierarchical_database()` not `create_db()`

**Variables:**
- Descriptive: `marker_sensitivity` not `ms`
- Boolean: `is_scaled`, `has_markers`

## Testing

### R Tests

We use `testthat` for unit testing:

```r
# tests/testthat/test-your-feature.R
context("Your feature")

test_that("function handles valid input", {
  result <- your_function(valid_input)
  expect_equal(result, expected_output)
})

test_that("function rejects invalid input", {
  expect_error(your_function(invalid_input), "error message")
})
```

**Running Tests:**
```r
testthat::test_dir("tests/testthat")
```

### Test Coverage

- Aim for >80% code coverage
- Test edge cases and error conditions
- Include integration tests for workflows

## Documentation

### Function Documentation (R)

Use roxygen2 format:

```r
#' Calculate ScType scores for cell type annotation
#'
#' Computes enrichment scores for each cell type in each cell based on
#' positive and negative marker gene expression.
#'
#' @param scRNAseqData Expression matrix (genes × cells)
#' @param scaled Boolean indicating if data is z-scaled (default: TRUE)
#' @param gs List of positive marker gene sets
#' @param gs2 List of negative marker gene sets (default: NULL)
#'
#' @return Matrix of cell types × cells with enrichment scores
#'
#' @examples
#' \dontrun{
#' scores <- sctype_score(expr_matrix, scaled = TRUE, gs = markers)
#' }
#'
#' @export
sctype_score <- function(scRNAseqData, scaled = TRUE, gs, gs2 = NULL) {
  # Implementation
}
```

### README Updates

Update README.md when adding:
- New features
- New dependencies
- Changes to installation process
- Breaking changes

## Submitting Changes

### Pull Request Process

1. **Update documentation** - All user-facing changes
2. **Add tests** - For new functionality
3. **Update CHANGELOG.md** - Describe your changes
4. **Run R CMD check** - Ensure package builds
5. **Describe your PR** - Clear title and description
6. **Link issues** - Reference related issues

### PR Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Breaking change
- [ ] Documentation update

## Testing
- [ ] Tests added/updated
- [ ] All tests pass
- [ ] R CMD check passes

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
- [ ] CHANGELOG.md updated
```

### Review Process

- Maintainers will review your PR
- Address feedback and make requested changes
- Once approved, a maintainer will merge

## Additional Resources

- **ScType Publication**: https://doi.org/10.1038/s41467-022-28803-w
- **Repository**: https://github.com/IanevskiAleksandr/sc-type
- **Web Portal**: http://sctype.app
- **Issue Tracker**: https://github.com/IanevskiAleksandr/sc-type/issues

## Questions?

- Open an issue for general questions
- Contact: aleksandr.ianevski@helsinki.fi

## License

By contributing, you agree that your contributions will be licensed under the GNU General Public License v3.0.

---

Thank you for contributing to ScType! Your efforts help make single-cell analysis more accessible to the research community.
