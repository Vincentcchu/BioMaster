# BioMaster - Quick Start for Benchmarking

This is a quick start guide specifically for running BioMaster on the cell type annotation benchmarking task. For full documentation, see the [main BioMaster README](README.md).

## Quick Setup

1. **Install dependencies**:
   ```bash
   cd agents/biomaster
   pip install -r requirements.txt
   ```

2. **Configure API keys**:
   ```bash
   cp config.yaml.template config.yaml
   # Edit config.yaml with your API keys
   ```

## Running for Cell Type Annotation

```bash
cd agents/biomaster
python run.py config.yaml
```

## Input Data

- **Location**: `../../data/dataset_restricted.h5ad`
- **Format**: AnnData h5ad file with gene expression data

## Output

- **Location**: `../../outputs/biomaster/`
- **Contains**: Cell type predictions and analysis logs

## Configuration

Edit `config.yaml` to customize:
- `biomaster.id` - IDs of cells/samples to process
- `models.main` - Main LLM model
- `biomaster.executor` - Execution mode

See the [full README](README.md) for complete documentation.
