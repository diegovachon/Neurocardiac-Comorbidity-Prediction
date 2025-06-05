# Neurocardiac Project - Alzheimer's Disease Data Processing Pipeline

This repository contains a data processing pipeline for Alzheimer's Disease gene expression datasets from various GEO (Gene Expression Omnibus) series. The pipeline handles multiple datasets and platforms, performing standardized preprocessing steps to prepare the data for downstream analysis.

## Overview

The project processes several key Alzheimer's Disease datasets:
- GSE84422 (MSBB) with multiple platforms (GPL570, GPL96, GPL97)
- GSE125050 (ROSMAP)
- GSE67333

## Features

- Automated detection and processing of GEO series matrix files
- Support for multiple microarray platforms
- Standardized preprocessing steps including:
  - Expression data normalization
  - Missing value imputation
  - Metadata extraction and processing
  - Label generation for disease status
- Results saved in structured format for downstream analysis

## Project Structure

```
.
├── working_pipeline.py    # Main processing pipeline
└── results/              # Processed data output directory
```

## Requirements

- Python 3.x
- Required Python packages:
  - pandas
  - numpy
  - scikit-learn
  - scipy
  - gzip

## Usage

1. Place your GEO series matrix files in the appropriate data directory
2. Run the pipeline:

```python
from working_pipeline import run_working_pipeline

# Specify your data directory
data_directory = "path/to/your/data"

# Run the pipeline
run_working_pipeline(data_directory)
```

## Data Processing Steps

1. File Detection: Automatically detects and validates GEO series matrix files
2. Metadata Extraction: Parses sample information and characteristics
3. Expression Data Processing:
   - Handles missing values
   - Performs standardization
   - Processes platform-specific data formats
4. Label Generation: Creates standardized disease status labels
5. Results Storage: Saves processed data in structured format

## Output

The pipeline generates processed data files in the `results/` directory, including:
- Processed expression matrices
- Metadata information
- Disease status labels

## Notes

- The pipeline prioritizes datasets based on their quality and completeness
- Metadata-only datasets (GSE125050, GSE67333) are processed with lower priority
- Platform-specific processing is implemented for different microarray types

## Contributing

Feel free to submit issues and enhancement requests.
