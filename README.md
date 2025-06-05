Neurocardiac Comorbidity Prediction
Overview
This project aims to develop an early detection system for the risk of simultaneously developing cardiovascular and neurodegenerative diseases (e.g., stroke + Alzheimer's) through analysis of blood transcriptomic profiles.
Objectives

Primary Goal: Early detection of neurocardiac comorbidity risk
Approach: Transcriptomic signature analysis to predict combined risks
Impact: Improve preventive patient care for at-risk individuals

Technical Architecture
Processing Pipeline

RNA-seq Preprocessing: Normalization, quality filtering, batch correction
Clustering: Identification of patient groups with similar profiles
Multi-label Classification: Simultaneous prediction of cardio + neuro risks
Multimodal Learning: Integration with brain imaging (when available)

Technologies Used

Python 3.8+: Primary language
Pandas/NumPy: Data manipulation
Scikit-learn: Machine learning algorithms
TensorFlow/PyTorch: Deep learning models
Seaborn/Matplotlib: Data visualization
Jupyter: Interactive development
Docker: Containerization

Repository Structure
neurocardiac-comorbidity-prediction/
├── README.md                 # This file
├── requirements.txt          # Python dependencies
├── environment.yml           # Conda environment
├── Dockerfile               # Container setup
├── data/
│   ├── raw/                 # Original datasets (not tracked)
│   ├── processed/           # Cleaned datasets
│   ├── synthetic/           # Example/demo data
│   └── README.md           # Data documentation
├── src/
│   ├── preprocessing/
│   │   ├── rna_seq_pipeline.py    # RNA-seq data processing
│   │   ├── quality_control.py     # QC metrics and filtering
│   │   └── clustering.py          # Patient clustering
│   ├── models/
│   │   ├── multi_label_classifier.py  # ML classification
│   │   ├── multimodal_fusion.py       # Image + RNA integration
│   │   └── feature_selection.py       # Gene selection methods
│   ├── visualization/
│   │   ├── plots.py               # Plotting functions
│   │   └── dashboard.py           # Interactive dashboard
│   └── utils/
│       ├── data_loader.py         # Data loading utilities
│       └── evaluation.py          # Model evaluation
├── notebooks/
│   ├── 01_exploratory_analysis.ipynb
│   ├── 02_preprocessing.ipynb
│   ├── 03_model_development.ipynb
│   └── 04_results_analysis.ipynb
├── results/
│   ├── figures/               # Generated plots
│   ├── models/               # Trained model files
│   └── reports/              # Analysis reports
├── tests/
│   ├── test_preprocessing.py
│   ├── test_models.py
│   └── test_utils.py
└── docs/
    ├── methodology.md        # Detailed methods
    └── api_reference.md      # Code documentation
Installation
Prerequisites

Python 3.8 or higher
Git
(Optional) Conda for environment management

Setup Instructions

Clone the repository

bashgit clone https://github.com/yourusername/neurocardiac-comorbidity-prediction.git
cd neurocardiac-comorbidity-prediction

Create virtual environment

bash
# Using conda (recommended)
conda env create -f environment.yml
conda activate neurocardiac-pred

# Or using pip
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt

Download sample data (if available)

bash
# Follow instructions in data/README.md

Data Requirements
Input Data Format

RNA-seq data: Gene expression matrix (genes × samples)
Clinical labels: Multi-label format (cardiovascular, neurological)
Optional: Brain imaging data for multimodal approach

Data Privacy Notice
Important: This repository does not contain any real patient data. All examples use synthetic or publicly available datasets. Always ensure compliance with data protection regulations (HIPAA, GDPR) when working with medical data.

License
This project is licensed under the MIT License - see the LICENSE file for details.

Contact

Author: Diego Vachon Galindo
Email: diego.vachongalindo@mail.mcgill.ca
GitHub: @diegovachon