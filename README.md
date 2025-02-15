# single-cell-integration-pipeline

This repository is a beginner-friendly demonstration of a **Single Cell Integration Pipeline**, primarily following the **Scanpy** tutorials.

## Overview

The project aims to explore basic preprocessing, quality control, and clustering methods for single-cell RNA sequencing (scRNA-seq) data using **Scanpy**. It serves as a learning exercise for understanding fundamental steps in single-cell analysis.

## Sources

The implementation is based on the following tutorials:

- **Quality Control & Preprocessing**: [SC Best Practices](https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html)  
- **Clustering & Basic Analysis**: [Scanpy Clustering Tutorial](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html#)  

## Requirements

To set up the environment, use the provided `environment.yml` file:

```bash
conda env create -f environment.yml
conda activate single_cell_env
python -m ipykernel install --user --name=single_cell_env --display-name "Python (single_cell_env)"
```