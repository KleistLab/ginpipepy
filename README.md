# ginpipepy - Package for temporal binning of dated sequences in BAM format and fast population size estimate computation

[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
[![DOI:10.1101/2021.05.14.21257234](http://img.shields.io/badge/DOI-10.1101/2021.05.14.21257234-blue.svg)](https://doi.org/10.1101/2021.05.14.21257234)

This Python package is used by GInPipe [[1]](#1) - Genome-based Incidence Estimation Pipeline. Functions implmented in *ginpipepy* are used to split aligned sequences in BAM format into chronologically sorted bins of predefined size and temporal length, and to calculate an estimate of population size. The package also includes a module for Variant Calling File reading that contains sites that have to be masked. 

## Dependencies

*ginpipepy* uses the following packages:

  - numpy
  - pysam
  - biopython
  - pandas
  - scipy
  - pyvcf

## Installation

*ginpipepy* can be installed using **pip**:

```
pip install git+https://github.com/trofimovamw/ginpipepy
```
 
## Reference
<a id="1">[1]</a> 
Smith, M. R. and Trofimova, M., et al. (2021). Rapid incidence estimation from SARS-CoV-2 genomes reveals decreased case detection in Europe during summer 2020. medRxiv, https://github.com/KleistLab/GInPipe

