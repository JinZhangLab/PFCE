# PFCE: A Parameter-Free Framework for Calibration Enhancement of Near-Infrared Spectroscopy

## Overview

PFCE is a versatile parameter-free calibration enhancement framework designed for near-infrared (NIR) spectra. This framework is capable of addressing the inconsistency issue in NIR spectra and maintaining the prediction ability of the calibration model under different conditions. 

## Available Methods

The following methods are included in the PFCE framework:
- Non-supervised PFCE (NS-PFCE)
- Semi-supervised PFCE (SS-PFCE)
- Full-supervised PFCE (FS-PFCE)
- Multi-temporal PFCE (MT-PFCE)

### Non-supervised PFCE (NS-PFCE)
NS-PFCE makes use of both the provided master and slave spectra of standard samples to construct a maintained calibration slave model by implementing a correlation constraint on the regression coefficients. This method is suitable for cases where both master and slave standard spectra are available, but no reference information of the standard samples is provided.

### Semi-supervised PFCE (SS-PFCE)
SS-PFCE integrates the slave spectra and reference information of standard samples into the slave spectral calibration. This method can be used when just slave spectra and reference information of a batch of samples are available. 

### Full-supervised PFCE (FS-PFCE)
FS-PFCE is similar to NS-PFCE, but the reference information of the standard samples is also included in the construction of the slave model. This method is appropriate when reference information is available for both the master and slave spectra.

### Multi-temporal PFCE (MT-PFCE)
MT-PFCE is designed for the most general scenario where multiple instruments require calibration enhancement. This method can be used to simultaneously enhance the calibration of both the master and slave models. MT-PFCE can be adapted to a wide range of scenarios, such as SS-PFCE or FS-PFCE, by adjusting the conditions for the master model. Additionally, this method has the potential to model data with inherent hierarchical structures, similar to pre-training and fine-tuning in deep learning. This is an area for future work.

## Required Data
- Non-supervised PFCE (NS-PFCE): Master spectra and slave spectra from the same batch of standard samples
- Semi-supervised PFCE (SS-PFCE): Slave spectra and the corresponding reference information from a batch of samples
- Full-supervised PFCE (FS-PFCE): Master spectra, slave spectra, and the corresponding reference information of batch of standard samples
- Multi-temporal PFCE (MT-PFCE):  Spectra and corresponding reference information for each task, not necessarily from the same batch of samples.

## Usage
This code is implemented in Matlab. To run the demo, simply execute the `demo_tablet.m` script in the Matlab environment.

```matlab
    demo_tablet.m
```
The PFCEs can be also used with python package at [pynir](https://pypi.org/project/pynir/), or online web server at [nir.chemoinfolab.com](nir.chemoinfolab.com).

## Successful applications
| Method | Application | Reference |
|--------|-------------|----------|
| SS-PFCE | Batch-to-batch variability of fruit samples in portable spectroscopy | [Mishra, 2021](https://www.sciencedirect.com/science/article/pii/S0003267021005973) |
| NS- and SS-PFCE | Calibration transfer from point spectrometers to visible/near-infrared spectral imaging machines | [Mishra, 2021](https://www.sciencedirect.com/science/article/pii/S0003267021009806)|
| SS-PFCE | Baseline for the deep transfer leanring of NIR spectral calibration of mango dry matter and melamine turbidity point prediction| [Mishra, 2021](https://www.sciencedirect.com/science/article/pii/S0169743921000514)|

## Citation
If you find this method useful in your research, we kindly ask that you cite the following paper:

[1] J Zhang, BY Li, Y Hu, et. al. A Parameter-Free Framework for Calibration Enhancement of Near-Infrared Spectroscopy Based on Correlation Constraint, analytica chimica acta, 2021, 1142: 169-178. https://doi.org/10.1016/j.aca.2020.11.006
