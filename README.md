# PFCE
A new parameter-free framework for calibration enhancement (PFCE) was proposed for dealing with the near-infrared (NIR) spectral inconsistency and maintaining the prediction ability of the calibration model under different conditions.  The calibration issues encountered in the maintenance with or without using standards, and even the enhancement between instruments have been thoroughly addressed.  The general calibration maintenance/enhancement cases were formulated into non-supervised PFCE (NS-PFCE), semi-supervised PFCE (SS-PFCE), and full-supervised PFCE (FS-PFCE).  The NS-PFCE made use of both the provided master and slave spectra of standard samples to construct a maintained calibration slave model by implementing a correlation constraint on the regression coefficients.  The SS-PFCE and FS-PFCE methods integrated the slave spectra and reference information of standard samples at the same time into the slave spectral calibration, and thus a maintenance or enhancement model could be achieved for the slave spectra, in particular measured on different instruments.

# Useage

---
Run demo_tablet.m script directly in Matlab environment
```matlab
    demo_tablet.m
```


# Graph Abstract

---
![Image](https://github.com/JinZhangLab/PFCE/blob/ccb9f7b9312c999593d3ae2357b53d1eb2f3b083/TOC.jpg)

---
# If the progrom is helpful for you, please cite the following paper:

[1] J Zhang, BY Li, Y Hu, et. al. A Parameter-Free Framework for Calibration Enhancement of Near-Infrared Spectroscopy Based on Correlation Constraint, analytica chimica acta, 2021, 1142: 169-178. https://doi.org/10.1016/j.aca.2020.11.006
