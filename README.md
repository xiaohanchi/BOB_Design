# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints

BOB design is a Bayesian optimal design proposed for biosimilar trials with co-primary endpoints. 

This repository contains R codes used to implement **numerical studies** and the **trial application** in the corresponding paper.

## Table of Contents

- [Numerical Studies](#numerical-studies)
  - [Fixed Designs](#Fixed-Designs)
  - [Bayesian Adaptive Designs](#Bayesian-Adaptive-Designs)
    - [Design Calibration](#Design-Calibration)
    - [power(type I error)](#power(type-I-error))
- [Trial Application](#trial-application)
- Simulation Settings

## Numerical Studies

### Fixed Designs

description xxx

* FE.R:
* FS.R:
* FES.R:

### Bayesian Adaptive Designs

description xxx

#### Design Calibration

We calibrate adaptive designs including BAE, BAS, $BOB_s$,  and $BOB_{avg}$.

* calibration_bae.R: R codes used to calibrate the design **BAE**, and output the optimal design parameters.
* calibration_bas.R: R codes used to calibrate the design **BAS**, and output the optimal design parameters.
* calibration_bobs.R: R codes used to calibrate the design **BOBs**, and output the optimal design parameters.

* BOBavg

  This folder contains 3 files used to implement the whole calibration procedure of the design **BOBavg** with the following settings: 

  * search_1.R: $\mu_T=\pm 0.223$, $p_T \sim unif(0.35,0.65)$
  * search_2.R: $\mu_T=0, p_T=0$ (i.e., power of the design)
  * search_3.R: $\mu_T \sim unif(-0.223,0.223)$, $p_T=0.5\pm 0.15$

  and the file output.R used to output the optimal design parameters.

#### **power(type I error)**



* simu_bae.R:
* simu_bas.R:
* simu_bob.R:





## Trial Application





## Authors and Reference

* Xiaohan Chi, Zhangsheng Yu, and Ruitao Lin

  
