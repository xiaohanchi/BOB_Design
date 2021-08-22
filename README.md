# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints

BOB design is a Bayesian optimal design proposed for biosimilar trials with co-primary endpoints. 

This repository contains R codes used to implement **numerical studies** and the **trial application** in the corresponding paper.

## Table of Contents

- [Simulation Settings](#Simulation-Settings)
- [Numerical Studies](#numerical-studies)
  - [Fixed Designs](#Fixed-Designs)
  - [Bayesian Adaptive Designs](#Bayesian-Adaptive-Designs)
    - [Design Calibration](#Design-Calibration)
    - [power(type I error)](#power(type-I-error))
- [Trial Application](#trial-application)

## Simulation Settings

Firstly, some important setting parameters and their meanings are explained in this document, as follows:

* maxnsample: $n_J$
* Tmax: $J$
* nsample: $n_j$
* tau2: $\tau^2$
* rho: $\rho$
* sn: number of replicated trials in this simulation study
* pR (pT): $p_R (p_T)$
* overallmuR (overallmuT) $\mu_R (\mu_T)$

## Numerical Studies

### Fixed Designs

This folder contains R codes used to implement three fixed-sample designs considered in the paper.

Fixed-sample designs includes two univariate designs such as FE and FS and a bivariate design FES. Among them, **FE** and **FS** are frequentist fixed-sample designs that adopt the frequentist two one-sided tests (TOST) procedure for both sides to test the respective efficacy and safety endpoints, and **FES** is a frequentist fixed-sample design that combines the FE and FS designs to test both efficacy and safety endpoints. 

### Bayesian Adaptive Designs

This folder contains R codes used to implement four Bayesian adaptive designs considered in the paper. The procedure of Bayesian designs requires two main steps: (1) design calibration and (2) design implementation.

#### Design Calibration

R codes in this folder help us to calibrate the design and return the optimal design parameters ($\lambda, \gamma$). We calibrate Bayesian adaptive designs including BAE, BAS, BOBs,  and BOBavg. In detail, BAE and BAS are Bayesian group-sequential designs that consider the respective efficacy and safety as a single primary endpoint. BOBs  and BOBavg are proposed BOB designs.

* calibration_bae.R: R codes used to calibrate the design **BAE**, and output the optimal design parameters.

* calibration_bas.R: R codes used to calibrate the design **BAS**, and output the optimal design parameters.

* calibration_bobs.R: R codes used to calibrate the design **BOBs**, and output the optimal design parameters.

* BOBavg

  This folder contains 3 files used to implement the whole calibration procedure of the design **BOBavg** with the following settings: 

  * search_1.R: $\mu_T=\pm 0.223$, $p_T \sim unif(0.35,0.65)$
  * search_2.R: $\mu_T=0, p_T=0$ (i.e., power of the design)
  * search_3.R: $\mu_T \sim unif(-0.223,0.223)$, $p_T=0.5\pm 0.15$

  and the file output.R used to output the optimal design parameters.

For example, for the proposed BOBs design, simply run the corresponding R script like

```shell
Rscript calibration_bobs.R
```

the optimal design parameters will be output as follows:

```R
#Optimal parameters for BOBs:(0.9408,1.06)
```



#### Design implementation: **power(type I error)**

Given the resulted optimal design parameters, simulation can be performed to obtain the operating characteristics such as the power (or the type I error rate) and the expected sample size.

For example, for the proposed BOBs design, 

```shell
Rscript simu_bob.R
```

```R
#The power (or type I error rate) of the design BOB:
#Expected Sample Size(EN):
```



## Trial Application





## Authors and Reference

* Xiaohan Chi, Zhangsheng Yu, and Ruitao Lin

  
