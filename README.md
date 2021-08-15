# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints

BOB design is a Bayesian optimal design proposed for biosimilar trials with co-primary endpoints. 

This repository contains R codes used to implement **numerical studies** and the **trial application** in the corresponding paper.

## Table of Contents

- [Numerical Studies](#numerical studies)
  - [Design Calibration](#Design Calibration)
  - [Type I error rate control](#Type I error rate Control)
  - [Empirical power](#Empirical power)
- [Trial Application](#trial application)

## Numerical Studies

### Design Calibration

We calibrate adaptive designs including BAE, BAS, $BOB_s$,  and $BOB_{avg}$.



* calibrate_1: $mu_T=\pm 0.223$, $p_T \sim unif(0.35,0.65)$
* calibrate_2: muT=0, pT=0 (i.e., power of the design)
* calibrate_3: $mu_T \sim unif(-0.223,0.223)$, $p_T=0.5\pm 0.15$

### Type I error rate Control



### Empirical power







## Trial Application

This module depends upon a knowledge of [Markdown]().

```
```



## Usage

```
```

Note: The `license` badge image link at the top of this file should be updated with the correct `:user` and `:repo`.





## Contributing

See [the contributing file](CONTRIBUTING.md)!

PRs accepted.

Small note: If editing the Readme, please conform to the [standard-readme](https://github.com/RichardLitt/standard-readme) specification.

### Any optional sections

## License

[MIT © Richard McRichface.](../LICENSE)
