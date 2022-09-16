# [FluTO: GRADED MULTISCALE FLUID TOPOLOGY OPTIMIZATION USING NEURAL NETWORKS]()

[Rahul Kumar Padhy](https://sites.google.com/view/rahulkp/home), [Aaditya Chandrasekhar*](https://aadityacs.github.io/), [Krishnan Suresh](https://directory.engr.wisc.edu/me/faculty/suresh_krishnan)  
University of Wisconsin-Madison


## Abstract

Fluid-flow devices with low dissipation, but high contact area, are of importance in many applications. A well-known strategy to design such devices is multi scale topology optimization (MTO), where optimal microstructures are designed within each cell of a discretized domain. Unfortunately, MTO is computationally very expensive since one must perform homogenization of the evolving microstructures, during each step of the homogenization process. As an alternate, we propose here a
graded multiscale topology optimization (GMTO) for designing fluid-flow devices. In the proposed method, several pre-selected but size-parameterized and orientable microstructures are used to fill the domain optimally. GMTO significantly reduces the computation while retaining many of the benefits of MTO.


In particular, GMTO is implemented here using a neural-network (NN) since: (1) homogenization can be performed off-line, and used by the NN during optimization, (2) it enables continuous switching between microstructures during optimization, (3) the number of design variables and computational effort is independent of number of microstructure used, and, (4) it supports automatic differentiation, thereby eliminating manual sensitivity analysis. Several numerical results are presented to illustrate the proposed framework.

## Citation

```
@article{Padhy2022FluTO,
author = {Padhy, Rahul Kumar and Chandrasekhar, Aaditya and Suresh, Krishnan},
title = {FluTO: GRADED MULTISCALE FLUID TOPOLOGY OPTIMIZATION USING NEURAL NETWORKS},
journal = {},
year={2022}
}
```

*contributed equally
