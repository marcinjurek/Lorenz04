VEnKF Introduction
------------------

- [Brief_Description](#Brief_Description)
  * [EnKF_Updates](#EnKF_Updates)
  * [VEcchia_Update](#Vecchia_Update)
  * [Lorenz_Model](#Lorenz_Model)
- [Examples_and_Simulations](#Examples_and_Simulations)
- [References](#References)



Brief Description
-----------------

This package is built to provide functionality to the method described
in my forthcoming paper with Dr. Matthias Katzfuss. The main purpose of
the package is the function vec.update, which calculates the update step
in an Ensemble Kalman Filter using Vecchia approximation to regularize
the covariance estimation.

Another purpose of this package is the other update functions for
Ensemble Kalman Filters. These are functions that provide basic update
steps to implement with other regularization methods for comparison
against our method as well as new methods in the future.

The final purpose of this package is the implementation of a Lorenz
chaos model for testing methods of data assimilation. The model is
described in Lorenz (2005) as Model II.

### EnKF Updates

The Ensemble Kalman Filter is common algorithm utilized in large data
assimilation applications based upon the state-space model where we have
a sample from a distribution
$\\hat{\\mathbf{x}}\_{t-1}^{(1)},\\dots,\\hat{\\mathbf{x}}\_{t-1}^{(N)}$.
The state-space model is typically broken up into two steps: an
Observation Model and an Evolution Model.

**y**<sub>*t*</sub> = **H**<sub>*t*</sub>**x**<sub>*t*</sub> + **v**<sub>*t*</sub>,   **v**<sub>*t*</sub> ∼ *N*<sub>*m*<sub>*t*</sub></sub>(0, **R**<sub>*t*</sub>)

**x**<sub>*t*</sub> = **M**<sub>*t*</sub>**x**<sub>*t* − 1</sub> + **w**<sub>*t*</sub>,   **w**<sub>*t*</sub> ∼ *N*<sub>*n*</sub>(0, **Q**<sub>*t*</sub>)

This model can also be thought of as a forecast step and an update step.
In most applications, the forecast step is considered a “black box”
algorithm for some process, so we do not have knowledge of the internal
mechanism and cannot with consistency improve upon the class of Ensemble
Kalman Filters by working on these “black box” forecast steps. Instead
we focus upon the update step:

**x**<sub>*i*</sub><sup>\*</sup> = *g*(**x**<sub>*i*</sub>|*Σ*) = *Σ*<sup>\*</sup>(*Σ*<sup> − 1</sup>**x**<sub>*i*</sub> + **H**′**R**<sup> − 1</sup>**y**<sub>*i*</sub>)

where

*Σ*<sup>\*</sup> = (*Σ*<sup> − 1</sup> + **H**′**R**<sup> − 1</sup>**H**)<sup> − 1</sup>

This provides the framework for various regularized update steps
including tapering and localization on the estimated covariance matrix.
It is also possible, but not recommended for large datasets, to utilize
the unregularized sample covariance matrix.

### Vecchia Update

### Lorenz Model

Examples and Simulations
------------------------

References
----------

Lorenz, Edward N. 2005. “Designing chaotic models.” *Journal of the
Atmospheric Sciences* 62 (5): 1574–87.
<https://doi.org/10.1175/JAS3430.1>.
