# Calibration of Cyclic Backbone Parameters and Hysteretic Energy Capacity for Ductile RC Members

Supplementary data for the paper referenced below. The repository collects the filtered
experimental database, the per-specimen backbone curve identifications, the calibration
figures of the `HystereticSM` material model used in the study, and an implementation of
the proposed equations.

## Reference

Hasanoglu, S. and O'Reilly, G. J. (2026). *Calibration of cyclic backbone parameters for
numerical modelling and hysteretic energy capacity of ductile reinforced concrete members.*
(UNDER REVIEW)

## Repository contents

| Path | Description |
| --- | --- |
| [database_filtered.xlsx](database_filtered.xlsx) | Filtered database of the 115 flexure-dominated reversed-cyclic tests on rectangular RC columns used for the calibration (geometry, reinforcement detailing, material properties, axial load and source reference of each specimen). |
| [Backbone_Curves/](Backbone_Curves/) | Identified cyclic backbone curve of each test specimen plotted over its experimental hysteresis (one PDF per specimen). |
| [HystereticSM_cal_figures/](HystereticSM_cal_figures/) | Calibration figures of the `HystereticSM` material against the experimental hysteresis of each specimen, collected in `All_Figures.pdf`. |
| [cyclic_backbone.py](cyclic_backbone.py) | Implementation of the proposed equations for the six cyclic backbone parameters and the hysteretic energy capacity. |

## Using the calibrated models

Material strengths are in MPa, lengths in metres, moments in kNm, and reinforcement and axial load ratios are dimensionless fractions rather than percentages.

```python
from cyclic_backbone import cyclic_backbone

parameters = cyclic_backbone(
    b=0.40, h=0.40, nbl_int=1, nbl_side=1, dbl=0.020, dbh=0.008,
    cover=0.03, st=0.10, fc=30.0, fyl=450.0, nu=0.2,
    rho_l=0.0157, rho_t=0.008, Ls=1.5)

parameters['My']          # yield moment [kNm]
parameters['theta_pl']    # plastic rotation prior to capping [rad]
parameters['theta_cap']   # rotation at the capping point [rad]
parameters['Eh_cap']      # hysteretic energy capacity [kNm]
```

### Use in SimDesign

The calibrated equations are also implemented in
[SimDesign](https://github.com/builtenvdata/simulated-design), the simulated design and
non-linear modelling framework for RC frame buildings, as the concentrated plasticity
modelling option `CP02` (`model="CP02"`). In that option, the beam and column end hinges
use the `HystereticSM` material model.
