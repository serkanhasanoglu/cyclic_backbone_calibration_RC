# Calibration of Cyclic Backbone Parameters and Hysteretic Energy Capacity for Ductile RC Members

Supplementary data for the paper referenced below. The repository collects the filtered
experimental database, the per-specimen backbone curve identifications, the calibration
figures of the `HystereticSM` material model used in the study, and an implementation of
the proposed equations.

## Reference

Hasanoglu, S. and O'Reilly, G. J. (2026). *Calibration of cyclic backbone parameters and
hysteretic energy capacity for numerical modelling of ductile reinforced concrete members.*
(UNDER REVIEW)

## Repository contents

| Path | Description |
| --- | --- |
| [database_filtered.xlsx](database_filtered.xlsx) | Filtered database of the 115 flexure-dominated reversed-cyclic tests on rectangular RC columns used for the calibration (geometry, reinforcement detailing, material properties, axial load and source reference of each specimen). |
| [Backbone_Curves/](Backbone_Curves/) | Identified cyclic backbone curve of each test specimen plotted over its experimental hysteresis (one PDF per specimen). |
| [HystereticSM_cal_figures/](HystereticSM_cal_figures/) | Calibration figures of the `HystereticSM` material against the experimental hysteresis of each specimen, collected in `All_Figures.pdf`. |
| [cyclic_backbone.py](cyclic_backbone.py) | Implementation of the proposed equations for the six cyclic backbone parameters and the hysteretic energy capacity. |

## Using the calibrated models

`cyclic_backbone.py` implements the proposed equations and requires `numpy`. Material
strengths are in MPa, lengths in metres, moments in kNm, and reinforcement and axial load
ratios are dimensionless fractions rather than percentages.

```python
from cyclic_backbone import cyclic_backbone

parameters = cyclic_backbone(b=0.40, h=0.40,          # section width and depth
                             nbl_int=1, nbl_side=1,   # intermediate bars, layer / side
                             dbl=0.020, dbh=0.008,    # longitudinal, stirrup diameters
                             cover=0.03, st=0.10,     # cover, stirrup spacing
                             fc=30.0, fyl=450.0,
                             nu=0.2,                  # axial load ratio, N / (fc * Ag)
                             rho_l=0.0157, rho_t=0.008,   # 8 bars of 20 mm
                             Ls=1.5)                  # shear span

parameters['My']          # yield moment [kNm]
parameters['theta_pl']    # plastic rotation prior to capping [rad]
parameters['theta_cap']   # rotation at the capping point [rad]
parameters['Eh_cap']      # hysteretic energy capacity [kNm]
```

## Range of applicability

The equations were calibrated on 115 flexure-dominated reversed-cyclic tests on rectangular
RC columns, selected from the database with the criteria below. Applying them outside these
ranges is not recommended.

| Parameter | Range |
| --- | --- |
| Concrete compressive strength, `fc` | < 80 MPa |
| Yield strength of the reinforcement, `fyl` and `fyt` | < 600 MPa |
| Axial load ratio, `nu` | <= 0.7 |
| Longitudinal reinforcement ratio, `rho_l` | <= 0.04 |
| Shear span to effective depth ratio, `Ls/d` | 2.0 to 7.34 |

The equations also assume flexure-dominated, ductile behaviour: specimens failing in shear
or in combined flexure-shear were excluded from the calibration.
