"""Predictive equations for the cyclic backbone parameters and the hysteretic
energy capacity of ductile RC beam-column members.

Implements the models calibrated in:

    Hasanoglu, S. and O'Reilly, G. J. (2026). Calibration of cyclic backbone
    parameters and hysteretic energy capacity for numerical modelling of
    ductile reinforced concrete members. (UNDER REVIEW)

Symbols
-------
fc        concrete compressive strength [MPa]
fyl       yield strength of the longitudinal reinforcement [MPa]
nu        axial load ratio, N / (fc * Ag), positive in compression
rho_l     longitudinal reinforcement ratio
rho_t     transverse reinforcement ratio
Ls        shear span length [m]
Ls_d      shear span length to effective depth ratio
d         effective depth [m]
h         section depth [m]
st        transverse reinforcement spacing in the confinement zone [m]
My        yield moment [kNm]
EIg       gross section flexural stiffness [kNm2]

Reinforcement and axial load ratios are dimensionless fractions, not
percentages (e.g. 0.04 instead of 4%). Rotations are in radians and the
hysteretic energy capacity in kNm.
"""

import math
import numpy as np


def elastic_modulus(fc):
    """Modulus of elasticity of the concrete [kPa].
    """
    Ecm_mpa = 22000.0 * (fc / 10.0) ** 0.3
    Ecm_kpa = 1000 * Ecm_mpa

    return Ecm_kpa


def effective_depth(h, cover, dbh, dbl):
    """Effective depth d [m] of a section with one tension layer.
    """

    return h - cover - dbh - 0.5 * dbl


def yield_moment(b, h, nbl_int, nbl_side, dbl, dbh, cover, fc, fyl, nu):
    """Yield moment [kNm] of a rectangular section, after Panagiotakos and
    Fardis (2001).

    b, h         section width and depth [m], h in the bending direction
    nbl_int      intermediate bars in the tension layer, excluding the two
                 corner bars of that layer
    nbl_side     intermediate bars along one side face
    dbl, dbh     longitudinal and transverse bar diameters [m]
    cover        cover to the transverse reinforcement [m]

    The section is assumed symmetrically reinforced, and yielding is taken to
    be controlled by the tension reinforcement or by the non-linearity of the
    concrete in compression, whichever governs.
    """

    fc_kpa = fc * 1000.0
    fyl_kpa = fyl * 1000.0
    axial_force = nu * fc_kpa * b * h

    # Depth of the equivalent rectangular stress block, ACI 318.
    if fc < 27.6:
        beta_c = 0.85
    elif fc > 55.17:
        beta_c = 0.65
    else:
        beta_c = 1.05 - 0.05 * fc / 6.9

    Ec = elastic_modulus(fc)
    Es = 2e8
    n_young = Es / Ec
    eps_sy = fyl_kpa / Es

    d = effective_depth(h, cover, dbh, dbl)
    d_prime = h - d
    delta = d_prime / d

    bar_area = 0.25 * math.pi * dbl ** 2

    # Two corner bars per layer, plus the intermediate ones. The section is
    # symmetric, so the tension and compression layers are identical.
    As_tens = (2 + nbl_int) * bar_area
    As_comp = As_tens
    As_side = 2 * nbl_side * bar_area

    rho_tens = As_tens / (b * d)
    rho_comp = As_comp / (b * d)
    rho_side = As_side / (b * d)

    # Neutral axis depth, against its value at balanced failure.
    c = (As_tens * fyl_kpa - As_comp * fyl_kpa + axial_force) \
        / (0.85 * fc_kpa * b * beta_c)
    c_balanced = 0.0035 * d / (0.0035 + eps_sy)

    concrete_controlled = c >= c_balanced

    if concrete_controlled:
        A = rho_tens + rho_comp + rho_side \
            - axial_force / (1.8 * n_young * b * d * fc_kpa)
        B = rho_tens + rho_comp * delta + 0.5 * rho_side * (1 + delta)
    else:
        A = rho_tens + rho_comp + rho_side \
            + axial_force / (b * d * fyl_kpa)
        B = rho_tens + rho_comp * delta + 0.5 * rho_side * (1 + delta) \
            + axial_force / (b * d * fyl_kpa)

    # Neutral axis depth at yielding, normalised by the effective depth.
    ky = math.sqrt(n_young ** 2 * A ** 2 + 2.0 * n_young * B) - n_young * A

    if concrete_controlled:
        phi_y = 1.8 * fc_kpa / (Ec * ky * d)
    else:
        phi_y = fyl_kpa / (Es * (1.0 - ky) * d)

    concrete = Ec * ky ** 2 / 2.0 * (0.5 * (1.0 + delta) - ky / 3.0)
    steel = Es / 2.0 * ((1.0 - ky) * rho_tens
                        + (ky - delta) * rho_comp
                        + rho_side / 6.0 * (1.0 - delta)) * (1.0 - delta)

    My = b * d ** 3 * phi_y * (concrete + steel)

    return My


def stiffness_ratio(nu, Ls_d):
    """Secant stiffness at yield over gross stiffness, Equation (1).
    Bounded to 0.2 <= EIy/EIg <= 0.8.
    """
    stiff_ratio = np.clip(0.08 * (1.25 ** (10*nu)) * (1.2 ** Ls_d), 0.2, 0.8)

    return stiff_ratio


def yield_rotation(My, Ls, EIg, EIy_EIg):
    """Yield rotation theta_y
    """
    EIy = EIy_EIg * EIg

    return My * Ls / (3.0 * EIy)


def plastic_rotation(fc, nu, rho_l):
    """Plastic rotation prior to capping, theta_pl, Equation (2).
    Valid for symmetric reinforcement; use asymmetric_correction otherwise.
    """
    return 0.009 * 0.25 ** (0.01 * fc) * 0.23 ** nu * 2.17 ** (100.0 * rho_l)


def asymmetric_correction(theta_pl, rho_comp, rho_tens, fyl, fc):
    """Correct theta_pl for asymmetric reinforcement, Equation (3).

    Adopted from Fardis and Biskinis (2003); rho_comp and rho_tens are the
    compressive and tensile reinforcement ratios.
    """
    compressive = rho_comp * fyl / fc
    tensile = rho_tens * fyl / fc

    return (max(0.01, compressive) / max(0.01, tensile)) ** 0.225 * theta_pl


def post_capping_rotation(nu, rho_t, fyl):
    """Capping to ultimate rotation, theta_pc, Equation (4).
    """

    return 0.01 * 0.07 ** nu * 1.78 ** (100.0 * rho_t) * 1.18 ** (0.01 * fyl)


def post_ultimate_rotation(fc, rho_t, theta_pc=None):
    """Ultimate to zero-moment rotation, theta_pu, Equation (5).

    Given theta_pc, the result is capped at 4 * theta_pc so that the final
    branch stays at least as steep as the post-capping one.
    """
    value = 0.0023 * 1.96 ** (0.1 * fc) * 2.74 ** (100.0 * rho_t)

    if theta_pc is None:
        return value

    return min(value, 4 * theta_pc)


def hardening_ratio(nu, rho_l, Ls_d):
    """Post-yield hardening ratio Mcap/My, Equation (6).
    """

    return 1.08 ** nu * 1.38 ** (10.0 * rho_l) * 1.14 ** (0.1 * Ls_d)


def softening_ratio(fc, Ls_d):
    """Post-capping softening ratio Mult/Mcap, Equation (7).
    """

    return 0.91 ** (0.01 * fc) * 0.80 ** (0.1 * Ls_d)


def energy_capacity_ratio(fc, fyl, nu, st, h):
    """Dimensionless energy capacity lambda, Equation (9).

    Defined as lambda = Eh_cap / (My * theta_pl_total).
    """

    return 52.0 * 0.41 ** (0.01 * (fc + 0.1 * fyl)) \
        * 0.34 ** nu * 0.21 ** (st / h)


def energy_capacity(My, theta_pl, theta_pc, lambda_):
    """Hysteretic energy capacity Eh_cap [kNm], Equation (8).
    """
    Eh_cap = My * (theta_pl + theta_pc) * lambda_

    return Eh_cap


def cyclic_backbone(b, h, nbl_int, nbl_side, dbl, dbh, cover, st,
                    fc, fyl, nu, rho_l, rho_t, Ls,
                    rho_comp=None, rho_tens=None):
    """Evaluate the cyclic backbone of a member from its section and loading.
    """

    # Effective section depth
    d = effective_depth(h, cover, dbh, dbl)

    # Gross section stiffness
    EIg = elastic_modulus(fc) * b * h ** 3 / 12.0

    # Shear span to effective depth ratio
    Ls_d = Ls / d

    # EIy_EIg
    EIy_EIg = stiffness_ratio(nu, Ls_d)

    # yield moment
    My = yield_moment(b, h, nbl_int, nbl_side, dbl,
                      dbh, cover, fc, fyl, nu)

    # Yield rotation
    theta_y = yield_rotation(My, Ls, EIg, EIy_EIg)

    # Plastic rotation
    theta_pl = plastic_rotation(fc, nu, rho_l)
    if rho_comp is not None and rho_tens is not None:
        theta_pl = asymmetric_correction(theta_pl, rho_comp, rho_tens, fyl, fc)

    # Post-capping rotation
    theta_pc = post_capping_rotation(nu, rho_t, fyl)

    # Post-ultimate rotation
    theta_pu = post_ultimate_rotation(fc, rho_t, theta_pc)

    # Post-yield hardening ratio
    Mcap_My = hardening_ratio(nu, rho_l, Ls_d)

    # Post-capping softening ratio
    Mult_Mcap = softening_ratio(fc, Ls_d)

    # Lambda for Eh,cap
    lambda_ = energy_capacity_ratio(fc, fyl, nu, st, h)

    # Rotation at capping
    theta_cap = theta_y + theta_pl

    # Rotation at ultimate point
    theta_ult = theta_cap + theta_pc

    # Moment at capping point
    Mcap = Mcap_My * My

    # Hysteretic energy capacity
    Eh_cap = energy_capacity(My, theta_pl, theta_pc, lambda_)

    return {'nu': nu,
            'EIg': EIg,
            'EIy_EIg': EIy_EIg,
            'theta_pl': theta_pl,
            'theta_pc': theta_pc,
            'theta_pu': theta_pu,
            'Mcap_My': Mcap_My,
            'Mult_Mcap': Mult_Mcap,
            'lambda': lambda_,
            'Eh_cap': Eh_cap,
            'theta_y': theta_y,
            'theta_cap': theta_cap,
            'theta_ult': theta_ult,
            'theta_zero': theta_ult + theta_pu,
            'My': My,
            'Mcap': Mcap,
            'Mult': Mult_Mcap * Mcap}


if __name__ == '__main__':
    # A 400 x 400 mm column with 3 + 3 longitudinal bars of 20 mm in the
    # bending direction and one bar on each side face, 1.5 m shear span,
    # under an axial load ratio of 0.2.

    parameters = cyclic_backbone(b=0.4, h=0.4, nbl_int=1, nbl_side=1,
                                 dbl=0.020, dbh=0.008, cover=0.03, st=0.10,
                                 fc=30, fyl=450.0, nu=0.2,
                                 rho_l=0.015, rho_t=0.008, Ls=1.5)

    print('\nBackbone parameters')
    for key, unit in (('My', 'kNm'), ('EIy_EIg', '-'),
                      ('theta_pl', 'rad'), ('theta_pc', 'rad'),
                      ('theta_pu', 'rad'), ('Mcap_My', '-'),
                      ('Mult_Mcap', '-'), ('lambda', '-'),
                      ('Eh_cap', 'kNm')):
        print('{:<10s} {:9.4f}  {}'.format(key, parameters[key], unit))

    print('\nBackbone points [rad, kNm]')
    for theta, moment in (('theta_y', 'My'), ('theta_cap', 'Mcap'),
                          ('theta_ult', 'Mult'), ('theta_zero', None)):
        print('  {:8.5f} {:8.1f}'.format(
            parameters[theta], parameters[moment] if moment else 0.0))
