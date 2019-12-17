"""
Fluid Inclusion Modelling with Kuster-Toksoz Method (simplification of Wang, 1997)
By: Yohanes Nuwara
"""

import numpy as np

def stuffs(Km, Kf, Gm, Gf, totalporo, compos, alpha):
    A = Gf / Gm - 1.0
    B = (Kf / Km - Gf / Gm) / 3.0
    R = Gm / (Km + (4.0 / 3.0) * Gm)
    Fm = (Gm / 6.0) * (9.0 * Km + 8.0 * Gm) / (Km + 2.0 * Gm) #zeta
    ci = compos*totalporo
    theta = alpha * (np.arccos(alpha) - alpha * np.sqrt(1.0 - alpha * alpha)) / (1.0 - alpha * alpha) ** (3.0 / 2.0)
    f = alpha * alpha * (3.0 * theta - 2.0) / (1.0 - alpha * alpha)
    return(A, B, R, Fm, ci, theta, f)

def PQ(A, B, R, theta, f):
    F1 = 1.0 + A * (1.5 * (f + theta) - R * (1.5 * f + 2.5 * theta - 4.0 / 3.0))
    F2 = 1.0 + A * (1.0 + 1.5 * (f + theta) - R * (1.5 * f + 2.5 * theta)) + B * (3.0 - 4.0 * R) + A * (A + 3.0 * B) * (
        1.5 - 2.0 * R) * (f + theta - R * (f - theta + 2.0 * theta * theta))
    F3 = 1.0 + A * (1.0 - f - 1.5 * theta + R * (f + theta))
    F4 = 1.0 + (A / 4.0) * (f + 3.0 * theta - R * (f - theta))
    F5 = A * (-f + R * (f + theta - 4.0 / 3.0)) + B * theta * (3.0 - 4.0 * R)
    F6 = 1.0 + A * (1.0 + f - R * (f + theta)) + B * (1.0 - theta) * (3.0 - 4.0 * R)
    F7 = 2.0 + (A / 4.0) * (3.0 * f + 9.0 * theta - R * (3.0 * f + 5.0 * theta)) + B * theta * (3.0 - 4.0 * R)
    F8 = A * (1.0 - 2.0 * R + (f / 2.0) * (R - 1.0) + (theta / 2.0) * (5.0 * R - 3.0)) + B * (1.0 - theta) * (
        3.0 - 4.0 * R)
    F9 = A * ((R - 1.0) * f - R * theta) + B * theta * (3.0 - 4.0 * R)
    P = 3.0 * F1 / F2
    Q = 2.0 / F3 + 1.0 / F4 + (F4 * F5 + F6 * F7 - F8 * F9) / (F2 * F4)
    return(P,Q)

def KusterToksoz(sigma_P, sigma_Q, Km, Gm, Kf, rhom, rhof):
    K_sat = ((3 * Km * (3 * Km + 4 * Gm)) + (4 * Gm * (Kf - Km) * sigma_P)) / ((3 * (3 * Km + 4 * Gm)) -
        (3 * (Kf - Km) * sigma_P))
    G_sat = ((25 * (Gm**2) * (3 * Km + 4 * Gm)) - ((Gm**2) * (9 * Km + 8 * Gm) * sigma_Q)) / ((25 * Gm * (3 * Km + 4 * Gm))
        + (6 * Gm * (Km + 2 * Gm) * sigma_Q))
    rho_sat = (rhom*(1-0.14))+(rhof*0.14)
    Vp_sat = np.sqrt((K_sat + 4/3*G_sat) / rho_sat)
    Vs_sat = np.sqrt(G_sat / rho_sat)
    return(K_sat, G_sat, rho_sat, Vp_sat, Vs_sat)
