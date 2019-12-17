"""
Rock Physics Modelling Program for CO2-EOR / CCS
By: Yohanes Nuwara
"""

import numpy as np
import matplotlib.pyplot as plt
from Batzle_and_Wang import BW_brine_density, BW_brine_bulk
from Kuster_Toksoz_2 import stuffs, PQ, KusterToksoz
from EOS_CO2 import EOS
from RockModel_maker import rockmod

"""
Input
"""

"Step 0. Matrix property calculation"
rhob_log = 2.48 #brine-saturated rock from RHOB log

# mineral elastic property database
rhocalc = 2.71; Kcalc = 76.8; Gcalc = 32
rhoclay = 2.58; Kclay = 20.9; Gclay = 6.9
rhodolo = 2.87; Kdolo = 94.9; Gdolo = 45
rhoqtz = 2.65; Kqtz = 36.6; Gqtz = 45
# mineral composition
calc = 0.4; clay = 0.24; dolo = 0.34; qtz = 0.02

"Step 1. Brine elastic property calculation"
temp = 158  # in celsius
Pp_baseline = 34  # pore pressure before injection in MPa / status 100% water
salinity = 0.011  # in percent

"Step 2. Brine-saturated rock elastic property calculation"

# CO2 properties
SG = 1.519
# porosity
totalporo = 0.14
alpha1 = 0.5  # aspect ratio of 1st pore
alpha2 = 0.1  # aspect ratio of 2nd pore
compos1 = 0.75  # composition of 1st pore relative to total poro
compos2 = 0.25  # composition of 2nd pore relative to total poro

"Step 4. Pore pressure calculation"
inj_rate = 800  # ton/day, assume constant
maxyear_monitor = 10 #maximum years of monitoring
period = np.linspace(1, maxyear_monitor, 10)  # years
bulkvol = 4.8E+11  # reservoir volume, m3
porevol = totalporo * bulkvol
Pp_prod = 4  # Pp decrease due to production

"Step 6. CO2-saturated rock elastic property calculation after max year"
#fluid saturations
Sw = 0.5
SCO2 = 0.5
e = 2 #default Brie exponent

"""
Modelling
"""

"Step 0. Matrix property calculation"
#density of brine should be calculated first to fit in the equation ...

rhobrine = BW_brine_density(temp, Pp_baseline, salinity)
rhom = (rhob_log - (totalporo*rhobrine)) / (1-totalporo)
# print(rhom)

# Voigt-Reuss-Hill
mincomposition = np.array([calc, clay, dolo, qtz])
rhomineral = np.array([rhocalc, rhoclay, rhodolo, rhoqtz])
Kmineral = np.array([Kcalc, Kclay, Kdolo, Kqtz])
Kmineral_inv = 1 / Kmineral
Gmineral = np.array([Gcalc, Gclay, Gdolo, Kqtz])
Gmineral_inv = 1 / Gmineral

Kv = sum(mincomposition * Kmineral)
Kr_inv = sum(mincomposition * Kmineral_inv)
Gv = sum(mincomposition * Gmineral)
Gr_inv = sum(mincomposition * Gmineral_inv)
Km = (Kv+(1/Kr_inv)) / 2
Gm = (Gv+(1/Gr_inv)) / 2
# print(Km, Gm, rhom)

"Step 1. Brine elastic property calculation"
rhobrine = BW_brine_density(temp, Pp_baseline, salinity)
Kbrine = BW_brine_bulk(temp, Pp_baseline, salinity, rhobrine)
Gbrine = 0

"Step 2. Brine-saturated rock elastic property calculation"
stuffs1_brine = stuffs(Km, Kbrine, Gm, Gbrine, totalporo, compos1, alpha1)
stuffs2_brine = stuffs(Km, Kbrine, Gm, Gbrine, totalporo, compos2, alpha2)

ci_1_brine = stuffs1_brine[4]
ci_2_brine = stuffs2_brine[4]

PQ_1_brine = PQ(stuffs1_brine[0], stuffs1_brine[1], stuffs1_brine[2], stuffs1_brine[5], stuffs1_brine[6])
PQ_2_brine = PQ(stuffs2_brine[0], stuffs2_brine[1], stuffs2_brine[2], stuffs2_brine[5], stuffs2_brine[6])

sigma_P_brine = (ci_1_brine * PQ_1_brine[0]) + (ci_2_brine * PQ_2_brine[0])
sigma_Q_brine = (ci_1_brine * PQ_1_brine[1]) + (ci_2_brine * PQ_2_brine[1])

KusTok_brine = KusterToksoz(sigma_P_brine, sigma_Q_brine, Km, Gm, Kbrine, rhom, rhobrine)
rho_baseline = KusTok_brine[2]; Vp_baseline = KusTok_brine[3]*1000; Vs_baseline = KusTok_brine[4]*1000
print(Vp_baseline)

"Step 3. Dry rock elastic property"
rhof_dry = 0  # when dry, no fluid, all equals to zero
Kf_dry = 0
Gf_dry = 0

stuffs1_dry = stuffs(Km, Kf_dry, Gm, Gf_dry, totalporo, compos1, alpha1)
stuffs2_dry = stuffs(Km, Kf_dry, Gm, Gf_dry, totalporo, compos2, alpha2)

ci_1_dry = stuffs1_dry[4]
ci_2_dry = stuffs2_dry[4]

PQ_1_dry = PQ(stuffs1_dry[0], stuffs1_dry[1], stuffs1_dry[2], stuffs1_dry[5], stuffs1_dry[6])
PQ_2_dry = PQ(stuffs2_dry[0], stuffs2_dry[1], stuffs2_dry[2], stuffs2_dry[5], stuffs2_dry[6])

sigma_P_dry = (ci_1_dry * PQ_1_dry[0]) + (ci_2_dry * PQ_2_dry[0])
sigma_Q_dry = (ci_1_dry * PQ_1_dry[1]) + (ci_2_dry * PQ_2_dry[1])

KusTok_dry = KusterToksoz(sigma_P_dry, sigma_Q_dry, Km, Gm, Kf_dry, rhom, rhof_dry)
K_dry = KusTok_dry[0]
K_pore = totalporo / ((1 / K_dry) - ((1 - totalporo) / Km))


"Step 4. Pore pressure calculation"
compressibility = 1 / (K_pore * 1E+09)

delta_Pp = 0
Pp_post_record = []
rhoCO2_record = []
KCO2_record = []

for i in period:
    Pp_post = 34 + delta_Pp
    rhoCO2_post = EOS(Pp_post, temp)[0]
    KCO2_post = EOS(Pp_post, temp)[1]
    vol_cum = (800 * i * 365) / rhoCO2_post
    delta_Pp = (1 / porevol) * (vol_cum / compressibility) / 1E+06

    "MAGIC CODE BELOW!!! To record all properties the iterated result. x can be anything; Pp_post, delta_Pp, vol_cum, rhoCO2_post, KCO2_post"
    Pp_post_record.append(float(Pp_post))
    rhoCO2_record.append(float(rhoCO2_post))
    KCO2_record.append(float(KCO2_post))
# for plotting Pp only; inactivate if don't want to plot
# plt.plot(period, Pp_post_record, '.-')
# plt.show()

# print(compressibility)

"Step 5. CO2 property calculation after x year"

# monitor after how many years?
howmany = 10 #max year of monitoring; if others, should be integers and less than max year
period_index = np.where(period == howmany) #find what the index of the year is
Pp_post_record = np.array(Pp_post_record)
rhoCO2_record = np.array(rhoCO2_record)
KCO2_record = np.array(KCO2_record)

Pp_post = Pp_post_record[period_index] #Pp post at x year from its index
rhoCO2 = rhoCO2_record[period_index]
KCO2 = KCO2_record[period_index]

"Step 6. CO2-saturated rock elastic property calculation after max year"
# mix brine + CO2; for carbonate use brie mixing
rho_inj = ((1-Sw)*rhoCO2) + (Sw*rhobrine)
Kf_inj = ((Kbrine - KCO2)*(Sw**e)) + KCO2
Gf_inj = 0
# print(rho_inj, Kf_inj)

stuffs1_inj = stuffs(Km, Kf_inj, Gm, Gf_inj, totalporo, compos1, alpha1)
stuffs2_inj = stuffs(Km, Kf_inj, Gm, Gf_inj, totalporo, compos2, alpha2)

ci_1_inj = stuffs1_inj[4]
ci_2_inj = stuffs2_inj[4]

PQ_1_inj = PQ(stuffs1_inj[0], stuffs1_inj[1], stuffs1_inj[2], stuffs1_inj[5], stuffs1_inj[6])
PQ_2_inj = PQ(stuffs2_inj[0], stuffs2_inj[1], stuffs2_inj[2], stuffs2_inj[5], stuffs2_inj[6])

sigma_P_inj = (ci_1_inj * PQ_1_inj[0]) + (ci_2_inj * PQ_2_inj[0])
sigma_Q_inj = (ci_1_inj * PQ_1_inj[1]) + (ci_2_inj * PQ_2_inj[1])

KusTok_inj = KusterToksoz(sigma_P_inj, sigma_Q_inj, Km, Gm, Kf_inj, rhom, rho_inj)
# print(KusTok_inj)
rho_post = KusTok_inj[2]; Vp_post = KusTok_inj[3]*1000; Vs_post = KusTok_inj[4]*1000
print(rho_post); print(Vp_post); print(Vs_post)

"Step 7. Rock elastic property analysis after CO2 injection"
rho_reduce = rho_baseline - rho_post; rho_reduce_percent = (rho_reduce / rho_baseline) * 100
Vp_reduce = Vp_baseline - Vp_post; Vp_reduce_percent = (Vp_reduce / Vp_baseline) * 100
Vs_increase = Vs_post - Vs_baseline; Vs_increase_percent = (Vs_increase / Vs_baseline) * 100

#change matrix to string
rho_reduce = " ".join(str(x) for x in rho_reduce); rho_reduce_percent = " ".join(str(x) for x in rho_reduce_percent)
Vp_reduce = " ".join(str(x) for x in Vp_reduce); Vp_reduce_percent = " ".join(str(x) for x in Vp_reduce_percent)
Vs_increase = " ".join(str(x) for x in Vs_increase); Vs_increase_percent = " ".join(str(x) for x in Vs_increase_percent)

"Step 8. Modelling result and visualization"

# visualize our rock model
rockmod = rockmod(compos1, compos2)
plt.suptitle('Nuwara Carbonate Rock Model', fontsize=14, fontweight='bold') #title

# results
fig = rockmod[2]
# a= "The reduce of Vp is ",Vp_reduce,"m/s or ", Vp_reduce_percent,"%"
plt.text(0.5, -0.12, "matplotlib", horizontalalignment='center', verticalalignment='center', transform=fig.transAxes)
plt.show()
