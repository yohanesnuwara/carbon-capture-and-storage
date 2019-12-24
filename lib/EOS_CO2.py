"""
CO2 Acoustic Properties from Equation of State database Span and Wagner (1996)
By: Yohanes Nuwara
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# load EOS database in a form of npy file
load_eos = np.load('/content/CO2Inject/lib/eostable.npy')
load_eos = pd.DataFrame(load_eos, columns=['Temperature (deg c)', 'Pressure (Pa)', 'Density (g/cc)',
                                           'Sound speed (m/s)'])

def EOS(Pp_post, temp):

# convert temperature to kelvin
    temp = temp + 273
# normalize pressure
    Pp_rounded = round(Pp_post/5)*5 #round to nearest 5 multiples, e.g. 32 to 30, 34 to 35, 36 to 35, 38 to 50
#normalize temp
    temp_rounded = round(temp/1)*1 #round to nearest number, e.g. 273.2 to 273, and 608.6 to 609
#find the normalized temp input in the dataframe
    findtemp = load_eos.loc[load_eos['Temperature (deg c)'] == temp_rounded]

# code for interpolation
    if Pp_post < Pp_rounded:
        Pp_rounded_low = (Pp_rounded - 5) * 1E+06
        Pp_rounded = Pp_rounded * 1E+06
        Pp_post = Pp_post * 1E+06

        #find density data from dataframe
        findPp_low = findtemp.loc[findtemp['Pressure (Pa)'] == Pp_rounded_low]
        findPp_rounded = findtemp.loc[findtemp['Pressure (Pa)'] == Pp_rounded]

        #for rho CO2
        rho_rounded_low = findPp_low.iloc[0]['Density (g/cc)']
        rho_rounded = findPp_rounded.iloc[0]['Density (g/cc)']
        rhoCO2 = (((Pp_post - Pp_rounded_low)*(rho_rounded - rho_rounded_low)) /
                  (Pp_rounded - Pp_rounded_low)) + rho_rounded_low

        #for bulk CO2
        vel_rounded_low = findPp_low.iloc[0]['Sound speed (m/s)']
        vel_rounded = findPp_rounded.iloc[0]['Sound speed (m/s)']
        velCO2 = (((Pp_post - Pp_rounded_low)*(vel_rounded - vel_rounded_low)) /
                  (Pp_rounded - Pp_rounded_low)) + vel_rounded_low
        KCO2 = (rhoCO2*(velCO2**2)) / 1E+06

    elif Pp_post > Pp_rounded:
        Pp_rounded_high = (Pp_rounded + 5) * 1E+06
        Pp_rounded = Pp_rounded * 1E+06
        Pp_post = Pp_post * 1E+06

        #find density data from dataframe
        findPp_high = findtemp.loc[findtemp['Pressure (Pa)'] == Pp_rounded_high]
        findPp_rounded = findtemp.loc[findtemp['Pressure (Pa)'] == Pp_rounded]

        # for rho CO2
        rho_rounded_high = findPp_high.iloc[0]['Density (g/cc)']
        rho_rounded = findPp_rounded.iloc[0]['Density (g/cc)']
        rhoCO2 = (((Pp_post - Pp_rounded) * (rho_rounded_high - rho_rounded)) /
                  (Pp_rounded_high - Pp_rounded)) + rho_rounded

        # for bulk CO2
        vel_rounded_high = findPp_high.iloc[0]['Sound speed (m/s)']
        vel_rounded = findPp_rounded.iloc[0]['Sound speed (m/s)']
        velCO2 = (((Pp_post - Pp_rounded) * (vel_rounded_high - vel_rounded)) /
                  (Pp_rounded_high - Pp_rounded)) + vel_rounded
        KCO2 = (rhoCO2*(velCO2**2)) / 1E+06

    elif Pp_post == Pp_rounded:
        Pp_post = Pp_post * 1E+06

        #find density data from dataframe
        findPp_exact = findtemp.loc[findtemp['Pressure (Pa)'] == Pp_post]

        #for rho CO2
        rho_exact = findPp_exact.iloc[0]['Density (g/cc)']
        rhoCO2 = rho_exact

        # for bulk CO2
        vel_exact = findPp_exact.iloc[0]['Sound speed (m/s)']
        velCO2 = vel_exact
        KCO2 = (rhoCO2 * (velCO2 ** 2)) / 1E+06

    return(rhoCO2, KCO2, velCO2)

# rho = EOS(34.046, 158)
# print(rho)
