def Ks(Kd, Km, Kf, phi):
    gamma = 1.0 - phi - Kd/Km
    return Kd + (gamma + phi)**2/(gamma/Km + phi/Kf)


def Kd(Ks, Km, Kf, phi):
    gamma = phi*(Km/Kf - 1.0)
    return (Ks*(gamma + 1.0) - Km)/(gamma - 1.0 + Ks/Km)
