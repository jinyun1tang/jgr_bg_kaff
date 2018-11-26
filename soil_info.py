"""python script for supporting the skopp test
Author: Jinyun Tang <jinyuntang@gmail.com>
"""

import numpy as np

def clapp_hornberg_par(pct_sand, pct_clay, fom=0.):
    """Calculate Clapp Honberg parameters
    chb: b constant
    sat: soil porosity
    psisat: mm
    ksat: mm/s
    """
    chb=2.91+0.195*pct_clay
    sat=0.489-0.00126*pct_sand
    psisat=-10.0**(1.88-0.0131*pct_sand)
    ksat=0.0070556*10.0**(-0.884+0.0153*pct_sand)
    if fom > 0.:
        chb=chb*(1.-fom)+fom*2.7
        psisat=psisat*(1.-fom)+fom*(-10.3)
        sat=sat*(1.-fom)+fom*0.9
    return [chb, sat,psisat,ksat]


def moldrup_tau(por,epsi,theta,kappa):
    "calculate tortuosity using Moldrup's (2003) equation"
    'equation (5)'
    taug=epsi*(epsi/por)**(3.0/kappa)
    'equation (3)'
    tauw=theta*((theta+1.e-20)/por)**(kappa/3.0-1.0)
    return [taug,tauw]

def cosby_psi(s_sat, psisat, sat, chb):
    """compute the soil water potential
    psisat: mm
    s_sat: level of saturation
    chb: b constant
    psi: Pa
    dpsidvsm: Pa
    """
    psi=np.fmax(psisat*(s_sat+1.e-10)**(-chb),-1.e8)*1.e-3*1.e4
    dpsidvsm=-chb*psi/(s_sat*sat+1.e-20)
    return [psi,dpsidvsm]

def cosby_hk(s_sat, ksat, chb):
    """compute hydraulic conductivity
    s_sat: level of saturation
    chb: b constant
    ksat: saturated hydraulic conductivity, mm/s
    hk: mm/s
    """
    hk=ksat*(s_sat)**(2.0*chb+3.0)
    return hk

def cosby_Dwpsi(s_sat, psisat, ksat, sat, chb):
    """compute the hydraulic diffusivity
    s_sat: level of saturation
    chb: b constant
    ksat: saturated hydraulic conductivity
    Dwpsi: m2/s
    """
    psi, dpsidvsm=cosby_psi(s_sat, psisat, sat, chb)
    hk=cosby_hk(s_sat, ksat, chb)
    Dwpsi=hk*1.e-3*dpsidvsm*1.e-4
    return Dwpsi
