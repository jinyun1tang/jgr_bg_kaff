"""python script for supporting the affinity calculation
Author: Jinyun Tang <jinyuntang@gmail.com>
"""

import numpy as np
import math
import numpy as np
from soil_info import clapp_hornberg_par, moldrup_tau, cosby_psi, cosby_Dwpsi
#soil property
#constant parameters
Npsite=3000.0
k2_0=100.0
rc=1.e-6
rp=1.e-9
Na=6.e23
Ratm=50 #atmospheric resistance, 50 s/m
#define temperature
temp=25+273.15
dzs=1.
def set_micpara(k20=100.,npsite=3000.):
    k2_0=k20
    Npsite=npsite
#assuming a spheric cell
def calc_cell_permolC(r_cell):
    #return mol number of cells per mol C
    #using allometric relationship from
    v_unit=1.e-18  #unit convertor for volume
    b_unit=1.e-15  #unit convertor for dry weight
    C_frac=0.5     #fraction of dry weight as carbon
    r_cell=1.15e-6      #cell radius
    v_cell=4./3.*math.pi*r_cell**3. #cell volume
    DW_C_cell=435*(v_cell/v_unit)**0.86*b_unit*C_frac
    DW_C_mol=DW_C_cell*Na
    return 12./DW_C_mol
def calc_Kaff_soil(pct_clay, pct_sand, fom=0.0):

    #basic cell relevant parameters

    #compute hydraulic properties
    chb,sat,psisat,ksat=clapp_hornberg_par(pct_sand, pct_clay,fom)
    nn=100
    s_sat=np.linspace(0,nn,nn)*0.01
    theta=s_sat*sat
    psi,dpsidvsm=cosby_psi(s_sat, psisat, sat, chb)
    epsi = sat-theta
    Dhk=cosby_Dwpsi(s_sat, psisat, ksat, sat, chb)

    #compute water film thickness
    ni=0
    film=np.zeros(nn)
    while (ni <nn):
        'film thickness is set to no greater than 1.e-6 m'
        film[ni]=np.fmax(np.exp(-13.65-0.857*np.log(-psi[ni]*1.e-6)),1.e-7)
        ni=ni+1
    ni=0
    taug=np.zeros(nn)
    tauw=np.zeros(nn)
    while (ni<nn):
        taug[ni],tauw[ni]=moldrup_tau(sat,epsi[ni],theta[ni],chb)
        ni=ni+1
    psi = psi*1.e-6
    return s_sat, theta, epsi, taug, tauw, film, psi,chb

def calc_diffus_O2(theta, epsi, taug, tauw):
    Dw_o2=2.4e-9*temp/298.0  #aqueous tracer diffusivity at 25
    Dg_o2=1.8e-5*(temp/273.0)**1.82 #oxygen diffusivity
    henry_o2=1.3e-3*math.exp(-1500.*(1./temp-1./298.15))
    bunsen_o2=henry_o2*temp/12.2
    #bulk aqueous molecular diffusivity
    Dwbo2=Dw_o2*tauw*theta+Dg_o2*taug*epsi/bunsen_o2
    return Dwbo2

def calc_diffus_solute(theta, tauw):
    Dw_s=1.4e-9
    Dwbs=Dw_s*tauw*theta
    return Dwbs

def calc_Kaff_O2(s_sat, theta, epsi, taug, tauw, film, DZ, Ncell, BT, alphaV, factw=1.0):
    Bdens=Ncell/Na
    rm=rc*(alphaV*Ncell)**(1/3.)
    vm=math.pi*4./3.*rm**3
    k2=k2_0*Npsite
    Nmsite=BT/Bdens*factw
    #compute the exact affinity
    fintc=Npsite*rp/(Npsite*rp+math.pi*rc) #interception rate
    #for O2, no temperature effect is considered
    Dw_o2=2.4e-9*temp/298.0  #aqueous tracer diffusivity at 25
    Dg_o2=1.8e-5*(temp/273.0)**1.82 #oxygen diffusivity
    henry_o2=1.3e-3*math.exp(-1500.*(1./temp-1./298.15))
    bunsen_o2=henry_o2*temp/12.2
    #bulk aqueous molecular diffusivity
    Dwbo2=Dw_o2*tauw*theta+Dg_o2*taug*epsi/bunsen_o2
    ksurf=Dwbo2/(DZ*DZ*0.5*dzs)
    rs_tops=1./ksurf+DZ*Ratm*bunsen_o2   #gaseous
    #using the DAMM model

    k1_o2=4.0*math.pi*Dw_o2*rc*fintc*Na
    Kaff_o2_0=k2/k1_o2   #reference affinity parameter, as number of molecules per m3, aqueous
    kapvoi=(film/(Dw_o2*rm*(rm+film))+1.0/(Dwbo2*(rm+film)))*vm/(4.0*math.pi)
    gamma_o2=1.0+k1_o2*Bdens*kapvoi/vm

    kapvoi_full=kapvoi/(Nmsite+1.e-12)+rs_tops
    gamma_o2_full=1.0+k1_o2*Bdens*kapvoi_full/vm
    k1_o2_full=k1_o2/gamma_o2_full*bunsen_o2     #gaseous delivery rate
    k1_o2=k1_o2/gamma_o2*bunsen_o2               #gaseous delivery rate
    Kaff_o2g=k2/k1_o2
    Kaff_o2g_full=k2/k1_o2_full
    kappa_tops=bunsen_o2/rs_tops
    return Kaff_o2g_full,Kaff_o2g,k2,k1_o2,k1_o2_full,kappa_tops

def calc_Kaff_SC(s_sat, theta, epsi, taug, tauw,  film, DZ, Ncell, alphaV):
    #compute affinity parameter for soluble carbon
    Bdens=Ncell/Na
    rm=rc*(alphaV*Ncell)**(1/3.)
    vm=math.pi*4./3.*rm**3
    k2=k2_0*Npsite
    #compute the exact affinity
    fintc=Npsite*rp/(Npsite*rp+math.pi*rc) #interception rate
    Dw_s=1.4e-9
    k1_s=4.0*math.pi*Dw_s*rc*fintc*Na
    Kaff_s_0=k2/k1_s   #reference affinity parameter
    kapvsi=(film/(Dw_s*rm*(rm+film))+1./(Dw_s*tauw*theta*(rm+film)+1.e-20))*vm/(4.0*math.pi)
    gamma_s=1.0+k1_s*Bdens*kapvsi/vm
    k1_s=k1_s/gamma_s
    Kaff_s=k2/k1_s
    return k1_s,Kaff_s,Kaff_s_0

def calc_Kaff_ch4(s_sat, theta, epsi, taug, tauw,  film, DZ, Ncell, BT, alphaV, factw=1.0):
    #computer affinity parameter for ch4
    Bdens=Ncell/Na
    rm=rc*(alphaV*Ncell)**(1/3.)
    vm=math.pi*4./3.*rm**3
    k2=k2_0*Npsite
    Nmsite=BT/Bdens*factw
    #compute the exact affinity
    fintc=Npsite*rp/(Npsite*rp+math.pi*rc) #interception rate
    Dw_ch4=1.5e-9*temp/298.0  #aqueous tracer diffusivity at 25
    Dg_ch4=1.9e-5*(temp/298.0)**1.82 #oxygen diffusivity
    henry_ch4=1.3e-3*math.exp(-1700.*(1./temp-1./298.15))
    bunsen_ch4=henry_ch4*temp/12.2
    Dwbch4=Dw_ch4*tauw*theta+Dg_ch4*taug*epsi/bunsen_ch4

    k1_ch4=4.0*math.pi*Dw_ch4*rc*fintc*Na
    Kaff_ch4_0=k2/k1_ch4   #reference affinity parameter
    kapvi=(film/(Dw_ch4*rm*(rm+film))+1.0/(Dwbch4*(rm+film)))*vm/(4.0*math.pi)
    gamma_ch4=1.0+k1_ch4*Bdens*kapvi/vm
    k1_ch4=k1_ch4/gamma_ch4*bunsen_ch4
    Kaff_ch4g=k2/k1_ch4

    return Kaff_ch4g, k1_ch4
