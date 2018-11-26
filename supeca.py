"""python script for supeca model
Author: Jinyun Tang <jinyuntang@gmail.com>
"""
import numpy as np

def supeca_model(o2a,S,mb, ms, Kb_o2, Kb_s, Km_s):

    ss=[o2a,S]
    ee=[mb,ms]

    Fra=ee[0]/Kb_s+ee[1]/Km_s
    Fca=ss[0]/Kb_o2

    Frb=ee[0]/Kb_s
    Fcb=ss[1]/Kb_s
    Fcabk=Fca+Fcb

    Gaik=Fca+Fra
    Gbjk=Fcb+Frb
    Gabijk=Gaik+Gbjk
    cplx=ee[0]*ss[0]*ss[1]/(Kb_o2*Kb_s)/(Gaik*Gbjk*Fcabk/Gabijk+Fcabk- \
        (Fca*Gbjk+Gaik*Fcb-Gaik*Gbjk)/Gabijk)
    return cplx


def supmic(o2a,ssc,mb,ms,factw,par):
    '''
    steady state microbial model based on supeca kinetics
    iterate to steady state"
    o2a: atmospheric o2, mol /m3
    ssc: total disolvable carbon, mol C/m3
    mb: total microbial biomass, mol C/m3
    ms: mineral surface, mol C/m3
    factw: volumetric moisture content
    par.kappa_s: atmospheric and soil resistance,
    '''
    o2gc=o2a
    dt=1800.0
    iter=0
    while(True):
        "obtain the microbe-substrate complex"
        mb1=mb#*factw
        Flx=supeca_model(o2gc,ssc,mb1, ms, par.Kb_o2, par.Kb_s, par.Km_s)*par.k2*factw
        o2gc1=par.kappa_s*o2a/(par.kappa_s+Flx/o2gc)
        if(abs(o2gc1-o2gc)<1.e-10):
            break
        o2gc=o2gc1
        iter=iter+1


    return Flx



class Para_supmic:
    Kb_o2=0.
    Kb_s=0.
    Km_s=0.0
    kappa_s=0.
    k2=0.


def supmic_flx(O2,S,factw, BT, ms, Km_s, Kaff_o2g, Kaff_s, kappa_tops, k2):
    nn=len(factw)
    ni=0
    Flx=np.zeros(nn)

    par=Para_supmic()

    while(ni<nn):
        par.Kb_o2=Kaff_o2g[ni]
        par.Kb_s=Kaff_s[ni]
        par.Km_s = Km_s
        par.kappa_s=kappa_tops[ni]
        par.k2=k2
        Flx[ni]=supmic(O2,S,BT,ms,factw[ni],par)
        ni=ni+1
    return Flx
