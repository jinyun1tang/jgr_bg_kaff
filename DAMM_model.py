"""python script for supeca model
Author: Jinyun Tang <jinyuntang@gmail.com>
"""
import numpy as np

def damm_model(o2a,S,par, factw):
    """
    DAMM model based on
    o2a: atmospheric o2, mol/m3
    S: soluble carbon, mol/m3
    par.Rmax: maximum respiration rate
    par.KO2: O2 affinity, mol /m3
    par.KS: soluble carbon affinity, mol /m3
    """
    o2gc=o2a
    dt=1800.0
    iter=0
    fS=S/(par.KS+S)
    while(True):
        fO2=o2gc/(par.KO2+o2gc)
        Flx=par.Rmax*fO2*fS*factw
        o2gc1=par.kappa_s*o2a/(par.kappa_s+Flx/o2gc)
        if(abs(o2gc1-o2gc)<1.e-10):
            break

        o2gc=o2gc1
        iter=iter+1
        #if(iter>200)break
    return Flx

class Para_damm:
    KO2=0.
    KS=0.
    kappa_s=0.
    Rmax=0.


def damm_flx(O2,S,factw,Kaff_s,Kaff_o2g,kappa_tops,Vmax):
    nn=len(factw)
    Flx=np.zeros(nn)
    par=Para_damm()
    ni=0
    while(ni<nn):
        par.KO2=Kaff_o2g[ni]
        par.KS=Kaff_s[ni]
        par.kappa_s=kappa_tops[ni]
        par.Rmax=Vmax
        Flx[ni]=damm_model(O2,S,par, factw[ni])
        ni=ni+1
    return Flx
