"""python script for supeca model
Author: Jinyun Tang <jinyuntang@gmail.com>
"""
import numpy as np

def su_model(o2a,S, par, factw):
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
    jS=par.k1_S*S
    while(True):
        jO2=par.k1_o2*o2gc
        Flx=par.BT/(1./par.k2+1./jO2+1./jS-1./(jO2+jS))*factw
        o2gc1=par.kappa_s*o2a/(par.kappa_s+Flx/o2gc)
        if(abs(o2gc1-o2gc)<1.e-10):
            break

        o2gc=o2gc1
        iter=iter+1

    return Flx

class Para_damm:
    k1_o2=0.
    k1_S=0.
    kappa_s=0.
    BT=0.
    k2=0.

def su_flx(O2,S, factw,BT, k2, k1_o2, k1_s,kappa_tops):
    nn=len(factw)
    ni=0
    Flx=np.zeros(nn)

    par=Para_damm()

    while(ni<nn):
        par.k1_o2=k1_o2[ni]
        par.k1_S=k1_s[ni]
        par.kappa_s=kappa_tops[ni]
        par.BT=BT
        par.k2=k2
        Flx[ni]=su_model(O2,S,par, factw[ni])
        ni=ni+1
    return Flx
