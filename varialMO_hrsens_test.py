from DAMM_model import damm_flx
import csv
from supeca import supmic_flx
from su_model import su_flx
import numpy as np
#set up model parameters
DZ=0.15  #topsoil thickness, 10 cm
alphaV=80.0  #volume a cell occupies
mb=2.   #mol of microbial C
Ncello=[10.0,20.,60.0,120.]  #number of cells per microsite, this is for oxygen
from Kaffapp import set_micpara
set_micpara(1.e2)

from Kaffapp import rc, calc_Kaff_soil, calc_Kaff_O2, calc_cell_permolC, calc_Kaff_SC,calc_Kaff_ch4
from supeca import supeca_model
#number of cells per mol carbon
cell_molC=calc_cell_permolC(rc)
BT=mb*cell_molC  #total number of cells, in mol /m3

O2=8.57 # mol /m3
S=[0.1,0.3,0.6,1.0] #mol/m3 ~ 3.6 mg/L

k=0
supflsMs=[]
supflsMsO=[]
df=2.52

o2scal=1.
pct_sand=65.
pct_clay=19.
kt=len(S)
k=0

k=0
ms=[0.0,50.,100.,200.]
supfmax0=0.
while k < kt:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(pct_clay, pct_sand)
    factw=s_sat**(1./df)
    Kaff_o2g_full,Kaff_o2g,k2,k1_o2,k1_o2_full,kappa_tops=calc_Kaff_O2(s_sat, theta, epsi, taug, tauw, film, DZ, Ncello[0], BT, alphaV,factw)

    Ncell=10
    k1_s,Kaff_s,Kaff_s_0=calc_Kaff_SC(s_sat, theta, epsi, taug, tauw,  film, DZ, Ncell, alphaV)
    #supeca model
    k2=k2*10.
    Km_s=2.0
    supf=supmic_flx(O2,S[2],factw, BT, ms[k], Km_s, Kaff_o2g, Kaff_s, kappa_tops, k2)
    supfmax=np.max(supf)
    if k==0:
        supfmax0=supfmax
    supf=supf/supfmax0
    supflsMs.append(supf)
    Kaff_o2g=Kaff_o2g*10.
    supf=supmic_flx(O2,S[2],factw, BT, ms[k], Km_s, Kaff_o2g, Kaff_s, kappa_tops, k2)
    supfmax=np.max(supf)
    supf=supf/supfmax0
    supflsMsO.append(supf)
    k=k+1


import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }


plt_to_file=True
if plt_to_file:
    pdf=PdfPages('figure/varMO_hrsens.pdf')
    fig=plt.figure()

plt.subplot(121)
for n in range(4):
    plt.plot(s_sat,np.transpose(supflsMs[n][:]),label=r'M=%.0f mol m$^{-3}$'%(ms[n],))
leg = plt.legend(loc='best', ncol=1)
leg.get_frame().set_alpha(0.1)
plt.yticks(np.arange(0, 1.2, step=0.25),[0.0,'',0.5,'',1.0])
plt.text(0.05, 0.9, '(a) ',fontdict=font)
plt.xticks(np.arange(0, 1.2, step=0.2))
plt.xlabel('Relative saturation',fontdict=font)
plt.subplot(122)
for n in range(4):
    plt.plot(s_sat,np.transpose(supflsMsO[n][:]),label=r'M=%.0f mol m$^{-3}$'%(ms[n],))
leg = plt.legend(loc='best', ncol=1)
leg.get_frame().set_alpha(0.1)
plt.yticks(np.arange(0, 1.2, step=0.25),[])
plt.xticks(np.arange(0, 1.2, step=0.2))
plt.text(0.05, 0.9, '(b) ',fontdict=font)
plt.xlabel('Relative saturation',fontdict=font)
plt.text(0.035, 0.75, 'Normalized respiration',fontdict=font,transform=plt.gcf().transFigure, rotation=90)
plt.text(0.15, 0.82, 'Mineral sorption effect on moisture-respiration relationship',fontdict=font,transform=plt.gcf().transFigure)
plt.subplots_adjust(top=0.80, bottom=0.35, left=0.12, right=0.95, hspace=0.15,
                    wspace=0.15)
if plt_to_file:
    pdf.savefig(fig)
    pdf.close()
else:
    plt.show()
