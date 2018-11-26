import csv
from DAMM_model import damm_flx
from supeca import supmic_flx
from su_model import su_flx
import numpy as np
#set up model parameters
DZ=0.15  #topsoil thickness, 10 cm
alphaV=80.0  #volume a cell occupies
mb=2.   #mol of microbial C
Ncello=10.0  #number of cells per microsite, this is for oxygen
from Kaffapp import set_micpara
set_micpara(1.e2)

from Kaffapp import rc, calc_Kaff_soil, calc_Kaff_O2, calc_cell_permolC, calc_Kaff_SC,calc_Kaff_ch4
from supeca import supeca_model
#number of cells per mol carbon
cell_molC=calc_cell_permolC(rc)
BT=mb*cell_molC  #total number of cells, in mol /m3

O2=8.57 # mol /m3
S=0.3 #mol/m3 ~ 3.6 mg/L

k=0
dammfls=[]
sufls=[]
supfls=[]
df=3.

o2scal=1.
pct_sand=20.
pct_clay=[10.,20.,40.,70.]
kt=len(pct_clay)
k=0

while k < kt:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(pct_clay[k], pct_sand)
    factw=s_sat**(1./df)
    Kaff_o2g_full,Kaff_o2g,k2,k1_o2,k1_o2_full,kappa_tops=calc_Kaff_O2(s_sat, theta, epsi, taug, tauw, film, DZ, Ncello, BT, alphaV,factw)
    Kaff_o2g_full=Kaff_o2g_full*o2scal
    k1_o2_full=k1_o2_full/o2scal
    k2=k2*10.
    Vmax=BT*k2
    Ncell=10
    k1_s,Kaff_s,Kaff_s_0=calc_Kaff_SC(s_sat, theta, epsi, taug, tauw,  film, DZ, Ncell, alphaV)
    #damm model
    dammf=damm_flx(O2,S,factw,Kaff_s*1.e0,Kaff_o2g,kappa_tops,Vmax)
    damfmax=np.max(dammf)
    dammf=dammf/damfmax
    #su model
    suf=su_flx(O2,S, factw,BT, k2*1.e0, k1_o2, k1_s,kappa_tops)
    sufmax=np.max(suf)
    suf=suf/sufmax
    #supeca model
    ms=0.0
    Km_s=2.0
    supf=supmic_flx(O2,S,factw, BT, ms, Km_s, Kaff_o2g, Kaff_s, kappa_tops, k2)
    supfmax=np.max(supf)
    supf=supf/supfmax
    dammfls.append(dammf)
    sufls.append(suf)
    supfls.append(supf)
    k=k+1

plt_to_file=True
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

if plt_to_file:
    pdf=PdfPages('figure/clay_hrsens.pdf')
    fig=plt.figure()
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

#plt.rcParams["figure.figsize"] = (11,4)
plt.subplot(221)
plt.plot(s_sat,np.transpose(dammfls))
plt.yticks(np.arange(0, 1.2, step=0.25),[0.0,'',0.5,'',1.0])
plt.text(0.05, 0.9, '(a) DM',fontdict=font)
plt.xticks(np.arange(0, 1.2, step=0.2),[])
plt.subplot(222)
plt.plot(s_sat,np.transpose(sufls))
plt.yticks(np.arange(0, 1.2, step=0.25),[])
plt.text(0.05, 0.9, '(b) SU',fontdict=font)
plt.xlabel('Relative saturation',fontdict=font)
plt.subplot(223)
for n in range(4):
    plt.plot(s_sat,np.transpose(supfls[n][:]),label="clay=%d%s"%(pct_clay[n],'%',))
plt.yticks(np.arange(0, 1.2, step=0.25),[])
plt.yticks(np.arange(0, 1.2, step=0.25),[0.0,'',0.5,'',1.0])
plt.text(0.05, 0.9, '(c) SUPECA',fontdict=font)
plt.xlabel('Relative saturation',fontdict=font)
plt.legend(loc=9, bbox_to_anchor=(1.7, 0.6))
plt.subplots_adjust(top=0.9, bottom=0.12, left=0.10, right=0.95, hspace=0.1,
                    wspace=0.15)
plt.text(0.15, 0.915, 'Clay content effect on moisture-respiration relationship',fontdict=font,transform=plt.gcf().transFigure)

plt.text(0.015, 0.75, 'Normalized respiration',fontdict=font,transform=plt.gcf().transFigure, rotation=90)
if plt_to_file:
    pdf.savefig(fig)
    pdf.close()
else:
    plt.show()
