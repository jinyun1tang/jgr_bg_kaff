

import numpy as np
import matplotlib.pyplot as plt
import math
from Kaffapp import rc, calc_Kaff_soil, calc_Kaff_O2, calc_cell_permolC, calc_Kaff_SC,calc_Kaff_ch4

#check clay content
pct_sand=20
pct_clay=[10,20,40,70]
fom=[0.01,0.1,0.3,0.5]

#determine hydraulic property
DZ=0.1   #topsoil thickness, 10 cm
alphaV=80.0  #volume a cell occupies
mb=2.0   #mol of microbial C
Ncell=[4.0,10.0,100.0,200.]  #number of cells per microsite
from Kaffapp import set_micpara
set_micpara(1.e2)

#number of cells per mol carbon
cell_molC=calc_cell_permolC()
BT=mb*cell_molC  #total number of cells, in mol /m3

plt_to_file=True

from matplotlib.backends.backend_pdf import PdfPages
if plt_to_file:
    pdf=PdfPages('figure/Figure3.pdf')
    fig=plt.figure()

plt.subplot(221)
#compute soil parameters
Kaff_o2g_full=[]
Kaff_o2g=[]
for clay in pct_clay:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(clay, pct_sand)
    Kaff_o2g_fulli,Kaff_o2gi,k2,k1_o2,k1_o2_full,kappa_tops=calc_Kaff_O2(s_sat, theta, epsi, taug, tauw, film, DZ, Ncell[1], BT, alphaV)
    Kaff_o2g_full.append(Kaff_o2g_fulli)
    Kaff_o2g.append(Kaff_o2gi)

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
plt.subplot(221)
for n in range(4):
    plt.plot(s_sat,np.transpose(Kaff_o2g[n][:]),label="%d%s"%(pct_clay[n],'%',))

leg = plt.legend(loc='best', ncol=1)
leg.get_frame().set_alpha(0.)

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.3, 0.95, r'Effective O$_2$ affinity (mol m$^{-3}$)',fontdict=font,transform=plt.gcf().transFigure)
plt.subplot(222)
#compute soil parameters
Kaff_o2g_full=[]
Kaff_o2g=[]
for om in fom:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(pct_clay[1], pct_sand,om)
    Kaff_o2g_fulli,Kaff_o2gi,k2,k1_o2,k1_o2_full,kappa_tops=calc_Kaff_O2(s_sat, theta, epsi, taug, tauw, film, DZ, Ncell[1], BT, alphaV)
    Kaff_o2g_full.append(Kaff_o2g_fulli)
    Kaff_o2g.append(Kaff_o2gi)
for n in range(4):
    plt.plot(s_sat,np.transpose(Kaff_o2g[n][:]),label="%d%s"%(fom[n]*100,'%',))
leg = plt.legend(loc='best', ncol=1)
leg.get_frame().set_alpha(0.)
plt.xlabel('Relative saturation')
plt.subplot(223)
#compute soil parameters
Kaff_o2g_full=[]
Kaff_o2g=[]
for ncell in Ncell:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(pct_clay[1], pct_sand)
    Kaff_o2g_fulli,Kaff_o2gi,k2,k1_o2,k1_o2_full,kappa_tops=calc_Kaff_O2(s_sat, theta, epsi, taug, tauw, film, DZ, ncell, BT, alphaV)
    Kaff_o2g_full.append(Kaff_o2g_fulli)
    Kaff_o2g.append(Kaff_o2gi)

for n in range(4):
    plt.plot(s_sat,np.transpose(Kaff_o2g[n][:]),label="%d"%(Ncell[n],))

leg = plt.legend(loc='best', ncol=1)
leg.get_frame().set_alpha(0.)


plt.xlabel('Relative saturation')

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

plt.text(0.35, 0.6, '(a) clay',fontdict=font,transform=plt.gcf().transFigure)
plt.text(0.6, 0.875, '(b) FOM',fontdict=font,transform=plt.gcf().transFigure)
plt.text(0.12, 0.15, r'(c) N$_{cell}$',fontdict=font,transform=plt.gcf().transFigure)


plt.subplots_adjust(top=0.92, bottom=0.10, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.35)

if plt_to_file:
    pdf.savefig(fig)
    pdf.close()

    pdf=PdfPages('figure/Figure2.pdf')
    fig=plt.figure()
else:
    plt.show()
#===============================================================================
plt.subplot(221)
#compute soil parameters
Kaff_ch4g=[]
for clay in pct_clay:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(clay, pct_sand)
    Kaff_ch4gi, k1_ch4=calc_Kaff_ch4(s_sat, theta, epsi, taug, tauw,  film, DZ, Ncell[1], BT, alphaV)
    Kaff_ch4g.append(Kaff_ch4gi)

for n in range(4):
    plt.plot(s_sat,np.transpose(Kaff_ch4g[n][:]),label="%d%s"%(pct_clay[n],'%',))
leg = plt.legend(loc='best', ncol=1)
leg.get_frame().set_alpha(0.)

plt.text(0.3, 0.95, r'Effective CH$_4$ affinity (mol m$^{-3}$)',fontdict=font,transform=plt.gcf().transFigure)

plt.subplot(222)
Kaff_ch4g=[]
for om in fom:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(pct_clay[1], pct_sand,om)
    Kaff_ch4gi, k1_ch4=calc_Kaff_ch4(s_sat, theta, epsi, taug, tauw,  film, DZ, Ncell[1], BT, alphaV)
    Kaff_ch4g.append(Kaff_ch4gi)
for n in range(4):
    plt.plot(s_sat,np.transpose(Kaff_ch4g[n][:]),label="%d%s"%(fom[n]*100,'%',))

leg = plt.legend(loc='best', ncol=1)
leg.get_frame().set_alpha(0.)
plt.xlabel('Relative saturation')
plt.subplot(223)
Kaff_ch4g=[]
for ncell in Ncell:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(pct_clay[1], pct_sand)
    Kaff_ch4gi, k1_ch4=calc_Kaff_ch4(s_sat, theta, epsi, taug, tauw,  film, DZ, ncell, BT, alphaV)
    Kaff_ch4g.append(Kaff_ch4gi)

for n in range(4):
    plt.plot(s_sat,np.transpose(Kaff_ch4g[n][:]),label="%d"%(Ncell[n],))

leg = plt.legend(loc='best', ncol=1)
leg.get_frame().set_alpha(0.)

plt.xlabel('Relative saturation')

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

plt.text(0.35, 0.6, '(a) clay',fontdict=font,transform=plt.gcf().transFigure)
plt.text(0.6, 0.875, '(b) FOM',fontdict=font,transform=plt.gcf().transFigure)
plt.text(0.12, 0.15, r'(c) N$_{cell}$',fontdict=font,transform=plt.gcf().transFigure)


plt.subplots_adjust(top=0.92, bottom=0.10, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.35)

if plt_to_file:
    pdf.savefig(fig)
    pdf.close()

    pdf=PdfPages('figure/Figure4.pdf')
    fig=plt.figure()
else:
    plt.show()


#================================================
#solute
Kaff_s=[]
for clay in pct_clay:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(clay, pct_sand)
    k1_s,Kaff_si,Kaff_s_0=calc_Kaff_SC(s_sat, theta, epsi, taug, tauw,  film, DZ, Ncell[1], alphaV)
    Kaff_s.append(Kaff_si)
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
plt.subplot(221)
for n in range(4):
    plt.semilogy(s_sat,np.transpose(Kaff_s[n][:]),label="%d%s"%(pct_clay[n],'%',))
leg = plt.legend(loc='best', ncol=1)
leg.get_frame().set_alpha(0.)

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.3, 0.95, r'Effective solute affinity (mol m$^{-3}$)',fontdict=font,transform=plt.gcf().transFigure)
plt.subplot(222)
Kaff_s=[]
for om in fom:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(pct_clay[1], pct_sand,om)
    k1_s,Kaff_si,Kaff_s_0=calc_Kaff_SC(s_sat, theta, epsi, taug, tauw,  film, DZ, Ncell[1], alphaV)
    Kaff_s.append(Kaff_si)

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
for n in range(4):
    plt.semilogy(s_sat,np.transpose(Kaff_s[n][:]),label="%d%s"%(fom[n]*100,'%',))
leg = plt.legend(loc='best', ncol=1)
leg.get_frame().set_alpha(0.)
plt.xlabel('Relative saturation')
plt.subplot(223)

Kaff_s=[]
for ncell in Ncell:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(pct_clay[1], pct_sand)
    k1_s,Kaff_si,Kaff_s_0=calc_Kaff_SC(s_sat, theta, epsi, taug, tauw,  film, DZ, ncell, alphaV)
    Kaff_s.append(Kaff_si)

Kaff_bm=[]

for n in range(4):
    plt.semilogy(s_sat,np.transpose(Kaff_s[n][:]),label="%d"%(Ncell[n],))
leg = plt.legend(loc='best', ncol=1)
leg.get_frame().set_alpha(0.)
plt.xlabel('Relative saturation')

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }
plt.text(0.12, 0.586, '(a) clay',fontdict=font,transform=plt.gcf().transFigure)
plt.text(0.6, 0.586, '(b) FOM',fontdict=font,transform=plt.gcf().transFigure)
plt.text(0.12, 0.15, r'(c) N$_{cell}$',fontdict=font,transform=plt.gcf().transFigure)


plt.subplots_adjust(top=0.92, bottom=0.10, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.35)

if plt_to_file:
    pdf.savefig(fig)
    pdf.close()
else:
    plt.show()
