import csv
from linearfit import linearfit
from scipy.interpolate import interp1d
r_sat1=[]
r_hr1=[]
with open('data/doran_1.csv', 'rb') as f:
    reader = csv.reader(f,delimiter=',')
    k=0
    for row in reader:
        #print row
        if k>5:
            r_sat1.append(float(row[0]))
            r_hr1.append(float(row[1]))
        k=k+1

r_sat2=[]
r_hr2=[]
with open('data/doran_2.csv', 'rb') as f:
    reader = csv.reader(f,delimiter=',')
    k=0
    for row in reader:
        #print row
        if k>5:
            r_sat2.append(float(row[0]))
            r_hr2.append(float(row[1]))
        k=k+1

r_sat3=[]
r_hr3=[]
with open('data/doran_3.csv', 'rb') as f:
    reader = csv.reader(f,delimiter=',')
    k=0
    for row in reader:
        #print row
        if k>5:
            r_sat3.append(float(row[0]))
            r_hr3.append(float(row[1]))
        k=k+1

r_sat4=[]
r_hr4=[]
with open('data/franz.csv', 'rb') as f:
    reader = csv.reader(f,delimiter=',')
    k=0
    for row in reader:
        #print row
        if k>0:
            r_sat4.append(float(row[0]))
            r_hr4.append(float(row[1]))
        k=k+1


o2scal=3.e-3
o2scal=1.
sand=[]
clay=[]
soil_type=[]
k1s=0
k2s=0
#data from Doran et al., 1990, table 1
with open('data/doran_soil.csv', 'rb') as f:
    reader = csv.reader(f,delimiter=',')
    k=0

    for row in reader:
        print row[0],row[1],row[2],row[4]
        if k>=1:
            sand.append(float(row[1]))
            clay.append(float(row[2]))
            if int(row[4])==1:
                k1s=k1s+1
            elif int(row[4])==2:
                k2s=k2s+1
        k=k+1

import numpy as np
k2s=k1s+k2s
pct_sand=np.array(sand)
pct_clay=np.array(clay)

kt=len(pct_sand)

#set up model parameters
DZ=0.15  #topsoil thickness, 10 cm
alphaV=80.0  #volume a cell occupies
mb=2.0   #mol of microbial C
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
df=2.52
k=0
dammfls=[]
sufls=[]
supfls=[]
kk=0
pct_claym=0.
pct_sandm=0.
while k < k1s:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(pct_clay[k], pct_sand[k])
    pct_claym=pct_claym+pct_clay[k]
    pct_sandm=pct_sandm+pct_sand[k]
    kk=kk+1
    factw=s_sat**(1./df)
    Kaff_o2g_full,Kaff_o2g,k2,k1_o2,k1_o2_full,kappa_tops=calc_Kaff_O2(s_sat, theta, epsi, taug, tauw, film, DZ, Ncello, BT, alphaV,factw)
    Kaff_o2g_full=Kaff_o2g_full*o2scal
    k1_o2_full=k1_o2_full/o2scal
    k2=k2*10.
    Vmax=BT*k2
    Ncell=10
    k1_s,Kaff_s,Kaff_s_0=calc_Kaff_SC(s_sat, theta, epsi, taug, tauw,  film, DZ, Ncell, alphaV)
    #damm model
    fO2=O2/(Kaff_o2g_full+O2)
    fS=S/(Kaff_s+S)
    dammf=fO2*fS*factw
    damfmax=np.max(dammf)
    dammf=dammf/damfmax
    #su model
    js=S*k1_s
    jo2=O2*k1_o2_full
    suf=1./(1./k2+1./jo2+1./js-1./(js+jo2))*factw
    sufmax=np.max(suf)
    suf=suf/sufmax
    #supeca model
    ms=0.0
    Km_s=2.0
    nn=len(jo2)
    ni=0
    supf=np.zeros(nn)
    while (ni<nn):
        supf[ni]=supeca_model(O2,S,BT, ms, Kaff_o2g_full[ni], Kaff_s[ni], Km_s)*factw[ni]*k2
        ni=ni+1
    supfmax=np.max(supf)
    supf=supf/supfmax
    dammfls.append(dammf)
    sufls.append(suf)
    supfls.append(supf)
    k=k+1


#compute the value for mean
def calc_mean_func(pct_claym, pct_sandm):
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(pct_claym, pct_sandm)
    factw=s_sat**(1./df)
    Kaff_o2g_full,Kaff_o2g,k2,k1_o2,k1_o2_full,kappa_tops=calc_Kaff_O2(s_sat, theta, epsi, taug, tauw, film, DZ, Ncello, BT, alphaV,factw)
    Kaff_o2g_full=Kaff_o2g_full*o2scal
    k1_o2_full=k1_o2_full/o2scal
    k2=k2*10.
    Vmax=BT*k2
    Ncell=10
    k1_s,Kaff_s,Kaff_s_0=calc_Kaff_SC(s_sat, theta, epsi, taug, tauw,  film, DZ, Ncell, alphaV)
    #damm model
    fO2=O2/(Kaff_o2g_full+O2)
    fS=S/(Kaff_s+S)
    dammf=fO2*fS*factw
    damfmax=np.max(dammf)
    dammf=dammf/damfmax
    #su model
    js=S*k1_s
    jo2=O2*k1_o2_full
    suf=1./(1./k2+1./jo2+1./js-1./(js+jo2))*factw
    sufmax=np.max(suf)
    suf=suf/sufmax
    #supeca model
    ms=0.0
    Km_s=2.0
    nn=len(jo2)
    ni=0
    supf=np.zeros(nn)
    while (ni<nn):
        supf[ni]=supeca_model(O2,S,BT, ms, Kaff_o2g_full[ni], Kaff_s[ni], Km_s)*factw[ni]*k2
        ni=ni+1
    supfmax=np.max(supf)
    supf=supf/supfmax
    return dammf, suf, supf
pct_claym=pct_claym/kk
pct_sandm=pct_sandm/kk
dammfm, sufm, supfm=calc_mean_func(pct_claym, pct_sandm)
plt_to_file=True
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
if plt_to_file:
    pdf=PdfPages('figure/topsoi1.pdf')
    fig=plt.figure()

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }
plt.rcParams["figure.figsize"] = (11,8)
if o2scal != 1.0:
    ff = open('model_tops1.txt', 'w')
meanf=dammfm
f = interp1d(s_sat, meanf)
r_hrnew=f(r_sat1)
w1=linearfit(r_hrnew,r_hr1)
print 'doran1 damm: y=x*%f+%f, R2=%f'%(w1[0],w1[2],w1[5])
if o2scal != 1.0:
    for jj in range(len(s_sat)):
        ff.write('%.2f,%.3f,%.3f,%.3f\n'%(s_sat[jj],dammfm[jj],sufm[jj],supfm[jj]))

plt.subplot(431)
plt.plot(s_sat,np.transpose(dammfls),'0.7')
plt.plot(s_sat,dammfm,'k')
plt.text(-0.04, 0.85, '(a1)',fontdict=font)
plt.plot(r_sat1,r_hr1,'bo',markersize=4,markerfacecolor='w')
plt.xlim([-0.04,1.04])
plt.xticks(np.arange(0, 1.2, step=0.2),[])
plt.yticks(np.arange(0, 1.2, step=0.25),[0.0,'',0.5,'',1.0])
plt.text(0.24, 0.05, 'Soil class 1',fontdict=font)

meanf=sufm
f = interp1d(s_sat, meanf)
r_hrnew=f(r_sat1)
w1=linearfit(r_hrnew,r_hr1)
print 'doran1 su: y=x*%f+%f, R2=%f'%(w1[0],w1[2],w1[5])
plt.subplot(432)
plt.plot(s_sat,np.transpose(sufls),'0.7')
plt.plot(s_sat,sufm,'k')
plt.plot(r_sat1,r_hr1,'bo',markersize=4,markerfacecolor='w')
plt.text(-0.04, 0.85, '(a2)',fontdict=font)
plt.xlim([-0.04,1.04])
plt.xticks(np.arange(0, 1.2, step=0.2),[])
plt.yticks(np.arange(0, 1.2, step=0.25),[])
plt.text(0.24, 0.05, 'Soil class 1',fontdict=font)
meanf=supfm
f = interp1d(s_sat, meanf)
r_hrnew=f(r_sat1)
w1=linearfit(r_hrnew,r_hr1)
print 'doran1 supeca: y=x*%f+%f, R2=%f'%(w1[0],w1[2],w1[5])
print '-----------------------------------------------------'
plt.subplot(433)
plt.plot(s_sat,np.transpose(supfls),'0.7')
plt.plot(s_sat,supfm,'k')
plt.plot(r_sat1,r_hr1,'bo',markersize=4,markerfacecolor='w')
plt.text(-0.04, 0.85, '(a3)',fontdict=font)
plt.xlim([-0.04,1.04])
plt.xticks(np.arange(0, 1.2, step=0.2),[])
plt.yticks(np.arange(0, 1.2, step=0.25),[])
plt.text(0.24, 0.05, 'Soil class 1',fontdict=font)
dammfls=[]
sufls=[]
supfls=[]
kk=0
pct_claym=0.
pct_sandm=0.
while k < k2s:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(pct_clay[k], pct_sand[k])
    pct_claym=pct_claym+pct_clay[k]
    pct_sandm=pct_sandm+pct_sand[k]
    kk=kk+1
    factw=s_sat**(1./df)
    Kaff_o2g_full,Kaff_o2g,k2,k1_o2,k1_o2_full,kappa_tops=calc_Kaff_O2(s_sat, theta, epsi, taug, tauw, film, DZ, Ncello, BT, alphaV,factw)
    Kaff_o2g_full=Kaff_o2g_full*o2scal
    k1_o2_full=k1_o2_full/o2scal
    k2=k2*10.
    Vmax=BT*k2
    Ncell=10
    k1_s,Kaff_s,Kaff_s_0=calc_Kaff_SC(s_sat, theta, epsi, taug, tauw,  film, DZ, Ncell, alphaV)
    #damm model
    fO2=O2/(Kaff_o2g_full+O2)
    fS=S/(Kaff_s+S)
    dammf=fO2*fS*factw
    damfmax=np.max(dammf)
    dammf=dammf/damfmax
    #su model
    js=S*k1_s
    jo2=O2*k1_o2_full
    suf=1./(1./k2+1./jo2+1./js-1./(js+jo2))*factw
    sufmax=np.max(suf)
    suf=suf/sufmax
    #supeca model
    Kaff_s=k2/k1_s
    ms=0.0
    Km_s=2.0
    nn=len(jo2)
    ni=0
    supf=np.zeros(nn)
    while (ni<nn):
        supf[ni]=supeca_model(O2,S,BT, ms, Kaff_o2g_full[ni], Kaff_s[ni], Km_s)*factw[ni]*k2
        ni=ni+1
    supfmax=np.max(supf)
    supf=supf/supfmax
    dammfls.append(dammf)
    sufls.append(suf)
    supfls.append(supf)
    k=k+1

pct_claym=pct_claym/kk
pct_sandm=pct_sandm/kk
dammfm, sufm, supfm=calc_mean_func(pct_claym, pct_sandm)
if o2scal != 1.0:
    for jj in range(len(s_sat)):
        ff.write('%.2f,%.3f,%.3f,%.3f\n'%(s_sat[jj],dammfm[jj],sufm[jj],supfm[jj]))

meanf=dammfm
f = interp1d(s_sat, meanf)
r_hrnew=f(r_sat2)
w1=linearfit(r_hrnew,r_hr2)
print 'doran2 damm: y=x*%f+%f, R2=%f'%(w1[0],w1[2],w1[5])
plt.subplot(434)
plt.plot(s_sat,np.transpose(dammfls),'0.7')
plt.plot(s_sat,dammfm,'k')
plt.plot(r_sat2,r_hr2,'bo',markersize=4,markerfacecolor='w')
plt.text(-0.04, 0.85, '(b1)',fontdict=font)
plt.xlim([-0.04,1.04])
plt.xticks(np.arange(0, 1.2, step=0.2),[])
plt.yticks(np.arange(0, 1.2, step=0.25),[0.0,'',0.5,'',1.0])
plt.text(0.24, 0.05, 'Soil class 2',fontdict=font)
meanf=sufm
f = interp1d(s_sat, meanf)
r_hrnew=f(r_sat2)
w1=linearfit(r_hrnew,r_hr2)
print 'doran2 su: y=x*%f+%f, R2=%f'%(w1[0],w1[2],w1[5])
plt.subplot(435)
plt.plot(s_sat,np.transpose(sufls),'0.7')
plt.plot(s_sat,sufm,'k')
plt.plot(r_sat2,r_hr2,'bo',markersize=4,markerfacecolor='w')
plt.text(-0.04, 0.85, '(b2)',fontdict=font)
plt.xlim([-0.04,1.04])
plt.xticks(np.arange(0, 1.2, step=0.2),[])
plt.yticks(np.arange(0, 1.2, step=0.25),[])
plt.text(0.24, 0.05, 'Soil class 2',fontdict=font)
meanf=supfm
f = interp1d(s_sat, meanf)
r_hrnew=f(r_sat2)
w1=linearfit(r_hrnew,r_hr2)
print 'doran2 supeca: y=x*%f+%f, R2=%f'%(w1[0],w1[2],w1[5])
print '-----------------------------------------------------'
plt.subplot(436)
plt.plot(s_sat,np.transpose(supfls),'0.7')
plt.plot(s_sat,supfm,'k')
plt.plot(r_sat2,r_hr2,'bo',markersize=4,markerfacecolor='w')
plt.text(-0.04, 0.85, '(b3)',fontdict=font)
plt.xlim([-0.04,1.04])
plt.xticks(np.arange(0, 1.2, step=0.2),[])
plt.yticks(np.arange(0, 1.2, step=0.25),[])
plt.text(0.24, 0.05, 'Soil class 2',fontdict=font)
dammfls=[]
sufls=[]
supfls=[]
kk=0
pct_claym=0.
pct_sandm=0.
while k < kt:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(pct_clay[k], pct_sand[k])
    pct_claym=pct_claym+pct_clay[k]
    pct_sandm=pct_sandm+pct_sand[k]
    kk=kk+1
    factw=s_sat**(1./df)
    Kaff_o2g_full,Kaff_o2g,k2,k1_o2,k1_o2_full,kappa_tops=calc_Kaff_O2(s_sat, theta, epsi, taug, tauw, film, DZ, Ncello, BT, alphaV,factw)
    k2=k2*10.
    Vmax=BT*k2
    Kaff_o2g_full=Kaff_o2g_full*o2scal
    k1_o2_full=k1_o2_full/o2scal
    Ncell=10
    k1_s,Kaff_s,Kaff_s_0=calc_Kaff_SC(s_sat, theta, epsi, taug, tauw,  film, DZ, Ncell, alphaV)
    #damm model
    fO2=O2/(Kaff_o2g_full+O2)
    fS=S/(Kaff_s+S)
    dammf=fO2*fS*factw
    damfmax=np.max(dammf)
    dammf=dammf/damfmax
    #su model
    js=S*k1_s
    jo2=O2*k1_o2_full
    suf=1./(1./k2+1./jo2+1./js-1./(js+jo2))*factw
    sufmax=np.max(suf)
    suf=suf/sufmax
    #supeca model
    Kaff_s=k2/k1_s
    ms=0.0
    Km_s=2.0
    nn=len(jo2)
    ni=0
    supf=np.zeros(nn)
    while (ni<nn):
        supf[ni]=supeca_model(O2,S,BT, ms, Kaff_o2g_full[ni], Kaff_s[ni], Km_s)*factw[ni]*k2
        ni=ni+1
    supfmax=np.max(supf)
    supf=supf/supfmax
    dammfls.append(dammf)
    sufls.append(suf)
    supfls.append(supf)
    k=k+1

pct_claym=pct_claym/kk
pct_sandm=pct_sandm/kk
dammfm, sufm, supfm=calc_mean_func(pct_claym, pct_sandm)
if o2scal != 1.0:
    for jj in range(len(s_sat)):
        ff.write('%.2f,%.3f,%.3f,%.3f\n'%(s_sat[jj],dammfm[jj],sufm[jj],supfm[jj]))


meanf=dammfm
f = interp1d(s_sat, meanf)
r_hrnew=f(r_sat3)
w1=linearfit(r_hrnew,r_hr3)
print 'doran3 damm: y=x*%f+%f, R2=%f'%(w1[0],w1[2],w1[5])
plt.subplot(437)
plt.plot(s_sat,np.transpose(dammfls),'0.7')
plt.plot(s_sat,dammfm,'k')
plt.plot(r_sat3,r_hr3,'bo',markersize=4,markerfacecolor='w')
plt.text(-0.04, 0.85, '(c1)',fontdict=font)
plt.xlim([-0.04,1.04])
plt.xticks(np.arange(0, 1.2, step=0.2),[])
plt.yticks(np.arange(0, 1.2, step=0.25),[0.0,'',0.5,'',1.0])
plt.text(0.24, 0.05, 'Soil class 3',fontdict=font)
meanf=sufm
f = interp1d(s_sat, meanf)
r_hrnew=f(r_sat3)
w1=linearfit(r_hrnew,r_hr3)
print 'doran3 su: y=x*%f+%f, R2=%f'%(w1[0],w1[2],w1[5])
plt.subplot(438)
plt.plot(s_sat,np.transpose(sufls),'0.7')
plt.plot(s_sat,sufm,'k')
plt.plot(r_sat3,r_hr3,'bo',markersize=4,markerfacecolor='w')
plt.text(-0.04, 0.85, '(c2)',fontdict=font)
plt.xlim([-0.04,1.04])
plt.xticks(np.arange(0, 1.2, step=0.2),[])
plt.yticks(np.arange(0, 1.2, step=0.25),[])
plt.text(0.24, 0.05, 'Soil class 3',fontdict=font)
meanf=supfm
f = interp1d(s_sat, meanf)
r_hrnew=f(r_sat3)
w1=linearfit(r_hrnew,r_hr3)
print 'doran3 supeca: y=x*%f+%f, R2=%f'%(w1[0],w1[2],w1[5])
print '-----------------------------------------------------'
plt.subplot(439)
plt.plot(s_sat,np.transpose(supfls),'0.7')
plt.plot(s_sat,supfm,'k')
plt.plot(r_sat3,r_hr3,'bo',markersize=4,markerfacecolor='w')
plt.text(-0.04, 0.85, '(c3)',fontdict=font)
plt.xlim([-0.04,1.04])
plt.xticks(np.arange(0, 1.2, step=0.2),[])
plt.yticks(np.arange(0, 1.2, step=0.25),[])
plt.text(0.24, 0.05, 'Soil class 3',fontdict=font)
sand=[]
clay=[]
soil_type=[]

with open('data/franz_soil.csv', 'rb') as f:
    reader = csv.reader(f,delimiter=',')
    k=0
    for row in reader:
        print row[0],row[1],row[2]
        if k>=1:
            sand.append(float(row[1]))
            clay.append(float(row[2]))
        k=k+1

pct_sand=np.array(sand)
pct_clay=np.array(clay)
kt=len(pct_sand)
dammfls=[]
sufls=[]
supfls=[]
k=0
kk=0
pct_claym=0.
pct_sandm=0.
while k < kt:
    s_sat, theta, epsi, taug, tauw, film,psi,chb=calc_Kaff_soil(pct_clay[k], pct_sand[k])
    pct_claym=pct_claym+pct_clay[k]
    pct_sandm=pct_sandm+pct_sand[k]
    kk=kk+1
    factw=s_sat**(1./df)
    Kaff_o2g_full,Kaff_o2g,k2,k1_o2,k1_o2_full,kappa_tops=calc_Kaff_O2(s_sat, theta, epsi, taug, tauw, film, DZ, Ncello, BT, alphaV,factw)
    Kaff_o2g_full=Kaff_o2g_full*o2scal
    k1_o2_full=k1_o2_full/o2scal
    k2=k2*10.
    Vmax=BT*k2
    Ncell=10
    k1_s,Kaff_s,Kaff_s_0=calc_Kaff_SC(s_sat, theta, epsi, taug, tauw,  film, DZ, Ncell, alphaV)
    #damm model
    fO2=O2/(Kaff_o2g_full+O2)
    fS=S/(Kaff_s+S)
    dammf=fO2*fS*factw
    damfmax=np.max(dammf)
    dammf=dammf/damfmax
    #su model
    js=S*k1_s
    jo2=O2*k1_o2_full
    suf=1./(1./k2+1./jo2+1./js-1./(js+jo2))*factw
    sufmax=np.max(suf)
    suf=suf/sufmax
    #supeca model
    Kaff_s=k2/k1_s
    ms=0.0
    Km_s=2.0
    nn=len(jo2)
    ni=0
    supf=np.zeros(nn)
    while (ni<nn):
        supf[ni]=supeca_model(O2,S,BT, ms, Kaff_o2g_full[ni], Kaff_s[ni], Km_s)*factw[ni]*k2
        ni=ni+1
    supfmax=np.max(supf)
    supf=supf/supfmax
    dammfls.append(dammf)
    sufls.append(suf)
    supfls.append(supf)
    k=k+1
pct_claym=pct_claym/kk
pct_sandm=pct_sandm/kk
dammfm, sufm, supfm=calc_mean_func(pct_claym, pct_sandm)
if o2scal != 1.0:
    for jj in range(len(s_sat)):
        ff.write('%.2f,%.3f,%.3f,%.3f\n'%(s_sat[jj],dammfm[jj],sufm[jj],supfm[jj]))
    ff.close()


plt.subplot(4,3,10)
#meanf=np.mean(dammfls,0)
meanf=dammfm
f = interp1d(s_sat, meanf)
r_hrnew=f(r_sat4)
w1=linearfit(r_hrnew,r_hr4)
print 'franz damm: y=x*%f+%f, R2=%f'%(w1[0],w1[2],w1[5])
plt.plot(s_sat,np.transpose(dammfls),'0.7')
plt.plot(s_sat,dammfm,'k')
plt.plot(r_sat4,r_hr4,'bo',markersize=4,markerfacecolor='w')
plt.text(-0.04, 0.85, '(d1)',fontdict=font)
plt.xlim([-0.04,1.04])
plt.xticks(np.arange(0, 1.2, step=0.2),[0.0,'',0.4,'',0.8,''])
plt.yticks(np.arange(0, 1.2, step=0.25),[0.0,'',0.5,'',1.0])
plt.text(0.24, 0.05, 'Soil class 4',fontdict=font)
plt.subplot(4,3,11)
#meanf=np.mean(sufls,0)
meanf=sufm
f = interp1d(s_sat, meanf)
r_hrnew=f(r_sat4)
w1=linearfit(r_hrnew,r_hr4)
print 'franz su: y=x*%f+%f, R2=%f'%(w1[0],w1[2],w1[5])
plt.plot(s_sat,np.transpose(sufls),'0.7')
plt.plot(s_sat,sufm,'k')
plt.plot(r_sat4,r_hr4,'bo',markersize=4,markerfacecolor='w')
plt.text(-0.04, 0.85, '(d2)',fontdict=font)
plt.xlim([-0.04,1.04])
plt.xticks(np.arange(0, 1.2, step=0.2),[0.0,'',0.4,'',0.8,''])
plt.yticks(np.arange(0, 1.2, step=0.25),[])
plt.text(0.24, 0.05, 'Soil class 4',fontdict=font)
plt.subplot(4,3,12)
#meanf=np.mean(supfls,0)
meanf=supfm
f = interp1d(s_sat, meanf)
r_hrnew=f(r_sat4)
w1=linearfit(r_hrnew,r_hr4)
print 'franz supeca: y=x*%f+%f, R2=%f'%(w1[0],w1[2],w1[5])
plt.plot(s_sat,np.transpose(supfls),'0.7')
plt.plot(s_sat,supfm,'k')
plt.plot(r_sat4,r_hr4,'bo',markersize=4,markerfacecolor='w')
plt.xticks(np.arange(0, 1.2, step=0.2),[0.0,'',0.4,'',0.8,''])
plt.yticks(np.arange(0, 1.2, step=0.25),[])
plt.xlim([-0.04,1.04])
plt.text(-0.04, 0.85, '(d3)',fontdict=font)
plt.text(0.24, 0.05, 'Soil class 4',fontdict=font)
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 14,
        }
plt.text(0.115, 0.95, 'DM model',fontdict=font,transform=plt.gcf().transFigure)
plt.text(0.42, 0.95, 'SU model',fontdict=font,transform=plt.gcf().transFigure)
plt.text(0.7, 0.95, 'SUPECA model',fontdict=font,transform=plt.gcf().transFigure)
plt.text(0.375, 0.015, 'Relative saturation',fontdict=font,transform=plt.gcf().transFigure)
plt.text(0.015, 0.75, 'Normalized respiration',fontdict=font,transform=plt.gcf().transFigure, rotation=90)


plt.subplots_adjust(top=0.92, bottom=0.10, left=0.10, right=0.95, hspace=0.1,
                    wspace=0.15)
if plt_to_file:
    pdf.savefig(fig)
    pdf.close()
else:
    plt.show()
