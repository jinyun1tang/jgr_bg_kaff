import os




py_scripts=['o2_ch4_sol_varKaff_plot.py','doranfranz_soil_topsoil_sol2_df252.py', \
	'doranfranz_soil_topsoil_sol2_df3.py','clay_hr_sens_moisture.py',\
	'varialMO_hrsens_test.py','BT_hr_sens_moisture.py','S_hr_sens_moisture.py']

os.system('mkdir -p figure')

for pys in py_scripts:
	os.system('python '+pys)
