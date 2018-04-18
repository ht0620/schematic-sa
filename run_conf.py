import os

tobs=100000000
Trng = [float(idx) / 20 for idx in range(2, 13)]

for Temp in Trng:
    print(Temp)
    File = 'Conf_T_%.2f.csv' %Temp
    os.system('julia conf_ensemble.jl %.4f %d > %s' %(Temp, tobs, File))
