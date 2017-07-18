import os

tobs=10000
kobs=10000

for it in range(0,11):
    temp = 0.10 + float(it) / 20
    os.system("julia conf_ensemble.jl %.4f %d %d > conf-%d.dat" %(temp, tobs, kobs, it))
