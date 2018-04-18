import os

trng = [pow(2, idx) for idx in range(7, 17)]

for it in range(10):
    tobs = trng[it]
    print(tobs)
    File = 'Temp_t_%d.csv' %(it + 7)
    os.system('julia temperature.jl %d > %s' %(tobs, File))
