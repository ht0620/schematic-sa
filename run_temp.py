import os

tobs = 100

for it in range(1,9):
    tobs *= 2
    os.system("julia temperature.jl %d > temp-%d.dat" %(tobs, it))
