#! /bin/bash
mkdir results
# problem 8 proving the AWBS model diffusive limit.
# This run corresponds to the Zbar = 1 run of "python M1_f1integrate.py".
#mpirun -np 8 nth -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 10001 -Tmin 10000 -Tgrad 1 -a0 1e15 -Z 4
# 1) linearity in mean free path
# 2) nonlinearity in temperature equal to 5/2.
# 3) linearity in tmeperature gradient
# Reference solution.
# Limiting case of working explicit Efield. Higher Kn destroys heat flux.
#mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 6e14

rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 1.395e12 | tee M1Philippe_results/M1Philippe_p51DKn0.00001.out
cp -r results M1Philippe_results/resultsp51DKn0.00001
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 1.395e11 | tee M1Philippe_results/M1Philippe_p51DKn0.0001.out
cp -r results M1Philippe_results/resultsp51DKn0.0001
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 1.395e10 | tee M1Philippe_results/M1Philippe_p51DKn0.001.out
cp -r results M1Philippe_results/resultsp51DKn0.001
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 1.395e9 | tee M1Philippe_results/M1Philippe_p51DKn0.01.out
cp -r results M1Philippe_results/resultsp51DKn0.01
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.6975e9 | tee M1Philippe_results/M1Philippe_p51DKn0.02.out
cp -r results M1Philippe_results/resultsp51DKn0.02
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.465e9 | tee M1Philippe_results/M1Philippe_p51DKn0.03.out
cp -r results M1Philippe_results/resultsp51DKn0.03
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.34875e9 | tee M1Philippe_results/M1Philippe_p51DKn0.04.out
cp -r results M1Philippe_results/resultsp51DKn0.04
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.279e9 | tee M1Philippe_results/M1Philippe_p51DKn0.05.out
cp -r results M1Philippe_results/resultsp51DKn0.05
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.2325e9 | tee M1Philippe_results/M1Philippe_p51DKn0.06.out
cp -r results M1Philippe_results/resultsp51DKn0.06
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.1992857142857143e9 | tee M1Philippe_results/M1Philippe_p51DKn0.07.out
cp -r results M1Philippe_results/resultsp51DKn0.07
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.174375e9 | tee M1Philippe_results/M1Philippe_p51DKn0.08.out
cp -r results M1Philippe_results/resultsp51DKn0.08
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.155e9 | tee M1Philippe_results/M1Philippe_p51DKn0.09.out
cp -r results M1Philippe_results/resultsp51DKn0.09

rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 1.395e8 | tee M1Philippe_results/M1Philippe_p51DKn0.1.out
cp -r results M1Philippe_results/resultsp51DKn0.1
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.6975e8 | tee M1Philippe_results/M1Philippe_p51DKn0.2.out
cp -r results M1Philippe_results/resultsp51DKn0.2
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.465e8 | tee M1Philippe_results/M1Philippe_p51DKn0.3.out
cp -r results M1Philippe_results/resultsp51DKn0.3
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.34875e8 | tee M1Philippe_results/M1Philippe_p51DKn0.4.out
cp -r results M1Philippe_results/resultsp51DKn0.4
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.279e8 | tee M1Philippe_results/M1Philippe_p51DKn0.5.out
cp -r results M1Philippe_results/resultsp51DKn0.5
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.2325e8 | tee M1Philippe_results/M1Philippe_p51DKn0.6.out
cp -r results M1Philippe_results/resultsp51DKn0.6
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.1992857142857143e8 | tee M1Philippe_results/M1Philippe_p51DKn0.7.out
cp -r results M1Philippe_results/resultsp51DKn0.7
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.174375e8 | tee M1Philippe_results/M1Philippe_p51DKn0.8.out
cp -r results M1Philippe_results/resultsp51DKn0.8
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 0.155e8 | tee M1Philippe_results/M1Philippe_p51DKn0.9.out
cp -r results M1Philippe_results/resultsp51DKn0.9
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 1.395e7 | tee M1Philippe_results/M1Philippe_p51DKn1.0.out
cp -r results M1Philippe_results/resultsp51DKn1.0

rm results/*
mpirun -np 8 nth -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tgrad 1 -Tmax 10001 -Tmin 10000 -a0 1e16 | tee M1Philippe_results/M1Philippe_p81Dreference.out
cp -r results M1Philippe_results/resultsp81DTreference
# Solution for 10x longer mean fee path, i.e. 10x smaller cross section.
rm results/*
mpirun -np 8 nth -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tgrad 1 -Tmax 10001 -Tmin 10000 -a0 1e15 | tee M1Philippe_results/M1Philippe_p81Dmfpx10.out
cp -r results M1Philippe_results/resultsp81Dmfpx10
# Solution for 2.511x temperature resulting in (2.511^2.5 = 10)x higher flux. 
rm results/*
mpirun -np 8 nth -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tgrad 1 -Tmax 25111 -Tmin 25110 -a0 1e16 | tee M1Philippe_results/M1Philippe_p81DTx2.511.out
cp -r results M1Philippe_results/resultsp81DTx2.511
# Solution for 10x higher temperature gradient.
rm results/*
mpirun -np 8 nth -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tgrad 1 -Tmax 10010 -Tmin 10000 -a0 1e16 | tee M1Philippe_results/M1Philippe_p81DgradTx10.out
cp -r results M1Philippe_results/resultsp81DgradTx10
