#! /bin/bash
rm results/*
mpirun -np 8 nth -p 4 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e13 | tee M1Philippe1e13.out
cp -r results results1e13
rm results/*
mpirun -np 8 nth -p 4 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e12 | tee M1Philippe1e12.out
cp -r results results1e12
rm results/*
mpirun -np 8 nth -p 4 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e11 | tee M1Philippe1e11.out
cp -r results results1e11
rm results/*
mpirun -np 8 nth -p 4 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e10 | tee M1Philippe1e10.out
cp -r results results1e10
rm results/*
mpirun -np 8 nth -p 4 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e9 | tee M1Philippe1e9.out
cp -r results results1e9
rm results/*
mpirun -np 8 nth -p 4 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e8 | tee M1Philippe1e8.out
cp -r results results1e8
rm results/*
mpirun -np 8 nth -p 4 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 5e7 | tee M1Philippe5e7.out
cp -r results results5e7
