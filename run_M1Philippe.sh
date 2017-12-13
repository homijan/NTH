#! /bin/bash
# problem 4 - step in density and constant temperature
rm results/*
mpirun -np 8 nth -p 4 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e6 | tee M1Philippe_p41D1e6.out
cp -r results M1Philippe_results/resultsp41D1e6
rm results/*
mpirun -np 8 nth -p 4 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e5 | tee M1Philippe_p41D1e5.out
cp -r results M1Philippe_results/resultsp41D1e5
rm results/*
mpirun -np 8 nth -p 4 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e4 | tee M1Philippe_p41D1e4.out
cp -r results M1Philippe_results/resultsp41D1e4
rm results/*
mpirun -np 8 nth -p 4 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e3 | tee M1Philippe_p41D1e3.out
cp -r results M1Philippe_results/resultsp41D1e3
rm results/*
mpirun -np 8 nth -p 4 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e2 | tee M1Philippe_p41D1e2.out
cp -r results M1Philippe_results/resultsp41D1e2
# problem 5 - constant density and step in temperature
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e14 | tee M1Philippe_p51D1e14.out
cp -r results M1Philippe_results/resultsp51D1e14
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e13 | tee M1Philippe_p51D1e13.out
cp -r results M1Philippe_results/resultsp51D1e13
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e12 | tee M1Philippe_p51D1e12.out
cp -r results M1Philippe_results/resultsp51D1e12
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e11 | tee M1Philippe_p51D1e11.out
cp -r results M1Philippe_results/resultsp51D1e11
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e10 | tee M1Philippe_p51D1e10.out
cp -r results M1Philippe_results/resultsp51D1e10
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e9 | tee M1Philippe_p51D1e9.out
cp -r results M1Philippe_results/resultsp51D1e9
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e8 | tee M1Philippe_p51D1e8.out
cp -r results M1Philippe_results/resultsp51D1e8
rm results/*
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 5e7 | tee M1Philippe_p51D5e7.out
cp -r results M1Philippe_results/resultsp51D5e7
# problem 6 - step in density and step in temperature as at critical density
rm results/*
mpirun -np 8 nth -p 6 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e14 | tee M1Philippe_p61D1e14.out
cp -r results M1Philippe_results/resultsp61D1e14
rm results/*
mpirun -np 8 nth -p 6 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e13 | tee M1Philippe_p61D1e13.out
cp -r results M1Philippe_results/resultsp61D1e13
rm results/*
mpirun -np 8 nth -p 6 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e12 | tee M1Philippe_p61D1e12.out
cp -r results M1Philippe_results/resultsp61D1e12
rm results/*
mpirun -np 8 nth -p 6 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e11 | tee M1Philippe_p61D1e11.out
cp -r results M1Philippe_results/resultsp61D1e11
rm results/*
mpirun -np 8 nth -p 6 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e10 | tee M1Philippe_p61D1e10.out
cp -r results M1Philippe_results/resultsp61D1e10
rm results/*
mpirun -np 8 nth -p 6 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e9 | tee M1Philippe_p61D1e9.out
cp -r results M1Philippe_results/resultsp61D1e9
rm results/*
mpirun -np 8 nth -p 6 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e8 | tee M1Philippe_p61D1e8.out
cp -r results M1Philippe_results/resultsp61D1e8
rm results/*
mpirun -np 8 nth -p 6 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 5e7 | tee M1Philippe_p61D5e7.out
cp -r results M1Philippe_results/resultsp61D5e7
# problem 7 - step in density and step in temperature as at shock
rm results/*
mpirun -np 8 nth -p 7 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e14 | tee M1Philippe_p71D1e14.out
cp -r results M1Philippe_results/resultsp71D1e14
rm results/*
mpirun -np 8 nth -p 7 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e13 | tee M1Philippe_p71D1e13.out
cp -r results M1Philippe_results/resultsp71D1e13
rm results/*
mpirun -np 8 nth -p 7 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e12 | tee M1Philippe_p71D1e12.out
cp -r results M1Philippe_results/resultsp71D1e12
rm results/*
mpirun -np 8 nth -p 7 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e11 | tee M1Philippe_p71D1e11.out
cp -r results M1Philippe_results/resultsp71D1e11
rm results/*
mpirun -np 8 nth -p 7 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e10 | tee M1Philippe_p71D1e10.out
cp -r results M1Philippe_results/resultsp71D1e10
rm results/*
mpirun -np 8 nth -p 7 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e9 | tee M1Philippe_p71D1e9.out
cp -r results M1Philippe_results/resultsp71D1e9
rm results/*
mpirun -np 8 nth -p 7 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e8 | tee M1Philippe_p71D1e8.out
cp -r results M1Philippe_results/resultsp71D1e8
rm results/*
mpirun -np 8 nth -p 7 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 5e7 | tee M1Philippe_p71D5e7.out
cp -r results M1Philippe_results/resultsp71D5e7
# problem 6 - step in density and step in temperature in 2D as at shock
rm results/*
mpirun -np 8 nth -p 7 -m data/square01_quad.mesh -rs 4 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e12 | tee M1Philippe_p72D1e12.out
cp -r results M1Philippe_results/resultsp72D1e12
rm results/*
mpirun -np 8 nth -p 7 -m data/square01_quad.mesh -rs 4 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e11 | tee M1Philippe_p72D1e11.out
cp -r results M1Philippe_results/resultsp72D1e11
rm results/*
mpirun -np 8 nth -p 7 -m data/square01_quad.mesh -rs 4 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e10 | tee M1Philippe_p72D1e10.out
cp -r results M1Philippe_results/resultsp72D1e10
rm results/*
mpirun -np 8 nth -p 7 -m data/square01_quad.mesh -rs 4 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e9 | tee M1Philippe_p72D1e9.out
cp -r results M1Philippe_results/resultsp72D1e9
rm results/*
mpirun -np 8 nth -p 7 -m data/square01_quad.mesh -rs 4 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -a0 1e8 | tee M1Philippe_p72D1e8.out
cp -r results M1Philippe_results/resultsp72D1e8
# problem 6 - step in density and step in temperature in 3D as at shock
rm results/*
mpirun -np 8 nth -p 7 -m data/cube01_hex.mesh -rs 4 -tf 0.0 -ok 1 -ot 0 -vis -fa -print -a0 1e12 | tee M1Philippe_p73D1e12.out
cp -r results M1Philippe_results/resultsp73D1e12
rm results/*
mpirun -np 8 nth -p 7 -m data/cube01_hex.mesh -rs 4 -tf 0.0 -ok 1 -ot 0 -vis -fa -print -a0 1e11 | tee M1Philippe_p73D1e11.out
cp -r results M1Philippe_results/resultsp73D1e11
rm results/*
mpirun -np 8 nth -p 7 -m data/cube01_hex.mesh -rs 4 -tf 0.0 -ok 1 -ot 0 -vis -fa -print -a0 1e10 | tee M1Philippe_p73D1e10.out
cp -r results M1Philippe_results/resultsp73D1e10
rm results/*
mpirun -np 8 nth -p 7 -m data/cube01_hex.mesh -rs 4 -tf 0.0 -ok 1 -ot 0 -vis -fa -print -a0 1e9 | tee M1Philippe_p73D1e9.out
cp -r results M1Philippe_results/resultsp73D1e9
rm results/*
mpirun -np 8 nth -p 7 -m data/cube01_hex.mesh -rs 4 -tf 0.0 -ok 1 -ot 0 -vis -fa -print -a0 1e8 | tee M1Philippe_p73D1e8.out
cp -r results M1Philippe_results/resultsp73D1e8
