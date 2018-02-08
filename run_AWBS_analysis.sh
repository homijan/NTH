#! /bin/bash
# Philippes tests, on the edge of locality.
# Explicit Efield.
#mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -S0 1.0 -E0 1.0 -a0 5e31
# Mimiced Efield.
#mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 5e31

# Pascal's setting for nonlocal test
# First, diffusive, case sets sigma 1000x higher, which assures SH solution.
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 3e33 -Tgrad 18 -Z 2 -ne 5e19
# Nonlocal solution very well corresponding to Pascal's solution with Aladin.
mpirun -np 8 nth -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 3e30 -Tgrad 18 -Z 2 -ne 5e19
 
mpirun -np 8 nth -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 2 -ot 1 -no-vis -fa -print -Tmax 800.5 -Tmin 799.5 -Tgrad 2.3 -S0 1.0 -E0 1.0 -a0 1e31 -Z 4 -ne 5e20
cp fe.txt docs/fe_data/fe_Ecorrect.txt

mpirun -np 8 nth -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 2 -ot 1 -no-vis -fa -print -Tmax 800.5 -Tmin 799.5 -Tgrad 2.3 -a0 1e31 -Z 4 -ne 5e20
cp fe.txt docs/fe_data/fe_Emimic.txt

cd docs
python AWBSf1_integrate.py -Z 4 -T 8e2 -g -2.3 -s 1e31 -n 5e20
