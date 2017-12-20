#!/bin/bash
# problem 5 in 1D
# generate figs locally
cd resultsp51DKn1e-9
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn1e-9Efield
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.0001
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.001
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.01
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.02
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.03
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.04
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.05
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.06
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.07
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.08
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.09
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.1
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.2
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.3
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.4
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.5
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.6
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.7
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.8
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn0.9
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ../resultsp51DKn1.0
../../../glvis-3.3/glvis -run ../create_parallel_M1Philippe_1Dfigures.glvs
cd ..
# copy figs
cp resultsp51DKn0.0001/density.png ../docs/figs/density_p51DNonlocAnal.png
cp resultsp51DKn0.0001/temperature.png ../docs/figs/temperature_p51DNonlocAnal.png
cp resultsp51DKn1e-9/Kn.png ../docs/figs/Kn_p51DKn1e-9.png
cp resultsp51DKn1e-9/hflux.png ../docs/figs/hflux_p51DKn1e-9.png
cp resultsp51DKn1e-9/jcurrent.png ../docs/figs/jcurrent_p51DKn1e-9.png
cp resultsp51DKn1e-9Efield/Kn.png ../docs/figs/Kn_p51DKn1e-9Efield.png
cp resultsp51DKn1e-9Efield/hflux.png ../docs/figs/hflux_p51DKn1e-9Efield.png
cp resultsp51DKn1e-9Efield/jcurrent.png ../docs/figs/jcurrent_p51DKn1e-9Efield.png
cp resultsp51DKn0.0001/Kn.png ../docs/figs/Kn_p51DKn0.0001.png
cp resultsp51DKn0.0001/hflux.png ../docs/figs/hflux_p51DKn0.0001.png
cp resultsp51DKn0.0001/jcurrent.png ../docs/figs/jcurrent_p51DKn0.0001.png
cp resultsp51DKn0.001/Kn.png ../docs/figs/Kn_p51DKn0.001.png
cp resultsp51DKn0.001/hflux.png ../docs/figs/hflux_p51DKn0.001.png
cp resultsp51DKn0.001/jcurrent.png ../docs/figs/jcurrent_p51DKn0.001.png
cp resultsp51DKn0.01/Kn.png ../docs/figs/Kn_p51DKn0.01.png
cp resultsp51DKn0.01/hflux.png ../docs/figs/hflux_p51DKn0.01.png
cp resultsp51DKn0.01/jcurrent.png ../docs/figs/jcurrent_p51DKn0.01.png
cp resultsp51DKn0.02/Kn.png ../docs/figs/Kn_p51DKn0.02.png
cp resultsp51DKn0.02/hflux.png ../docs/figs/hflux_p51DKn0.02.png
cp resultsp51DKn0.02/jcurrent.png ../docs/figs/jcurrent_p51DKn0.02.png
cp resultsp51DKn0.03/Kn.png ../docs/figs/Kn_p51DKn0.03.png
cp resultsp51DKn0.03/hflux.png ../docs/figs/hflux_p51DKn0.03.png
cp resultsp51DKn0.03/jcurrent.png ../docs/figs/jcurrent_p51DKn0.03.png
cp resultsp51DKn0.04/Kn.png ../docs/figs/Kn_p51DKn0.04.png
cp resultsp51DKn0.04/hflux.png ../docs/figs/hflux_p51DKn0.04.png
cp resultsp51DKn0.04/jcurrent.png ../docs/figs/jcurrent_p51DKn0.04.png
cp resultsp51DKn0.05/Kn.png ../docs/figs/Kn_p51DKn0.05.png
cp resultsp51DKn0.05/hflux.png ../docs/figs/hflux_p51DKn0.05.png
cp resultsp51DKn0.05/jcurrent.png ../docs/figs/jcurrent_p51DKn0.05.png
cp resultsp51DKn0.06/Kn.png ../docs/figs/Kn_p51DKn0.06.png
cp resultsp51DKn0.06/hflux.png ../docs/figs/hflux_p51DKn0.06.png
cp resultsp51DKn0.06/jcurrent.png ../docs/figs/jcurrent_p51DKn0.06.png
cp resultsp51DKn0.07/Kn.png ../docs/figs/Kn_p51DKn0.07.png
cp resultsp51DKn0.07/hflux.png ../docs/figs/hflux_p51DKn0.07.png
cp resultsp51DKn0.07/jcurrent.png ../docs/figs/jcurrent_p51DKn0.07.png
cp resultsp51DKn0.08/Kn.png ../docs/figs/Kn_p51DKn0.08.png
cp resultsp51DKn0.08/hflux.png ../docs/figs/hflux_p51DKn0.08.png
cp resultsp51DKn0.08/jcurrent.png ../docs/figs/jcurrent_p51DKn0.08.png
cp resultsp51DKn0.09/Kn.png ../docs/figs/Kn_p51DKn0.09.png
cp resultsp51DKn0.09/hflux.png ../docs/figs/hflux_p51DKn0.09.png
cp resultsp51DKn0.09/jcurrent.png ../docs/figs/jcurrent_p51DKn0.09.png
cp resultsp51DKn0.1/Kn.png ../docs/figs/Kn_p51DKn0.1.png
cp resultsp51DKn0.1/hflux.png ../docs/figs/hflux_p51DKn0.1.png
cp resultsp51DKn0.1/jcurrent.png ../docs/figs/jcurrent_p51DKn0.1.png
cp resultsp51DKn0.2/Kn.png ../docs/figs/Kn_p51DKn0.2.png
cp resultsp51DKn0.2/hflux.png ../docs/figs/hflux_p51DKn0.2.png
cp resultsp51DKn0.2/jcurrent.png ../docs/figs/jcurrent_p51DKn0.2.png
cp resultsp51DKn0.3/Kn.png ../docs/figs/Kn_p51DKn0.3.png
cp resultsp51DKn0.3/hflux.png ../docs/figs/hflux_p51DKn0.3.png
cp resultsp51DKn0.3/jcurrent.png ../docs/figs/jcurrent_p51DKn0.3.png
cp resultsp51DKn0.4/Kn.png ../docs/figs/Kn_p51DKn0.4.png
cp resultsp51DKn0.4/hflux.png ../docs/figs/hflux_p51DKn0.4.png
cp resultsp51DKn0.4/jcurrent.png ../docs/figs/jcurrent_p51DKn0.4.png
cp resultsp51DKn0.5/Kn.png ../docs/figs/Kn_p51DKn0.5.png
cp resultsp51DKn0.5/hflux.png ../docs/figs/hflux_p51DKn0.5.png
cp resultsp51DKn0.5/jcurrent.png ../docs/figs/jcurrent_p51DKn0.5.png
cp resultsp51DKn0.6/Kn.png ../docs/figs/Kn_p51DKn0.6.png
cp resultsp51DKn0.6/hflux.png ../docs/figs/hflux_p51DKn0.6.png
cp resultsp51DKn0.6/jcurrent.png ../docs/figs/jcurrent_p51DKn0.6.png
cp resultsp51DKn0.7/Kn.png ../docs/figs/Kn_p51DKn0.7.png
cp resultsp51DKn0.7/hflux.png ../docs/figs/hflux_p51DKn0.7.png
cp resultsp51DKn0.7/jcurrent.png ../docs/figs/jcurrent_p51DKn0.7.png
cp resultsp51DKn0.8/Kn.png ../docs/figs/Kn_p51DKn0.8.png
cp resultsp51DKn0.8/hflux.png ../docs/figs/hflux_p51DKn0.8.png
cp resultsp51DKn0.8/jcurrent.png ../docs/figs/jcurrent_p51DKn0.8.png
cp resultsp51DKn0.9/Kn.png ../docs/figs/Kn_p51DKn0.9.png
cp resultsp51DKn0.9/hflux.png ../docs/figs/hflux_p51DKn0.9.png
cp resultsp51DKn0.9/jcurrent.png ../docs/figs/jcurrent_p51DKn0.9.png
cp resultsp51DKn1.0/Kn.png ../docs/figs/Kn_p51DKn1.0.png
cp resultsp51DKn1.0/hflux.png ../docs/figs/hflux_p51DKn1.0.png
cp resultsp51DKn1.0/jcurrent.png ../docs/figs/jcurrent_p51DKn1.0.png

# Now create a .gif from the hflux.png snapshots.
rm *.png
cp resultsp51DKn0.001/hflux.png snapshot-00.png
cp resultsp51DKn0.01/hflux.png snapshot-01.png
cp resultsp51DKn0.02/hflux.png snapshot-02.png
cp resultsp51DKn0.03/hflux.png snapshot-03.png
cp resultsp51DKn0.04/hflux.png snapshot-04.png
cp resultsp51DKn0.05/hflux.png snapshot-05.png
cp resultsp51DKn0.06/hflux.png snapshot-06.png
cp resultsp51DKn0.07/hflux.png snapshot-07.png
cp resultsp51DKn0.08/hflux.png snapshot-08.png
cp resultsp51DKn0.09/hflux.png snapshot-09.png
cp resultsp51DKn0.1/hflux.png snapshot-10.png
cp resultsp51DKn0.2/hflux.png snapshot-11.png
cp resultsp51DKn0.3/hflux.png snapshot-12.png
cp resultsp51DKn0.4/hflux.png snapshot-13.png
cp resultsp51DKn0.5/hflux.png snapshot-14.png
cp resultsp51DKn0.6/hflux.png snapshot-15.png
cp resultsp51DKn0.7/hflux.png snapshot-16.png
cp resultsp51DKn0.8/hflux.png snapshot-17.png
cp resultsp51DKn0.9/hflux.png snapshot-18.png
cp resultsp51DKn1.0/hflux.png snapshot-19.png

convert -delay 100 -loop 0 *.png Nonlocality_analysis_woutE.gif
rm *.png
cp Nonlocality_analysis_woutE.gif ../docs/ 
