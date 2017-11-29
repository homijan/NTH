#!/bin/sh
cd triplepoint
~/MFEM_2017/glvis-3.3/glvis -run ../create_parallel_temperature_triplepoint_figures.glvs
~/MFEM_2017/glvis-3.3/glvis -run ../create_parallel_density_triplepoint_figures.glvs
~/MFEM_2017/glvis-3.3/glvis -run ../create_parallel_nonlocalI0_triplepoint_figures.glvs
~/MFEM_2017/glvis-3.3/glvis -run ../create_parallel_nonlocalI1z_triplepoint_figures.glvs
~/MFEM_2017/glvis-3.3/glvis -run ../create_parallel_nonlocalI1x_triplepoint_figures.glvs
cd ../OmegaLaser
~/MFEM_2017/glvis-3.3/glvis -run ../create_parallel_material_OmegaLaser_figures.glvs
~/MFEM_2017/glvis-3.3/glvis -run ../create_parallel_density_OmegaLaser_figures.glvs
~/MFEM_2017/glvis-3.3/glvis -run ../create_parallel_temperature_OmegaLaser_figures.glvs
~/MFEM_2017/glvis-3.3/glvis -run ../create_parallel_nonlocalI0_OmegaLaser_figures.glvs
~/MFEM_2017/glvis-3.3/glvis -run ../create_parallel_nonlocalI1z_OmegaLaser_figures.glvs
~/MFEM_2017/glvis-3.3/glvis -run ../create_parallel_nonlocalI1x_OmegaLaser_figures.glvs
