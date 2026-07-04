prefix=/glade/derecho/scratch/taoyuan/Proj_RSL_3Dvisc/case_B_ICE1_VE300_VertAve1/case_B_ICE1_VE300_VertAve1
nstep=1
prefix_out=./case_B_ICE1_VE300_VertAve1/out
nproc_surf=192
nprocz=2
noxy=33
noz=33
./process_visc_all_layers $prefix 1 ${prefix_out} ${nproc_surf} ${nprocz} ${noxy} ${noz}
