#!/bin/bash
# run with bash
# modify the parameters in this file to control the steps to run
# may need to modify the code to change some default settings (e.g. gmt path, gmt interpolation parameters)

# set to 0 to skip the step
run_step_A=0
run_step_B=0
run_step_C=0


# Step A: get timestep for each **epoch** (same epoch as input ice model), and get RSL_c
    # the suffix of saved files is based on epoch (0,1,.. nepoch-1), not timestep. Thats why we need to convert epoch to timestep.

# parameters to set:
stage_time_file=../ice6g_122ka_stages_timeB.dat  
                        # the stage_time file used in SVE calculation

saved_stage_timestep_file=ice6g_122ka_stages_timeB.epoch_to_timestep.txt   
                        # filename, to save file for epoch to timestep conversion

one_topo_file_without_suffix=../../link_scratch_BM_SLE/case_test/case_test.topo_s.1  
                        # one topo file without suffix (suffix is timestep)

position_of_RSL_c_in_Header=8   # the position of RSL_c in the header file 
output_RSL_c_file=case_test.RSL_c    # the output RSL_c file, only for epoch, not every timestep output
# end of parameters

a=$stage_time_file
b=$saved_stage_timestep_file
c=$one_topo_file_without_suffix
d=$position_of_RSL_c_in_Header
e=$output_RSL_c_file

if [ $run_step_A -eq 1 ]; then
    ./A1_epoch_to_timestep  $a  $b
    ./A2_get_RSL_c  $b  $c  $d  $e
fi


# Step B: get uplift and geoid on regular grid (same grid as input ice model)

# parameters to set:
SVE_output_prefix=../../link_scratch_BM_SLE/case_test/case_test
save_to_dir=./combined_files/
save_to_file_prefix=${save_to_dir}/case_test
ncpu_surface=48        #48    # number of cpus in the surface
ncpu_z=2            # number of cpus in the z direction
node_x_y=41         #41      # number of nodes in x and y directions (for each cpu)
nepochs=122
resolution=1    # resolution of regular grid (in degree)
# end of parameters

if [ $run_step_B -eq 1 ]; then
    echo "***** Step B *****"
    mkdir -p $save_to_dir
    ./combine_u_g_files ${SVE_output_prefix} \
                        ${saved_stage_timestep_file} \
                        ${save_to_file_prefix} \
                        ${ncpu_surface}  ${ncpu_z}  ${node_x_y} 

    sh ./interp_ug_regular_grid.sh      ${save_to_file_prefix}      $nepochs    $resolution

fi

wait
# Step C: update ocean function and ice model (check for floating ice). Following [Kendall 2005]

# parameters: 
# !!! modify refine_ocnfunc.input !!!

if [ $run_step_C -eq 1 ]; then
    ./A3_refine_ocnfunc refine_ocnfunc.input
fi

# finally, remove intermediate files

rm -r ./combined_files