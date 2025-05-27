#!/bin/bash
# run with bash
# modify the parameters in this file to control the steps to run
# may need to modify the code to change some default settings (e.g. gmt path, gmt interpolation parameters)




# parameters to set:
load_stages_time_file=/glade/derecho/scratch/taoyuan/Proj_Nonlinear_GIA/ice6g_40ka_stages_timeB.dat  
                        # the stage_time file used in SVE calculation
                        # filename, to save file for epoch to timestep conversion
SVE_output_prefix=/glade/derecho/scratch/taoyuan/Proj_Nonlinear_GIA/case_a0/case_a0

postp_id=case_a0

ncpu_surface=192        #48    # number of cpus in the surface
ncpu_z=2            # number of cpus in the z direction
node_x_y=33         #41      # number of nodes in x and y directions (for each cpu)
nepochs=81
resolution=1    # resolution of regular grid (in degree)


###################### end of parameters ##########################
# Step A: get timestep for each **epoch** (same epoch as input ice model), and get RSL_c
    # the suffix of saved files is based on epoch (0,1,.. nepoch-1), not timestep. Thats why we need to convert epoch to timestep.

# hardcoded parameters
position_of_RSL_c_in_Header=9   # the position of RSL_c in the header file 
saved_stage_timestep_file=${postp_id}.epoch_to_timestep.txt 
output_RSL_c_file=${postp_id}.RSL_c    # the output RSL_c file, only for epoch, not every timestep output
save_to_dir=./combined_files_${postp_id}/
save_to_file_prefix=${save_to_dir}/${postp_id}




a=$load_stages_time_file
b=$saved_stage_timestep_file
c=${SVE_output_prefix}.time_dep
d=$position_of_RSL_c_in_Header
e=$output_RSL_c_file


run_step_A=1
if [ $run_step_A -eq 1 ]; then
    ./A1_epoch_to_timestep  $a  $b
    ./A2_get_RSL_c  $b  $c  $d  $e
fi


# Step B: get uplift and geoid on regular grid (same grid as input ice model)


run_step_B=1
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
run_step_C=1
if [ $run_step_C -eq 1 ]; then
    ./A3_refine_ocnfunc input.refine_ocnfunc
fi

# finally, remove intermediate files

# rm -r ./combined_files
