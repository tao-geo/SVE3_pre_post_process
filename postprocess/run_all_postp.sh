# run with bash run_all_postp.sh

## control which steps to run: 0 means skip ##
run_step_1=1

run_step_2=0
run_step_2_MultiThreads=1

run_step_3=1

run_step_4=1



#######################    parameters needed to be specified    #######################

SVE_output_prefix=../BM_ICE6G/case6   # prefix of the SVE generated output files
fn_header=case6.header     # the header file name that will be generated in the first step, and used in the following steps
                                # it constains the timesteps, time, delta_t, RSL_c and other information.

start_timestep=392    # the starting timestep to do post process. (not necessarily the first timestep of SVE generated output files)
end_timestep=592    # the ending timestep to do post process.

save_to_dir=./combined_files_case_v3/        # the directory to save the combined files (and other generated files)
save_to_file_prefix=${save_to_dir}/case_v3   # the prefix of the combined files
ncpu_surface=192 #48    # number of cpus in the surface
ncpu_z=2            # number of cpus in the z direction
node_x_y=33   #41      # number of nodes in x and y directions (for each cpu)

nthreads=12         # used in the multi-thread version

resolution=0.5    # resolution of regular grid (in degree)
                # the regular grid's resolution created by GMT nearneighbor method,
                # on a regular grid with lon from 0 to 360, lat from -90 to 90. (Caution: not from resolution/2 to 360-resolution/2, -90+resolution/2 to 90-resolution/2)
nlat=361   # number of latitudes in the regular grid
nlon=721   # number of longitudes in the regular grid  (although this could be calculated from $resolution)


fn_RSL_c=case_v3.RSL_c    # the filename for generated RSL_c (timestep, time, RSL_c)

fn_RSL_sites=RSL_sites_groupA.posi  # (USER PREPARED FILE) the file that contains the lon, lat (and may also name of sites)
# !!! fn_RSL_sites's first line should be the number of sites

fn_RSL_out=RSL_sites_groupA.out.case_v3  # the output file that contains the calculated RSL for fn_RSL_sites




################  no need to modify the following part, unless need to change some default setting ################

# default settings that may need to be changed
# 1. gmt path: need to modify A2_interp_to_regular_grid.sh if gmt could not be found
# 2. gmt interpolation parameters: modify A2_interp_to_regular_grid.sh
# 3. other settings: modify the .c files if needed




################   Step 1: generate header file (timestep, time, dt, etc)  ################

# First, generate header file from output files, 
    # including the timesteps and time that will be used in the following process.


# parameters

a=$start_timestep          # starting_timestep
b=$end_timestep            # ending_timestep
c=${SVE_output_prefix}     # prefix for topo files, need this input to get at least one surface cpu
d=${fn_header}             # output file name

# run
if [ $run_step_1 -eq 1 ]; then
    sh ./A0_get_timedep.sh $a $b $c $d
fi





################  Step 2, Combine SVE generated files  ################
# combine output files from multiple cpus into one file for each timestep.
# it is possible to remove some timesteps (lines) in fn_header if they are not needed in the following process.
# 
# There is a multi-thread version of this step, which is useful when the number of timesteps is large.
# The multi-thread version splits the header file (each line is one timestep) into multiple files, and each file is processed by one thread.


# parameters

a=${SVE_output_prefix}   # prefix of the SVE generated output files
b=${fn_header}            # header file
c=${save_to_file_prefix}  # save files using this prefix
d=$ncpu_surface           # number of cpus in the surface
e=$ncpu_z                # number of cpus in the z direction
f=$node_x_y              # number of nodes in x and y directions

# run
mkdir -p $save_to_dir

# single thread version
if [ $run_step_2 -eq 1 ]; then
    ./A1_combine_u_g_files $a $b $c $d $e $f 
fi

# multi-thread version
# use split to split the header file into multiple files, and then run the multi-thread version
# use nthreads defined above

if [ $run_step_2_MultiThreads -eq 1 ]; then
    
    # first make sure no temp files exist
    rm -f ${fn_header}_split_head ${fn_header}_split_body ${fn_header}_split_body_*

    # split the header file (body) into multiple files (nfiles ~= nthreads)
    # run A1_combine_u_g_files for each file on a separate thread (background process)

    head -n 1 $fn_header > ${fn_header}_split_head
    tail -n +2 $fn_header > ${fn_header}_split_body
    nrows=$(wc -l ${fn_header}_split_body | awk '{print $1}')
    nrows_per_thread=$((nrows/nthreads))
    split -l $nrows_per_thread ${fn_header}_split_body ${fn_header}_split_body_
    all_files=$(ls ${fn_header}_split_body_*)

    for file in $all_files; do
        cat ${fn_header}_split_head $file > ${file}_temp
        ./A1_combine_u_g_files $a \
                        ${file}_temp \
                        $c $d $e $f &
    done
    wait # wait for all background processes to

    rm ${fn_header}_split_head ${fn_header}_split_body ${fn_header}_split_body_*
fi


################  Step 3, interpolate disp and geoid to regular grid  ################

# fn_header=case_v1_t2.header_split_body_ak_temp # debug

# parameters

a=${save_to_file_prefix}  # prefix of the combined files
b=${fn_header}            # header file
c=$resolution             # resolution
d=$nthreads               # number of threads

if [ $run_step_3 -eq 1 ]; then
    echo ""
    echo ""
    echo "**** interpolate disp and geoid to regular grid (including cumulative and rate) ****"
    echo file_prefix: ${save_to_file_prefix}
    echo fn_header: $fn_header
    echo resolution: $resolution
    echo nthreads: $nthreads

    ./A2_interp_to_regular_grid.sh $a $b $c $d
fi



################  Step 4: RSL calculation  ################

# 4.1 get RSL_c


a=$fn_header
b=$fn_RSL_c

# run
if [ $run_step_4 -eq 1 ]; then
    echo getting RSL_c from fn_header: $a
    echo saving to fn_RSL_c: $b

    awk '(NR>1){printf "%d %f %f\n", $1, $2, $9}' $a > $b
    echo RSL_c saved in $b
fi


# 4.2 get RSL for sites
# parameters


a=$fn_RSL_sites
b=$fn_RSL_out
c=$fn_RSL_c
d=${save_to_file_prefix}
e=$nlat
f=$nlon

# ./A3_get_RSL_sites
# run
if [ $run_step_4 -eq 1 ]; then
    echo "**** A3_get_RSL_sites ****"

    ./A3_get_RSL_sites $a $b $c $d $e $f
fi




# # 5 plot RSL for sites
# nsites=$(head -n 1 $fn_RSL_sites)
# echo nsites: $nsites



# if [ $run_step_5 -eq 1 ]; then
#     echo "**** plot RSL for sites ****"
#     for i in $(seq 1 $nsites); do
#         echo "site $i"
#         ./plot_RSL_sites $fn_RSL_out $i
#     done
# fi
