# interpolate the disp and geoid to regular grid (same grid as the input ice model)
    # requires a prefix, interpolate all the steps
    # save file to same directory as the input file
# usage: bash interp_ug_regular_grid.sh prefix fn_header resolution
# Modify gmt nearneighbor parameters if needed
    # -S200k: search radius
# Modify MAX_CPU if needed: the maximum number of threads

#########################
####### Check GMT #######
#########################

# for derecho
# module load ncarenv-basic/23.09  gcc/12.2.0 gmt/6.5.0

# for casper
module load ncarenv/23.10  gcc/12.2.0  gmt/6.5.0

# check if gmt can be found, if not, load module above
if ! command -v gmt &> /dev/null
then
    # try to load module
    echo "gmt could not be found, try module spider gmt"
    module spider gmt
    exit
fi

#########################
####### Parameters #######
#########################

data_prefix=$1
fn_header=$2   # header
resolution=$3  # resolution of the regular grid (in degree)

MAX_CPU=$4

#########################




# half_res=$(awk -v res="$resolution" 'BEGIN {printf "%.2f", res/2}')
# echo "half_res: ${half_res}"
# range="${half_res}/$(echo "360-${half_res}" | bc)/$(echo "-90+${half_res}" | bc)/$(echo "90-${half_res}" | bc)"
range="0/360/-90/90"
echo "range: ${range}"


export GMT_END_SHOW=off


one_step_saving() {
    # $1: TIMESTEP -> epoch

    # TIMESTEP=$1
    SVE_data_path=$1

    gmt nearneighbor ${SVE_data_path} -R${range} -I${resolution} -E-999999 -G${SVE_data_path}.regular.grd -S200k -N4/2 #-Lg
    gmt grd2xyz ${SVE_data_path}.regular.grd > ${SVE_data_path}.regular.xyz --FORMAT_FLOAT_OUT=%11.5E

    # rm ${SVE_data_path}.regular.grd
    # rm ${SVE_data_path}
}


# find all the steps
steps=$(awk '(NR>1){printf "%d ", $1}' ${fn_header})
echo steps: ${steps}
step_array=($steps)


######### Interpolate for Cumulative Disp/Geoid #########
echo ""
echo "Interpolating for Cumulative Disp/Geoid: "

cpu_id=0
for step in "${step_array[@]}"; do
    echo -n " ${step} "

    SVE_data_path=${data_prefix}.map_disp_and_geoid.${step}

    ((cpu_id=cpu_id%MAX_CPU)); ((cpu_id++==0)) && wait
    SVE_uplift_path=${data_prefix}.map_uplift.${step}
    awk '{print $1, $2, $3}' ${SVE_data_path} > ${SVE_uplift_path}
    one_step_saving   $SVE_uplift_path &

    ((cpu_id=cpu_id%MAX_CPU)); ((cpu_id++==0)) && wait
    SVE_south_path=${data_prefix}.map_south.${step}
    awk '{print $1, $2, $4}' ${SVE_data_path} > ${SVE_south_path}
    one_step_saving   $SVE_south_path &

    ((cpu_id=cpu_id%MAX_CPU)); ((cpu_id++==0)) && wait
    SVE_east_path=${data_prefix}.map_east.${step}
    awk '{print $1, $2, $5}' ${SVE_data_path} > ${SVE_east_path}
    one_step_saving   $SVE_east_path &

    ((cpu_id=cpu_id%MAX_CPU)); ((cpu_id++==0)) && wait
    SVE_geoid_path=${data_prefix}.map_geoid.${step}
    awk '{print $1, $2, $6}' ${SVE_data_path} > ${SVE_geoid_path}
    one_step_saving   $SVE_geoid_path &

done

wait

######### Interpolate for Disp/Geoid Rate #########
echo "Interpolating for Disp/Geoid Rate: "
cpu_id=0
for step in "${step_array[@]}"; do
    echo -n " ${step} "

    SVE_data_path=${data_prefix}.map_rate_disp_and_geoid.${step}

    ((cpu_id=cpu_id%MAX_CPU)); ((cpu_id++==0)) && wait
    SVE_uplift_path=${data_prefix}.map_rate_uplift.${step}
    awk '{print $1, $2, $3}' ${SVE_data_path} > ${SVE_uplift_path}
    one_step_saving   $SVE_uplift_path &

    ((cpu_id=cpu_id%MAX_CPU)); ((cpu_id++==0)) && wait
    SVE_south_path=${data_prefix}.map_rate_south.${step}
    awk '{print $1, $2, $4}' ${SVE_data_path} > ${SVE_south_path}
    one_step_saving   $SVE_south_path &

    ((cpu_id=cpu_id%MAX_CPU)); ((cpu_id++==0)) && wait
    SVE_east_path=${data_prefix}.map_rate_east.${step}
    awk '{print $1, $2, $5}' ${SVE_data_path} > ${SVE_east_path}
    one_step_saving   $SVE_east_path &

    ((cpu_id=cpu_id%MAX_CPU)); ((cpu_id++==0)) && wait
    SVE_geoid_path=${data_prefix}.map_rate_geoid.${step}
    awk '{print $1, $2, $6}' ${SVE_data_path} > ${SVE_geoid_path}
    one_step_saving   $SVE_geoid_path &

done

wait


##### delete some files #####
cpu_id=0
for step in "${step_array[@]}"; do
    ((cpu_id=cpu_id%MAX_CPU)); ((cpu_id++==0)) && wait

    (rm ${data_prefix}.map_uplift.${step}
    rm ${data_prefix}.map_south.${step}
    rm ${data_prefix}.map_east.${step}
    rm ${data_prefix}.map_geoid.${step}
    rm ${data_prefix}.map_rate_uplift.${step}
    rm ${data_prefix}.map_rate_south.${step}
    rm ${data_prefix}.map_rate_east.${step}
    rm ${data_prefix}.map_rate_geoid.${step}) &
done

echo ""