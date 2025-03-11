# interpolate the uplift and geoid to regular grid (same grid as the input ice model)
    # requires a prefix, interpolate all the steps
    # save file to same directory as the input file
# usage: bash interp_ug_regular_grid.sh prefix nepoch resolution
# save to prefix.map_uplift.epoch${TIMESTEP}.regular.grd and prefix.map_geoid.epoch${TIMESTEP}.regular.grd


# module load ncarenv-basic/23.09  gcc/12.2.0 gmt/6.5.0

module load ncarenv/23.10  gcc/12.2.0  gmt/6.5.0

# check if gmt can be found, if not, load module above
if ! command -v gmt &> /dev/null
then
    # try to load module
    echo "gmt could not be found, try module spider gmt"
    module spider gmt
    exit
fi

data_prefix=$1
nepoch=$2   # number of epochs
resolution=$3  # resolution of the regular grid (in degree)

half_res=$(awk -v res="$resolution" 'BEGIN {printf "%.2f", res/2}')
echo "half_res: ${half_res}"
range="${half_res}/$(echo "360-${half_res}" | bc)/$(echo "-90+${half_res}" | bc)/$(echo "90-${half_res}" | bc)"
echo "range: ${range}"


export GMT_END_SHOW=off

# each call do for all models, just one step
one_step_saving() {
    # $1: TIMESTEP -> epoch

    TIMESTEP=$1
    model_prefix=$2

    echo "epoch: ${TIMESTEP}"

    # cumulative uplift
    SVE_data_path=${model_prefix}.map_uplift.epoch${TIMESTEP}
    gmt nearneighbor ${SVE_data_path} -R${range} -I${resolution} -E-999999 -G${SVE_data_path}.regular.grd -S200k -N4/2 -Lg
    gmt grd2xyz ${SVE_data_path}.regular.grd > ${SVE_data_path}.regular.xyz --FORMAT_FLOAT_OUT=%11.5E

    # geoid
    SVE_data_path=${model_prefix}.map_geoid.epoch${TIMESTEP}
    gmt nearneighbor ${SVE_data_path} -R${range} -I${resolution} -E-999999 -G${SVE_data_path}.regular.grd -S200k -N4/2 -Lg
    gmt grd2xyz ${SVE_data_path}.regular.grd > ${SVE_data_path}.regular.xyz --FORMAT_FLOAT_OUT=%11.5E
}


# find all the steps
# steps=$(ls ${data_prefix}.map_geoid.epoch* | grep -oP '\.\K\d+' | sort -n -u )
# echo steps: ${steps}
# step_array=($steps)
# for step in "${step_array[@]}"; do ...


# run all the steps
MAX_CPU=12
cpu_id=0

for step in $(seq 0 $((${nepoch}-1)) ); do
    ((cpu_id=cpu_id%MAX_CPU)); ((cpu_id++==0)) && wait
    one_step_saving $step $data_prefix &
done
