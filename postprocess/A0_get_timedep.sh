# get time dependent quantities from SVE output files (*.time_dep)
# a=$start_timestep          # starting_timestep
# b=$end_timestep            # ending_timestep
# c=${SVE_output_prefix}     # prefix for topo files, need this input to get at least one surface cpu
# d=${fn_header}             # output file name

echo "starting timestep $1 " 
echo "ending timestep $2 "
echo "from file ${3}.time_dep" 
echo "writting file $4 "

input_file=${3}.time_dep
output_file=$4

# get header
head -n 1 ${input_file} > ${output_file}

# get time dependent quantities
awk -v start_timestep=$1 -v end_timestep=$2 '{
    if (NR>1 && $1>=start_timestep && $1<=end_timestep) {
        print $0
    }
}' ${input_file} >> ${output_file}