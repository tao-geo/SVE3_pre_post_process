# gen ascii file from nc (grd) files
# for IceT.I6F_C.131QB_VM5a_1deg.nc, transform lat direction: in nc file it is -89.5 to 89.5
#       Now we need from 89.5 to -89.5
# No need.

nepoch=122
outdir=./ICE6G_original/
mkdir -p $outdir

for (( epoch=0; epoch<$nepoch; epoch++ ))
do
    echo "epoch: $epoch"
    infile=IceT.I6F_C.131QB_VM5a_1deg.nc?stgit\[${epoch}\]  #

    # transform lat direction: in nc file it is -89.5 to 89.5
    # temp1=${outdir}/temp1.nc
    # gmt grdedit $infile -Ev -G$temp1

    temp1=$infile
    # convert to ascii
    gmt grd2xyz $temp1 > ${outdir}/ice180x360.${epoch} --FORMAT_FLOAT_OUT=%11.5E

    # rm $temp1
    # exit
done

