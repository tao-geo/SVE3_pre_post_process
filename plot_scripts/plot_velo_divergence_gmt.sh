# Purpose: plot divergence and horizontal velocity vector on the same map

# need to be combined with other scripts to make the final figure:
#    mk_figA_s2_combine_panels.sh

########################### Required files  #############################

# NCAR_DATA_v4/download_from_ncar/dilitation_rate_regulargrid/${case_id}/${case_id}.map_dilitation.${step}.regular.CmCorrect.grd
# NCAR_DATA_v4/download_from_ncar/velocity_regulargrid_allsteps/${case_id}_rm_CM/${case_id}.map_incr_south.${step}.regular.grd
# NCAR_DATA_v4/download_from_ncar/velocity_regulargrid_allsteps/${case_id}_rm_CM/${case_id}.map_incr_east.${step}.regular.grd

# case_id: case6, case10_3Dcrust (case10)
# step: 648, 584, 536, 480, 104

############################# Functions #############################
FN_RIDGE=/home/tao/Prj_GIA_PlateTectonic/Data/viscosity_model/PlateBoundaryInter_0Seton

create_north_from_south() {
    # $1: input south file, $2: output north file
    gmt grdmath $1 -1 MUL = $2
}

main () {
    # $1: case id , $2, step  (descrip purpose)
    local case_id="$1"
    local step=$2
    # mkdir -p ./fig_${case_id}
    local FILE_div=$3
    local FILE_east=$4
    local File_north=$5
    local output_file=$6    # fig (without extension)
    local title="$7"
    local scale=$8
    local cpt_range=$9

    gmt set MAP_TITLE_OFFSET -0.25c
    
    ## transform to mm/yr
    # gmt grdmath $FILE_div 1000 MUL = ${output_file}_div_.grd
    gmt grdmath $FILE_east 1000 MUL = ${output_file}_east_mm.grd
    gmt grdmath $File_north 1000 MUL = ${output_file}_north_mm.grd

    ## plot divergence rate
    cpt_step=$(awk -v cpt_range=${cpt_range} 'BEGIN{print cpt_range/5}')
    echo $cpt_step
    # echo $cpt_step
    gmt makecpt -Croma -T-${cpt_range}/${cpt_range}  -Do -G-0.7/0.7 -I  #roma -G-0.8/0.8 -I
    # gmt makecpt -C123/50/148,194/165/207,166/219/160,0/136/55 -T-${cpt_range}/${cpt_range}  -Do
    # 123/50/148
    # 194/165/207
    # 166/219/160
    # 0/136/55
    # gmt makecpt -Cvik -T-${cpt_range}/${cpt_range}  -Do -G-0.5/0.5
    # exit
    # gmt makecpt -Chawaii -T-${cpt_range}/${cpt_range}  -Do -I # -G0.2/0.8
    
    # if global, use -R-180/180/-90/90 -JQ25c; if regional, use -R-135/40/20/85 -JQ25c
    # gmt grdimage ${FILE_div} -R-135/40/20/85 -JQ25c -Baf -B+t"${title}" --FONT_ANNOT_PRIMARY=20p,Helvetica-Bold,black
    gmt grdimage ${FILE_div} 

    gmt plot $FN_RIDGE -Sc0.1c -Gred -Wred

    #### draw colorbar ####
    # gmt colorbar  -Bx+l"divergence rate (1e-8/yr)" -DJBC+w80\%/2\%+e -W1e8 --FONT_ANNOT_PRIMARY=10p,Helvetica-Bold,black
    interval=$(awk -v cpt_range=${cpt_range} 'BEGIN{print cpt_range*1e8/5}')
    gmt colorbar -Ba${interval} -Bx+l"Divergence (1e-8/yr)" -DJBC+w9c/0.25c+e+h -W1e8 --FONT_ANNOT_PRIMARY=12p,Helvetica,black

    #####
    # gmt coast -W1/thin,100 -W2/0p,200 -Dc -A50000 #-Wthinnest,black
    gmt coast -Wthin,black -Dc -A50000 #-Wthinnest,black

    #################################### calculate horizontal velocity (magnitude and degree)
    #### downsample
    gmt grdsample ${output_file}_east_mm.grd -G${output_file}_east_mm_downsample.grd -I2/2
    gmt grdsample ${output_file}_north_mm.grd -G${output_file}_north_mm_downsample.grd -I2/2
    local FILE_east=${output_file}_east_mm_downsample.grd
    local File_north=${output_file}_north_mm_downsample.grd
    gmt grdmath $FILE_east $File_north HYPOT = ${output_file}_vecmag.grd
    # gmt grdmath $File_north $FILE_east ATAN2D = ${output_file}_vecazimuth.grd
    gmt grdmath $File_north $FILE_east ATAN2D = ${output_file}_vecdeg.grd

    #### to xyz
    gmt grd2xyz ${output_file}_vecmag.grd > ${output_file}_vecmag.xyz
    # gmt grd2xyz ${output_file}_vecazimuth.grd > ${output_file}_vecazimuth.xyz
    gmt grd2xyz ${output_file}_vecdeg.grd > ${output_file}_vecdeg.xyz
    #### paste deg and mag
    awk '{print $3}' ${output_file}_vecmag.xyz > temp1.txt
    paste -d' ' ${output_file}_vecdeg.xyz temp1.txt > ${output_file}_vecdeg_mag.xyz
    rm temp1.txt

    #################### plot vector, using -Sv (degree and magnitude) to keep the direction
    ### method:
        ### use a 'scale' to scale both vector field and legend vector.
        ### for example, for scale=0.25, 1mm/yr is 1*scale cm on plot
        ### for the legend, the length of vector is fixed as LEGEND_LENGTH (e.g, 0.5 cm)
        ### then, calculate the corresponding vector magnitude as label
        ### that is, label = LEGEND_LENGTH / scale (mm/yr), will be used to label the legend of vector.

    gmt plot ${output_file}_vecdeg_mag.xyz -Sv0.3c+e+v${scale}c+n0.5c -W1.1p,red -Gred

    # plot vector legend, 1mm/yr is 1*scale cm on plot, so 1cm is 1/scale mm/yr
    # here plot 1cm vector as legend (scale), and magnitude it represents depends on scale: 1/scale mm/yr
    # 0 0 0 1 v0.2c+e+vl0.5c : 0 0 0 1 are starting point (0,0) and direction and polar length; However -Sv0.2c+e+vl0.5c, if +vl is prepended then it is taken as a fixed length to override input lengths.
    gmt inset begin -X7.5c  -DjBL+w2.5c/0.5c -F+gwhite # -Y0.1c
    # gmt basemap -X0.3c -Y0.3c  #-F+gblack

    ### BEGIN legend vector
    LEGEND_LENGTH=0.5
    gmt plot -Baf -Btblr -R0/2/-0.004/0.004 -JX2.5c/0.5c -W1.1p,red -S -Gred << EOF
0.1 0 0 1 v0.3c+e+vl${LEGEND_LENGTH}c
EOF
    # magnitude of velocity the a LEGEND_LENGTH cm vector represents: LEGEND_LENGTH/scale mm/yr
    LEGEND_unit=$(awk -v scale=${scale} -v LEGEND_LENGTH=${LEGEND_LENGTH} 'BEGIN{print LEGEND_LENGTH/scale}')
    gmt text -F+f10p --FONT_ANNOT_PRIMARY=10p,Helvetica,black << EOF
1.3 0 ${LEGEND_unit} mm/yr
EOF
    gmt inset end

    ##### END of vector legend

    ################################## add time label for the figure, if title is not NONE
    if [ "$title" != "NONE" ]; then
        gmt inset begin -X-7.5c  -DjBL+w2.1c/0.5c -F+gwhite
        gmt text -Baf -Btblr -R0/2/-0.004/0.004 -JX2.1c/0.5c -F+f10p --FONT_ANNOT_PRIMARY=10p,Helvetica,black << EOF
1 0 ${title}
EOF
        gmt inset end
    fi


    # echo saved ${output_file}_NA.png


    # # draw colorbar separately
    # # echo "cpt_range is:" $cpt_range
    # interval=$(awk -v cpt_range=${cpt_range} 'BEGIN{print cpt_range*1e8/5}')
    # # echo "interval is:" $interval
    # gmt begin ${output_file}_colorbar png
    #     # gmt makecpt -Chaxby -T-${cpt_range}/${cpt_range}  -Do -G0.2/0.8  #roma -G-0.8/0.8 -I
    #     gmt makecpt -Croma -T-${cpt_range}/${cpt_range}  -Do -G-0.7/0.7 -I
    #     # gmt makecpt -Cvik -T-${cpt_range}/${cpt_range}  -Do -G-0.5/0.5
    #     # gmt makecpt -Chawaii -T-${cpt_range}/${cpt_range}  -Do -I 
    #     # gmt makecpt -C123/50/148,194/165/207,166/219/160,0/136/55 -T-${cpt_range}/${cpt_range}  -Do
    #     gmt colorbar -Ba${interval} -By+l"Divergence (1e-8/yr)" -D+w18c/0.5c+e+h -W1e8 --FONT_ANNOT_PRIMARY=16p,Helvetica,black
    # gmt end

    # rm ${output_file}_uplift_mm.grd 
    rm ${output_file}_east_mm.grd ${output_file}_north_mm.grd
    rm ${output_file}_east_mm_downsample.grd ${output_file}_north_mm_downsample.grd
    rm ${output_file}_vecdeg.grd ${output_file}_vecmag.grd
    rm ${output_file}_vecdeg.xyz ${output_file}_vecmag.xyz
    rm ${output_file}_vecdeg_mag.xyz
}


plot_one_step(){
    case_id=$1 #case6
    step=$2 #648
    title="$3" #"test"
    scale=$4 #0.05
    cpt_range=$5 #5e-8
    PROJECT_BASE=../../../../
    FILE_div=${PROJECT_BASE}/NCAR_DATA_v4/download_from_ncar/dilitation_rate_regulargrid/${case_id}/${case_id}.map_dilitation.${step}.regular.CmCorrect.grd
    FILE_south=${PROJECT_BASE}/NCAR_DATA_v4/download_from_ncar/velocity_regulargrid_allsteps/${case_id}_rm_CM/${case_id}.map_incr_south.${step}.regular.grd
    FILE_east=${PROJECT_BASE}/NCAR_DATA_v4/download_from_ncar/velocity_regulargrid_allsteps/${case_id}_rm_CM/${case_id}.map_incr_east.${step}.regular.grd
    File_north=${PROJECT_BASE}/NCAR_DATA_v4/download_from_ncar/velocity_regulargrid_allsteps/${case_id}_rm_CM/${case_id}.map_incr_north.${step}.regular.grd

    Output_file=./temp/${case_id}_div_vel_map_${step}
    mkdir -p ./temp
    create_north_from_south $FILE_south $File_north
    main $case_id $step $FILE_div $FILE_east $File_north $Output_file "$title" $scale $cpt_range

    # copy files to data/main_fig1
    # cp $FILE_div $FILE_east $File_south $File_north ./data/main_fig1/

}


plot() {

cb_scale=3e-8
NCOL=5

gmt begin ${output_file}_NA_smaller png

ipanel=0
for ((step=648; step>=536; step-=8)); do
# for ((step=648; step>=600; step-=8)); do
    gmt basemap -R-135/-110/20/60 -JQ10c -Baf  --FONT_ANNOT_PRIMARY=10p,Helvetica,black
    plot_one_step case10_3Dcrust $step $step 0.1 1e-8
    ipanel=$((ipanel+1))
    # if if ipanel % NCOL == 0, shift y and x
    if ((ipanel % NCOL == 0)); then
        gmt plot /dev/null -X-50c -Y-20c #-R-135/-110/40/60 -JQ20c -Baf --FONT_ANNOT_PRIMARY=16p,Helvetica,black
    else
        gmt plot /dev/null -X12.5c #-R-135/-110/40/60 -JQ20c -Baf --FONT_ANNOT_PRIMARY=16p,Helvetica,black
    fi

done

gmt end


}

plot

exit







