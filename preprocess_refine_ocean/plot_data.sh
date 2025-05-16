# module load ncarenv-basic/23.10  gcc/12.2.0 gmt/6.5.0
# module load gmt/6.5.0
# ascii data
# file=./ICE6G_1x1_processed/ocn180x360.1


for epoch in $(seq 0 1 80); do
    echo "epoch: $epoch"
    file=./ICE6G_1x1_refine_ocean_iter2/ocn.${epoch}
    # file_ice=./ICE6G_1x1_processed_iter2/ice180x360.${epoch}
    file_ice=../preprocess_v2/ICE6G_original/ice.${epoch}
    gmt begin $file png
        gmt basemap -R0/360/-90/90 -JQ180/9i -Bxa60f30g30 -Bya30f15g15 -BWSne
        # gmt makecpt -Cgray -T0/1/0.5
        # gmt contour $file -C -W0.5p -I

        gmt makecpt -Ccool -T0/3000/500
        gmt plot $file_ice -Sc0.1c -C
        gmt colorbar -DJBC+w10c/0.5c+o0c/1c 

        # ocean contour at 0.5
        gmt xyz2grd $file -Gtmp.grd -I1 -R0.5/359.5/-89.5/89.5
        gmt grdcontour tmp.grd -A0.5, -W0.5p
        # gmt coast -W0.5p
    gmt end

    rm tmp.grd
done

echo "output: $file.png"