module load ncarenv-basic/23.10  gcc/12.2.0 gmt/6.5.0
# module load gmt/6.5.0
# ascii data
# file=./ICE6G_1x1_processed/ocn180x360.1

file=./combined_files_case_test/case_test.map_rate_uplift.976.regular.grd
gmt begin $file png
    gmt basemap -R0/360/-90/90 -JQ180/9i -Bxa60f30g30 -Bya30f15g15 -BWSne

    # gmt makecpt -Cwysiwyg -T-1000/1000/50
    # gmt plot $file -Sc0.1c -Cwysiwyg
    gmt grdimage $file -Cwysiwyg
    gmt colorbar -DJBC+w10c/0.5c+o0c/1c  -Bx+l"uplift rate (m/yr)"

    gmt coast -W0.5p -Dc -A1000 
gmt end



echo "output: $file.png"