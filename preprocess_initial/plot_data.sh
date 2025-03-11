# module load ncarenv-basic/23.09  gcc/12.2.0
# module load gmt/6.5.0
# ascii data
# file=./ICE6G_1x1_processed/ocn180x360.1
file=./ICE6G_1x1_122_iter1/ocn180x360.1

gmt begin $file png
    gmt basemap -R0/360/-90/90 -JQ180/9i -Bxa60f30g30 -Bya30f15g15 -BWSne
    gmt makecpt -Cgray -T0/1/0.1
    gmt contour $file -C -W0.5p -I
    gmt coast -W0.5p
gmt end

echo "output: $file.png"