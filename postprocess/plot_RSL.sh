module load ncarenv-basic/23.10  gcc/12.2.0 gmt/6.5.0

RSL_file=RSL_sites_groupA.out.case_test


gmt begin ${RSL_file} png
    # Plot the RSL data
    awk '(NR>1){print $2, $3}' $RSL_file | \
        gmt plot -R-120/0/-200/900 -JX15c/7c -BWSne -Bxa30f10 -Bya200f100 -W0.5p,black 
    
    awk '(NR>1){print $2, $4}' $RSL_file | \
        gmt plot -W0.5p,red
    
gmt end

echo "output: ${RSL_file}.png"