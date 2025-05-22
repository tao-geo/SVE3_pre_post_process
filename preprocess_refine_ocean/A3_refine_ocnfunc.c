/*
refine ocean function (and ice load) by using the u,g from the previous calculation.

Input: 



Generated files:
1) grounded ice load and ocean function for each epoch, saved to files.
2) the calculated (paleo-)topography for the initial epoch. (for the next possible iteration)
    (named paleotopo_epoch0.dat)

Notes:

1) order of grid nodes: first lon changes, then lat changes.
    for(int i=0; i<nlat; i++)
        for(int j=0; j<nlon; j++)
            node = i*nlon + j;
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>

/*
    Function prototypes
*/
void construct_variables(char* inputfilename);
void read_grid_file(char * filename, double * griddata);
void read_grid_file_int(char * filename, int * griddata, int nskip);
void read_ice_height(const char * ice_prefix, int epoch, double * ice_height);
double get_surface_integral_int(const int * mask, const double * lat_grid, const double * lon_grid, int nlat, int nlon);
double get_surface_integral(const double * data, const double * lat_grid, const double * lon_grid, int nlat, int nlon);
void find_groundice_ocean_mask(const double *topo, const double *ice_thickness,
                               int *groundIce_mask, int *ocean_function);
void find_topo_ocean_and_total_load(int epoch, double * A00, double * I00);
void write_ice_and_ocean_func();
void write_topo_initialEpoch();

void helper_get_Delta_SL(int epoch, double * Delta_SL);

// Define the struct
struct {
    // Grid parameters
    int nlat;
    int nlon;
    char ice_prefix[250];
    char topo_filename[250];
    int nepoch;
    int current_epoch; // current epoch. To make sure when a function is called with **epoch**, it is consistent with the current epoch.
                        // if not, there will be problems when XX_currentEpoch is used in the function because it is not consistent with the epoch passed as argument.
    
    /* TODO */
    char prevIter_ocean_prefix[250]; // prefix for ocean function files from previous iteration
    char prevIter_ug_prefix[250]; // prefix for u,g files from previous iteration
    char fn_prevIter_topo_initialEpoch[250]; // filename for the topography for the initial epoch used in the previous iteration

    char out_dir[250];
            // output files (directory). files name are ice${nlat}x${nlon}.${epoch} and ocn${nlat}x${nlon}.${epoch}

    // Data arrays, make sure every one is assigned with values (not just allocated memory)
    // naming: initial/current means initial epoch (0) and current epoch (during loop); presentday means present-day (last epoch)
            // ocn mean ocean function
    double *lat_grid;
    double *lon_grid;
    double *topo_presentday;
    double *topo_initial;    
    double *ice_height_current;  // used in finding floating ice
    int *groundIce_mask_current;   // used in finding floating ice

    int * ocn_current;         // saved to files
    int * ocn_initial;    // needed for topo correction
    double * groundIce_height_current;  // saved to files
    double * groundIce_height_initial;  // needed to subtract for following epoches

    // previous iteration data
    double * prevIter_groundIce_height_current;  // used to get delta I.
    double * prevIter_groundIce_height_initial;
    double * prevIter_topo_initial; // used to get delta SL
    
    int * prevIter_ocn_current;
    int * prevIter_ocn_initial;  // ocean function of initial epoch for the previous iteration 
    // double * Delta_SL_presentday;
    // int *ocean_function_presentday; // NOT USED IN THIS CODE
    // int *groundIce_mask_presentday;  // NOT USED IN THIS CODE
    // double *ice_height_presentday;  // NOT USED IN THIS CODE

} all_data;


double RHO_water = 1000.0;
double RHO_ice = 917.4;


/*
-----------------------------------------------------
*/


int main(int argc, char ** argv){

    FILE * fp_log;
    if((fp_log = fopen("log.txt", "w")) == NULL){
        fprintf(stderr, "Error: cannot open file log.txt\n");
        return 1;
    }



    if(argc != 2){
        fprintf(stderr, "Usage: %s inputfile\n", argv[0]);
        return 1;
    }

    

    construct_variables(argv[1]);

    int nlat = all_data.nlat;
    int nlon = all_data.nlon;
    double A00, I00; // ocean area and ice load (surface integral)

    /*
        Need to calculate Delta_SL_presentday. Only need to calculate once.
        Make sure helper_get_Delta_SL() will be status independent (will not change whenever it is called with the same epoch)
    */



    all_data.current_epoch = -1;



    for(int epoch=0; epoch < all_data.nepoch; epoch++){
        all_data.current_epoch = epoch;
        
        /*
           Read ocean function from previous iteration;
           Read ice height (before removing floating ice) using original ice model;
           Construct grounded ice load used in previous iteration;
           prevIter_groundIce_height = ice_height * (1 - ocean_function);
        */

        char prev_ocean_filename[250];
        sprintf(prev_ocean_filename, "%s.%d", all_data.prevIter_ocean_prefix, epoch);
        read_grid_file_int(prev_ocean_filename, all_data.prevIter_ocn_current, 1); 
                // temoprary use ocn_current to store previous iteration ocean function


        // get previous iteration grounded ice load (total ice height)
        read_ice_height(all_data.ice_prefix, epoch, all_data.ice_height_current);
        for(int node=0; node<nlat*nlon; node++){
            all_data.prevIter_groundIce_height_current[node] = 
                    all_data.ice_height_current[node] * (1.0 - all_data.prevIter_ocn_current[node]);
        }

        // save groundIce_height for epoch 0 for future use
        if(epoch==0){   
            for(int node=0; node<nlat*nlon; node++){
                all_data.prevIter_groundIce_height_initial[node] = 
                        all_data.prevIter_groundIce_height_current[node];
                all_data.prevIter_ocn_initial[node] = all_data.prevIter_ocn_current[node];
            }
        }

        /*
            Get new ocean function and grounded ice load for current epoch
            find_topo_ocean_and_total_load() constructs: 
                groundIce_height_current
                ocn_current
                ocn_initial (used for topo correction)
        */

        
        find_topo_ocean_and_total_load(epoch, &A00, &I00);

        printf("ocean area: %f\n", A00/(4.0*M_PI));

        fprintf(fp_log, "epoch: %d, ocean area: %f\n", epoch, A00/(4.0*M_PI));

        // write out ice and ocean function
        write_ice_and_ocean_func();

        if(epoch==0){
            write_topo_initialEpoch();
        }
    }

    fclose(fp_log);
    return 0;
}


/*
    write out files as SVE input:
        1. ice load (grounded: groundIce_height_current) + topo correction. Relative to first epoch.
        2. ocean function
        3. years.reverse

*/

void write_ice_and_ocean_func(){

    int epoch = all_data.current_epoch;
    int nlat = all_data.nlat;
    int nlon = all_data.nlon;
    const double * lat_grid = all_data.lat_grid;
    const double * lon_grid = all_data.lon_grid;

    char ice_filename[250];
    char ocean_filename[250];

    sprintf(ice_filename, "%s/ice.%d", all_data.out_dir, epoch);
    sprintf(ocean_filename, "%s/ocn.%d", all_data.out_dir, epoch);

    printf("Writing to files: %s, %s\n", ice_filename, ocean_filename);

    FILE * fp_icefile;
    FILE * fp_oceanfile;

    if((fp_icefile = fopen(ice_filename, "w")) == NULL){
        fprintf(stderr, "Error: cannot open file %s\n", ice_filename);
        return;
    }

    if((fp_oceanfile = fopen(ocean_filename, "w")) == NULL){
        fprintf(stderr, "Error: cannot open file %s\n", ocean_filename);
        return;
    }
    // write the header (nlon, nlat)
    fprintf(fp_icefile, "%d %d\n", nlon, nlat);
    fprintf(fp_oceanfile, "%d %d\n", nlon, nlat);

    double * topo_correction = (double *) malloc(nlat * nlon * sizeof(double));

    for(int node=0; node<nlat*nlon; node++){
        topo_correction[node] = all_data.topo_initial[node] *
                                (all_data.ocn_current[node] -
                                all_data.ocn_initial[node]) *
                                (RHO_water / RHO_ice);
    }

    // first lon increases, then lat increases
    int node;
    for(int i=0; i<nlat; i++){
        for(int j=0; j<nlon; j++){
            node = i*nlon + j;
            // write ice load
            fprintf(fp_icefile, "%10.4E %10.4E %11.5E\n", lon_grid[j], lat_grid[i], 
                    all_data.groundIce_height_current[node] - topo_correction[node]);

            // write ocean function
            fprintf(fp_oceanfile, "%10.4E %10.4E %d\n", lon_grid[j], lat_grid[i], 
                    all_data.ocn_current[node]);
        }
    }

    // close the files
    fclose(fp_icefile);
    fclose(fp_oceanfile);
    
    // free memory
    free(topo_correction);

    return;
}

void write_topo_initialEpoch(){

    int nlat = all_data.nlat;
    int nlon = all_data.nlon;
    const double * lat_grid = all_data.lat_grid;
    const double * lon_grid = all_data.lon_grid;

    char topo0_filename[250];

    sprintf(topo0_filename, "%s/paleotopo_epoch0.dat", all_data.out_dir);

    printf("Writing to files: %s\n", topo0_filename);

    FILE * fp_topo;

    if((fp_topo = fopen(topo0_filename, "w")) == NULL){
        fprintf(stderr, "Error: cannot open file %s\n", topo0_filename);
        return;
    }

    

    // first lon increases, then lat increases
    int node;
    for(int i=0; i<nlat; i++){
        for(int j=0; j<nlon; j++){
            node = i*nlon + j;
            // write ice load
            fprintf(fp_topo, "%10.4E %10.4E %11.5E\n", lon_grid[j], lat_grid[i], 
                    all_data.topo_initial[node]);
        }
    }

    // close the files
    fclose(fp_topo);

    return;
}




/*
get Delta_SL;

    Delta_SL: geoid - uplift + c
    where c = -1/A00[epoch] * (RHO_ice/RHO_water * I00[epoch] + dyn00[epoch] )

    A00 is ocean area (current epoch)
    I is change of grounded ice load (compared to the initial epoch)
    I00 is surface intergral of I

    dyn = (g - u) * ocean_mask[epoch] - Topo[initial_epoch] * (ocean_mask[epoch] - ocean_mask[initial_epoch])
    dyn00 is surface integral of dyn

A potential problem for this function:
    It is not status independent - here if not only depends on epoch,
    but also **implicitly** depends on prevIter_ocn_current and prevIter_groundIce_height_current,
    which could be inconsistent with the **epoch** passed as argument.

*/

void helper_get_Delta_SL(int epoch, double * Delta_SL){
    int nlat = all_data.nlat;
    int nlon = all_data.nlon;

    if(epoch != all_data.current_epoch){
        fprintf(stderr, "Error: epoch is not consistent with the current epoch but X_currentEpoch array is used\n");
        exit(1);
    }

    /* 
        read g,u from previous iteration 
    */


    char temp_char[250];
    // FILE * fp_u, * fp_g;
    double * g = (double *) malloc(nlat * nlon * sizeof(double));
    double * u = (double *) malloc(nlat * nlon * sizeof(double));

    sprintf(temp_char, "%s.map_geoid.epoch%d.regular.xyz", all_data.prevIter_ug_prefix, epoch);
    read_grid_file(temp_char, g);

    sprintf(temp_char, "%s.map_uplift.epoch%d.regular.xyz", all_data.prevIter_ug_prefix, epoch);
    read_grid_file(temp_char, u);


    double A00, I00, dyn00;
    double c;

    A00 = get_surface_integral_int(all_data.prevIter_ocn_current, 
                                all_data.lat_grid, all_data.lon_grid, 
                                    nlat, nlon);

    I00 = get_surface_integral(all_data.prevIter_groundIce_height_current,
                                all_data.lat_grid, all_data.lon_grid,
                                    nlat, nlon);
    
    I00 -= get_surface_integral(all_data.prevIter_groundIce_height_initial,
                                all_data.lat_grid, all_data.lon_grid,
                                    nlat, nlon);
    
    double * dyn = (double *) malloc(nlat * nlon * sizeof(double));
    
    for(int node=0; node<nlat*nlon; node++){
        dyn[node] = (g[node] - u[node]) * all_data.prevIter_ocn_current[node] -
                    all_data.prevIter_topo_initial[node] * 
                    (all_data.prevIter_ocn_current[node] - all_data.prevIter_ocn_initial[node]);
    }

    dyn00 = get_surface_integral(dyn, all_data.lat_grid, all_data.lon_grid, nlat, nlon);

    c = -1.0/A00 * (RHO_ice/RHO_water * I00 + dyn00);

    for(int node=0; node<nlat*nlon; node++){
        Delta_SL[node] = g[node] - u[node] + c;
    }

    // fclose(fp_g);
    // fclose(fp_u);

    free(g);
    free(u);
    free(dyn);
}


/*
Why cann't I just call helper_get_Delta_SL() for presentday?
    Because helper_get_Delta_SL() needs some variables named with "current", 
    and the loop starts with epoch 0, all the *current* variables (like ice height current)
    are with data for epoch 0, not the data for the present day (last epoch).
    so I can not just call helper_get_Delta_SL(presentday_epoch, Delta_SL_presentday) for present day.

    Unless I start the loop backward, from epoch n-1 to 0, then the first loop will be present day,
    and all the current variables will be with data for present day.
*/
void helper_get_Delta_SL_presentday(double * Delta_SL_presentday)
{
    int nlat = all_data.nlat;
    int nlon = all_data.nlon;

    int presentday_epoch = all_data.nepoch - 1;

    char prev_ocean_filename[250];

    int * prev_ocean_function_presentday = (int *) malloc(nlat * nlon * sizeof(int));
    int * prevIter_ocn_initial = (int *) malloc(nlat * nlon * sizeof(int));
    double * ice_height_presentday = (double *) malloc(nlat * nlon * sizeof(double));
    double * ice_height_initialEpoch = (double *) malloc(nlat * nlon * sizeof(double));
    double * prevIter_groundIce_height_presentday = (double *) malloc(nlat * nlon * sizeof(double));
    double * prevIter_groundIce_height_initial = (double *) malloc(nlat * nlon * sizeof(double));

    sprintf(prev_ocean_filename, "%s.%d", all_data.prevIter_ocean_prefix, presentday_epoch);
    read_grid_file_int(prev_ocean_filename, prev_ocean_function_presentday, 1); 

    sprintf(prev_ocean_filename, "%s.%d", all_data.prevIter_ocean_prefix, 0);
    read_grid_file_int(prev_ocean_filename, prevIter_ocn_initial, 1); 



    // get previous iteration grounded ice load
    read_ice_height(all_data.ice_prefix, 0, ice_height_initialEpoch);
    read_ice_height(all_data.ice_prefix, presentday_epoch, ice_height_presentday);

    for(int node=0; node<nlat*nlon; node++){
        prevIter_groundIce_height_presentday[node] = 
                ice_height_presentday[node] * (1.0 - prev_ocean_function_presentday[node]);
        prevIter_groundIce_height_initial[node] = 
                ice_height_initialEpoch[node] * (1.0 - prevIter_ocn_initial[node]);
    }

    
    char temp_char[250];
    // FILE * fp_u, * fp_g;
    double * g = (double *) malloc(nlat * nlon * sizeof(double));
    double * u = (double *) malloc(nlat * nlon * sizeof(double));

    sprintf(temp_char, "%s.map_geoid.epoch%d.regular.xyz", all_data.prevIter_ug_prefix, presentday_epoch);
    read_grid_file(temp_char, g);

    sprintf(temp_char, "%s.map_uplift.epoch%d.regular.xyz", all_data.prevIter_ug_prefix, presentday_epoch);
    read_grid_file(temp_char, u);

    double A00, I00, dyn00;
    double c;

    A00 = get_surface_integral_int(prev_ocean_function_presentday, 
                                all_data.lat_grid, all_data.lon_grid, 
                                    nlat, nlon);

    I00 = get_surface_integral(prevIter_groundIce_height_presentday,
                                all_data.lat_grid, all_data.lon_grid,
                                    nlat, nlon);
    
    I00 -= get_surface_integral(prevIter_groundIce_height_initial,
                                all_data.lat_grid, all_data.lon_grid,
                                    nlat, nlon);
    
    double * dyn = (double *) malloc(nlat * nlon * sizeof(double));
    
    for(int node=0; node<nlat*nlon; node++){
        dyn[node] = (g[node] - u[node]) * prev_ocean_function_presentday[node] -
                    all_data.prevIter_topo_initial[node] * 
                    (prev_ocean_function_presentday[node] - prevIter_ocn_initial[node]);
    }

    dyn00 = get_surface_integral(dyn, all_data.lat_grid, all_data.lon_grid, nlat, nlon);

    c = -1.0/A00 * (RHO_ice/RHO_water * I00 + dyn00);

    for(int node=0; node<nlat*nlon; node++){
        Delta_SL_presentday[node] = g[node] - u[node] + c;
    }

    // fclose(fp_g);
    // fclose(fp_u);

    free(g);
    free(u);
    free(dyn);
    free(prev_ocean_function_presentday);
    free(prevIter_ocn_initial);
    free(ice_height_presentday);
    free(ice_height_initialEpoch);
    free(prevIter_groundIce_height_presentday);
    free(prevIter_groundIce_height_initial);

    return;

}

/*
Purpose:
        Find topo, ocean, and total load for one epoch

Construct: 
    groundIce_height_current
    ocn_current
    ocn_initial
*/
void find_topo_ocean_and_total_load(int epoch, double * A00, double * I00){
    int nlat = all_data.nlat;
    int nlon = all_data.nlon;

    static int been_here = 0;
    static double * Delta_SL_presentday;
    double * Delta_SL_currentEpoch;
    double * topo_currentEpoch;

    if(epoch != all_data.current_epoch){
        fprintf(stderr, "Error: epoch is not consistent with the current epoch\n");
        exit(1);
    }

    if(been_here == 0){
        Delta_SL_presentday = (double *) malloc(nlat * nlon * sizeof(double));
        
        // calculate Delta_SL_presentday
        helper_get_Delta_SL_presentday(Delta_SL_presentday);
        been_here = 1;
    }


    Delta_SL_currentEpoch = (double *) malloc(nlat * nlon * sizeof(double));
    topo_currentEpoch = (double *) malloc(nlat * nlon * sizeof(double));

    helper_get_Delta_SL(epoch, Delta_SL_currentEpoch);

    /*
        calculate topo

        topo = presentday_topo + Delta_SL_presentday - Delta_SL_currentEpoch
    */

    for(int node=0; node<nlat*nlon; node++){
        topo_currentEpoch[node] = all_data.topo_presentday[node] +
                                    Delta_SL_presentday[node] - Delta_SL_currentEpoch[node];
    }

    if (epoch == 0) {
        for (int node = 0; node < nlat * nlon; node++) {
            all_data.topo_initial[node] = topo_currentEpoch[node];
        }
    }

    if (epoch == all_data.nepoch - 1) {
        for (int node=0; node<nlat*nlon; node++){
            assert(fabs(Delta_SL_presentday[node] - Delta_SL_currentEpoch[node]) < 1e-6);
        }
    }

    // find ground ice and ocean mask
    find_groundice_ocean_mask(topo_currentEpoch, all_data.ice_height_current,
                                all_data.groundIce_mask_current, all_data.ocn_current);

    // find ground ice height
    for(int node=0; node<nlat*nlon; node++){
        all_data.groundIce_height_current[node] = all_data.ice_height_current[node] * 
                                                        all_data.groundIce_mask_current[node];
    }

    // subtract initial epoch -> get relative change in ground ice height
    // groundIce_height_current now is relative to the first epoch
    if (epoch == 0) {
        for (int node = 0; node < nlat * nlon; node++) {
            all_data.groundIce_height_initial[node] =
                all_data.groundIce_height_current[node];
            
            all_data.groundIce_height_current[node] = 0.0;

            all_data.ocn_initial[node] = all_data.ocn_current[node];
        }
    } else {
        for (int node = 0; node < nlat * nlon; node++) {
            // subtract initial epoch
            all_data.groundIce_height_current[node] =
                all_data.groundIce_height_current[node] -
                all_data.groundIce_height_initial[node];
        }
    }

    // double I00, A00;
    I00[0] = get_surface_integral(all_data.groundIce_height_current,
                                all_data.lat_grid, all_data.lon_grid,
                                    nlat, nlon);
    
    // or *A00
    A00[0] = get_surface_integral_int(all_data.ocn_current,
                                    all_data.lat_grid, all_data.lon_grid,
                                        nlat, nlon);

    free(Delta_SL_currentEpoch);
    free(topo_currentEpoch);

    return;

}



/*
Check and remove floating ice

Reading in: topo, ice thickness (before checking floating ice)
Write out: ground ice mask, ocean mask

Note, the topo is the current topo (relative to current sea level), NOT (necessarily) present-day topo.
*/
void find_groundice_ocean_mask(const double *topo, const double *ice_thickness,
                               int *groundIce_mask, int *ocean_function) 
{
    int nlat = all_data.nlat;
    int nlon = all_data.nlon;
    // int land_mask[nlat * nlon];

    // int * land_mask = (int *) malloc(nlat * nlon * sizeof(int));
    // int * floating_ice_mask = (int *) malloc(nlat * nlon * sizeof(int));

    int land_mask = 0;
    int floating_ice_mask = 0;
    int ice_mask = 0;
    

    // check floating ice: floating ice: ice where topo < 0 and ice thickness < abs(topo) * RHO_water / RHO_ice

    for (int i=0; i<nlat*nlon; i++){
        land_mask = topo[i] > 0 ? 1 : 0; // 1 for land, 0 for ocean
        ice_mask = ice_thickness[i] > 0 ? 1 : 0; // 1 for ice, 0 for no ice

        floating_ice_mask = 0;
        if( (land_mask == 0) && ice_mask && 
            ( ice_thickness[i] < (fabs(topo[i]) * RHO_water / RHO_ice) ))
            {
                floating_ice_mask = 1;
            }
        
        // ground ice mask
        groundIce_mask[i] = 0;
        if((ice_mask == 1) && (floating_ice_mask == 0)){
            groundIce_mask[i] = 1;
        }

        // ocean function
        ocean_function[i] = 0;
        if( (land_mask == 0) && (groundIce_mask[i] == 0) ) {
            ocean_function[i] = 1;
        }
    }

    // free(land_mask);
    // free(floating_ice_mask);
    return;
}


// integral on unit sphere.
double get_surface_integral(const double * data, const double * lat_grid, const double * lon_grid, int nlat, int nlon){
    double integral = 0.0;
    double dlat = M_PI / nlat;
    double dlon = 2.0 * M_PI/ nlon;
    for(int i=0; i<nlat; i++){
        for(int j=0; j<nlon; j++){
            // area = cos(lat) * dlat * dlon
            integral += cos(lat_grid[i] * M_PI / 180.0) * dlat * dlon * data[i*nlon + j];
        }
    }
    return integral;
}

// this is for integer data
double get_surface_integral_int(const int * mask, const double * lat_grid, const double * lon_grid, int nlat, int nlon){
    double integral = 0.0;
    double dlat = M_PI / nlat;
    double dlon = 2.0 * M_PI/ nlon;
    for(int i=0; i<nlat; i++){
        for(int j=0; j<nlon; j++){
            // area = cos(lat) * dlat * dlon
            integral += cos(lat_grid[i] * M_PI / 180.0) * dlat * dlon * mask[i*nlon + j];
        }
    }
    return integral;
}




/*
Construct non-changing variables:

- input parameters from input file

    // First line: prefix for ice model files, for example: ICE6G_1x1_122/ice6g180x360
    // Second line: nlat, nlon: two numbers for grid: 180 360
    // third line: nepoch: number of epoches: for example 122, then it reads file ICE6G_1x1_122/ice6g180x360.X, where X is from 0 to 121 (including).
    // forth line: filename,   present-day topography file,  
    // fifth line: output directory
    // sixth line: prefix for ocean function files from previous iteration
    // seventh line: prefix for u,g files from previous iteration (on regular grid, ascii format)
    // eighth line: filename for the topography for the initial epoch used in the previous iteration

- lat, lon grid
- present-day topo
- present-day ice height
- allocate memory for other arrays
*/
void construct_variables(char* inputfilename){

    FILE * fp_inputfile; // file pointer for input file
    int nlat, nlon, nepoch;
    char buffer[250];

    // open file
    if((fp_inputfile = fopen(inputfilename, "r")) == NULL){
        fprintf(stderr, "Error: cannot open file %s\n", inputfilename);
        return;
    }

    printf("Reading input file: %s\n", inputfilename);

    // read parameters from input file
    // First line: prefix for ice model files, for example: ICE6G_1x1_122/ice6g180x360
    // Second line: nlat, nlon: two numbers for grid: 180 360
    // third line: nepoch: number of epoches: for example 122, then it reads file ICE6G_1x1_122/ice6g180x360.X, where X is from 0 to 121 (including).
    // forth line: filename,   present-day topography file,  

    // read the lines: printout for checking, instead of using return values to check
    fgets(buffer, 250, fp_inputfile);
    sscanf(buffer, "%s", all_data.ice_prefix);
    printf("ice_prefix: %s\n", all_data.ice_prefix);

    fgets(buffer, 250, fp_inputfile);
    sscanf(buffer, "%d %d", &nlat, &nlon);
    printf("nlat: %d, nlon: %d\n", nlat, nlon);

    fgets(buffer, 250, fp_inputfile);
    sscanf(buffer, "%d", &nepoch);
    printf("nepoch: %d\n", nepoch);

    fgets(buffer, 250, fp_inputfile);
    sscanf(buffer, "%s", all_data.topo_filename);
    printf("topo_filename: %s\n", all_data.topo_filename);

    // get the output file prefix for ice
    fgets(buffer, 250, fp_inputfile);
    sscanf(buffer, "%s", all_data.out_dir);
    printf("out_dir: %s\n", all_data.out_dir);

    // sixth line: prefix for ocean function files from previous iteration
    fgets(buffer, 250, fp_inputfile);
    sscanf(buffer, "%s", all_data.prevIter_ocean_prefix);
    printf("prevIter_ocean_prefix: %s\n", all_data.prevIter_ocean_prefix);

    // seventh line: prefix for u,g files from previous iteration
    fgets(buffer, 250, fp_inputfile);
    sscanf(buffer, "%s", all_data.prevIter_ug_prefix);
    printf("prevIter_ug_prefix: %s\n", all_data.prevIter_ug_prefix);

    // eighth line: filename for the topography for the initial epoch used in the previous iteration
    fgets(buffer, 250, fp_inputfile);
    sscanf(buffer, "%s", all_data.fn_prevIter_topo_initialEpoch);
    printf("fn_prevIter_topo_initialEpoch: %s\n", all_data.fn_prevIter_topo_initialEpoch);

    // close
    fclose(fp_inputfile);

    // check if out_dir exists, if not, create it
    // struct stat st = {0};
    // if (stat(all_data.out_dir, &st) == -1) {

    int err = mkdir(all_data.out_dir, 0744);
    if (err && errno != EEXIST) {
        fprintf(stderr, "Error: cannot create directory %s\n", all_data.out_dir);
        return;
    }

    // }

    // store
    all_data.nlat = nlat;
    all_data.nlon = nlon;
    all_data.nepoch = nepoch;


    // construct arrays
    all_data.lat_grid = (double *) malloc(nlat * sizeof(double));
    all_data.lon_grid = (double *) malloc(nlon * sizeof(double));
    all_data.topo_presentday = (double *) malloc(nlat * nlon * sizeof(double));
    // all_data.ice_height_presentday = (double *) malloc(nlat * nlon * sizeof(double));
    all_data.ice_height_current = (double *) malloc(nlat * nlon * sizeof(double));
    // all_data.groundIce_mask_presentday = (int *) malloc(nlat * nlon * sizeof(int));
    all_data.groundIce_mask_current = (int *) malloc(nlat * nlon * sizeof(int));
    // all_data.ocean_function_presentday = (int *) malloc(nlat * nlon * sizeof(int));
    all_data.ocn_current = (int *) malloc(nlat * nlon * sizeof(int));
    all_data.groundIce_height_current = (double *) malloc(nlat * nlon * sizeof(double));
    all_data.groundIce_height_initial = (double *) malloc(nlat * nlon * sizeof(double));
    all_data.ocn_initial = (int *) malloc(nlat * nlon * sizeof(int));
    all_data.prevIter_groundIce_height_current = (double *) malloc(nlat * nlon * sizeof(double));
    all_data.prevIter_groundIce_height_initial = (double *) malloc(nlat * nlon * sizeof(double));
    all_data.topo_initial = (double *) malloc(nlat * nlon * sizeof(double));
    all_data.prevIter_topo_initial = (double *) malloc(nlat * nlon * sizeof(double));
    all_data.prevIter_ocn_current = (int *) malloc(nlat * nlon * sizeof(int));
    all_data.prevIter_ocn_initial = (int *) malloc(nlat * nlon * sizeof(int));

    double dlat = 180.0 / nlat; // latitude grid size
    double dlon = 360.0 / nlon; // longitude grid size

    // create the grids
    for (int i=0; i<nlat; i++){
        all_data.lat_grid[i] = 90.0 - dlat/2 - i*dlat;
    }

    for (int i=0; i<nlon; i++){
        all_data.lon_grid[i] = dlon/2 + i*dlon;
    }

    // print out the grids
    // printf("lat_grid: ");
    // for (int i=0; i<nlat; i++){
    //     printf("%f ", lat_grid[i]);
    // }
    // printf("\n");
    // printf("lon_grid: ");
    // for (int i=0; i<nlon; i++){
    //     printf("%f ", lon_grid[i]);
    // }

    // load_present_topo
    read_grid_file(all_data.topo_filename, all_data.topo_presentday);

    // load the initial epoch topo used in the previous iteration
    read_grid_file(all_data.fn_prevIter_topo_initialEpoch, all_data.prevIter_topo_initial);

    // get the present-day ice model data (why? This will be used in every other epoch)
    // all_data.ice_height_presentday already allocated
    // read_ice_height(all_data.ice_prefix, nepoch-1, all_data.ice_height_presentday);

    return;
}

/*
Read grid files, including ice and topo.
The grid data should be (lon, lat, data) format, where lon is increasing (between 0 to 360), lat is decreasing (between 90 to -90).
The grid order is, first lon changes, then lat changes.
*/
void read_grid_file(char * filename, double * griddata){
    FILE * fp_gridfile; // file pointer for grid file
    char buffer[250];
    double lon, lat, data;
    int node;
    if((fp_gridfile = fopen(filename, "r")) == NULL){
        fprintf(stderr, "Error: cannot open file %s\n", filename);
        exit(1);
        // return;
    }

    int nlat = all_data.nlat;
    int nlon = all_data.nlon;
    const double * lat_grid = all_data.lat_grid;
    const double * lon_grid = all_data.lon_grid;

    // read the grid data
    for(int i=0; i<nlat; i++){
        for(int j=0; j<nlon; j++){
            node = i*nlon + j;
            fgets(buffer, 250, fp_gridfile);
            if(sscanf(buffer, "%lf %lf %lf", &lon, &lat, &data) != 3){
                fprintf(stderr, "Error: reading grid file %s\n", filename);
                return;
            }
            // check lon and lat
            if(lon != lon_grid[j] || lat != lat_grid[i]){
                fprintf(stderr, "Error: reading grid file %s, lon or lat does not match\n", filename);
                fprintf(stderr, "expected lon: %f, lat: %f\n", lon_grid[j], lat_grid[i]);
                fprintf(stderr, "read lon: %f, lat: %f\n", lon, lat);
                return;
            }
            griddata[node] = data;
        }
    }

    fclose(fp_gridfile);
}

/* for int type */
void read_grid_file_int(char * filename, int * griddata, int nskips){
    FILE * fp_gridfile; // file pointer for grid file
    char buffer[250];
    double lon, lat;
    int data;
    int node;
    if((fp_gridfile = fopen(filename, "r")) == NULL){
        fprintf(stderr, "Error: cannot open file %s\n", filename);
        exit(1);
        // return;
    }

    int nlat = all_data.nlat;
    int nlon = all_data.nlon;
    const double * lat_grid = all_data.lat_grid;
    const double * lon_grid = all_data.lon_grid;

    if (nskips > 0){
        for(int i=0; i<nskips; i++){
            fgets(buffer, 250, fp_gridfile);
        }
    }

    // read the grid data
    for(int i=0; i<nlat; i++){
        for(int j=0; j<nlon; j++){
            node = i*nlon + j;
            fgets(buffer, 250, fp_gridfile);
            if(sscanf(buffer, "%lf %lf %d", &lon, &lat, &data) != 3){
                fprintf(stderr, "Error: reading grid file %s\n", filename);
                return;
            }
            // check lon and lat
            if(lon != lon_grid[j] || lat != lat_grid[i]){
                fprintf(stderr, "Error: reading grid file %s, lon or lat does not match\n", filename);
                fprintf(stderr, "expected lon: %f, lat: %f\n", lon_grid[j], lat_grid[i]);
                fprintf(stderr, "read lon: %f, lat: %f\n", lon, lat);
                return;
            }
            griddata[node] = data;
        }
    }

    fclose(fp_gridfile);
}


/*
get ice height for certain epoch, a wrapper function for reading ice model data
input:
- ice_prefix: prefix for ice model files, for example: ICE6G_1x1_122/ice6g180x360
- epoch: epoch number, from 0 to nepoch-1
- array to write ice height data
- nlat, nlon: number of lat, number of lon
- lat_grid, lon_grid: lat, lon grid data
*/

void read_ice_height(const char * ice_prefix, int epoch, double * ice_height){

    // first get file name
    char ice_filename[250];
    sprintf(ice_filename, "%s.%d", ice_prefix, epoch);
    printf("Reading ice model file: %s\n", ice_filename);

    // read ice model data
    read_grid_file(ice_filename, ice_height);

    return;

}