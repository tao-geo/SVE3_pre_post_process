/*

input parameters are read from input file.


Change log:

2025-5-17: NEW initial topography estimation
    - before the initial topography was approximated by the present-day topography
    - Now, the correct the initial topography by using 1) isostasy approximation and 2) eustatic sea level change.
    That is, assuming the difference between initial topo and present day topo are from:
        1. isostasy
        2. eustatic sea level change 
    T_end - T_initial = - \Delta I * (RHO_ice/RHO_rock) + sealevel_relative2_presentday
    T_end = topo_presentday
    T_initial = topo_presentday + \Delta I * (RHO_ice/RHO_rock) - sealevel_relative2_presentday
    where \Delta I is the change in ice thickness from initial to present day.


Notes:

1) order of grid nodes: first lon changes, then lat changes.
    for(int i=0; i<nlat; i++)
        for(int j=0; j<nlon; j++)
            node = i*nlon + j;
2) If one array is constructed and used only in one function, 
    I just allocate it in that function (static pointer with been_here)
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>


// Define the struct
struct {
    // Grid parameters
    int nlat;
    int nlon;
    char ice_prefix[250];
    char topo_filename[250];
    int nepoch;

    // output files (directory). files name are ice${nlat}x${nlon}.${epoch} and ocn${nlat}x${nlon}.${epoch}
    char out_dir[250];

    // Data arrays
    double *lat_grid;
    double *lon_grid;
    double *topo_presentday;
    double *topo_initialEpoch; // initial epoch topo
    double *eustatic_sealevel_relative2presentday; // note this is relative to present day (change backward)
    double *ice_height_presentday;
    double *ice_height_currentEpoch;
    int *groundIce_mask_presentday;
    int *groundIce_mask_currentEpoch;
    int *ocean_function_presentday;
    int *ocean_function_currentEpoch;         // saved to files
    double * groundIce_height_currentEpoch;  // saved to files

    double * groundIce_height_initialEpoch;  // needed to subtract for following epoches
    int * ocean_function_initialEpoch;    // needed for topo correction

} all_data;

double RHO_water = 1000.0;
double RHO_ice = 917.4;
double RHO_rock = 3300.0;

// function prototypes
void construct_variables(char* inputfilename);
void read_grid_file(char * filename, double * griddata);
void read_ice_height(const char * ice_prefix, int epoch, double * ice_height);
void find_groundice_ocean_mask(const double *topo, const double *ice_thickness,
                                 int *groundIce_mask, int *ocean_function);

double get_surface_integral(const double * data, const double * lat_grid, const double * lon_grid, int nlat, int nlon);
double get_surface_integral_int(const int * mask, const double * lat_grid, const double * lon_grid, int nlat, int nlon);
double solve_for_ocean_func(int epoch);
void write_ice_ocean_func(int epoch);   

void get_initial_topo();


int main(int argc, char **argv)
{
    FILE * fp_log;
    if((fp_log = fopen("log.txt", "a")) == NULL){
        fprintf(stderr, "Error: cannot open file log.txt\n");
        return 1;
    }
    
    char inputfilename[250];
    sprintf(inputfilename, "%s", argv[1]);

    // construct variables by reading input file
    construct_variables(inputfilename);

    int nlat = all_data.nlat;
    int nlon = all_data.nlon;
    

    // get groundIce_mask_presentday and ocean_function_presentday
    
    find_groundice_ocean_mask(all_data.topo_presentday, all_data.ice_height_presentday, 
                                all_data.groundIce_mask_presentday, 
                                all_data.ocean_function_presentday);
    
    // loop through all epoches
    // get ocean_function_currentEpoch, groundIce_mask_currentEpoch
    // and write to file
    // 
    
    double ocean_area;
 

    for(int epoch=0; epoch<all_data.nepoch; epoch++){

        printf("Calculating epoch: %d\n", epoch);

        // read ice height for current epoch
        read_ice_height(all_data.ice_prefix, epoch, all_data.ice_height_currentEpoch);

        // get groundIce_mask_currentEpoch, ocean_function_currentEpoch
        ocean_area = solve_for_ocean_func(epoch);

        // subtract initial epoch -> get relative change in ground ice height
        // groundIce_height_currentEpoch now is relative to the first epoch
        if (epoch == 0) {
            for (int node = 0; node < nlat * nlon; node++) {
                all_data.groundIce_height_initialEpoch[node] =
                    all_data.groundIce_height_currentEpoch[node];
                
                all_data.groundIce_height_currentEpoch[node] = 0.0;

                all_data.ocean_function_initialEpoch[node] = all_data.ocean_function_currentEpoch[node];
            }
        } else {
            for (int node = 0; node < nlat * nlon; node++) {
                // subtract initial epoch
                all_data.groundIce_height_currentEpoch[node] =
                    all_data.groundIce_height_currentEpoch[node] -
                    all_data.groundIce_height_initialEpoch[node];
            }
        }

        printf("Ocean area: %f\n", ocean_area/(4*M_PI));
        fprintf(fp_log, "epoch: %d, ocean area: %f\n", epoch, ocean_area/(4.0*M_PI));

        // get the initial topography, only for the first epoch
        if(epoch == 0){
            get_initial_topo();
        }
        
        // write to file
        write_ice_ocean_func(epoch);
    }

    fclose(fp_log);
    return 0;
}

/*
get initial topography

Need know the total ice height change from initial to present day,
    IceThickness_presentday_grounded[i] = 
        all_data.ice_height_presentday[i] * all_data.groundIce_mask_presentday[i]
    groundIce_height_initialEpoch: absolute height of ground ice at initial epoch
*/
void get_initial_topo(){
    int nlat = all_data.nlat;
    int nlon = all_data.nlon;
    const double * lat_grid = all_data.lat_grid;
    const double * lon_grid = all_data.lon_grid;
    
    // eustatic correction
    for (int i=0; i<nlat*nlon; i++){
        all_data.topo_initialEpoch[i] = all_data.topo_presentday[i] - all_data.eustatic_sealevel_relative2presentday[0];
    }
    // isostatic correction
    for (int i=0; i<nlat*nlon; i++){
        all_data.topo_initialEpoch[i] = all_data.topo_initialEpoch[i] + (all_data.ice_height_presentday[i] * all_data.groundIce_mask_presentday[i]
        - all_data.groundIce_height_initialEpoch[i]) * (RHO_ice/RHO_rock);
    }

    // write to file
    char topo_filename[250];
    sprintf(topo_filename, "%s/topo_initialEpoch", all_data.out_dir);
    FILE * fp_topofile;
    if((fp_topofile = fopen(topo_filename, "w")) == NULL){
        fprintf(stderr, "Error: cannot open file %s\n", topo_filename);
        return;
    }

    // write (lon, lat, topo) to file
    int node;
    for(int i=0; i<nlat; i++){
        for(int j=0; j<nlon; j++){
            node = i*nlon + j;
            // write ice load
            fprintf(fp_topofile, "%10.4E %10.4E %11.5E\n", lon_grid[j], lat_grid[i], 
                    all_data.topo_initialEpoch[node]);
        }
    }
    fclose(fp_topofile);
    printf("Initial topography written to file: %s\n", topo_filename);
    return;
}

/*
solve_for_ocean_func:
    write to ocean_function_currentEpoch, groundIce_height_currentEpoch
    return ocean area
*/
double solve_for_ocean_func(int epoch){
    
    int nlat = all_data.nlat;
    int nlon = all_data.nlon;
    const double * lat_grid = all_data.lat_grid;
    const double * lon_grid = all_data.lon_grid;

    const double * IceThickness_absolute = all_data.ice_height_currentEpoch;  // current epoch ice height
    const double * IceThickness_presentday = all_data.ice_height_presentday;
    static double * IceThickness_presentday_grounded;

    static int been_here = 0;

    int N_ITER = 3; // number of iterations
    
    // The iteration is the find sealevel_relative2_presentday: the sea level for current epoch
    // Then, use that to find the ocean mask for current epoch (and ground ice mask)

    double sealevel_relative2_presentday = 0.0;
    double ocean_area = 0.0;
    double IceVolume_relative2_presentday = 0.0;
    // double * IceThickness_grounded; 
    static double * IceThickness_grounded_relative2_presentday;  // grounded ice thickness of current epoch relative to present day
                // used to get IceVolume_relative2_presentday
    static int * ocean_function; // ocean mask for current epoch
    static int * groundIce_mask; // ground ice mask for current epoch
    static double * topo_currentEpoch; // topo for current epoch
    static double * topo_correction_grid; // topo correction (grid) for current epoch
    
    double topo_correction_term; // topo correction (volume) for current epoch
    
        
    /*
    */
    if (been_here == 0){
        
        // If this is the first time, calculate IceThickness_presentday_grounded, do not release it.
        IceThickness_presentday_grounded = (double *) malloc(nlat * nlon * sizeof(double));
        
        for (int i=0; i<nlat*nlon; i++){
            IceThickness_presentday_grounded[i] = 
                IceThickness_presentday[i] * all_data.groundIce_mask_presentday[i];
        }

        // other arrays
        // [ ] release memory for this! or only do this in been_here == 0
        // better release memory for these arrays after each epoch, and allocate again for next epoch
        IceThickness_grounded_relative2_presentday = (double *) malloc(nlat * nlon * sizeof(double));
        ocean_function = (int *) malloc(nlat * nlon * sizeof(int));
        groundIce_mask = (int *) malloc(nlat * nlon * sizeof(int));
        topo_currentEpoch = (double *) malloc(nlat * nlon * sizeof(double));
        topo_correction_grid = (double *) malloc(nlat * nlon * sizeof(double));

        been_here = 1;
    }


    /*
        The first initial guess of sealevel_relative2_presentday: ice change / present day ocean area
    */

    for (int i=0; i<nlat*nlon; i++){
        IceThickness_grounded_relative2_presentday[i] = 
            IceThickness_absolute[i] * all_data.groundIce_mask_presentday[i] - 
            IceThickness_presentday_grounded[i];
    }


    IceVolume_relative2_presentday = get_surface_integral(
        IceThickness_grounded_relative2_presentday, lat_grid, lon_grid, nlat, nlon
    );

    ocean_area = get_surface_integral_int(all_data.ocean_function_presentday, lat_grid, lon_grid, nlat, nlon);

    sealevel_relative2_presentday = -1.0 * RHO_ice/RHO_water * IceVolume_relative2_presentday / ocean_area;

    // Done with first initial guess

    /* do following iterations */
    for(int iter=0; iter<N_ITER; iter++){
        // get current topo relative to sealevel, using sealevel_relative2_presentday
        for (int i=0; i<nlat*nlon; i++){
            topo_currentEpoch[i] = all_data.topo_presentday[i] - sealevel_relative2_presentday;
        }
        // get ground ice mask and ocean mask for current epoch
        find_groundice_ocean_mask(topo_currentEpoch, IceThickness_absolute, 
                                    groundIce_mask, ocean_function);
        
        // get grounded ice thickness relative to present day
        for (int i=0; i<nlat*nlon; i++){
            IceThickness_grounded_relative2_presentday[i] = 
                IceThickness_absolute[i] * groundIce_mask[i] - 
                IceThickness_presentday_grounded[i];
        }
        IceVolume_relative2_presentday = get_surface_integral(
            IceThickness_grounded_relative2_presentday, lat_grid, lon_grid, nlat, nlon
        );

        ocean_area = get_surface_integral_int(ocean_function, lat_grid, lon_grid, nlat, nlon);


        // note: correct the wedge area. always is negative
        for (int i=0; i<nlat*nlon; i++){
            topo_correction_grid[i] = all_data.topo_presentday[i] * 
                                        ( all_data.ocean_function_presentday[i] -
                                            ocean_function[i]);
        }  
        topo_correction_term = get_surface_integral(
            topo_correction_grid, lat_grid, lon_grid, nlat, nlon
        );

        // if sea level increase (or decrease), add a positive topo correction term, that's why it's subtracted
        sealevel_relative2_presentday = (-1.0 * RHO_ice/RHO_water * IceVolume_relative2_presentday
                                        - topo_correction_term) / ocean_area;
    }

    // Done with iterations

    all_data.eustatic_sealevel_relative2presentday[epoch] = sealevel_relative2_presentday;

    // get current topo relative to sealevel, using sealevel_relative2_presentday
    for (int i=0; i<nlat*nlon; i++){
        topo_currentEpoch[i] = all_data.topo_presentday[i] - sealevel_relative2_presentday;
    }
    // get ground ice mask and ocean mask for current epoch
    find_groundice_ocean_mask(topo_currentEpoch, IceThickness_absolute, 
                                groundIce_mask, ocean_function);

    // write to arrays 
    for (int i=0; i<nlat*nlon; i++){
        all_data.ocean_function_currentEpoch[i] = ocean_function[i];
        all_data.groundIce_mask_currentEpoch[i] = groundIce_mask[i];
        all_data.groundIce_height_currentEpoch[i] = IceThickness_absolute[i] * groundIce_mask[i];
    }

    ocean_area = get_surface_integral_int(ocean_function, lat_grid, lon_grid, nlat, nlon);
    return ocean_area; 
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
    write out files as SVE input:
        1. ice load (grounded: groundIce_height_currentEpoch) + topo correction. Relative to first epoch.
        2. ocean function
        3. years.reverse

*/

void write_ice_ocean_func(int epoch){

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

    // for topo correction, use the initial epoch topo
    for(int node=0; node<nlat*nlon; node++){
        topo_correction[node] = all_data.topo_initialEpoch[node] *
                                (all_data.ocean_function_currentEpoch[node] -
                                all_data.ocean_function_initialEpoch[node]) *
                                (RHO_water / RHO_ice);
    }

    // first lon increases, then lat increases
    int node;
    for(int i=0; i<nlat; i++){
        for(int j=0; j<nlon; j++){
            node = i*nlon + j;
            // write ice load
            fprintf(fp_icefile, "%10.4E %10.4E %11.5E\n", lon_grid[j], lat_grid[i], 
                    all_data.groundIce_height_currentEpoch[node] - topo_correction[node]);

            // write ocean function
            fprintf(fp_oceanfile, "%10.4E %10.4E %d\n", lon_grid[j], lat_grid[i], 
                    all_data.ocean_function_currentEpoch[node]);
        }
    }

    // close the files
    fclose(fp_icefile);
    fclose(fp_oceanfile);
    
    // free memory
    free(topo_correction);

    return;
}



/*
Construct non-changing variables:
- input parameters from input file
    (ice_prefix, nlat, nlon, nepoch, topo_filename)
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
    all_data.topo_initialEpoch = (double *) malloc(nlat * nlon * sizeof(double));
    all_data.eustatic_sealevel_relative2presentday = (double *) malloc(nepoch * sizeof(double));
    all_data.ice_height_presentday = (double *) malloc(nlat * nlon * sizeof(double));
    all_data.ice_height_currentEpoch = (double *) malloc(nlat * nlon * sizeof(double));
    all_data.groundIce_mask_presentday = (int *) malloc(nlat * nlon * sizeof(int));
    all_data.groundIce_mask_currentEpoch = (int *) malloc(nlat * nlon * sizeof(int));
    all_data.ocean_function_presentday = (int *) malloc(nlat * nlon * sizeof(int));
    all_data.ocean_function_currentEpoch = (int *) malloc(nlat * nlon * sizeof(int));
    all_data.groundIce_height_currentEpoch = (double *) malloc(nlat * nlon * sizeof(double));
    all_data.groundIce_height_initialEpoch = (double *) malloc(nlat * nlon * sizeof(double));
    all_data.ocean_function_initialEpoch = (int *) malloc(nlat * nlon * sizeof(int));

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

    // get the present-day ice model data (why? This will be used in every other epoch)
    // all_data.ice_height_presentday already allocated
    read_ice_height(all_data.ice_prefix, nepoch-1, all_data.ice_height_presentday);

    return;
}


/*
Check and remove floating ice

Reading in: topo, ice thickness (before checking floating ice)
Write out: ground ice mask, ocean mask

Note, the topo is the current topo (relative to current sea level), NOT (necessarily) present-day topo.

TODO: Note, here topo, ice_thickness depends on input (presentday or current epoch)
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
