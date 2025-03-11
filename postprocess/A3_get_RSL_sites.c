/*
input parameters:
fn_RSL_sites=RSL_sites.info  # the file that contains the lon, lat and name of sites
    - first line is the number of sites
fn_RSL_out=RSL_sites.out  # the output file that contains the RSL for sites
fn_RSL_c => get timestep, time, RSL_c
prefix_data for u and g.
nlat, nlon: number of latitudes and longitudes in the regular grid. The regular grid is assumed to be from 0 to 360 for lon and 90 to -90 for lat
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>

int DEBUG = 0;

int get_number_of_lines(const char * filename, int n_skip);
void interp_u_g_to_sites(char * fn_u, char * fn_g, int nlon, int nlat, float * sites_lon, float * sites_lat, int n_sites, float * RSL_array, double RSL_c);

int main(int argc, char** argv){
    

    if (argc != 7){
        printf("Usage: %s fn_RSL_sites fn_RSL_out fn_RSL_c prefix_data nlat nlon\n", argv[0]);
        return 1;
    }

    // read in parameters

    char *fn_RSL_sites = argv[1];
    char *fn_RSL_out = argv[2];
    char *fn_RSL_c = argv[3];
    char *prefix_data = argv[4];
    int nlat = atoi(argv[5]);
    int nlon = atoi(argv[6]);

    // read in sites

    FILE * fp_sites = fopen(fn_RSL_sites, "r");
    if(fp_sites == NULL){ printf("Error: cannot open file %s\n", fn_RSL_sites); return 1; }

    int n_sites;
    char buffer[300];
    fgets(buffer, 300, fp_sites);
    sscanf(buffer, "%d", &n_sites);

    printf("Number of sites: %d\n", n_sites);

    float * sites_lon = (float *)malloc(n_sites * sizeof(float));
    float * sites_lat = (float *)malloc(n_sites * sizeof(float));
    
    for(int i = 0; i < n_sites; i++){
        fgets(buffer, 300, fp_sites);
        sscanf(buffer, "%f %f", &sites_lon[i], &sites_lat[i]);
        // transform lon to 0 - 360
        if(sites_lon[i] < 0){
            sites_lon[i] += 360;
        }
        if(sites_lon[i] > 360 || sites_lon[i] < 0 || sites_lat[i] > 90 || sites_lat[i] < -90){
            printf("Error: invalid lon or lat for site %d\n", i);
            return 1;
        }
    }
    fclose(fp_sites);

    // get number of timesteps

    int n_timesteps;
    n_timesteps = get_number_of_lines(fn_RSL_c, 0); // the last parameter is the number of lines to skip

    float * RSL_array = (float *)malloc(n_timesteps * n_sites * sizeof(float));

    // get timestep, time, RSL_c

    int * timestep_arr = (int *)malloc(n_timesteps * sizeof(int));
    float * time_arr = (float *)malloc(n_timesteps * sizeof(float));
    float * RSL_c_arr = (float *)malloc(n_timesteps * sizeof(float));

    FILE * fp_header = fopen(fn_RSL_c, "r");
    if(fp_header == NULL){ printf("Error: cannot open file %s\n", fn_RSL_c); return 1; }

    for(int i = 0; i < n_timesteps; i++){
        fgets(buffer, 300, fp_header);
        sscanf(buffer, "%d %f %f", &timestep_arr[i], &time_arr[i], &RSL_c_arr[i]);
    }
    fclose(fp_header);

    /* convert time */
    for(int i = 0; i < n_timesteps; i++){
        time_arr[i] = (time_arr[i] - time_arr[n_timesteps - 1])/1000; // convert to kyrs
    }


    /*
    Get RSL for sites:
     1. read in u and g file (regular grid)
     2. interpolate u and g to sites
     3. RSL = g - u + RSL_c

     RSL_array is a 2D array with n_timesteps rows and n_sites columns
    */
    for (int step = 0; step < n_timesteps; step++){
        int current_timestep;
        char fn_u[300], fn_g[300];
        current_timestep = timestep_arr[step];
        sprintf(fn_u, "%s.map_uplift.%d.regular.xyz", prefix_data, current_timestep);
        sprintf(fn_g, "%s.map_geoid.%d.regular.xyz", prefix_data, current_timestep);

        double RSL_c = RSL_c_arr[step];

        interp_u_g_to_sites(fn_u, fn_g, nlon, nlat, sites_lon, sites_lat, n_sites, RSL_array + step * n_sites, RSL_c); // RSL_array + step * n_sites is the pointer to the RSL values at the current timestep
    }

    /*
       substract RSL at present day from RSL at all timesteps
    */

    for (int i = 0; i < n_sites; i++){
        double RSL_present = RSL_array[(n_timesteps - 1) * n_sites + i];
        for (int j = 0; j < n_timesteps; j++){
            RSL_array[j * n_sites + i] -= RSL_present;
        }
    }

    /*
        save to file
    */

    FILE * fp_out = fopen(fn_RSL_out, "w");
    if(fp_out == NULL){ printf("Error: cannot open file %s\n", fn_RSL_out); return 1; }

    fprintf(fp_out, "RSL for sites in %s: timestep, year, site1, site2, .... \n", fn_RSL_sites);
    for (int i = 0; i < n_timesteps; i++){
        fprintf(fp_out, "%d %f ", timestep_arr[i], time_arr[i]);
        for (int j = 0; j < n_sites; j++){
            fprintf(fp_out, "%f ", RSL_array[i * n_sites + j]);
        }
        fprintf(fp_out, "\n");
    }

    fclose(fp_out);

    printf("RSL for sites saved to %s\n", fn_RSL_out);

    free(sites_lon);
    free(sites_lat);
    free(RSL_array);
    free(timestep_arr);
    free(time_arr);
    free(RSL_c_arr);

    return 0;

}


/*
Here the interpolation is done by assuming
1. data on regular grid, lon from 0 to 360, lat from 90 to -90. lon changes faster than lat

*/
void interp_u_g_to_sites(char * fn_u, char * fn_g, int nlon, int nlat, float * sites_lon, float * sites_lat, int n_sites, float * RSL_array, double RSL_c)
{
    FILE * fp_u = fopen(fn_u, "r");
    FILE * fp_g = fopen(fn_g, "r");

    if(fp_u == NULL || fp_g == NULL){
        printf("Error: cannot open file %s or %s\n", fn_u, fn_g);
        exit(1);
    }

    // printf("reading u and g from %s and %s\n", fn_u, fn_g);

    float grid_lon[nlon], grid_lat[nlat];
    float grid_u[nlat*nlon], grid_g[nlat*nlon];

    for (int i = 0; i < nlon; i++){
        grid_lon[i] = i * 360.0 / (nlon - 1);
    }

    for (int i = 0; i < nlat; i++){
        grid_lat[i] = 90 - i * 180.0 / (nlat - 1);
    }

    for (int i = 0; i < nlat; i++){
        for (int j = 0; j < nlon; j++){
            float lon, lat;
            fscanf(fp_u, "%f %f %f", &lon, &lat, &grid_u[i*nlon + j]);
            fscanf(fp_g, "%f %f %f", &lon, &lat, &grid_g[i*nlon + j]);
            if( fabs(lon - grid_lon[j]) > 1e-6 || fabs(lat - grid_lat[i]) > 1e-6){
                printf("Error: lon or lat mismatch\n");
                exit(1);
            }
        }
    }

    fclose(fp_u);
    fclose(fp_g);

    /*
    Interpolate u and g to sites
    */
    float dlon = 360.0 / (nlon - 1);
    float dlat = 180.0 / (nlat - 1);

    for (int i=0; i<n_sites; i++){
        float lon = sites_lon[i];
        float lat = sites_lat[i];


        
        int lon_idx = (int)(lon / dlon);    
        int lat_idx = (int)((90 - lat) / dlat);     

        float dist_left = (lon - grid_lon[lon_idx]) / dlon;
        float dist_right = 1 - dist_left;
        float dist_top = fabs(lat - grid_lat[lat_idx]) / dlat;
        float dist_bottom = 1 - dist_top;

        float u = grid_u[lat_idx * nlon + lon_idx] * dist_right * dist_bottom +
                    grid_u[lat_idx * nlon + lon_idx + 1] * dist_left * dist_bottom +
                    grid_u[(lat_idx + 1) * nlon + lon_idx] * dist_right * dist_top +
                    grid_u[(lat_idx + 1) * nlon + lon_idx + 1] * dist_left * dist_top;
        
        float g = grid_g[lat_idx * nlon + lon_idx] * dist_right * dist_bottom +
                    grid_g[lat_idx * nlon + lon_idx + 1] * dist_left * dist_bottom +
                    grid_g[(lat_idx + 1) * nlon + lon_idx] * dist_right * dist_top +
                    grid_g[(lat_idx + 1) * nlon + lon_idx + 1] * dist_left * dist_top;
        
        if(DEBUG){
            printf("site %d:  ", i );
            printf("u %f g %f\n", u, g);
        }

        RSL_array[i] = g - u + RSL_c;
    }

    return;

}
    

int get_number_of_lines(const char * filename, int n_skip){
    FILE * fp = fopen(filename, "r");
    int n_lines = 0;
    char buffer[300];
    while(fgets(buffer, 300, fp) != NULL){
        n_lines++;
    }
    n_lines -= n_skip;
    fclose(fp);
    return n_lines;
}