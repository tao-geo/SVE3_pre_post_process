/*
usage: ./A0_timestep_and_time start_timestep end_timestep prefix outputfile ncpuz
*/

#include <stdio.h>
#include <stdlib.h>

int DEBUG = 0;

int main(int argc, char ** argv)
{
    int start_timestep, end_timestep;
    char SVE_prefix[250];
    int ncpuz;
    char topo_prefix[250], outputfile[250];
    static int been_here = 0;

    /* input parameters */

    if(argc != 6)
    {
        printf("Error! usage: ./A0_timestep_and_time start_timestep end_timestep prefix outputfile ncpuz\n");
        exit(1);
    }
    
    sscanf(argv[1], "%d", &start_timestep);
    sscanf(argv[2], "%d", &end_timestep);
    sscanf(argv[3], "%s", SVE_prefix);
    sscanf(argv[4], "%s", outputfile);
    sscanf(argv[5], "%d", &ncpuz);

    sprintf(topo_prefix, "%s.topo_s.%d", SVE_prefix, ncpuz-1);

    char temp_char[400];

    FILE * fp_topo, * fp_out;

    fp_out = fopen(outputfile, "w"); 
    if(fp_out == NULL) {printf("Error! opening file %s", outputfile); exit(1);}

    /*
        loop over time steps (start_timestep to end_timestep, step 1),
        check if the file exists, if yes, write the time step and time to the output file.
    */
    float time, dt;
    int nsteps=0;
    for(int step = start_timestep; step <= end_timestep; step++)
    {
        sprintf(temp_char, "%s.%d", topo_prefix, step);
        fp_topo = fopen(temp_char, "r");
        if(fp_topo == NULL) { continue;} // file does not exist, skip to next time step


        fgets(temp_char, 400, fp_topo); // first line is header description, only write it once

        if(been_here == 0)
        {
            fprintf(fp_out, "%s", temp_char);
            been_here = 1;
        }

        fgets(temp_char, 400, fp_topo); // the headers

        fprintf(fp_out, "%s", temp_char); 

        fclose(fp_topo);
        nsteps++;
    }

    fclose(fp_out);
    printf("file %s written\n", outputfile);
    printf("number of time steps written: %d\n", nsteps);
    return 0;
}