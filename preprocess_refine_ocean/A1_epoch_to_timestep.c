/*
Find time steps for each epoch.
Input parameter: 1) The load_stages_time_file used in SVE inputfile. 2) the filename to write the output.
*/

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char ** argv)
{
    char inputfile[200], outputfile[200];
    sprintf(inputfile, "%s", argv[1]);
    sprintf(outputfile, "%s", argv[2]);
    FILE * fp_in, * fp_out;

    if((fp_in = fopen(inputfile, "r")) == NULL)
    {
        printf("Error! opening file");
        exit(1);
    }

    if((fp_out = fopen(outputfile, "w")) == NULL)
    {
        printf("Error! opening file");
        exit(1);
    }


    /* Read stage file */
    char buffer[200];
    int n_stages;
    int * stage2timestep; // stage2timestep[stage] = timestep
    fgets(buffer, 200, fp_in);
    sscanf(buffer, "%d", &n_stages); // number of stages
    n_stages = n_stages + 1; // add 1 for the initial stage (or the last stage)

    stage2timestep = (int *)malloc(n_stages * sizeof(int));
    stage2timestep[0] = 0; // initial stage
    int cumulative_timesteps = 0;
    int steps_current_stage;
    float temp_float;
    for (int i = 1; i < n_stages; i++){
        fgets(buffer, 200, fp_in);
        sscanf(buffer, "%f %d", &temp_float, &steps_current_stage);
        cumulative_timesteps += steps_current_stage;
        stage2timestep[i] = cumulative_timesteps;
    }

    /* Write to output file */
    fprintf(fp_out, "%d\n", n_stages);
    for (int i = 0; i < n_stages; i++){
        fprintf(fp_out, "%d\t%d\n", i, stage2timestep[i]);
    }

    fclose(fp_in);
    fclose(fp_out);
    free(stage2timestep);

    return 0;
}