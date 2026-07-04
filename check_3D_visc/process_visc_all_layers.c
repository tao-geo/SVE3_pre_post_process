/*
load 3D distributed visc datafile, create 2D map files of a depth slice;
assume the depth interested is on top processors.

argv list:
1st: prefix of CitcomSVE: the x in 'x.coord_s.y'
2nd: nstep
3rd: output files prefix
4th: nproc of surface
5th: nprocz
6th: nox (=noy) per process
7th: noz per process
8th: depth slice id, starting from surface, which is 1
*/

#include <math.h>
#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <sys/types.h>
#include <stdlib.h>

void main(int argc, char **argv)
 {
 char filename[250],outputfile[250],inputfile[250],input_s[300];
 FILE *fp0,*fp1,*fp2,*fp3,*fp4,*fp5;
 float time_s,hf_s,rad1,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,timea[10001],qt[10001],qb[10001];
 float viscosity,rr[200],tt[200],visc[200];

 int nproc_surf,i,ip,ipp,j,nstep,compute_node,nno,nox,noy;

 float scale_visc = 1e21;
 
  nproc_surf = atoi(argv[4]);    // number of the processes for the surface
  int nprocz = atoi(argv[5]);
  nox = noy = atoi(argv[6]); 
  int noz = atoi(argv[7]);

  int procz; // loop throught procs in z direction [0.. nprocz-1]
  nno = nox*noy;
  
  nstep=atoi(argv[2]);

  int depth_id_local; // the depth id in local cpu [1, noz]
  int depth_id_global_min; // minimal depth id for a local cpu, could be 1, or noz, ...

  // loop through procz [0, nprocz-1]
  for(procz=0; procz< nprocz; procz++) {
    depth_id_global_min = procz*(noz-1) + 1;  

    // loop through each depths in one procz, depth_id_local order is bottom up (bottom is 0, surface is noz-1)
    for(depth_id_local=0; depth_id_local<noz; depth_id_local++) {  

      sprintf(filename,"%s.map_visc.%d.%d",argv[3],depth_id_local + depth_id_global_min,nstep);  //global depth id for this layer.
      fp2=fopen(filename,"w");

      for (ipp=0;ipp<nproc_surf;ipp++)  {
        ip = nprocz*ipp+procz;  
        sprintf(filename,"%s.coord_s.%d",argv[1], nprocz*ipp + nprocz - 1); // only top processors
        fprintf(stderr,"%s\n",filename);
        fp0=fopen(filename,"r");
        fgets(input_s,300,fp0);

        sprintf(filename,"%s.stress_visc.%d.%d",argv[1],ip,nstep);
        fprintf(stderr,"%s\n",filename);
        fp1=fopen(filename,"r");
        fgets(input_s,300,fp1);  // first line is additional information, ignore it

        rad1 = 180.0/M_PI;

        for (i=1;i<=nno;i++) {
          fgets(input_s,300,fp0);  // read the coordinates
          sscanf(input_s,"%f %f",&temp1,&temp2);
          for(j=0;j<depth_id_local;j++)  // skip the data
              fgets(input_s,300,fp1);
            if (j!=depth_id_local)
            {
              fprintf(stderr,"Error: depth_id is not correct!\n");
              exit(1);
            }
          fgets(input_s,300,fp1);
          sscanf(input_s,"%f",&temp3);
          for(j=depth_id_local+1;j<noz;j++)  // skip the data
              fgets(input_s,300,fp1);
          temp1=90-temp1*rad1;
          temp2=temp2*rad1;
          temp3 = temp3*scale_visc;
          fprintf(fp2,"%g %g %g\n",temp2,temp1,temp3); 
          }
        fclose(fp0);
        fclose(fp1);
        }

      fclose(fp2);
    }
  }
  return;
 }
