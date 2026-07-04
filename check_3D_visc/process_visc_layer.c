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

 main(argc,argv)
  int argc;
  char **argv;
 {
 char filename[250],outputfile[250],inputfile[250],input_s[300];
 FILE *fp0,*fp1,*fp2,*fp3,*fp4,*fp5;
 float time_s,hf_s,rad1,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,timea[10001],qt[10001],qb[10001];
 float scale_length,scale_potential,viscosity,rr[200],tt[200],visc[200];

 int nproc,i,ip,ipp,j,nstep,compute_node,nno,nox,noy;

 scale_length = 6371000;
 scale_potential = 4.0*M_PI*6.6742e-11*4400*scale_length*scale_length;
 float scale_visc = 1e21;
 
  nproc = atoi(argv[4]);    // number of the processes for the surface
  int nprocz = atoi(argv[5]);
  nox = noy = atoi(argv[6]); //case0AA:25; case0AB:65; case0A/case0C/case0E:41 
  int noz = atoi(argv[7]);
  int depth_id = atoi(argv[8]); //depth slice id
  depth_id = noz - depth_id;  // depth_id is the id of the depth slice from the surface, starting from 1
  nno = nox*noy;
  
  nstep=atoi(argv[2]);

  sprintf(filename,"%s.map_visc.%d.%d",argv[3],noz - depth_id,nstep);  //TODO: depth information needed
  fp2=fopen(filename,"w");


  for (ipp=1;ipp<=nproc;ipp++)  {
    ip = nprocz*ipp-1;  // only top processors
    sprintf(filename,"%s.coord_s.%d",argv[1],ip);
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
       for(j=0;j<depth_id;j++)  // skip the data
          fgets(input_s,300,fp1);
        if (j!=depth_id)
        {
          fprintf(stderr,"Error: depth_id is not correct!\n");
          exit(1);
        }
       fgets(input_s,300,fp1);
       sscanf(input_s,"%f",&temp3);
       for(j=depth_id+1;j<noz;j++)  // skip the data
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
