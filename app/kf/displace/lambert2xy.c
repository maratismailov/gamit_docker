#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXLINE 100

void calctrans_(float *lon0, float *lat0, float *maxlon, float *minlon,
		float *maxlat,float *minlat, float *p1, float *p2,
		float *lon, float *lat, int *numpnts)
{
  float theta;
  int n;

/*explain which map projection are using and why
  printf("\nTranslating Longitudes and latittudes to x and y coordinates\n");
  printf("using a Lambert comformable conical projection\n");
  printf("and GMT program mapproject to minimize the distortion.\n");
  printf("If you want to use a different projection,\n");
  printf("replace the subroutines lamberttrans and calctrans.\n\n"); */
  printf("Lambert Projection used from long/lat to XY\n");
  
/* find minimum and maximum lons and lats */
  *minlon = *maxlon = lon[0];
  *minlat = *maxlat = lat[0];

  for (n=0; n<*numpnts; n++) {
    *maxlon = (lon[n] > *maxlon) ? lon[n] : *maxlon;
    *maxlat = (lat[n] > *maxlat) ? lat[n] : *maxlat;
    *minlon = (lon[n] < *minlon) ? lon[n] : *minlon;
    *minlat = (lat[n] < *minlat) ? lat[n] : *minlat;
  }
  
/* Set up standard parallels.  The best fit is generally */
/* 15% of the total difference between the upper and lower */
/* latitudes of the area in question */
  theta = (*maxlat - *minlat)*0.135;
  *p1 = *minlat + theta;
  *p2 = *maxlat - theta;
  *lon0 = (*maxlon + *minlon)/2.0;
  *lat0 = (*maxlat + *minlat)/2.0;

/*
  printf("minlon:\t%f\tmaxlon:\t%f\n",*minlon,*maxlon);
  printf("minlat:\t%f\tmaxlat:\t%f\n",*minlat,*maxlat);
  printf("Standard Parallels:\t%f\t%f\n",*p1,*p2);
  printf("Center of projection:\t%f\t%f\n",*lon0,*lat0);
*/

  return;
}

void lamberttrans_(float *lon0, float *lat0, float *maxlon, float *minlon,
		   float *maxlat, float *minlat, float *p1, float *p2,
		   int *numpnts, int *invert,
		   float *a, float *b, float *c, float *d)
{
  FILE *ab_file, *cd_file;
  char *ab_name = "ccab.tmp";
  char *cd_name = "cccd.tmp";
  char line[80];
  char *runstring;
  int n;

  runstring = (char *) malloc(sizeof(char) * 1000);
  ab_file = (FILE *) malloc(sizeof(FILE));
  cd_file = (FILE *) malloc(sizeof(FILE));

  /* write out the longitudes and latitudes */
  if ((ab_file=fopen(ab_name,"w"))==NULL){
    printf("lamberttrans_: could not open file %s\n",ab_name);
    exit(1);
  } 
  
//  printf("Translating longitudes and latitudes:\n");
  for (n=0; n<*numpnts; n++) {
    fprintf(ab_file,"%.4f %.4f\n",a[n],b[n]);
    // printf("\t%.4f\t%.4f\n",a[n],b[n]);
  }
   if (fclose(ab_file)!=0){
     printf("file %s did not close correctly!\n",ab_name);
     exit(3);
   }
  /*put the command string together */
  if (*invert){
    sprintf(runstring,"gmt mapproject %s -Jl%.2f/%.2f/%.2f/%.2f/1:1 "
	    "-R%.2f/%.2f/%.2f/%.2f -I -C  -F > %s",
	    ab_name,*lon0,*lat0,*p1,*p2,
	    *minlon,*maxlon,*minlat,*maxlat,cd_name);
  } else {
    sprintf(runstring,"gmt mapproject %s -Jl%.2f/%.2f/%.2f/%.2f/1:1 "
	    "-R%.2f/%.2f/%.2f/%.2f -C -F > %s",
	    ab_name,*lon0,*lat0,*p1,*p2,
	    *minlon,*maxlon,*minlat,*maxlat,cd_name);
  }
  /* printf("Run %s\n",runstring); */
  system(runstring);

/*read cd files*/
  if ((cd_file=fopen(cd_name,"r"))==NULL)
    exit(2);
  
//  printf("\nx and y coordinates:\n");
  for (n=0; n<*numpnts; n++) {
    fscanf(cd_file,"%f %f",&(c[n]),&(d[n]));
    c[n]/=1000.0;
    d[n]/=1000.0;
    // printf("\t%.4f\t%.4f\n",c[n],d[n]);
  }
  if (fclose(cd_file)!=0){
    printf("file %s did not close correctly!\n",cd_name);
    exit(4);
  }
  return;
}




