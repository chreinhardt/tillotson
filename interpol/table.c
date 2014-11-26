#include	<math.h>#include	<stdio.h>#include	<stdlib.h>#include	"coeff.h"#include	"interpol.h"float brent(float (*func)(float,float,float), float a, float b, float rho, float u, float tol); int Ntable,SplineDegree=3;float *raster,vmax,rhomax;void init(){ FILE *fil;  fil = fopen("lookup.txt","r");  fscanf(fil,"%i %e %e",&Ntable,&rhomax,&vmax);  int i,j;  raster = malloc(Ntable*Ntable*sizeof(float));  for (i=0; i<Ntable; i++)    for (j=0; j<Ntable; j++)      fscanf(fil,"%e",raster+Ntable*j+i);  fclose(fil);  SamplesToCoefficients(raster, Ntable, Ntable, SplineDegree);}float uenergy(float v, float rho){ float iv,irho,u;  iv = Ntable*v/vmax;  irho = Ntable*rho/rhomax;  printf("indices of v,rho %f %f\n",iv,irho);  u = InterpolatedValue(raster, Ntable, Ntable, iv, irho, SplineDegree);  return u;}float denergy(float v, float rho, float u){ return uenergy(v,rho) - u;}float slike(float rho, float u){ float tol=1e-6;  return brent(denergy,0,vmax,rho,u,tol);}int main (int argc, const char *argv[]){ float v,rho,u;  init();  sscanf(argv[1],"%f",&rho);  sscanf(argv[2],"%f",&u);  v = slike(rho,u);  printf("v=%f u=%f\n",v,uenergy(v,rho));}