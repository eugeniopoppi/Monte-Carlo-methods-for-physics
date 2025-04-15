#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<complex.h>
#include"/home/eugenio/Documenti/c/header/random.h"
#define PI 3.14159265358979323846
#define Nx 60
#define Nt 15

int metropolis(double complex lattice[Nx][Nt], int x, int t, double complex R){
  double complex Ut,Uold,Unew,S;
  int acc = 0;

  Uold = lattice[x][t];
  Ut = R*Uold;
  //
  S = lattice

  return acc;
}


int main(){
  double complex lattice[Nx][Nt],R;
  double r,rho = 0.2,theta = PI;
  int n,x,t,acc = 0,numbofmet = 0,iter;
  const unsigned long int seed1=(unsigned long int) time(NULL);
  const unsigned long int seed2=seed1+127;
  FILE *fp;


//initialize the lattice
for(i=0;i<Nx;i++){
  for(j=0;j<Nt;j++){
    lattice[i][j] = 1+0*I;
  }
}

//iter updates of the entire lattice 1/5 micro/overelaxation
  for(n=0;n<iter;n++){
    for(x=0;x<Ns;x++){
      for(t=0;t<Nt;t++){
        r = myrand();
        if(r<rho){
          numbofmet+=1;
          r = (myrand()*2-1)*theta;
          R = (1+ r*I)/sqrt(1+theta*theta);
          acc+=metropolis(lattice,x,t,R);
          w=wilson();

          }
        microcanonic(lattice);
        w=wilson();

      }
    }
  }




  return EXIT_SUCCESS;
}
