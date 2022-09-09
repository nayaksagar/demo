#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

#define N 64

int main(){
  size(1);
  origin (-0.5,-0.5);
  init_grid(N);
  const face vector muc[]={1e-3,1e-3};
  mu=muc;
  DT=0.1;
  CFL=0.8;
  run();
}

u.t[top] = dirichlet(1);

/**
For the other no-slip boundaries this gives */

u.t[bottom] = dirichlet(0);
u.t[left]   = dirichlet(0);
u.t[right]  = dirichlet(0);

#if !MAC
uf.n[left]   = 0;
uf.n[right]  = 0;
uf.n[top]    = 0;
uf.n[bottom] = 0;
#endif

event init (t = 0){
  solid(cs,fs,sq(x)+sq(y)-sq(0.125));
  foreach()
   u.x[]=cs[];
}
scalar un[];

//error file and stopping when attained steady state
event logfile (t += 0.1; i <= 10000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-5)
    return 1; /* stop */
  //fprintf (stderr, "%f %.9f %g\n", t, energy(), du);
}

double dx=1.0/N;
double dy=1.0/N;
//profile printing
event profile(t=end){
  FILE*fp=fopen("velocity.dat","w");
  fprintf(fp,"VARIABLES=X,Y,U,V\n");
  fprintf(fp,"ZONE T=\"0\" i=%d, j=%d ZONETYPE=ORDERED, DATAPACKING=POINT\n",N,N);
  for(int i=1;i<=N;i++){
    for(int j=1;j<=N;j++){
      double x=-0.5+(i-1)*dx;
      double y=-0.5+(j-1)*dy;
      fprintf(fp,"%f\t%f\t%g\t%g\n", x, y, interpolate(u.x,x,y), interpolate(u.y,x,y));
    }
  }
  fclose(fp);
}
