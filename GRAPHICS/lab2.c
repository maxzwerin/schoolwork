#include "FPToolkit.c"
#include "M2d_matrix_tools.c"


#define MAXOBJS 20
#define MAXPTS 1000
#define MAXPOLYS 10000

int numobjects ;
int numpoints[MAXOBJS] ;
double x[MAXOBJS][MAXPTS],y[MAXOBJS][MAXPTS] ;
int numpolys[MAXOBJS] ;
int psize[MAXOBJS][MAXPOLYS] ;
int con[MAXOBJS][MAXPOLYS][20] ;
double red[MAXOBJS][MAXPOLYS], grn[MAXOBJS][MAXPOLYS], blu[MAXOBJS][MAXPOLYS] ;



int read_object(int onum, char **fname) {

  FILE *f;
  f = fopen(fname[onum+1], "r") ;
  if (f == NULL) {
    printf("File not found.\n") ;
    exit(1) ;
  }
  
  fscanf(f, "%d", &numpoints[onum]) ;
  if (numpoints[onum] > MAXPTS) {
    printf("Too many points.\n") ;
    exit(1) ;
  }
  for (int i=0 ; i<numpoints[onum] ; i++) {
    fscanf(f, "%lf %lf", &x[onum][i], &y[onum][i]) ;
  }
  
  fscanf(f, "%d", &numpolys[onum]) ;
  if (numpolys[onum] > MAXPOLYS) {
    printf("Too many polygons.\n") ;
    exit(1) ;
  }
  for (int i=0 ; i<numpolys[onum] ; i++) {
    fscanf(f, "%d", &psize[onum][i]) ;
    for (int j=0 ; j<psize[onum][i] ; j++) {
      fscanf(f, "%d", &con[onum][i][j]) ;
    }
  }
  
  for (int i=0 ; i<numpolys[onum] ; i++) {
    fscanf(f, "%lf", &red[onum][i]) ;
    fscanf(f, "%lf", &grn[onum][i]) ;
    fscanf(f, "%lf", &blu[onum][i]) ;
  }
  
}



int draw_object(int onum) {

  G_rgb(1,1,1);
  G_clear();

  double xp[20], yp[20];
  int p;

  for (int i=0 ; i<numpolys[onum] ; i++) {
    for (int j=0 ; j<psize[onum][i] ; j++) {
      p = con[onum][i][j] ;
      xp[j] = x[onum][p] ;
      yp[j] = y[onum][p] ;
    }

    G_rgb(red[onum][i], grn[onum][i], blu[onum][i]) ;
    G_fill_polygon(xp, yp, psize[onum][i]) ;

  }

}



int rotate_object(int onum, double degrees) {

  double t = degrees*M_PI/180;
  double A[3][3], B[3][3], C[3][3];
  M2d_make_translation (A, -400, -400);
  M2d_make_rotation    (B, t);
  M2d_make_translation (C, 400, 400);
  
  M2d_mat_mult(B, B,A); M2d_mat_mult(C, C,B);
  M2d_mat_mult_points(x[onum], y[onum], C, x[onum], y[onum], numpoints[onum]);
  
}



int main(int argc, char **argv) {


  // scan object(s) in
  numobjects = argc-1;
  if (numobjects > MAXOBJS) {
    printf("Too many objects.\n") ;
    exit(1) ;
  }
  printf("displaying %d objects :\n", numobjects - 1);
  for (int i = 1; i < numobjects; i++) {
    printf("  [%d] %s\n", i - 1, argv[i]);
  }
  for (int onum=0 ; onum<numobjects ; onum++) {read_object(onum, argv);}
  
  
  // shift n scale
  for (int onum=0 ; onum<numobjects ; onum++) {
    
    double minx=x[onum][0], maxx=x[onum][0], miny=y[onum][0], maxy=y[onum][0];
    for (int i=1 ; i<numpoints[onum] ; i++) {
      if (minx > x[onum][i]) minx = x[onum][i];
      if (maxx < x[onum][i]) maxx = x[onum][i];
      if (miny > y[onum][i]) miny = y[onum][i];
      if (maxy < y[onum][i]) maxy = y[onum][i];
    }
    double centerx = (minx + maxx)/2, centery = (miny + maxy)/2;
    double width = fabs(minx - maxx), height = fabs(miny - maxy);
    double square = (width > height) ? width : height;
    double mag = 700/square;
    
    
    // matrix mults
    double A[3][3], B[3][3], C[3][3];
    M2d_make_translation (A, -centerx, -centery);
    M2d_make_scaling     (B, mag, mag);
    M2d_make_translation (C, 400, 400);
    
    // final mult happens once
    M2d_mat_mult(B, B,A); M2d_mat_mult(C, C,B); 
    M2d_mat_mult_points(x[onum], y[onum], C, x[onum], y[onum], numpoints[onum]);
    
  }

  G_init_graphics(800,800) ;
  
  // draw object(s)
  int onum=0, key ;
  while (1) {
    if (48 <= key && key < numobjects+48) onum = key - 48 ;
    if (key == 32) rotate_object(onum, 5) ;
    if (key == 113) exit(1) ;
    draw_object(onum) ;
    key = G_wait_key() ;
  }

}
