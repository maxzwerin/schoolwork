#include "FPToolkit.c"
#include "M3d_matrix_tools.c"


#define MAXOBJS 10
#define MAXPTS 100000
#define MAXPOLYS 100000

// XYZ data
int numobjects ;
int numpoints[MAXOBJS] ;
double x[MAXOBJS][MAXPTS],y[MAXOBJS][MAXPTS],z[MAXOBJS][MAXPTS] ;
int numpolys[MAXOBJS] ;
int psize[MAXOBJS][MAXPOLYS] ;
int con[MAXOBJS][MAXPOLYS][20] ;

// light stuff
double lx, ly, lz ;
double ambient, diffuse, specular, specpow ;

double colors[MAXOBJS][3] = {
    {1.0, 0.2, 0.2},   // red
    {0.2, 1.0, 1.0},   // cyan
    {0.2, 1.0, 0.2},   // green
    {0.3, 0.6, 0.9},   // light blue
    {1.0, 0.2, 1.0},   // magenta
    {0.2, 0.2, 1.0},   // blue
    {1.0, 1.0, 0.2},   // yellow
    {0.5, 0.5, 0.5},   // gray
    {0.8, 0.4, 0.2},   // brown
    {0.9, 0.3, 0.7}    // pink
} ;

typedef
struct {
  double a ; // Nx
  double b ; // Ny
  double c ; // Nz
  double d ;
}
PLANE ; // ax + by + cz + d = 0

// clip stuff
double halfangle ;
double hither, yon ;
PLANE vv[6] ;	// view volume





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
    fscanf(f, "%lf %lf %lf", &x[onum][i], &y[onum][i], &z[onum][i]) ;
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
  
  
}




int center_object(int onum) {

  double minx=x[onum][0], maxx=x[onum][0], 
         miny=y[onum][0], maxy=y[onum][0],
         minz=z[onum][0], maxz=z[onum][0];
    for (int i=1 ; i<numpoints[onum] ; i++) {
      minx = (minx < x[onum][i]) ? minx : x[onum][i];
      maxx = (maxx > x[onum][i]) ? maxx : x[onum][i];
      miny = (miny < y[onum][i]) ? miny : y[onum][i];
      maxy = (maxy > y[onum][i]) ? maxy : y[onum][i];
      minz = (minz < z[onum][i]) ? minz : z[onum][i];
      maxz = (maxz > z[onum][i]) ? maxz : z[onum][i];
    }
  double X = (minx + maxx)/2, Y = (miny + maxy)/2, Z = (minz + maxz)/2 ;

  double M[4][4];
  M3d_make_translation(M, -X, -Y, -Z);
  M3d_mat_mult_points(x[onum], y[onum], z[onum], M, x[onum], y[onum], z[onum], numpoints[onum]);
  
  x[onum][numpoints[onum]] = 0;
  y[onum][numpoints[onum]] = 0;
  z[onum][numpoints[onum]] = 0;
  
}





int translate_object(int onum, double distance, char axis, double sign) {
   
   double M[4][4];

   // x
   if (axis == 'x') M3d_make_translation(M, distance*sign,0,0);
   // y
   if (axis == 'y') M3d_make_translation(M, 0,distance*sign,0);
   // z
   if (axis == 'z') M3d_make_translation(M, 0,0,distance*sign);
   
   M3d_mat_mult_points(x[onum], y[onum], z[onum], M, x[onum], y[onum], z[onum], numpoints[onum]+1);
   
}





int rotate_object(int onum, double degrees, char axis, double sign) {

  double t = degrees*M_PI/180;
  int np = numpoints[onum];
  double A[4][4], B[4][4], C[4][4];
  
  // x
  if (axis == 'x') {
    M3d_make_translation   (A, 0, -y[onum][np], -z[onum][np]);
    M3d_make_x_rotation_cs (B, cos(sign*t), sin(sign*t));
    M3d_make_translation   (C, 0,  y[onum][np],  z[onum][np]);
  }
  
  // y
  if (axis == 'y') {
    M3d_make_translation   (A, -x[onum][np], 0, -z[onum][np]);
    M3d_make_y_rotation_cs (B, cos(sign*t), sin(sign*t));
    M3d_make_translation   (C,  x[onum][np], 0,  z[onum][np]);
  }
  
  // z
  if (axis == 'z') {
    M3d_make_translation   (A, -x[onum][np], -y[onum][np], 0);
    M3d_make_z_rotation_cs (B, cos(sign*t), sin(sign*t));
    M3d_make_translation   (C,  x[onum][np],  y[onum][np], 0);
  }
  
  M3d_mat_mult(A, B,A); M3d_mat_mult(A, C,A);
  M3d_mat_mult_points(x[onum], y[onum], z[onum], A, x[onum], y[onum], z[onum], numpoints[onum]+1);
  
}



int move_light_source(char axis, double sign)
{

  if (axis == 'x') lx += 5*sign ;
  if (axis == 'y') ly += 5*sign ;
  if (axis == 'z') lz += 5*sign ;

}




typedef
struct {
  int objnum ;
  int polynum ;
  double dist ;
  double cx ;
  double cy ;
  double cz ;
}
THING ;

int compare(const void *p, const void *q)
{
  THING *a, *b ;

  a = (THING*)p ;
  b = (THING*)q ;

  if  (((*a).dist) < ((*b).dist)) return -1 ;
  else if (((*a).dist) > ((*b).dist)) return 1 ;
  else return 0 ;
}



int unitize(double *v, int n)
{
  double mag = 0 ;
  for (int i=0 ; i<n ; i++) mag += v[i]*v[i] ;
  mag = sqrt(mag) ;
  for (int i=0 ; i<n ; i++) v[i] /= mag;
}

double light_model(THING thing)
{

  // get 2 points from polygon
  double a[3];
  a[0] = x[thing.objnum][con[thing.objnum][thing.polynum][0]];
  a[1] = y[thing.objnum][con[thing.objnum][thing.polynum][0]];
  a[2] = z[thing.objnum][con[thing.objnum][thing.polynum][0]];
  
  double b[3];
  b[0] = x[thing.objnum][con[thing.objnum][thing.polynum][1]];
  b[1] = y[thing.objnum][con[thing.objnum][thing.polynum][1]];
  b[2] = z[thing.objnum][con[thing.objnum][thing.polynum][1]];
  
  // get center
  double c[3];
  c[0] = thing.cx ;
  c[1] = thing.cy ;
  c[2] = thing.cz ;
  
  // find two vectors on the polygon
  double AC[3], BC[3];
  for (int i=0 ; i<3 ; i++) {
    AC[i] = c[i] - a[i];
    BC[i] = c[i] - b[i];
  }
  
  // find normal vector
  double N[3];
  M3d_x_product(N, AC, BC);
  unitize(N,3);
  
  // find light vector
  double L[3];
  L[0] = c[0] - lx ;
  L[1] = c[1] - ly ;
  L[2] = c[2] - lz ;
  unitize(L,3);
  
  // find eye vector
  double E[3];
  for (int i=0 ; i<3 ; i++) E[i] = -c[i] ;
  unitize(E,3);
  
  // find N dot L (= cos(theta))
  double N_L = 0 ;
  for (int i=0 ; i<3 ; i++) N_L += N[i]*L[i] ;
  
  // find reflect vector
  double R[3];
  for (int i=0 ; i<3 ; i++) R[i] = (2*N_L)*N[i] - L[i] ;
  
  // find E dot R (= cos(beta))
  double E_R = 0 ;
  for (int i=0 ; i<3 ; i++) E_R += E[i]*R[i] ;
  
  // find N dot E (for detecting degeneracy)
  double N_E = 0 ;
  for (int i=0 ; i<3 ; i++) N_E += N[i]*E[i] ;
  
  // find light intensity (with degeneracies)
  double intensity ;
  if (N_L/fabs(N_L) != N_E/fabs(N_E)) {
    intensity = ambient ;
  } else if (N_L < 0 && N_E < 0) {
    N_L = -N_L ;
    E_R = -E_R ;
    intensity = ambient + diffuse*N_L + specular*pow(E_R,specpow) ;
  } else {
    intensity = ambient + diffuse*N_L + specular*pow(E_R,specpow) ;
  }
  return intensity ;

}



double intersect_line_plane(double P[3], double Q[3], int v, double intersection[3])
{
  
  // setup
  double a = vv[v].a ; double b = vv[v].b ; double c = vv[v].c ; double d = vv[v].d ;
  
  // find intersection using parametric form of a line
  double t = - (a*P[0] + b*P[1] + c*P[2] + d) / (a*(Q[0] - P[0]) + b*(Q[1] - P[1]) + c*(Q[2] - P[2])) ;
  for (int i=0 ; i<3 ; i++) intersection[i] = P[i] + t*(Q[i] - P[i]) ;
  
}



double clip_to_volume(THING thing, double *X, double *Y, double *Z)
{

  // init temporary polygon
  int i = thing.objnum  ;	// i = object#
  int j = thing.polynum ;	// j = polygon# 
  
  double xp[100], yp[100], zp[100] ;
  int np = psize[i][j] ;
  for (int k=0 ; k<np ; k++) {	// k = point#
    
    xp[k] = x[i][con[i][j][k]] ;
    yp[k] = y[i][con[i][j][k]] ;
    zp[k] = z[i][con[i][j][k]] ;
    
  }

  int N ;
  int h ;
  
  // plane by plane in the view volume // v = plane#
  for (int v=0 ; v<6 ; v++) {
  
    // line by line in the polygon // k = line/point#
    N = 0 ;
    for (int k=0 ; k<np ; k++) {
    
      // get line endpoints
      double P[3] = {xp[k],   yp[k],   zp[k]} ;
      h = k+1 ; if (h == np) h = 0 ;
      double Q[3] = {xp[h], yp[h], zp[h]} ;
      
      // are the points in or out? // in <= 0 < out
      double Pside = vv[v].a*P[0] + vv[v].b*P[1] + vv[v].c*P[2] + vv[v].d ;
      double Qside = vv[v].a*Q[0] + vv[v].b*Q[1] + vv[v].c*Q[2] + vv[v].d ; 
      
      // find intersection
      double intersection[3] ;
      intersect_line_plane(P, Q, v, intersection) ;
      
      // in to in: keep point
      if (Pside <= 0 && Qside <= 0) {
        X[N] = xp[h] ; Y[N] = yp[h] ; Z[N] = zp[h] ; N++ ;
      }
      
      // in to out: keep point and add intersection
      if (Pside <= 0 && 0 < Qside) {
        X[N] = intersection[0] ; Y[N] = intersection[1] ; Z[N] = intersection[2] ; N++ ;
      }
      
      // out to out: do nothing
      if (0 < Pside && 0 < Qside) {
        
      }
      
      // out to in: add intersection
      if (0 < Pside && Qside <= 0) {
        X[N] = intersection[0] ; Y[N] = intersection[1] ; Z[N] = intersection[2] ; N++ ;
        X[N] = xp[h] ; Y[N] = yp[h] ; Z[N] = zp[h] ; N++ ;

      }
      
    } // end for k
    
    // REPLACE OLD POLYGON WITH NEW POLYGON
    for (int k=0 ; k<N ; k++) {
      xp[k] = X[k] ; yp[k] = Y[k] ; zp[k] = Z[k] ;
    }
    
    np = N ; 
  } // end for v
  
  return N ;

}



// PAINTERS ALGORYTHM
int draw_all_objects()
{
    // make polygon soup
    THING thing[99999];
    int t = 0;
    for (int i = 0; i < numobjects; i++) {
        for (int j = 0; j < numpolys[i]; j++) {
            double cx = 0, cy = 0, cz = 0;
            for (int k = 0; k < psize[i][j]; k++) {
                int p = con[i][j][k];
                cx += x[i][p];
                cy += y[i][p];
                cz += z[i][p];
            }
            cx /= psize[i][j];
            cy /= psize[i][j];
            cz /= psize[i][j];
            double distance = sqrt(cx * cx + cy * cy + cz * cz);

            thing[t].objnum = i;
            thing[t].polynum = j;
            thing[t].dist = distance;
            thing[t].cx = cx;
            thing[t].cy = cy;
            thing[t].cz = cz;
            t++;
        }
    }

    // sort polygon soup
    qsort(thing, t, sizeof(THING), compare);

    // light, clip, project, and draw
    double H = tan(halfangle);
    for (t--; t >= 0; t--) {
        int objnum = thing[t].objnum;
        int polynum = thing[t].polynum;

        // light calculation
        double light = light_model(thing[t]);

        // apply color from array
        double r = colors[objnum % MAXOBJS][0] * light;
        double g = colors[objnum % MAXOBJS][1] * light;
        double b = colors[objnum % MAXOBJS][2] * light;

        // clip
        double x1[100], y1[100], z1[100];
        int N = clip_to_volume(thing[t], x1, y1, z1);

        // project and draw
        if (N > 0) {
            double x2[100], y2[100];
            for (int k = 0; k < N; k++) {
                x2[k] = (400 / H) * (x1[k] / z1[k]) + 400;
                y2[k] = (400 / H) * (y1[k] / z1[k]) + 400;
            }

            G_rgb(r, g, b);
            G_fill_polygon(x2, y2, N);
        }
    }
}






	////////////////////////////////////////////////////////////////////////////////





// 0, 1, 2... to choose an object
// x, y, z to specify an axis
// t to translate
// r to rotate
// l to move light source
// c to change the sign
// -, = to zoom in or out
// [, ] to move hither clipping plane
// q to quit

int main(int argc, char **argv) {

  // light stuff
  lx = 100 ; ly = 200 ; lz = 0 ;
  ambient  = 0.2 ;
  diffuse  = 0.4 ;
  specular = 1 - ambient - diffuse ;
  specpow  = 50 ;
  
  // clip stuff
  halfangle = M_PI/4 ; double H = tan(halfangle);
  hither = 0.1 ; yon = 100 ;



  // init window
  G_init_graphics(800,800);
  G_rgb(0,0,0);
  G_clear();
  
  // display objects
  numobjects = argc-1;
  if (numobjects > MAXOBJS) {
    printf("Too many objects.\n") ;
    exit(1) ;
  }
  printf("Displaying %d objects.\n\n", numobjects);
  for (int onum=0 ; onum<numobjects ; onum++) {
    read_object(onum, argv);
    center_object(onum);
    translate_object(onum, 10, 'z', 1);
  }
  
  // init variables
  double zoom = 45 ;
  int  sign = 1 ;
  int  action = 't' ;  
  int  onum = 0 ; // onum marks the current object
  int  q,k ;
  double V[4][4] ;
  
  
  // CONTROLS
  while (1) {

  for (int i=0 ; i<6 ; i++) {
    vv[i].a = 0 ; vv[i].b = 0 ; vv[i].c = 0 ; vv[i].d = 0 ;
  }
					// plane:
  vv[0].c = -1 ; vv[0].d = hither ;     // hither  0x + 0y - 1z + 0.1 = 0
  vv[1].c =  1 ; vv[1].d = -yon   ;     // yon     0x + 0y + 1z - 100 = 0
  vv[2].b =  1 ; vv[2].c = -H     ;     // top     0x + 1y - Hz + 0 = 0
  vv[3].b = -1 ; vv[3].c = -H     ;     // bottom  0x - 1y - Hz + 0 = 0
  vv[4].a =  1 ; vv[4].c = -H     ;     // right   1x + 0y - Hz + 0 = 0
  vv[5].a = -1 ; vv[5].c = -H     ;     // left   -1x + 0y - Hz + 0 = 0
  
  


    
    
    G_rgb(0,0,0) ; 
    G_clear() ;
    draw_all_objects() ;

    q = G_wait_key() ;
    
    // quit
    if (q == 'q') exit(0) ;
    
    // change sign
    else if (q == 'c') sign = -sign ;

    // translate mode
    else if (q == 't') action = q ;
    
    // rotate mode
    else if (q == 'r') action = q ;
    
    // light mode
    else if (q == 'l') action = q ;

    // choose object
    else if (('0' <= q) && (q <= '9')) {
      k = q - '0' ;  
      if (k < numobjects) { onum = k ; printf("object %d\n", onum) ;}
      else printf("not an object\n");
    }
    
    // translate & rotate
    else if ((q == 'x') || (q == 'y') || (q == 'z')) {
      if (action == 't') translate_object(onum, .5, q, sign);
      if (action == 'r') rotate_object(onum, 5, q, sign);
      if (action == 'l') move_light_source(q, sign);
    }
    
    // zoom in
    else if (q == '=') {
      zoom-- ;
      halfangle = zoom*M_PI/180 ;
    }
    
    // zoom out
    else if (q == '-') {
      zoom++ ;
      halfangle = zoom*M_PI/180 ;
    }
    
    // hither in
    else if (q == '[') {
      if (0.1 < hither) hither -= 0.1 ;
    }
    
    // hither out
    else if (q == ']') {
      if (hither < yon) hither += 0.1 ;
    }
    
    // nothing

  }



}