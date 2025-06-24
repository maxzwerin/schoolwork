#include "FPToolkit.c"
#include "M3d_matrix_tools.c"


#define MAXOBJS 10
#define MAXPTS 100000
#define MAXPOLYS 100000

int numobjects ;
int numpoints[MAXOBJS] ;
double x[MAXOBJS][MAXPTS],y[MAXOBJS][MAXPTS],z[MAXOBJS][MAXPTS] ;
int numpolys[MAXOBJS] ;
int psize[MAXOBJS][MAXPOLYS] ;
int con[MAXOBJS][MAXPOLYS][20] ;

int window_size = 800 ;

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

double light_model(THING thing) {
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
    c[0] = thing.cx;
    c[1] = thing.cy;
    c[2] = thing.cz;
  
    // find two vectors on the polygon
    double AC[3], BC[3];
    for (int i = 0; i < 3; i++) {
        AC[i] = c[i] - a[i];
        BC[i] = c[i] - b[i];
    }
  
    // find normal vector
    double N[3];
    M3d_x_product(N, AC, BC);
    unitize(N, 3);
  
    // find light vector
    double L[3];
    L[0] = c[0] - lx;
    L[1] = c[1] - ly;
    L[2] = c[2] - lz;
    unitize(L, 3);
  
    // find eye vector
    double E[3];
    for (int i = 0; i < 3; i++) E[i] = -c[i];
    unitize(E, 3);
  
    // find N dot L (= cos(theta))
    double N_L = 0;
    for (int i = 0; i < 3; i++) N_L += N[i] * L[i];
  
    // find reflect vector
    double R[3];
    for (int i = 0; i < 3; i++) R[i] = (2 * N_L) * N[i] - L[i];
  
    // find E dot R (= cos(beta))
    double E_R = 0;
    for (int i = 0; i < 3; i++) E_R += E[i] * R[i];
  
    // find N dot E (for detecting degeneracy)
    double N_E = 0;
    for (int i = 0; i < 3; i++) N_E += N[i] * E[i];
  
    // find light intensity (with degeneracies)
    double intensity;
    if (N_L / fabs(N_L) != N_E / fabs(N_E)) {
        intensity = ambient;
    } else if (N_L < 0 && N_E < 0) {
        N_L = -N_L;
        E_R = -E_R;
        intensity = ambient + diffuse * N_L + specular * pow(E_R, specpow);
    } else {
        intensity = ambient + diffuse * N_L + specular * pow(E_R, specpow);
    }

    return intensity;
}






// for 3D clipping: write subroutine that finds intersect between plane and line




// PAINTERS ALGORYTHM
int draw_all_objects(double halfangle) {
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
      
            thing[t].objnum  = i;
            thing[t].polynum = j;
            thing[t].dist    = distance;
            thing[t].cx = cx;
            thing[t].cy = cy;
            thing[t].cz = cz;
            t++;
        }
    }

    // sort polygon soup
    qsort(thing, t, sizeof(THING), compare);

    // project 3d soup to 2d soup and draw
    double H = tan(halfangle);
    for (t--; t >= 0; t--) {
        int i = thing[t].objnum;
        int j = thing[t].polynum;

        // project
        double x2[100], y2[100];
        for (int k = 0; k < psize[i][j]; k++) {
            int p = con[i][j][k];
            x2[k] = (400 / H) * (x[i][p] / z[i][p]) + 400;
            y2[k] = (400 / H) * (y[i][p] / z[i][p]) + 400;
        }

        // light intensity
        double light = light_model(thing[t]);
        
        // get color of object (use the object color defined earlier)
        double r = colors[i][0] * light;
        double g = colors[i][1] * light;
        double b = colors[i][2] * light;

        // apply color and fill the polygon
        G_rgb(r, g, b);
        G_fill_polygon(x2, y2, psize[i][j]);
    }
}

void draw_light() {
    double x, y, z;
    x = 400 - lx ;
    y = 400 - ly ;
    z = 50 * (lz+100)/1000 ;
    G_rgb(1, 1, 0.5) ;
    G_fill_circle(x, y, z) ;
}


void display(char text_disp_object[100], char text_disp_num[100], int action, int start, double zoom) {
    G_rgb(0, 0, 0);
    G_clear();
    draw_all_objects(zoom*M_PI/180) ;
    char text_light[100] ;

    G_rgb(0.1, 0.6, 0.3);
    G_fill_rectangle(10, window_size - 70, window_size / 5, 60);
    G_rgb(1, 1, 1);
    
    if (action == 'r') {
        G_draw_string(text_disp_object, 15, window_size - 24);
        G_draw_string(text_disp_num, 15, window_size - 44);
        G_draw_string("action : rotate", 15, window_size - 64);
    } else if (action == 't') {
        G_draw_string(text_disp_object, 15, window_size - 24);
        G_draw_string(text_disp_num, 15, window_size - 44);
        G_draw_string("action : translate", 15, window_size - 64);
    } else if (action == 'l') {
        draw_light() ;
        G_fill_rectangle(10, window_size - 70, window_size / 5, 60);
        G_rgb(0, 0, 0);
        G_draw_string("shape : light", 15, window_size - 24);
        snprintf(text_light, sizeof(text_light), "location : %.0f, %.0f, %.0f", lx, ly, lz);
        G_draw_string(text_light, 15, window_size - 44);
        G_draw_string("action : x, y, z", 15, window_size - 64);
    }

    if (start) {
        G_rgb(0.3, 0.6, 0.6);
        G_fill_rectangle(10, window_size - 280, window_size / 5, 200);

        G_rgb(1, 1, 1);
        G_draw_string("to use program : ", 15, window_size - 94);
        G_draw_string("[0-9] > change objects", 15, window_size - 114);
        G_draw_string("T > translate", 15, window_size - 134);
        G_draw_string("R > rotate", 15, window_size - 154);
        G_draw_string("L > control light", 15, window_size - 174);
        G_draw_string("X, Y, X > perform action", 15, window_size - 194);
        G_draw_string("C > reverse direction", 15, window_size - 214);
        G_draw_string("+, - > zoom in/out", 15, window_size - 234);
        G_draw_string("M > close this menu", 15, window_size - 254);
        G_draw_string("Q > quit program", 15, window_size - 274);
    } else {
        G_rgb(0.3, 0.6, 0.6);
        G_fill_rectangle(10, window_size - 100, window_size / 5, 20);
        G_rgb(1, 1, 1);
        G_draw_string("press M to see menu", 15, window_size - 94);
    }
}




// 0, 1, 2... to choose an object
// x, y, z to specify an axis
// t to translate
// r to rotate
// c to change the sign
// q to quit

int main(int argc, char **argv) {

  lx = 100 ; ly = 200 ; lz = 5 ;
  ambient  = 0.2 ;
  diffuse  = 0.4 ;
  specular = 1 - ambient - diffuse ;
  specpow  = 50 ;


  G_init_graphics(800,800);

  numobjects = argc-1;
  if (numobjects > MAXOBJS) {
    printf("Too many objects.\n") ;
    exit(1) ;
  }
  for (int onum=0 ; onum<numobjects ; onum++) {
    read_object(onum, argv);
    center_object(onum);
    translate_object(onum, 10, 'z', 1);
  }

  // Initial setup
  double zoom = 45 ;
  int  flip = 1 ;
  int  sign = 1 ;
  int  action = 't' ;  
  int  onum = 0 ;
  int  q, k ;
  double V[4][4] ;
  int start = 1 ;

  char text_disp_object[100];
  char text_disp_num[100];

  snprintf(text_disp_object, sizeof(text_disp_object), "shape : %s", 4 + argv[1 + onum]);
  snprintf(text_disp_num, sizeof(text_disp_num), "location : %i", onum);
  
  display(text_disp_object, text_disp_num, action, start, zoom);

  while (1) {

    q = G_wait_key() ;

    if (q == 'q') exit(0) ;
    else if (q == 'c') sign = -sign ;
    else if (q == 't' || q == 'r' || q == 'l') action = q ;

    else if (('0' <= q) && (q <= '9')) {
      k = q - '0' ;
      if (k < numobjects) onum = k ;
      if (action == 'l') action = 't' ;
    }

    snprintf(text_disp_object, sizeof(text_disp_object), "shape : %s", 4 + argv[1 + onum]);
    snprintf(text_disp_num, sizeof(text_disp_num), "location : %i", onum);

    if (q == 'm') start = (start > 0) ? 0 : 1 ;

    else if ((q == 'x') || (q == 'y') || (q == 'z')) {
      if (action == 't') translate_object(onum, .5, q, sign);
      if (action == 'r') rotate_object(onum, 5, q, sign);
      if (action == 'l') move_light_source(q, sign);
    }
    else if (q == '=') zoom-- ;
    else if (q == '-') zoom++ ;
    
    display(text_disp_object, text_disp_num, action, start, zoom);

  }

}

