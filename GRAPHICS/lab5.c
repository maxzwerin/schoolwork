#include "FPToolkit.c"
#include "M3d_matrix_tools.c"  

#define MAXPTS 59000
#define MAXPOLYS 57500
#define MAXOBJS 10

int numobjects ;                   
int numpoints[MAXOBJS] ;  
int numpolys[MAXOBJS] ;  
double x[MAXOBJS][MAXPTS] ; 
double y[MAXOBJS][MAXPTS] ;
double z[MAXOBJS][MAXPTS] ;  
int psize[MAXOBJS][MAXPOLYS] ;    
int con[MAXOBJS][MAXPOLYS][20] ;  

int halfangle = 45 ; 
double window_size = 800 ;

int cull_flag[MAXOBJS] ;

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
};

typedef struct {
    int object ;
    int vertices ; 
    double x[MAXPTS] ;
    double y[MAXPTS] ;
    double zavg ;
} Polygon3D[MAXPOLYS*MAXOBJS] ;

int poly_count = 0 ; 

// Read object data from a file
int read_object(FILE *f, int object) {
    int i, j;

    // Read the number of points for the object
    fscanf(f, "%d", &numpoints[object]);

    // Check if the number of points exceeds the maximum allowed
    if (numpoints[object] >= MAXPTS) {
        printf("MAXPTS = %d : exceeded.\n", MAXPTS);
        exit(1); // Exit if too many points
    }

    // Read point coordinates
    for (i = 0; i < numpoints[object]; i++) {
        fscanf(f, "%lf %lf %lf", &x[object][i], &y[object][i], &z[object][i]);
    }

    // Read the number of polygons
    fscanf(f, "%d", &numpolys[object]);
    
    // Check if the number of polygons exceeds the maximum allowed
    if (numpolys[object] > MAXPOLYS) {
        printf("MAXPOLYS = %d : exceeded.\n", MAXPOLYS);
        exit(1); // Exit if too many polygons
    }

    // Read polygon data
    for (i = 0; i < numpolys[object]; i++) {
        fscanf(f, "%d", &psize[object][i]); // Size of the polygon
        for (j = 0; j < psize[object][i]; j++) {
            fscanf(f, "%d", &con[object][i][j]); // Indices of points for the polygon
        }
    }

    return 1 ; 
}




int draw_3d_polygon(int object, double xp[], double yp[], double zp[], int np) {

    double xbb[100], ybb[100] ;
    double H = tan(halfangle * (M_PI / 180)) ;
    double h_wind = window_size / 2 ; 
    int i ;

    // Project the 3D points into 2D space
    for (i = 0 ; i < np ; i++) {
        xbb[i] = ((h_wind / H) * (xp[i] / zp[i]) + h_wind) ;
        ybb[i] = ((h_wind / H) * (yp[i] / zp[i]) + h_wind) ; 
    }

    G_rgb(colors[object][0], colors[object][1], colors[object][2]);
    G_fill_polygon(xbb, ybb, np) ;
    G_rgb(0,0,0) ; 
    G_polygon(xbb, ybb, np) ; 

    return 1 ; 
}



int cull_polygon(int object, double xp[], double yp[], double zp[], int np) {
    double A[3], B[3], normal[3], eye[3] ;
    double dot_product ;

    double a[3] = {xp[0], yp[0], zp[0]} ;
    double b[3] = {xp[1], yp[1], zp[1]} ;
    double c[3] = {xp[2], yp[2], zp[2]} ;

    A[0] = b[0] - a[0] ;
    A[1] = b[1] - a[1] ;
    A[2] = b[2] - a[2] ;

    B[0] = c[0] - a[0] ;
    B[1] = c[1] - a[1] ;
    B[2] = c[2] - a[2] ;

    M3d_x_product(normal, A,B);

    eye[0] = -a[0] ;
    eye[1] = -a[1] ;
    eye[2] = -a[2] ;

    dot_product = (eye[0] * normal[0]) + (eye[1] * normal[1]) + (eye[2] * normal[2]) ;

    if (cull_flag[object] > 0) return dot_product < 0 ;
    return dot_product > 0 ; 
}




int draw_single_object(int object) {
    int h, i, j; 
    double xp[100], yp[100], zp[100], index[100];
    double zavg = 0 ; 
    int np;
    int count = 0;

    for (i = 0; i < numpolys[object]; i++) {
        np = psize[object][i];

        for (j = 0; j < np; j++) {
            h = con[object][i][j];
            xp[j] = x[object][h];
            yp[j] = y[object][h];
            zp[j] = z[object][h];
        }

        if (cull_polygon(object, xp, yp, zp, np)) continue;  

        draw_3d_polygon(object, xp, yp, zp, np);
        count++ ; 
    }

    return 1;
}





int display(char text_disp_object[100], char text_disp_num[100], int action, int start) {
    G_rgb(0, 0, 0) ;
    G_clear() ;

    for (int i = 0; i < numobjects; i++) {
        draw_single_object(i) ;
    }

    G_rgb(0.1, 0.6, 0.3) ;
    G_fill_rectangle(10, window_size-70, window_size/5, 60) ;

    G_rgb(1,1,1) ;
    G_draw_string(text_disp_object, 15, window_size-24) ;
    G_draw_string(text_disp_num, 15, window_size-44) ;
    if (action == 'r') {
        G_draw_string("action : rotate", 15, window_size-64) ;
    } else if (action == 't') {
        G_draw_string("action : translate", 15, window_size-64) ;
    }
    
    if (start) {
      G_rgb(0.3, 0.6, 0.6) ;
      G_fill_rectangle(10, window_size-280, window_size/5, 200) ;
      
      G_rgb(1,1,1) ;
      G_draw_string("to use program : ", 15, window_size-94) ;
      G_draw_string("[0-9] > change objects", 15, window_size-114) ;
      G_draw_string("t > translate", 15, window_size-134) ;
      G_draw_string("r > rotate", 15, window_size-154) ;
      G_draw_string("v > reverse backface elim", 15, window_size-174) ;
      G_draw_string("x, y, z > perform action", 15, window_size-194) ;
      G_draw_string("c > reverse direction", 15, window_size-214) ;
      G_draw_string("+, - > rotation speed", 15, window_size-234) ;
      G_draw_string("m > see this menu", 15, window_size-254) ;
      G_draw_string("q > quit program", 15, window_size-274) ; 
    } else {
      G_rgb(0.3, 0.6, 0.6) ;
      G_fill_rectangle(10, window_size-100, window_size/5, 20) ;
      G_rgb(1,1,1) ;
      G_draw_string("press m to see menu", 15, window_size-94) ;
    }
    
    return 1 ; 
}




int main(int numFiles, char **file) 
{
    FILE *fin ; 
    int sign = 1 ; 
    int action = 't' ; 
    int onum = 0 ; 
    int q, k ; 
    double V[4][4], M1[4][4], M2[4][4] ; 
    double XT, YT, ZT, COS, SIN ;
    double radians ;  
    int rotation_speed = 8 ;
    int initial_depth = 10 ; 

    char text_disp_object[100] ;
    char text_disp_num[100] ;
    char text_disp_action[100] ; 

    numobjects = (numFiles - 2 > MAXOBJS) ? MAXOBJS : numFiles - 1 ;

    for (int i = 1; i < numFiles; i++) {
        fin = fopen(file[i], "r") ; 
        if (fin == NULL) {
            printf("Error: invalid file '%s'\n", file[i]) ;
            exit(0) ; 
        }
        read_object(fin, i - 1) ; 

        cull_flag[i - 1] = 1 ;
    }

    G_init_graphics(window_size, window_size) ;

    int start = 1 ; 

    while (1) {

        M3d_make_identity(V) ;

        if (start) {
            q = 't' ;
            for (int i = 0 ; i < numobjects ; i++) {
                for (int j = 0 ; j < initial_depth ; j++) {
                    M3d_make_translation(V, 0, 0, sign * 1) ;
                    M3d_mat_mult_points(x[i], y[i], z[i], V, x[i], y[i], z[i], numpoints[i] + 1) ;
                }
            }
        } else q = G_wait_key() ; 

        radians = (M_PI / 180) * rotation_speed ;

        if (q == 'q') exit(0) ; 
        else if (q == 'c') sign = -sign ;
        else if (q == 't' || q == 'r') action = q ;

        else if (q == '=') rotation_speed++ ;
        else if (q == '-') rotation_speed = (rotation_speed > 1) ? rotation_speed - 1 : 1 ;

        else if ('0' <= q && q <= '9') {
            k = q - '0' ;  
            if (k < numobjects) { onum = k ; } 
        }

        // put file name and number into char strings
        snprintf(text_disp_object, sizeof(text_disp_object), "shape : %s", 4 + file[1 + onum]) ;
        snprintf(text_disp_num, sizeof(text_disp_num), "location : %i", onum) ;

        XT = x[onum][numpoints[onum]] ;
        YT = y[onum][numpoints[onum]] ;
        ZT = z[onum][numpoints[onum]] ;
        COS = cos(sign * radians) ;
        SIN = sin(sign * radians) ;

        if (q == 'p') printf("rotation speed = %d\n", rotation_speed) ;
        else if (q == 'v') cull_flag[onum] = -cull_flag[onum] ;
        else if (q == 'm') start = 1 ; 

        if ((q == 'x' || q == 'y' || q == 'z') && action == 't') {
            if (q == 'x') M3d_make_translation(V, sign * 1, 0, 0) ;
            if (q == 'y') M3d_make_translation(V, 0, sign * 1, 0) ;
            if (q == 'z') M3d_make_translation(V, 0, 0, sign * 1) ;
        }
        
        else if ((q == 'x' || q == 'y' || q == 'z') && action == 'r') {
            M3d_make_translation(V, -XT, -YT, -ZT) ;
            if (q == 'x') M3d_make_x_rotation_cs(M1, COS, SIN) ;
            if (q == 'y') M3d_make_y_rotation_cs(M1, COS, SIN) ;
            if (q == 'z') M3d_make_z_rotation_cs(M1, COS, SIN) ;
            M3d_mat_mult(M1, M1, V) ;
            M3d_make_translation(M2, XT, YT, ZT) ;
            M3d_mat_mult(V, M2, M1) ;
        }

        M3d_mat_mult_points(x[onum], y[onum], z[onum], V, x[onum], y[onum], z[onum], numpoints[onum] + 1) ;

        display(text_disp_object, text_disp_num, action, start) ; 
        start = 0 ;
    }
}
