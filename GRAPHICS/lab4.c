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

double red = 1.0 ; 
double grn = 0.2 ; 
double blu = 0.2 ; 

int halfangle = 45 ; 
double window_size = 800 ;

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
}



int draw_3d_polygon(double xp[], double yp[], double zp[], int np) {

    double xbb[100], ybb[100] ;
    double H = tan(halfangle * (M_PI / 180)) ;
    double h_wind = window_size / 2 ; 
    int i ;

    // Project the 3D points into 2D space
    for (i = 0; i < np; i++) {
        xbb[i] = ((h_wind / H) * (xp[i] / zp[i]) + h_wind) ;
        ybb[i] = ((h_wind / H) * (yp[i] / zp[i]) + h_wind) ; 
    }

    G_rgb(red, grn, blu) ; 
    G_polygon(xbb, ybb, np) ; 

    return 1 ; 

}



int draw_single_object(int object) {
    int h, i, j ; 
    double xp[100], yp[100], zp[100] ;
    int np;

    for (i = 0; i < numpolys[object]; i++) {

        np = psize[object][i] ; 

        for (j = 0; j < np; j++) {
            h = con[object][i][j] ;
            xp[j] = x[object][h] ; 
            yp[j] = y[object][h] ; 
            zp[j] = z[object][h] ; 

        }

        draw_3d_polygon(xp, yp, zp, np) ;

    }

    return 1 ; 
}



int translate(int onum, double dx, double dy, double dz) {

    int i ;
    double m[4][4] ; 
    M3d_make_translation(m, dx, dy, dz) ;

    return 1 ; 

}



int main(int numFiles, char **file) 
{
    FILE *fin ; 
    char fname[100] ; 
    int sign = 1 ; 
    int action = 't' ; 
    int onum = 0 ; 
    int q, k ; 
    double V[4][4], M1[4][4], M2[4][4] ; 
    double XT, YT, ZT, COS, SIN ;
    double radians ; 
    int rotation_speed = 10 ;

    

    if (numFiles - 2 > 2) {
        numobjects = numFiles - 2 ; 
    } else {
        numobjects = numFiles ; 
    }

    if (numobjects > MAXOBJS) {
        exit(0) ;
    }

    for (int i = 1; i < numFiles; i++) {
        fin = fopen(file[i], "r") ; 
        if (fin == NULL) {
            printf("error: invalid file '%s\n'", file[i]) ;
            exit(0) ; 
        }
        read_object(fin, i - 1) ; 
    }


    /* ######################################################### */
    /* ######################## GRAPHCS ######################## */
    /* ######################################################### */


    G_init_graphics(window_size, window_size) ;
    G_rgb(0, 0, 0) ;  
    G_clear() ;
    draw_single_object(0) ; 

    while (1) {

        M3d_make_identity(V) ;
        q = G_wait_key() ; 
        radians = (M_PI / 180) * rotation_speed ;

        // initiation keys
        if (q == 'q') exit(0) ; 
        else if (q == 'c') sign = -sign ;
        else if (q == 't') action = q ; 
        else if (q == 'r') action = q ;

        else if (q == '=') {

            rotation_speed ++ ;

        } else if (q == '-') {

            rotation_speed -- ; 
            if (rotation_speed <= 0) rotation_speed = 1 ; 

        } else if (('0' <= q) && (q <= '9')) {

            k = q - '0' ;  
            if (k < numobjects) { onum = k ; }

        }

        XT = x[onum][numpoints[onum]] ; 
        YT = y[onum][numpoints[onum]] ; 
        ZT = z[onum][numpoints[onum]] ;
        COS = cos(sign * radians) ;
        SIN = sin(sign * radians) ;

        
        if (q == 'p') printf("rotation speed = %d\n", rotation_speed) ;

        // translate
        else if ((q == 'x') && (action == 't')) M3d_make_translation(V, sign * 1, 0, 0) ;
        else if ((q == 'y') && (action == 't')) M3d_make_translation(V, 0, sign * 1, 0) ;
        else if ((q == 'z') && (action == 't')) M3d_make_translation(V, 0, 0, sign * 1) ;
        
        // rotate
        else if ((q == 'x') && (action == 'r')) {

            M3d_make_translation(V, -XT, -YT, -ZT) ;
            M3d_make_x_rotation_cs(M1, COS, SIN) ;
            M3d_mat_mult(M1, M1, V) ;
            M3d_make_translation(M2, XT, YT, ZT) ;
            M3d_mat_mult(V, M2, M1) ;

        } else if ((q == 'y') && (action == 'r')) {

            M3d_make_translation(V, -XT, -YT, -ZT) ;
            M3d_make_y_rotation_cs(M1, COS, SIN) ;
            M3d_mat_mult(M1, M1, V) ;
            M3d_make_translation(M2, XT, YT, ZT) ;
            M3d_mat_mult(V, M2, M1) ;

        } else if ((q == 'z') && (action == 'r')) {

            M3d_make_translation(V, -XT, -YT, -ZT) ;
            M3d_make_z_rotation_cs(M1, COS, SIN) ;
            M3d_mat_mult(M1, M1, V) ;
            M3d_make_translation(M2, XT, YT, ZT) ;
            M3d_mat_mult(V, M2, M1) ;

        } 

        M3d_mat_mult_points (x[onum],y[onum],z[onum],  V,
           x[onum],y[onum],z[onum],numpoints[onum]+1) ;
          // the numpoints[onum]+1 is because we have stored the center
          // of the object at the arrays' end

        G_rgb(0, 0, 0) ;
        G_clear() ;
        draw_single_object(onum) ;

        G_display_image() ; 
    }

}