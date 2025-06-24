#include "FPToolkit.c"
#include "M3d_matrix_tools.c"

#define MAXPTS 59000
#define MAXPOLYS 57500
#define MAXOBJS 10

int numobjects;
int numpoints[MAXOBJS];
int numpolys[MAXOBJS];
double x[MAXOBJS][MAXPTS];
double y[MAXOBJS][MAXPTS];
double z[MAXOBJS][MAXPTS];
int psize[MAXOBJS][MAXPOLYS];
int con[MAXOBJS][MAXPOLYS][20];

int halfangle = 45 ; 
double window_size = 800 ;

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
    int objnum;
    int polynum;
    double dist;
} POLY;

POLY Polys[MAXPOLYS];
int PolyCount;

// Function to read a single object from file
int readObject(FILE *file, int objectIndex) {
    int i, j;

    // Read points
    fscanf(file, "%d", &numpoints[objectIndex]);
    if (numpoints[objectIndex] >= MAXPTS) {
        printf("Error: Exceeded maximum points limit of %d.\n", MAXPTS);
        exit(1);
    }

    for (i = 0; i < numpoints[objectIndex]; i++) {
        fscanf(file, "%lf %lf %lf", &x[objectIndex][i], &y[objectIndex][i], &z[objectIndex][i]);
    }

    // Read polygons
    fscanf(file, "%d", &numpolys[objectIndex]);
    if (numpolys[objectIndex] > MAXPOLYS) {
        printf("Error: Exceeded maximum polygons limit of %d.\n", MAXPOLYS);
        exit(1);
    }

    for (i = 0; i < numpolys[objectIndex]; i++) {
        fscanf(file, "%d", &psize[objectIndex][i]);
        for (j = 0; j < psize[objectIndex][i]; j++) {
            fscanf(file, "%d", &con[objectIndex][i][j]);
        }
    }

    return 0;
}


void init_poly_struct() {
    PolyCount = 0;
    for (int obj = 0; obj < numobjects; obj++) {
        for (int poly = 0; poly < numpolys[obj]; poly++) {
            // Compute the average depth
            double avgZ = 0.0;
            for (int j = 0; j < psize[obj][poly]; j++) {
                int pointIndex = con[obj][poly][j];
                avgZ += z[obj][pointIndex];
            }
            avgZ /= psize[obj][poly];
            Polys[PolyCount++] = (POLY){obj, poly, avgZ};
        }
    }
}



// Function to draw a single 3D polygon
void draw3DPolygon(double xp[], double yp[], double zp[], int numpoints, int object) {
    double screenX[100], screenY[100];
    double halfangle = 60;
    double scale = 400.0 / tan(halfangle * M_PI / 180.0);
    int i;

    for (i = 0; i < numpoints; i++) {
        screenX[i] = scale * (xp[i] / zp[i]) + 400;
        screenY[i] = scale * (yp[i] / zp[i]) + 400;
    }

    G_rgb(colors[object][0], colors[object][1], colors[object][2]);
    G_fill_polygon(screenX, screenY, numpoints);

    G_rgb(0, 0, 0);
    G_polygon(screenX, screenY, numpoints);
}


// Function to draw a single object with back-face culling
void drawObject(int objectIndex, int cullBackFaces) {
    int i, j, h;
    double xp[1000], yp[1000], zp[1000];
    int np;

    for (i = 0; i < numpolys[objectIndex]; i++) {

        np = psize[objectIndex][i];

        for (j = 0; j < np; j++) {
            h = con[objectIndex][i][j];
            xp[j] = x[objectIndex][h];
            yp[j] = y[objectIndex][h];
            zp[j] = z[objectIndex][h];
        }

        double vector1[3] = { xp[1] - xp[0], yp[1] - yp[0], zp[1] - zp[0] };
        double vector2[3] = { xp[2] - xp[0], yp[2] - yp[0], zp[2] - zp[0] };
        double normal[3];
        M3d_x_product(normal, vector1, vector2);

        //double dotProduct = normal[0] * xp[0] + normal[1] * yp[0] + normal[2] * zp[0];
        //int visible = (dotProduct > 0) ^ cullBackFaces;

        //if (visible) {
            draw3DPolygon(xp, yp, zp, np, objectIndex);
        //}
    }
}


int compare (const void *p, const void *q)
{
  POLY *a, *b ;

  a = (POLY*)p ;
  b = (POLY*)q ;

  if  (((*a).dist) < ((*b).dist)) return -1 ;
  else if (((*a).dist) > ((*b).dist)) return 1 ;
  else return 0 ;
}


// Function to draw all objects using the painter's algorithm
void drawAllObjects() {
    int object, i, j, np, h;
    double xp[1000], yp[1000], zp[1000];
    PolyCount = 0;

    // Compute distances and populate Poly array
    for (object = 0; object < numobjects; object++) {
        for (i = 0; i < numpolys[object]; i++) {
            np = psize[object][i];
            for (j = 0; j < np; j++) {
                h = con[object][i][j];
                xp[j] = x[object][h];
                yp[j] = y[object][h];
                zp[j] = z[object][h];
            }
            Polys[PolyCount++] = (POLY){object, i, zp[0]};
        }
    }

    // Sort polygons by distance
    qsort(Polys, PolyCount, sizeof(POLY), compare);

    // Draw polygons in sorted order
    for (int k = PolyCount - 1; k >= 0; k--) {
        object = Polys[k].objnum;
        i = Polys[k].polynum;
        np = psize[object][i];
        for (j = 0; j < np; j++) {
            h = con[object][i][j];
            xp[j] = x[object][h];
            yp[j] = y[object][h];
            zp[j] = z[object][h];
        }
        draw3DPolygon(xp, yp, zp, np, object);
    }
}


int display(char text_disp_object[100], char text_disp_num[100], int action, int start) {
    G_rgb(0, 0, 0) ;
    G_clear() ;
    drawAllObjects();

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
    FILE *fin;
    int sign = 1;
    int action = 't';
    int onum = 0;
    int q, k;
    double V[4][4], M1[4][4], M2[4][4];
    double XT, YT, ZT, COS, SIN;
    double radians;
    int rotation_speed = 14;
    int initial_depth = 8;

    char text_disp_object[100];
    char text_disp_num[100];

    numobjects = (numFiles - 2 > MAXOBJS) ? MAXOBJS : numFiles - 1;

    for (int i = 1; i < numFiles; i++) {
        fin = fopen(file[i], "r");
        if (fin == NULL) {
            printf("Error: invalid file '%s'\n", file[i]);
            exit(0);
        }
        readObject(fin, i - 1);
    }

    G_init_graphics(window_size, window_size);

    int start = -1;

    while (1) {

        init_poly_struct(); 

        M3d_make_identity(V);

        if (start == -1) {
            q = 't';
            for (int i = 0; i < numobjects; i++) {
                for (int j = 0; j < initial_depth; j++) {
                    M3d_make_translation(V, 0, 0, sign * 1);
                    M3d_mat_mult_points(x[i], y[i], z[i], V, x[i], y[i], z[i], numpoints[i] + 1);
                }
            }
            start = 0;
        } else {
            q = G_wait_key();
        }

        radians = (M_PI / 180) * rotation_speed;

        if (q == 'q') exit(0);
        else if (q == 'c') sign = -sign;
        else if (q == 't' || q == 'r') action = q;

        else if (q == '=') rotation_speed++;
        else if (q == '-') rotation_speed = (rotation_speed > 1) ? rotation_speed - 1 : 1;

        else if ('0' <= q && q <= '9') {
            k = q - '0';
            if (k < numobjects) { onum = k; }
        }

        // Update display information
        snprintf(text_disp_object, sizeof(text_disp_object), "shape : %s", 4 + file[1 + onum]);
        snprintf(text_disp_num, sizeof(text_disp_num), "location : %i", onum);

        XT = x[onum][numpoints[onum]];
        YT = y[onum][numpoints[onum]];
        ZT = z[onum][numpoints[onum]];
        COS = cos(sign * radians);
        SIN = sin(sign * radians);

        if (q == 'p') printf("rotation speed = %d\n", rotation_speed);
        else if (q == 'm') start = (start > 0) ? 0 : 1;

        if ((q == 'x' || q == 'y' || q == 'z') && action == 't') {
            if (q == 'x') M3d_make_translation(V, sign * 1, 0, 0);
            if (q == 'y') M3d_make_translation(V, 0, sign * 1, 0);
            if (q == 'z') M3d_make_translation(V, 0, 0, sign * 1);
        }
        
        else if ((q == 'x' || q == 'y' || q == 'z') && action == 'r') {
            M3d_make_translation(V, -XT, -YT, -ZT);
            if (q == 'x') M3d_make_x_rotation_cs(M1, COS, SIN);
            if (q == 'y') M3d_make_y_rotation_cs(M1, COS, SIN);
            if (q == 'z') M3d_make_z_rotation_cs(M1, COS, SIN);
            M3d_mat_mult(M1, M1, V);
            M3d_make_translation(M2, XT, YT, ZT);
            M3d_mat_mult(V, M2, M1);
        }

        M3d_mat_mult_points(x[onum], y[onum], z[onum], V, x[onum], y[onum], z[onum], numpoints[onum] + 1);

        display(text_disp_object, text_disp_num, action, start) ; 
    }
}
