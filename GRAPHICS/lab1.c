#include "FPToolkit.c"

void find_min_max_y(double y[], int n, double *min_y, double *max_y) {
    *min_y = y[0];
    *max_y = y[0];

    for (int i = 1; i < n; i++) {
        if (y[i] < *min_y) {
            *min_y = y[i];
        }
        if (y[i] > *max_y) {
            *max_y = y[i];
        }
    }
}

void my_polygon(double x[], double y[], int n) {
    int i;
    G_rgb(1, 1, 1);

    for (i = 0; i < n - 1; i++) {
        G_line(x[i], y[i], x[i + 1], y[i + 1]);
    }
    G_line(x[0], y[0], x[n - 1], y[n - 1]);
}

int click_and_save(double x[], double y[]) {
    int n = 0;
    double p[2];

    while (1) {
        G_wait_click(p);
        if (p[0] < 40 && p[1] < 40) break;
        x[n] = p[0];
        y[n] = p[1];
        n++;
       
        G_rgb(1, 1, 1);
        G_fill_circle(p[0], p[1], 2);
    }

    return n;
}

int find_intersections(double a[], double b[], int n, double yy, double array[]) {
    int i;
    double miny, maxy;
    int count = 0;

    for (i = 0; i < n - 1; i++) {
        if (b[i] < b[i + 1]) {
            miny = b[i];
            maxy = b[i + 1];
        } else {
            miny = b[i + 1];
            maxy = b[i];
        }
       
        if (miny < yy && yy < maxy) {
            array[count] = a[i] + (yy - b[i]) * (a[i + 1] - a[i]) / (b[i + 1] - b[i]);
            G_rgb(1, 1, 1);
            G_point(array[count], yy);
            count++;
        }
    }

    if (b[n - 1] < b[0]) {
        miny = b[n - 1];
        maxy = b[0];
    } else {
        miny = b[0];
        maxy = b[n - 1];
    }

    if (miny < yy && yy < maxy) {
        array[count] = a[n - 1] + (yy - b[n - 1]) * (a[0] - a[n - 1]) / (b[0] - b[n - 1]);
        G_rgb(1, 1, 1);
        G_point(array[count], yy);
        count++;
    }
   
    return count;
}

void selection_sort(double x[], int n) {
    int i, s, j;
    double tmp;

    for (i = 0; i < n; i++) {
        s = i;
        for (j = i + 1; j < n; j++) {
            if (x[j] < x[s]) {
                s = j;
            }
        }
        tmp = x[i];
        x[i] = x[s];
        x[s] = tmp;
    }
}

void print_array(double x[], int n) {
    printf("n = %d\n", n);
    for (int i = 0; i < n; i++) {
        printf("x[%d] = %lf\n", i, x[i]);
    }
    printf("\n");
}

void drawlines(double array[], double yy, int count) {
    selection_sort(array, count);

    for (int i = 0; i < count - 1; i += 2) {
        G_rgb(1, 0, 0);
        //printf("[%lf, %lf]\n", array[i], array[i + 1]);
        G_line(array[i], yy, array[i + 1], yy);
    }
}

void homescreen() {
    G_rgb(0, 0, 0);
    G_clear();
    G_rgb(1, 0, 0);
    G_fill_rectangle(0, 0, 40, 40);
}

void clear_array(double array[], int size) {
    for (int i = 0 ; i < size ; i++) {
      array[i] = 0 ;
    }
}

int main()
{
    G_init_graphics(800, 800) ;
    homescreen() ;

    double a[1000], b[1000], array[1000] ;
    int nab, count ;
    double yy, LOWERY, UPPERY ;

    nab = click_and_save(a, b) ;

    find_min_max_y(b, nab, &LOWERY, &UPPERY) ;

    for (yy = LOWERY + 0.1; yy < UPPERY; yy++) {
        clear_array(array, 1000) ;
        count = find_intersections(a, b, nab, yy, array) ;
        drawlines(array, yy, count) ;
        G_wait_key() ;
    }

    my_polygon(a, b, nab) ;

    while(1) {
      G_wait_key() ;
    }
}

