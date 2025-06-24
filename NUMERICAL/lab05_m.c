#include "FPToolkit.c"

#define MAXN 1000
#define MAXDIST 7
#define SIZE 800
#define EXIT 30
#define STEP 0.01


int click_and_save(double x[], double y[]) {
    int n = 0;
    double p[2];

    while (1) {
		G_wait_click(p);
		if (p[1] < EXIT) break;

		x[n] = p[0];
		y[n] = p[1];

		G_rgb(1,1,1);
		G_circle(p[0], p[1], 4);
        n++;
	}
	return n;
}


int tridiagonal(double L[], double D[],
                double R[], double Q[],
                int n, double x[]) {
    int i;
    double m;

    for (i = 1; i < n; i++) {
        if (D[i - 1] == 0) return 0;

        m = L[i] / D[i - 1];

        D[i] -= m * R[i - 1];
        Q[i] -= m * Q[i - 1];
    }

    if (D[n - 1] == 0) return 0;

    x[n - 1] = Q[n - 1] / D[n - 1];

    for (i = n - 2; i >= 0; i--) {
        if (D[i] == 0) return 0;
        x[i] = (Q[i] - R[i] * x[i + 1]) / D[i];
    }
    return 1;
}


int natural_spline(double x[], double y[], int n,
                   double M[], double A[], double B[]) {
    double L[MAXN], D[MAXN], R[MAXN], Q[MAXN], AB[MAXN];
    int i,j;
    double a,b,c,d,e,f,g;
    int s;

    for (i = 0; i <= n-2; i++) {
        M[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]);
    }

    j = 0;
    D[j] = 1; R[j] = -(x[1] - x[0]); Q[j] = 0; j++;

    for (i = 1; i <= n-2; i++) {
        a = x[i] - x[i-1];
        b = a*a;
        c = x[i+1] - x[i];
        d = b;
        e = x[i-1] - x[i+1];
        f = a*c;
        g = M[i] - M[i-1];

        L[j] = a; D[j] = b; R[j] = c; Q[j] =  g; j++;
        L[j] = d; D[j] = e; R[j] = f; Q[j] = -g; j++;
    }

    L[j] = 1; D[j] = 2*(x[i] - x[i-1]); Q[j] = 0; j++;

    s = tridiagonal(L,D,R,Q,j,AB);

    j = 0;
    for (i = 0; i <= n-2; i++) {
        A[i] = AB[j++];
        B[i] = AB[j++];
    }
    return s;
}


void graph_natural_spline(double x[], double y[], int n) {
    double A[MAXN], B[MAXN], M[MAXN];
    double xp,yp;
    int i,s;

    s = natural_spline(x, y, n, M, A, B);
    if (s == 0) return;

    G_rgb(0.5,0.5,1.0);
    for (i = 0; i <= n-2; i++) {
        if (x[i] < x[i+1]) {
            for (xp = x[i]; xp < x[i+1]; xp += STEP) {

                yp = y[i] + M[i]*(xp - x[i]) +
                    (xp - x[i])*(xp - x[i+1])*
                    (A[i] + B[i]*(xp - x[i]));
                G_point(xp, yp);
            }
        } else {
            for (xp = x[i]; xp >= x[i+1]; xp -= STEP) {
                yp = y[i] + M[i]*(xp - x[i]) +
                    (xp - x[i])*(xp - x[i+1])*
                    (A[i] + B[i]*(xp - x[i]));
                G_point(xp, yp);
            }
        }
    }
}


void clear() {
    double a,b,c;

    a = SIZE-EXIT;
    b = SIZE;
    c = EXIT;

	G_rgb(0,0,0);
	G_clear();
	G_rgb(1,0.3,0.3);
	G_fill_rectangle(0,0,b,c);

    G_rgb(1,0,0);
    G_fill_rectangle(a,a,c,c);
    G_rgb(1,1,1);
    G_line(a,a,b,b);
    G_line(a,b,b,a);
}


int bubblesort(double x[], double y[], int n, int I) {
    double tempx, tempy;

    // Bubble right
    while (I < n - 1 && x[I] > x[I+1]) {
        tempx = x[I];       tempy = y[I];
        x[I] = x[I+1];      y[I] = y[I+1];
        x[I+1] = tempx;     y[I+1] = tempy;
        I++;
    }

    // Bubble left
    while (I > 0 && x[I] < x[I-1]) {
        tempx = x[I];       tempy = y[I];
        x[I] = x[I-1];      y[I] = y[I-1];
        x[I-1] = tempx;     y[I-1] = tempy;
        I--;
    }

    return I;
}

void redraw_pts(double x[], double y[], int n) {
    clear();
    G_rgb(1,1,1);
    for (int i = 0; i < n; i++) {
        G_circle(x[i], y[i], 4);
    }
}


int is_on_point(int p[], double x[], double y[], int n) {
    for (int i = 0; i < n; i++) {
        double dx = x[i] - p[0];
        double dy = y[i] - p[1];
        double dist2 = dx * dx + dy * dy;
        if (dist2 < MAXDIST * MAXDIST) return i;
    }
    return -1;
}


void move_point(double x[], double y[], int i, int n) {
    int p[2];

    while (S_mouse_coord_window(p) == 1) {
        p[1] = SIZE - p[1]; // y coord is inversed

        x[i] = p[0];
        y[i] = p[1];

        i = bubblesort(x,y,n,i);
        redraw_pts(x,y,n);
        graph_natural_spline(x,y,n);
        G_display_image(); // needed for animation
    }
}


int main() {
    int key, n;
    double x[MAXN], y[MAXN];
    int p[2];

	G_init_graphics(SIZE, SIZE);

    clear();

    n = click_and_save(x, y);

    // still haven't found a good way to close the program

    while (1) {
        int mouse_state = S_mouse_coord_window(p);
        p[1] = SIZE - p[1]; // y coord is inversed

        redraw_pts(x,y,n);
        graph_natural_spline(x,y,n);

        if (mouse_state == 1) {
            if (p[0] > SIZE-EXIT && p[1] > SIZE-EXIT) exit(0);
            int I = is_on_point(p, x,y,n);
            if (I >= 0) move_point(x,y,I,n);
        }

        G_display_image();
    }
}
