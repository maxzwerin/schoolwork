#include "FPToolkit.c"

/* LAB06 :: B SPLINES

   points 0123 make a small cubic between 1 & 2. 
   0 and 3 are control points for C0123

   C0123 = a(t)(x0y0) + b(t)(x1y1) + c(t)(x2y2) + d(t)(x3y3)

   equations:
   C 0123(1) = C 1234(0)
   C'0123(1) = C'1234(0)
   C"0123(1) = C"1234(0)
 
   a(t) = a0 + a1t + a2t^2 + a3t^3
   b(t) = b0 + b1t + b2t^2 + b3t^3
   c(t) = c0 + c1t + c2t^2 + c3t^3
   d(t) = d0 + d1t + d2t^2 + d3t^3
   -------------------------------
   a'(t) = a1 + 2a2t + 3a3t^2
   b'(t) = b1 + 2b2t + 3b3t^2
   c'(t) = c1 + 2c2t + 3c3t^2
   d'(t) = d1 + 2d2t + 3d3t^2
   -------------------------------
   a"(t) = 2a2 + 6a3t
   b"(t) = 2b2 + 6b3t
   c"(t) = 2c2 + 6c3t
   d"(t) = 2d2 + 6d3t

   --- 16 unknown #s

   a(1)P0 + b(1)P1 + c(1)P2 + d(1)P3 = a(0)P1 + b(0)P2 + c(0)P3 + d(0)P4
  *a(1) = 0     | a'(1) = 0      | a"(1) = 0
  *b(1) = a(0)  | b'(1) = a'(0)  | b"(1) = a"(0)
  *c(1) = b(0)  | c'(1) = b'(0)  | c"(1) = b"(0)
  *d(1) = c(0)  | d'(1) = c'(0)  | d"(1) = c"(0)
  *   0 = d(0)  |     0 = d'(0)  |     0 = d"(0)
 
   a0+a1+a2+a3 = 0  | a1 + 2a2 + 3a3 = 0  | 2a2 + 6a3 = 0
   b0+b1+b2+b3 = a0 | b1 + 2b2 + 3b3 = a1 | 2b2 + 6b3 = 2a2
   c0+c1+c2+c3 = b0 | c1 + 2c2 + 3c3 = b1 | 2c2 + 6c3 = 2b2
   d0+d1+d2+d3 = c0 | d1 + 2d2 + 3d3 = c1 | 2d2 + 6d3 = 2c2
             0 = d0 |              0 = d1 |         0 = 2d2

   a0 + b0 + c0 + d0 = 1

    a0 = 1/6    b0 = 2/3    c0 = 1/6    d0 = 0
    a1 = -1/2   b1 = 0      c1 = 1/2    d1 = 0
    a2 = 1/2    b2 = -1     c2 = 1/2    d2 = 0
    a3 = -1/6   b3 = 1/2    c3 = -1/2   d3 = 1/6
*/

#define MAXN 1000
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

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

int tridiagonal (double L[], double D[],
                 double R[], double Q[],
                 int n, double x[]) {
    int i;
    double m;

    // forward elimination
    for (i = 1; i < n; i++) {
        if (D[i - 1] == 0) return 0;

        m = L[i] / D[i - 1];

        D[i] -= m * R[i - 1];
        Q[i] -= m * Q[i - 1];
    }

    if (D[n - 1] == 0) return 0;

    // back substitution
    x[n - 1] = Q[n - 1] / D[n - 1];

    for (i = n - 2; i >= 0; i--) {
        if (D[i] == 0) return 0;
        x[i] = (Q[i] - R[i] * x[i + 1]) / D[i];
    }

    return 1;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void graph_b_spline(double x[], double y[], int n) {
    if (n < 4) return;

    double a[4],b[4],c[4],d[4],M[4];
    double t,t2,t3;
    int i,j,k;
    double X,Y;

    a[0] = 1.0/6;    b[0] = 2.0/3;    c[0] = 1.0/6;    d[0] = 0.0;
    a[1] = -0.5;     b[1] = 0.0;      c[1] = 0.5;      d[1] = 0.0;
    a[2] = 0.5;      b[2] = -1.0;     c[2] = 0.5;      d[2] = 0.0;
    a[3] = -1.0/6;   b[3] = 0.5;      c[3] = -0.5;     d[3] = 1.0/6;

    G_rgb(1,1,0.5);
    for (i = 0; i < n - 3; i++) {
        for (t = 0.0; t <= 1.0; t += STEP) {
            t2 = t * t;
            t3 = t2 * t;

            M[0] = a[0] + a[1]*t + a[2]*t2 + a[3]*t3;
            M[1] = b[0] + b[1]*t + b[2]*t2 + b[3]*t3;
            M[2] = c[0] + c[1]*t + c[2]*t2 + c[3]*t3;
            M[3] = d[0] + d[1]*t + d[2]*t2 + d[3]*t3;

            X = Y = 0;
            for (j = 0; j < 4; j++) {
                X += M[j]*x[i + j];
                Y += M[j]*y[i + j];
            }
            G_point(X, Y);
        }
    }
    G_draw_string("B SPLINE",200,SIZE-80);
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

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

    s = tridiagonal (L,D,R,Q,j,AB);

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

    s = natural_spline(x,y,n,M,A,B);
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
                    (xp - x[i])*(xp - x[i+1])*(A[i] + B[i]*(xp - x[i]));
                G_point(xp, yp);
            }
        }
    }
    G_draw_string("N SPLINE",200,SIZE-60);
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void best_fit_line (double x[], double y[], int n, 
                    double *A, double *B) {
    double a, b;
    double sumX = 0, sumY = 0, sumXX = 0, sumXY = 0;
    int i;

    for (i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumXX += x[i] * x[i];
        sumXY += x[i] * y[i];
    }

    // n     +  sumX    =  sumY
    // sumX  +  sumYY   =  sumXY

    double denom = (n * sumXX - sumX * sumX);

    *A = (sumY * sumXX - sumX * sumXY) / denom;
    *B = (n * sumXY - sumX * sumY)     / denom;
}


void graph_best_fit_line (double x[], double y[], int n) {
    double p, A, B;
    double x0, y0, x1, y1;

    G_rgb(0.5,1.0,1.0);

    best_fit_line (x, y, n, &A, &B);

    x0 = 0;
    y0 = A;
    x1 = SIZE;
    y1 = A + (B * SIZE);

    G_line(x0, y0, x1, y1);
    G_draw_string("BEST FIT",200,SIZE-40);
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

double lagrange (double X, double x[], double y[], int n) {
    double E, D, Y = 0, sum;
    int i, j;

    for (i = 0; i < n; i++) {
        E = 1; D = 1;
        for (j = 0; j < n; j++) {
            if (j != i) {
                E *= (X - x[j]);
                D *= (x[i] - x[j]);
            }
        }
        sum = y[i] * E / D;
        Y += sum;
    }
    return Y;
}

void graph_lagrange (double x[], double y[], int n) {
    int i;
    double xmin, xmax, X, Y;
    double ymin = EXIT, ymax = SIZE;
    xmin = xmax = x[0] ;

    for (i = 1 ; i < n ; i++) {
        if (x[i] < xmin) { xmin = x[i] ; }
        if (x[i] > xmax) { xmax = x[i] ; }
    }

	G_rgb(0.5,1.0,0.5);
    for (X = xmin; X <= xmax; X += STEP) {
        Y = lagrange (X, x,y,n);
        if (Y < ymin || Y > ymax) continue;
        G_point(X,Y);
    }
    G_draw_string("LAGRANGE",200,SIZE-20);
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void menu() {
    int x = 20, y = SIZE-x, dy = 20;

    G_rgb(0.7,0.3,0.7);
    G_rectangle(195,SIZE-170,65,165);
    G_rgb(0.3,0.3,0.3);
    G_draw_string("LAGRANGE",200,SIZE-20);
    G_draw_string("BEST FIT",200,SIZE-40);
    G_draw_string("N SPLINE",200,SIZE-60);
    G_draw_string("B SPLINE",200,SIZE-80);


	G_rgb(0.7,0.3,0.7);
    G_fill_rectangle(5,SIZE-170,180,165);
    G_rgb(1,1,1);
    G_draw_string("1 to display lagrange",      x,y); y-=dy;
    G_draw_string("2 to dsiplay best fit line", x,y); y-=dy;
    G_draw_string("3 to display natural spline",x,y); y-=dy;
    G_draw_string("4 to display b spline",      x,y); y-=dy;
    G_draw_string("A to display all",           x,y); y-=dy;
    G_draw_string("C to display none",          x,y); y-=dy;
    G_draw_string("R to use new points",        x,y); y-=dy;
    G_draw_string("Q to quit",                  x,y);
}


void clear() {
	G_rgb(0,0,0);
	G_clear();
	G_rgb(1,0.3,0.3);
	G_fill_rectangle(0,0,SIZE,EXIT);
}


void redraw_pts (double x[], double y[], int n) {
    clear();
    G_rgb(1,1,1);
    for (int i = 0; i < n; i++) {
        G_circle(x[i], y[i], 4);
    }
}


void graph_all(double x[], double y[], int n) {
    graph_lagrange(x,y,n);
    graph_best_fit_line(x,y,n);
    graph_natural_spline(x,y,n);
    graph_b_spline(x,y,n);
}


int main() {
    int key, n;
    double x[1000], y[1000];

	G_init_graphics(SIZE, SIZE);

	while (1) {
        clear();
        n = click_and_save(x, y);
        menu();
        do {
            key = G_wait_key();
            if (key == '1')      { redraw_pts(x,y,n); menu(); graph_lagrange(x,y,n); }
            else if (key == '2') { redraw_pts(x,y,n); menu(); graph_best_fit_line(x,y,n); }
            else if (key == '3') { redraw_pts(x,y,n); menu(); graph_natural_spline(x,y,n); }
            else if (key == '4') { redraw_pts(x,y,n); menu(); graph_b_spline(x,y,n); }
            else if (key == 'a') { redraw_pts(x,y,n); menu(); graph_all(x,y,n); }
            else if (key == 'c') { redraw_pts(x,y,n); menu(); }
		    else if (key == 'q') exit(0);
        } while (key != 'r');
	}
}

