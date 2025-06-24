#include "FPToolkit.c"

/* LAB05 :: SPLINES

   a series of connected cubics

   idea:
   s1(x) = q0 + (q10/p10)(x-p0)  A1(x-p0)(x-p1) + B1(x-p0)^2(x-p1)
   s2(x) = q1 + (q21/p21)(x-p1)  A2(x-p1)(x-p2) + B2(x-p1)^2(x-p2)
   s3(x) = q2 + (q32/p32)(x-p2)  A3(x-p2)(x-p3) + B2(x-p2)^2(x-p3)
 
   s1'(p1) should equal s2'(p1) to ensure line is smooth
   s1"(p1) should equal s2"(p1) to ensure line is smooth
   s"(endpoints) should equal 0
 
   tri-diagonal system - allows for O(n) vs gaussian O(n^3)
 
   A1 --- B1 --- A2 --- B2 --- A3 --- B3 --- A4 --- B4 --- A5 --- B5     |
   ----------------------------------------------------------------------------------------------
   1      -p10                                                           | 0
   p10    p10^2  p21                                                     | q21/p21 - q10/p10
          p10^2  -p20   p10p21                                           | - (...)
                 p21    p21^2  p32                                       | q32/p32 - q21/p21
                        p21^2  -p31   p21p32                             | - (...)
                               p32    p32^2  p43                         | q43/p43 - q32/p32
                                      p32^2  -p42   p32p43               | - (...)
                                             p43    p43^2  p54           | q54/p54 - q43/p43
                                                    p43^2  -p53   p43p54 | - (...)
                                                           1      2p54   | 0

   note : For the interpolating polynomial, the order of the x[i]'s  is immaterial but
   for the natural spline, the order of x[i]'s IS important and one would normally want
   x[0] < x[1] < x[2] < x[3] < .... although this program will work with them as given
   and one can use the program to see why one would want the x[i]'s ordered.

   The use of the natural spline is a bad idea for designing non-functional shapes.
   When the direction of x switches, in order to make the derivatives match, the graph
   will arch back the way it came.
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

    s = natural_spline(x, y, n,
                       M, A, B);
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
    G_draw_string("SPLINE",200,SIZE-60);
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

    // y = a + bx

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
    G_draw_string("SPLINE",200,SIZE-60);

	G_rgb(0.7,0.3,0.7);
    G_fill_rectangle(5,SIZE-170,180,165);
    G_rgb(1,1,1);
    G_draw_string("1 to display lagrange",      x,y); y-=dy;
    G_draw_string("2 to dsiplay best fit line", x,y); y-=dy;
    G_draw_string("3 to display spline",        x,y); y-=dy;
    G_draw_string("A to display all",           x,y); y-=dy;
    G_draw_string("C to display none",          x,y); y-=dy;
                                                      y-=dy;
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
    graph_lagrange(x, y, n);
    graph_best_fit_line(x, y, n);
    graph_natural_spline(x, y, n);
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
            else if (key == 'a') { redraw_pts(x,y,n); menu(); graph_all(x,y,n); }
            else if (key == 'c') { redraw_pts(x,y,n); menu(); }
		    else if (key == 'q') exit(0);
        } while (key != 'r');
	}
}
