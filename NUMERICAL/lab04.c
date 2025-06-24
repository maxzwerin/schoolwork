#include "FPToolkit.c"

#define SIZE 800
#define EXIT 30

/* LAB04 :: BEST FIT LINE 
idea:
e0   (a + bx0) - y0
e1   (a + bx1) - y1
e2   (a + bx2) - y2
e3   (a + bx3) - y3
total error = e0^2 + e1^2 + e2^2 + e3^2  *standard*
total error = |e0| + |e1| + |e2| + |e3|  *alternative*
*/

double x[100], y[100];


void pder(int n, double *A, double *B)
// partial derivative
// line a + bx
{
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

    return;
}


int click_and_save()
{
	int n;
	double p[2];

	for (n = 0; ; n++) {
		G_wait_click(p);
		if (p[1] < EXIT) break;

		x[n] = p[0];
		y[n] = p[1];

		G_rgb(1,1,1);
		G_circle(p[0], p[1], 4);
	}
	return n;
}


void clear() 
{
	G_rgb(0,0,0);
	G_clear();
	G_rgb(1,0.3,0.3);
	G_fill_rectangle(0,0,SIZE,EXIT);
}


void draw()
{
	double p, A, B;

	clear();
	int n = click_and_save();

	G_rgb(0.5,0.5,1.0);
	
	pder(n, &A, &B);

	// y = a + bx
	double x0, y0, x1, y1;

	x0 = 0;
	y0 = A;
	x1 = SIZE;
	y1 = A + (B * SIZE);

	G_line(x0, y0, x1, y1); 
	
}


int main()
{
    int key;
	G_init_graphics(SIZE, SIZE);

	draw();

	while (1) {
		key = G_wait_key();
		if (key == 'q') exit(0);
		else draw();
	}
}
