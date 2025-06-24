#include "FPToolkit.c"

#define SIZE 800
#define EXIT 30
#define STEP 0.001

/* LAB03 :: LAGRANGE FORM
p3(x) = y0(x-x1)(x-x2)(x-x3)     y1(x-x0)(x-x2)(x-x3)
        --------------------  +  --------------------  + ...
        (x0-x1)(x0-x2)(x0-x3)    (x1-x0)(x1-x2)(x1-x3)
*/

double x[100], y[100];

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


double lagrange(int n, double X)
{
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

void clear() 
{
	G_rgb(0,0,0);
	G_clear();
	G_rgb(1,0.3,0.3);
	G_fill_rectangle(0,0,SIZE,EXIT);
}

void draw()
{
	clear();
	int n = click_and_save();

	double p, X, Y;

	G_rgb(0.5,0.5,1.0);
	for (X = 0; X < SIZE; X += STEP) {
		Y = lagrange(n, X);
		G_point(X, Y);
	}
}


int main()
{
    int key;
	G_init_graphics(SIZE, SIZE);

	draw();

	while (1) {
		key = G_wait_key();
		if (key == 'q') exit(0);
		else if (key == 'c') draw();
	}
}
