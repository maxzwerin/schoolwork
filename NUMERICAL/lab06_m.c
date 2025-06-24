#include "FPToolkit.c"

#define MAXDIST 7
#define SIZE 800
#define EXIT 30
#define STEP 0.001


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

    G_rgb(0.2,1.0,1.0);
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
}


void clear() {
	G_rgb(0,0,0);
	G_clear();
	G_rgb(1,0.3,0.3);
	G_fill_rectangle(0,0,SIZE,EXIT);
}


void redraw_pts(double x[], double y[], int n) {
    clear();
    for (int i = 0; i < n; i++) {
        G_rgb(1,1,1);
        G_circle(x[i], y[i], 4);
        if (i == n - 1) continue;
        G_rgb(0.2,0.2,0.6);
        G_line(x[i],y[i],x[i+1],y[i+1]);
    }
}


void menu() {
    int x = 30, y = SIZE-x, dy = 20;

    G_rgb(1,1,1);
    G_draw_string("click points to move them", x,y); y-=dy;
    G_draw_string("C to use new points",       x,y); y-=dy;
    G_draw_string("Q to quit",                 x,y);
}


void update(double x[], double y[], int n) {
    redraw_pts(x,y,n);
    menu();
    graph_b_spline(x,y,n);
    G_display_image();
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
        p[1] = SIZE - p[1];

        x[i] = p[0];
        y[i] = p[1];

        update(x,y,n);
    }
}

int move(double x[], double y[], int n) {
    int key, p[2],i,has_moved = 0;
    double pi[2];
    while (1) {
        key = G_wait_event(pi);
        if (key == 'c') return key;
		else if (key == 'q') exit(0);

        int mouse_state = S_mouse_coord_window(p);
        p[1] = SIZE - p[1]; // y coord is inversed

        update(x,y,n);

        if (mouse_state == 0 && has_moved) return key;

        else if (mouse_state == 1) {
            i = is_on_point(p, x,y,n);
            if (i >= 0) move_point(x,y,i,n);
            has_moved = 1;
        }
    }
    return key;
}



int main() {
    int key, n;
    double x[1000], y[1000], p[2];

	G_init_graphics(SIZE, SIZE);

	while (1) {
        clear();
        n = click_and_save(x, y);

        do {
            update(x,y,n);
            key = move(x,y,n);
		    if (key == 'q') exit(0);
        } while (key != 'c');
	}
}
