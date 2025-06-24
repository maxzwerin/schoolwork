// Create an object of revolution
// Outputs object to file XYZ/output.xyz

#include "FPToolkit.c"
#include "M3d_matrix_tools.c"


int window_size = 800;


void clear_graphics() {
   G_rgb(0, 0.2, 0.2);
   G_clear();
   G_rgb(1, 0, 0);
   G_fill_rectangle(0, window_size - 40, window_size, 40);
}


int click_and_save(double x[], double y[]) {
   int n = 0;
   double p[2];
   while (1) {
       G_wait_click(p);
       if (p[1] > window_size - 40) break;
       x[n] = p[0];
       y[n] = p[1];
       n++;
       G_rgb(1, 1, 1);
       G_fill_circle(p[0], p[1], 2);
   }
   return n;
}


int main() {
   FILE *g;
   double x[1000], y[1000], z[1000];
   double xt[1000], yt[1000], zt[1000];
   double M[4][4], V[4][4];
   int i, j, n, slices;


   g = fopen("XYZ/output.xyz", "w");
   if (g == NULL) {
       printf("Failed to open output file.\n");
       exit(0);
   }


   G_init_graphics(window_size, window_size);
   clear_graphics();


   n = click_and_save(x, y);


   for (i = 0; i < n; i++) {
       z[i] = 0; // initial z points at 0
   }


   printf("\nEnter number of slices: ");
   scanf("%d", &slices);
   double radians = (2 * M_PI) / slices;


   fprintf(g, "%d\n", n * slices);


   for (i = 0; i < slices; i++) {
       double cos_theta = cos(i * radians);
       double sin_theta = sin(i * radians);


       M3d_make_x_rotation_cs(M, cos_theta, sin_theta);


       for (j = 0; j < n; j++) {
           xt[j] = x[j];
           yt[j] = y[j];
           zt[j] = z[j];
       }


       M3d_mat_mult_points(xt, yt, zt, M, x, y, z, n);


       for (j = 0; j < n; j++) {
           fprintf(g, "%lf %lf %lf\n", xt[j] / 100, yt[j] / 100, zt[j] / 100);
       }
   }


   int num_polys = (n - 1) * slices;
   fprintf(g, "%d\n", num_polys);


   for (i = 0; i < slices; i++) {
       for (j = 0; j < n - 1; j++) {
           int a = i * n + j;
           int b = ((i + 1) % slices) * n + j;
           int c = ((i + 1) % slices) * n + (j + 1);
           int d = i * n + (j + 1);
           fprintf(g, "4  %d %d %d %d\n", a, b, c, d);
       }
   }


   fclose(g);
   printf("Model saved to output.xyz\n");


   return 0;
}
