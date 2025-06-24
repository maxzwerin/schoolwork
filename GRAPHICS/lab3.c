#include "FPToolkit.c"

int window_size = 800 ; 

void clear_graphics()
{
  G_rgb(0,0,0) ; 
  G_clear() ; 
  G_rgb(1,0,0) ; 
  G_fill_rectangle(0,0,800,40) ;
}


int click_and_save(double x[], double y[]) {
    int n = 0;
    double p[2];

    while (1) {
        G_wait_click(p);
        if (p[1] < 40) break;
        x[n] = p[0];
        y[n] = p[1];
        n++;
       
        G_rgb(1, 1, 1);
        G_fill_circle(p[0], p[1], 2);
    }

    return n;
}

int in_out (double x[], double y[], int n, double P[2], int I, int J)
// return 1 if point P is on the correct side of the line
// using the convex polygon to find center mass
// else return 0
{
  int i ;
  double MX = 0, MY = 0 ;
  for (i = 0 ; i < n ; i++) {
    MX += x[i] ; MY += y[i] ;
  }
  
  MX /= n ; MY /= n ;

  double a = (y[J] - y[I]);
  double b = (x[J] - x[I]);
  double c = (x[J] - x[I]) * y[I] - (y[J] - y[I]) * x[I] ;
    
  double signM = a * MX - b * MY + c ;
  double signP = a * P[0] - b * P[1] + c ;
      
  if (signM > 0 && signP < 0) return 0 ; 
  if (signM < 0 && signP > 0) return 0 ; 
  
  return 1 ; 
}


int intersect_2_lines (double A[2], double B[2],
                       double C[2], double D[2],
                       double intersection[2])
// return 0 if lines do NOT intersect
// return 1 if they do  
{
  double m1, m2, b1, b2 ; 
  
  if (A[0] == B[0] && C[0] == D[0]) return 0; // both lines are vertical
  
  if (A[0] == B[0]) { // line AB is vertical
    
    intersection[0] = A[0];
    m2 = (D[1] - C[1]) / (D[0] - C[0]);
    b2 = C[1] - m2 * C[0];

    intersection[1] = m2 * intersection[0] + b2;
 
  } else if (C[0] == D[0]) { // line CD is vertical

    intersection[0] = C[0];
    m1 = (B[1] - A[1]) / (B[0] - A[0]);
    b1 = A[1] - m1 * A[0];

    intersection[1] = m1 * intersection[0] + b1;
    
  } else { // normal lines

    m1 = (B[1] - A[1]) / (B[0] - A[0]);
    m2 = (D[1] - C[1]) / (D[0] - C[0]);

    if (m1 == m2) return 0; // lines are parallel

    b1 = A[1] - m1 * A[0];
    b2 = C[1] - m2 * C[0];

    intersection[0] = (b2 - b1) / (m1 - m2);
    intersection[1] = m1 * intersection[0] + b1;
  }
  
  return 1;
}


int clip_line(double A[], double B[], int NAB, 
	            double X[], double Y[], int NXY, 
	            int a, int b) 
// AB is the clipping polygon
// XY is the original polygon
// clip XY using line A[a],B[a] to A[b],B[b]
{
  double P[2], intersection[2] ;
  double aa[2], bb[2], xx[2], yy[2] ; 
  int inout1, inout2, i, j ; 

  double tempx[1000], tempy[1000] ;
  int ntemp = 0 ; 

  aa[0] = A[a] ; aa[1] = B[a] ;
  bb[0] = A[b] ; bb[1] = B[b] ; 

  G_rgb(0,1,1) ; 
  G_polygon(A, B, NAB) ; 
  
  G_rgb(1,0,1) ;
  G_line(A[a],B[a], A[b],B[b]) ;
  
  for (i = 0 ; i < NXY; i++) {

    j = (i + 1) % NXY ; 

    xx[0] = X[i] ; xx[1] = Y[i] ;
    yy[0] = X[j] ; yy[1] = Y[j] ;

    P[0] = X[i] ; P[1] = Y[i] ; 
    inout1 = in_out(A, B, NAB, P, a, b) ; 

    P[0] = X[j] ; P[1] = Y[j] ; 
    inout2 = in_out(A, B, NAB, P, a, b) ; 

    if (inout1 == 1 && inout2 == 1) { 
      // good to good
      tempx[ntemp] = X[j] ;
      tempy[ntemp] = Y[j] ;
      ntemp++ ;

    } else if (inout1 == 1 && inout2 == 0) {
      // good to bad
      intersect_2_lines(aa,bb, xx,yy, intersection) ; 
      G_rgb(0,1,0) ; 
      G_fill_circle(intersection[0], intersection[1], 3) ;
      tempx[ntemp] = intersection[0] ;
      tempy[ntemp] = intersection[1] ;
      ntemp++ ;

    } else if (inout1 == 0 && inout2 == 1) {
      // bad to good
      intersect_2_lines(aa,bb, xx,yy, intersection) ; 
      G_rgb(0,1,0) ; 
      G_fill_circle(intersection[0], intersection[1], 3) ;
      tempx[ntemp] = intersection[0] ;
      tempy[ntemp] = intersection[1] ;
      ntemp++ ;
      tempx[ntemp] = X[j] ;
      tempy[ntemp] = Y[j] ;
      ntemp++ ; 

    }

  }

  G_rgb(1,0,0) ; 
  G_polygon(X, Y, NXY) ; 

  G_rgb(0,1,1) ; 
  G_polygon(A, B, NAB) ; 

  G_rgb(0,1,0) ;
  G_polygon(tempx, tempy, ntemp) ;

  for (i = 0 ; i < ntemp ; i++) {
    X[i] = tempx[i] ;
    Y[i] = tempy[i] ;
  }

  return ntemp ; 
}

int main()
{
  double A[1000], B[1000] ;   // clipping polygon
  double X[1000], Y[1000] ;   // clipped  polygon
  double cx[1000], cy[1000] ; // original polygon
  int NAB, NXY, cxy, i, j ; 
  char q ;
  
  G_init_graphics(window_size, window_size) ;
  clear_graphics() ; 
  
  cxy = click_and_save(cx, cy); 

  while (1) {
    clear_graphics();
    G_rgb(1, 0, 0); 
    G_polygon(cx, cy, cxy);
    
    NAB = click_and_save(A, B); 
    G_rgb(0, 1, 1); 
    G_polygon(A, B, NAB);
    
    NXY = cxy;
    for (i = 0; i < cxy; i++) {
      X[i] = cx[i];
      Y[i] = cy[i];
    }
    
    for (i = 0; i < NAB; i++) {
      j = (i + 1) % NAB;
      NXY = clip_line(A, B, NAB, X, Y, NXY, i, j);
    }

    clear_graphics();
    G_rgb(0.3, 0, 0);
    G_polygon(cx, cy, cxy);
    G_rgb(0, 0.3, 0.3); 
    G_polygon(A, B, NAB); 
    G_rgb(0, 1, 0); 
    G_fill_polygon(X, Y, NXY);   

    q = G_wait_key();
    if (q == 113) break; 
  }

  return 0; 
}