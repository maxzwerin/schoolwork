#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

// LAB01 :: FIND ALL ROOTS FROM COMPLEX FUNCTION

#define TOLERANCE 1e-10

void print_complex_number(complex c)
{
  double a,b ;
  a = creal(c) ;
  b = cimag(c) ;
  if (b >= 0) {
    printf("%22.16lf + %22.16lf I ",a,fabs(b)) ;
  } else {
    printf("%22.16lf - %22.16lf I ",a,fabs(b)) ;
  }
  printf("\n");
}

void evaluate_poly(complex poly[], int degree, complex x, complex *p, complex *dp) 
// Evaluate polynomial (p) and derivative (dp) with synthetic division
{
    complex b = poly[0];
    complex db = 0.0;
    for (int i = 1; i <= degree; i++) {
        db = db * x + b;
        b = b * x + poly[i];
    }
    *p = b;
    *dp = db;
}

complex newton_method(complex poly[], int degree, complex guess) 
// Newton's method to find a root
{
    complex x = guess;
    for (int i = 0; i < 50; i++) {
        complex p, dp;
        evaluate_poly(poly, degree, x, &p, &dp);
        if (cabs(p) < TOLERANCE) break; // Root found (close enough to zero)
        if (cabs(dp) < TOLERANCE) { // Dividing by zero bad
            x += 1.0 + 1.0 * I;  
            continue;
        }
        x -= p / dp;
    }
    return x;
}

void synthetic_division(complex poly[], int degree, complex root, complex new_poly[]) 
// Synthetic division to deflate polynomial
{
    new_poly[0] = poly[0];
    for (int i = 1; i < degree; i++) {
        new_poly[i] = new_poly[i - 1] * root + poly[i];
    }
}

int compare_complex(const void *a, const void *b) 
// Compare two complex numbers by magnitude
{
    complex c1 = *(complex*)a;
    complex c2 = *(complex*)b;
    
    double mag1 = cabs(c1); 
    double mag2 = cabs(c2);
    
    if (mag1 < mag2) return -1;
    if (mag1 > mag2) return 1;
    return 0;  
}

void print_roots(complex mag[], int n) 
{
    for (int i = 0; i < n; i++) {
        if (i + 1 == 1) printf("%4dst root found: ", i + 1); 
        else if (i + 1 == 2) printf("%4dnd root found: ", i + 1); 
        else if (i + 1 == 3) printf("%4drd root found: ", i + 1); 
        else printf("%4dth root found: ", i + 1); 
        print_complex_number(mag[i]);       
    } 
    printf("\n");
}

int main() {
    int degree;
    scanf("%d", &degree);
    
    complex poly[degree + 1];
    for (int i = 0; i <= degree; i++) {
        double real;
        scanf("%lf", &real);
        poly[i] = real;
    }
    
    printf("\n");
    complex mag[degree];
    for (int i = degree; i > 0; i--) {
        complex guess = 0.4 + 0.9 * I; // Must have a non-zero imaginary guess
        complex root = newton_method(poly, i, guess);
        mag[degree - i] = root;

        // Deflate polynomial
        complex new_poly[i];
        synthetic_division(poly, i, root, new_poly);
        for (int j = 0; j < i; j++) {
            poly[j] = new_poly[j];
        }
    }
    
    int n = sizeof(mag) / sizeof(mag[0]);
    qsort(mag, n, sizeof(complex), compare_complex);
    print_roots(mag, n);
}
