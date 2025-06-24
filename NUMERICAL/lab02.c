#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/** LAB02 :: GAUSSIAN ELIMINATION
  * 
  * Solve systems of equations by converting
  * a matrix to upper triangular form, and then
  * using back substitution to get results for 
  * each x.
  * 
  * Compile        : gcc lab02.c -lm 
  * Run            : ./a.out < lab02_tests/p03.in 
  * Test against   :           lab02_tests/p03.res
  *
  */

double m[20][20];
int degree;

void swap_rows(int row1, int row2) 
{
  double temp[20];
  for (int i = 0; i <= degree; i++) {
    temp[i] = m[row1][i];
    m[row1][i] = m[row2][i];
    m[row2][i] = temp[i];
  }
}

void reshuffler(int x)
// find the largest number in the column
// for row x. swap row x with largest number
{
  int max_row = x;
  double max_value = fabs(m[x][x]);

  for (int i = x + 1; i < degree; i++) {
    if (fabs(m[i][x]) > max_value) {
      max_value = fabs(m[i][x]);
      max_row = i;
    }
  }

  if (max_value == 0.0) {
    printf("no unique solution\n");
    exit(0);
  }

  if (max_row != x) swap_rows(x, max_row);
}

int main()
{
  scanf("%d", &degree);
  if (degree > 20) {
    printf("error: too many unknowns\n");
    exit(0);
  }

  for (int i = 0; i < degree; i++) {
    for (int j = 0; j < degree + 1; j++) {
      scanf("%lf", &m[i][j]);
    }
  }

  printf("\n");
  printf("initial m (with degree %d) : \n", degree);
  for (int i = 0; i < degree; i++) {
    for (int j = 0; j < degree + 1; j++) {
      if (j == degree) printf("|%10.2lf\n", m[i][j]);
      else printf("%10.2lf ", m[i][j]);
    }
  }

  // upper triangular form
  int count = 1;
  for (int i = count - 1; i < degree; i++) {
    reshuffler(i);
    for (int j = count; j < degree; j++) {

      double x = m[j][i] / m[i][i];

      for (int k = i; k < degree + 1; k++) {
          m[j][k] = m[j][k] - m[i][k] * x;
      }
    }
    count++;
  }
  printf("\n");

  printf("upper triangular form\n");
  for (int i = 0; i < degree; i++) {
    for (int j = 0; j < degree + 1; j++) {
      if (j == degree) printf("|%10.2lf\n", m[i][j]);
      else printf("%10.2lf ", m[i][j]);
    }
  }

  // back substitution
  double ans[20];
  for (int i = degree - 1; i >= 0; i--) {
    ans[i] = m[i][degree];
    for (int j = i + 1; j < degree; j++) {
      ans[i] -= m[i][j] * ans[j];
    }
    ans[i] /= m[i][i];
  }

  printf("\n");
  for (int i = 0; i < degree; i++) {
    printf("+++ x%d : %20.16lf\n", i + 1, ans[i]);
  }
  printf("\n");

  return 0;
}
