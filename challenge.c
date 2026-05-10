#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
// #include <string.h> // Uncomment if 'string.h' is included in system.

// =========================================================
// FUNCTION PROTOTYPE
// =========================================================
void sparse__multiply(
    int rows,
    int cols,
    const double *A,
    const double *x,
    int *out_nnz,
    double *values,
    int *col_indices,
    int *row_ptrs,
    double *y);

// =========================================================
// TODO: USER IMPLEMENTATION
// =========================================================
void sparse_multiply(
    int rows, int cols, const double *A, const double *x,
    int *out_nnz, double *values, int *col_indices, int *row_ptrs,
    double *y)
{
  int i, j, k = 0;
  double current_value = 0.0;

  // Cleaning 'y'
  // memset(y, 0 , rows * sizeof(double)); // Uncomment if 'string.h' is included in system.
  // Comment cycle if 'string.h' is included in system, it could be benefitial.
  for (i = 0; i < rows; i++)
  {
    y[i] = 0.0;
  }

  // Compressed Sparse Row
  for (i = 0; i < rows; i++)
  {
    row_ptrs[i] = k;

    for (j = 0; j < cols; j++)
    {
      current_value = A[i * cols + j];
      if (current_value != 0)
      {
        values[k] = current_value;
        col_indices[k] = j;
        k += 1;
      }
    }
  }
  row_ptrs[rows] = k;

  *out_nnz = row_ptrs[rows];

  // Matrix-Vector product
  for (i = 0; i < rows; i++)
  {
    for (j = row_ptrs[i]; j < row_ptrs[i + 1]; j++)
    {
      y[i] += values[j] * x[col_indices[j]];
    }
  }
}

// =========================================================
// TEST HARNESS
// =========================================================
int main(void)
{
  srand(time(NULL));

  const int num_iterations = 100;
  int passed_count = 0;

  for (int iter = 0; iter < num_iterations; ++iter)
  {
    int rows = rand() % 41 + 5;
    int cols = rand() % 41 + 5;
    double density = 0.05 + (rand() / (double)RAND_MAX) * 0.35;

    size_t mat_sz = (size_t)rows * cols;

    double *A = calloc(mat_sz, sizeof(double));
    for (size_t i = 0; i < mat_sz; ++i)
    {
      if (((double)rand() / RAND_MAX) < density)
      {
        A[i] = ((double)rand() / RAND_MAX) * 20.0 - 10.0;
      }
    }

    double *values = malloc(mat_sz * sizeof(double));
    int *col_indices = malloc(mat_sz * sizeof(int));
    int *row_ptrs = malloc((rows + 1) * sizeof(int));
    double *x = malloc(cols * sizeof(double));
    double *y_user = malloc(rows * sizeof(double));
    double *y_ref = calloc(rows, sizeof(double));
    int out_nnz = 0;

    for (int i = 0; i < cols; ++i)
    {
      x[i] = ((double)rand() / RAND_MAX) * 20.0 - 10.0;
    }

    for (int i = 0; i < rows; ++i)
    {
      double sum = 0.0;
      for (int j = 0; j < cols; ++j)
      {
        sum += A[i * cols + j] * x[j];
      }
      y_ref[i] = sum;
    }

    sparse_multiply(rows, cols, A, x, &out_nnz, values, col_indices, row_ptrs, y_user);

    double max_err = 0.0;
    int passed = 1;
    for (int i = 0; i < rows; ++i)
    {
      double diff = fabs(y_user[i] - y_ref[i]);
      double tol = 1e-7 + 1e-7 * fabs(y_ref[i]); // Mixed absolute/relative tolerance
      if (diff > tol)
      {
        max_err = fmax(max_err, diff);
        passed = 0;
      }
    }

    if (passed)
    {
      passed_count++;
    }

    printf(
        "Iter %2d [%3dx%3d, density=%.2f, nnz=%4d]: %s (Max error: %.2e)\n",
        iter, rows, cols, density, out_nnz, passed ? "PASS" : "FAIL", max_err);

    free(A);
    free(values);
    free(col_indices);
    free(row_ptrs);
    free(x);
    free(y_user);
    free(y_ref);
  }

  printf(
      "\n%s (%d/%d iterations passed)\n",
      passed_count == num_iterations ? "All tests passed!" : "Some tests failed.",
      passed_count, num_iterations);

  return passed_count == num_iterations ? 0 : 1;
}