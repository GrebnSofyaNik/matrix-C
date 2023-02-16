#include "s21_matrix.h"

int is_valid(matrix_t *matrix) {
  int res = FAILURE;
  if (matrix && matrix->rows > 0 && matrix->columns > 0 && matrix->matrix)
    res = SUCCESS;
  return res;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int res = OK;
  result->rows = rows;
  result->columns = columns;
  if (rows < 1 || columns < 1) {
    result->matrix = NULL;
    res = INCORRECT_MATRIX;
  } else {
    result->matrix = (double **)calloc(rows, sizeof(double *));
    if (result->matrix) {
      for (int i = 0; i < rows; i++) {
        result->matrix[i] = (double *)calloc(columns, sizeof(double));
      }
    } else {
      res = INCORRECT_MATRIX;
    }
  }
  return res;
}

void s21_remove_matrix(matrix_t *A) {
  if (A != NULL) {
    if (A->matrix) {
      for (int i = 0; i < A->rows; i++)
        if (A->matrix[i]) free(A->matrix[i]);
      free(A->matrix);
    }
    A->columns = 0;
    A->rows = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int res = SUCCESS;
  if (is_valid(A) && is_valid(B)) {
    if (A->rows == B->rows && A->columns == B->columns) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          if (fabs(A->matrix[i][j] - B->matrix[i][j]) >= EPS) res = FAILURE;
        }
      }
    } else {
      res = FAILURE;
    }
  } else {
    res = FAILURE;
  }
  return res;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  s21_create_matrix(A->rows, A->columns, result);
  if (is_valid(A) && is_valid(B)) {
    if (A->rows == B->rows && A->columns == B->columns) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    } else {
      res = CALCULATION_ERROR;
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  s21_create_matrix(A->rows, A->columns, result);
  if (is_valid(A) && is_valid(B)) {
    if (A->rows == B->rows && A->columns == B->columns) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    } else {
      res = CALCULATION_ERROR;
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int res = OK;
  s21_create_matrix(A->rows, A->columns, result);
  if (is_valid(A)) {
    if (is_valid(result)) {
      for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->columns; j++) {
          if (number == 0)
            result->matrix[i][j] = 0;
          else
            result->matrix[i][j] = number * A->matrix[i][j];
        }
      }
    } else {
      res = CALCULATION_ERROR;
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  s21_create_matrix(A->rows, B->columns, result);
  if (is_valid(A) && is_valid(B)) {
    if (A->columns == B->rows) {
      for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->columns; j++) {
          result->matrix[i][j] = 0;
          for (int k = 0; k < A->columns; k++) {
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
          }
        }
      }
    } else {
      res = CALCULATION_ERROR;
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int res = OK;
  s21_create_matrix(A->columns, A->rows, result);
  if (is_valid(A)) {
    if (is_valid(result)) {
      for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->columns; j++) {
          result->matrix[i][j] = A->matrix[j][i];
        }
      }
    } else {
      res = CALCULATION_ERROR;
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int res = OK;
  s21_create_matrix(A->rows, A->columns, result);
  if (is_valid(A)) {
    if (A->rows == A->columns) {
      if (A->rows == 1) {
        result->matrix[0][0] = 1 / A->matrix[0][0];
      } else {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            double d = 0;
            if (A->rows - 1 == 1) {
              d = A->matrix[!i][!j];
            } else {
              matrix_t minor;
              s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
              res |= create_minor(A, i, j, &minor);
              res |= s21_determinant(&minor, &d);
              s21_remove_matrix(&minor);
            }
            result->matrix[i][j] = ((i + j) % 2 == 1) ? d * (-1) : d;
          }
        }
      }
    } else {
      res = CALCULATION_ERROR;
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

int create_minor(matrix_t *A, int i, int j, matrix_t *result) {
  int res = OK;
  if (is_valid(A)) {
    int m = 0, n = 0;
    for (int k = 0; k < A->rows; k++) {
      if (k == i) continue;
      for (int l = 0; l < A->columns; l++) {
        if (l == j) continue;
        result->matrix[m][n] = A->matrix[k][l];
        n++;
      }
      m++;
      n = 0;
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

int s21_determinant(matrix_t *A, double *result) {
  int res = OK;
  if (is_valid(A) && result) {
    *result = 0;
    if (A->rows == A->columns) {
      if (A->rows == 1) {
        *result = A->matrix[0][0];
      } else if (A->rows == 2) {
        *result = A->matrix[0][0] * A->matrix[1][1] -
                  A->matrix[0][1] * A->matrix[1][0];
      } else {
        for (int i = 0; i < A->rows; i++) {
          int a = (i % 2 == 0) ? 1 : -1;
          matrix_t minor;
          double d = 0;
          s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
          res |= create_minor(A, 0, i, &minor);
          res |= s21_determinant(&minor, &d);
          *result += a * A->matrix[0][i] * d;
          s21_remove_matrix(&minor);
        }
      }
    } else {
      res = CALCULATION_ERROR;
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int res = OK;
  s21_create_matrix(A->rows, A->columns, result);
  if (is_valid(A)) {
    if (A->rows == A->columns) {
      double d = 0;
      res |= s21_determinant(A, &d);
      if (d == 0) {
        res = CALCULATION_ERROR;
      } else {
        d = 1 / d;
        matrix_t temp, temp1, temp2;
        res |= s21_calc_complements(A, &temp);
        res |= s21_transpose(&temp, &temp1);
        res |= s21_mult_number(&temp1, d, &temp2);
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            result->matrix[i][j] = temp2.matrix[i][j];
          }
        }
        s21_remove_matrix(&temp);
        s21_remove_matrix(&temp1);
        s21_remove_matrix(&temp2);
      }
    } else {
      res = CALCULATION_ERROR;
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}