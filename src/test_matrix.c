#include <check.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "s21_matrix.h"

START_TEST(test_create_remove_1) {
  matrix_t A;
  ck_assert(0 == s21_create_matrix(3, 3, &A));
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < A.columns; j++) {
      ck_assert(A.matrix[i][j] == 0);
    }
  }
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_create_remove_2) {
  matrix_t A;
  ck_assert(INCORRECT_MATRIX == s21_create_matrix(0, 2, &A));
}
END_TEST

START_TEST(test_determinant_1) {
  double det1, det2;
  matrix_t A;
  ck_assert(OK == s21_create_matrix(2, 2, &A));
  A.matrix[0][0] = 1.6;
  A.matrix[0][1] = -1;
  A.matrix[1][0] = -6745;
  A.matrix[1][1] = 23;
  det1 = -6708.2;
  ck_assert(OK == s21_determinant(&A, &det2));
  ck_assert(SUCCESS == (det1 == det2));
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_2) {
  double det1, det2;
  matrix_t A;
  ck_assert(OK == s21_create_matrix(2, 2, &A));
  A.matrix[0][0] = 0;
  A.matrix[0][1] = 0;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 0;
  det1 = 0;
  ck_assert(OK == s21_determinant(&A, &det2));
  ck_assert(SUCCESS == (det1 == det2));
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_3) {
  double det1, det2;
  matrix_t A;
  ck_assert(OK == s21_create_matrix(2, 2, &A));
  A.matrix[0][0] = 1.1;
  A.matrix[0][1] = 1.1;
  A.matrix[1][0] = 1.1;
  A.matrix[1][1] = 1.1;
  det1 = 0;
  ck_assert(OK == s21_determinant(&A, &det2));
  ck_assert(SUCCESS == (det1 == det2));
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_4) {
  double det1, det2;
  matrix_t A;
  ck_assert(OK == s21_create_matrix(3, 3, &A));
  A.matrix[0][0] = 1.1;
  A.matrix[0][1] = 1.1;
  A.matrix[0][2] = 5;
  A.matrix[1][0] = 1.1;
  A.matrix[1][1] = 1.1;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 1;
  A.matrix[2][1] = 9;
  A.matrix[2][2] = 3;
  det1 = 26.4;
  ck_assert(OK == s21_determinant(&A, &det2));
  ck_assert(SUCCESS == (det1 == det2));
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_5) {
  double det1, det2;
  matrix_t A;
  ck_assert(OK == s21_create_matrix(3, 3, &A));
  A.matrix[0][0] = 0;
  A.matrix[0][1] = 3;
  A.matrix[0][2] = 23;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 3;
  A.matrix[2][0] = 4;
  A.matrix[2][1] = 123;
  A.matrix[2][2] = 5;
  det1 = 16644;
  ck_assert(OK == s21_determinant(&A, &det2));
  ck_assert(SUCCESS == (det1 == det2));
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_6) {
  double det1, det2;
  matrix_t A;
  ck_assert(OK == s21_create_matrix(3, 3, &A));
  A.matrix[0][0] = 0;
  A.matrix[0][1] = 0;
  A.matrix[0][2] = 0;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 0;
  A.matrix[1][2] = 0;
  A.matrix[2][0] = 0;
  A.matrix[2][1] = 0;
  A.matrix[2][2] = 0;
  det1 = 0;
  ck_assert(OK == s21_determinant(&A, &det2));
  ck_assert(SUCCESS == (det1 == det2));
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_7) {
  double det1, det2;
  matrix_t A;
  ck_assert(OK == s21_create_matrix(4, 4, &A));
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[0][3] = 4;
  A.matrix[1][0] = 1;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 9;
  A.matrix[1][3] = 25;
  A.matrix[2][0] = 4;
  A.matrix[2][1] = 3;
  A.matrix[2][2] = 2;
  A.matrix[2][3] = 1;
  A.matrix[3][0] = 25;
  A.matrix[3][1] = 9;
  A.matrix[3][2] = 4;
  A.matrix[3][3] = 1;
  det1 = -585;
  ck_assert(OK == s21_determinant(&A, &det2));
  ck_assert(SUCCESS == (det1 == det2));
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_determinant_8) {
  ck_assert(INCORRECT_MATRIX == s21_determinant(NULL, NULL));
}
END_TEST

START_TEST(test_determinant_9) {
  double det;
  matrix_t A;
  ck_assert(OK == s21_create_matrix(2, 3, &A));
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  A.matrix[1][2] = 2;
  ck_assert(CALCULATION_ERROR == s21_determinant(&A, &det));
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_calc_complements_1) {
  matrix_t A, B, C;
  ck_assert(OK == s21_create_matrix(2, 2, &A));
  A.matrix[0][0] = 4;
  A.matrix[0][1] = 3;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;

  ck_assert(OK == s21_create_matrix(2, 2, &B));
  B.matrix[0][0] = 3;
  B.matrix[0][1] = -6;
  B.matrix[1][0] = -3;
  B.matrix[1][1] = 4;

  ck_assert(OK == s21_calc_complements(&A, &C));
  ck_assert(SUCCESS == s21_eq_matrix(&B, &C));

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST

START_TEST(test_calc_complements_2) {
  matrix_t A, B, C;
  ck_assert(OK == s21_create_matrix(3, 3, &A));
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = 2;
  A.matrix[2][2] = 1;

  ck_assert(OK == s21_create_matrix(3, 3, &B));
  B.matrix[0][0] = 0;
  B.matrix[0][1] = 10;
  B.matrix[0][2] = -20;
  B.matrix[1][0] = 4;
  B.matrix[1][1] = -14;
  B.matrix[1][2] = 8;
  B.matrix[2][0] = -8;
  B.matrix[2][1] = -2;
  B.matrix[2][2] = 4;

  ck_assert(OK == s21_calc_complements(&A, &C));
  ck_assert(SUCCESS == s21_eq_matrix(&B, &C));

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST

START_TEST(test_calc_complements_3) {
  matrix_t A, result, res_origin;
  s21_create_matrix(2, 4, &A);
  s21_create_matrix(2, 2, &res_origin);
  int check = s21_calc_complements(&A, &result);
  int check_origin = 2;
  ck_assert_int_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(test_calc_complements_4) {
  matrix_t A, result, res_origin;
  s21_create_matrix(-2, 2, &A);
  s21_create_matrix(2, 2, &res_origin);
  int check = s21_calc_complements(&A, &result);
  int check_origin = 1;
  ck_assert_uint_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(test_inverse_1) {
  matrix_t A, A1, A2;
  ck_assert(OK == s21_create_matrix(2, 2, &A));
  ck_assert(OK == s21_create_matrix(2, 2, &A1));

  A.matrix[0][0] = 4;
  A.matrix[0][1] = 3;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;

  A1.matrix[0][0] = -0.5;
  A1.matrix[0][1] = 0.5;
  A1.matrix[1][0] = 1;
  A1.matrix[1][1] = -2.0 / 3;

  ck_assert(OK == s21_inverse_matrix(&A, &A2));
  ck_assert(SUCCESS == s21_eq_matrix(&A1, &A2));

  s21_remove_matrix(&A);
  s21_remove_matrix(&A1);
  s21_remove_matrix(&A2);
}
END_TEST

START_TEST(test_inverse_2) {
  matrix_t A, A1, A2;
  ck_assert(OK == s21_create_matrix(2, 2, &A));
  ck_assert(OK == s21_create_matrix(2, 2, &A1));

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3;
  A.matrix[1][1] = 4;

  A1.matrix[0][0] = -2;
  A1.matrix[0][1] = 1;
  A1.matrix[1][0] = 1.5;
  A1.matrix[1][1] = -0.5;

  ck_assert(OK == s21_inverse_matrix(&A, &A2));
  ck_assert(SUCCESS == s21_eq_matrix(&A1, &A2));

  s21_remove_matrix(&A);
  s21_remove_matrix(&A1);
  s21_remove_matrix(&A2);
}
END_TEST

START_TEST(test_inverse_3) {
  matrix_t A, B;
  ck_assert(OK == s21_create_matrix(2, 2, &A));

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 0;

  ck_assert(CALCULATION_ERROR == s21_inverse_matrix(&A, &B));

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_inverse_4) {
  matrix_t A, B;

  ck_assert(OK == s21_create_matrix(2, 3, &A));
  ck_assert(CALCULATION_ERROR == s21_inverse_matrix(&A, &B));

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_inverse_5) {
  matrix_t A, B;

  ck_assert(OK == s21_create_matrix(3, 3, &A));
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 4;
  A.matrix[1][1] = 5;
  A.matrix[1][2] = 6;
  A.matrix[2][0] = 7;
  A.matrix[2][1] = 8;
  A.matrix[2][2] = 9;
  ck_assert(CALCULATION_ERROR == s21_inverse_matrix(&A, &B));

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_inverse_6) {
  matrix_t A, result;
  s21_create_matrix(-5, 5, &A);
  int check = s21_inverse_matrix(&A, &result);
  int check_origin = 1;
  ck_assert_uint_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_eq_matrix) {
  matrix_t A, B, C;

  ck_assert(OK == s21_create_matrix(3, 3, &A));
  ck_assert(OK == s21_create_matrix(3, 3, &B));
  ck_assert(OK == s21_create_matrix(3, 3, &C));

  ck_assert(SUCCESS == s21_eq_matrix(&A, &B));

  for (int i = 0; i < A.rows; ++i) {
    for (int j = 0; j < A.columns; ++j) {
      A.matrix[i][j] = 1;
      B.matrix[i][j] = 1;
      C.matrix[i][j] = 2;
    }
  }

  ck_assert(SUCCESS == s21_eq_matrix(&A, &B));
  ck_assert(FAILURE == s21_eq_matrix(&A, &C));
  ck_assert(FAILURE == s21_eq_matrix(&B, &C));

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST

START_TEST(test_mult_number_1) {
  matrix_t A, B, result;
  ck_assert(OK == s21_create_matrix(3, 3, &A));
  ck_assert(OK == s21_create_matrix(3, 3, &B));

  for (int i = 0; i < A.rows; ++i) {
    for (int j = 0; j < A.columns; ++j) {
      A.matrix[i][j] = 1;
      B.matrix[i][j] = 2;
    }
  }
  ck_assert(OK == s21_mult_number(&A, 2, &result));
  ck_assert(SUCCESS == s21_eq_matrix(&result, &B));

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_mult_matrix_1) {
  matrix_t A, B, C1, C2;

  ck_assert(OK == s21_create_matrix(2, 2, &A));
  A.matrix[0][0] = 1;
  A.matrix[1][0] = 1;
  ck_assert(OK == s21_create_matrix(2, 2, &B));
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 3;
  B.matrix[1][0] = 2;
  B.matrix[1][1] = 4;
  ck_assert(OK == s21_mult_matrix(&A, &B, &C1));

  ck_assert(OK == s21_create_matrix(2, 2, &C2));
  C2.matrix[0][0] = 1;
  C2.matrix[0][1] = 3;
  C2.matrix[1][0] = 1;
  C2.matrix[1][1] = 3;
  ck_assert(SUCCESS == s21_eq_matrix(&C1, &C2));

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C1);
  s21_remove_matrix(&C2);
}
END_TEST

START_TEST(test_mult_matrix_2) {
  matrix_t A, B, result;
  int check_origin = 1;
  s21_create_matrix(-12, 3, &A);
  s21_create_matrix(3, 12, &B);
  int check = s21_mult_matrix(&A, &B, &result);
  ck_assert_int_eq(check, check_origin);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_transpose_1) {
  matrix_t A, AT1, AT2;

  ck_assert(OK == s21_create_matrix(2, 2, &A));
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 3;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = 4;
  ck_assert(OK == s21_create_matrix(2, 2, &AT1));
  AT1.matrix[0][0] = 1;
  AT1.matrix[0][1] = 2;
  AT1.matrix[1][0] = 3;
  AT1.matrix[1][1] = 4;
  ck_assert(OK == s21_transpose(&A, &AT2));
  ck_assert(SUCCESS == s21_eq_matrix(&AT1, &AT2));

  s21_remove_matrix(&A);
  s21_remove_matrix(&AT1);
  s21_remove_matrix(&AT2);
}
END_TEST

START_TEST(test_transpose_2) {
  matrix_t A, result;
  s21_create_matrix(4, -7, &A);
  int check = s21_transpose(&A, &result);
  int check_origin = 1;
  ck_assert_uint_eq(check, check_origin);
  s21_remove_matrix(&result);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_sum_1) {
  matrix_t A, B, C, C1;
  ck_assert(OK == s21_create_matrix(2, 2, &A));
  ck_assert(OK == s21_create_matrix(2, 2, &B));
  ck_assert(OK == s21_create_matrix(2, 2, &C));

  A.matrix[0][0] = 1.6;
  A.matrix[0][1] = -1;
  A.matrix[1][0] = -67;
  A.matrix[1][1] = 23;

  B.matrix[0][0] = 1.6;
  B.matrix[0][1] = -1;
  B.matrix[1][0] = -67;
  B.matrix[1][1] = 23;

  C.matrix[0][0] = 3.2;
  C.matrix[0][1] = -2;
  C.matrix[1][0] = -134;
  C.matrix[1][1] = 46;

  ck_assert(OK == s21_sum_matrix(&A, &B, &C1));
  ck_assert(SUCCESS == s21_eq_matrix(&C, &C1));
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
  s21_remove_matrix(&C1);
}
END_TEST

START_TEST(test_sum_2) {
  matrix_t A, B, C;
  ck_assert(OK == s21_create_matrix(2, 2, &A));
  ck_assert(OK == s21_create_matrix(2, 3, &B));

  ck_assert(CALCULATION_ERROR == s21_sum_matrix(&A, &B, &C));
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST

START_TEST(test_sub_1) {
  matrix_t A, B, C, C1;
  ck_assert(OK == s21_create_matrix(2, 2, &A));
  ck_assert(OK == s21_create_matrix(2, 2, &B));
  ck_assert(OK == s21_create_matrix(2, 2, &C));

  A.matrix[0][0] = 1.6;
  A.matrix[0][1] = -1;
  A.matrix[1][0] = -67;
  A.matrix[1][1] = 23;

  B.matrix[0][0] = 1.6;
  B.matrix[0][1] = -1;
  B.matrix[1][0] = -67;
  B.matrix[1][1] = 23;

  ck_assert(OK == s21_sub_matrix(&A, &B, &C1));
  ck_assert(SUCCESS == s21_eq_matrix(&C, &C1));
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
  s21_remove_matrix(&C1);
}
END_TEST

START_TEST(test_sub_2) {
  matrix_t A, B, C;
  ck_assert(OK == s21_create_matrix(2, 2, &A));
  ck_assert(OK == s21_create_matrix(2, 3, &B));

  ck_assert(CALCULATION_ERROR == s21_sub_matrix(&A, &B, &C));
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST

Suite *s21_matrix_suite(void) {
  Suite *s = suite_create("s21_matrix");
  TCase *tc = tcase_create("tests");

  tcase_add_test(tc, test_create_remove_1);
  tcase_add_test(tc, test_create_remove_2);

  tcase_add_test(tc, test_determinant_1);
  tcase_add_test(tc, test_determinant_2);
  tcase_add_test(tc, test_determinant_3);
  tcase_add_test(tc, test_determinant_4);
  tcase_add_test(tc, test_determinant_5);
  tcase_add_test(tc, test_determinant_6);
  tcase_add_test(tc, test_determinant_7);
  tcase_add_test(tc, test_determinant_8);
  tcase_add_test(tc, test_determinant_9);

  tcase_add_test(tc, test_calc_complements_1);
  tcase_add_test(tc, test_calc_complements_2);
  tcase_add_test(tc, test_calc_complements_3);
  tcase_add_test(tc, test_calc_complements_4);

  tcase_add_test(tc, test_inverse_1);
  tcase_add_test(tc, test_inverse_2);
  tcase_add_test(tc, test_inverse_3);
  tcase_add_test(tc, test_inverse_4);
  tcase_add_test(tc, test_inverse_5);
  tcase_add_test(tc, test_inverse_6);

  tcase_add_test(tc, test_eq_matrix);

  tcase_add_test(tc, test_mult_number_1);
  tcase_add_test(tc, test_mult_matrix_1);
  tcase_add_test(tc, test_mult_matrix_2);

  tcase_add_test(tc, test_transpose_1);
  tcase_add_test(tc, test_transpose_2);

  tcase_add_test(tc, test_sum_1);
  tcase_add_test(tc, test_sum_2);

  tcase_add_test(tc, test_sub_1);
  tcase_add_test(tc, test_sub_2);

  suite_add_tcase(s, tc);
  return s;
}

int main() {
  int number_failed;
  Suite *s = s21_matrix_suite();
  SRunner *sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
