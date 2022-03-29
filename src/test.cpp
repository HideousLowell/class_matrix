#include "gtest/gtest.h"
#include "s21_matrix_oop.h"

void S21Matrix::set_number(double number, int row, int col) {
  _matrix[row][col] = number;
}

TEST(test_sum_matrix, 1) {
  S21Matrix matr1(3, 3);
  matr1.fill();
  S21Matrix matr2(matr1);
  S21Matrix matr3 = matr2;
  matr1.sum_matrix(matr2);
  matr2 += matr2;
  matr3 = matr3 + matr3;
  S21Matrix res(3, 3);
  res.set_number(0, 0, 0);
  res.set_number(2., 0, 1);
  res.set_number(4., 0, 2);
  res.set_number(2., 1, 0);
  res.set_number(4., 1, 1);
  res.set_number(6., 1, 2);
  res.set_number(4., 2, 0);
  res.set_number(6., 2, 1);
  res.set_number(8., 2, 2);
  ASSERT_TRUE(matr1 == res && matr2 == res && matr3 == res);
}

TEST(test_sub_matrix, 1) {
  S21Matrix matr1(4, 4);
  matr1.fill();
  S21Matrix matr2(matr1);
  S21Matrix matr3(4, 4);
  S21Matrix res(matr1);
  matr2 += matr2;
  matr3 = matr2 - matr1;
  matr2 -= matr1;
  ASSERT_TRUE(matr2 == res && matr3 == res);
}

TEST(test_mul_matrix, 1) {
  S21Matrix matr1(4, 5);
  matr1.fill();
  S21Matrix matr2(5, 4);
  matr2.fill();
  S21Matrix matr3 = matr1 * matr2;
  matr1 *= matr2;
  S21Matrix res(4, 4);
  res.set_number(30, 0, 0);
  res.set_number(40., 0, 1);
  res.set_number(50, 0, 2);
  res.set_number(60, 0, 3);
  res.set_number(40, 1, 0);
  res.set_number(55, 1, 1);
  res.set_number(70, 1, 2);
  res.set_number(85, 1, 3);
  res.set_number(50, 2, 0);
  res.set_number(70, 2, 1);
  res.set_number(90, 2, 2);
  res.set_number(110, 2, 3);
  res.set_number(60, 3, 0);
  res.set_number(85, 3, 1);
  res.set_number(110, 3, 2);
  res.set_number(135, 3, 3);
  ASSERT_TRUE(matr1 == res && matr3 == res);
}

TEST(test_transpose_matrix, 1) {
  S21Matrix matr1;
  matr1.fill();
  matr1.set_number(-3, 0, 1);
  matr1.set_number(50, 0, 2);
  matr1.set_number(99, 1, 0);
  matr1.set_number(88, 2, 0);
  S21Matrix matr2 = matr1.transpose();
  S21Matrix res(matr1);
  res.set_number(-3, 1, 0);
  res.set_number(50, 2, 0);
  res.set_number(99, 0, 1);
  res.set_number(88, 0, 2);
  ASSERT_TRUE(matr2 == res);
}

TEST(test_determinant_matrix, 1) {
  S21Matrix matr1(5, 5);
  matr1.fill();
  matr1.set_number(-3, 0, 1);
  matr1.set_number(5.7, 0, 2);
  matr1.set_number(9.9, 1, 0);
  matr1.set_number(8.8, 2, 0);
  matr1.set_number(4, 3, 3);
  matr1.set_number(2.4, 3, 3);
  double res = -178.776;
  double func_res = matr1.determinant();
  ASSERT_NEAR(res, func_res, 0.0001);
}

TEST(test_braces_matrix, 1) {
  S21Matrix matr1(5, 5);
  matr1.set_number(-3, 0, 1);
  matr1.set_number(5.7, 0, 2);
  matr1.set_number(9.9, 1, 0);
  double func_res1 = matr1(0, 1);
  double res1 = -3.;
  double func_res2 = matr1(0, 2);
  double res2 = 5.7;
  double func_res3 = matr1(1, 0);
  double res3 = 9.9;
  ASSERT_DOUBLE_EQ(res1, func_res1);
  ASSERT_DOUBLE_EQ(res2, func_res2);
  ASSERT_DOUBLE_EQ(res3, func_res3);
}

TEST(test_inverse_matrix, 3) {
  S21Matrix matr1(5, 5);
  matr1.fill();
  matr1.set_number(5.7, 0, 2);
  matr1.set_number(9.9, 1, 0);
  matr1.set_number(8.8, 2, 0);
  matr1.set_number(4, 3, 3);
  matr1.set_number(2.4, 3, 3);
  S21Matrix matr2 = matr1.inverse_matrix();
  S21Matrix res(5, 5);
  res.set_number(1, 0, 0);
  res.set_number(1, 1, 1);
  res.set_number(1, 2, 2);
  res.set_number(1, 3, 3);
  res.set_number(1, 4, 4);
  matr2 *= matr1;
  ASSERT_TRUE(matr2 == res);
}

TEST(test_mul_number_matrix, 3) {
  S21Matrix matr1(2, 3);
  matr1.fill();
  S21Matrix matr2 = 3 * matr1;
  matr2 = matr2 * -1;
  matr1 *= -3;
  S21Matrix res(2, 3);
  res.set_number(-0.0, 0, 0);
  res.set_number(-3.0, 0, 1);
  res.set_number(-6.0, 0, 2);
  res.set_number(-3.0, 1, 0);
  res.set_number(-6.0, 1, 1);
  res.set_number(-9.0, 1, 2);
  ASSERT_TRUE(matr1 == res && matr2 == res);
}

TEST(test_move_matrix, 3) {
  S21Matrix matr1(5, 3);
  matr1.fill();
  S21Matrix matr2 = std::move(matr1);
  int row1 = matr1.get_rows();
  int col1 = matr1.get_cols();
  int row2 = matr2.get_rows();
  int col2 = matr2.get_cols();
  ASSERT_EQ(row1, 0);
  ASSERT_EQ(col1, 0);
  ASSERT_EQ(row2, 5);
  ASSERT_EQ(col2, 3);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
