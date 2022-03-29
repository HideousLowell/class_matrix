#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_

#include <cmath>
#include <iostream>

class S21Matrix {
  int _rows, _cols;
  double** _matrix;

  // overload multiply for numbers
  friend S21Matrix operator*(const double number, const S21Matrix& matr);
  friend S21Matrix operator*(const S21Matrix& matr, const double number);

 protected:
  void create_new_matrix(int rows, int cols);

 public:
  // constructors and destruktor
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();

  // accessors and mutators
  int get_rows() const { return _rows; }
  int get_cols() const { return _cols; }
  void set_rows(int number) { _rows = number; }
  void set_cols(int number) { _cols = number; }

  // overloaded operators
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const double num);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  double operator()(int row_num, int col_num);
  bool operator==(const S21Matrix& other);

  // mathematical operations with matrices
  bool eq_matrix(const S21Matrix& other);
  void sum_matrix(const S21Matrix& other);
  void sub_matrix(const S21Matrix& other);
  void mul_number(const double num);
  void mul_matrix(const S21Matrix& other);
  S21Matrix transpose();
  S21Matrix calc_complements();
  double determinant();
  S21Matrix inverse_matrix();
  S21Matrix minor_matrix(int row_num, int col_num);

// for tests
  void fill();
  void set_number(double number, int row, int col);
  void minus_zero();
};

#endif  // SRC_S21_MATRIX_OOP_H_
