#include "s21_matrix_oop.h"

void S21Matrix::create_new_matrix(int rows, int cols) {
  _rows = rows;
  _cols = cols;
  _matrix = new double *[rows];
  for (int i = 0; i < _rows; i++) _matrix[i] = new double[_cols]{};
}

S21Matrix::S21Matrix() { create_new_matrix(3, 3); }

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows < 1 || cols < 1)
    throw std::out_of_range(
        "The number of rows and cols must be greater than 0");
  create_new_matrix(rows, cols);
}

S21Matrix::~S21Matrix() {
  if (_matrix) {
    for (int i = 0; i < _rows; i++) delete[](_matrix[i]);
    delete[](_matrix);
  }
  _rows = 0;
  _cols = 0;
}

S21Matrix::S21Matrix(const S21Matrix &other) {
  if (this != &other) {
    create_new_matrix(other._rows, other._cols);
    *this = other;
  }
}

S21Matrix::S21Matrix(S21Matrix &&other) : _rows(0), _cols(0), _matrix(nullptr) {
  std::swap(_rows, other._rows);
  std::swap(_cols, other._cols);
  std::swap(_matrix, other._matrix);
}

bool is_equal(double x, double y) {
  return std::fabs(x - y) < 0.0000000000001;
}

bool S21Matrix::eq_matrix(const S21Matrix &other) {
  bool result = true;
  if (_rows == other._rows && _cols == other._cols) {
    for (int i = 0; i < _rows && result; i++)
      for (int j = 0; j < _cols && result; j++)
        if (_matrix[i][j] != other._matrix[i][j]) {
          result = is_equal(_matrix[i][j], other._matrix[i][j]);
        }
  } else {
    result = false;
  }
  return result;
}

void S21Matrix::sum_matrix(const S21Matrix &other) {
  if (_rows != other._rows || _cols != other._cols)
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _cols; j++) _matrix[i][j] += other._matrix[i][j];
}

void S21Matrix::sub_matrix(const S21Matrix &other) {
  if (_rows != other._rows || _cols != other._cols)
    throw std::out_of_range("Matrices should have the same size");
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _cols; j++) _matrix[i][j] -= other._matrix[i][j];
}

void S21Matrix::mul_number(const double num) {
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _cols; j++) _matrix[i][j] *= num;
}

void S21Matrix::mul_matrix(const S21Matrix &other) {
  if (_cols != other._rows) {
    throw std::out_of_range(
        "Differrent number of 1st matrix cols and 2nd matrix rows");
  }
  S21Matrix res(_rows, other._cols);
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < other._cols; j++)
      for (int k = 0; k < _cols && k < other._rows; k++)
        res._matrix[i][j] += _matrix[i][k] * other._matrix[k][j];
  *this = res;
}

S21Matrix S21Matrix::transpose() {
  S21Matrix transponed(_rows, _cols);
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _cols; j++) transponed._matrix[j][i] = _matrix[i][j];
  return transponed;
}

S21Matrix S21Matrix::calc_complements() {
  if (_cols != _rows) throw std::out_of_range("Matrix must be square");
  S21Matrix res(_cols, _rows);
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      S21Matrix buf = minor_matrix(i, j);
      res._matrix[i][j] = pow(-1, i + j) * buf.determinant();
    }
  }
  res.minus_zero();
  return res;
}

S21Matrix S21Matrix::minor_matrix(int row_num, int col_num) {
  S21Matrix res(_rows - 1, _rows - 1);
  int shift_row = 0;
  for (int i = 0; i < _rows - 1; i++) {
    if (i == row_num) shift_row = 1;
    int shift_col = 0;
    for (int j = 0; j < _rows - 1; j++) {
      if (j == col_num) shift_col = 1;
      res._matrix[i][j] = _matrix[i + shift_row][j + shift_col];
    }
  }
  res.minus_zero();
  return res;
}

double S21Matrix::determinant() {
  if (_rows != _cols) return NAN;
  double result = 0;
  int degree = 1;
  if (_rows == 1) {
    result = _matrix[0][0];
  } else if (_rows == 2) {
    result = _matrix[0][0] * _matrix[1][1] - _matrix[0][1] * _matrix[1][0];
  } else {
    for (int j = 0; j < _rows; j++) {
      S21Matrix less = minor_matrix(0, j);
      result += pow(-1.0, j) * _matrix[0][j] * less.determinant();
      degree = -degree;
    }
  }
  return result;
}

void S21Matrix::fill() {
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _cols; j++) _matrix[i][j] = (i + j);
}

S21Matrix S21Matrix::inverse_matrix() {
  double det = determinant();
  if (std::isnan(det) || det == 0)
    throw std::out_of_range("Can't inverse, determinant is zero or nan");
  S21Matrix buf = calc_complements();
  S21Matrix res = buf.transpose();
  for (int i = 0; i < res._rows; i++)
    for (int j = 0; j < res._cols; j++)
      res._matrix[i][j] /= det * pow(1, i + j);
  res.minus_zero();
  return res;
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (_rows != other._rows || _cols != other._cols) {
    this->~S21Matrix();
    create_new_matrix(other._rows, other._cols);
  }
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _cols; j++)
    _matrix[i][j] = other._matrix[i][j];
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix result(*this);
  result.sum_matrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix result(*this);
  result.mul_matrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix result(*this);
  result.sub_matrix(other);
  return result;
}

double S21Matrix::operator()(int row_num, int col_num) {
  if (row_num >= _rows || col_num >= _cols)
    throw std::out_of_range("Out of matrix");
  return _matrix[row_num][col_num];
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  sum_matrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  sub_matrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  mul_matrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const double number) {
  mul_number(number);
  return *this;
}

bool S21Matrix::operator==(const S21Matrix &other) { return eq_matrix(other); }

S21Matrix operator*(const S21Matrix &matr, const double number) {
  S21Matrix result(matr);
  result.mul_number(number);
  return result;
}

S21Matrix operator*(const double number, const S21Matrix &matr) {
  return matr * number;
}

void S21Matrix::minus_zero() {
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _cols; j++)
      if (_matrix[i][j] > -0.0000001 && _matrix[i][j] < 0)
        _matrix[i][j] = 0;
}
