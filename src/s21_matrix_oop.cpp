#include "s21_matrix_oop.h"

#include <gtest/gtest.h>

void S21Matrix::NewMatrix() {
  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; ++i) {
    matrix_[i] = new double[cols_]();
  }
}

void S21Matrix::DeleteMatrix() {
  for (int i = 0; i < rows_; ++i) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
  matrix_ = NULL;
}

S21Matrix::S21Matrix() {
  rows_ = 0;
  cols_ = 0;
  matrix_ = NULL;
}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ < 1 || cols_ < 1) {
    throw std::invalid_argument("Invalid matrix size");
  }
  NewMatrix();
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
  NewMatrix();

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix &&other) {
  rows_ = other.rows_;
  cols_ = other.cols_;
  NewMatrix();
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }

  for (int i = 0; i < rows_; ++i) {
    if (other.matrix_[i] != NULL) delete[] other.matrix_[i];
  }
  if (other.matrix_ != NULL) delete[] other.matrix_;
  other.matrix_ = NULL;
  other.rows_ = 0;
  other.cols_ = 0;
}

S21Matrix::~S21Matrix() {
  if (matrix_ != NULL) {
    DeleteMatrix();
  }
}

bool S21Matrix::EqMatrix(const S21Matrix &other) const {
  bool e = true;
  if (rows_ != other.rows_ || cols_ != other.cols_)
    e = false;

  else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (matrix_[i][j] != other.matrix_[i][j]) e = false;
      }
    }
  }
  return e;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    if (matrix_ == NULL || other.matrix_ == NULL || rows_ <= 0 ||
        other.rows_ <= 0 || cols_ <= 0 || other.cols_ <= 0) {
      throw std::out_of_range(
          "Incorrect input, matrices should have the same size");
    }
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    if (matrix_ == NULL || other.matrix_ == NULL || rows_ <= 0 ||
        other.rows_ <= 0 || cols_ <= 0 || other.cols_ <= 0) {
      throw std::out_of_range(
          "Incorrect input, matrices should have the same size");
    }
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  if (matrix_ == NULL || rows_ <= 0 || cols_ <= 0) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] * num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    if (matrix_ == NULL || other.matrix_ == NULL || rows_ <= 0 ||
        other.rows_ <= 0 || cols_ <= 0 || other.cols_ <= 0) {
      if (rows_ != other.cols_) {
        throw std::out_of_range(
            "Incorrect input, matrices should have the same size");
      }
    }
  }

  S21Matrix result(rows_, cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++)
      for (int add = 0; add < cols_; add++) {
        result.matrix_[i][j] += (matrix_[i][add] * other.matrix_[add][j]);
      }
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = result.matrix_[i][j];
    }
  }
}

S21Matrix S21Matrix::Transpose() {
  if (matrix_ == NULL || rows_ <= 0 || cols_ <= 0) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }
  S21Matrix result = S21Matrix(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      result.matrix_[i][j] = matrix_[j][i];
    }
  }
  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  if (matrix_ == NULL || rows_ <= 0 || cols_ <= 0) {
    if (rows_ != cols_) {
      throw std::out_of_range(
          "Incorrect input, matrices should have the same size");
    }
  }

  S21Matrix result(rows_, cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      S21Matrix minor = s21_minor(i + 1, j + 1);

      result.matrix_[i][j] = std::pow(-1, i + j) * minor.Determinant();
    }
  }
  return result;
}

double S21Matrix::Determinant() {
  if (matrix_ == NULL || rows_ <= 0 || cols_ <= 0) {
    if (rows_ != cols_) {
      throw std::out_of_range(
          "Incorrect input, matrices should have the same size");
    }
  }

  double result = s21_recursive_det();
  if (fabs(result) < 1e-6) result = 0;
  return result;
}

double S21Matrix::s21_recursive_det() {
  double res = 0;
  if (rows_ == 1) {
    res = matrix_[0][0];
  }
  if (rows_ == 2) {
    res = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  }

  if (rows_ > 2) {
    for (int i = 0; i < rows_; i++) {
      S21Matrix minor = s21_minor(1, i + 1);

      res += std::pow(-1, i) * matrix_[0][i] * minor.s21_recursive_det();

      minor.DeleteMatrix();
    }
  }
  return res;
}

S21Matrix S21Matrix::s21_minor(int delete_rows, int delete_columns) {
  if (matrix_ == NULL)
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");

  int minor_rows, minor_columns;

  S21Matrix result = S21Matrix((rows_ - 1), (cols_ - 1));

  for (int i = 0; i < rows_; i++) {
    minor_rows = i;

    if (i > delete_rows - 1) minor_rows--;

    for (int j = 0; j < cols_; j++) {
      minor_columns = j;

      if (j > delete_columns - 1) minor_columns--;

      {
        if (i != delete_rows - 1 && j != delete_columns - 1) {
          result.matrix_[minor_rows][minor_columns] = matrix_[i][j];
        }
      }
    }
  }

  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  if (matrix_ == NULL || rows_ <= 0 || cols_ <= 0) {
    if (rows_ != cols_) {
      throw std::out_of_range(
          "Incorrect input, matrices should have the same size");
    }
  }
  double determ = S21Matrix::Determinant();
  if (determ == 0) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }

  S21Matrix data = S21Matrix::CalcComplements();
  S21Matrix datatranspose = data.Transpose();
  S21Matrix result = datatranspose * (1.0 / determ);
  return result;
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  if (cols_ != other.cols_ || rows_ != other.rows_) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  // creating result matrix
  S21Matrix res(rows_, cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      res.matrix_[i][j] = matrix_[i][j];
    }
  }
  res.SumMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  if (cols_ != other.cols_ || rows_ != other.rows_) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  // creating result matrix
  S21Matrix res(rows_, cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      res.matrix_[i][j] = matrix_[i][j];
    }
  }
  res.SubMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  // creating result matrix
  S21Matrix res(rows_, cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      res.matrix_[i][j] = matrix_[i][j];
    }
  }

  res.MulMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix res(rows_, cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      res.matrix_[i][j] = matrix_[i][j];
    }
  }
  res.MulNumber(num);
  return res;
}

bool S21Matrix::operator==(const S21Matrix &other) const {
  S21Matrix res(rows_, cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      res.matrix_[i][j] = matrix_[i][j];
    }
  }
  bool value = res.EqMatrix(other);

  return value;
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (this != &other) {
    DeleteMatrix();

    rows_ = other.rows_;
    cols_ = other.cols_;
    NewMatrix();

    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
  return *this;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  *this = *this + other;
  return *this;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  *this = *this - other;
  return *this;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  *this = *this * other;
  return *this;
}

S21Matrix &S21Matrix::operator*=(double num) {
  *this = *this * num;
  return *this;
}

// index operator overload
double &S21Matrix::operator()(int i, int j) {
  if (i < 0 || i >= rows_ || j < 0 || j >= cols_) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  return matrix_[i][j];
}

// публичные методы доступа (accessor)

int S21Matrix::GetRows() const { return rows_; }

int S21Matrix::GetCols() const { return cols_; }

// публичные методы изменения (mutator)

void S21Matrix::SetRows(const int rows) { this->rows_ = rows; }

void S21Matrix::SetCols(int cols) { this->cols_ = cols; }

void S21Matrix::s21_print_matrix() {
  std::cout << "-------------------------------\n";
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      std::cout << " " << matrix_[i][j];
    }
    std::cout << '\n';
  }
  std::cout << "-------------------------------\n";
}
