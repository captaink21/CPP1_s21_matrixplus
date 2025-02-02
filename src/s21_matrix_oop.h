
#ifndef __S21MATRIX_H__
#define __S21MATRIX_H__

#include <cmath>
#include <cstring>
#include <iostream>

class S21Matrix {
 private:
  // attributes
  int rows_, cols_;  // rows and columns attributes
  double** matrix_;  // pointer to the memory where the matrix will be allocated

 public:
  S21Matrix();                        // default constructor
  S21Matrix(int rows, int cols);      // parameterized constructor
  S21Matrix(const S21Matrix& other);  // copy cnstructor
  S21Matrix(S21Matrix&& other);       // move cnstructor
  ~S21Matrix();                       // destructor

  // public methods

  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  // some operators overloads
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(const double num);
  bool operator==(const S21Matrix& other) const;
  S21Matrix& operator=(const S21Matrix& other);  // assignment operator overload
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(double num);
  double& operator()(int row, int col);  // index operator overload

  // публичные методы доступа (accessor)
  int GetRows() const;
  int GetCols() const;

  // публичные методы изменения (mutator)

  void SetRows(const int rows);
  void SetCols(int cols);

  // Вспомогательные

  void NewMatrix();
  void DeleteMatrix();
  double s21_recursive_det();
  void s21_print_matrix();
  S21Matrix s21_minor(int delete_rows, int delete_columns);
};

#endif
