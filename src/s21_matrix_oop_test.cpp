

#include "s21_matrix_oop.h"

#include <gtest/gtest.h>

#include <cstdio>
#include <iostream>

void print(S21Matrix &test) {
  std::cout << "-------------------------------\n";
  for (int i = 0; i < test.GetRows(); i++) {
    for (int j = 0; j < test.GetCols(); j++) {
      std::cout << test(i, j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "-------------------------------\n";
}

namespace {
TEST(Constructor, Default) {
  auto test = S21Matrix();
  EXPECT_EQ(test.GetRows(), 0);
  EXPECT_EQ(test.GetCols(), 0);
}

TEST(Constructor, arg2) {
  auto test = S21Matrix(7, 7);
  EXPECT_EQ(test.GetRows(), 7);
  EXPECT_EQ(test.GetCols(), 7);
}

TEST(Constructor, copy) {
  auto test = S21Matrix(7, 7);
  test(0, 0) = 7;
  auto data = S21Matrix(test);
  EXPECT_EQ(test, data);
}

TEST(Constructor, move) {
  auto test1 = S21Matrix(7, 7);
  test1(1, 1) = 7;
  auto test2 = S21Matrix(test1);
  auto test3 = S21Matrix(std::move(test1));
  EXPECT_EQ(test2, test3);
}

TEST(Operation, EqMatrix) {
  auto test1 = S21Matrix(6, 6);

  test1(3, 3) = 7;
  auto test2 = test1;
  EXPECT_TRUE(test1.EqMatrix(test2));

  EXPECT_TRUE(test1.EqMatrix(test2) == (test1 == test2));
  test2(4, 4) = 7;
  EXPECT_FALSE(test1.EqMatrix(test2));
  EXPECT_TRUE(test1.EqMatrix(test2) == (test1 == test2));
  auto test3 = S21Matrix(5, 5);

  EXPECT_FALSE(test1.EqMatrix(test3));

  EXPECT_TRUE(test1.EqMatrix(test3) == (test1 == test3));
}

TEST(Operation, SumMatrix) {
  auto test1 = S21Matrix(7, 7);
  test1(3, 3) = 7;
  auto test2 = test1;
  test1.SumMatrix(test2);
  EXPECT_EQ(test1(3, 3), 14);
}

TEST(Operation, SubMatrix) {
  auto test1 = S21Matrix(7, 7);
  test1(3, 3) = 7;
  auto test2 = test1;
  test2(3, 3) = 2;
  test1.SubMatrix(test2);
  EXPECT_EQ(test1(3, 3), 5);
}

TEST(Operation, MulNumber) {
  auto test1 = S21Matrix(7, 7);
  test1(3, 3) = 3;
  test1.MulNumber(2);
  EXPECT_EQ(test1(3, 3), 6);
}

TEST(Operation, MulMatrix) {
  auto test1 = S21Matrix(3, 3);
  for (int i = 0; i < test1.GetRows(); i++) {
    for (int j = 0; j < test1.GetCols(); j++) {
      test1(i, j) = i * 3 + j;
    }
  }

  auto test2 = test1;
  test1.MulMatrix(test2);
  test2(0, 0) = 15;
  test2(0, 1) = 18;
  test2(0, 2) = 21;
  test2(1, 0) = 42;
  test2(1, 1) = 54;
  test2(1, 2) = 66;
  test2(2, 0) = 69;
  test2(2, 1) = 90;
  test2(2, 2) = 111;

  EXPECT_EQ(test1, test2);
}

TEST(LinearOperations, Transpose) {
  auto test1 = S21Matrix(3, 3);
  for (int i = 0; i < test1.GetRows(); i++) {
    for (int j = 0; j < test1.GetCols(); j++) {
      test1(i, j) = i * 3 + j;
    }
  }
  test1 = test1.Transpose();
  auto test2 = S21Matrix(3, 3);
  test2(0, 0) = 0;
  test2(0, 1) = 3;
  test2(0, 2) = 6;
  test2(1, 0) = 1;
  test2(1, 1) = 4;
  test2(1, 2) = 7;
  test2(2, 0) = 2;
  test2(2, 1) = 5;
  test2(2, 2) = 8;
  EXPECT_EQ(test1, test2);
}

TEST(LinearOperations, CalcComplements) {
  auto test1 = S21Matrix(3, 3);
  test1(0, 0) = 1;
  test1(0, 1) = 2;
  test1(0, 2) = 3;
  test1(1, 0) = 4;
  test1(1, 1) = 5;
  test1(1, 2) = 6;
  test1(2, 0) = 7;
  test1(2, 1) = 8;
  test1(2, 2) = 10;

  test1 = test1.CalcComplements();

  auto test2 = S21Matrix(3, 3);
  test2(0, 0) = 2;
  test2(0, 1) = 2;
  test2(0, 2) = -3;
  test2(1, 0) = 4;
  test2(1, 1) = -11;
  test2(1, 2) = 6;
  test2(2, 0) = -3;
  test2(2, 1) = 6;
  test2(2, 2) = -3;
  EXPECT_EQ(test1, test2);
}

TEST(LinearOperations, Determinant) {
  auto test1 = S21Matrix(3, 3);

  test1(0, 0) = 7;
  test1(0, 1) = 2;
  test1(0, 2) = 3;
  test1(1, 0) = 4;
  test1(1, 1) = 7;
  test1(1, 2) = 6;
  test1(2, 0) = 9;
  test1(2, 1) = 8;
  test1(2, 2) = 7;

  EXPECT_DOUBLE_EQ(test1.Determinant(), -34);
}

TEST(LinearOperations, InverseMatrix) {
  auto test1 = S21Matrix(3, 3);
  test1(0, 0) = 2;
  test1(0, 1) = 5;
  test1(0, 2) = 7;
  test1(1, 0) = 6;
  test1(1, 1) = 3;
  test1(1, 2) = 4;
  test1(2, 0) = 5;
  test1(2, 1) = -2;
  test1(2, 2) = -3;
  test1 = test1.InverseMatrix();
  auto test2 = S21Matrix(3, 3);
  test2(0, 0) = 1;
  test2(0, 1) = -1;
  test2(0, 2) = 1;
  test2(1, 0) = -38;
  test2(1, 1) = 41;
  test2(1, 2) = -34;
  test2(2, 0) = 27;
  test2(2, 1) = -29;
  test2(2, 2) = 24;
  EXPECT_EQ(test1, test2);
}

TEST(accessor, GetRows) {
  auto test1 = S21Matrix();
  EXPECT_EQ(test1.GetRows(), 0);
  auto test2 = S21Matrix(7, 7);
  EXPECT_EQ(test2.GetRows(), 7);
}

TEST(accessor, GetCols) {
  auto test1 = S21Matrix();
  EXPECT_EQ(test1.GetCols(), 0);
  auto test2 = S21Matrix(7, 7);
  EXPECT_EQ(test2.GetCols(), 7);
}

TEST(mutator, SetRows) {
  auto test1 = S21Matrix(7, 7);
  test1.SetRows(4);
  EXPECT_EQ(test1.GetRows(), 4);
}

TEST(mutator, SetCols) {
  auto test1 = S21Matrix(7, 7);
  test1.SetCols(4);
  EXPECT_EQ(test1.GetCols(), 4);
}

}  // namespace
