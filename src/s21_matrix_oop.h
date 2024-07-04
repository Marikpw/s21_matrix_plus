#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H
#include <algorithm>
#include <cmath>

class S21Matrix {
 private:
  int rows_, cols_;
  double **matrix_;

 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other);
  ~S21Matrix();
  int getRows();
  int getCols();
  bool EqMatrix(const S21Matrix &other);
  double &operator()(int rows, int cols);
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  S21Matrix Transpose();
  bool operator==(const S21Matrix &other);
  S21Matrix operator+(const S21Matrix &other);
  S21Matrix operator-(const S21Matrix &other);
  S21Matrix &operator=(S21Matrix &&other);
  S21Matrix &operator=(const S21Matrix &other);
  S21Matrix &operator+=(const S21Matrix &other);
  S21Matrix &operator-=(const S21Matrix &other);
  void MulMatrix(const S21Matrix &other);
  S21Matrix operator*(const S21Matrix &other);
  S21Matrix &operator*=(const S21Matrix &other);
  S21Matrix &operator*=(const double other);
  friend S21Matrix operator*(double, S21Matrix &);
  friend S21Matrix operator*(S21Matrix &matrix, double num);
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();
  S21Matrix Minor(int row, int col);
  void setRows(int newRows);
  void setCols(int newCols);
  void createMatrix();
  void free();
};

#endif