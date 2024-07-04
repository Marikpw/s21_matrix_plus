#include <iostream>

#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr) {}

void S21Matrix::createMatrix() {
  matrix_ = new double *[rows_];

  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = 0.0;
    }
  }
}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ <= 0 || cols_ <= 0)
    throw std::length_error("Rows and cols must be > 0");

  createMatrix();
}

S21Matrix::~S21Matrix() {
  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
    matrix_[i] = nullptr;
  }

  delete[] matrix_;
  matrix_ = nullptr;
}

int S21Matrix::getRows() { return rows_; }

int S21Matrix::getCols() { return cols_; }

double &S21Matrix::operator()(int rows, int cols) {
  if (rows >= rows_ || cols >= cols_ || rows < 0 || cols < 0) {
    throw std::out_of_range("Index is out of range");
  }

  return matrix_[rows][cols];
}

bool S21Matrix::EqMatrix(const S21Matrix &other) {
  bool flag = true;

  if (rows_ == other.rows_ && cols_ == other.cols_) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-07) {
          flag = false;
        }
      }
    }
  } else {
    flag = false;
  }
  return flag;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
    throw std::out_of_range("Matrices are not the same size");
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] += other.matrix_[i][j];
      }
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_) {
    throw std::out_of_range("Matrices are not the same size");
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] -= other.matrix_[i][j];
      }
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result(cols_, rows_);
  for (int i = 0; i < cols_; i++) {
    for (int j = 0; j < rows_; j++) {
      result.matrix_[i][j] = matrix_[j][i];
    }
  }
  return result;
}

S21Matrix::S21Matrix(const S21Matrix &other) {
  rows_ = other.rows_;
  cols_ = other.cols_;

  if (rows_ <= 0 || cols_ <= 0)
    throw std::length_error("Rows and cols must be > 0");

  createMatrix();

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix &&other) {
  matrix_ = other.matrix_;
  rows_ = other.rows_;
  cols_ = other.cols_;

  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

bool S21Matrix::operator==(const S21Matrix &other) { return EqMatrix(other); }

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

S21Matrix &S21Matrix::operator=(S21Matrix &&other) {
  if (&other != this) {
    free();

    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = other.matrix_;

    other.rows_ = 0;
    other.cols_ = 0;
    other.matrix_ = nullptr;
  }
  return *this;
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (&other != this) {
    free();

    rows_ = other.rows_;
    cols_ = other.cols_;

    createMatrix();

    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }

  return *this;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  SumMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  SubMatrix(other);
  return *this;
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (rows_ != other.cols_ || cols_ != other.rows_) {
    throw std::out_of_range("Number of rows and cols must be even");
  } else {
    S21Matrix result(rows_, other.cols_);
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < other.cols_; j++) {
        for (int n = 0; n < cols_; n++) {
          result.matrix_[i][j] += matrix_[i][n] * other.matrix_[n][j];
        }
      }
    }
    *this = result;
  }
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  MulMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const double other) {
  MulNumber(other);
  return *this;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw std::invalid_argument("Rows must be equal to cols");
  }

  S21Matrix result(rows_, cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j != cols_; j++) {
      S21Matrix minorMatrix = Minor(i, j);
      result.matrix_[i][j] = pow((-1), i + j) * minorMatrix.Determinant();
    }
  }

  return result;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) {
    throw std::invalid_argument("Rows must be equal to cols");
  } else {
    double result = 0.0;
    if (rows_ == 1) {
      result = matrix_[0][0];
    } else if (rows_ == 2) {
      result = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
    } else {
      for (int i = 0; i < cols_; i++) {
        S21Matrix new_matrix = Minor(0, i);
        result += matrix_[0][i] * pow(-1, i) * new_matrix.Determinant();
      }
    }
    return result;
  }
}

S21Matrix S21Matrix::InverseMatrix() {
  double thisDet = Determinant();
  if (fabs(thisDet) < 1e-07) {
    throw std::invalid_argument("Determinant must be >0");
  }
  S21Matrix result(rows_, cols_);

  if (rows_ == 1) {
    result.matrix_[0][0] = 1.0 / matrix_[0][0];
  } else {
    S21Matrix temp = CalcComplements();
    result = temp.Transpose();
    result.MulNumber(1.0 / thisDet);
  }
  return result;
}

S21Matrix S21Matrix::Minor(int row, int col) {
  S21Matrix result(rows_ - 1, cols_ - 1);

  for (int i = 0, min_i = 0; min_i < result.rows_; min_i++) {
    if (row == i) i++;
    for (int j = 0, min_j = 0; min_j < result.cols_; min_j++) {
      if (col == j) j++;
      result.matrix_[min_i][min_j] = matrix_[i][j];
      j++;
    }
    i++;
  }
  return result;
}

S21Matrix operator*(double num, S21Matrix &matrix) {
  S21Matrix result(matrix);
  result.MulNumber(num);
  return result;
}

S21Matrix operator*(S21Matrix &matrix, double num) {
  S21Matrix result(matrix);
  result.MulNumber(num);
  return result;
}

void S21Matrix::setRows(int newRows) {
  if (newRows <= 0) {
    throw std::out_of_range("Amount of rows must be above 0");
  }

  S21Matrix result(newRows, cols_);

  for (int i = 0; i < newRows; i++) {
    for (int j = 0; j < cols_; j++) {
      if (i < rows_) result.matrix_[i][j] = matrix_[i][j];
    }
  }
  *this = result;
}

void S21Matrix::setCols(int newCols) {
  if (newCols <= 0) {
    throw std::out_of_range("Amount of cols must be above 0");
  }
  S21Matrix result(rows_, newCols);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < newCols; j++) {
      if (j < cols_) result.matrix_[i][j] = matrix_[i][j];
    }
  }
  *this = result;
}

void S21Matrix::free() {
  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
}