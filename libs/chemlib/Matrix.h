#ifndef TNT_CMAT_H
#define TNT_CMAT_H

#include <cstdlib>
#include <cassert>
#include <sstream>
#include <vector>
#include <functional>
#include <algorithm>
#include "nonguilib.h"

namespace TNT {

typedef int Subscript;

template <class T>
class Matrix {


public:

  typedef Subscript size_type;
  typedef         T value_type;
  typedef         T element_type;
  typedef         T*  pointer;
  typedef         T*  iterator;
  typedef         T&  reference;
  typedef const T*  const_iterator;
  typedef const T&  const_reference;

  Subscript lbound() const {
    return 1;
  }

protected:
  Subscript m_;
  Subscript n_;
  Subscript mn_;        // total size
  T* v_;
  T** row_;
  T* vm1_ ;         // these point to the same data, but are 1-based
  T** rowm1_;

  // internal helper function to create the array
  // of row pointers

  void initialize(Subscript M, Subscript N) {
    mn_ = M*N;
    m_ = M;
    n_ = N;

    v_ = new T[mn_];
    row_ = new T*[M];
    rowm1_ = new T*[M];

    assert(v_  != NULL);
    assert(row_  != NULL);
    assert(rowm1_ != NULL);

    T* p = v_;
    vm1_ = v_ - 1;
    for (Subscript i = 0; i < M; i++) {
      row_[i] = p;
      rowm1_[i] = p-1;
      p += N ;

    }

    rowm1_-- ;          // compensate for 1-based offset
  }

  void copy(const T*  v) {
    Subscript N = m_ * n_;
    Subscript i;

#ifdef TNT_UNROLL_LOOPS
    Subscript Nmod4 = N & 3;
    Subscript N4 = N - Nmod4;

    for (i = 0; i < N4; i += 4) {
      v_[i] = v[i];
      v_[i+1] = v[i+1];
      v_[i+2] = v[i+2];
      v_[i+3] = v[i+3];
    }

    for (i = N4; i < N; i++) {
      v_[i] = v[i];
    }
#else

    for (i = 0; i < N; i++) {
      v_[i] = v[i];
    }
#endif
  }

  void set(const T& val) {
    Subscript N = m_ * n_;
    Subscript i;

#ifdef TNT_UNROLL_LOOPS
    Subscript Nmod4 = N & 3;
    Subscript N4 = N - Nmod4;

    for (i = 0; i < N4; i += 4) {
      v_[i] = val;
      v_[i+1] = val;
      v_[i+2] = val;
      v_[i+3] = val;
    }

    for (i = N4; i < N; i++) {
      v_[i] = val;
    }
#else

    for (i = 0; i < N; i++) {
      v_[i] = val;
    }

#endif
  }

  void destroy() {
    /* do nothing, if no memory has been previously allocated */
    if (v_ == NULL) {
      return ;
    }

    /* if we are here, then matrix was previously allocated */
    if (v_ != NULL) {
      delete[] (v_);
    }
    if (row_ != NULL) {
      delete[] (row_);
    }

    /* return rowm1_ back to original value */
    rowm1_++;
    if (rowm1_ != NULL) {
      delete[] (rowm1_);
    }
  }

public:

  operator T**() {
    return row_;
  }

  operator T**() const {
    return row_;
  }

  Subscript size() const {
    return mn_;
  }

  // constructors

  Matrix() : m_(0), n_(0), mn_(0), v_(0), row_(0), vm1_(0), rowm1_(0) {
  };

  Matrix(const Matrix<T>& A) {
    initialize(A.m_, A.n_);
    copy(A.v_);
  }

  Matrix(Subscript M, Subscript N, const T& value = T()) {
    initialize(M, N);
    set(value);
  }

  Matrix(Subscript M, Subscript N, const T* v) {
    initialize(M, N);
    copy(v);
  }

  Matrix(Subscript M, Subscript N, const char* s) {
    initialize(M, N);
    //std::istrstream ins(s);
    std::istringstream ins(s);

    Subscript i, j;

    for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++) {
        ins >> row_[i][j];
      }
    }
  }

  // destructor
  //
  ~Matrix() {
    destroy();
  }

  // reallocating
  //
  Matrix<T>& newsize(Subscript M, Subscript N) {
    if (num_rows() == M && num_cols() == N) {
      return *this;
    }

    destroy();
    initialize(M, N);

    return *this;
  }

  // assignments
  //
  Matrix<T>& operator=(const Matrix<T>& A) {
    if (v_ == A.v_) {
      return *this;
    }

    if (m_ == A.m_  && n_ == A.n_) {        // no need to re-alloc
      copy(A.v_);
    } else {
      destroy();
      initialize(A.m_, A.n_);
      copy(A.v_);
    }

    return *this;
  }

  Matrix<T>& operator=(const T& scalar) {
    set(scalar);
    return *this;
  }

  Subscript dim(Subscript d) const {
#ifdef TNT_BOUNDS_CHECK
    assert(d >= 1);
    assert(d <= 2);
#endif
    return (d == 1) ? m_ : ((d == 2) ? n_ : 0);
  }

  Subscript num_rows() const {
    return m_;
  }

  Subscript num_cols() const {
    return n_;
  }

  inline T* operator[](Subscript i) {
#ifdef TNT_BOUNDS_CHECK
    assert(0 <= i);
    assert(i < m_) ;
#endif
    return row_[i];
  }

  inline const T* operator[](Subscript i) const {
#ifdef TNT_BOUNDS_CHECK
    assert(0 <= i);
    assert(i < m_) ;
#endif
    return row_[i];
  }

  inline reference operator()(Subscript i) {
#ifdef TNT_BOUNDS_CHECK
    assert(1 <= i);
    assert(i <= mn_) ;
#endif
    return vm1_[i];
  }

  inline const_reference operator()(Subscript i) const {
#ifdef TNT_BOUNDS_CHECK
    assert(1 <= i);
    assert(i <= mn_) ;
#endif
    return vm1_[i];
  }

  inline reference operator()(Subscript i, Subscript j) {
#ifdef TNT_BOUNDS_CHECK
    assert(1 <= i);
    assert(i <= m_) ;
    assert(1 <= j);
    assert(j <= n_);
#endif
    return rowm1_[i][j];
  }

  inline const_reference operator()(Subscript i, Subscript j) const {
#ifdef TNT_BOUNDS_CHECK
    assert(1 <= i);
    assert(i <= m_) ;
    assert(1 <= j);
    assert(j <= n_);
#endif
    return rowm1_[i][j];
  }

  void SwapRows(Subscript M1, Subscript M2) {
#ifdef TNT_BOUNDS_CHECK
    assert(M1 < m_);
    assert(M2 < m_);
#endif

    if (M1 == M2) {
      return;
    }

    T* tmp = new T[n_];

    memcpy(tmp, row_[M1], n_ * sizeof(T));
    memcpy(row_[M1], row_[M2], n_ * sizeof(T));
    memcpy(row_[M2], tmp, n_ * sizeof(T));

    delete tmp;
  }

  void SwapColumns(Subscript N1, Subscript N2) {
#ifdef TNT_BOUNDS_CHECK
    assert(N1 < n_);
    assert(N2 < n_);
#endif
    if (N1 == N2) {
      return;
    }

    T tmp;
    int i;
    for (i = 0; i < m_; ++i) {
      tmp = row_[i][N1];
      row_[i][N1] = row_[i][N2];
      row_[i][N2] = tmp;
    }
  }

};      //End of Matrix Class

/****************************  I/O  ********************************/
/*
   template <class T>
   std::ostream& operator<<(std::ostream &s, const Matrix<T> &A)
   {
    Subscript M=A.num_rows();
    Subscript N=A.num_cols();

    s << M << " " << N << "\n";

    for (Subscript i=0; i<M; i++)
    {
        for (Subscript j=0; j<N; j++)
        {
            s << A[i][j] << " ";
        }
        s << "\n";
    }


    return s;
   }

   template <class T>
   std::istream& operator>>(std::istream &s, Matrix<T> &A)
   {

    Subscript M, N;

    s >> M >> N;

    if ( !(M == A.num_rows() && N == A.num_cols() ))
    {
        A.newsize(M,N);
    }


    for (Subscript i=0; i<M; i++)
        for (Subscript j=0; j<N; j++)
        {
            s >>  A[i][j];
        }


    return s;
   }
 */
// *******************[ basic matrix algorithms ]***************************

template <class T>
Matrix<T> operator+(const Matrix<T>& A,
                    const Matrix<T>& B) {
  Subscript M = A.num_rows();
  Subscript N = A.num_cols();

  assert(M == B.num_rows());
  assert(N == B.num_cols());

  Matrix<T> tmp(M, N);
  Subscript i, j;

  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      tmp[i][j] = A[i][j] + B[i][j];
    }
  }

  return tmp;
}

template <class T>
Matrix<T> operator-(const Matrix<T>& A,
                    const Matrix<T>& B) {
  Subscript M = A.num_rows();
  Subscript N = A.num_cols();

  assert(M == B.num_rows());
  assert(N == B.num_cols());

  Matrix<T> tmp(M, N);
  Subscript i, j;

  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      tmp[i][j] = A[i][j] - B[i][j];
    }
  }

  return tmp;
}

template <class T>
Matrix<T> mult_element(const Matrix<T>& A,
                       const Matrix<T>& B) {
  Subscript M = A.num_rows();
  Subscript N = A.num_cols();

  assert(M == B.num_rows());
  assert(N == B.num_cols());

  Matrix<T> tmp(M, N);
  Subscript i, j;

  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      tmp[i][j] = A[i][j] * B[i][j];
    }
  }

  return tmp;
}

template <class T>
Matrix<T> transpose(const Matrix<T>& A) {
  Subscript M = A.num_rows();
  Subscript N = A.num_cols();

  Matrix<T> S(N, M);
  Subscript i, j;

  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      S[j][i] = A[i][j];
    }
  }

  return S;
}

template <class T>
inline Matrix<T> matmult(const Matrix<T>& A,
                         const Matrix<T>& B) {

#ifdef TNT_BOUNDS_CHECK
  assert(A.num_cols() == B.num_rows());
#endif

  Subscript M = A.num_rows();
  Subscript N = A.num_cols();
  Subscript K = B.num_cols();

  Matrix<T> tmp(M, K);
  T sum;

  for (Subscript i = 0; i < M; i++) {
    for (Subscript k = 0; k < K; k++) {
      sum = 0;
      for (Subscript j = 0; j < N; j++) {
        sum = sum +  A[i][j] * B[j][k];
      }

      tmp[i][k] = sum;
    }
  }

  return tmp;
}

template <class T>
inline Matrix<T> operator*(const Matrix<T>& A,
                           const Matrix<T>& B) {
  return matmult(A, B);
}

template <class T>
inline int matmult(Matrix<T>& C, const Matrix<T>& A,
                   const Matrix<T>& B) {

  assert(A.num_cols() == B.num_rows());

  Subscript M = A.num_rows();
  Subscript N = A.num_cols();
  Subscript K = B.num_cols();

  C.newsize(M, K);

  T sum;

  const T* row_i;
  const T* col_k;

  for (Subscript i = 0; i < M; i++) {
    for (Subscript k = 0; k < K; k++) {
      row_i  = &(A[i][0]);
      col_k  = &(B[0][k]);
      sum = 0;
      for (Subscript j = 0; j < N; j++) {
        sum  += *row_i * *col_k;
        row_i++;
        col_k += K;
      }
      C[i][k] = sum;
    }
  }

  return 0;
}

template <class T>
void AugmentMatrix(const Matrix<T>& A, Matrix<T>& B) {
  Subscript i, j, m = A.num_rows(), n = A.num_cols();
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) {
      B[i][j] = A[i][j];
      B[i][j+n] = (i == j) ? 1 : 0;
    }
  }
}

template <class T>
bool InvertMatrix(Matrix<T>& A, double& det) {

  if (A.num_rows() != A.num_cols()) {
    det = 0.0;
    return false;
  }

  int d = A.num_rows();
  Matrix<T> augmented(d, 2 * d);
  AugmentMatrix(A, augmented);
  det = 1.0;

  Subscript row, col, i, j;

  for (i = 0; i < d; i++) {
    col = row = i;
    T max = augmented[col][col];
    for (j = col+1; j < d; ++j) {
      if (max*max < augmented[j][col] * augmented[j][col]) {
        row = j;
        max = augmented[j][col];
      }
    }

    //If there are no non-zero values in the column, then A was not invertible
    if (max == 0) {
      det = 0.0;
      return false;
    }

    det *= (row == col) ? augmented[row][col] : -augmented[row][col];

    Pivot(augmented, row, col);
  }

  for (i = 0; i < d; ++i) {
    for (j = 0; j < d; ++j) {
      A[i][j] = augmented[i][j+d];
    }
  }

  return true;
}

template <class T>
void Pivot(Matrix<T>& A, Subscript M, Subscript N) {
  int i;
  Subscript n = A.num_rows();
  Subscript m = A.num_cols();

  if (M != N) {
    A.SwapRows(M, N);
  }

  T pivot = A[N][N];

  for (i = 0; i < m; ++i) {
    A[N][i] /= pivot;
  }
  A[N][N] = 1;

  T tmp;
  for (i = 0; i < n; ++i) {
    if (i != N) {
      tmp = A[i][N];
      A[i][N] = 0;
      for (int j = 0; j < m; j++) {
        A[i][j] -= A[N][j] * tmp;
      }
    }
  }
}

template <class T>
T Determinant(Matrix<T>& A) {

  if (A.num_rows() != A.num_cols()) {
    return T(0);
  }

  int d = A.num_rows();
  Matrix<T> augmented(d, 2 * d);
  AugmentMatrix(A, augmented);
  Subscript row, col, i, j;

  T det = T();
  for (i = 0; i < d; i++) {
    col = row = i;
    T max = augmented[col][col];
    for (j = col+1; j < d; ++j) {
      if (max*max < augmented[j][col] * augmented[j][col]) {
        row = j;
        max = augmented[j][col];
      }
    }

    //If there are no non-zero values in the column, then A was not invertible
    if (max == 0) {
      return T(0);
    }
    if (i == 0) {
      det = (row == col) ? augmented[row][col] : -augmented[row][col];
    } else {
      det *= (row == col) ? augmented[row][col] : -augmented[row][col];
    }

    Pivot(augmented, row, col);
  }

  return det;
}

template <class T>
Matrix<T> Minor(Matrix<T>& A, Subscript skip_row, Subscript skip_col) {
  int m = A.num_rows();
  int n = A.num_cols();
  Matrix<T> minorm(m - 1, n - 1);

  int i, j, xRow, xCol;
  for (i = 0, xRow = 0; i < m; ++i) {
    if (i == skip_row) {
      continue;
    }
    for (j = 0, xCol = 0; j < n; ++j) {
      if (j == skip_col) {
        continue;
      }
      minorm[xRow][xCol] = A[i][j];
      xCol++;
    }
    xRow++;
  }

  return minorm;
}

template <class T>
T DeterminantToo(Matrix<T>& A) {

  if (A.num_rows() == 1) {
    return A[0][0];
  }

  int d = A.num_rows();
  int sign = 1;
  T det = 0;

  for (int i = 0; i < d; ++i) {
    det += sign * A[0][i] * TNT::DeterminantToo(Minor(A, 0, i));
    sign *= -1;
  }

  return det;
}

inline int DeterminantSign(const std::vector<int>& perm) {
  int i = perm.size(), sign = 1, offset;
  while (i > 0) {
    offset = std::count_if(perm.begin() + i, perm.end(), std::bind1st(std::greater<int>(), perm[i]) );
    if ((i + perm[i] - offset) % 2 == 1) {
      sign *= -1;
    }
    --i;
  }

  return sign;
}

template <class T>
T DeterminantThree(Matrix<T>& A) {
  int i, d = A.num_rows();

  //Generate the indices in order
  std::vector<int> columnIndices;
  columnIndices.reserve(d);
  for (i = 0; i < d; ++i) {
    columnIndices.push_back(i);
  }

  T det = 0;
  do {
    T term = DeterminantSign(columnIndices);
    for (i = 0; i < d; ++i) {
      term *= A[i][columnIndices[i]];
    }
    det += term;
  } while (next_permutation(columnIndices.begin(), columnIndices.end()));

  return det;
}

} // namespace TNT

#endif
// CMAT_H

