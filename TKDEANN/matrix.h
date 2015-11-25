// change the extension of this file to .H

#ifndef H__MATRIX
#define H__MATRIX 1

class matrix
{
  private:
    double *cells;
  public:
    matrix(int _height, int _width);
    matrix(const matrix& A);
    void operator = (matrix);
    ~matrix();
    int height,width;
    matrix transpose();
    matrix inverse();
    double& operator () (int, int) const;
    void print();
};

void matrix::print() {
  for (int i = 1; i <= height; i++) {
    for (int j = 1; j <= width; j++)
		printf("%f\t",(*this)(i,j));
    printf("\n");
  }
}

matrix::matrix(int _height, int _width) {
  height = _height;
  width = _width;
  cells = new double[height*width];
  for (int i=0;i<height*width;i++)
  	cells[i]=0;
}

matrix::matrix(const matrix &A) {
  height = A.height;
  width = A.width;
  cells = new double[height*width];
  memcpy(cells, A.cells, height*width*sizeof(double));
}

void matrix::operator = (matrix A) {
  assert( (height != A.height || width != A.width)==false);
    //throw matrix::exception("Incompatible matrix dimensions");
  height = A.height;
  width = A.width;
  memcpy(cells, A.cells, height*width*sizeof(double));
}

matrix::~matrix() {
  delete[] cells;
}

// column-wise storage
double& matrix::operator () (int row, int column) const {
  	return cells[(row-1)*width + column - 1];
}

matrix operator * (matrix A, matrix B) {
	assert( (A.width != B.height)==false);
		//throw matrix::exception("Incompatible matrix dimensions");
	matrix AB(A.height, B.width);
	for (int i = 1; i <= AB.height; i++) {
		for (int j = 1; j <= AB.width; j++) {
			AB(i,j) = 0.0;
			for (int k = 1; k <= A.width; k++)
				AB(i,j) = AB(i,j) + A(i,k) * B(k,j);
		}
	}
	return AB;
}

matrix matrix::transpose() {
  matrix T(width, height);
  for (int i = 1; i <= T.height; i++)
    for (int j = 1; j <= T.width; j++)
      T(i,j) = (*this)(j,i);
  return T;
}

matrix matrix::inverse() {
  assert( (height != width)==false);
    //throw exception("Inversion of non-square matrix");
  const int N = height;
  matrix B(*this), R(N,N);
  
  // assume 0 init.
  for (int i = 1; i <= N; i++)
    R(i,i) = 1.0;
  
  for (int k = 1; k <= N; k++) {
    int pivot_row = 0;
    double pivot_value,abs_pivot_value;
    for (int i = k; i <= N; i++)
    {
      double x = B(i,k);
      double absx = x;
      if (absx < 0)
        absx = -absx;
      if (pivot_row == 0 || absx > abs_pivot_value) {
        pivot_row = i;
        abs_pivot_value = absx;
        pivot_value = x;
      }
      assert( (pivot_value == 0.0)==false);
        //throw exception("Inversion of singular matrix");
    }
    
    if (pivot_row == k) {
      for (int j = k; j <= N; j++)
        B(k,j) = B(k,j) / pivot_value;
      for (int j = 1; j <= N; j++)
        R(k,j) = R(k,j) / pivot_value;
    } else {
      for (int j = k; j <= N; j++)
      {
        double x = B(k,j);
        B(k,j) = B(pivot_row,j) / pivot_value;
        B(pivot_row,j) = x;
      }
      for (int j = 1; j <= N; j++)
      {
        double x = R(k,j);
        R(k,j) = R(pivot_row,j) / pivot_value;
        R(pivot_row,j) = x;
      }
    }
    for (int i = 1; i <= N; i++) {
      if (i != k) {
        for (int j = k+1; j <= N; j++)
          B(i,j) = B(i,j) - B(i,k) * B(k,j);
        for (int j = 1; j <= N; j++)
          R(i,j) = R(i,j) - B(i,k) * R(k,j);
      }
    }
  }
  return R;
}

#endif
