using System;
using System.Collections.Generic;


namespace BaseMatrixCalculation
{
	public class Matrix //: IEnumerable<IEnumerable<double>>
    {
        private static Random _rnd = new Random();
        private List<List<double>> _matrix;

        
        #region Constructors

        public Matrix()
        {
            _matrix = new List<List<double>>();
        }

        public Matrix(int m, int n)
        {
            _matrix = new List<List<double>>();
            for (int i = 0; i < m; i++)
            {
                _matrix.Add(new List<double>());
				for (int j = 0; j < n; j++)
				{
					_matrix[i].Add(0);
				}
            }
        }

		public Matrix(int dim) : this(dim, dim)
		{
		}

		public Matrix(double[] source, MatrixType type)
		{
			if(type == MatrixType.Row)
			{
				_matrix = new List<List<double>>();
				_matrix.Add(new List<double>());
				for (int j = 0; j < source.Length; j++)
				{
					_matrix[0].Add(source[j]);
				}
			}
			else
			{
				_matrix = new List<List<double>>();
				for (int i = 0; i < source.Length; i++)
				{
					_matrix.Add(new List<double>());
					_matrix[i].Add(source[i]);
				}
			}
		}

		public Matrix(List<double> source, MatrixType type) : this(source.ToArray(), type)
		{
		}

		public Matrix(double[,] matrix) : this(matrix.GetLength(0), matrix.GetLength(1))
        {
            for (int i = 0; i < RowCount; i++)
            {
                for (int j = 0; j < ColumnCount; j++)
                {
                    _matrix[i][j] = matrix[i, j];
                }
            }
        }

        public Matrix(double[][] matrix) : this(matrix.Length, matrix[0].Length)
        {
            for (int i = 0; i < matrix.Length - 1; i++)
            {
				if (matrix[i].Length != matrix[i + 1].Length)
				{
					throw new ArgumentException("Matrix must be a rectangle");
				}
            }

            for (int i = 0; i < RowCount; i++)
            {
                for (int j = 0; j < ColumnCount; j++)
                {
                    _matrix[i][j] = matrix[i][j];
                }
            }
        }

        public Matrix(List<double[]> matrix) : this(matrix.Count, matrix[0].Length)
        {
            for (int i = 0; i < matrix.Count - 1; i++)
            {
                if (matrix[i].Length != matrix[i + 1].Length)
                    throw new ArgumentException("Matrix must be a rectangle");
            }

            for (int i = 0; i < RowCount; i++)
            {
                for (int j = 0; j < ColumnCount; j++)
                {
                    _matrix[i][j] = matrix[i][j];
                }
            }
        }

        public Matrix(List<List<double>> matrix) : this(matrix.Count, matrix[0].Count)
        {
            for (int i = 0; i < matrix.Count - 1; i++)
            {
                if (matrix[i].Count != matrix[i + 1].Count)
                    throw new ArgumentException("Matrix must be a rectangle");
            }

            for (int i = 0; i < RowCount; i++)
            {
                for (int j = 0; j < ColumnCount; j++)
                {
                    _matrix[i][j] = matrix[i][j];
                }
            }
        }

        #endregion


        public double this[int i, int j]
        {
            get
            {
                if (i < RowCount && i >= 0 && j < ColumnCount && j >= 0)
                {
                    return _matrix[i][j];
                }
                else
                {
                    throw new ArgumentOutOfRangeException("Attempt to adress out of matrix boundary");
                }
            }
            set
            {
                if (i < RowCount && i >= 0 && j < ColumnCount && j >= 0)
                {
                    _matrix[i][j] = value;
                }
                else
                {
                    throw new ArgumentOutOfRangeException("Attempt to adress out of matrix boundary");
                }
            }
        }


		#region Public Properties

		public int RowCount
		{
			get
			{
				return _matrix.Count;
			}
		}

		public int ColumnCount
		{
			get
			{
				return RowCount > 0 ? _matrix[0].Count : 0;
			}
		}

		public List<List<double>> InternalList
		{
			get
			{
				return _matrix;
			}
		}

		public List<List<double>> InternalListCopy
		{
			get
			{
				var res = new List<List<double>>();
				for (int i = 0; i < RowCount; i++)
				{
					res.Add(new List<double>());
					for (int j = 0; j < ColumnCount; j++)
					{
						res[i].Add(_matrix[i][j]);
					}
				}
				return res;
			}
		}

		public double[][] InternalArrayCopy
		{
			get
			{
				var res = new double[RowCount][];
				for (int i = 0; i < RowCount; i++)
				{
					res[i] = new double[ColumnCount];
					for (int j = 0; j < ColumnCount; j++)
					{
						res[i][j] = _matrix[i][j];
					}
				}
				return res;
			}
		}

		public double Determinant
		{
			get
			{
				return new LUDecomposition(this).Determinant();
			}
		}

		#endregion


		#region Public Methods

		/// <summary>
		/// Returns submatrix with specified indexes: [i0,j0] - top-left point, [i1,j1] - bot-right point
		/// </summary>
		/// <param name="i0"></param>
		/// <param name="j0"></param>
		/// <param name="i1"></param>
		/// <param name="j1"></param>
		/// <returns></returns>
		public Matrix GetSubMatrix(int i0, int j0, int i1, int j1)
		{
			var subMatrix = new Matrix(i1 - i0 + 1, j1 - j0 + 1);
			for (int i = i0; i <= i1; i++)
			{
				for (int j = j0; j <= j1; j++)
				{
					subMatrix[i - i0, j - j0] = _matrix[i][j];
				}
			}
			return subMatrix;
		}

		/// <summary>
		/// Rank of matrix
		/// </summary>
		/// <returns></returns>
		public int Rank()
		{
			return new SingularValueDecomposition(this).Rank();
		}

		/// <summary>
		/// Returns new matrix without row and column with specified index
		/// </summary>
		/// <param name="ind"></param>
		/// <returns></returns>
		public Matrix DeleteRowAndColumn(int ind)
		{
			var res = new Matrix(RowCount - 1, ColumnCount - 1);
			var resRowInd = 0;
			var resColInd = 0;
			for (int i = 0; i < RowCount; i++)
			{
				for (int j = 0; j < ColumnCount; j++)
				{
					if(j != ind)
					{
						res[resRowInd, resColInd] = this[i, j];
						resColInd++;
					}
				}
				if(i != ind)
				{
					resRowInd++;
				}
			}
			return res;
		}

		/// <summary>
		/// Returns matrix with random integer values from minValue to maxValue
		/// </summary>
		/// <param name="rowCount"></param>
		/// <param name="columnCount"></param>
		/// <param name="minValue"></param>
		/// <param name="maxValue"></param>
		/// <returns></returns>
		public static Matrix Random(int rowCount, int columnCount, int minValue, int maxValue)
        {
            Matrix res = new Matrix(rowCount, columnCount);
            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < columnCount; j++)
                {
					res[i, j] = _rnd.Next(minValue, maxValue + 1);
                }
            }
            return res;
        }

		/// <summary>
		/// Returns matrix with random double values from 0 to 1
		/// </summary>
		/// <param name="rowCount"></param>
		/// <param name="columnCount"></param>
		/// <returns></returns>
		public static Matrix Random(int rowCount, int columnCount)
		{
			Matrix res = new Matrix(rowCount, columnCount);
			for (int i = 0; i < rowCount; i++)
			{
				for (int j = 0; j < columnCount; j++)
				{
					res[i, j] = _rnd.NextDouble();
				}
			}
			return res;
		}

		/// <summary>
		/// Return transposed matrix
		/// </summary>
		/// <returns></returns>
        public Matrix GetTranspose()
        {
			var res = new Matrix(ColumnCount, RowCount);
            for (int i = 0; i < ColumnCount; i++)
            {
                for (int j = 0; j < RowCount; j++)
                {
                    res[i, j] = this[j, i];
                }
            }
			return res;
        }


		/// <summary>Solve A*X = B</summary>
		/// <param name="B">   right hand side
		/// </param>
		/// <returns>     solution if A is square, least squares solution otherwise
		/// </returns>
		public Matrix Solve(Matrix B)
		{
			return (RowCount == ColumnCount ? (new LUDecomposition(this)).Solve(B) : (new QRDecomposition(this)).Solve(B));
		}

		/// <summary>
		/// Returns inverse matrix
		/// </summary>
		/// <param name="A"></param>
		/// <returns></returns>
		public Matrix GetInverse()
        {
			return Solve(E(RowCount));

			//    List<List<double>> dopoln = new List<List<double>>();
			//    int n = RowCount;
			//    for (int i = 0; i < n; i++)
			//    {
			//        dopoln.Add(new List<double>());
			//        for (int j = 0; j < n; j++)
			//        {
			//            dopoln[i].Add(this[i, j]);
			//        }
			//    }
			//    for (int i = 0; i < n; i++)
			//    {
			//        for (int j = 0; j < n; j++)
			//        {
			//            dopoln[i][j] = 0;
			//        }
			//    }
			//    for (int i = 0; i < n; i++)
			//    {

			//        for (int j = 0; j < n; j++)
			//        {
			//            List<List<double>> minor = new List<List<double>>();
			//            for (int k = 0; k < n; k++)
			//            {
			//                minor.Add(new List<double>());
			//                for (int count = 0; count < n; count++)
			//                {
			//                    minor[k].Add(this[k, count]);
			//                }
			//            }
			//            minor.RemoveAt(i);
			//            for (int k = 0; k < minor.Count; k++)
			//            {
			//                minor[k].RemoveAt(j);
			//            }
			//            double s11 = Math.Pow(-1, i + 1 + j + 1);
			//            Matrix tmpminor = new Matrix(minor.Count, minor.Count);
			//            for (int count = 0; count < minor.Count; count++)
			//            {
			//                for (int count1 = 0; count1 < minor.Count; count1++)
			//                {
			//                    tmpminor[count, count1] = minor[count][count1];
			//                }
			//            }
			//            double s22 = tmpminor.Determinant;
			//            dopoln[i][j] = s11 * s22;
			//            ;
			//        }
			//    }
			//    var tmp = new Matrix(dopoln);
			//    tmp = tmp.GetTranspose();
			//    double[,] inversed = new double[n, n];
			//    double mn = Determinant;
			//    for (int i = 0; i < n; i++)
			//    {
			//        for (int j = 0; j < n; j++)
			//        {
			//            inversed[i, j] = tmp[i, j] / mn;
			//        }
			//    }

			//    return new Matrix(inversed);
		}

		/// <summary>
		/// Returns identity matrix of specified dimension
		/// </summary>
		/// <param name="dimension"></param>
		/// <returns></returns>
		public static Matrix E(int dimension)
		{
			var res = new Matrix(dimension);
			for (int i = 0; i < dimension; i++)
			{
				for (int j = 0; j < dimension; j++)
				{
					if(i == j)
					{
						res[i, j] = 1;
					}
				}
			}
			return res;
		}

		/// <summary>
		/// Returns empty matrix (filled by zeros) of specified dimension
		/// </summary>
		/// <param name="dimension"></param>
		/// <returns></returns>
		public static Matrix Empty(int dimension)
		{
			return new Matrix(dimension);
		}

		/// <summary>
		/// Matrix multiplication
		/// </summary>
		/// <param name="A"></param>
		/// <param name="B"></param>
		/// <returns></returns>
        private static Matrix Multi(Matrix A, Matrix B)
        {
			if (A.ColumnCount != B.RowCount)
			{
				throw new ArgumentException("Matrixes inner dimensions must agree.");
			}

			int m = A.RowCount;
            int n = A.ColumnCount;
            int q = B.ColumnCount;
			
			var result = new Matrix(m, q);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < q; j++)
                {
                    double c = 0;
                    for (int r = 0; r < n; r++)
                    {
                        c += A[i, r] * B[r, j];
                    }
                    result[i, j] = c;
                }
            }


            return result;
        }

		/// <summary>
		/// Add each element of first matrix to second. Matrixes must be equal size
		/// </summary>
		/// <param name="A"></param>
		/// <param name="B"></param>
		/// <returns></returns>
		public static Matrix Sum(Matrix A, Matrix B)
        {
			int n = A.ColumnCount;
            var result = new Matrix(n, n);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    result[i, j] = A[i, j] + B[i, j];
                }
            }

            return result;
        }

		/// <summary>
		/// Substitute all elements of second matrix from first
		/// </summary>
		/// <param name="A"></param>
		/// <param name="B"></param>
		/// <returns></returns>
        private static Matrix Substr(Matrix A, Matrix B)
        {
            int n = A.ColumnCount;
            var result = new Matrix(n, n);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    result[i, j] = A[i, j] - B[i, j];
                }
            }

            return result;
        }

		public Matrix Copy()
		{
			var clone = new Matrix(RowCount, ColumnCount);
			for (int i = 0; i < RowCount; i++)
			{
				for (int j = 0; j < ColumnCount; j++)
				{
					clone[i, j] = this[i, j];
				}
			}
			return clone;
		}

		public void AddRow(List<double> list)
		{
			if (list.Count != ColumnCount)
			{
				throw new ArgumentException();
			}

			_matrix.Add(list);
		}

		public void AddColumn(List<double> list)
		{
			if (list.Count != RowCount)
			{
				throw new ArgumentException();
			}

			for (int j = 0; j < RowCount; j++)
			{
				_matrix[j].Add(list[j]);
			}
		}

		//public IEnumerator<IEnumerable<double>> GetEnumerator()
		//{
		//	return _matrix.GetEnumerator();
		//}

		//IEnumerator IEnumerable.GetEnumerator()
		//{
		//	return GetEnumerator();
		//}

		#endregion


		#region Operators

		public static bool operator ==(Matrix left, Matrix right)
		{
			if (left.RowCount != right.RowCount || left.ColumnCount != right.ColumnCount)
			{
				return false;
			}

			for (int i = 0; i < left.RowCount; i++)
			{
				for (int j = 0; j < left.ColumnCount; j++)
				{
					if (left[i, j] != right[i, j])
					{
						return false;
					}
				}
			}
			return true;
		}

		public static bool operator !=(Matrix left, Matrix right)
		{
			return !(left == right);
		}

		public static Matrix operator *(Matrix left, Matrix right)
		{
			return Multi(left, right);
		}

		public static Matrix operator +(Matrix left, Matrix right)
		{
			if(left.RowCount != right.RowCount || left.ColumnCount != right.ColumnCount)
			{
				throw new ArgumentException("Matrix dimensions must be equal");
			}

			var res = new Matrix(left.RowCount, left.ColumnCount);
			for (int i = 0; i < res.RowCount; i++)
			{
				for (int j = 0; j < res.ColumnCount; j++)
				{
					res[i, j] = left[i, j] + right[i, j];
				}
			}
			return res;
		}

		public static Matrix operator -(Matrix left, Matrix right)
		{
			return left + -right;
		}

		public static Matrix operator -(Matrix left)
		{
			var res = left.Copy();
			for (int i = 0; i < res.RowCount; i++)
			{
				for (int j = 0; j < res.ColumnCount; j++)
				{
					res[i, j] = -res[i, j];
				}
			}
			return res;
		}
		
		public static Matrix operator *(Matrix left, double number)
		{
			var res = left.Copy();
			for (int i = 0; i < res.RowCount; i++)
			{
				for (int j = 0; j < res.ColumnCount; j++)
				{
					res[i, j] = left[i, j] * number;
				}
			}
			return res;
		}

		public static Matrix operator /(Matrix left, double number)
		{
			var res = left.Copy();
			for (int i = 0; i < res.RowCount; i++)
			{
				for (int j = 0; j < res.ColumnCount; j++)
				{
					res[i, j] /= number;
				}
			}
			return res;
		}

		public static Matrix operator +(Matrix left, double number)
		{
			var res = left.Copy();
			for (int i = 0; i < res.RowCount; i++)
			{
				for (int j = 0; j < res.ColumnCount; j++)
				{
					res[i, j] += number;
				}
			}
			return res;
		}

		public static Matrix operator -(Matrix left, double number)
		{
			var res = left.Copy();
			for (int i = 0; i < res.RowCount; i++)
			{
				for (int j = 0; j < res.ColumnCount; j++)
				{
					res[i, j] -= number;
				}
			}
			return res;
		}

		public static Matrix operator *(double number, Matrix left)
		{
			return left * number;
		}

		public static Matrix operator ^(Matrix left, int pow)
		{
			if(pow == 0)
			{
				return E(left.RowCount);
			}

			var res = left.Copy();
			for (int i = 0; i < pow; i++)
			{
				res *= left;
			}
			return res;
		}

		#endregion
	}
}
