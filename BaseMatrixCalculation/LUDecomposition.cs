namespace BaseMatrixCalculation
{
	class LUDecomposition
	{
		#region Class variables

		/// <summary>Array for internal storage of decomposition.
		/// @serial internal array storage.
		/// </summary>
		private double[][] LU;

		/// <summary>Row and column dimensions, and pivot sign.
		/// @serial column dimension.
		/// @serial row dimension.
		/// @serial pivot sign.
		/// </summary>
		private int m, n, pivsign;

		/// <summary>Internal storage of pivot vector.
		/// @serial pivot vector.
		/// </summary>
		private int[] piv;

		#endregion //  Class variables

		#region Constructor

		/// <summary>LU Decomposition</summary>
		/// <param name="A">  Rectangular matrix
		/// </param>
		/// <returns>     Structure to access L, U and piv.
		/// </returns>

		public LUDecomposition(Matrix A)
		{
			// Use a "left-looking", dot-product, Crout/Doolittle algorithm.

			LU = A.InternalArrayCopy;
			m = A.RowCount;
			n = A.ColumnCount;
			piv = new int[m];
			for (int i = 0; i < m; i++)
			{
				piv[i] = i;
			}
			pivsign = 1;
			double[] LUrowi;
			double[] LUcolj = new double[m];

			// Outer loop.

			for (int j = 0; j < n; j++)
			{

				// Make a copy of the j-th column to localize references.

				for (int i = 0; i < m; i++)
				{
					LUcolj[i] = LU[i][j];
				}

				// Apply previous transformations.

				for (int i = 0; i < m; i++)
				{
					LUrowi = LU[i];

					// Most of the time is spent in the following dot product.

					int kmax = System.Math.Min(i, j);
					double s = 0.0;
					for (int k = 0; k < kmax; k++)
					{
						s += LUrowi[k] * LUcolj[k];
					}

					LUrowi[j] = LUcolj[i] -= s;
				}

				// Find pivot and exchange if necessary.

				int p = j;
				for (int i = j + 1; i < m; i++)
				{
					if (System.Math.Abs(LUcolj[i]) > System.Math.Abs(LUcolj[p]))
					{
						p = i;
					}
				}
				if (p != j)
				{
					for (int k = 0; k < n; k++)
					{
						double t = LU[p][k]; LU[p][k] = LU[j][k]; LU[j][k] = t;
					}
					int k2 = piv[p]; piv[p] = piv[j]; piv[j] = k2;
					pivsign = -pivsign;
				}

				// Compute multipliers.

				if (j < m & LU[j][j] != 0.0)
				{
					for (int i = j + 1; i < m; i++)
					{
						LU[i][j] /= LU[j][j];
					}
				}
			}
		}
		#endregion //  Constructor

		#region Public Properties
		/// <summary>Is the matrix nonsingular?</summary>
		/// <returns>     true if U, and hence A, is nonsingular.
		/// </returns>
		virtual public bool IsNonSingular
		{
			get
			{
				for (int j = 0; j < n; j++)
				{
					if (LU[j][j] == 0)
						return false;
				}
				return true;
			}
		}

		/// <summary>Return lower triangular factor</summary>
		/// <returns>     L
		/// </returns>
		virtual public Matrix L
		{
			get
			{
				Matrix X = new Matrix(m, n);
				var L = X.InternalList;
				for (int i = 0; i < m; i++)
				{
					for (int j = 0; j < n; j++)
					{
						if (i > j)
						{
							L[i][j] = LU[i][j];
						}
						else if (i == j)
						{
							L[i][j] = 1.0;
						}
						else
						{
							L[i][j] = 0.0;
						}
					}
				}
				return X;
			}
		}

		/// <summary>Return upper triangular factor</summary>
		/// <returns>     U
		/// </returns>
		virtual public Matrix U
		{
			get
			{
				Matrix X = new Matrix(n, n);
				var U = X.InternalList;
				for (int i = 0; i < n; i++)
				{
					for (int j = 0; j < n; j++)
					{
						if (i <= j)
						{
							U[i][j] = LU[i][j];
						}
						else
						{
							U[i][j] = 0.0;
						}
					}
				}
				return X;
			}
		}

		/// <summary>Return pivot permutation vector</summary>
		/// <returns>     piv
		/// </returns>
		virtual public int[] Pivot
		{
			get
			{
				int[] p = new int[m];
				for (int i = 0; i < m; i++)
				{
					p[i] = piv[i];
				}
				return p;
			}
		}

		/// <summary>Return pivot permutation vector as a one-dimensional double array</summary>
		/// <returns>     (double) piv
		/// </returns>
		virtual public double[] DoublePivot
		{
			get
			{
				double[] vals = new double[m];
				for (int i = 0; i < m; i++)
				{
					vals[i] = (double)piv[i];
				}
				return vals;
			}
		}

		#endregion //  Public Properties

		#region Public Methods

		/// <summary>Determinant</summary>
		/// <returns>     det(A)
		/// </returns>
		/// <exception cref="System.ArgumentException">  Matrix must be square
		/// </exception>

		public virtual double Determinant()
		{
			if (m != n)
			{
				throw new System.ArgumentException("Matrix must be square.");
			}
			double d = (double)pivsign;
			for (int j = 0; j < n; j++)
			{
				d *= LU[j][j];
			}
			return d;
		}

		/// <summary>Solve A*X = B</summary>
		/// <param name="B">  A Matrix with as many rows as A and any number of columns.
		/// </param>
		/// <returns>     X so that L*U*X = B(piv,:)
		/// </returns>
		/// <exception cref="System.ArgumentException"> Matrix row dimensions must agree.
		/// </exception>
		/// <exception cref="System.SystemException"> Matrix is singular.
		/// </exception>

		public virtual Matrix Solve(Matrix B)
		{
			if (B.RowCount != m)
			{
				throw new System.ArgumentException("Matrix row dimensions must agree.");
			}
			if (!this.IsNonSingular)
			{
				throw new System.SystemException("Matrix is singular.");
			}

			// Copy right hand side with pivoting
			int nx = B.ColumnCount;
			Matrix Xmat = GetMatrix(B, piv, 0, nx - 1);
			var X = Xmat.InternalList;

			// Solve L*Y = B(piv,:)
			for (int k = 0; k < n; k++)
			{
				for (int i = k + 1; i < n; i++)
				{
					for (int j = 0; j < nx; j++)
					{
						X[i][j] -= X[k][j] * LU[i][k];
					}
				}
			}
			// Solve U*X = Y;
			for (int k = n - 1; k >= 0; k--)
			{
				for (int j = 0; j < nx; j++)
				{
					X[k][j] /= LU[k][k];
				}
				for (int i = 0; i < k; i++)
				{
					for (int j = 0; j < nx; j++)
					{
						X[i][j] -= X[k][j] * LU[i][k];
					}
				}
			}
			return Xmat;
		}

		#endregion //  Public Methods

		private Matrix GetMatrix(Matrix A, int[] r, int j0, int j1)
		{
			Matrix X = new Matrix(r.Length, j1 - j0 + 1);
			var B = X.InternalList;
			try
			{
				for (int i = 0; i < r.Length; i++)
				{
					for (int j = j0; j <= j1; j++)
					{
						B[i][j - j0] = A[r[i], j];
					}
				}
			}
			catch (System.IndexOutOfRangeException e)
			{
				throw new System.IndexOutOfRangeException("Submatrix indices", e);
			}
			return X;
		}
	}
}
