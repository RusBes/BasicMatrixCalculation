using Microsoft.VisualStudio.TestTools.UnitTesting;
using BaseMatrixCalculation;

namespace MatrixLibraryTest
{
	[TestClass]
	public class MatrixTest
	{
		private Matrix _matrix1 = new Matrix(new double[,] { { 1, 2, 3 }, { 1, 2, 2 }, { 2, 4, 5 } });

		private Matrix _matrix2 = new Matrix(new double[,] { { 1, 2, 3 }, { 5, 6, 2 }, { 5, 2, 6 } });

		[TestMethod]
		public void TestMultiplyByMatrix()
		{
			var expected = new Matrix(new double[,] { { 26, 20, 25 }, { 21, 18, 19 }, { 47, 38, 44 } });
			var actual = _matrix1 * _matrix2;

			Assert.IsTrue(expected == actual);
		}

		[TestMethod]
		public void TestMultiplyByNumber()
		{
			var expected = new Matrix(new double[,] { { 2, 4, 6 }, { 2, 4, 4 }, { 4, 8, 10 } });
			var actual1 = _matrix1 * 2;
			var actual2 = 2 * _matrix1;

			Assert.IsTrue(expected == actual1);
			Assert.IsTrue(expected == actual2);
		}

		[TestMethod]
		public void TestSummurizeByMatrix()
		{
			var expected = new Matrix(new double[,] { { 2, 4, 6 }, { 6, 8, 4 }, { 7, 6, 11 } });
			var actual1 = _matrix1 + _matrix2;
			var actual2 = _matrix2 + _matrix1;

			Assert.IsTrue(expected == actual1);
			Assert.IsTrue(expected == actual2);
		}

		[TestMethod]
		public void TestSummurizeByNumber()
		{
			var expected = new Matrix(new double[,] { { 11, 12, 13 }, { 11, 12, 12 }, { 12, 14, 15 } });
			var actual = _matrix1 + 10;

			Assert.IsTrue(expected == actual);
		}

		[TestMethod]
		public void TestDeterminant()
		{
			var expected = 8;
			// Almost equal to _matrix1, but first element is 5
			var matrix = new Matrix(new double[,] { { 5, 2, 3 }, { 1, 2, 2 }, { 2, 4, 5 } });
			var actual = matrix.Determinant;

			Assert.AreEqual(expected, actual);
		}

		[TestMethod]
		public void TestInverse()
		{
			var expected = new Matrix(new double[,] { { -2, 1 }, { 1.5, -0.5 } });
			var matrix = new Matrix(new double[,] { { 1, 2 }, { 3, 4 } });
			var actual = matrix.GetInverse();
			// set epsilon, because of computational mistakes. The result must be in range [expected - eps, expected + eps]
			var eps = 0.0001;

			for (int i = 0; i < expected.RowCount; i++)
			{
				for (int j = 0; j < expected.ColumnCount; j++)
				{
					Assert.IsTrue(expected[i, j] + eps >= actual[i, j]
						&& expected[i, j] - eps <= actual[i, j]);
				}
			}
		}

		//[TestMethod]
		//[ExpectedException(typeof(ArgumentException), "Haven't exception when accessing incorrect element")]
		//public void TestElementAccessing()
		//{
		//	var expected = new Matrix(new double[,] { { 26, 20, 25 }, { 21, 18, 19 }, { 47, 38, 44 } });
		//	var actual = _matrix1 * _matrix2;

		//	Assert.AreEqual(expected, actual);
		//}
	}
}
