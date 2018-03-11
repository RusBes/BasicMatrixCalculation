Very simple library for matrix calculations

Just a convinient way for me to use matrixes. It contains very basic operations, but soon it will be more.
Some functions I took from this project: https://www.codeproject.com/Articles/5835/DotNetMatrix-Simple-Matrix-Library-for-NET
and some is mine. I include some constructors and operators, so you can use it lik this:

var matrix = new Matrix(new double[,] { { 1, 2 }, { 3, 4 } });
var result = ((2 * matrix.GetInverse() * matrix.GetTranspose()) ^ 2).Determinant;
