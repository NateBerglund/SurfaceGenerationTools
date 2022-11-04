using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace SurfaceGenerator
{
    /// <summary>
    /// Set of utility functions for generating STL data as a 3n x 3 matrix of vertex coordinates,
    /// 3 per triangle.
    /// </summary>
    public static class SurfaceGenerator
    {
        #region Parametric Surface Generation

        /// <summary>
        /// Generates the vertices of a quad surface, 3 vertices per triangle.
        /// </summary>
        /// <param name="f">Function describing the (x,y,z) position of the surface as a function of two input parameters, t and u.</param>
        /// <param name="tStart">Starting value of t</param>
        /// <param name="tEnd">Ending value of t</param>
        /// <param name="nT">Number of quadrilaterals to generate in the t-direction</param>
        /// <param name="uStart">Starting value of u</param>
        /// <param name="uEnd">Ending value of u</param>
        /// <param name="nU">Number of quadrilaterals to generate in the u-direction</param>
        /// <returns>Matrix of vertex coordinates, which can be used to generate STL data</returns>
        public static Matrix<double> GenerateQuadSurface(
            Func<double[], double[], Matrix<double>> f,
            double tStart, double tEnd, int nT,
            double uStart, double uEnd, int nU)
        {
            List<double> tParams = Generate.LinearSpaced(nT + 1, tStart, tEnd).ToList();
            List<double> tParamsShifted = new List<double>(tParams);
            tParamsShifted.RemoveAt(0);
            tParams.RemoveAt(tParams.Count - 1);

            List<double> uParams = Generate.LinearSpaced(nU + 1, uStart, uEnd).ToList();
            List<double> uParamsShifted = new List<double>(uParams);
            uParamsShifted.RemoveAt(0);
            uParams.RemoveAt(uParams.Count - 1);

            Grid(tParams.ToArray(), uParams.ToArray(), out Matrix<double> tGrid1, out Matrix<double> uGrid1);
            Grid(tParams.ToArray(), uParamsShifted.ToArray(), out Matrix<double> tGrid2, out Matrix<double> uGrid2);
            Grid(tParamsShifted.ToArray(), uParamsShifted.ToArray(), out Matrix<double> tGrid3, out Matrix<double> uGrid3);
            Grid(tParamsShifted.ToArray(), uParams.ToArray(), out Matrix<double> tGrid4, out Matrix<double> uGrid4);

            int nQuads = nT * nU;
            Matrix<double> vertices = Matrix<double>.Build.DenseOfArray(new double[6 * nQuads, 3]);
            SetSubMatrix(vertices, 0, 6, nQuads, 0, 1, 3, f(Flatten(tGrid1), Flatten(uGrid1)));
            SetSubMatrix(vertices, 1, 6, nQuads, 0, 1, 3, f(Flatten(tGrid2), Flatten(uGrid2)));
            SetSubMatrix(vertices, 2, 6, nQuads, 0, 1, 3, f(Flatten(tGrid3), Flatten(uGrid3)));
            SetSubMatrix(vertices, 3, 6, nQuads, 0, 1, 3, f(Flatten(tGrid3), Flatten(uGrid3)));
            SetSubMatrix(vertices, 4, 6, nQuads, 0, 1, 3, f(Flatten(tGrid4), Flatten(uGrid4)));
            SetSubMatrix(vertices, 5, 6, nQuads, 0, 1, 3, f(Flatten(tGrid1), Flatten(uGrid1)));
            return vertices;
        }

        /// <summary>
        /// Generates the vertices of a fan surface, 3 vertices per triangle. A fan surface differs from a quad surface
        /// in that we assume there is a singularity in the function at u = uStart, such that f(t, uStart) is always the
        /// same point, regarless of what the value of t is.
        /// </summary>
        /// <param name="f">Function describing the (x,y,z) position of the surface as a function of two input parameters, t and u.</param>
        /// <param name="tStart">Starting value of t</param>
        /// <param name="tEnd">Ending value of t</param>
        /// <param name="nT">Number of quadrilaterals to generate in the t-direction</param>
        /// <param name="uStart">Starting value of u</param>
        /// <param name="uEnd">Ending value of u</param>
        /// <param name="nU">Number of quadrilaterals to generate in the u-direction</param>
        /// <returns>Matrix of vertex coordinates, which can be used to generate STL data</returns>
        public static Matrix<double> GenerateFanSurface(
            Func<double[], double[], Matrix<double>> f,
            double tStart, double tEnd, int nT,
            double uStart, double uEnd, int nU)
        {
            double uNewStart = uStart + (uEnd - uStart) / nU;
            Matrix<double> quadPart = GenerateQuadSurface(
                f,
                tStart, tEnd, nT,
                uNewStart, uEnd, nU - 1);

            Matrix<double> centerMatrix = f(new double[] { tStart }, new double[] { uStart });
            double[] centerPoint = new double[] { centerMatrix[0, 0], centerMatrix[0, 1], centerMatrix[0, 2] };

            Matrix<double> fanPart = GenerateBasicFan(
                tValues => f(tValues, Enumerable.Repeat(uNewStart, tValues.Length).ToArray()),
                tStart, tEnd, nT,
                centerPoint);

            return Matrix<double>.Build.DenseOfMatrixArray(new Matrix<double>[,] { { fanPart }, { quadPart } });
        }

        /// <summary>
        /// Generates the vertices of a "crescent moon" shaped surface, 3 vertices per triangle. A crescent surface differs from a quad surface
        /// in that we assume there are two singularities in the function, one at u = uStart and one at u = uEnd, such that f(t, uStart) is always the
        /// same point, regarless of what the value of t is, and such that f(t, uEnd) is always the same point, regarless of what the value of t is.
        /// </summary>
        /// <param name="f">Function describing the (x,y,z) position of the surface as a function of two input parameters, t and u.</param>
        /// <param name="tStart">Starting value of t</param>
        /// <param name="tEnd">Ending value of t</param>
        /// <param name="nT">Number of quadrilaterals to generate in the t-direction</param>
        /// <param name="uStart">Starting value of u</param>
        /// <param name="uEnd">Ending value of u</param>
        /// <param name="nU">Number of quadrilaterals to generate in the u-direction</param>
        /// <returns>Matrix of vertex coordinates, which can be used to generate STL data</returns>
        public static Matrix<double> GenerateCrescentSurface(
            Func<double[], double[], Matrix<double>> f,
            double tStart, double tEnd, int nT,
            double uStart, double uEnd, int nU)
        {
            double uNewStart = uStart + (uEnd - uStart) / nU;
            double uNewEnd = uEnd - (uEnd - uStart) / nU;
            Matrix<double> quadPart = GenerateQuadSurface(
                f,
                tStart, tEnd, nT,
                uNewStart, uNewEnd, nU - 2);

            Matrix<double> centerMatrix1 = f(new double[] { tStart }, new double[] { uStart });
            double[] centerPoint1 = new double[] { centerMatrix1[0, 0], centerMatrix1[0, 1], centerMatrix1[0, 2] };
            Matrix<double> centerMatrix2 = f(new double[] { tStart }, new double[] { uEnd });
            double[] centerPoint2 = new double[] { centerMatrix2[0, 0], centerMatrix2[0, 1], centerMatrix2[0, 2] };

            Matrix<double> fanPart1 = GenerateBasicFan(
                tValues => f(tValues, Enumerable.Repeat(uNewStart, tValues.Length).ToArray()),
                tStart, tEnd, nT,
                centerPoint1);

            Matrix<double> fanPart2 = GenerateBasicFan(
                tValues => f(tValues, Enumerable.Repeat(uNewEnd, tValues.Length).ToArray()),
                tEnd, tStart, nT,
                centerPoint2);

            return Matrix<double>.Build.DenseOfMatrixArray(new Matrix<double>[,] { { fanPart1 }, { quadPart }, { fanPart2 } });
        }

        /// <summary>
        /// Generates the vertices of a simple triangle fan, 3 vertices per triangle. A triangle fan as a single
        /// center point, and "fans out" to connect to the curve described by the 1-parameter input function f.
        /// </summary>
        /// <param name="f">Function describing the (x,y,z) position of the fan boundary as a function of one input parameter t.</param>
        /// <param name="tStart">Starting value of t</param>
        /// <param name="tEnd">Ending value of t</param>
        /// <param name="nT">Number of triangles to generate</param>
        /// <returns>Matrix of vertex coordinates, which can be used to generate STL data</returns>
        public static Matrix<double> GenerateBasicFan(
            Func<double[], Matrix<double>> f,
            double tStart, double tEnd, int nT,
            double[] centerPoint)
        {
            List<double> tParams = Generate.LinearSpaced(nT + 1, tStart, tEnd).ToList();
            List<double> tParamsShifted = new List<double>(tParams);
            tParamsShifted.RemoveAt(0);
            tParams.RemoveAt(tParams.Count - 1);

            Matrix<double> vertices = Matrix<double>.Build.DenseOfArray(new double[3 * nT, 3]);
            SetSubMatrix(vertices, 0, 3, nT, 0, 1, 3, Matrix<double>.Build.Dense(nT, 3, (i, j) => centerPoint[j]));
            SetSubMatrix(vertices, 1, 3, nT, 0, 1, 3, f(tParams.ToArray()));
            SetSubMatrix(vertices, 2, 3, nT, 0, 1, 3, f(tParamsShifted.ToArray()));
            return vertices;
        }

        #endregion Parametric Surface Generation

        #region Polyhedral Surface Generation

        /// <summary>
        /// Generates the vertices of an STL surface in the shape of a disdyakis dodecahedron
        /// </summary>
        /// <param name="radius">Desired radius of the vertices (from the origin)</param>
        /// <returns>Matrix of vertex coordinates, which can be used to generate STL data</returns>
        public static Matrix<double> GenerateDisdyakisDodecahedron(double radius)
        {
            double a = 4 / (1 + 2 * Math.Sqrt(2));
            double b = 4 / (2 + 3 * Math.Sqrt(2));
            double c = 4 / Math.Sqrt(27 + 18 * Math.Sqrt(2));

            double[,] array = new double[144, 3];
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        for (int m = 0; m < 3; m++)
                        {
                            for (int n = 0; n < 2; n++)
                            {
                                int sx = 2 * i - 1;
                                int sy = 2 * j - 1;
                                int sz = 2 * k - 1;

                                int s = ((sx * sy * sz * (2 * n - 1)) + 1) / 2;

                                int t = (24 * i) + (12 * j) + (6 * k) + (2 * m) + n;

                                array[(3 * t) + 0, 0] = sx * c * radius;
                                array[(3 * t) + 0, 1] = sy * c * radius;
                                array[(3 * t) + 0, 2] = sz * c * radius;

                                array[(3 * t) + 1 + s, 0] = sx * b * radius;
                                array[(3 * t) + 1 + s, 1] = sy * b * radius;
                                array[(3 * t) + 1 + s, 2] = sz * b * radius;
                                array[(3 * t) + 1 + s, m] = 0;

                                int p = (m + n + 1) % 3;
                                array[(3 * t) + 2 - s, 0] = sx * a * radius;
                                array[(3 * t) + 2 - s, 1] = sy * a * radius;
                                array[(3 * t) + 2 - s, 2] = sz * a * radius;
                                array[(3 * t) + 2 - s, m] = 0;
                                array[(3 * t) + 2 - s, p] = 0;
                            }
                        }
                    }
                }
            }

            return Matrix<double>.Build.DenseOfArray(array);
        }

        #endregion Polyhedral Surface Generation

        #region Surface Manipulation

        /// <summary>
        /// Flips the orientation of all triangles in a surface
        /// </summary>
        /// <param name="vertices">Input surface</param>
        /// <returns>Surface with the orientation of each triangle flipped</returns>
        public static Matrix<double> OrientationFlip(Matrix<double> vertices)
        {
            double[,] array = vertices.ToArray();
            int nT = array.GetLength(0) / 3;
            for (int t = 0; t < nT; t++)
            {
                // Swap rows 3*t and 3*t+1
                for (int j = 0; j < 3; j++)
                {
                    double temp = array[3 * t, j];
                    array[3 * t, j] = array[3 * t + 1, j];
                    array[3 * t + 1, j] = temp;
                }
            }
            return Matrix<double>.Build.DenseOfArray(array);
        }

        /// <summary>
        /// Normalizes a row of a matrix to have a certain magnitude as a vector.
        /// </summary>
        /// <param name="array">Array containing the data</param>
        /// <param name="rowIndex">Which row to normalize</param>
        /// <param name="scale">Scale factor to apply after normalizing</param>
        public static void NormalizeRow(double[,] array, int rowIndex, double scale)
        {
            double sf = scale / Math.Sqrt(
                (array[rowIndex, 0] * array[rowIndex, 0])
                + (array[rowIndex, 1] * array[rowIndex, 1])
                + (array[rowIndex, 2] * array[rowIndex, 2]));
            array[rowIndex, 0] *= sf;
            array[rowIndex, 1] *= sf;
            array[rowIndex, 2] *= sf;
        }

        #endregion Surface Manipulation

        #region Functions to mimic Matlab

        /// <summary>
        /// Implementation of Matlab's ndGrid for 2-dimensional data.
        /// </summary>
        /// <param name="x1">first input vector (1-dimensional)</param>
        /// <param name="x2">second input vector (1-dimensional)</param>
        /// <param name="y1">first output matrix, consists of replicas of x1 interpreted as columns</param>
        /// <param name="y2">second output matrix, consists of replicas of x2 interpreted as rows</param>
        public static void Grid(double[] x1, double[] x2, out Matrix<double> y1, out Matrix<double> y2)
        {
            int m = x1.Length;
            int n = x2.Length;
            double[,] array1 = new double[m, n];
            double[,] array2 = new double[m, n];
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    array1[i, j] = x1[i];
                    array2[i, j] = x2[j];
                }
            }
            y1 = Matrix<double>.Build.DenseOfArray(array1);
            y2 = Matrix<double>.Build.DenseOfArray(array2);
        }

        /// <summary>
        /// Meant to mimic the colon operator in Matlab, flattening a matrix column-by-column.
        /// </summary>
        /// <param name="input">Input matrix</param>
        /// <returns>Flattened matrix as an array of doubles</returns>
        public static double[] Flatten(Matrix<double> input)
        {
            double[,] array = input.ToArray();
            int m = array.GetLength(0);
            int n = array.GetLength(1);
            double[] output = new double[m * n];
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    output[m * j + i] = array[i, j];
                }
            }
            return output;
        }

        /// <summary>
        /// Meant to mimic the ability in Matlab to set a linearly spaced submatrix of one matrix to another matrix.
        /// </summary>
        /// <param name="destination">Destination matrix</param>
        /// <param name="iStart">Starting i index in the destination matrix</param>
        /// <param name="iStep">Step in i in the destination matrix</param>
        /// <param name="iCount">Total count of i-values to set in the destination matrix</param>
        /// <param name="jStart">Starting j index in the destination matrix</param>
        /// <param name="jStep">Step in j in the destination matrix</param>
        /// <param name="jCount">Total count of j-values to set in the destination matrix</param>
        /// <param name="source">Source matrix providing the data to copy to the submatrix of the destination matrix</param>
        public static void SetSubMatrix(Matrix<double> destination,
            int iStart, int iStep, int iCount,
            int jStart, int jStep, int jCount,
            Matrix<double> source)
        {
            for (int i = 0; i < iCount; i++)
            {
                for (int j = 0; j < jCount; j++)
                {
                    destination[iStart + (iStep * i), jStart + (jStep * j)] = source[i, j];
                }
            }
        }

        #endregion Functions to mimic Matlab

        /// <summary>
        /// Converts a matrix to a string suitable for writing to a (text) STL file.
        /// </summary>
        /// <param name="vertices">Matrix of vertices, 3 per triangle</param>
        /// <param name="name">Nane to give the solid in the resulting STL</param>
        /// <returns>String containing the contents to be written to an STL file</returns>
        public static string SurfaceToSTLData(Matrix<double> vertices, string name)
        {
            int nTriangles = vertices.RowCount / 3;
            StringBuilder output = new StringBuilder();
            output.AppendLine("solid " + name);
            for (int facet = 0; facet < nTriangles; facet++)
            {
                output.AppendLine("facet normal 0.0 0.0 1.0");
                output.AppendLine("    outer loop");
                output.AppendLine(string.Format("        vertex {0:F3} {1:F3} {2:F3}", vertices[3 * facet, 0], vertices[3 * facet, 1], vertices[3 * facet, 2]));
                output.AppendLine(string.Format("        vertex {0:F3} {1:F3} {2:F3}", vertices[3 * facet + 1, 0], vertices[3 * facet + 1, 1], vertices[3 * facet + 1, 2]));
                output.AppendLine(string.Format("        vertex {0:F3} {1:F3} {2:F3}", vertices[3 * facet + 2, 0], vertices[3 * facet + 2, 1], vertices[3 * facet + 2, 2]));
                output.AppendLine("    endloop");
                output.AppendLine("endfacet");
            }
            output.AppendLine("endsolid " + name);
            return output.ToString();
        }
    }
}
