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
        /// <summary>
        /// Generates the vertices of a quad surface, 3 vertices per triangle.
        /// </summary>
        /// <param name="f">Function describing the (x,y,z) position of the surface as a function of two input parameters, t and u.</param>
        /// <param name="tMin">Minimum value of t</param>
        /// <param name="tMax">Maximum value of t</param>
        /// <param name="nT">Number of quadrilaterals to generate in the t-direction</param>
        /// <param name="uMin">Minimum value of u</param>
        /// <param name="uMax">Maximum value of u</param>
        /// <param name="nU">Number of quadrilaterals to generate in the u-direction</param>
        /// <returns>Matrix of vertex coordinates, which can be used to generate STL data</returns>
        public static Matrix<double> GenerateQuadSurface(
            Func<double[], double[], Matrix<double>> f,
            double tMin, double tMax, int nT,
            double uMin, double uMax, int nU)
        {
            List<double> tParams = Generate.LinearSpaced(nT + 1, tMin, tMax).ToList();
            List<double> tParamsShifted = new List<double>(tParams);
            tParamsShifted.RemoveAt(0);
            tParams.RemoveAt(tParams.Count - 1);

            List<double> uParams = Generate.LinearSpaced(nU + 1, uMin, uMax).ToList();
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
