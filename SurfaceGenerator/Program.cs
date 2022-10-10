using MathNet.Numerics.LinearAlgebra;
using System;
using System.IO;
using System.Linq;
using static SurfaceGenerator.SurfaceGenerator;

namespace SurfaceGenerator
{
    /// <summary>
    /// Class containing the main entry point for the application.
    /// </summary>
    class Program
    {
        /// <summary>
        /// Main entry point for the application
        /// </summary>
        /// <param name="args">Command line arguments (if any)</param>
        static void Main(string[] args)
        {
            // Simple example to test that our functions are working as expected.
            const string filename = "wavy_cylinder.stl";

            const double cylinderRadius = 1.0;
            const double waveAmplitude = 0.3;
            const double thetaMin = 0;
            const double thetaMax = 1.5 * Math.PI;
            const int nTheta = 10;
            const double yMin = -1.0;
            const double yMax = 1.0;
            const int nY = 8;

            Func<double[], double[], Matrix<double>> wavyCylinder = 
                (theta, y) =>
                    Matrix<double>.Build.DenseOfColumnArrays(
                        theta.Select((t, i) => (cylinderRadius + waveAmplitude * Math.Cos(4 * y[i])) * System.Math.Cos(t)).ToArray(),
                        y,
                        theta.Select((t, i) => (cylinderRadius + waveAmplitude * Math.Cos(4 * y[i])) * System.Math.Sin(t)).ToArray());

            Matrix<double> cylinderSurface = GenerateQuadSurface(
                wavyCylinder,
                thetaMin, thetaMax, nTheta,
                yMin, yMax, nY);

            File.WriteAllText(filename, SurfaceToSTLData(cylinderSurface, "wavy_cylinder"));
        }
    }
}
