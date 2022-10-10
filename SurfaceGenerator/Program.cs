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
            const string filename = "wavy_sphere.stl";

            const double sphereRadius = 1.0;
            const double waveAmplitude = 0.1;
            const double thetaStart = 0;
            const double thetaEnd = 1.5 * Math.PI;
            const int nTheta = 10;
            const double phiStart = 0;
            const double phiEnd = Math.PI;
            const int nPhi = 32;

            Func<double[], double[], Matrix<double>> wavySphere = 
                (theta, phi) =>
                    Matrix<double>.Build.DenseOfColumnArrays(
                        theta.Select((t, i) => (sphereRadius + waveAmplitude * Math.Sin(8 * phi[i])) * Math.Cos(t) * Math.Sin(phi[i])).ToArray(),
                        phi.Select(p => -(sphereRadius + waveAmplitude * Math.Sin(8 * p)) * Math.Cos(p)).ToArray(),
                        theta.Select((t, i) => (sphereRadius + waveAmplitude * Math.Sin(8 * phi[i])) * Math.Sin(t) * Math.Sin(phi[i])).ToArray());

            Matrix<double> sphereSurface = GenerateCrescentSurface(
                wavySphere,
                thetaStart, thetaEnd, nTheta,
                phiStart, phiEnd, nPhi);

            File.WriteAllText(filename, SurfaceToSTLData(sphereSurface, "wavy_sphere"));
        }
    }
}
