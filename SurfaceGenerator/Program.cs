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
            File.WriteAllText("skeleton_long_piece.stl",
                SurfaceToSTLData(
                    GenerateSkeletonLongPiece(),
                    "skeleton_long_piece"));

            File.WriteAllText("skeleton_short_piece.stl",
                SurfaceToSTLData(
                    GenerateSkeletonShortPiece(),
                    "skeleton_short_piece"));
        }

        /// <summary>
        /// Generates the long piece for the skeleton joints.
        /// </summary>
        /// <returns>Matrix of vertex coordinates, which can be used to generate STL data</returns>
        public static Matrix<double> GenerateSkeletonLongPiece()
        {
            const double innerGearRadius = 7.0;
            const double outerRadius = 10.0;
            const double thetaStart = 0;
            const double thetaEnd = 2 * Math.PI;
            const int nTeeth = 15;
            const int nTheta = 8 * nTeeth;
            const int nR = 5;
            const double totalHeight = 28.0;
            const double gearBaseHeight = 0.5;
            const double baseThickness = 3.5;
            const double innerHoleRadius = 3.0;
            const double outerHoleRadius = 4.0;
            const double endThickness = 3.0;
            const double boreLength = 10.5;
            const double wallThickness = 1.5;

            Matrix<double> teeth = GenerateQuadSurface(
                (r, theta) =>
                    Matrix<double>.Build.DenseOfColumnArrays(
                        theta.Select((t, i) => r[i] * Math.Cos(t)).ToArray(),
                        theta.Select((t, i) => r[i] * Math.Sin(t)).ToArray(),
                        theta.Select((t, i) => totalHeight - (outerRadius * 0.5 * Math.PI / nTeeth) + r[i] * (0.5 * Math.PI / nTeeth) * SawTooth(nTeeth * t)).ToArray()),
                outerRadius, innerGearRadius, nR,
                thetaStart, thetaEnd, nTheta);

            Matrix<double> teethInnerWall = GenerateQuadSurface(
                (z, theta) =>
                    Matrix<double>.Build.DenseOfColumnArrays(
                        theta.Select(t => innerGearRadius * Math.Cos(t)).ToArray(),
                        theta.Select(t => innerGearRadius * Math.Sin(t)).ToArray(),
                        theta.Select((t, i) => totalHeight - (outerRadius * 0.5 * Math.PI / nTeeth) - (innerGearRadius * 0.5 * Math.PI / nTeeth + gearBaseHeight)
                            + z[i] * ((innerGearRadius * 0.5 * Math.PI / nTeeth + gearBaseHeight) + innerGearRadius * (0.5 * Math.PI / nTeeth) * SawTooth(nTeeth * t))).ToArray()),
                1, 0, 1,
                thetaStart, thetaEnd, nTheta);

            Matrix<double> outerWall = GenerateQuadSurface(
               (z, theta) =>
                   Matrix<double>.Build.DenseOfColumnArrays(
                       theta.Select(t => outerRadius * Math.Cos(t)).ToArray(),
                       theta.Select(t => outerRadius * Math.Sin(t)).ToArray(),
                       theta.Select((t, i) => 
                             (1 - z[i]) * baseThickness
                           + z[i] * (totalHeight - outerRadius * 0.5 * Math.PI / nTeeth + outerRadius * (0.5 * Math.PI / nTeeth) * SawTooth(nTeeth * t))).ToArray()),
               0, 1, 1,
               thetaStart, thetaEnd, nTheta);

            double zBase = totalHeight - (outerRadius * 0.5 * Math.PI / nTeeth) - (innerGearRadius * 0.5 * Math.PI / nTeeth + gearBaseHeight);
            Matrix<double> gearBase = GenerateQuadSurface(
                (r, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => r[i] * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => r[i] * Math.Sin(t)).ToArray(),
                    theta.Select(t => zBase).ToArray()),
                innerGearRadius, innerHoleRadius, 1,
                thetaStart, thetaEnd, nTheta);

            double zBoreEnd = zBase - endThickness - boreLength;
            Matrix<double> innerBore = GenerateQuadSurface(
                (z, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => innerHoleRadius * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => innerHoleRadius * Math.Sin(t)).ToArray(),
                    z),
                zBase, zBoreEnd, 1,
                thetaStart, thetaEnd, nTheta);

            Matrix<double> boreEndFlat = GenerateQuadSurface(
                (r, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => r[i] * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => r[i] * Math.Sin(t)).ToArray(),
                    theta.Select(t => zBoreEnd).ToArray()),
                innerHoleRadius, outerHoleRadius, 1,
                thetaStart, thetaEnd, nTheta);

            Matrix<double> outerBore = GenerateQuadSurface(
                (z, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => outerHoleRadius * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => outerHoleRadius * Math.Sin(t)).ToArray(),
                    z),
                zBoreEnd, zBase - endThickness, 1,
                thetaStart, thetaEnd, nTheta);

            Matrix<double> boreBaseFlat = GenerateQuadSurface(
                (r, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => r[i] * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => r[i] * Math.Sin(t)).ToArray(),
                    theta.Select(t => zBase - endThickness).ToArray()),
                outerHoleRadius, outerRadius - wallThickness, 1,
                thetaStart, thetaEnd, nTheta);

            Matrix<double> innerWall = GenerateQuadSurface(
                (z, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => (outerRadius - wallThickness) * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => (outerRadius - wallThickness) * Math.Sin(t)).ToArray(),
                    z),
                zBase - endThickness, 0, 1,
                thetaStart, thetaEnd, nTheta);

            Matrix<double> bottom = GenerateQuadSurface(
                (r, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => r[i] * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => r[i] * Math.Sin(t)).ToArray(),
                    theta.Select(t => 0.0).ToArray()),
                outerRadius - wallThickness, outerRadius + 10, 1,
                thetaStart, thetaEnd, nTheta);

            Matrix<double> topOfBottom = GenerateQuadSurface(
                (r, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => r[i] * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => r[i] * Math.Sin(t)).ToArray(),
                    theta.Select(t => baseThickness).ToArray()),
                outerRadius + 10, outerRadius, 1,
                thetaStart, thetaEnd, nTheta);

            Matrix<double> outerEdgeOfBottom = GenerateQuadSurface(
                (z, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => (outerRadius + 10) * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => (outerRadius + 10) * Math.Sin(t)).ToArray(),
                    z),
                0, baseThickness, 1,
                thetaStart, thetaEnd, nTheta);

            return Matrix<double>.Build.DenseOfMatrixArray(new Matrix<double>[,] {
                { teeth }, { teethInnerWall }, { outerWall }, { gearBase }, { innerBore }, { boreEndFlat },
                { outerBore }, { boreBaseFlat }, { innerWall }, { bottom }, { topOfBottom }, { outerEdgeOfBottom } });

            // TODO: After generating this STL, construct a rectangular solid in Blender that is 26 mm in x by 20 mm in y by 30 mm in z (centered on z = 15 mm),
            // then intersect this surface with that rectangular solid.

            // For the simple washer piece, this can be made entirely in Blender:
            // Cylinder 6.5 mm height, 4 mm radius, z-offset 3.25
            // Union with : cylinder 2.5 mm height, 8.3 mm radius, z-offset: 1.25
            // Subtract with: cylinder >6.5 mm height, 3 mm radius
        }

        /// <summary>
        /// Generates the long piece for the skeleton joints.
        /// </summary>
        /// <returns>Matrix of vertex coordinates, which can be used to generate STL data</returns>
        public static Matrix<double> GenerateSkeletonShortPiece()
        {
            const double innerGearRadius = 7.0;
            const double outerRadius = 10.0;
            const double thetaStart = 0;
            const double thetaEnd = 2 * Math.PI;
            const int nTeeth = 15;
            const int nTheta = 8 * nTeeth;
            const int nR = 5;
            const double totalHeight = 7.0;
            const double gearBaseHeight = 0.5;
            const double innerHoleRadius = 3.0;
            const double endThickness = 3.0;
            const double wallThickness = 1.5;

            Matrix<double> teeth = GenerateQuadSurface(
                (r, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => r[i] * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => r[i] * Math.Sin(t)).ToArray(),
                    theta.Select((t, i) => totalHeight - (outerRadius * 0.5 * Math.PI / nTeeth) + r[i] * (0.5 * Math.PI / nTeeth) * SawTooth(nTeeth * t)).ToArray()),
                outerRadius, innerGearRadius, nR,
                thetaStart, thetaEnd, nTheta);

            Matrix<double> teethInnerWall = GenerateQuadSurface(
                (z, theta) =>
                    Matrix<double>.Build.DenseOfColumnArrays(
                        theta.Select(t => innerGearRadius * Math.Cos(t)).ToArray(),
                        theta.Select(t => innerGearRadius * Math.Sin(t)).ToArray(),
                        theta.Select((t, i) => totalHeight - (outerRadius * 0.5 * Math.PI / nTeeth) - (innerGearRadius * 0.5 * Math.PI / nTeeth + gearBaseHeight)
                            + z[i] * ((innerGearRadius * 0.5 * Math.PI / nTeeth + gearBaseHeight) + innerGearRadius * (0.5 * Math.PI / nTeeth) * SawTooth(nTeeth * t))).ToArray()),
                1, 0, 1,
                thetaStart, thetaEnd, nTheta);

            Matrix<double> outerWall = GenerateQuadSurface(
               (z, theta) =>
                   Matrix<double>.Build.DenseOfColumnArrays(
                       theta.Select(t => outerRadius * Math.Cos(t)).ToArray(),
                       theta.Select(t => outerRadius * Math.Sin(t)).ToArray(),
                       theta.Select((t, i) =>
                           z[i] * (totalHeight - outerRadius * 0.5 * Math.PI / nTeeth + outerRadius * (0.5 * Math.PI / nTeeth) * SawTooth(nTeeth * t))).ToArray()),
               0, 1, 1,
               thetaStart, thetaEnd, nTheta);

            double zBase = totalHeight - (outerRadius * 0.5 * Math.PI / nTeeth) - (innerGearRadius * 0.5 * Math.PI / nTeeth + gearBaseHeight);
            Matrix<double> gearBase = GenerateQuadSurface(
                (r, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => r[i] * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => r[i] * Math.Sin(t)).ToArray(),
                    theta.Select(t => zBase).ToArray()),
                innerGearRadius, innerHoleRadius, 1,
                thetaStart, thetaEnd, nTheta);

            double zBoreEnd = zBase - endThickness;
            Matrix<double> innerBore = GenerateQuadSurface(
                (z, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => innerHoleRadius * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => innerHoleRadius * Math.Sin(t)).ToArray(),
                    z),
                zBase, zBoreEnd, 1,
                thetaStart, thetaEnd, nTheta);

            Matrix<double> boreBaseFlat = GenerateQuadSurface(
                (r, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => r[i] * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => r[i] * Math.Sin(t)).ToArray(),
                    theta.Select(t => zBoreEnd).ToArray()),
                innerHoleRadius, outerRadius - wallThickness, 1,
                thetaStart, thetaEnd, nTheta);

            Matrix<double> innerWall = GenerateQuadSurface(
                (z, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => (outerRadius - wallThickness) * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => (outerRadius - wallThickness) * Math.Sin(t)).ToArray(),
                    z),
                zBase - endThickness, 0, 1,
                thetaStart, thetaEnd, nTheta);

            Matrix<double> bottom = GenerateQuadSurface(
                (r, theta) =>
                Matrix<double>.Build.DenseOfColumnArrays(
                    theta.Select((t, i) => r[i] * Math.Cos(t)).ToArray(),
                    theta.Select((t, i) => r[i] * Math.Sin(t)).ToArray(),
                    theta.Select(t => 0.0).ToArray()),
                outerRadius - wallThickness, outerRadius, 1,
                thetaStart, thetaEnd, nTheta);

            Console.WriteLine("zBoreEnd = " + zBoreEnd);
            Console.WriteLine("zBoreEnd - 3.5 = " + (zBoreEnd - 3.5).ToString());

            return Matrix<double>.Build.DenseOfMatrixArray(new Matrix<double>[,] {
                { teeth }, { teethInnerWall }, { outerWall }, { gearBase }, { innerBore },
                { boreBaseFlat }, { innerWall }, { bottom } });

            // TODO: After generating this STL, construct four cylinders in Blender that are 5 mm diameter by 7 mm tall,
            // with z-offset: zBoreEnd - 3.5, x/y offsets +/- 5 mm. Union with these cylinders.
        }

        /// <summary>
        /// Sawtooth wave function
        /// </summary>
        /// <param name="theta">Input theta</param>
        /// <returns>Value of the sawtooth function at theta</returns>
        public static double SawTooth(double theta)
        {
            return Math.Asin(Math.Sin(theta)) * 2 / Math.PI;
        }
    }
}
