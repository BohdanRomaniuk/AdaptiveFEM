using MathNet.Numerics.LinearAlgebra.Double;
using System.Collections.Generic;

namespace AdaptiveFEM.Models
{
    public class Iteration
    {
        public List<Element> Elements { get; set; }

        public int N => Elements?.Count ?? 0;

        public Vector SolutionCenter { get; set; }

        public Vector SolutionCenterDeriv { get; set; }

        public Vector Errors { get; set; }

        public Vector Solution { get; set; }

        public Vector ErrorsNormsV { get; set; }

        public double UNormV2 { get; set; }

        public double ENormV2 { get; set; }

        public double UNorm { get; set; }

        public double ENorm { get; set; }

        public double OrderOfConvergence { get; set; }

        public double MaxRelativeError { get; set; }

        public double GlobalError => (UNorm != 0.0) ? (ENorm / UNorm) * 100 : 0;
    }
}
