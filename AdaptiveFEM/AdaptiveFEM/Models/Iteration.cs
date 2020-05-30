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

        public double SolutionNormV { get; set; }

        public double ErrorNormV { get; set; }

        public double uNv { get; set; }

        public double eNv { get; set; }

        public double OrderOfConvergence { get; set; }

        public double MaxRelativeError { get; set; }

        public double GlobalError => (uNv != 0.0) ? (eNv / uNv) * 100 : 0;
    }
}
