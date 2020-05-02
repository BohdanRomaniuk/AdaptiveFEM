using MathNet.Numerics.LinearAlgebra.Double;
using System.Collections.Generic;

namespace AdaptiveFEM.Models
{
    public class Iteration
    {
        public List<Element> Elements { get; set; }

        public Vector SolutionCenter { get; set; }

        public Vector SolutionCenterDeriv { get; set; }

        public Vector Errors { get; set; }

        public Vector Solution { get; set; }

        public Vector ErrorsNormsV { get; set; }

        public double SolutionNormV { get; set; }

        public double ErrorNormV { get; set; }
    }
}
