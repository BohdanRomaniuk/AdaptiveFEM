using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra.Double;

namespace AdaptiveFEM
{
    public class IterationData
    {
        //N
        public List<Element> Elements { get; set; }
        //N
        public Vector SolutionCenter { get; set; }
        //N
        public Vector SolutionCenterDeriv { get; set; }
        //N
        public Vector Errors { get; set; }
        //N+1
        public Vector Solution { get; set; }

        public Vector ErrorsNormsV { get; set; }

        public double SolutionNormV { get; set; }

        public double ErrorNormV { get; set; }
    }
}
