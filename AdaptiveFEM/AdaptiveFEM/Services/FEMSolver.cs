using AdaptiveFEM.Models;
using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra.Double;
using AdaptiveFEM.Services.Interfaces;
using NCalc;
using AdaptiveFEM.Helpers;

namespace AdaptiveFEM.Services
{
    public class FEMSolver : IFEMSolver
    {
        public Expression Mu { get; set; }
        public Expression Beta { get; set; }
        public Expression Sigma { get; set; }
        public Expression F { get; set; }

        public double A { get; set; }
        public double B { get; set; }
        public double Alpha { get; set; }
        public double Gamma { get; set; }
        public double Ua { get; set; }
        public double Ub { get; set; }
        public double Error { get; set; }
        public int InitialN { get; set; }
        public List<Iteration> Iterations { get; set; }
        private double startEn = 0.0;

        public FEMSolver(Expression mu, Expression beta, Expression sigma, Expression f, double a, double b, double alpha, double gamma, double ua, double ub, double error, int n)
        {
            Iterations = new List<Iteration>();
            Mu = mu;
            Beta = beta;
            Sigma = sigma;
            F = f;
            A = a;
            B = b;
            Alpha = alpha;
            Gamma = gamma;
            Ua = ua;
            Ub = ub;
            Error = error;
            InitialN = n;
        }

        public void Solve()
        {
            FirstIteration();
            //Indexes of elements with the largest Error
            var indexes = new List<int>();
            int iterationNumber = 1;

            do
            {
                indexes.Clear();
                //Finding indexes of elements with Error, larger then needed
                for (int i = 0; i < Iterations[iterationNumber - 1].Elements.Count; ++i)
                {
                    if (Iterations[iterationNumber - 1].Errors[i] > Error)
                    {
                        indexes.Add(i);
                    }
                }

                if (indexes.Count == 0)
                {
                    continue;
                }
                //Divide elements with unappropriate Error
                var elements = new List<Element>();
                for (int i = 0; i < Iterations[iterationNumber - 1].Elements.Count; ++i)
                {
                    if (!indexes.Contains(i))
                    {
                        elements.Add(Iterations[iterationNumber - 1].Elements[i]);
                    }
                    else
                    {
                        var elem = Iterations[iterationNumber - 1].Elements[i];
                        var leftElem = new Element(elem.Begin, elem.MidPoint, (elem.Begin + elem.MidPoint) / 2, elem.MidPoint - elem.Begin);
                        var rightElem = new Element(elem.MidPoint, elem.End, (elem.MidPoint + elem.End) / 2, elem.End - elem.MidPoint);
                        elements.Add(leftElem);
                        elements.Add(rightElem);
                    }
                }

                var n = elements.Count;
                var matrix = GenereteMatrix(elements);
                var vector = GenereteVector(elements);
                var solution = Solve(matrix, vector);
                var solutionCenter = new DenseVector(n);
                var solutionCenterDeriv = new DenseVector(n);
                for (int i = 1; i < solution.Count; ++i)
                {
                    solutionCenter[i - 1] = (solution[i] + solution[i - 1]) * 0.5;
                    solutionCenterDeriv[i - 1] = (solution[i] - solution[i - 1]) / elements[i - 1].H;
                }

                var iteration = new Iteration();
                iteration.Elements = elements;
                iteration.Solution = solution;
                iteration.SolutionCenter = solutionCenter;
                iteration.SolutionCenterDeriv = solutionCenterDeriv;
                CalculateErrors(ref iteration);
                Iterations.Add(iteration);
                ++iterationNumber;
            }
            while (indexes.Count != 0);
        }

        private void FirstIteration()
        {
            var elements = FirstIterationElements();
            var matrix = GenereteMatrix(elements);
            var vector = GenereteVector(elements);
            var solution = Solve(matrix, vector);

            //Solution in midPoint
            var solutionCenter = new DenseVector(InitialN);
            var solutionCenterDeriv = new DenseVector(InitialN);
            for (int i = 1; i < solution.Count; ++i)
            {
                solutionCenter[i - 1] = (solution[i] + solution[i - 1]) * 0.5;
                solutionCenterDeriv[i - 1] = (solution[i] - solution[i - 1]) / elements[i - 1].H;
            }

            var iteration = new Iteration();
            iteration.Elements = elements;
            iteration.Solution = solution;
            iteration.SolutionCenter = solutionCenter;
            iteration.SolutionCenterDeriv = solutionCenterDeriv;
            CalculateErrors(ref iteration);
            Iterations.Add(iteration);
        }

        private List<Element> FirstIterationElements()
        {
            var elements = new List<Element>();
            var step = Math.Abs(B - A) / InitialN;
            var xi = A;
            while (xi < B)
            {
                elements.Add(new Element(xi, Math.Round(xi + step, 10), Math.Round(xi + step / 2, 10), step));
                xi = Math.Round(xi + step, 10);
            }
            return elements;
        }

        private Matrix GenereteMatrix(List<Element> elem)
        {
            var n = elem.Count;
            var matrix = new DenseMatrix(n + 1, n + 1);
            matrix[0, 0] = Mu.Evaluate(elem[0].MidPoint) / elem[0].H - Beta.Evaluate(elem[0].MidPoint) / 2 + Sigma.Evaluate(elem[0].MidPoint) * elem[0].H / 3 + Alpha;
            matrix[0, 1] = -Mu.Evaluate(elem[0].MidPoint) / elem[0].H + Beta.Evaluate(elem[0].MidPoint) / 2 + Sigma.Evaluate(elem[0].MidPoint) * elem[0].H / 6;
            for (int i = 1; i < n; ++i)
            {
                matrix[i, i - 1] = -Mu.Evaluate(elem[i - 1].MidPoint) / elem[i - 1].H - Beta.Evaluate(elem[i - 1].MidPoint) / 2 + Sigma.Evaluate(elem[i - 1].MidPoint) * elem[i - 1].H / 6;
                matrix[i, i] = Mu.Evaluate(elem[i - 1].MidPoint) / elem[i - 1].H + Beta.Evaluate(elem[i - 1].MidPoint) / 2 + Sigma.Evaluate(elem[i - 1].MidPoint) * elem[i - 1].H / 3 + Mu.Evaluate(elem[i].MidPoint) / elem[i].H - Beta.Evaluate(elem[i].MidPoint) / 2 + Sigma.Evaluate(elem[i].MidPoint) * elem[i].H / 3;
                matrix[i, i + 1] = -Mu.Evaluate(elem[i].MidPoint) / elem[i].H + Beta.Evaluate(elem[i].MidPoint) / 2 + Sigma.Evaluate(elem[i].MidPoint) * elem[i].H / 6;
            }
            matrix[n, n - 1] = -Mu.Evaluate(elem[n - 1].MidPoint) / elem[n - 1].H - Beta.Evaluate(elem[n - 1].MidPoint) / 2 + Sigma.Evaluate(elem[n - 1].MidPoint) * elem[n - 1].H / 6;
            matrix[n, n] = Mu.Evaluate(elem[n - 1].MidPoint) / elem[n - 1].H + Beta.Evaluate(elem[n - 1].MidPoint) / 2 + Sigma.Evaluate(elem[n - 1].MidPoint) * elem[n - 1].H / 3 + Gamma;
            
            return matrix;
        }

        private Vector GenereteVector(List<Element> elem)
        {
            var n = elem.Count;
            Vector vector = new DenseVector(n + 1);
            vector[0] = 0.5 * elem[0].H * F.Evaluate(elem[0].MidPoint) + Alpha * Ua;
            for (int i = 1; i < n; ++i)
            {
                vector[i] = 0.5 * (elem[i - 1].H * F.Evaluate(elem[i - 1].MidPoint) + elem[i].H * F.Evaluate(elem[i].MidPoint));
            }
            vector[n] = 0.5 * elem[n - 1].H * F.Evaluate(elem[n - 1].MidPoint) + Gamma * Ub;
            
            return vector;
        }

        private void CalculateErrors(ref Iteration iter)
        {
            var n = iter.Elements.Count;
            var errors = new DenseVector(n);
            var errorsNorms = new DenseVector(n);

            var eNormV2 = 0.0;
            var uNormV2 = 0.0;
            for (int i = 0; i < n; ++i)
            {
                var m = iter.Elements[i].H * iter.Elements[i].H * iter.Elements[i].H / Mu.Evaluate(iter.Elements[i].Begin);
                var B = F.Evaluate(iter.Elements[i].MidPoint) - Beta.Evaluate(iter.Elements[i].MidPoint) * iter.SolutionCenterDeriv[i] - Sigma.Evaluate(iter.Elements[i].MidPoint) * iter.SolutionCenter[i];
                //var d = 10 + ((iter.Elements[i].H * Beta.Evaluate(iter.Elements[i].MidPoint)) / Mu.Evaluate(iter.Elements[i].MidPoint))*((iter.Elements[i].H * iter.Elements[i].H * Sigma.Evaluate(iter.Elements[i].MidPoint)) / Mu.Evaluate(iter.Elements[i].MidPoint));
                var d = 10 + ((iter.Elements[i].H * Beta.Evaluate(iter.Elements[i].MidPoint)) / Mu.Evaluate(iter.Elements[i].MidPoint) * ((iter.Elements[i].H * iter.Elements[i].H * Sigma.Evaluate(iter.Elements[i].MidPoint)) / Mu.Evaluate(iter.Elements[i].MidPoint)));
                var e_h2 = (5.0 / 6) * m * (B * B / d);
                errorsNorms[i] = Math.Sqrt(Math.Abs(e_h2));
                eNormV2 += Math.Abs(e_h2);
                //var qa = iter.SolutionCenter[i];
                //eNormV2 += iter.Elements[i].H * qa * qa;

                var q = iter.SolutionCenterDeriv[i];
                uNormV2 += iter.Elements[i].H * q * q;
            }

            var matrxi = GenereteMatrix(iter.Elements);
            double uN = iter.Solution * matrxi * iter.Solution;

            iter.SolutionNormV = uN;
            iter.ErrorNormV = eNormV2;
            iter.ErrorsNormsV = errorsNorms;

            iter.uNv = Math.Sqrt(uN);
            iter.eNv = Math.Sqrt(eNormV2);

            for (int i = 0; i < n; ++i)
            {
                errors[i] = (errorsNorms[i] * Math.Sqrt(n) * 100) / Math.Sqrt(uN + eNormV2);
            }
            iter.Errors = errors;
            if(iter.N == InitialN)
            {
                startEn = iter.eNv;
            }

            iter.OrderOfConvergence = (Math.Log(startEn) - Math.Log(iter.eNv)) / (Math.Log(iter.N) - Math.Log(iter.uNv));
            iter.MaxRelativeError = iter.Errors.Maximum();
        }

        private Vector Solve(Matrix matrix, Vector vector)
        {
            if (matrix.ColumnCount != matrix.RowCount || matrix.ColumnCount != vector.Count)
            {
                System.Windows.MessageBox.Show("Sizes of matrix and vector are not the same");
                return null;
            }

            return (Vector)matrix.Solve(vector);
        }
    }
}
