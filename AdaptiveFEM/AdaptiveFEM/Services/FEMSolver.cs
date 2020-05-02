using AdaptiveFEM.Models;
using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra.Double;

namespace AdaptiveFEM.Services
{
    public class FEMSolver
    {
        //Functions
        Function mu;
        Function beta;
        Function sigma;
        Function f;

        //Constants
        double a;
        double b;
        double alpha;
        double gamma;
        double u_a;
        double u_b;
        double error;
        int initialN;

        //Result
        List<Iteration> iterations;

        //Constructor with parameters
        public FEMSolver(Function mu, Function beta, Function sigma, Function f,
                             double a, double b, double alpha, double gamma, double u_a, double u_b, double error, int N)
        {
            this.mu = mu;
            this.beta = beta;
            this.sigma = sigma;
            this.f = f;
            this.a = a;
            this.b = b;
            this.alpha = alpha;
            this.gamma = gamma;
            this.u_a = u_a;
            this.u_b = u_b;
            this.error = error;
            this.initialN = N;

            iterations = new List<Iteration>();
        }

        public List<Iteration> Iterations { get { return this.iterations; } }
        public double A { get { return this.a; } }
        public double B { get { return this.b; } }

        public List<Element> InitialFiniteElements()
        {
            List<Element> e = new List<Element>();
            double step = Math.Abs(b - a) / initialN;
            double i = a;
            while (i < b)
            {
                e.Add(new Element(i, Math.Round(i + step, 10), Math.Round(i + step / 2, 10), step));
                i = Math.Round(i + step, 10);
            }
            return e;
        }

        public void InitialIteration()
        {
            //Set elements
            List<Element> e = InitialFiniteElements();

            //Solve system
            Matrix matrix = GetMatrix(e);
            Vector vector = GetVector(e);
            Vector solution = Solve(matrix, vector);

            //Solution in midPoint
            Vector solutionCenter = new DenseVector(initialN);
            Vector solutionCenterDeriv = new DenseVector(initialN);
            for (int i = 1; i < solution.Count; ++i)
            {
                solutionCenter[i - 1] = (solution[i] + solution[i - 1]) * 0.5;
                solutionCenterDeriv[i - 1] = (solution[i] - solution[i - 1]) / e[i - 1].H;
            }

            //Setting iteration data
            Iteration data = new Iteration();
            data.Elements = e;
            data.Solution = solution;
            data.SolutionCenter = solutionCenter;
            data.SolutionCenterDeriv = solutionCenterDeriv;

            //Calculating errors
            CalculateErrors(ref data);

            //Save iteration data
            iterations.Add(data);
        }

        public Matrix GetMatrix(List<Element> e)
        {
            int N = e.Count;
            //Matrix
            Matrix matrix = new DenseMatrix(N + 1, N + 1);
            matrix[0, 0] = mu.Evaluate(e[0].MidPoint) / e[0].H - beta.Evaluate(e[0].MidPoint) / 2 + sigma.Evaluate(e[0].MidPoint) * e[0].H / 3 + alpha;
            matrix[0, 1] = -mu.Evaluate(e[0].MidPoint) / e[0].H + beta.Evaluate(e[0].MidPoint) / 2 + sigma.Evaluate(e[0].MidPoint) * e[0].H / 6;
            for (int i = 1; i < N; ++i)
            {
                matrix[i, i - 1] = -mu.Evaluate(e[i - 1].MidPoint) / e[i - 1].H - beta.Evaluate(e[i - 1].MidPoint) / 2 + sigma.Evaluate(e[i - 1].MidPoint) * e[i - 1].H / 6;

                matrix[i, i] = mu.Evaluate(e[i - 1].MidPoint) / e[i - 1].H + beta.Evaluate(e[i - 1].MidPoint) / 2 + sigma.Evaluate(e[i - 1].MidPoint) * e[i - 1].H / 3 +
                    mu.Evaluate(e[i].MidPoint) / e[i].H - beta.Evaluate(e[i].MidPoint) / 2 + sigma.Evaluate(e[i].MidPoint) * e[i].H / 3;

                matrix[i, i + 1] = -mu.Evaluate(e[i].MidPoint) / e[i].H + beta.Evaluate(e[i].MidPoint) / 2 + sigma.Evaluate(e[i].MidPoint) * e[i].H / 6;
            }
            matrix[N, N - 1] = -mu.Evaluate(e[N - 1].MidPoint) / e[N - 1].H - beta.Evaluate(e[N - 1].MidPoint) / 2 + sigma.Evaluate(e[N - 1].MidPoint) * e[N - 1].H / 6;
            matrix[N, N] = mu.Evaluate(e[N - 1].MidPoint) / e[N - 1].H + beta.Evaluate(e[N - 1].MidPoint) / 2 + sigma.Evaluate(e[N - 1].MidPoint) * e[N - 1].H / 3 + gamma;
            return matrix;
        }

        public Vector GetVector(List<Element> e)
        {
            int N = e.Count;
            //Vector
            Vector vector = new DenseVector(N + 1);
            vector[0] = 0.5 * e[0].H * f.Evaluate(e[0].MidPoint) + alpha * u_a;
            for (int i = 1; i < N; ++i)
            {
                vector[i] = 0.5 * (e[i - 1].H * f.Evaluate(e[i - 1].MidPoint) + e[i].H * f.Evaluate(e[i].MidPoint));
            }
            vector[N] = 0.5 * e[N - 1].H * f.Evaluate(e[N - 1].MidPoint) + gamma * u_b;
            return vector;
        }

        public void Run()
        {
            //zero iteration
            InitialIteration();

            //Indexes of elements with the largest error
            List<int> indexes = new List<int>();

            //Setting iteration number
            int iterationNumber = 1;


            do
            {
                indexes.Clear();

                //Finding indexes of elements with error, larger then needed
                for (int i = 0; i < iterations[iterationNumber - 1].Elements.Count; ++i)
                {
                    if (iterations[iterationNumber - 1].Errors[i] > error)
                    {
                        indexes.Add(i);
                    }
                }

                if (indexes.Count != 0)
                {
                    //Divide elements with unappropriate error
                    List<Element> elements = new List<Element>();
                    for (int i = 0; i < iterations[iterationNumber - 1].Elements.Count; ++i)
                    {
                        if (!indexes.Contains(i))
                        {
                            elements.Add(iterations[iterationNumber - 1].Elements[i]);
                        }
                        else
                        {
                            Element e = iterations[iterationNumber - 1].Elements[i];
                            Element e1 = new Element(e.Begin, e.MidPoint, (e.Begin + e.MidPoint) / 2, e.MidPoint - e.Begin);
                            Element e2 = new Element(e.MidPoint, e.End, (e.MidPoint + e.End) / 2, e.End - e.MidPoint);
                            elements.Add(e1);
                            elements.Add(e2);
                        }
                    }

                    int N = elements.Count;
                    Matrix matrix = GetMatrix(elements);
                    Vector vector = GetVector(elements);
                    Vector solution = Solve(matrix, vector);
                    Vector solutionCenter = new DenseVector(N);
                    Vector solutionCenterDeriv = new DenseVector(N);
                    for (int i = 1; i < solution.Count; ++i)
                    {
                        solutionCenter[i - 1] = (solution[i] + solution[i - 1]) * 0.5;
                        solutionCenterDeriv[i - 1] = (solution[i] - solution[i - 1]) / elements[i - 1].H;
                    }

                    Iteration data = new Iteration();
                    data.Elements = elements;
                    data.Solution = solution;
                    data.SolutionCenter = solutionCenter;
                    data.SolutionCenterDeriv = solutionCenterDeriv;
                    CalculateErrors(ref data);
                    iterations.Add(data);
                    ++iterationNumber;
                }
            }
            while (indexes.Count != 0);

        }

        public void CalculateErrors(ref Iteration data)
        {
            int N = data.Elements.Count;
            Vector errors = new DenseVector(N);
            Vector errorsNorms = new DenseVector(N);

            double eNormV2 = 0;
            double uNormV2 = 0;
            for (int i = 0; i < N; ++i)
            {
                //eNormV2_i = m * (b * b / d)
                double m = data.Elements[i].H * data.Elements[i].H * data.Elements[i].H / mu.Evaluate(data.Elements[i].MidPoint);
                double b = f.Evaluate(data.Elements[i].MidPoint) - beta.Evaluate(data.Elements[i].MidPoint) * data.SolutionCenterDeriv[i] - sigma.Evaluate(data.Elements[i].MidPoint) * data.SolutionCenter[i];
                double d = 10 + data.Elements[i].H * data.Elements[i].H * sigma.Evaluate(data.Elements[i].MidPoint) / mu.Evaluate(data.Elements[i].MidPoint);
                double e_h2 = (5.0 / 6) * m * (b * b / d);
                //double a = (2.0 / 3) * data.Elements[i].H * data.Elements[i].H * (f.Evaluate(data.Elements[i].MidPoint) - beta.Evaluate(data.Elements[i].MidPoint) * data.SolutionCenterDeriv[i] - sigma.Evaluate(data.Elements[i].MidPoint) * data.SolutionCenter[i]);
                //double b = (16.0 / 3) * mu.Evaluate(data.Elements[i].MidPoint) / data.Elements[i].H + (8.0 / 15) * sigma.Evaluate(data.Elements[i].MidPoint) * data.Elements[i].H;
                //double e_h2 = a*a/b;
                errorsNorms[i] = Math.Sqrt(Math.Abs(e_h2));
                eNormV2 += Math.Abs(e_h2);

                double q = data.SolutionCenterDeriv[i];
                //uNormV2_i
                uNormV2 += data.Elements[i].H * q * q;
            }

            data.SolutionNormV = uNormV2;
            data.ErrorNormV = eNormV2;
            data.ErrorsNormsV = errorsNorms;

            double sqrtN = Math.Sqrt(N);
            double sqrtNorms = Math.Sqrt(uNormV2 + eNormV2);

            for (int i = 0; i < N; ++i)
            {
                errors[i] = (errorsNorms[i] * sqrtN * 100) / sqrtNorms;
            }
            data.Errors = errors;
        }

        public Vector Solve(Matrix matrix, Vector vector)
        {
            if (matrix.ColumnCount != matrix.RowCount || matrix.ColumnCount != vector.Count)
            {
                throw new Exception("Sizes of matrix and vector are not the same");
            }

            return (Vector)matrix.Solve(vector);
        }
    }
}
