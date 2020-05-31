using System;
using System.Collections.Generic;
using System.Linq;

namespace AdaptiveFEM.Services
{
    class LinealAproksimation
    {
        private Functions f = new Functions();

        private const double x1 = -1.0;
        private const double xn = 1;
        private double alpha = Math.Pow(10, 28);
        private double u0 = 0;

        public int N;

        public int prevN = 0;
        private int startN;
        private double startEn = 0.0;
        private double uNv = 0.0;
        private double eNv = 0.0;

        private double[] xi;
        private double[] new_xi;

        private double[] q;
        private List<double> allCeh = new List<double>();

        private double accuracy = 5;
        double[] qder;
        public LinealAproksimation(int _N)
        {
            N = _N;
            xi = new double[N + 1];
            startN = N;

            q = new double[N + 1];
            qder = new double[N];

            aproksimation();
            calcU();
        }

        private void aproksimation()
        {
            double h = (xn - x1) / (N);
            xi[0] = x1;
            for (int i = 1; i <= N; i++)
            {
                xi[i] = x1 + i * h;
            }
        }


        public double[,] GetMatrix()
        {
            double[,] matrix = new double[N + 1, N + 1];
            double h = xi[1] - xi[0];
            double hp = 0;

            matrix[0, 0] = f.calcMu(xi[0] + h) / h - f.calcBeta(xi[0] + h) / 2 + f.calcSigma(xi[0] + h) * h / 3 + alpha;
            matrix[0, 1] = -f.calcMu(xi[0] + h) / h + f.calcBeta(xi[0] + h) / 2 + f.calcSigma(xi[0] + h) * h / 6;

            for (int i = 1; i < N; ++i)
            {
                h = xi[i + 1] - xi[i];
                hp = xi[i] - xi[i - 1];

                matrix[i, i - 1] = -f.calcMu(xi[i - 1] + hp) / hp - f.calcBeta(xi[i - 1] + hp) / 2 + f.calcSigma(xi[i - 1] + hp) * hp / 6;

                matrix[i, i] = f.calcMu(xi[i - 1] + hp) / hp + f.calcBeta(xi[i - 1] + hp) / 2 + f.calcSigma(xi[i - 1] + hp) * hp / 3 +
                    f.calcMu(xi[i] + h) / h - f.calcBeta(xi[i] + h) / 2 + f.calcBeta(xi[i] + h) * h / 3;

                matrix[i, i + 1] = -f.calcMu(xi[i] + h) / h + f.calcBeta(xi[i] + h) / 2 + f.calcSigma(xi[i] + h) * h / 6;
            }

            hp = xi[N] - xi[N - 1];

            matrix[N, N - 1] = -f.calcMu(xi[N - 1] + hp) / hp - f.calcBeta(xi[N - 1] + hp) / 2 + f.calcSigma(xi[N - 1] + hp) * hp / 6;
            matrix[N, N] = f.calcMu(xi[N - 1] + hp) / hp + f.calcBeta(xi[N - 1] + hp) / 2 + f.calcSigma(xi[N - 1] + hp) * hp / 3 + alpha;
            return matrix;
        }

        public double[] GetVector()
        {
            double[] vector = new double[N + 1];

            double h = xi[1] - xi[0];
            //vector[0] = 0.5 * e[0].H * f.Evaluate(e[0]+ h) + alpha * u_a;
            vector[0] = 0.5 * h * f.calcF(xi[0] + h) + alpha * 0;//u_a;
            double hp = 0;
            for (int i = 1; i < N; ++i)
            {
                hp = xi[i] - xi[i - 1];
                h = xi[i + 1] - xi[i];
                vector[i] = 0.5 * (hp * f.calcF(xi[i - 1] + hp) + h * f.calcF(xi[i] + h));
                //vector[i] = 0.5 * (h * f.Evaluate(xi[i-1]+ h) + xi[i].H * f.Evaluate(xi[i]+ h));
            }
            //vector[N] = 0.5 * e[N - 1].H * f.Evaluate(e[N - 1]+ h) + gamma * u_b;
            hp = xi[N] - xi[N - 1];
            vector[N] = 0.5 * hp * f.calcF(xi[N - 1] + hp) + alpha * 0;// u_b;
            return vector;
        }



        private void calcU()
        {
            q = Accord.Math.Matrix.Solve(GetMatrix(), GetVector());

            for (int i = 1; i < q.Length; ++i)
            {
                qder[i - 1] = (q[i] - q[i - 1]) / (xi[i] - xi[i - 1]);
            }
        }

        private double calcE(int i)
        {
            double h = xi[i + 1] - xi[i];
            double m = (h * h * h) / f.calcMu(xi[i] + h);

            double b = (f.calcF(xi[i] + h)) - (f.calcBeta(xi[i] + h) * qder[i]) - (f.calcSigma(xi[i] + h) * q[i]);
            double d = 10 + (((h * f.calcBeta(xi[i] + h)) / f.calcMu(xi[i] + h)) * ((h * h * f.calcSigma(xi[i] + h)) / f.calcMu(xi[i] + h)));

            return Math.Abs((5.0 / 6.0) * m * (b * b / d));
        }

        private double mistakeIndicator(int j)
        {
            double uN = Accord.Math.Matrix.DotAndDot(q, GetMatrix(), q);
            double eN = 0.0;
            var result = 0.0;

            for (int i = 0; i < xi.Length - 1; ++i)
            {
                double h = xi[i + 1] - xi[i];
                double m = (h * h * h) / f.calcMu(xi[i] + h);

                double b = (f.calcF(xi[i] + h)) - (f.calcBeta(xi[i] + h) * qder[i]) - (f.calcSigma(xi[i] + h) * q[i]);
                double d = 10 + ((h * f.calcBeta(xi[i] + h)) / f.calcMu(xi[i] + h) * ((h * h * f.calcSigma(xi[i] + h)) / f.calcMu(xi[i] + h)));

                eN += Math.Abs((5.0 / 6.0) * m * (b * b / d));
            }
            double eNaEl = Math.Sqrt(calcE(j));
            uNv = Math.Sqrt(uN);
            eNv = Math.Sqrt(eN);
            result = (eNaEl * Math.Sqrt(N) * 100) / Math.Sqrt(uN + eN);
            if (N == startN)
            {
                startEn = eNv;
            }
            allCeh.Add(result);
            return result;
        }

        public void nextstep()
        {
            prevN = N;
            new_xi = xi;
            var newN = N;
            List<double> nx = new List<double>();

            var xil = new List<double>();
            for (int i = 0; i < N; i++)
            {
                if (mistakeIndicator(i) > accuracy)
                {
                    nx.Add((new_xi[i] + new_xi[i + 1]) / 2);
                    newN++;
                }
            }

            if (newN != N)
            {
                nx.AddRange(xi);
                nx.Sort();
                new_xi = nx.ToArray();
                N = newN;

                xi = new_xi;
                q = new double[N + 1];
                qder = new double[N + 1];
                calcU();
            }
        }


        public double[] getXi()
        {
            return xi;
        }

        public double[] getQ()
        {
            return q;
        }

        public Tuple<double[], double[]> getAllCeh()
        {
            var ceh = allCeh.ToArray();
            List<double> x = new List<double>();
            for (int i = 0; i < ceh.Length; i++)
            {
                x.Add(i);
            }

            var ret = Tuple.Create(ceh, x.ToArray());
            allCeh.Clear();
            return ret;
        }

        public Tuple<double, double> getNorms()
        {
            return Tuple.Create(eNv, uNv);
        }

        public double getGlobalMistake()
        {
            return (eNv / uNv) * 100;
        }

        public Tuple<double, double> getVars()
        {
            double p;
            if (N == startN)
            {
                p = 0;
            }
            else
            {
                p = (Math.Log(startEn) - Math.Log(eNv)) / (Math.Log(N) - Math.Log(startN));
            }

            return Tuple.Create(p, allCeh.Max());
        }
    }
}
