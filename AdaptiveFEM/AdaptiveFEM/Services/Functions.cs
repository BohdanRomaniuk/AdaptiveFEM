using System;

namespace AdaptiveFEM.Services
{
    class Functions
    {
        public double calcBeta(double x = 0)
        {
            return 1500.0 * Math.Pow(x, 8); 
        }

        public double calcSigma(double x = 0)
        {
            return 80.0 + 2 * x * x;//0.0;
        }

        public double calcMu(double x = 0)
        {
            return 1.0;
        }

        public double calcF(double x = 0)
        {
            return 100 * Math.Exp(Math.Pow(x - 0.15, 7)) * x;
        }
    }
}
