using NCalc;

namespace AdaptiveFEM.Helpers
{
    public static class MathHelper
    {
        public static double Evaluate(this Expression expression, double x)
        {
            expression.Parameters["X"] = x;
            return (double)expression.Evaluate();
        }
    }
}
