using PoohMathParser;

namespace AdaptiveFEM.Models
{
    public class Function
    {
        private string function;

        public Function()
        {
            function = "0";
        }

        public Function(string newFunction)
        {
            function = newFunction;
        }

        public double Evaluate(double x)
        {
            MathExpression expr;
            if (function[0] == '-')
            {
                expr = new MathExpression(function.Substring(1));
                return -1 * expr.Calculate(x);
            }
            else
            {
                expr = new MathExpression(function);
                return expr.Calculate(x);
            }
        }

        public override string ToString()
        {
            return function;
        }
    }
}
