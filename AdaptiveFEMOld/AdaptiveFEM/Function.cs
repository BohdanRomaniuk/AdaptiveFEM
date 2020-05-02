using System;
//using NCalc;
using PoohMathParser;
using org.mariuszgromada.math.mxparser;

namespace AdaptiveFEM
{
    public class Function
    {
        private string function;

        public Function()
        {
            this.function = "0";
        }

        public Function(string newFunction)
        {
            this.function = newFunction;
        }

        public double Evaluate(double x)
        {
            //Expression e = new Expression(function);
            //e.Parameters["x"] = x;
            //return Convert.ToDouble(e.Evaluate());
            //switch (function)
            //{
            //    case "0": return 0;
            //    case "1": return 1;
            //    case "100": return 100;
            //    case "-100": return -100;
            //    default: return 0;
            //}

            //pooh
            MathExpression expr;
            if (function[0] == '-')
            {
                expr = new MathExpression(function.Substring(1));
                return -1 * expr.Calculate(x);
            }
            else
            {
                expr = new MathExpression(function);
               return  expr.Calculate(x);
            }

            //Function f = new Function("f(x) = " + function);
            //Expression e = new Expression("f(" + x.ToString() + ")", f);
            //return e.calculate();
        }

        public override string ToString()
        {
            return this.function;
        }
    }
}