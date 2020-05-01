using PoohMathParser;
using System;
using System.Collections.Generic;
using System.Text;

namespace AdaptiveFEM.Models
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
            return this.function;
        }
    }
}
