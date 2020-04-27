using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using PoohMathParser;

namespace Tests
{
    [TestClass]
    public class Tests
    {
        MathExpression expr;

        [TestMethod]
        public void PlusTest1()
        {
            expr = new MathExpression("3+2");
            Assert.AreEqual(5, expr.Calculate());
        }

        [TestMethod]
        public void PlusTest2()
        {
            expr = new MathExpression("x+2");
            Assert.AreEqual(5, expr.Calculate(3));
        }

        [TestMethod]
        public void MinusTest1()
        {
            expr = new MathExpression("3-2");
            Assert.AreEqual(1, expr.Calculate());
        }

        [TestMethod]
        public void MinusTest2()
        {
            expr = new MathExpression("x-2");
            Assert.AreEqual(1, expr.Calculate(3));
        }

        [TestMethod]
        public void MinusTest3()
        {
            expr = new MathExpression("sin(x)-2");
            Assert.AreEqual(-1, expr.Calculate(Math.PI / 2));
        }

        [TestMethod]
        public void MultiplicationTest1()
        {
            expr = new MathExpression("3*2");
            Assert.AreEqual(6, expr.Calculate());
        }

        [TestMethod]
        public void MultiplicationTest2()
        {
            expr = new MathExpression("x*2");
            Assert.AreEqual(6, expr.Calculate(3));
        }

        [TestMethod]
        public void DivisionTest1()
        {
            expr = new MathExpression("10/2");
            Assert.AreEqual(5, expr.Calculate());
        }

        [TestMethod]
        public void DivisionTest2()
        {
            expr = new MathExpression("x/2");
            Assert.AreEqual(5, expr.Calculate(10));
        }

        [TestMethod]
        public void ComplexExpressionTest1()
        {
            expr = new MathExpression("x^2+2*x+1");
            Assert.AreEqual(9, expr.Calculate(2));
        }

        [TestMethod]
        public void ComplexExpressionTest2()
        {
            expr = new MathExpression("sin(pi/2)+(x+2)*5");
            Assert.AreEqual(21, expr.Calculate(2));
        }

        [TestMethod]
        public void ComplexExpressionTest3()
        {
            expr = new MathExpression("abs(2-4)+tg(pi)");
            Assert.AreEqual(2, expr.Calculate());
        }

        [TestMethod]
        public void ComplexExpressionTest4()
        {
            expr = new MathExpression("x^2 + 2,5");
            Assert.AreEqual(6.5, expr.Calculate(2));
        }

        [TestMethod]
        public void ComplexExpressionTest5()
        {
            expr = new MathExpression("x^ 2 + 2.5");
            Assert.AreEqual(6.5, expr.Calculate(2));
        }

        [TestMethod]
        public void ComplexExpressionTest6()
        {
            expr = new MathExpression("sin(pi/3.8)*ln(100500)+arctg(2)");
            Console.WriteLine(expr.Calculate().ToString());
            Assert.AreEqual(9.58115, expr.Calculate(), 0.00001);
        }

        [TestMethod]
        public void UnaryMinusTest1()
        {
            expr = new MathExpression("-2+3");
            Console.WriteLine(expr.ToString());
            Assert.AreEqual(1, expr.Calculate());
        }

        [TestMethod]
        public void DerivativeTest()
        {
            MathExpression expr = new MathExpression("e^x+sin(x)");
            Assert.AreEqual(2, expr.Derivative(0), 0.000001);
        }

        [TestMethod]
        public void SeveralVariablesTest1()
        {
            MathExpression expr = new MathExpression("x+y");
            Assert.AreEqual(5, expr.Calculate(new Var("x", 2), new Var("y", 3)));
        }

        [TestMethod]
        public void SeveralVariablesTest2()
        {
            MathExpression expr = new MathExpression("1/(x*y)");
            Assert.AreEqual(0.1, expr.Calculate(new Var("x", 2), new Var("y", 5)));
        }

        [TestMethod]
        public void SeveralVariablesTest3()
        {
            MathExpression expr = new MathExpression("x+y");
            Assert.AreEqual(5, expr.Calculate(2, 3));
        }

        [TestMethod]
        public void SeveralVariablesTest4()
        {
            MathExpression expr = new MathExpression("1/(x*y)");
            Assert.AreEqual(0.1, expr.Calculate(2, 5));
        }

        [TestMethod]
        public void SeveralVariablesTest5()
        {
            MathExpression expr = new MathExpression("x*y+y*z+z*x");
            Assert.AreEqual(11, expr.Calculate(1, 2, 3));
        }

        [TestMethod]
        public void SimpleTest1()
        {
            MathExpression expr = new MathExpression("x");
            Assert.AreEqual(0.5, expr.Calculate(0.5));
        }

        [TestMethod]
        public void SimpleTest2()
        {
            MathExpression expr = new MathExpression("x");
            Assert.AreEqual(2, expr.Calculate(2));
        }

        [TestMethod]
        public void SimpleTest3()
        {
            MathExpression expr = new MathExpression("x");
            Assert.AreEqual(0.1, expr.Calculate(0.1));
        }
    }
}
