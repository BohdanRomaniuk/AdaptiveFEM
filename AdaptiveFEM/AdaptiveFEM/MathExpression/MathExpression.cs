using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Globalization;

namespace PoohMathParser
{
    /// <summary>
    /// Class for storing mathematical expressions. 
    /// </summary>
    public class MathExpression
    {
        private string expression;
        private List<Token> tokens;
        private List<Token> reversePolishNotation;

        public List<Token> Tokens
        {
            get
            {
                return tokens;
            }
        }

        public List<Token> ReversePolishNotation
        {
            get
            {
                return reversePolishNotation;
            }
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="expression">String expression to parse</param>
        public MathExpression(string expression)
        {
            this.expression = PrepareString(expression);
            this.tokens = GetTokens(this.expression);
            this.reversePolishNotation = ConvertToReversePolishNotation(this.tokens);
        }

        #region Dictionaries with operators, functions etc. and their delegates

        private Dictionary<string, Func<double, double, double>> operators =
            new Dictionary<string, Func<double, double, double>>
            {
                { "+", (operand1, operand2) => Plus(operand1, operand2) },
                { "-", (operand1, operand2) => Minus(operand1, operand2) },
                { "*", (operand1, operand2) => Multiplication(operand1, operand2) },
                { "/", (operand1, operand2) => Division(operand1, operand2) },
                { "^", (operand1, operand2) => Pow(operand1, operand2) },
                { "(", null },
                { ")", null },
            };

        private Dictionary<string, Func<double, double>> functions = new Dictionary<string, Func<double, double>>
            {
                { "sin", (operand) => Sin(operand) },
                { "cos", (operand) => Cos(operand) },
                { "tg", (operand) => Tg(operand) },
                { "ctg", (operand) => Ctg(operand) },
                { "arcsin", (operand) => Arcsin(operand) },
                { "arccos", (operand) => Arccos(operand) },
                { "arctg", (operand) => Arctg(operand) },
                { "arcctg", (operand) => Arcctg(operand) },
                { "sinh", (operand) => Sinh(operand) },
                { "cosh", (operand) => Cosh(operand) },
                { "tgh", (operand) => Tgh(operand) },
                { "ctgh", (operand) => Ctgh(operand) },
                { "ln", (operand) => Ln(operand) },
                { "lg", (operand) => Lg(operand) },
                { "sqrt", (operand) => Sqrt(operand) },
                { "abs", (operand) => Abs(operand) },
                { "sign", (operand) => Sign(operand) },
            };

        private Dictionary<string, double> constants = new Dictionary<string, double>
            {
                { "e", Math.E },
                { "pi", Math.PI },
            };

        private List<char> variables = new List<char>()
            {
                'a', 'b', 'c', 'd', 'i', 'j', 'k', 'l', 'm', 'n', 'p', 'q', 'r', 's', 't', 'u', 'v', 'x', 'y', 'z'
            };

        #endregion

        #region Operators implementation

        private static double Plus(double operand1, double operand2)
        {
            return operand1 + operand2;
        }
        private static double Minus(double operand1, double operand2)
        {
            return operand2 - operand1;
        }
        private static double Multiplication(double operand1, double operand2)
        {
            return operand1 * operand2;
        }
        private static double Division(double operand1, double operand2)
        {
            return operand2 / operand1;
        }
        private static double Pow(double operand1, double operand2)
        {
            return Math.Pow(operand2, operand1);
        }

        #endregion

        #region Functions implementation

        private static double Sin(double operand)
        {
            return Math.Sin(operand);
        }
        private static double Cos(double operand)
        {
            return Math.Cos(operand);
        }
        private static double Tg(double operand)
        {
            return Math.Tan(operand);
        }
        private static double Ctg(double operand)
        {
            return 1 / Math.Tan(operand);
        }
        private static double Arcsin(double operand)
        {
            return Math.Asin(operand);
        }
        private static double Arccos(double operand)
        {
            return Math.Acos(operand);
        }
        private static double Arctg(double operand)
        {
            return Math.Atan(operand);
        }
        private static double Arcctg(double operand)
        {
            return 1 / Math.Atan(operand);
        }
        private static double Sinh(double operand)
        {
            return Math.Sinh(operand);
        }
        private static double Cosh(double operand)
        {
            return Math.Cosh(operand);
        }
        private static double Tgh(double operand)
        {
            return Math.Tanh(operand);
        }
        private static double Ctgh(double operand)
        {
            return 1 / Math.Tanh(operand);
        }
        private static double Ln(double operand)
        {
            return Math.Log(operand);
        }
        private static double Lg(double operand)
        {
            return Math.Log10(operand);
        }
        private static double Sqrt(double operand)
        {
            return Math.Sqrt(operand);
        }
        private static double Abs(double operand)
        {
            return Math.Abs(operand);
        }
        private static double Sign(double operand)
        {
            return Math.Sign(operand);
        }

        #endregion

        #region Methods for string preparation, splitting string to tokens and so on.

        /// <summary>
        /// Preparing string for it's further convertion to the sequence of tokens. 
        /// This method is deleting whitespaces etc.
        /// </summary>
        /// <param name="expression">String expression to prepare</param>
        /// <returns>String prepared for converting to the sequence of tokens</returns>
        private string PrepareString(string expression)
        {
            //Deleting all the spaces from the string
            expression = expression.Replace(" ", string.Empty);
            expression = expression.Replace(',', '.');
            return expression;
        }

        /// <summary>
        /// Parses a string of characters and converts it to sequence of tokens.
        /// </summary>
        /// <returns>Sequence of tokens</returns>
        private List<Token> GetTokens(string expression)
        {
            List<Token> tokens = new List<Token>();

            bool isError;
            bool success;
            for (int i = 0; i < expression.Length; ++i)
            {
                isError = true;
                success = false;

                string number = "";
                double num;
                while (double.TryParse(expression[i].ToString(CultureInfo.InvariantCulture), System.Globalization.NumberStyles.Number, CultureInfo.InvariantCulture, out num) || expression[i] == '.')
                {
                    if (expression[i] == '.')
                    {
                        number += ".";
                    }
                    else
                    {
                        number += num.ToString(CultureInfo.InvariantCulture);
                    }

                    if (i < expression.Length - 1)
                    {
                        ++i;
                    }
                    else break;
                }
                if (number != "")
                {
                    Token t = new Token(number, TokenType.Number);
                    tokens.Add(t);
                    success = true;
                    isError = false;
                    number = "";
                }

                foreach (string s in operators.Keys)
                {
                    if (expression[i].ToString(CultureInfo.InvariantCulture) == s)
                    {
                        if (isUnaryOparator(tokens))
                        {
                            expression = expression.Insert(i, "0");
                            --i;
                        }
                        else
                        {
                            Token t = new Token(expression[i].ToString(CultureInfo.InvariantCulture), TokenType.Operator);
                            tokens.Add(t);
                        }
                        success = true;
                        isError = false;
                        break;
                    }
                }

                foreach (string s in functions.Keys)
                {
                    if (success)
                    {
                        break;
                    }

                    string function = "";

                    int k = i;
                    for (int j = 0; j < s.Length; ++j)
                    {
                        if (s[j] == expression[k])
                        {
                            function += expression[k];
                        }
                        if (k < expression.Length - 1)
                        {
                            ++k;
                        }
                        else break;
                    }
                    if (function == s)
                    {
                        Token t = new Token(function, TokenType.Function);
                        tokens.Add(t);
                        isError = false;
                        success = true;
                        i += function.Length - 1;
                        break;
                    }
                }

                foreach (string s in constants.Keys)
                {
                    if (success)
                    {
                        break;
                    }

                    string constant = "";

                    int k = i;
                    for (int j = 0; j < s.Length; ++j)
                    {
                        if (s[j] == expression[k])
                        {
                            constant += expression[k];
                        }
                        if (k < expression.Length - 1)
                        {
                            ++k;
                        }
                        else break;
                    }
                    if (constant == s)
                    {
                        Token t = new Token(constant, TokenType.Constant);
                        tokens.Add(t);
                        isError = false;
                        success = true;
                        i += constant.Length - 1;
                        break;
                    }
                }

                foreach (char c in variables)
                {
                    if (success)
                    {
                        break;
                    }

                    if (expression[i] == c)
                    {
                        Token t = new Token(expression[i].ToString(CultureInfo.InvariantCulture), TokenType.Variable);
                        tokens.Add(t);
                        isError = false;
                        success = true;
                        break;
                    }
                }

                if (isError)
                {
                    throw new ArgumentException("Bad expression string.");
                }
            }

            return tokens;
        }

        /// <summary>
        /// Converts sequence of tokens from infix notation to postfix (reverse polish notation). 
        /// See Shunting-yard algorithm for details.
        /// </summary>
        /// <param name="tokens">Sequence of tokens in infix notation</param>
        /// <returns>Sequence of tokens in reverse polish notation</returns>
        private List<Token> ConvertToReversePolishNotation(List<Token> tokens)
        {
            List<Token> reversePolishNotation = new List<Token>();
            Stack stack = new Stack();

            foreach (Token t in tokens)
            {
                if (t.Type == TokenType.Number || t.Type == TokenType.Variable || t.Type == TokenType.Constant)
                {
                    reversePolishNotation.Add(t);
                }
                else if (t.Type == TokenType.Function)
                {
                    stack.Push(t);
                }
                else if (t.Lexeme == "(")
                {
                    stack.Push(t);
                }
                else if (t.Type == TokenType.Operator && t.Lexeme != ")")
                {
                    if (!stack.Empty())
                    {
                        while (stack.Top().Type == TokenType.Operator && t.Priority <= stack.Top().Priority)
                        {
                            reversePolishNotation.Add(stack.Pop());
                            if (stack.Empty())
                            {
                                break;
                            }
                        }
                        stack.Push(t);
                    }
                    else
                    {
                        stack.Push(t);
                    }
                }
                else if (t.Type == TokenType.Operator && t.Lexeme == ")")
                {
                    if (!stack.Empty())
                    {
                        while (stack.Top().Lexeme != "(")
                        {
                            reversePolishNotation.Add(stack.Pop());
                            if (stack.Empty())
                            {
                                break;
                            }
                        }
                        stack.Pop();
                        if (!stack.Empty())
                        {
                            if (stack.Top().Type == TokenType.Function)
                            {
                                reversePolishNotation.Add(stack.Pop());
                            }
                        }
                    }
                }
            }

            while (!stack.Empty())
            {
                reversePolishNotation.Add(stack.Pop());
            }

            return reversePolishNotation;
        }

        #endregion

        /// <summary>
        /// Calculates value of the expression without variables.
        /// </summary>
        /// <returns>Value of the expression</returns>
        public double Calculate()
        {
            return this.Calculate(0);
        }

        /// <summary>
        /// Calculates value of the expression of one variable.
        /// </summary>
        /// <param name="variableValue">Value of the variable</param>
        /// <returns>Value of the expression</returns>
        public double Calculate(double variableValue)
        {
            string variableName = "";
            foreach (Token t in this.tokens)
            {
                if (t.Type == TokenType.Variable)
                {
                    variableName = t.Lexeme;
                    break;
                }
            }

            return this.Calculate(new Var(variableName, variableValue));
        }

        /// <summary>
        /// Calculates value of the expression of many variables.
        /// </summary>
        /// <param name="variables">Variables</param>
        /// <returns>Value of the expression</returns>
        public double Calculate(params Var[] variables)
        {
            Stack stack = new Stack();

            foreach (Token t in reversePolishNotation)
            {
                if (t.Type == TokenType.Number)
                {
                    stack.Push(t);
                }
                else if (t.Type == TokenType.Operator)
                {
                    string operand1Str = stack.Pop().Lexeme.Replace(',', '.');
                    string operand2Str = stack.Pop().Lexeme.Replace(',', '.');
                    double operand1 = double.Parse(operand1Str, CultureInfo.InvariantCulture);
                    double operand2 = double.Parse(operand2Str, CultureInfo.InvariantCulture);
                    double smallResult = operators[t.Lexeme](operand1, operand2);
                    stack.Push(new Token(smallResult.ToString(CultureInfo.InvariantCulture), TokenType.Number));
                }
                else if (t.Type == TokenType.Function)
                {
                    string operandStr = stack.Pop().Lexeme.Replace(',', '.');
                    double operand = double.Parse(operandStr, CultureInfo.InvariantCulture);
                    double smallResult = functions[t.Lexeme](operand);
                    stack.Push(new Token(smallResult.ToString(CultureInfo.InvariantCulture), TokenType.Number));
                }
                else if (t.Type == TokenType.Constant)
                {
                    double constant = constants[t.Lexeme];
                    stack.Push(new Token(constant.ToString(CultureInfo.InvariantCulture), TokenType.Number));
                }
                else if (t.Type == TokenType.Variable)
                {
                    Token smallToken = new Token();

                    foreach (Var variable in variables)
                    {
                        if (t.Lexeme == variable.Name)
                        {
                            smallToken = new Token(variable.Value.ToString(CultureInfo.InvariantCulture), TokenType.Number);
                            break;
                        }
                    }

                    stack.Push(smallToken);
                }
                else
                {
                    throw new ArgumentException("Bad expression string.");
                }
            }

            double result = double.Parse(stack.Pop().Lexeme, CultureInfo.InvariantCulture);

            return result;
        }

        /// <summary>
        /// Calculates value of the expression of many variables.
        /// </summary>
        /// <param name="variables">Values of variables</param>
        /// <returns>Value of the expression</returns>
        public double Calculate(params double[] variables)
        {
            Dictionary<string, double> varNamesAndValues = new Dictionary<string, double>();

            for (int i = 0; i < variables.Count(); ++i)
            {
                foreach (Token t in tokens)
                {
                    if (t.Type == TokenType.Variable)
                    {
                        if (!varNamesAndValues.ContainsKey(t.Lexeme))
                        {
                            varNamesAndValues.Add(t.Lexeme, variables[i]);
                            break;
                        }
                    }
                }
            }

            Var[] varNamesAndValuesArr = new Var[varNamesAndValues.Count];

            int j = 0;
            foreach (KeyValuePair<string, double> kv in varNamesAndValues)
            {
                varNamesAndValuesArr[j] = new Var(kv.Key, kv.Value);
                ++j;
            }

            return this.Calculate(varNamesAndValuesArr);
        }

        /// <summary>
        /// Calculates the first derivative of the expression at specified point.
        /// </summary>
        /// <param name="x0">Point derivative is calculated at</param>
        /// <returns></returns>
        public double Derivative(double x0)
        {
            double result;

            double dx = 10e-10;
            result = (this.Calculate(x0 + dx) - this.Calculate(x0)) / dx;

            return result;
        }

        /// <summary>
        /// Calculates the second derivative of the expression at specified point.
        /// </summary>
        /// <param name="x0">Point derivative is calculated at</param>
        /// <returns></returns>
        public double SecondDerivative(double x0)
        {
            double result;

            double dx = Math.Pow(10, -5);
            result = (2 / Math.Pow(dx, 2)) * ((this.Calculate(x0 + dx)
                + this.Calculate(x0 - dx)) / 2 - this.Calculate(x0));

            return result;
        }

        /// <summary>
        /// Converts expression to string. 
        /// </summary>
        /// <returns>String representation of the expression</returns>
        public override string ToString()
        {
            return expression;
        }

        private bool isUnaryOparator(List<Token> pTokens)
        {
            if (pTokens == null || pTokens.Capacity == 0)
            {
                return true;
            }
            else
            {
                Token token = pTokens[pTokens.Count - 1];
                return token.Lexeme == "(";
            }
        }
    }
}

/// <summary>
/// Represents a variable in MathExpression
/// </summary>
public class Var
{
    /// <summary>
    /// Name of the variable
    /// </summary>
    public string Name
    {
        get;
        set;
    }

    /// <summary>
    /// Value of the variable
    /// </summary>
    public double Value
    {
        get;
        set;
    }

    /// <summary>
    /// Creates new variable
    /// </summary>
    /// <param name="name">Name of the variable</param>
    /// <param name="value">Value of the variable</param>
    public Var(string name, double value)
    {
        this.Name = name;
        this.Value = value;
    }
}
