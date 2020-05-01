using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PoohMathParser
{
    /// <summary>
    /// Token is string of one or more characters which represents
    /// either a number, or an operator, a function, a variable or a constant.
    /// </summary>
    public class Token
    {
        #region Fields

        /// <summary>
        /// A string of characters.
        /// </summary>
        private string lexeme;

        /// <summary>
        /// Type of the token (number, function etc.)
        /// </summary>
        private TokenType type;

        /// <summary>
        /// Priority of the token.
        /// </summary>
        private int priority;

        #endregion

        #region Constructors

        /// <summary>
        /// Default constructor.
        /// </summary>
        public Token()
        {
            lexeme = "";
        }

        /// <summary>
        /// Constructor with parameters.
        /// </summary>
        /// <param name="lexeme">Lexeme of the token</param>
        /// <param name="type">Type of the token</param>
        public Token(string lexeme, TokenType type)
        {
            this.lexeme = lexeme;
            this.type = type;
            switch (lexeme)
            {
                case "^": this.priority = 3; break;
                case "*": this.priority = 2; break;
                case "/": this.priority = 2; break;
                case "+": this.priority = 1; break;
                case "-": this.priority = 1; break;
                default: this.priority = 0; break;
            }
        }

        #endregion

        #region Properties

        /// <summary>
        /// Gets lexeme of the token.
        /// </summary>
        public string Lexeme
        {
            get { return lexeme; }
        }

        /// <summary>
        /// Gets type of the token.
        /// </summary>
        public TokenType Type
        {
            get { return type; }
        }

        /// <summary>
        /// Gets priority of the token.
        /// </summary>
        public int Priority
        {
            get { return priority; }
        }

        #endregion

        /// <summary>
        /// Converts token to it's string representation.
        /// </summary>
        /// <returns>String representation of the token</returns>
        public override string ToString()
        {
            string result = String.Format("Lexeme: {0}, Token type: {1}", lexeme, type);
            return result;
        }
    }

    /// <summary>
    /// Type of the token
    /// </summary>
    public enum TokenType { Number, Operator, Function, Variable, Constant };
}
