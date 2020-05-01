using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PoohMathParser
{
    /// <summary>
    /// LIFO stack containing tokens. Based on System.Collections.Generic.List<T>
    /// </summary>
    public class Stack
    {
        /// <summary>
        /// List containing tokens.
        /// </summary>
        List<Token> tokens;

        /// <summary>
        /// Default constructor.
        /// </summary>
        public Stack()
        {
            tokens = new List<Token>();
        }

        /// <summary>
        /// Adds a token to the top of stack.
        /// </summary>
        /// <param name="t">Token to add</param>
        public void Push(Token t)
        {
            tokens.Add(t);
        }

        /// <summary>
        /// Removes a token which is currently on the top of stack.
        /// </summary>
        /// <returns>Removed token</returns>
        public Token Pop()
        {
            if (this.Empty() == true)
            {
                throw new IndexOutOfRangeException();
            }

            Token result = tokens[tokens.Count - 1];
            tokens.RemoveAt(tokens.Count - 1);
            return result;
        }

        /// <summary>
        /// Gets a token which is currently on the top of stack.
        /// </summary>
        /// <returns>Token which is currently on the top of stack</returns>
        public Token Top()
        {
            if (this.Empty() == true)
            {
                throw new IndexOutOfRangeException();
            }

            return tokens[tokens.Count - 1];
        }

        /// <summary>
        /// Checks whether stack is empty.
        /// </summary>
        /// <returns>True is stack is empty; else false</returns>
        public bool Empty()
        {
            return (tokens.Count == 0 ? true : false);
        }
        
        /// <summary>
        /// Gets the number of tokens in stack.
        /// </summary>
        /// <returns>Numver of tokens in stack</returns>
        public int Size()
        {
            return tokens.Count;
        }
    }
}
