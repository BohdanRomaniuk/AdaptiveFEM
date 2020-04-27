using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using PoohMathParser;

namespace VisualTester
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private void buttonCalculateExpr_Click(object sender, RoutedEventArgs e)
        {
            listBoxTokens.Items.Clear();
            textBoxReversePolishNotation.Clear();
            textBoxResult.Clear();

            string input = textBoxInput.Text;
            MathExpression expr = new MathExpression(input);

            foreach (Token t in expr.Tokens)
            {
                listBoxTokens.Items.Add(t);
            }

            foreach (Token t in expr.ReversePolishNotation)
            {
                textBoxReversePolishNotation.Text += t.Lexeme;
            }

            double result = 0;
            if (textBoxVar.Text != "")
            {
                double var = double.Parse(textBoxVar.Text);
                result = expr.Calculate(var);
            }
            else
            {
                result = expr.Calculate();
            }
            textBoxResult.Text = result.ToString();

            MathExpression expr2 = new MathExpression("1/(x*y)");
            MessageBox.Show(expr2.Calculate(new Var("x", 2), new Var("y", 5)).ToString());
        }
    }
}
