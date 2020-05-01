using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Globalization;
using ZedGraph;

namespace AdaptiveFEM
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        GraphPane pane;
        AdaptiveFEM_H solver;

        public MainWindow()
        {
            InitializeComponent();
            InitializeGraphPane();
            //this.fTxt.Text = "cos(pi*x)*cos(pi*x)+0,005*pi*pi*cos(2*pi*x)";
            //this.muTxt.Text = "0.0025";
            //this.betaTxt.Text = "0";
            //this.sigmaTxt.Text = "1";
            //this.bTxt.Text = "1";
        }

        public void InitializeGraphPane()
        {
            canvas.GraphPane.XAxis.Scale.Min = 0;
            canvas.GraphPane.XAxis.Scale.Max = 10;
            canvas.GraphPane.YAxis.Scale.Min = 0;
            canvas.GraphPane.YAxis.Scale.Max = 10;
            canvas.GraphPane.Title.Text = "";
            canvas.GraphPane.XAxis.Title.Text = "x";
            canvas.GraphPane.YAxis.Title.Text = "u(x)";
            canvas.AxisChange();
            canvas.Invalidate();
            pane = canvas.GraphPane;
        }

        private void DrawGraph(IterationData data)
        {
            canvas.GraphPane.CurveList.Clear();
            PointPairList numList = new PointPairList();
            canvas.GraphPane.XAxis.Scale.Min = solver.A;
            canvas.GraphPane.XAxis.Scale.Max = solver.B;
            double min = data.Solution.Min();
            double max = data.Solution.Max();
            canvas.GraphPane.XAxis.Title.Text = "x";
            canvas.GraphPane.YAxis.Title.Text = "U";
            canvas.GraphPane.YAxis.Scale.Min = min - 0.25 * Math.Abs(min);
            canvas.GraphPane.YAxis.Scale.Max = max + 0.25 * Math.Abs(max);
            canvas.GraphPane.Title.Text = data.Elements.Count.ToString() + " Finite Elements";
            canvas.AxisChange();
            canvas.Invalidate();
            for (int i = 0; i < data.Elements.Count; ++i)
            {
                numList.Add(data.Elements[i].Begin, data.Solution[i]);
            }
            numList.Add(data.Elements[data.Elements.Count - 1].End, data.Solution[data.Elements.Count]);
            pane.AddCurve("", numList, System.Drawing.Color.BlueViolet, SymbolType.Star);
            canvas.AxisChange();
            canvas.Invalidate();
        }

        private void DrawErrorNorm(IterationData data)
        {
            canvas.GraphPane.CurveList.Clear();
            PointPairList numList = new PointPairList();
            canvas.GraphPane.XAxis.Scale.Min = solver.A;
            canvas.GraphPane.XAxis.Scale.Max = solver.B;
            canvas.GraphPane.XAxis.Title.Text = "x";
            canvas.GraphPane.YAxis.Title.Text = "Norm of error";
            double min = data.ErrorsNormsV.Min();
            double max = data.ErrorsNormsV.Max();
            canvas.GraphPane.YAxis.Scale.Min = min - 0.25 * Math.Abs(min);
            canvas.GraphPane.YAxis.Scale.Max = max + 0.25 * Math.Abs(max);
            canvas.GraphPane.Title.Text = data.Elements.Count.ToString() + " Finite Elements";
            canvas.AxisChange();
            canvas.Invalidate();
            for (int i = 0; i < data.Elements.Count; ++i)
            {
                numList.Add(data.Elements[i].MidPoint, data.ErrorsNormsV[i]);
            }
            pane.AddCurve("", numList, System.Drawing.Color.BlueViolet, SymbolType.Star);
            canvas.AxisChange();
            canvas.Invalidate();
        }

        private void DrawRelError(IterationData data)
        {
            canvas.GraphPane.CurveList.Clear();
            PointPairList numList = new PointPairList();
            canvas.GraphPane.XAxis.Scale.Min = solver.A;
            canvas.GraphPane.XAxis.Scale.Max = solver.B;
            canvas.GraphPane.XAxis.Title.Text = "x";
            canvas.GraphPane.YAxis.Title.Text = "Relative error";
            double min = data.Errors.Min();
            double max = data.Errors.Max();
            canvas.GraphPane.YAxis.Scale.Min = min - 0.25 * Math.Abs(min);
            canvas.GraphPane.YAxis.Scale.Max = max + 0.25 * Math.Abs(max);
            canvas.GraphPane.Title.Text = data.Elements.Count.ToString() + " Finite Elements";
            canvas.AxisChange();
            canvas.Invalidate();
            for (int i = 0; i < data.Elements.Count; ++i)
            {
                numList.Add(data.Elements[i].MidPoint, data.Errors[i]);
            }
            pane.AddCurve("", numList, System.Drawing.Color.BlueViolet, SymbolType.Star);
            canvas.AxisChange();
            canvas.Invalidate();
        }

        private void DrawDeriv(IterationData data)
        {
            canvas.GraphPane.CurveList.Clear();
            PointPairList numList = new PointPairList();
            canvas.GraphPane.XAxis.Scale.Min = solver.A;
            canvas.GraphPane.XAxis.Scale.Max = solver.B;
            canvas.GraphPane.XAxis.Title.Text = "x";
            canvas.GraphPane.YAxis.Title.Text = "du / dx";
            double min = data.SolutionCenterDeriv.Min();
            double max = data.SolutionCenterDeriv.Max();
            canvas.GraphPane.YAxis.Scale.Min = min - 0.25 * Math.Abs(min);
            canvas.GraphPane.YAxis.Scale.Max = max + 0.25 * Math.Abs(max);
            canvas.GraphPane.Title.Text = data.Elements.Count.ToString() + " Finite Elements";
            canvas.AxisChange();
            canvas.Invalidate();
            for (int i = 0; i < data.Elements.Count; ++i)
            {
                numList.Add(data.Elements[i].MidPoint, data.SolutionCenterDeriv[i]);
            }
            pane.AddCurve("", numList, System.Drawing.Color.BlueViolet, SymbolType.Star);
            canvas.AxisChange();
            canvas.Invalidate();
        }

        private void ShowTable(IterationData data)
        {
            currentN.Content = data.Elements.Count.ToString();
            List<NumericalResultItem> result = new List<NumericalResultItem>(3);
            for (int i = 0; i < data.Elements.Count; ++i)
            {
                result.Add(new NumericalResultItem(data.Elements[i].Begin, data.Solution[i]));
            }
            result.Add(new NumericalResultItem(data.Elements[data.Elements.Count - 1].End, data.Solution[data.Elements.Count]));
            table.ItemsSource = result;

            List<NumericalResultErrorItem> errors = new List<NumericalResultErrorItem>();
            for (int i = 0; i < data.Elements.Count; ++i)
            {
                errors.Add(new NumericalResultErrorItem(data.Elements[i].MidPoint, data.SolutionCenter[i], data.SolutionCenterDeriv[i], data.ErrorsNormsV[i], data.Errors[i]));
            }
            tableErrors.ItemsSource = errors;
        }

        private void solveBtn_Click(object sender, RoutedEventArgs e)
        {
            Function mu = new Function(muTxt.Text);
            Function beta = new Function(betaTxt.Text);
            Function sigma = new Function(sigmaTxt.Text);
            Function f = new Function(fTxt.Text);

            double a = Convert.ToDouble(aTxt.Text);
            double b = Convert.ToDouble(bTxt.Text);
            double alpha = Convert.ToDouble(alphaTxt.Text);
            double gamma = Convert.ToDouble(gammaTxt.Text);
            double u_a = Convert.ToDouble(u_aTxt.Text);
            double u_b = Convert.ToDouble(u_bTxt.Text);
            double error = Convert.ToDouble(errorTxt.Text);
            int N = Convert.ToInt32(NTxt.Text);

            solver = new AdaptiveFEM_H(mu, beta, sigma, f, a, b, alpha, gamma, u_a, u_b, error, N);

            DateTime begin = DateTime.Now;
            solver.Run();
            DateTime end = DateTime.Now;

            iterationLabel.Content = solver.Iterations.Count.ToString();
            timeLabel.Content = (end - begin).ToString(@"hh\:mm\:ss\.ff");

            listBox.Items.Clear();
            for (int i = 0; i < solver.Iterations.Count; ++i)
            {
                Button btn = new Button();
                btn.Content = solver.Iterations[i].Elements.Count.ToString();
                btn.Height = 25;
                btn.Width = 50;
                btn.HorizontalAlignment = System.Windows.HorizontalAlignment.Left;
                btn.Click += iterationBtnClick;
                listBox.Items.Add(btn);
            }
            DrawGraph(solver.Iterations[solver.Iterations.Count - 1]);
            ShowTable(solver.Iterations[solver.Iterations.Count - 1]);
            //errorNormLabel.Content = solver.Iterations[solver.Iterations.Count - 1].ErrorNormV.ToString();
            drawSolutionBtn.IsEnabled = true;
            drawErrorNormBtn.IsEnabled = true;
            drawRellErrorBtn.IsEnabled = true;
            drawU_xDerivBtn.IsEnabled = true;
        }

        void iterationBtnClick(object sender, RoutedEventArgs e)
        {
            Button btn = (Button)sender;
            int N = Convert.ToUInt16(btn.Content.ToString());
            int numberOfIteration = 0;
            for (int i = 0; i < solver.Iterations.Count; ++i)
            {
                if (solver.Iterations[i].Elements.Count == N)
                {
                    numberOfIteration = i;
                    break;
                }
            }
            DrawGraph(solver.Iterations[numberOfIteration]);
            ShowTable(solver.Iterations[numberOfIteration]);
            //errorNormLabel.Content = solver.Iterations[numberOfIteration].ErrorNormV.ToString();
        }

        private void drawSolutionBtn_Click(object sender, RoutedEventArgs e)
        {
            int N = Convert.ToUInt16(currentN.Content.ToString());
            int numberOfIteration = 0;
            for (int i = 0; i < solver.Iterations.Count; ++i)
            {
                if (solver.Iterations[i].Elements.Count == N)
                {
                    numberOfIteration = i;
                    break;
                }
            }
            DrawGraph(solver.Iterations[numberOfIteration]);
        }

        private void drawErrorNormBtn_Click(object sender, RoutedEventArgs e)
        {
            int N = Convert.ToUInt16(currentN.Content.ToString());
            int numberOfIteration = 0;
            for (int i = 0; i < solver.Iterations.Count; ++i)
            {
                if (solver.Iterations[i].Elements.Count == N)
                {
                    numberOfIteration = i;
                    break;
                }
            }
            DrawErrorNorm(solver.Iterations[numberOfIteration]);
        }

        private void drawRellErrorBtn_Click(object sender, RoutedEventArgs e)
        {
            int N = Convert.ToUInt16(currentN.Content.ToString());
            int numberOfIteration = 0;
            for (int i = 0; i < solver.Iterations.Count; ++i)
            {
                if (solver.Iterations[i].Elements.Count == N)
                {
                    numberOfIteration = i;
                    break;
                }
            }
            DrawRelError(solver.Iterations[numberOfIteration]);
        }

        private void drawU_xDeriv_Click(object sender, RoutedEventArgs e)
        {
            int N = Convert.ToUInt16(currentN.Content.ToString());
            int numberOfIteration = 0;
            for (int i = 0; i < solver.Iterations.Count; ++i)
            {
                if (solver.Iterations[i].Elements.Count == N)
                {
                    numberOfIteration = i;
                    break;
                }
            }
            DrawDeriv(solver.Iterations[numberOfIteration]);
        }
    }
}