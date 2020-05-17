using AdaptiveFEM.Helpers;
using AdaptiveFEM.Models;
using AdaptiveFEM.Services;
using LiveCharts;
using LiveCharts.Defaults;
using LiveCharts.Wpf;
using NCalc;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Windows.Input;
using System.Windows.Media;

namespace AdaptiveFEM.ViewModels
{
    public class MainViewModel : BaseViewModel
    {
        #region props
        private string mu;
        private string beta;
        private string sigma;
        private string f;

        private double a;
        private double b;
        private double alpha;
        private double gamma;
        private double ua;
        private double ub;
        private double error;
        private int n;

        private FEMSolver solver;
        private int iterationCount;
        private int elementsCount;
        private int currentIteration;

        public string Mu
        {
            get => mu;
            set
            {
                mu = value;
                OnPropertyChanged(nameof(Mu));
            }
        }
        public string Beta
        {
            get => beta;
            set
            {
                beta = value;
                OnPropertyChanged(nameof(Beta));
            }
        }
        public string Sigma
        {
            get => sigma;
            set
            {
                sigma = value;
                OnPropertyChanged(nameof(Sigma));
            }
        }
        public string F
        {
            get => f;
            set
            {
                f = value;
                OnPropertyChanged(nameof(F));
            }
        }

        public double A
        {
            get => a;
            set
            {
                a = value;
                OnPropertyChanged(nameof(A));
            }
        }
        public double B
        {
            get => b;
            set
            {
                b = value;
                OnPropertyChanged(nameof(B));
            }
        }
        public double Alpha
        {
            get => alpha;
            set
            {
                alpha = value;
                OnPropertyChanged(nameof(Alpha));
            }
        }
        public double Gamma
        {
            get => gamma;
            set
            {
                gamma = value;
                OnPropertyChanged(nameof(Gamma));
            }
        }
        public double Ua
        {
            get => ua;
            set
            {
                ua = value;
                OnPropertyChanged(nameof(Ua));
            }
        }
        public double Ub
        {
            get => ub;
            set
            {
                ub = value;
                OnPropertyChanged(nameof(Ub));
            }
        }
        public double Error
        {
            get => error;
            set
            {
                error = value;
                OnPropertyChanged(nameof(Error));
            }
        }
        public int N
        {
            get => n;
            set
            {
                n = value;
                OnPropertyChanged(nameof(N));
            }
        }

        public FEMSolver Solver
        {
            get => solver;
            set
            {
                solver = value;
                OnPropertyChanged(nameof(Solver));
            }
        }
        public int IterationCount
        {
            get => iterationCount;
            set
            {
                iterationCount = value;
                OnPropertyChanged(nameof(IterationCount));
            }
        }
        public int ElementsCount
        {
            get => elementsCount;
            set
            {
                elementsCount = value;
                OnPropertyChanged(nameof(ElementsCount));
            }
        }
        public int CurrentIteration
        {
            get => currentIteration;
            set
            {
                if (currentIteration != value && value != 0 && value <= IterationCount)
                {
                    currentIteration = value;
                    if (Solver.Iterations?.Count != 0)
                    {
                        var lastIter = Solver.Iterations[currentIteration - 1];
                        ElementsCount = lastIter.Elements.Count;
                        var result = PopulateResult(lastIter);
                        ShowTable(result);
                        DrawGraphic(result);
                    }
                    OnPropertyChanged(nameof(CurrentIteration));
                }
            }
        }
        #endregion

        private ObservableCollection<KeyValuePair<double, double>> dots;

        public ObservableCollection<KeyValuePair<double, double>> Dots
        {
            get => dots;
            set
            {
                dots = value;
                OnPropertyChanged(nameof(SeriesCollection));
            }
        }

        private string expectedFunction;
        public string ExpectedFunction
        {
            get => expectedFunction;
            set
            {
                expectedFunction = value;
                OnPropertyChanged(nameof(ExpectedFunction));
            }
        }

        public ObservableCollection<Solution> NumResults { get; set; }
        public ICommand SolveCommand { get; private set; }
        public ICommand ClearCommand { get; private set; }

        public MainViewModel()
        {
            NumResults = new ObservableCollection<Solution>();
            Dots = new ObservableCollection<KeyValuePair<double, double>>();

            SolveCommand = new Command(Solve);
            ClearCommand = new Command(Clear);

            Mu = "-1.0";
            Beta = "1500*Pow([X],8)";
            Sigma = "80*Pow([X],2)";
            F = "100*Exp(Pow([X]-0.15,17))*[X]";
            A = -1;
            B = 1;
            N = 4;
            Alpha = 100000;
            Gamma = 100000;
            Error = 10;
            Ua = 0;
            Ub = 0;
        }

        private void Solve(object parameter)
        {
            currentIteration = 0;

            var mu = new Expression(Mu);
            var beta = new Expression(Beta);
            var sigma = new Expression(Sigma);
            var f = new Expression(F);

            Solver = new FEMSolver(mu, beta, sigma, f, A, B, Alpha, Gamma, Ua, Ub, Error, N);
            Solver.Solve();

            IterationCount = Solver.Iterations.Count;
            CurrentIteration = Solver.Iterations.Count;
        }

        private void Clear(object parameter)
        {
            currentIteration = 0;
            Dots.Clear();
            NumResults.Clear();
        }

        private List<Solution> PopulateResult(Iteration data)
        {
            var result = new List<Solution>();
            for (int i = 0; i < data.Elements.Count; ++i)
            {
                result.Add(new Solution(data.Elements[i].Begin, data.Solution[i]));
            }
            result.Add(new Solution(data.Elements[data.Elements.Count - 1].End, data.Solution[data.Elements.Count]));
            return result;
        }

        private void ShowTable(List<Solution> numList)
        {
            NumResults.Clear();
            foreach (var num in numList)
            {
                NumResults.Add(num);
            }
        }

        private void DrawGraphic(List<Solution> numList)
        {
            Dots.Clear();
            foreach (var elem in numList)
            {
                Dots.Add(new KeyValuePair<double, double>(elem.X, elem.Ux));
            }
            OnPropertyChanged(nameof(Dots));
        }

        private string Subscript(int number)
        {
            var strNumber = number.ToString();
            if (strNumber.Length == 1)
            {
                switch (number)
                {
                    case 0: return "\u2080";
                    case 1: return "\u2081";
                    case 2: return "\u2082";
                    case 3: return "\u2083";
                    case 4: return "\u2084";
                    case 5: return "\u2085";
                    case 6: return "\u2086";
                    case 7: return "\u2087";
                    case 8: return "\u2088";
                    case 9: return "\u2089";
                }
            }

            var result = string.Empty;
            foreach (var symbol in strNumber)
            {
                result += Subscript(Convert.ToInt32(symbol.ToString()));
            }
            return result;
        }
    }
}
