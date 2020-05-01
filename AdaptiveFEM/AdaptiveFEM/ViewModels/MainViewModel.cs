using AdaptiveFEM.Models;
using AdaptiveFEM.Services;
using LiveCharts;
using LiveCharts.Wpf;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Text;
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

        #endregion

        public ObservableCollection<NumericalResultItem> NumResults { get; set; }
        public ICommand SolveCommand { get; private set; }


        public SeriesCollection SeriesCollection { get; set; }
        public string[] Labels { get; set; }
        public Func<double, string> YFormatter { get; set; }

        public MainViewModel()
        {
            NumResults = new ObservableCollection<NumericalResultItem>();
            SeriesCollection = new SeriesCollection();
            SolveCommand = new Command(Solve);
            Mu = "1";
            Beta = "100";
            Sigma = "0";
            F = "100";
            A = 0;
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
            NumResults.Clear();
            var mu = new Function(Mu);
            var beta = new Function(Beta);
            var sigma = new Function(Sigma);
            var f = new Function(F);

            var solver = new AdaptiveFEM_H(mu, beta, sigma, f, A, B, Alpha, Gamma, Ua, Ub, Error, N);
            solver.Run();
            ShowTable(solver.Iterations.LastOrDefault());
            var data = solver.Iterations.LastOrDefault();
            var result = new List<NumericalResultItem>();
            for (int i = 0; i < data.Elements.Count; ++i)
            {
                result.Add(new NumericalResultItem(data.Elements[i].Begin, data.Solution[i]));
            }
            result.Add(new NumericalResultItem(data.Elements[data.Elements.Count - 1].End, data.Solution[data.Elements.Count]));

            SeriesCollection.Add(new LineSeries
            {
                Title = "u(x)",
                Values = new ChartValues<double>(result.Select(x => x.Ux).ToList())
            }
            );
            Labels = result.Select(x=>x.X.ToString()).ToArray();
            YFormatter = value => value.ToString("C");
        }

        private void ShowTable(IterationData data)
        {
            for (int i = 0; i < data.Elements.Count; ++i)
            {
                NumResults.Add(new NumericalResultItem(data.Elements[i].Begin, data.Solution[i]));
            }
            NumResults.Add(new NumericalResultItem(data.Elements[data.Elements.Count - 1].End, data.Solution[data.Elements.Count]));
        }
    }
}
