using System;
using System.Collections.Generic;
using System.Text;

namespace AdaptiveFEM.Models
{
    public class Error
    {
        public double XCenter { get; set; }
        public double UxCenter { get; set; }
        public double UxDeriv { get; set; }
        public double ErrorNorm { get; set; }
        public double RelError { get; set; }

        public Error(double xc, double uxc, double uxd, double ern, double rele)
        {
            XCenter = xc;
            UxCenter = uxc;
            UxDeriv = uxd;
            ErrorNorm = ern;
            RelError = rele;
        }
    }
}
