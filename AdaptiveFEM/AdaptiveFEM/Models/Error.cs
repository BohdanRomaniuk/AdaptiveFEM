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
        public double UNorm { get; set; }
        public double ENormI { get; set; }
        public double ENormV { get; set; }
        public double EthaError { get; set; }

        public Error(double xc, double uxc, double uxd, double enormi, double ethaerr)
        {
            XCenter = xc;
            UxCenter = uxc;
            UxDeriv = uxd;
            ENormI = enormi;
            EthaError = ethaerr;
        }
    }
}
