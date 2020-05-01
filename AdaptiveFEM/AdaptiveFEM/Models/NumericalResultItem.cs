namespace AdaptiveFEM.Models
{
    public class NumericalResultItem
    {
        public NumericalResultItem(double x, double u_x)
        {
            this.X = x;
            this.U_x = u_x;
        }
        public double X { get; set; }
        public double U_x { get; set; }
    }

    public class NumericalResultErrorItem
    {
        public NumericalResultErrorItem(double xCenter, double u_xCenter, double u_xDeriv, double errorNorm, double relError)
        {
            this.XCenter = xCenter;
            this.U_xCenter = u_xCenter;
            this.U_xDeriv = u_xDeriv;
            this.ErrorNorm = errorNorm;
            this.RelError = relError;
        }

        public double XCenter { get; set; }
        public double U_xCenter { get; set; }
        public double U_xDeriv { get; set; }
        public double ErrorNorm { get; set; }
        public double RelError { get; set; }
    }
}
