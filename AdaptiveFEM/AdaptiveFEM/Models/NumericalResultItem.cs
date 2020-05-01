namespace AdaptiveFEM.Models
{
    public class NumericalResultItem
    {
        public double X { get; set; }
        public double Ux { get; set; }

        public NumericalResultItem(double x, double ux)
        {
            X = x;
            Ux = ux;
        }
    }
}
