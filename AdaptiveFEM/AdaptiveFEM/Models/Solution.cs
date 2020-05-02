namespace AdaptiveFEM.Models
{
    public class Solution
    {
        public double X { get; set; }
        public double Ux { get; set; }

        public Solution(double x, double ux)
        {
            X = x;
            Ux = ux;
        }
    }
}
