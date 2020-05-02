namespace AdaptiveFEM.Models
{
    public class Element
    {
        public double Begin { get; set; }
        public double End { get; set; }
        public double MidPoint { get; set; }
        public double H { get; set; }

        public Element(double begin, double end, double midPoint, double h)
        {
            Begin = begin;
            End = end;
            MidPoint = midPoint;
            H = h;
        }
    }
}
