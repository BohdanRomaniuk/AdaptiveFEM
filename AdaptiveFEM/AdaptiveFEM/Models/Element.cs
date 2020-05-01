namespace AdaptiveFEM.Models
{
    public class Element
    {
        private double begin;
        private double end;
        private double midPoint;
        private double h;

        public Element(double begin, double end, double midPoint, double h)
        {
            this.begin = begin;
            this.end = end;
            this.midPoint = midPoint;
            this.h = h;
        }

        public double Begin { get { return this.begin; } }
        public double End { get { return this.end; } }
        public double MidPoint { get { return this.midPoint; } }
        public double H { get { return this.h; } }
    }
}
