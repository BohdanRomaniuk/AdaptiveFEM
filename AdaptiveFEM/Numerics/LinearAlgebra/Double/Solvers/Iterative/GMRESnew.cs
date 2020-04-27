using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;
using MathNet.Numerics.LinearAlgebra.Generic.Solvers.Status;
using MathNet.Numerics.LinearAlgebra.Double.Solvers.Preconditioners;

namespace GMRES
{
    /// <summary>
    /// Generalized Minimal Resilual iterative matrix solver.
    /// The gmres solver can be used to non-symetric matrices.
    /// Much of the success depends on structure of matrix and proper
    /// preconditioner.
    /// The gmres algorithm was taken from:
    /// GMRES(spaceDimention) : Parallel scientific computing in C++ and MPI
    /// George Em Karniadakis
    /// Robert space.Kirby II
    /// page 590-597 
    /// </summary>
    public sealed class GMRES : IIterativeSolver
    {
        /// <summary>
        /// Dimension of Krylov space,
        /// By default is 3
        /// </summary>
        private int spaceDimention = 50;

        /// <summary>
        /// The iterative process controller.
        /// </summary>
        private IIterator _iterator;

        /// <summary>
        /// The preconditioner that will be used. Can be set to <c>null</c>, in which case the default
        /// pre-conditioner will be used.
        /// </summary>
        private IPreConditioner _preconditioner;

        /// <summary>
        /// Indicates if the user has stopped the solver.
        /// </summary>
        private bool _hasBeenStopped = false;

        /// <summary>
        /// The status used if there is no status, i.e. the solver hasn't run yet and there is no
        /// iterator.
        /// </summary>
        private static readonly ICalculationStatus DefaultStatus = new CalculationIndetermined();

        /// <summary>
        /// Initializes a new instance of the <see cref="GMRES"/> class.
        /// </summary>
        /// <remarks>
        /// When using this constructor the solver will use the <see cref="IIterator"/> with
        /// the standard settings and a default preconditioner.
        /// </remarks>
        public GMRES()
            : this(null, null)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="GMRES"/> class.
        /// </summary>
        /// <remarks>
        /// <para>
        /// When using this constructor the solver will use a default preconditioner.
        /// </para>
        /// <para>
        /// The main advantages of using a user defined <see cref="IIterator"/> are:
        /// <list type="number">
        /// <item>It is possible to set the desired convergence limits.</item>
        /// <item>
        /// It is possible to check the reason for which the solver finished 
        /// the iterative procedure by calling the <see cref="IIterator.Status"/> property.
        /// </item>
        /// </list>
        /// </para>
        /// </remarks>
        /// <param name="iterator">The <see cref="IIterator"/> that will be used to monitor the iterative process. </param>
        public GMRES(IIterator iterator)
            : this(null, iterator)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="GMRES"/> class.
        /// </summary>
        /// <remarks>
        /// When using this constructor the solver will use the <see cref="IIterator"/> with
        /// the standard settings.
        /// </remarks>
        /// <param name="preconditioner">The <see cref="IPreConditioner"/> that will be used to precondition the matrix equation.</param>
        public GMRES(IPreConditioner preconditioner)
            : this(preconditioner, null)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="GMRES"/> class.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The main advantages of using a user defined <see cref="IIterator"/> are:
        /// <list type="number">
        /// <item>It is possible to set the desired convergence limits.</item>
        /// <item>
        /// It is possible to check the reason for which the solver finished 
        /// the iterative procedure by calling the <see cref="IIterator.Status"/> property.
        /// </item>
        /// </list>
        /// </para>
        /// </remarks>
        /// <param name="preconditioner">The <see cref="IPreConditioner"/> that will be used to precondition the matrix equation.</param>
        /// <param name="iterator">The <see cref="IIterator"/> that will be used to monitor the iterative process.</param>
        public GMRES(IPreConditioner preconditioner, IIterator iterator)
        {
            _iterator = iterator;
            _preconditioner = preconditioner;
        }

        /// <summary>
        /// Gets the status of the iteration once the calculation is finished.
        /// </summary>
        public ICalculationStatus IterationResult
        {
            get
            {
                return (_iterator != null) ? _iterator.Status : DefaultStatus;
            }
        }

        /// <summary>
        /// Sets the <see cref="IPreConditioner"/> that will be used to precondition the iterative process.
        /// </summary>
        /// <param name="preconditioner">The preconditioner.</param>
        public void SetPreconditioner(IPreConditioner preconditioner)
        {
            _preconditioner = preconditioner;
        }

        /// <summary>
        /// Sets the <see cref="IIterator"/> that will be used to track the iterative process.
        /// </summary>
        /// <param name="iterator">The iterator.</param>
        public void SetIterator(IIterator iterator)
        {
            _iterator = iterator;
        }

        /// <summary>
        /// Stops the solve process. 
        /// </summary>
        /// <remarks>
        /// Note that it may take an indetermined amount of time for the solver to actually
        /// stop the process.
        /// </remarks>
        public void StopSolve()
        {
            _hasBeenStopped = true;
        }

        /// <summary>
        /// Solves the matrix equation Ax = b, where A is the coefficient matrix, b is the
        /// solution vector and x is the unknown vector.
        /// </summary>
        /// <param name="matrix">The coefficient matrix, <c>A</c>.</param>
        /// <param name="vector">The solution vector, <c>b</c>.</param>
        /// <returns>The result vector, <c>x</c>.</returns>
        public Vector Solve(Matrix matrix, Vector vector)
        {
            return solve(matrix, vector);
        }

        /// <summary>
        /// Solves the matrix equation AX = B, where A is the coefficient matrix, B is the
        /// solution matrix and X is the unknown matrix.
        /// </summary>
        /// <param name="matrix">The coefficient matrix, <c>A</c>.</param>
        /// <param name="input">The solution matrix, <c>B</c>.</param>
        /// <returns>The result matrix, <c>X</c>.</returns>
        public Matrix Solve(Matrix matrix, Matrix input)
        {

            Matrix resultMatrix = new DenseMatrix(input.RowCount, input.ColumnCount);
            for (int i = 0; i < input.ColumnCount; i++)
            {
                resultMatrix.SetColumn(i, solve(matrix, (Vector)input.Column(i)));
            }
            return resultMatrix;
        }

        /// <summary>
        /// Solves the matrix equation Ax = b, where A is the coefficient matrix, b is the
        /// solution vector and x is the unknown vector.
        /// </summary>
        /// <param name="matrix">The coefficient matrix, <c>A</c>.</param>
        /// <param name="input">The solution vector, <c>b</c></param>
        /// <param name="result">The result vector, <c>x</c></param>
        public void Solve(Matrix matrix, Vector input, Vector result)
        {
            result = solve(matrix, input);

        }

        /// <summary>
        /// Solves the matrix equation AX = B, where A is the coefficient matrix, B is the
        /// solution matrix and X is the unknown matrix.
        /// </summary>
        /// <param name="matrix">The coefficient matrix, <c>A</c>.</param>
        /// <param name="input">The solution matrix, <c>B</c>.</param>
        /// <param name="result">The result matrix, <c>X</c></param>
        public void Solve(Matrix matrix, Matrix input, Matrix result)
        {
            for (int i = 0; i < input.ColumnCount; i++)
            {
                result.SetColumn(i, solve(matrix, (Vector)(input.Column(i))));
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="iterationNumber"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        private bool ShouldContinue(int iterationNumber, Vector result, Vector source, Vector residuals)
        {
            if (_hasBeenStopped)
            {
                _iterator.IterationCancelled();
                return true;
            }
            _iterator.DetermineStatus(iterationNumber, result, source, residuals);

            var status = _iterator.Status;

            // We stop if either:
            // - the user has stopped the calculation
            // - the calculation needs to be stopped from a numerical point of view (divergence, convergence etc.)
            return (!status.TerminatesCalculation) && (!_hasBeenStopped);
        }
        /// <summary>
        /// method to solve linear system A x = b
        /// </summary>
        /// <param name="A">matrix of system</param>
        /// <param name="b">right hand side vector</param>
        /// <param name="x">left hand side vector</param>
        /// <param name="eps">precision</param>
        /// <param name="spaceDimention">size of Krylov subspace</param>
        /// <param name="nsteps">maximum number of steps</param>
        /// <returns>vector - the result</returns>
        private Vector solve(Matrix A, Vector b)
        {
            // Error checks
            if (A == null)
            {
                throw new ArgumentNullException("matrix");
            }

            if (A.RowCount != A.ColumnCount)
            {
                throw new ArgumentException("matrix is not square");
            }

            if (b == null)
            {
                throw new ArgumentNullException("input");
            }

            if (b.Count != A.RowCount)
            {
                throw new ArgumentException("matrix row count is not equal to input dimension");
            }

            _hasBeenStopped = false;

            // Initialize the solver fields
            // Set the convergence monitor
            if (_iterator == null)
            {
                _iterator = Iterator.CreateDefault();
            }

            if (_preconditioner == null)
            {
                _preconditioner = new UnitPreconditioner();
                    //IncompleteLU();
                    //// new Diagonal();
                
            }

            _preconditioner.Initialize(A);

            Vector x = new DenseVector(A.RowCount);
            int iterationsNumber = 1;
            int n = A.RowCount;

            //the residual vector
            Vector r = new DenseVector(b);
            //matrix with orthonormal vectors
            Matrix V = new DenseMatrix(n, spaceDimention + 1);
            //Hessenberg matrix
            Matrix H = new DenseMatrix(spaceDimention + 1, spaceDimention);
            //coefficients for Givens rotation
            Vector cs = new DenseVector(spaceDimention + 1);
            Vector sn = new DenseVector(spaceDimention + 1);

            Vector s = new DenseVector(spaceDimention + 1);

            Vector xh = new DenseVector(spaceDimention + 1);
            Vector unity = new DenseVector(1);
            unity[0] = 1.0;
            Vector r1 = new DenseVector(1);
            double normb = b.Norm(2);

            //r -= Ax            
            //r.Add(A.Multiply(x.Negate()));
            r.Subtract(A.Multiply(x));
            // beta = norm(r)
            double beta = r.Norm(2);

            if (!ShouldContinue(iterationsNumber, x, b, r))
            {
                //CalculateResidual(beta, normb, resid);
                return x;
            }

            //to do : perform iterations of algorithm 
            while (ShouldContinue(iterationsNumber, x, b, r)) //j <= iterationsNumber
            {
                //V0 - first orthonormal vector
                //V = (DenseMatrix)V.InsertColumn(0, r.Normalize(2));
                V.SetColumn(0, r.Normalize(2));
                s.Clear();
                s[0] = beta;
                //to do : perform Arnoldi iterations
                for (int i = 0; i < spaceDimention; i++, iterationsNumber++)
                {
                    //to do :  w = A * v[i];
                    xh = (Vector)V.Column(i);
                    //precondition expected
                    //...
                    xh = _preconditioner.Approximate(xh);
                    //to do : V[i+1] = A xh                                     
                    V.SetColumn(i + 1, A.Multiply(xh));

                    //to do : generate matrix H
                    for (int k = 0; k <= i; k++)
                    {
                        //to do : H[k,i] = V[i+1]*V[i]
                        H[k, i] = V.Column(i + 1).DotProduct(V.Column(k));
                        //to do : V[i+1] -= H[k,i]*V[i]                    
                        V.SetColumn(i + 1, V.Column(i + 1) + V.Column(k).Multiply(-H[k, i]));
                    }
                    //to do : H[i+1,i] = ||V[i+1]||
                    H[i + 1, i] = V.Column(i + 1).Norm(2);
                    //to do : normalize V[i+1]
                    V.SetColumn(i + 1, V.Column(i + 1).Normalize(2));

                    //to do : apply old Givens rotations to the last column in H
                    for (int k = 0; k < i; k++)
                    {
                        double h1 = H[k, i], h2 = H[k + 1, i];
                        ApplyPlRot(ref h1, ref h2, cs[k], sn[k]);
                        H[k, i] = h1;
                        H[k + 1, i] = h2;
                    }
                    //to do :  generate new Givens rotation which eleminates H[i+1,i]                    
                    double c1 = cs[i], s1 = sn[i];
                    GenPlRot(H[i, i], H[i + 1, i], ref c1, ref s1);
                    cs[i] = c1;
                    sn[i] = s1;
                    // apply it to H and s                   
                    double h3 = H[i, i], h4 = H[i + 1, i];
                    ApplyPlRot(ref h3, ref h4, cs[i], sn[i]);
                    H[i, i] = h3;
                    H[i + 1, i] = h4;
                    double s3 = s[i], s4 = s[i + 1];
                    ApplyPlRot(ref s3, ref s4, cs[i], sn[i]);
                    s[i] = s3;
                    s[i + 1] = s4;

                    // b = {1,0,0,0,..};
                    // r = {fabs(s[i+1]/normb) ,0,0,0,....};
                    //if ((resid=fabs(s[i+1]/normb)) < eps)\\

                    r1[0] = Math.Abs(s[i + 1]) / normb;
                    if (!ShouldContinue(iterationsNumber, r1, unity, r1))
                    {
                        //update result x with size of Krylov subspace as i + 1
                        Update(A, i + 1, H, s, V, ref x);

                        return x;
                    }
                }

                //update result x with size of Krylov subspace as spaceDimention
                Update(A, spaceDimention, H, s, V, ref x);
                //to do : r = b - Ax
                r = (Vector)b.Add(A.Multiply(x.Negate()));
                beta = r.Norm(2);

                if (!ShouldContinue(iterationsNumber, x, b, r))
                {
                    //CalculateResidual(beta, normb, resid);
                    return x;
                }
            }
            return x;
        }

        /// <summary>
        /// method to generate coefficients for Givens rotation
        /// using first and second parameters
        /// and store the result in first and second coefficient
        /// </summary>
        /// <param name="dx">first parameter</param>
        /// <param name="dy">second parameter</param>
        /// <param name="cs">first coefficient</param>
        /// <param name="sn">second coefficient</param>
        private void GenPlRot(double dx, double dy, ref double cs, ref double sn)
        {
            if (dy == 0.0)
            {
                cs = 1.0;
                sn = 0.0;
            }
            else if (Math.Abs(dy) > Math.Abs(dx))
            {
                double tmp = dx / dy;
                sn = 1.0 / Math.Sqrt(1.0 + tmp * tmp);
                cs = tmp * sn;
            }
            else
            {
                double tmp = dy / dx;
                cs = 1.0 / Math.Sqrt(1.0 + tmp * tmp);
                sn = tmp * cs;
            }
        }

        /// <summary>
        /// method to apply Givens rotation
        /// using coefficients and parameters to apply
        /// </summary>
        /// <param name="dx">first param</param>
        /// <param name="dy">second param</param>
        /// <param name="cs">first coef</param>
        /// <param name="sn">second coef</param>
        private void ApplyPlRot(ref double dx, ref double dy, double cs, double sn)
        {
            double tmp = cs * dx + sn * dy;
            dy = cs * dy - sn * dx;
            dx = tmp;
        }

        /// <summary>
        /// method to update solution of linear system A x = b
        /// </summary>
        /// <param name="A">matrix of system</param>
        /// <param name="k">size of Krylov subspace</param>
        /// <param name="H">Hessenberg matrix</param>
        /// <param name="s">right hand side vector</param>
        /// <param name="V">matrix with orthonormal vectors</param>
        /// <param name="x">vector to update</param>
        private void Update(Matrix A, int k, Matrix H, Vector s, Matrix V, ref Vector x)
        {
            //to solve system H y = s
            Vector y = new DenseVector(k);
            s.CopyTo(y, 0, 0, k);
            y = BackSolve(H, y);

            Vector y1 = new DenseVector(V.ColumnCount);
            y.CopyTo(y1, 0, 0, k);

            Vector xh = new DenseVector(A.RowCount);
            //to do xh += Vy            
            xh = (Vector)V.Multiply(y1);

            //preconditioner
            xh = _preconditioner.Approximate(xh);
            //to do x += xh
            x = (Vector)x.Add(xh);
        }

        /// <summary>
        /// Method to solve linear system H x = r.
        /// </summary>
        /// <param name="H">Upper triangular matrix</param>
        /// <param name="r">Right hand side vector</param>
        /// <returns>Vector x - solution of linear system</returns>
        private Vector BackSolve(Matrix H, Vector r)
        {
            Vector x = new DenseVector(r);
            for (int i = x.Count - 1; i >= 0; i--)
            {
                for (int j = x.Count - 1; j > i; j--)
                {
                    x[i] -= H[i, j] * x[j];
                }
                x[i] /= H[i, i];
            }
            return x;
        }

        /// <summary>
        /// Property for restart parameter.
        /// </summary>
        public int RestartParameter
        {
            set
            {
                spaceDimention = value;
            }
        }
    }
}
