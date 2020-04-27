﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra.Double.Solvers.Preconditioners;


namespace MathNet.Numerics.LinearAlgebra.Double.Solvers.Iterative
{
    /*
    #region Copyright ©2004 Joannes Vermorel
    //http://code.google.com/p/paz-search/source/browse/MathNet/src/LinearAlgebra/Sparse/Linear/GMRESSolver.cs
    // MathNet Numerics, part of MathNet
    //
    // Copyright (c) 2004,﻿  Joannes Vermorel, http://www.vermorel.com
    // Based on JMP , Copyright (c) 2003 Bjørn-Ove Heimsund
    //
    // This program is free software; you can redistribute it and/or modify
    // it under the terms of the GNU Lesser General Public License as published
    // by the Free Software Foundation; either version 2 of the License, or
    // (at your option) any later version.
    //
    // This program is distributed in the hope that it will be useful,
    // but WITHOUT ANY WARRANTY; without even the implied warranty of
    // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    // GNU Lesser General Public License for more details.
    //
    // You should have received a copy of the GNU Lesser General Public
    // License along with this program; if not, write to the Free Software
    // Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    #endregion


  /// <summary> GMRES solver.
﻿  /// GMRES solves the unsymmetric linear system <c>Ax = b</c> using the
﻿  /// Generalized Minimum Residual method. The GMRES iteration is restarted after
﻿  /// a given number of iterations. By default it is restarted after 30 iterations.
﻿  /// </summary>
﻿  /// <author>  Templates </author>
   
﻿  public class GMRESSolver : IIterativeSolver
﻿  {
       
       /// <summary>
       /// The preconditioner that will be used. Can be set to <see langword="null" />, in which case the default
       /// pre-conditioner will be used.
       /// </summary>
       private IPreConditioner _preconditioner;
       /// <summary>
       /// Sets the <see cref="IPreConditioner"/> that will be used to precondition the iterative process.
       /// </summary>
       /// <param name="preconditioner">The preconditioner.</param>
       public void SetPreconditioner(IPreConditioner preconditioner)
       {
           _preconditioner = preconditioner;
       }


       /// <summary>
       /// The iterative process controller.
       /// </summary>
       private IIterator _iterator;

       /// <summary>
       /// Sets the <see cref="IIterator"/> that will be used to track the iterative process.
       /// </summary>
       /// <param name="iterator">The iterator.</param>
       public void SetIterator(IIterator iterator)
       {
           _iterator = iterator;
       }


﻿  ﻿  /// <summary> Sets restart to use</summary>
﻿  ﻿  /// <remark>GMRES iteration is restarted after this number of
﻿  ﻿  /// iterations </remark>
﻿  ﻿  public virtual int Restart
﻿  ﻿  {
﻿  ﻿  ﻿  set { this.restart = value; }

﻿  ﻿  }

﻿  ﻿  /// <summary> After this many iterations, the GMRES will be restarted.</summary>
﻿  ﻿  private int restart;

﻿  ﻿  /// <summary> Default restart of 30</summary>
﻿  ﻿  public GMRESSolver()
          : base()
﻿  ﻿  {
﻿  ﻿  ﻿  restart = 30;
﻿  ﻿  }
       
﻿  ﻿  public void Solve(Matrix A, Vector b, Vector x)
﻿  ﻿  {
  ﻿  ﻿  Vector w = new SparseVector(b.Count), temp = new SparseVector(b.Count), r = new SparseVector(b.Count);

﻿  ﻿  ﻿  double[] s = new double[restart + 1], cs = new double[restart + 1], sn = new double[restart + 1];

﻿  ﻿  ﻿  double[,] H = new double[restart + 1, restart];
﻿  ﻿  ﻿  for (int i = 0; i < restart + 1; i++)
             //﻿  ﻿  ﻿  {
             //﻿  ﻿  ﻿  ﻿  H[i] = new double[restart];
             //﻿  ﻿  ﻿  }

             //x.Add(v[j].Multiply(yy[j]));  

﻿  ﻿  ﻿  //temp = Blas.Default.MultAdd(-1.0, A, x, b, temp);
             //z = alpha*A*x + y
  ﻿  ﻿  temp = A.Multiply(x).Multiply(-1.0).Add(b) as SparseVector;
          //﻿  ﻿  ﻿  r = M.Apply(A, temp, r);
         _preconditioner.Initialize(A);          
         _preconditioner.Approximate(temp, r);

       //Blas.Default.Norm(r, NORMS.NORM2);
﻿  ﻿  ﻿  double normr = r.Norm(2); 
//﻿  ﻿  ﻿  temp = M.Apply(A, b, temp);
         _preconditioner.Approximate(b, temp);

﻿  ﻿  ﻿  Vector[] v = createVectors(b, restart + 1);

//﻿  ﻿  ﻿  for (iter.Reset(); !iter.Converged(r, x); iter.MoveNext())
﻿  ﻿  ﻿  {
//﻿  ﻿  ﻿  ﻿  v[0] = Blas.Default.ScaleCopy(1.0 / normr, r, v[0]);
//﻿  ﻿  ﻿  ﻿  SupportClass.ArraySupport.Fill(s, 0.0);
//﻿  ﻿  ﻿  ﻿  s[0] = normr;

//﻿  ﻿  ﻿  ﻿  for (int i = 0; i < restart; i++, iter.MoveNext())
//﻿  ﻿  ﻿  ﻿  {
//﻿  ﻿  ﻿  ﻿  ﻿  temp = Blas.Default.Mult(A, v[i], temp);
//﻿  ﻿  ﻿  ﻿  ﻿  w = M.Apply(A, temp, w);
//﻿  ﻿  ﻿  ﻿  ﻿  for (int k = 0; k <= i; k++)
//﻿  ﻿  ﻿  ﻿  ﻿  {
//﻿  ﻿  ﻿  ﻿  ﻿  ﻿  H[k, i] = Blas.Default.Dot(w, v[k]);
//﻿  ﻿  ﻿  ﻿  ﻿  ﻿  w = Blas.Default.Add(-H[k, i], v[k], w);
//﻿  ﻿  ﻿  ﻿  ﻿  }
//﻿  ﻿  ﻿  ﻿  ﻿  H[i + 1, i] = Blas.Default.Norm(w, NORMS.NORM2);
//﻿  ﻿  ﻿  ﻿  ﻿  v[i + 1] = Blas.Default.ScaleCopy(1.0 / H[i + 1, i], w, v[i + 1]);

//﻿  ﻿  ﻿  ﻿  ﻿  for (int k = 0; k < i; k++)
//﻿  ﻿  ﻿  ﻿  ﻿  {
//﻿  ﻿  ﻿  ﻿  ﻿  ﻿  RotationData RD = new RotationData(this, H[k, i], H[k + 1, i], cs[k], sn[k]);
//﻿  ﻿  ﻿  ﻿  ﻿  ﻿  applyPlaneRotation(RD);
//﻿  ﻿  ﻿  ﻿  ﻿  ﻿  H[k, i] = RD.dx;
//﻿  ﻿  ﻿  ﻿  ﻿  ﻿  H[k + 1, i] = RD.dy;
//﻿  ﻿  ﻿  ﻿  ﻿  ﻿  cs[k] = RD.cs;
//﻿  ﻿  ﻿  ﻿  ﻿  ﻿  sn[k] = RD.sn;
//﻿  ﻿  ﻿  ﻿  ﻿  }

//﻿  ﻿  ﻿  ﻿  ﻿  RotationData RD2 = new RotationData(this, H[i, i], H[i + 1, i], cs[i], sn[i]);
//﻿  ﻿  ﻿  ﻿  ﻿  generatePlaneRotation(RD2);
//﻿  ﻿  ﻿  ﻿  ﻿  H[i, i] = RD2.dx;
//﻿  ﻿  ﻿  ﻿  ﻿  H[i + 1, i] = RD2.dy;
//﻿  ﻿  ﻿  ﻿  ﻿  cs[i] = RD2.cs;
//﻿  ﻿  ﻿  ﻿  ﻿  sn[i] = RD2.sn;

//﻿  ﻿  ﻿  ﻿  ﻿  applyPlaneRotation(RD2);
//﻿  ﻿  ﻿  ﻿  ﻿  H[i, i] = RD2.dx;
//﻿  ﻿  ﻿  ﻿  ﻿  H[i + 1, i] = RD2.dy;
//﻿  ﻿  ﻿  ﻿  ﻿  cs[i] = RD2.cs;
//﻿  ﻿  ﻿  ﻿  ﻿  sn[i] = RD2.sn;

//﻿  ﻿  ﻿  ﻿  ﻿  RD2.dx = s[i];
//﻿  ﻿  ﻿  ﻿  ﻿  RD2.dy = s[i + 1];
//﻿  ﻿  ﻿  ﻿  ﻿  applyPlaneRotation(RD2);
//﻿  ﻿  ﻿  ﻿  ﻿  s[i] = RD2.dx;
//﻿  ﻿  ﻿  ﻿  ﻿  s[i + 1] = RD2.dy;
//﻿  ﻿  ﻿  ﻿  ﻿  cs[i] = RD2.cs;
//﻿  ﻿  ﻿  ﻿  ﻿  sn[i] = RD2.sn;

//﻿  ﻿  ﻿  ﻿  ﻿  if (iter.Converged(Math.Abs(s[i + 1]), x))
//﻿  ﻿  ﻿  ﻿  ﻿  {
//﻿  ﻿  ﻿  ﻿  ﻿  ﻿  update(x, i, s, v, H);
//﻿  ﻿  ﻿  ﻿  ﻿  ﻿  return;
//﻿  ﻿  ﻿  ﻿  ﻿  }
//﻿  ﻿  ﻿  ﻿  }

﻿  ﻿  ﻿  ﻿  update(x, restart - 1, s, v, H);
﻿  ﻿  ﻿  ﻿  temp = Blas.Default.MultAdd(-1.0, A, x, b, temp);
﻿  ﻿  ﻿  ﻿  r = M.Apply(A, temp, r);
﻿  ﻿  ﻿  ﻿  normr = Blas.Default.Norm(r, NORMS.NORM2);
﻿  ﻿  ﻿  }
﻿  ﻿  }

      private Vector[] createVectors(Vector template, int quantity)
      {
          List<Vector> entitys = new List<Vector>();
          for (int i = 0; i < quantity; i++)
          {
              entitys.Add(new Vector(template.Count));
          }
          return entitys.ToArray();
      }

﻿  ﻿  private void update(Vector x, int k, double[] s, Vector[] v, double[,] H)
﻿  ﻿  {
﻿  ﻿  ﻿  double[] yy = new double[s.Length];
﻿  ﻿  ﻿  Array.Copy(s, 0, yy, 0, s.Length);

﻿  ﻿  ﻿  // Backsolve
﻿  ﻿  ﻿  for (int i = k; i >= 0; i--)
﻿  ﻿  ﻿  {
﻿  ﻿  ﻿  ﻿  yy[i] /= H[i, i];
﻿  ﻿  ﻿  ﻿  for (int j = i - 1; j >= 0; j--)
﻿  ﻿  ﻿  ﻿  ﻿  yy[j] -= H[j, i] * yy[i];
﻿  ﻿  ﻿  }

﻿  ﻿  ﻿  for (int j = 0; j <= k; j++)
      {
          //y = alpha*x + y
//﻿  ﻿  ﻿  ﻿  x = Blas.Default.Add(yy[j], v[j], x);
          x.Add(v[j].Multiply(yy[j]));               
      }
﻿  ﻿  }

﻿  ﻿  /// <summary> Constructs a Givens rotation, to be used by applyPlaneRotation</summary>
﻿  ﻿  private void generatePlaneRotation(RotationData RD)
﻿  ﻿  {
﻿  ﻿  ﻿  if (RD.dy == 0.0)
﻿  ﻿  ﻿  {
﻿  ﻿  ﻿  ﻿  RD.cs = 1.0;
﻿  ﻿  ﻿  ﻿  RD.sn = 0.0;
﻿  ﻿  ﻿  }
﻿  ﻿  ﻿  else if (Math.Abs(RD.dy) > Math.Abs(RD.dx))
﻿  ﻿  ﻿  {
﻿  ﻿  ﻿  ﻿  double temp = RD.dx / RD.dy;
﻿  ﻿  ﻿  ﻿  RD.sn = 1.0 / Math.Sqrt(1.0 + temp * temp);
﻿  ﻿  ﻿  ﻿  RD.cs = temp * RD.sn;
﻿  ﻿  ﻿  }
﻿  ﻿  ﻿  else
﻿  ﻿  ﻿  {
﻿  ﻿  ﻿  ﻿  double temp = RD.dy / RD.dx;
﻿  ﻿  ﻿  ﻿  RD.cs = 1.0 / Math.Sqrt(1.0 + temp * temp);
﻿  ﻿  ﻿  ﻿  RD.sn = temp * RD.cs;
﻿  ﻿  ﻿  }
﻿  ﻿  }

﻿  ﻿  /// <summary> Applies a Givens plane rotation</summary>
﻿  ﻿  private void applyPlaneRotation(RotationData RD)
﻿  ﻿  {
﻿  ﻿  ﻿  double temp = RD.cs * RD.dx + RD.sn * RD.dy;
﻿  ﻿  ﻿  RD.dy = (-RD.sn) * RD.dx + RD.cs * RD.dy;
﻿  ﻿  ﻿  RD.dx = temp;
﻿  ﻿  }

﻿  ﻿  /// <summary> Data for a Givens rotation in GMRES.</summary>
﻿  ﻿  private class RotationData
﻿  ﻿  {
﻿  ﻿  ﻿  private GMRESSolver gmresSolver;

﻿  ﻿  ﻿  public double dx, dy, cs, sn;

﻿  ﻿  ﻿  private void InitBlock(GMRESSolver gmresSolver)
﻿  ﻿  ﻿  {
﻿  ﻿  ﻿  ﻿  this.gmresSolver = gmresSolver;
﻿  ﻿  ﻿  }

﻿  ﻿  ﻿  public GMRESSolver GmresSolver
﻿  ﻿  ﻿  {
﻿  ﻿  ﻿  ﻿  get { return gmresSolver; }

﻿  ﻿  ﻿  }

﻿  ﻿  ﻿  public RotationData(GMRESSolver enclosingInstance, double dx, double dy, double cs, double sn)
﻿  ﻿  ﻿  {
﻿  ﻿  ﻿  ﻿  InitBlock(enclosingInstance);
﻿  ﻿  ﻿  ﻿  this.dx = dx;
﻿  ﻿  ﻿  ﻿  this.dy = dy;
﻿  ﻿  ﻿  ﻿  this.cs = cs;
﻿  ﻿  ﻿  ﻿  this.sn = sn;
﻿  ﻿  ﻿  }

﻿  ﻿  ﻿  public RotationData(GMRESSolver enclosingInstance)
﻿  ﻿  ﻿  {
﻿  ﻿  ﻿  ﻿  InitBlock(enclosingInstance);
﻿  ﻿  ﻿  ﻿  this.dx = 0.0;
﻿  ﻿  ﻿  ﻿  this.dy = 0.0;
﻿  ﻿  ﻿  ﻿  this.cs = 0.0;
﻿  ﻿  ﻿  ﻿  this.sn = 0.0;
﻿  ﻿  ﻿  }
﻿  ﻿  }
﻿  
public void StopSolve()
      {
          throw new System.NotImplementedException();
      }

//public void SetIterator(IIterator iterator)
//{
//    throw new System.NotImplementedException();
//}

public Generic.Solvers.Status.ICalculationStatus IterationResult
{
    get { throw new System.NotImplementedException(); }
}

public Vector Solve(Matrix matrix, Vector vector)
{
    throw new System.NotImplementedException();
}

public Matrix Solve(Matrix matrix, Matrix input)
{
    throw new System.NotImplementedException();
}

public void Solve(Matrix matrix, Matrix input, Matrix result)
{
    throw new System.NotImplementedException();
}
   }
*/
}

