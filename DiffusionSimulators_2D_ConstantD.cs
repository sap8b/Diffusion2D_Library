using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using static Diffusion2D_Library.BoundaryCondition;

namespace Diffusion2D_Library
{
    /// <summary>
    /// Solves the parabolic partial differential equation: ∂u/∂t-(∂^2 u)/(∂x^2) - (∂^2 u)/(∂y^2)=f(x,y,t) in two dimensions
    /// </summary>
    public class DiffusionSimulators_2D_ConstantD : DiffusionSimulator2D
    {        
        // Constructors
        public DiffusionSimulators_2D_ConstantD(double D, double dx, double dy, int nx, int ny, double dt, int nt,
            string[] Boundary_Conditions, Del_BC_xy[] bc_s, Del_IC_xy I0, Del_Source_xy g)
        {
            d_coeff = D;
            this.dx = dx;
            this.dy = dy;
            this.dt = dt;
            cf_2D = new CompositionField2D(ny, nx);

            for (int i = 0; i < ny; i++)
            {
                for (int j = 0; j < nx; j++)
                {
                    cf_2D.XPositionValues[j, i] = i * dx;
                    cf_2D.YPositionValues[j, i] = j * dy;
                }
            }

            C_Initial = I0(cf_2D.XPositionValues, cf_2D.YPositionValues);
            this.I0 = I0;

            int num_bounds = Boundary_Conditions.Length;
            border_with_function = new BoundaryWithFunction[num_bounds];
            RVector C0;
            for (int i = 0; i < num_bounds; i++)
            {
                border_with_function[i] = new BoundaryWithFunction
                {
                    BoundaryLocation = i switch
                    {
                        0 => BoundingBox.top,
                        1 => BoundingBox.right,
                        2 => BoundingBox.left,
                        3 => BoundingBox.bottom,
                        _ => BoundingBox.bottom,
                    },
                    TypeBC = ConvertStringToEnumBC(Boundary_Conditions[i]),
                    BoundaryFunction = bc_s[i]
                };
                switch (border_with_function[i].BoundaryLocation)
                {
                    case BoundingBox.top:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[X.GetnRows - 1, 0];
                            C0 = border_with_function[i].BoundaryFunction(0, dt, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(ny - 1) + C0;
                            C_Initial.ReplaceRow(C0, ny - 1);
                        }
                        break;
                    case BoundingBox.right:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, Y.GetnCols - 1];
                            C0 = border_with_function[i].BoundaryFunction(0, dt, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(nx - 1) + C0;
                            C_Initial.ReplaceCol(C0, nx - 1);
                        }
                        break;
                    case BoundingBox.left:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0, dt, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(0) + C0;
                            C_Initial.ReplaceCol(C0, 0);
                        }
                        break;
                    case BoundingBox.bottom:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0, dt, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(0) + C0;
                            C_Initial.ReplaceRow(C0, 0);
                        }
                        break;
                }
            }

            gxt_xy = g;
            Chat_mode = Mode.quiet;
        }

        /// <summary>
        /// Constructor for the 2D diffusion simulator with a constant diffusivity
        /// </summary>
        /// <param name="D"></param>
        /// <param name="dx"></param>
        /// <param name="dy"></param>
        /// <param name="nx"></param>
        /// <param name="ny"></param>
        /// <param name="dt"></param>
        /// <param name="nt"></param>
        /// <param name="Boundary_Conditions"></param>
        /// <param name="bc_s"></param>
        /// <param name="I0"></param>
        /// <param name="g"></param>
        /// <param name="chat">Flag for whether or not the time-steps are displayed on the screen</param>
        /// <param name="base_filename">Base filename for output filenaming</param>
        public DiffusionSimulators_2D_ConstantD(double D, double dx, double dy, int nx, int ny, double dt, int nt,
            string[] Boundary_Conditions, Del_BC_xy[] bc_s, Del_IC_xy I0, Del_Source_xy g, Mode chat, string base_filename)
        {
            d_coeff = D;
            this.dx = dx;
            this.dy = dy;
            this.dt = dt;
            cf_2D = new CompositionField2D(ny, nx);

            for (int i = 0; i < ny; i++)
            {
                for (int j = 0; j < nx; j++)
                {
                    cf_2D.XPositionValues[j, i] = i * dx;
                    cf_2D.YPositionValues[j, i] = j * dy;
                }
            }

            C_Initial = I0(cf_2D.XPositionValues, cf_2D.YPositionValues);

            int num_bounds = Boundary_Conditions.Length;
            border_with_function = new BoundaryWithFunction[num_bounds];
            RVector C0;
            for (int i = 0; i < num_bounds; i++)
            {
                border_with_function[i] = new BoundaryWithFunction
                {
                    BoundaryLocation = i switch
                    {
                        0 => BoundingBox.top,
                        1 => BoundingBox.right,
                        2 => BoundingBox.left,
                        3 => BoundingBox.bottom,
                        _ => BoundingBox.bottom,
                    },
                    TypeBC = ConvertStringToEnumBC(Boundary_Conditions[i]),
                    BoundaryFunction = bc_s[i]
                };
                switch (border_with_function[i].BoundaryLocation)
                {
                    case BoundingBox.top:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[X.GetnRows - 1, 0];
                            C0 = border_with_function[i].BoundaryFunction(0, dt, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(ny - 1) + C0;
                            C_Initial.ReplaceRow(C0, ny - 1);
                        }
                        break;
                    case BoundingBox.right:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, Y.GetnCols - 1];
                            C0 = border_with_function[i].BoundaryFunction(0, dt, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(nx - 1) + C0;
                            C_Initial.ReplaceCol(C0, nx - 1);
                        }
                        break;
                    case BoundingBox.left:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0, dt, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(0) + C0;
                            C_Initial.ReplaceCol(C0, 0);
                        }
                        break;
                    case BoundingBox.bottom:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0, dt, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(0) + C0;
                            C_Initial.ReplaceRow(C0, 0);
                        }
                        break;
                }
            }
            this.I0 = I0;
            gxt_xy = g;
            Chat_mode = chat;
            Errors = new string[nt];
        }

        public DiffusionSimulators_2D_ConstantD(double[] coeffs, int[] n, Del_IC_xy I0)  //string Boundary_Conditions,
        {
            if (coeffs.Length >= 3)
            {
                d_coeff = coeffs[0];
                dx = coeffs[1];
                dt = coeffs[2];
            }
            if (n.Length == 2)
            {
                cf_2D = new CompositionField2D(n[0], n[1]);
            }

            //this.Boundary_Conditions = Boundary_Conditions;
            C_Initial = I0(cf_2D.XPositionValues, cf_2D.YPositionValues);
        }

        // Solvers
        /// <summary>
        /// Method for solving the 2D diffusion equation using the Alternating Direction Implicit algorithm
        /// </summary>
        /// <param name="n_time_steps">Integer specifying the number of time-steps to take during the simulation</param>
        public void Solve(int n_time_steps, int output_interval)
        {
            int nrows = cf_2D.InitialCompositionValues.GetnRows;
            int ncols = cf_2D.InitialCompositionValues.GetnCols;

            RMatrix C_Ex = new(nrows, ncols);
            RMatrix C_Im1 = new(nrows, ncols);
            RMatrix C_Im2 = new(nrows, ncols);

            RMatrix fn = new(nrows, ncols);
            RMatrix f0 = new(nrows, ncols);
            RMatrix f12 = new(nrows, ncols);

            RVector CL = new(nrows);
            RVector CR = new(nrows);
            RVector CT = new(ncols);
            RVector CB = new(ncols);

            double nu = d_coeff * dt / (2 * Math.Pow(dx, 2));

            // Define the A matrix for the explicit steps
            double off_d_val = nu;
            double diag_val = 1 - (2 * nu);
            TridiagonalMatrix A = new(ncols, diag_val, off_d_val, off_d_val);

            // Define the B matrices for the implicit steps
            off_d_val = -nu;
            diag_val = 1 + (2 * nu);
            TridiagonalMatrix B = new(nrows, diag_val, off_d_val, off_d_val);

            // Define the BC Matrices for the implicit half time-step
            double tau = dt; // / 2.0;
            off_d_val = -nu; // -tau;
            diag_val = 1.0 - (2.0 * nu); // 1 - (2 * tau);
            TridiagonalMatrix Bm = new(nrows, diag_val, off_d_val, off_d_val);

            off_d_val = nu; // tau;
            diag_val = 1.0 + (2.0 * nu);// 1 + (2 * tau);
            TridiagonalMatrix Bp = new(nrows, diag_val, off_d_val, off_d_val);

            string full_file_name;
            // Time evolution           
            for (int t = 0; t < n_time_steps; t++)
            {
                try
                {
                    if (Chat_mode == Mode.verbose && output_interval > 0 && Base_filename != null)
                    {
                        if (t % output_interval == 0)
                        {
                            decimal time = (decimal)(t * dt);
                            Console.ForegroundColor = ConsoleColor.Green;
                            Console.WriteLine("{0}s have been simulated", time);
                            full_file_name = Base_filename + time.ToString() + suffix;
                            if (File.Exists(full_file_name)) { File.Delete(full_file_name); }
                            if (t == 0) { FileWriteConcentrationData_CSV(full_file_name, X, Y, C_Initial); }
                            else { FileWriteConcentrationData_CSV(full_file_name, X, Y, C_Im2); }
                        }
                    }
                    // ===================
                    // Explicit time-step
                    // ===================
                    // 0 = top, 1 = right, 2 = left, 3 = bottom
                    switch (BCs_Functions[1].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            CR = BCs_Functions[1].BoundaryFunction(t, dt, BCs_Functions[1].PositionVaries, BCs_Functions[1].PositionFixed);
                            break;
                        case ABoundaryCondition.neumann:
                            RVector C1, C2; //,C3;
                            double cUpper1, cUpper2, cLower1, cLower2;
                            if (t == 0)
                            {
                                C1 = C_Initial.GetColVector(ncols - 1);
                                C2 = C_Initial.GetColVector(ncols - 2);

                                cLower1 = C_Initial[0, ncols - 1];
                                cLower2 = C_Initial[1, ncols - 2];
                                cUpper1 = C_Initial[nrows - 1, ncols - 1];
                                cUpper2 = C_Initial[nrows - 2, ncols - 2];
                            }
                            else
                            {
                                C1 = C_Im2.GetColVector(ncols - 1);
                                C2 = C_Im2.GetColVector(ncols - 2);

                                cLower1 = C_Im2[0, ncols - 1];
                                cLower2 = C_Im2[1, ncols - 2];
                                cUpper1 = C_Im2[nrows - 1, ncols - 1];
                                cUpper2 = C_Im2[nrows - 2, ncols - 2];
                            }
                            //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                            CR = (-2 * nu * C1) + (2 * nu * C2); // + C3
                            CR[0] = (-2 * nu * cLower1) + (2 * nu * cLower2);
                            CR[ncols - 1] = (-2 * nu * cUpper1) + (2 * nu * cUpper2);
                            break;
                        default:
                            break;
                    }
                    switch (BCs_Functions[2].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            CL = BCs_Functions[2].BoundaryFunction(t, dt, BCs_Functions[2].PositionVaries, BCs_Functions[2].PositionFixed);
                            break;
                        case ABoundaryCondition.neumann:
                            RVector C1, C2; //,C3;
                            double cUpper1, cUpper2, cLower1, cLower2;
                            if (t == 0)
                            {
                                C1 = C_Initial.GetColVector(0);
                                C2 = C_Initial.GetColVector(1);

                                cLower1 = C_Initial[0, 0];
                                cLower2 = C_Initial[1, 1];
                                cUpper1 = C_Initial[nrows - 1, 0];
                                cUpper2 = C_Initial[nrows - 2, 1];
                            }
                            else
                            {
                                C1 = C_Im2.GetColVector(0);
                                C2 = C_Im2.GetColVector(1);

                                cLower1 = C_Im2[0, 0];
                                cLower2 = C_Im2[1, 1];
                                cUpper1 = C_Im2[nrows - 1, 0];
                                cUpper2 = C_Im2[nrows - 2, 1];
                            }
                            //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                            CL = (-2 * nu * C1) + (2 * nu * C2); // + C3
                            CL[0] = (-2 * nu * cLower1) + (2 * nu * cLower2);
                            CL[ncols - 1] = (-2 * nu * cUpper1) + (2 * nu * cUpper2);
                            break;
                        default:
                            break;
                    }
                    for (int i = 1; i < nrows - 1; i++)
                    {
                        RVector xold;
                        if (t == 0) { xold = C_Initial.GetRowVector(i); }
                        else { xold = C_Im2.GetRowVector(i); } // C_Im2.GetRowVector(i); }

                        RVector v1 = A.Dot(xold);
                        C_Ex.ReplaceRow(v1, i);
                        C_Ex[i, 0] = CL[i]; //nu *
                        C_Ex[i, ncols - 1] = CR[i]; //nu *  
                    }

                    // ===================
                    // ===================
                    // One-half implicit time-step
                    // ===================
                    // Source terms
                    double t12 = (t + 0.5) * dt;
                    if (t == 0)
                    {
                        //f0 = gxt(X, Y, t * dt, D, C_Initial);
                        //fn = gxt(X, Y, (t + 1) * dt, D, C_Initial);
                        fn = gxt_xy(X, Y, t12, d_coeff, C_Initial);
                    }
                    else
                    {
                        //f0 = gxt(X, Y, t * dt, D, C_Im2);
                        //fn = gxt(X, Y, (t + 1) * dt, D, C_Im2);
                        fn = gxt_xy(X, Y, t12, d_coeff, C_Im2);
                    }

                    //f12 = (tau / 2.0) * (fn + f0);
                    f12 = tau / 2.0 * fn;

                    // BCs
                    switch (BCs_Functions[0].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            RVector gn = BCs_Functions[0].BoundaryFunction((t + 1), dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
                            RVector g0 = BCs_Functions[0].BoundaryFunction(t, dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
                            CT = (0.5 * B.Dot(gn)) + (0.5 * A.Dot(g0)); ////((0.5 * Bm.Dot(gn)) + (0.5 * Bp.Dot(g0))); //nu * 
                            break;
                        case ABoundaryCondition.neumann:
                            RVector C1, C2; //,C3;
                            if (t == 0)
                            {
                                C1 = C_Initial.GetRowVector(nrows - 1);
                                C2 = C_Initial.GetRowVector(nrows - 2);
                            }
                            else
                            {
                                C1 = C_Im2.GetRowVector(nrows - 1);
                                C2 = C_Im2.GetRowVector(nrows - 2);
                            }
                            //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                            CT = (-2 * nu * C1) + (2 * nu * C2); // + C3
                            break;
                        default:
                            break;
                    }
                    switch (BCs_Functions[3].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            RVector gn = BCs_Functions[3].BoundaryFunction((t + 1), dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
                            RVector g0 = BCs_Functions[3].BoundaryFunction(t, dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
                            CB = (0.5 * B.Dot(gn)) + (0.5 * A.Dot(g0)); //((0.5 * Bm.Dot(gn)) + (0.5 * Bp.Dot(g0))); //nu * 
                            break;
                        case ABoundaryCondition.neumann:
                            RVector C1, C2; //,C3;
                            if (t == 0)
                            {
                                C1 = C_Initial.GetRowVector(0);
                                C2 = C_Initial.GetRowVector(1);
                            }
                            else
                            {
                                C1 = C_Im2.GetRowVector(0);
                                C2 = C_Im2.GetRowVector(1);
                            }
                            //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                            CB = (-2 * nu * C1) + (2 * nu * C2); // + C3
                            break;
                        default:
                            break;
                    }
                    for (int j = 1; j < ncols - 1; j++)
                    {
                        RVector v1 = C_Ex.GetColVector(j);
                        v1[0] = CB[j]; //nu * 
                        v1[ncols - 1] = CT[j]; //nu * 

                        RVector f12s = f12.GetColVector(j); // 

                        RVector u12 = TridiagonalMatrix.Thomas_Algorithm(B, v1 + f12s);
                        //u12[0] = CB[j]; //nu * 
                        //u12[ncols - 1] = CT[j]; //nu * 
                        C_Im1.ReplaceCol(u12, j);
                    }

                    // ===================
                    // ===================
                    // Full implicit time-step
                    // ===================
                    switch (BCs_Functions[1].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            //CR = BCs_Functions[1].BoundaryFunction((t + 1) * dt, BCs_Functions[1].PositionVaries, BCs_Functions[1].PositionFixed);
                            CT = BCs_Functions[1].BoundaryFunction((t + 1), dt, BCs_Functions[1].PositionVaries, BCs_Functions[0].PositionFixed);
                            RVector ctab = C_Initial.GetRowVector(0);
                            C_Im2.ReplaceRow(CT, nrows - 1);
                            break;
                        case ABoundaryCondition.neumann:
                            RVector C1, C2; //,C3;
                            if (t == 0)
                            {
                                C1 = C_Initial.GetColVector(ncols - 1);
                                C2 = C_Initial.GetColVector(ncols - 2);
                            }
                            else
                            {
                                C1 = C_Im2.GetColVector(ncols - 1);
                                C2 = C_Im2.GetColVector(ncols - 2);
                            }
                            //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                            CR = (-2 * nu * C1) + (2 * nu * C2); // + C3
                            break;
                        default:
                            break;
                    }
                    switch (BCs_Functions[2].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            //CL = BCs_Functions[2].BoundaryFunction((t + 1) * dt, BCs_Functions[2].PositionVaries, BCs_Functions[2].PositionFixed);
                            CB = BCs_Functions[2].BoundaryFunction((t + 1), dt, BCs_Functions[2].PositionVaries, BCs_Functions[3].PositionFixed);

                            RVector ctab = C_Initial.GetRowVector(nrows - 1);

                            C_Im2.ReplaceRow(CB, 0);
                            break;
                        case ABoundaryCondition.neumann:
                            RVector C1, C2; //,C3;
                            if (t == 0)
                            {
                                C1 = C_Initial.GetColVector(0);
                                C2 = C_Initial.GetColVector(1);
                            }
                            else
                            {
                                C1 = C_Im2.GetColVector(0);
                                C2 = C_Im2.GetColVector(1);
                            }
                            //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                            CL = (-2 * nu * C1) + (2 * nu * C2); // + C3
                            break;
                        default:
                            break;
                    }
                    for (int k = 1; k < nrows - 1; k++)
                    {
                        RVector v1 = C_Ex.GetRowVector(k);
                        RVector u12 = C_Im1.GetRowVector(k);
                        RVector b = (2 * u12) - v1;
                        b[0] = CL[k]; //nu * 
                        b[ncols - 1] = CR[k]; //nu * 

                        RVector u1 = TridiagonalMatrix.Thomas_Algorithm(B, b);
                        u1[0] = CL[k];  //nu * 
                        u1[ncols - 1] = CR[k]; //nu * 
                        C_Im2.ReplaceRow(u1, k);
                    }
                    // ===================
                }
                catch (Exception e)
                {
                    Errors[t] = e.Message;
                    error_flag = true;
                    Console.ForegroundColor = ConsoleColor.Red;
                    Console.WriteLine("Fault!");
                }
            }

            // Setup the return composition field
            if (error_flag == false) { C_Final = C_Im2; }
            else { C_Final = C_Initial; }
        }


        //
        // =======================================================================
        // =======================================================================
        // =======================================================================
    }
}
