using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using static Diffusion2D_Library.BoundaryCondition;
using static Diffusion2D_Library.DiffusionSimulators_1D;

namespace Diffusion2D_Library
{
    /// <summary>
    /// This class simulates diffusion in a 2D region with a non-constant diffusion coefficient 
    /// </summary>
   public class DiffusionSimulators_2D_MatrixD : DiffusionSimulator2D 
    {
        // Fields
        private readonly TridiagonalMatrix[] A_explicit;
        private readonly TridiagonalMatrix[] A_implicit;
        private readonly TridiagonalMatrix[] B_row;
        private readonly TridiagonalMatrix[] B_col;
        private double[] NeumannBCs_L_A;
        private double[] NeumannBCs_R_A;


        // Properties
        public RMatrix D
        {
            get => cf_2D.DValues;
            set { if (value.GetnCols > 0 && value.GetnRows > 0) { cf_2D.DValues = value; } }
        }

        // ====================================================================
        // Constructors
        // ====================================================================
        /// <summary>
        /// Constructs a 2D solver for the diffusion equation with a matrix diffusivity
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
        /// <param name="Tmode"></param>
        /// <param name="base_filename"></param>
        public DiffusionSimulators_2D_MatrixD(
            RMatrix D,
            double dx,
            double dy,
            int nx,
            int ny,
            double dt,
            int nt,
            ABoundaryCondition[] Boundary_Conditions,
            Del_BC_xy[] bc_s,
            Del_IC_xy I0,
            Del_Source_MatrixD g,
            Mode Tmode,
            string base_filename)
        {
            NX = nx;
            NY = ny;
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
                    cf_2D.DValues[j, i] = D[j, i];
                }
            }

            // ====================================
            // Set-up the initial condition function and start establishing the initial composition field
            // ===================================
            this.I0 = I0;
            C_Initial = I0(cf_2D.XPositionValues, cf_2D.YPositionValues);
            // ====================================
            // Set the boundary conditions and apply them to the initial composition field, if needed
            // ====================================
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
                    TypeBC = Boundary_Conditions[i],
                    BoundaryFunction = bc_s[i]
                };
                switch (border_with_function[i].BoundaryLocation)
                {
                    case BoundingBox.top:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[X.GetnRows - 1, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(ny - 1) + C0;
                            C_Initial.ReplaceRow(Ctab, ny - 1);
                        }
                        else if (border_with_function[i].TypeBC == ABoundaryCondition.neumann)
                        {
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[X.GetnRows - 1, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(ny - 1) + C0;
                            C_Initial.ReplaceRow(Ctab, ny - 1);
                        }
                        break;
                    case BoundingBox.right:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, Y.GetnCols - 1];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(nx - 1) + C0;
                            C_Initial.ReplaceCol(Ctab, nx - 1);
                        }
                        else if (border_with_function[i].TypeBC == ABoundaryCondition.neumann)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, Y.GetnCols - 1];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(nx - 1) + C0;
                            C_Initial.ReplaceCol(Ctab, nx - 1);
                        }
                        break;
                    case BoundingBox.left:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(0) + C0;
                            C_Initial.ReplaceCol(Ctab, 0);
                        }
                        else if (border_with_function[i].TypeBC == ABoundaryCondition.neumann)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(0) + C0;
                            C_Initial.ReplaceCol(Ctab, 0);
                        }
                        break;
                    case BoundingBox.bottom:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(0) + C0;
                            C_Initial.ReplaceRow(Ctab, 0);
                        }
                        else if (border_with_function[i].TypeBC == ABoundaryCondition.neumann)
                        {
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(0) + C0;
                            C_Initial.ReplaceRow(Ctab, 0);
                        }
                        break;
                }
            }
            // ====================================
            // Create the solver method matrices
            // ====================================
            nx_less1 = NX - 1;
            // Total points minus the left/right or up/down boundaries
            nx_less2 = NX - 2;
            // nx_less1 - one point because the total entries are one fewer because these are on the off-diagonals
            nx_less3 = nx_less2 - 1;

            end_idx1 = nx_less3;
            end_idx2 = nx_less3 - 1;
            start_idx1 = 0;
            start_idx2 = 1;

            A_explicit = new TridiagonalMatrix[nx];
            A_implicit = new TridiagonalMatrix[nx];
            B_row = new TridiagonalMatrix[nx];
            B_col = new TridiagonalMatrix[nx];

            NeumannBCs_L_A = new double[nx];
            NeumannBCs_R_A = new double[nx];

            double nu0 = dt / (2 * Math.Pow(dx, 2));

            //RVector nu = new(nx_less2);
            RVector off_d_val_l0 = new(NX-1);
            RVector off_d_val_u0 = new(NX - 1);
            RVector main0 = new(NX);
            RVector off_d_val_l1 = new(NX - 1);
            RVector off_d_val_u1 = new(NX - 1);
            RVector main1 = new(NX);
            RVector off_d_val_l2 = new(NX - 1);
            RVector off_d_val_u2 = new(NX - 1);
            RVector main2 = new(NX);

            double D_minus, D_plus;
            for (int i = 0; i < NX; i++)
            {
                RVector D_row = D.GetRowVector(i);
                double D0 = D_row[0];
                double D1 = D_row[1];
                double Dend = D_row[NX-1];
                double Dend1 = D_row[NX-2];

                // Working this for diffusion heterogeneous media
                //================================================
                // A Matrix first
                //================================================
                #region A_explicit
                D_plus = 0.5 * (D1 + D0);
                //main0[1] = 1.0 - 2.0 * D_plus * nu0;
                NeumannBCs_L_A[i] = 1.0 - D_plus * nu0;

                D_minus = 0.5 * (Dend1 + Dend);
                //main0[NX - 2] = 1.0 - 2.0 * D_minus * nu0;
                NeumannBCs_R_A[i] = 1.0 - D_minus * nu0;

                for (int j = 1; j < NX-1; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
                    main0[j] = 1.0 - ((D_minus * nu0) + (D_plus * nu0));
                }
                for (int j = 1; j < NX; j++)
                {
                    D_minus = 0.5 * (D_row[j] + D_row[j - 1]);
                    off_d_val_l0[j - 1] = nu0 * D_minus;
                }
                for (int j = 0; j < NX-1; j++)
                {
                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
                    off_d_val_u0[j] = nu0 * D_plus;
                }

                A_explicit[i] = new TridiagonalMatrix(main0, off_d_val_u0, off_d_val_l0);
                
                if (border_with_function[2].TypeBC == ABoundaryCondition.dirichlet)
                {
                    //D_plus = 0.5 * (D1 + D0);
                    A_explicit[i][0, 0] = 1.0; // D_plus * nu0;
                    A_explicit[i][0, 1] = 0.0;
                }
                else
                {
                    D_plus = 0.5 * (D1 + D0);
                    A_explicit[i][0, 0] = 1.0; // dx; //1.0 - dx * D_plus * nu0; // C0[0]
                    A_explicit[i][0, 1] = 0.0; // 1.0; // dx * D_plus * nu0; // C0[0]
                    A_explicit[i][1, 1] = 1.0 - D_plus * nu0; // C0[0]
                    //A_explicit[i][1, 0] = dx * D_plus * nu0; // C0[0]

                }
                if (border_with_function[1].TypeBC == ABoundaryCondition.dirichlet)
                {
                    //D_minus = 0.5 * (Dend1 + Dend);
                    A_explicit[i][NX - 1, NX - 1] = 1.0; // D_minus * nu0;
                    A_explicit[i][NX - 1, NX - 2] = 0.0;
                }
                else
                {
                    A_explicit[i][NX - 1, NX - 1] = 1.0; // dx; //1.0 - dx * D_plus * nu0; // C0[0]
                    A_explicit[i][NX - 1, NX - 2] = 0.0; // 1.0; // dx * D_plus * nu0; // C0[0]
                    A_explicit[i][NX - 2, NX - 2] = 1.0 - D_minus * nu0; // C0[0]
                    //A_explicit[i][NX - 2, NX - 1] = dx * D_minus * nu0; // C0[0]
                }
                #endregion
                //
                // A Implicit for the one-half time-step
                //
                #region A_implicit
                D_plus = 0.5 * (D1 + D0);
                //main0[1] = 1.0 - 2.0 * D_plus * nu0;
                NeumannBCs_L_A[i] = 1.0 - D_plus * nu0;

                D_minus = 0.5 * (Dend1 + Dend);
                //main0[NX - 2] = 1.0 - 2.0 * D_minus * nu0;
                NeumannBCs_R_A[i] = 1.0 - D_minus * nu0;

                for (int j = 1; j < NX - 1; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
                    main0[j] = 1.0 - ((D_minus * nu0) + (D_plus * nu0));
                }
                for (int j = 1; j < NX; j++)
                {
                    D_minus = 0.5 * (D_row[j] + D_row[j - 1]);
                    off_d_val_l0[j - 1] = nu0 * D_minus;
                }
                for (int j = 0; j < NX - 1; j++)
                {
                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
                    off_d_val_u0[j] = nu0 * D_plus;
                }

                A_implicit[i] = new TridiagonalMatrix(main0, off_d_val_u0, off_d_val_l0);
                A_implicit[i][0, 0] = 1.0;
                A_implicit[i][0, 1] = 0.0;
                A_implicit[i][NX - 1, NX - 1] = 1.0;
                A_implicit[i][NX - 1, NX - 2] = 0.0;
                #endregion
                //================================================
                // The B-row Matrix
                //================================================
                #region B_row
                D_plus = 0.5 * (D1 + D0);
                //main1[start_idx1] = 1.0 + 2.0 * D_plus * nu0;

                D_minus = 0.5 * (Dend1 + Dend);
                //main1[end_idx1] = 1.0 + 2.0 * D_minus * nu0;

                for (int j = 1; j < NX - 1; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
                    main1[j] = 1.0 + ((D_minus * nu0) + (D_plus * nu0));
                }
                for (int j = 1; j < NX; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
                    off_d_val_l1[j - 1] = -nu0 * D_minus;
                }
                for (int j = 0; j < NX-1; j++)
                {
                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
                    off_d_val_u1[j] = -nu0 * D_plus;
                }

                B_row[i] = new TridiagonalMatrix(main1, off_d_val_u1, off_d_val_l1);

                if (border_with_function[2].TypeBC == ABoundaryCondition.dirichlet)
                {
                    //D_plus = 0.5 * (D1 + D0);
                    B_row[i][0, 0] = 1.0; // D_plus * nu0;
                    B_row[i][0, 1] = 0.0;
                }
                else
                {
                    D_plus = 0.5 * (D1 + D0);
                    B_row[i][0, 0] = 1.0 + 2.0 * D_plus * nu0; // C0[0]
                    B_row[i][0, 1] = -2.0 * D_plus * nu0; // C0[0]
                    //B_row[i][1, 1] = 1.0 - dx * D_plus * nu0; // C0[0]

                }

                if (border_with_function[1].TypeBC == ABoundaryCondition.dirichlet)
                {
                    //
                    B_row[i][NX - 1, NX - 1] = 1.0; // D_minus * nu0;
                    B_row[i][NX - 1, NX - 2] = 0.0;
                }
                else
                {
                    D_minus = 0.5 * (Dend1 + Dend);
                    B_row[i][NX - 1, NX - 1] = 1.0 + 2.0 * D_minus * nu0; // dx; //1.0 - dx * D_plus * nu0; // C0[0]
                    B_row[i][NX - 1, NX - 2] = -2.0 * D_minus * nu0; //1.0; // dx * D_plus * nu0; // C0[0]
                    //B_row[i][NX - 2, NX - 2] = 1.0 - dx * D_plus * nu0; // C0[0]
                }
                #endregion
                ////================================================
                //// The B-column Matrix
                ////================================================
                #region B_col
                RVector D_col = D.GetColVector(i);
                D0 = D_col[0];
                D1 = D_col[1];
                Dend = D_col[nx - 1];
                Dend1 = D_col[nx - 2];

                D_plus = 0.5 * (D1 + D0);
                //main2[start_idx1] = 1.0 + 2.0 * D_plus * nu0;

                D_minus = 0.5 * (Dend1 + Dend);
                //main2[end_idx1] = 1.0 + 2.0 * D_minus * nu0;

                for (int j = 1; j < NX-1; j++)
                {
                    D_minus = 0.5 * (D_col[j - 1] + D_col[j]);
                    D_plus = 0.5 * (D_col[j + 1] + D_col[j]);
                    main2[j] = 1.0 + ((D_minus * nu0) + (D_plus * nu0));
                }
                for (int j = 1; j < NX; j++)
                {
                    D_minus = 0.5 * (D_col[j - 1] + D_col[j]);
                    off_d_val_l2[j - 1] = -nu0 * D_minus;
                }
                for (int j = 0; j < NX-1; j++)
                {
                    D_plus = 0.5 * (D_col[j + 1] + D_col[j]);
                    off_d_val_u2[j] = -nu0 * D_plus;
                }

                B_col[i] = new TridiagonalMatrix(main2, off_d_val_l2, off_d_val_u2);

                if (border_with_function[3].TypeBC == ABoundaryCondition.dirichlet)
                {
                    //D_plus = 0.5 * (D1 + D0);
                    B_col[i][0, 0] = 1.0; // D_plus * nu0;
                    B_col[i][0, 1] = 0.0;
                }
                else
                {
                    D_plus = 0.5 * (D1 + D0);
                    B_col[i][0, 0] = 1.0 + 2.0 * D_plus * nu0; // C0[0]
                    B_col[i][0, 1] = -2.0 * D_plus * nu0; // C0[0]

                }
                if (border_with_function[0].TypeBC == ABoundaryCondition.dirichlet)
                {
                    //
                    B_col[i][NX - 1, NX - 1] = 1.0; // D_minus * nu0;
                    B_col[i][NX - 1, NX - 2] = 0.0;
                }
                else
                {
                    D_minus = 0.5 * (Dend1 + Dend);
                    B_col[i][NX - 1, NX - 1] = 1.0 + 2.0 * D_minus * nu0; // dx; //1.0 - dx * D_plus * nu0; // C0[0]
                    B_col[i][NX - 1, NX - 2] = -2.0 * D_minus * nu0; //1.0; // dx * D_plus * nu0; // C0[0]
                    //B_row[i][NX - 2, NX - 2] = 1.0 - dx * D_plus * nu0; // C0[0]
                }
                #endregion
            }

            gxt_function = g;
            Chat_mode = Tmode;
            b_filename = base_filename;
            Errors = new string[nt];
        }

        // ====================================================================
        // Solvers
        // ====================================================================
        /// <summary>
        /// Method for solving the 2D diffusion equation using the 2D Alternating Direction Implicit algorithm.  Will output time-steps elapsed.
        /// </summary>
        /// <param name="n_time_steps"></param>
        /// <param name="output_interval"></param>
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
            RVector CT0 = new(ncols);
            RVector CT1 = new(ncols);
            RVector CB0 = new(ncols);
            RVector CB1 = new(ncols);

            double nu0 = dt / (2 * Math.Pow(dx, 2));
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
                            if (t == 0) { FileWriteData_CSV(full_file_name, X, Y, C_Initial); }
                            else { FileWriteData_CSV(full_file_name, X, Y, C_Im2); }
                        }
                    }

                    // 0 = top, 1 = right, 2 = left, 3 = bottom
                    CT = BCs_Functions[0].BoundaryFunction(t * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
                    CR = BCs_Functions[1].BoundaryFunction(t * dt, BCs_Functions[1].PositionVaries, BCs_Functions[1].PositionFixed);
                    CL = BCs_Functions[2].BoundaryFunction(t * dt, BCs_Functions[2].PositionVaries, BCs_Functions[2].PositionFixed);
                    CB = BCs_Functions[3].BoundaryFunction(t * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);

                    // ===================
                    // Explicit time-step
                    // ===================
                    for (int i = 0; i < NX; i++)
                    {
                        RVector xold;
                        if (t == 0) { xold = C_Initial.GetRowVector(i); } //xold = C_Initial.GetRowVector(i, start_idx2, nx_less1); 
                        else { xold = C_Im2.GetRowVector(i); } //, start_idx2, nx_less1
                        RVector v1 = A_explicit[i].Dot(xold);

                        // Left BC
                        if (BCs_Functions[2].TypeBC == ABoundaryCondition.dirichlet) { v1[0] = v1[0] + CL[i]; } //A_explicit[i][0, 1] *
                        // Right BCs
                        if (BCs_Functions[1].TypeBC == ABoundaryCondition.dirichlet) { v1[NX - 1] = v1[NX - 1] + CR[i]; } //A_explicit[i][0, 1] * 
                        C_Ex.ReplaceRow(v1, i); //, start_idx2, nx_less1
                    }

                    // Top BC
                    if (BCs_Functions[0].TypeBC == ABoundaryCondition.dirichlet) { C_Ex.ReplaceRow(CT, NY - 1); }
                    // Bottom BCs
                    if (BCs_Functions[3].TypeBC == ABoundaryCondition.dirichlet) { C_Ex.ReplaceRow(CB, 0); }

                    // ===================
                    // ===================
                    // One-half implicit time-step
                    // ===================
                    // Source terms
                    double t12 = (t + 0.5) * dt;
                    if (t == 0) { fn = gxt_function(X, Y, t12, cf_2D.DValues, C_Initial); }
                    else { fn = gxt_function(X, Y, t12, cf_2D.DValues, C_Im2); }

                    f12 = dt / 2.0 * fn;
                    RVector gn, g0;
                    for (int j = 0; j < NX; j++)
                    {
                        // BCs                        
                        gn = BCs_Functions[0].BoundaryFunction((t + 1) * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
                        g0 = BCs_Functions[0].BoundaryFunction(t * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
                        CT = (0.5 * B_row[j].Dot(gn)) + (0.5 * A_implicit[j].Dot(g0));

                        gn = BCs_Functions[3].BoundaryFunction((t + 1) * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
                        g0 = BCs_Functions[3].BoundaryFunction(t * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
                        CB = (0.5 * B_row[j].Dot(gn)) + (0.5 * A_implicit[j].Dot(g0));

                        RVector v1 = C_Ex.GetColVector(j);
                        RVector f12s = f12.GetColVector(j);
                        RVector v2, u12;
                        double s = -B_col[j][2, 1];
                        v2 = v1 + f12s;

                        // Bottom BC
                        if (BCs_Functions[3].TypeBC == ABoundaryCondition.dirichlet) { v2[1] = v2[1] + s * CB[j]; }
                        else if (BCs_Functions[3].TypeBC == ABoundaryCondition.neumann) { v2[0] = v2[0] + 2.0 * s * dy * CB[j]; }
                        // Top BCs
                        if (BCs_Functions[0].TypeBC == ABoundaryCondition.dirichlet) { v2[NX - 2] = v2[NX - 2] + s * CT[j]; }
                        else if (BCs_Functions[0].TypeBC == ABoundaryCondition.neumann) { v2[NX - 1] = v2[NX - 1] + 2.0 * s * dy * CT[j]; }

                        u12 = TridiagonalMatrix.Thomas_Algorithm(B_col[j], v2);
                        C_Im1.ReplaceCol(u12, j); //, ns, u12.GetRVectorSize
                    }

                    // Left BC
                    switch (BCs_Functions[2].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            C_Im1.ReplaceCol(CL, 0);
                            break;
                        case ABoundaryCondition.neumann:
                            //RVector c1 = C_Im1.GetColVector(1);
                            //RVector avalue = (dx * CL) + c1;
                            //C_Im1.ReplaceCol(avalue, 0);
                            break;
                    }
                    // Right BCs
                    switch (BCs_Functions[1].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            C_Im1.ReplaceCol(CR, NX - 1);
                            break;
                        case ABoundaryCondition.neumann:
                            //RVector c1 = 48 * C_Im1.GetColVector(NX - 2);
                            //RVector c2 = 36 * C_Im1.GetColVector(NX - 3);
                            //RVector c3 = 16 * C_Im1.GetColVector(NX - 4);
                            //RVector c4 = 3 * C_Im1.GetColVector(NX - 5);
                            //RVector avalue = 12 * dx * CR - c1 + c2 - c3 + c4;
                            //C_Im1.ReplaceCol(-1.0 / 25.0 * avalue, NX - 1);

                            //RVector c1 = C_Im1.GetColVector(NX - 2);
                            //RVector avalue = (dx * CR) + c1;
                            //C_Im1.ReplaceCol(avalue, NX - 1);
                            break;
                    }
                    // ===================
                    // ===================
                    // Full implicit time-step
                    // ===================
                    CR = BCs_Functions[1].BoundaryFunction((t + 1) * dt, BCs_Functions[1].PositionVaries, BCs_Functions[1].PositionFixed);
                    CL = BCs_Functions[2].BoundaryFunction((t + 1) * dt, BCs_Functions[2].PositionVaries, BCs_Functions[2].PositionFixed);
                    
                    //TridiagonalMatrix B_row2;
                    for (int k = 0; k < NY; k++)
                    {
                        RVector u1;
                        RVector v1 = C_Ex.GetRowVector(k);
                        RVector u12 = C_Im1.GetRowVector(k);
                        double s = -B_row[k][2, 1];
                        RVector b = (2 * u12) - v1;

                        // Left BC
                        if (BCs_Functions[2].TypeBC == ABoundaryCondition.dirichlet) { b[1] = b[1] + s * CL[k]; }
                        else if (BCs_Functions[2].TypeBC == ABoundaryCondition.neumann) { b[0] = b[0] + 2.0 * s * dx * CL[k]; }
                        // Right BC
                        if (BCs_Functions[1].TypeBC == ABoundaryCondition.dirichlet) { b[NX - 2] = b[NX - 2] + s * CR[k]; }
                        else if (BCs_Functions[1].TypeBC == ABoundaryCondition.neumann) { b[NX - 1] = b[NX - 1] + 2.0 * s * dx * CR[k]; }

                        u1 = TridiagonalMatrix.Thomas_Algorithm(B_row[k], b);
                        C_Im2.ReplaceRow(u1, k); //, ns, u12.GetRVectorSize + 1
                    }
                    // ===================
                    // Set top and bottom boundary conditions on the C_Im2 matrix
                    // ===================
                    switch (BCs_Functions[0].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            CT = BCs_Functions[0].BoundaryFunction(t * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
                            C_Im2.ReplaceRow(CT, nrows - 1);
                            break;
                        case ABoundaryCondition.neumann:
                            //CT = BCs_Functions[0].BoundaryFunction(t * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);

                            //RVector c1 = 48 * C_Im2.GetRowVector(NY - 2);
                            //RVector c2 = 36 * C_Im2.GetRowVector(NY - 3);
                            //RVector c3 = 16 * C_Im2.GetRowVector(NY - 4);
                            //RVector c4 = 3 * C_Im2.GetRowVector(NY - 5);
                            //RVector avalue = (12 * dy * CT) - c1 + c2 - c3 + c4;
                            //C_Im2.ReplaceRow((-1.0 / 25.0) * avalue, 0);

                            //RVector c1 = C_Im2.GetRowVector(NY - 2);
                            //RVector avalue = (2 * dy * CT) + c1;
                            //C_Im2.ReplaceRow(avalue, 0);

                            //CB = C_Im2.GetRowVector(nrows - 2);
                            //C_Im2.ReplaceRow(CT * 2 * dy + CB, nrows - 1);
                            break;
                    }
                    switch (BCs_Functions[3].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            CB = BCs_Functions[3].BoundaryFunction(t * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
                            C_Im2.ReplaceRow(CB, 0);
                            break;
                        case ABoundaryCondition.neumann:
                            //CT = C_Im2.GetRowVector(1);
                            //CB = BCs_Functions[3].BoundaryFunction(t * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);

                            //RVector c1 = 48 * C_Im2.GetRowVector(1);
                            //RVector c2 = 36 * C_Im2.GetRowVector(2);
                            //RVector c3 = 16 * C_Im2.GetRowVector(3);
                            //RVector c4 = 3 * C_Im2.GetRowVector(4);

                            //RVector c1 = C_Im2.GetRowVector(1);
                            //RVector avalue = (2 * dy * CB) + c1; // (12 * dy * CB) - c1 + c2 - c3 + c4;
                            //C_Im2.ReplaceRow(avalue, 0);

                            //RVector c1 = 2 * C_Im2.GetRowVector(1);
                            //RVector c2 = 0.5 * C_Im2.GetRowVector(2);
                            //RVector avalue = (dy * CB) - c1 + c2;
                            //C_Im2.ReplaceRow((-2.0 / 3.0) * avalue, 0);

                            //RVector c1 = C_Im2.GetRowVector(1) * 1.0e-4;
                            //C_Im2.ReplaceRow(c1, 0);

                            break;
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
            // Setup the final composition field
            if (error_flag == false) { C_Final = C_Im2; }
            else { C_Final = C_Initial; }
        }
    }

    // Hold this code.  This version of obtaining the derivative doesn't seem to work because there is a huge magnitude change over these points...
    //RVector c1 = 48 * C_Ex.GetRowVector(NY - 2);
    //RVector c2 = 36 * C_Ex.GetRowVector(NY - 3);
    //RVector c3 = 16 * C_Ex.GetRowVector(NY - 4);
    //RVector c4 = 3 * C_Ex.GetRowVector(NY - 5);
    //RVector avalue = 12 * dx * CT - c1 + c2 - c3 + c4;
    //C_Ex.ReplaceRow(-1.0 / 25.0 * avalue, NY - 1);
}
