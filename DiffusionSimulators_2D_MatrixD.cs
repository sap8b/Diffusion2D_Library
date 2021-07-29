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
        private readonly TridiagonalMatrix[] A;
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

            A = new TridiagonalMatrix[nx];
            B_row = new TridiagonalMatrix[nx];
            B_col = new TridiagonalMatrix[nx];
            NeumannBCs_L_A = new double[nx];
            NeumannBCs_R_A = new double[nx];

            double nu0 = dt / (2 * Math.Pow(dx, 2));

            RVector nu = new(nx_less2);
            RVector off_d_val_l0 = new(nx_less3);
            RVector off_d_val_u0 = new(nx_less3);
            RVector main0 = new(nx_less2);
            RVector off_d_val_l1 = new(nx_less3);
            RVector off_d_val_u1 = new(nx_less3);
            RVector main1 = new(nx_less2);
            RVector off_d_val_l2 = new(nx_less3);
            RVector off_d_val_u2 = new(nx_less3);
            RVector main2 = new(nx_less2);

            double D_minus, D_plus;
            for (int i = start_idx1; i < NX; i++)
            {
                RVector D_row = D.GetRowVector(i);
                double D0 = D_row[start_idx1];
                double D1 = D_row[start_idx2];
                double Dend = D_row[nx_less1];
                double Dend1 = D_row[nx_less2];

                // Working this for diffusion heterogeneous media
                //================================================
                // A Matrix first
                //================================================
                #region A
                D_plus = 0.5 * (D1 + D0);
                main0[start_idx1] = 1.0 - 2.0 * D_plus * nu0;
                NeumannBCs_L_A[i] = 1.0 - D_plus * nu0;

                D_minus = 0.5 * (Dend1 + Dend);
                main0[end_idx1] = 1.0 - 2.0 * D_minus * nu0;
                NeumannBCs_R_A[i] = 1.0 - D_minus * nu0;

                for (int j = start_idx2; j < end_idx1; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
                    main0[j] = 1.0 - ((D_minus * nu0) + (D_plus * nu0));
                }

                for (int j = start_idx2; j < end_idx1 + 1; j++)
                {
                    D_minus = 0.5 * (D_row[j] + D_row[j - 1]);
                    off_d_val_l0[j - 1] = nu0 * D_minus;
                }

                for (int j = start_idx1; j < end_idx1; j++)
                {
                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
                    off_d_val_u0[j] = nu0 * D_plus;
                }

                A[i] = new TridiagonalMatrix(main0, off_d_val_l0, off_d_val_u0);
                #endregion
                ////================================================
                //// The B-row Matrix next
                ////================================================
                #region B_row
                D_plus = 0.5 * (D1 + D0);
                main1[start_idx1] = 1.0 + 2.0 * D_plus * nu0;

                D_minus = 0.5 * (Dend1 + Dend);
                main1[end_idx1] = 1.0 + 2.0 * D_minus * nu0;

                for (int j = start_idx2; j < end_idx1; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
                    main1[j] = 1.0 + ((D_minus * nu0) + (D_plus * nu0));
                }

                for (int j = start_idx2; j < end_idx1 + 1; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
                    off_d_val_l1[j - 1] = -nu0 * D_minus;
                }

                for (int j = start_idx1; j < end_idx1; j++)
                {
                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
                    off_d_val_u1[j] = -nu0 * D_plus;
                }

                B_row[i] = new TridiagonalMatrix(main1, off_d_val_l1, off_d_val_u1);
                #endregion
                ////================================================
                //// The B-column Matrix last
                ////================================================
                #region B_col
                RVector D_col = D.GetColVector(i);
                D0 = D_col[0];
                D1 = D_col[1];
                Dend = D_col[nx - 1];
                Dend1 = D_col[nx - 2];

                D_plus = 0.5 * (D1 + D0);
                main2[start_idx1] = 1.0 + 2.0 * D_plus * nu0;

                D_minus = 0.5 * (Dend1 + Dend);
                main2[end_idx1] = 1.0 + 2.0 * D_minus * nu0;

                for (int j = start_idx2; j < end_idx1; j++)
                {
                    D_minus = 0.5 * (D_col[j - 1] + D_col[j]);
                    D_plus = 0.5 * (D_col[j + 1] + D_col[j]);
                    main2[j] = 1.0 + ((D_minus * nu0) + (D_plus * nu0));
                }

                for (int j = start_idx2; j < end_idx1 + 1; j++)
                {
                    D_minus = 0.5 * (D_col[j - 1] + D_col[j]);
                    off_d_val_l2[j - 1] = -nu0 * D_minus;
                }

                for (int j = start_idx1; j < end_idx1; j++)
                {
                    D_plus = 0.5 * (D_col[j + 1] + D_col[j]);
                    off_d_val_u2[j] = -nu0 * D_plus;
                }

                B_col[i] = new TridiagonalMatrix(main2, off_d_val_l2, off_d_val_u2);
                #endregion
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

                            RVector Ctab = C_Initial.GetRowVector(ny - 2) + C0;
                            C_Initial.ReplaceRow(Ctab, ny - 1);
                        }
                        else if (border_with_function[i].TypeBC == ABoundaryCondition.neumann)
                        {
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[X.GetnRows - 1, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(ny - 2) + C0;
                            C_Initial.ReplaceRow(Ctab, ny - 1);
                        }
                        break;
                    case BoundingBox.right:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, Y.GetnCols - 1];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(nx - 2) + C0;
                            C_Initial.ReplaceCol(Ctab, nx - 1);
                        }
                        else if (border_with_function[i].TypeBC == ABoundaryCondition.neumann)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, Y.GetnCols - 1];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(nx - 2) + C0;
                            C_Initial.ReplaceCol(Ctab, nx - 1);
                        }
                        break;
                    case BoundingBox.left:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(1) + C0;
                            C_Initial.ReplaceCol(Ctab, 0);
                        }
                        else if (border_with_function[i].TypeBC == ABoundaryCondition.neumann)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(1) + C0;
                            C_Initial.ReplaceCol(C0, 0);
                        }
                        break;
                    case BoundingBox.bottom:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(1) + C0;
                            C_Initial.ReplaceRow(Ctab, 0);
                        }
                        else if (border_with_function[i].TypeBC == ABoundaryCondition.neumann)
                        {
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(1) + C0;
                            C_Initial.ReplaceRow(Ctab, 0);
                        }
                        break;
                }
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

                    // ===================
                    // Explicit time-step
                    // ===================
                    // 0 = top, 1 = right, 2 = left, 3 = bottom
                    CR = BCs_Functions[1].BoundaryFunction(t * dt, BCs_Functions[1].PositionVaries, BCs_Functions[1].PositionFixed);
                    CL = BCs_Functions[2].BoundaryFunction(t * dt, BCs_Functions[2].PositionVaries, BCs_Functions[2].PositionFixed);
                    double holdA0 = 0.0, holdA1 = 0.0;
                    for (int i = start_idx2; i < nx_less1; i++)
                    {
                        RVector xold;
                        if (t == 0) { xold = C_Initial.GetRowVector(i, start_idx2, nx_less1); }
                        else { xold = C_Im2.GetRowVector(i, start_idx2, nx_less1); }
                        if (BCs_Functions[2].TypeBC == ABoundaryCondition.neumann) { holdA0 = A[i][start_idx1, start_idx1]; A[i][start_idx1, start_idx1] = NeumannBCs_L_A[i]; }
                        if (BCs_Functions[1].TypeBC == ABoundaryCondition.neumann) { holdA1 = A[i][end_idx1, end_idx1]; A[i][end_idx1, end_idx1] = NeumannBCs_R_A[i]; }

                        RVector v1 = A[i].Dot(xold);
                        C_Ex.ReplaceRow(v1, i, start_idx2, nx_less1);

                        // Left BCs
                        switch (BCs_Functions[2].TypeBC)
                        {
                            case ABoundaryCondition.dirichlet:
                                C_Ex[i, start_idx1] = CL[i];
                                C_Ex[i, start_idx2] = C_Ex[i, start_idx2] + (-(NeumannBCs_L_A[i] - 1)) * CL[i];
                                break;
                            case ABoundaryCondition.neumann:
                                //C_Ex[i, start_idx2] = C_Ex[i, start_idx2] + (-(NeumannBCs_L_A[i] - 1)) * dx * CL[i];
                                //C_Ex[i, start_idx1] = C_Ex[i, start_idx2];
                                double c1 = 48 * C_Ex[i, 1];
                                double c2 = 36 * C_Ex[i, 2];
                                double c3 = 16 * C_Ex[i, 3];
                                double c4 = 3 * C_Ex[i, 4];
                                double avalue = 12 * dx * CL[i] - c1 + c2 + c3 - c4;
                                C_Ex[i, start_idx1] = -1 / 25 * avalue;
                                break;
                        }
                        // Right BCs
                        switch (BCs_Functions[1].TypeBC)
                        {
                            case ABoundaryCondition.dirichlet:
                                C_Ex[i, nx_less1] = CR[i];
                                C_Ex[i, nx_less2] = C_Ex[i, nx_less2] + (-(NeumannBCs_R_A[i] - 1)) * CR[i];
                                break;
                            case ABoundaryCondition.neumann:
                                //C_Ex[i, nx_less2] = C_Ex[i, nx_less2] + (-(NeumannBCs_R_A[i] - 1)) * dx * CR[i];
                                //C_Ex[i, nx_less1] = C_Ex[i, nx_less2];
                                double c1 = 48 * C_Ex[i, C_Ex.GetnCols - 1 - 1];
                                double c2 = 36 * C_Ex[i, C_Ex.GetnCols - 1 - 2];
                                double c3 = 16 * C_Ex[i, C_Ex.GetnCols - 1 - 3];
                                double c4 = 3 * C_Ex[i, C_Ex.GetnCols - 1 - 4];
                                double avalue = 12 * dx * CR[i] - c1 + c2 + c3 - c4;
                                C_Ex[i, C_Ex.GetnCols - 1] = -1 / 25 * avalue;
                                break;
                        }

                        if (BCs_Functions[2].TypeBC == ABoundaryCondition.neumann) { A[i][start_idx1, start_idx1] = holdA0; }
                        if (BCs_Functions[1].TypeBC == ABoundaryCondition.neumann) { A[i][end_idx1, end_idx1] = holdA1; }
                    }
                    // ===================
                    // ===================
                    // One-half implicit time-step
                    // ===================
                    // Source terms
                    double t12 = (t + 0.5) * dt;
                    if (t == 0) { fn = gxt_function(X, Y, t12, cf_2D.DValues, C_Initial); }
                    else { fn = gxt_function(X, Y, t12, cf_2D.DValues, C_Im2); }

                    f12 = dt / 2.0 * fn;

                    // BCs
                    RVector gn, g0;
                    gn = BCs_Functions[0].BoundaryFunction((t + 1) * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
                    g0 = BCs_Functions[0].BoundaryFunction(t * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
                    CT = (0.5 * B_row[nx_less1].Dot(gn.Section(start_idx2, nx_less1))) + (0.5 * A[nx_less1].Dot(g0.Section(start_idx2, nx_less1)));

                    gn = BCs_Functions[3].BoundaryFunction((t + 1) * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
                    g0 = BCs_Functions[3].BoundaryFunction(t * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
                    CB = (0.5 * B_row[start_idx1].Dot(gn.Section(start_idx2, nx_less1))) + (0.5 * A[start_idx1].Dot(g0.Section(start_idx2, nx_less1)));

                    TridiagonalMatrix B_col2;
                    for (int j = start_idx2; j < nx_less1; j++)
                    {
                        RVector v1 = C_Ex.GetColVector(j);
                        RVector f12s = f12.GetColVector(j);
                        RVector v2, u12;
                        double s = -B_col[j][1, 2];

                        // Bottom dirichlet, Top dirichlet
                        if (BCs_Functions[3].TypeBC == ABoundaryCondition.dirichlet && BCs_Functions[0].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            v2 = (v1 + f12s).Section(start_idx2, nx_less1);
                            v2[0] = s * (v1[0] + f12s[0]) + v1[1] + f12s[1];
                            v2[v2.GetRVectorSize - 1] = s * (v1[nx_less1] + f12s[nx_less1]) + v1[nx_less2] + f12s[nx_less2];
                            B_col2 = B_col[j];

                            u12 = TridiagonalMatrix.Thomas_Algorithm(B_col2, v2);
                            int ns = (nrows - 1) - u12.GetRVectorSize;
                            C_Im1.ReplaceCol(u12, j, ns, u12.GetRVectorSize + 1);

                            if (j < CB.GetRVectorSize) { C_Im1[0, j] = CB[j]; C_Im1[nrows - 1, j] = CT[j]; }
                            else { C_Im1[0, j] = CB[1]; C_Im1[nrows - 1, j] = CT[CT.GetRVectorSize - 1]; }
                        }
                        // Bottom neumann, Top dirichlet
                        else if (BCs_Functions[3].TypeBC == ABoundaryCondition.neumann && BCs_Functions[0].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            if (j < CB.GetRVectorSize)
                            {
                                v1[0] = (2 * dy * CB[j] / -3) - (4 / 3 * v1[1]) + (1 / 3 * v1[2]); //nu * 
                                v1[ncols - 1] = CT[j]; //nu * 
                            }
                            else
                            {
                                v1[0] = (2 * dy * CB[1] / -3) - (4 / 3 * v1[1]) + (1 / 3 * v1[2]); //nu * 
                                v1[ncols - 1] = CT[CT.GetRVectorSize - 1]; //nu * 
                            }
                            v2 = (v1 + f12s).Section(start_idx1, nx_less1);
                            v2[0] = 2 * s * dx * (v1[0] + f12s[0]) + v2[0]; //v1[0] + f12s[0]
                            v2[v2.GetRVectorSize - 1] = s * (v1[nx_less1] + f12s[nx_less1]) + v1[nx_less2] + f12s[nx_less2];

                            B_col2 = new(v2.GetRVectorSize, v2.GetRVectorSize);
                            B_col2.FitMainUpperLower(start_idx2, nx_less1, start_idx2, nx_less1 - 1, start_idx2, nx_less1 - 1, B_col[j].GetMain(), B_col[j].GetUpper(), B_col[j].GetLower());
                            B_col2[0, 1] = -2 * s;
                            B_col2[0, 0] = B_col2[1, 1];
                            B_col2[1, 0] = B_col2[1, 2];

                            u12 = TridiagonalMatrix.Thomas_Algorithm(B_col2, v2);

                            int ns = start_idx1;
                            C_Im1.ReplaceCol(u12, j, ns, u12.GetRVectorSize);
                            if (j < CB.GetRVectorSize) { C_Im1[nrows - 1, j] = CT[j]; }
                            else { C_Im1[nrows - 1, j] = CT[CT.GetRVectorSize - 1]; }
                        }
                        // Bottom dirichlet, Top neumann
                        else if (BCs_Functions[3].TypeBC == ABoundaryCondition.dirichlet && BCs_Functions[0].TypeBC == ABoundaryCondition.neumann)
                        {
                            v2 = (v1 + f12s).Section(start_idx2, NX);
                            v2[0] = s * (v1[0] + f12s[0]) + v1[1] + f12s[1];
                            v2[v2.GetRVectorSize - 1] = 2 * s * dx * (v1[nx_less1] + f12s[nx_less1]) + v1[nx_less2] + f12s[nx_less2];

                            B_col2 = new(v2.GetRVectorSize, v2.GetRVectorSize);
                            B_col2.FitMainUpperLower(start_idx1, nx_less2, start_idx1, nx_less2 - 1, start_idx1, nx_less2 - 1, B_col[j].GetMain(), B_col[j].GetUpper(), B_col[j].GetLower());
                            B_col2[nx_less2, nx_less3] = -2 * s;
                            B_col2[nx_less2, nx_less2] = B_col2[nx_less3, nx_less3];
                            B_col2[nx_less3, nx_less2] = B_col2[nx_less3, nx_less3 - 1];

                            u12 = TridiagonalMatrix.Thomas_Algorithm(B_col2, v2);

                            int ns = start_idx2;
                            C_Im1.ReplaceCol(u12, j, ns, u12.GetRVectorSize);
                            if (j < CB.GetRVectorSize) { C_Im1[0, j] = CB[j]; }
                            else { C_Im1[0, j] = CB[1]; }
                        }
                        // Bottom neumann, top neumann
                        else
                        {
                            v2 = v1 + f12s;
                            v2[0] = 2 * s * dx * (v1[0] + f12s[0]) + v1[1] + f12s[1];
                            v2[v2.GetRVectorSize - 1] = 2 * s * dx * (v1[nx_less1] + f12s[nx_less1]) + v1[nx_less2] + f12s[nx_less2];

                            B_col2 = new(v2.GetRVectorSize, v2.GetRVectorSize);
                            B_col2.FitMainUpperLower(start_idx2, nx_less1, start_idx2, nx_less1 - 1, start_idx2, nx_less1 - 1, B_col[j].GetMain(), B_col[j].GetUpper(), B_col[j].GetLower());

                            B_col2[0, 1] = -2 * s;
                            B_col2[0, 0] = B_col2[1, 1];
                            B_col2[1, 0] = B_col2[1, 2];

                            B_col2[nx_less1, nx_less2] = -2 * s;
                            B_col2[nx_less1, nx_less1] = B_col2[nx_less3, nx_less3];
                            B_col2[nx_less2, nx_less1] = B_col2[nx_less3, nx_less3 - 1];

                            u12 = TridiagonalMatrix.Thomas_Algorithm(B_col2, v2);

                            int ns = start_idx1;
                            C_Im1.ReplaceCol(u12, j, ns, u12.GetRVectorSize);
                        }
                    }
                    // ===================
                    // ===================
                    // Full implicit time-step
                    // ===================
                    CR = BCs_Functions[1].BoundaryFunction((t + 1) * dt, BCs_Functions[1].PositionVaries, BCs_Functions[1].PositionFixed);
                    CL = BCs_Functions[2].BoundaryFunction((t + 1) * dt, BCs_Functions[2].PositionVaries, BCs_Functions[2].PositionFixed);
                    
                    TridiagonalMatrix B_row2;
                    for (int k = 1; k < nrows - 1; k++)
                    {
                        RVector v1 = C_Ex.GetRowVector(k);
                        RVector u12 = C_Im1.GetRowVector(k);
                        RVector b = (2 * u12) - v1;
                        //if (k < CL.GetRVectorSize) { b[0] = CL[k]; b[nrows - 1] = CR[k]; }
                        //else { b[0] = CL[1]; b[nrows - 1] = CR[CR.GetRVectorSize - 1]; }

                        RVector b2, u1;
                        double s = -B_row[k][1, 2];
                        // Left dirichlet, Right dirichlet
                        if (BCs_Functions[2].TypeBC == ABoundaryCondition.dirichlet && BCs_Functions[1].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            b2 = b.Section(start_idx2, nx_less1);
                            //b2[0] = s * b[0];
                            //b2[b2.GetRVectorSize - 1] = s * b[nx_less1];
                            if (k < CL.GetRVectorSize)
                            {
                                b2[0] = s * b[0] + b[0];
                                b2[b2.GetRVectorSize - 1] = s * CR[k] + b[b.GetRVectorSize - 1];
                            }
                            else
                            {
                                b2[0] = s * b[1] + b[0];
                                b2[b2.GetRVectorSize - 1] = s * CR[CR.GetRVectorSize - 1] + b[b.GetRVectorSize - 1];
                            }

                            B_row2 = B_row[k];

                            u1 = TridiagonalMatrix.Thomas_Algorithm(B_row2, b2);
                            int ns = (nrows - 1) - u1.GetRVectorSize;
                            C_Im2.ReplaceRow(u1, k, ns, u12.GetRVectorSize + 1);
                            if (k < CL.GetRVectorSize) { C_Im2[0, k] = CL[k]; C_Im2[nrows - 1, k] = CR[k]; }
                            else { C_Im2[0, k] = CL[1]; C_Im2[nrows - 1, k] = CR[CR.GetRVectorSize - 1]; }
                        }
                        // Left neumann, Right dirichlet
                        else if (BCs_Functions[2].TypeBC == ABoundaryCondition.neumann && BCs_Functions[1].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            b2 = b.Section(start_idx1, nx_less1);
                            //b2[0] = 2 * s * dx * b[0];
                            //b2[b2.GetRVectorSize - 1] = s * b[nx_less1];
                            if (k < CL.GetRVectorSize)
                            {
                                b2[0] = 2 * s * dx * CL[k] + b[0];
                                b2[b2.GetRVectorSize - 1] = s * CR[k] + b[b.GetRVectorSize - 1];
                            }
                            else
                            {
                                b2[0] = 2 * s * dx * CL[1] + b[0];
                                b2[b2.GetRVectorSize - 1] = s * CR[CR.GetRVectorSize - 1] + b[b.GetRVectorSize - 1];
                            }

                            B_row2 = new(b2.GetRVectorSize, b2.GetRVectorSize);
                            B_row2.FitMainUpperLower(start_idx2, nx_less1, start_idx2, nx_less1 - 1, start_idx2, nx_less1 - 1, B_row[k].GetMain(), B_row[k].GetUpper(), B_row[k].GetLower());
                            B_row2[0, 1] = -2 * s;
                            B_row2[0, 0] = B_row2[1, 1];
                            B_row2[1, 0] = B_row2[1, 2];

                            u1 = TridiagonalMatrix.Thomas_Algorithm(B_row2, b2);

                            int ns = start_idx1;
                            C_Im2.ReplaceRow(u1, k, ns, u1.GetRVectorSize);
                            if (k < CL.GetRVectorSize) { C_Im2[nrows - 1, k] = CR[k]; }
                            else { C_Im2[nrows - 1, k] = CR[CR.GetRVectorSize - 1]; }
                        }
                        // Left dirichlet, Right neumann
                        else if (BCs_Functions[3].TypeBC == ABoundaryCondition.dirichlet && BCs_Functions[0].TypeBC == ABoundaryCondition.neumann)
                        {
                            b2 = b.Section(start_idx2, NX);
                            //b2[0] = s * b[0];
                            //b2[b2.GetRVectorSize - 1] = 2 * s * dx * b[nx_less1];
                            if (k < CL.GetRVectorSize)
                            {
                                b2[0] = s * b[0] + b[0];
                                b2[b2.GetRVectorSize - 1] = 2 * s * dx * CR[k] + b[b.GetRVectorSize - 1];
                            }
                            else
                            {
                                b2[0] = s * b[1] + b[0];
                                b2[b2.GetRVectorSize - 1] = 2 * s * dx * CR[CR.GetRVectorSize - 1] + b[b.GetRVectorSize - 1];
                            }

                            B_row2 = new(b2.GetRVectorSize, b2.GetRVectorSize);
                            B_row2.FitMainUpperLower(start_idx1, nx_less2, start_idx1, nx_less2 - 1, start_idx1, nx_less2 - 1, B_row[k].GetMain(), B_row[k].GetUpper(), B_row[k].GetLower());
                            B_row2[nx_less2, nx_less3] = -2 * s;
                            B_row2[nx_less2, nx_less2] = B_row2[nx_less3, nx_less3];
                            B_row2[nx_less3, nx_less2] = B_row2[nx_less3, nx_less3 - 1];

                            u1 = TridiagonalMatrix.Thomas_Algorithm(B_row2, b2);

                            int ns = start_idx2;
                            C_Im2.ReplaceRow(u1, k, ns, u1.GetRVectorSize);
                            if (k < CL.GetRVectorSize) { C_Im2[0, k] = CL[k]; }
                            else { C_Im2[0, k] = CL[1]; }
                        }
                        // Left neumann, Right neumann
                        else
                        {
                            b2 = new RVector(b.GetRVectorSize);
                            for (int kk = 0; kk < b.GetRVectorSize-1; kk++) { b2[kk] = b[kk]; }
                            
                            if (k < CL.GetRVectorSize)
                            {
                                b[0] = (2 * dx * CL[k] / -3) - (4 / 3 * b[1]) + (1 / 3 * b[2]);
                                b[b.GetRVectorSize - 1] = (2 * dx * CR[k] / -3) - (4 / 3 * b[b.GetRVectorSize - 2]) + (1 / 3 * b[b.GetRVectorSize - 3]);
                                b2[0] = s * b[0]; // 2 * s * dx * CL[k] + b[0]; // 
                                b2[b2.GetRVectorSize - 1] = s * b[b.GetRVectorSize - 1];  //2 * s * dx * CR[k] + b[b.GetRVectorSize - 1]; // 
                            }
                            else
                            {
                                b[0] = (2 * dx * CL[1] / -3) - (4 / 3 * b[1]) + (1 / 3 * b[2]);
                                b[b.GetRVectorSize - 1] = (2 * dx * CR[CR.GetRVectorSize - 1] / -3) - (4 / 3 * b[b.GetRVectorSize - 2]) + (1 / 3 * b[b.GetRVectorSize - 3]);
                                b2[0] = s * b[0]; // 2 * s * dx * CL[k] + b[0]; // 
                                b2[b2.GetRVectorSize - 1] = s * b[b.GetRVectorSize - 1];  //2 * s * dx * CR[k] + b[b.GetRVectorSize - 1]; //

                                //b2[0] = 2 * s * dx * CL[1] + b[1]; // 
                                //b2[b2.GetRVectorSize - 1] = 2 * s * dx * CR[CR.GetRVectorSize - 1] + b[b.GetRVectorSize - 1]; // 
                            }


                            B_row2 = new(b2.GetRVectorSize, b2.GetRVectorSize);
                            B_row2.FitMainUpperLower(start_idx2, nx_less1, start_idx2, nx_less1 - 1, start_idx2, nx_less1 - 1, B_row[k].GetMain(), B_row[k].GetUpper(), B_row[k].GetLower());

                            B_row2[0, 1] = -2 * s;
                            B_row2[0, 0] = B_row2[1, 1];
                            B_row2[1, 0] = B_row2[1, 2];

                            B_row2[nx_less1, nx_less2] = -2 * s;
                            B_row2[nx_less1, nx_less1] = B_row2[nx_less2, nx_less2];
                            B_row2[nx_less2, nx_less1] = B_row2[nx_less2, nx_less2 - 1];

                            u1 = TridiagonalMatrix.Thomas_Algorithm(B_row2, b2);

                            int ns = start_idx1;
                            C_Im2.ReplaceRow(u1, k, ns, u1.GetRVectorSize);

                            // Apply the BCs to the left and right sides
                            double c1 = 48 * C_Im2[k, 1];
                            double c2 = 36 * C_Im2[k, 2];
                            double c3 = 16 * C_Im2[k, 3];
                            double c4 = 3 * C_Im2[k, 4];
                            double avalue = 12 * dx * CL[k] - c1 + c2 + c3 - c4;
                            C_Im2[k, start_idx1] = -1 / 25 * avalue;

                            c1 = 48 * C_Im2[k, C_Im2.GetnCols - 1 - 1];
                            c2 = 36 * C_Im2[k, C_Im2.GetnCols - 1 - 2];
                            c3 = 16 * C_Im2[k, C_Im2.GetnCols - 1 - 3];
                            c4 = 3 * C_Im2[k, C_Im2.GetnCols - 1 - 4];
                            avalue = 12 * dx * CR[k] - c1 + c2 + c3 - c4;
                            C_Im2[k, C_Im2.GetnCols - 1] = -1 / 25 * avalue;
                        }
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
                            CT = BCs_Functions[0].BoundaryFunction(t * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
                            CB = C_Im2.GetRowVector(nrows - 2);
                            C_Im2.ReplaceRow(CT * 2 * dy + CB, nrows - 1);
                            break;
                    }
                    switch (BCs_Functions[3].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            CB = BCs_Functions[3].BoundaryFunction(t * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
                            C_Im2.ReplaceRow(CB, 0);
                            break;
                        case ABoundaryCondition.neumann:
                            CT = C_Im2.GetRowVector(1);
                            CB = BCs_Functions[3].BoundaryFunction(t * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);

                            C_Im2.ReplaceRow(CB * 2 * dy + CT, 0);
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
}
