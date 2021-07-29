using System;
using System.IO;

using static Diffusion2D_Library.BoundaryCondition;

namespace Diffusion2D_Library
{
    // ====================================================================
    // Method delegates used by the 2D diffusion equation solvers
    // ====================================================================

    /// <summary>
    /// Delegate for methods for boundary conditions that specify both constant and variable locations.  Used if there is a position-dependence to the boundary condition
    /// </summary>
    /// <param name="time">time, in seconds</param>
    /// <returns>A vector storing the composition on a boundary of a 2D field</returns>
    //public delegate RVector BoundaryCondition_Del(double time, int n);
    public delegate RVector Del_BC_xy(double t, RVector SlidingSide, double FixedValue);
    /// <summary>
    /// Delegate for assigning the initial composition field in a 2D region
    /// </summary>
    /// <param name="c"></param>
    /// <returns>A matrix of composition values</returns>
    public delegate RMatrix Del_IC_xy(RMatrix x, RMatrix y);
    /// <summary>
    /// Delegate for calculating sources/sinks for reactions occurring in a 2D composition field
    /// </summary>
    /// <param name="position"></param>
    /// <param name="time"></param>
    /// <param name="composition"></param>
    /// <returns>A matrix of composition values</returns>
    public delegate RMatrix Del_Source_xy(RMatrix xposition, RMatrix yposition, double time, double DiffCoeff, RMatrix composition);
    /// <summary>
    /// Delegate for calculating sources/sinks for reactions occurring in a 2D composition field with diffusivity as a matrix
    /// </summary>
    /// <param name="xposition"></param>
    /// <param name="yposition"></param>
    /// <param name="time"></param>
    /// <param name="DiffCoeff"></param>
    /// <param name="composition"></param>
    /// <returns></returns>
    public delegate RMatrix Del_Source_MatrixD(RMatrix xposition, RMatrix yposition, double time, RMatrix DiffCoeff, RMatrix composition);
    public delegate RVector Del_BC_Constant(int n, double v);
    public delegate RMatrix Del_IC_Constant(int nr, int nc, double v);
    public delegate RMatrix Del_Source_Constant(int nr, int nc, double v);
    // ====================================================================

    public class DiffusionSimulator2D
    {
        /// <summary>
        /// Defines a structure to enclose the matrices that define a 2D composition
        /// </summary>
        public struct CompositionField2D
        {
            public RMatrix InitialCompositionValues;
            public RMatrix FinalCompositionValues;
            public RMatrix XPositionValues;
            public RMatrix YPositionValues;
            public RMatrix DValues;

            /// <summary>
            /// Conctructor for CompositionField2D
            /// </summary>
            /// <param name="nx">number of x positions</param>
            /// <param name="ny">number of y positions</param>
            public CompositionField2D(int nx, int ny)
            {
                InitialCompositionValues = new(nx, ny);
                FinalCompositionValues = new(nx, ny);
                XPositionValues = new(nx, ny);
                YPositionValues = new(nx, ny);
                DValues = new(nx, ny);
            }
        }

        /// <summary>
        /// Enumeration that specifies the type of boundary condition to be applied to a boundary
        /// </summary>
        public enum BoundingBox { left, right, top, bottom }
        /// <summary>
        /// Enumeration that specifies whether text is output to the output stream
        /// </summary>
        public enum Mode { quiet, verbose }
        /// <summary>
        /// 
        /// </summary>
        public struct BoundaryWithFunction
        {
            public BoundingBox BoundaryLocation;
            public ABoundaryCondition TypeBC;
            public Del_BC_xy BoundaryFunction;
            public RVector PositionVaries;
            public double PositionFixed;
        }
        public struct BoundaryWithVector
        {
            public BoundingBox BoundaryLocation;
            public ABoundaryCondition TypeBC;
            public RVector FunctionValues;
        }

        // Constants
        public const string suffix = ".csv";

        // Fields
        internal string b_filename;
        internal string[] errors;
        internal double dt;
        internal double dx;
        internal double dy;
        internal double d_coeff;
        internal CompositionField2D cf_2D;
        internal BoundaryWithFunction[] border_with_function;
        internal BoundaryWithVector[] border_with_vector;
        internal Mode text_mode;
        internal bool error_flag = false;
        internal int nx, ny, nx_less1, nx_less2, nx_less3, end_idx1, end_idx2, start_idx1, start_idx2;
        public Del_IC_xy I0;
        public Del_Source_MatrixD gxt_function;
        public Del_Source_xy gxt_xy;
        public RMatrix gxt_matrix;

        // Properties
        public int NX
        {
            get => nx;
            set => nx = value;
        }
        public int NY
        {
            get => ny;
            set => ny = value;
        }
        public RMatrix C_Initial
        {            
            get => cf_2D.InitialCompositionValues;
            set
            {
                if (value.GetnCols > 0 && value.GetnRows > 0) { cf_2D.InitialCompositionValues = value; }
            }
        }
        public RMatrix C_Final
        {
            get => cf_2D.FinalCompositionValues;
            set { if (value.GetnCols > 0 && value.GetnRows > 0) { cf_2D.FinalCompositionValues = value; } }
        }
        public RMatrix X => cf_2D.XPositionValues;
        public RMatrix Y => cf_2D.YPositionValues;
        public BoundaryWithFunction[] BCs_Functions
        {
            get => border_with_function;
            set => border_with_function = value;
        }
        public BoundaryWithVector[] BCs_Vector
        {
            get => border_with_vector;
            set => border_with_vector = value;
        }
        public Mode Chat_mode
        {
            get => text_mode;
            set => text_mode = value;
        }
        public string Base_filename
        {
            get => b_filename;
            set => b_filename = value;
        }
        public string[] Errors
        {
            get => errors;
            set => errors = value;
        }


        // ====================================================================
        // Helpful methods
        // ====================================================================
        /// <summary>
        /// Method for outputting composition and position data to a csv file
        /// </summary>
        /// <param name="of">Filename</param>
        /// <param name="x">Position</param>
        /// <param name="c">Composition</param>
        public static void FileWriteData_CSV(string of, RMatrix x, RMatrix y, RMatrix c)
        {
            string check_dir = Path.GetDirectoryName(of);
            string owd;

            if (check_dir == null)
            {
                string cwd = Directory.GetCurrentDirectory();
                owd = Path.Combine(cwd, of);
            }
            else
            {
                owd = of;
            }
            if (File.Exists(owd)) { File.Delete(owd); }

            FileStream fS = new(owd, FileMode.OpenOrCreate);
            StreamWriter sW = new(fS);

            string header = "x,y,c";
            sW.WriteLine(header);

            int nrows = x.GetnRows;
            int ncols = x.GetnCols;

            if (nrows > 0 && ncols > 0)
            {
                for (int i = 0; i < nrows; i++)
                {
                    for (int j = 0; j < ncols; j++)
                    {
                        string line = x[i, j].ToString() + "," + y[i, j].ToString() + "," + c[i, j].ToString();
                        sW.WriteLine(line);
                    }
                }

            }
            else
            {
                throw new Exception("No data available to write to the file!");
            }

            sW.Close();
            fS.Close();

        }

        // =======================================================================
        // Private methods
        // =======================================================================
        /// <summary>
        /// Converts the enumerated boundary condition to a string
        /// </summary>
        /// <param name="bc">variable of type BoundaryCondition</param>
        /// <returns></returns>
        private static string ConvertEnumBCToString(ABoundaryCondition bc)
        {
            string result = bc switch
            {
                ABoundaryCondition.dirichlet => "Dirichlet",
                ABoundaryCondition.neumann => "Neumann",
                _ => "Dirichlet",
            };
            return result;
        }
        /// <summary>
        /// Method to convert a string value to one of the enumerated boundary conditions
        /// </summary>
        /// <param name="value">string variable to be converted to the boundary condition enum</param>
        /// <returns></returns>
        private static ABoundaryCondition ConvertStringToEnumBC(string value)
        {
            ABoundaryCondition bc = value switch
            {
                "Dirichlet" => ABoundaryCondition.dirichlet,
                "Neumann" => ABoundaryCondition.neumann,
                _ => ABoundaryCondition.dirichlet,
            };
            return bc;
        }
        private static RVector ShapeMatrixToVector(RMatrix rm)
        {
            int nrows = rm.GetnRows;
            int ncols = rm.GetnCols;

            int nvals = nrows * ncols;

            RVector rv = new(nvals);

            int counter = 0;
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    rv[counter] = rm[i, j];
                    counter++;
                }
            }
            return rv;
        }
        private static RMatrix ShapeVectorToMatrix(RVector rv, int nrows, int ncols)
        {
            RMatrix rm = new(nrows, ncols);

            int counter = 0;
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    rm[i, j] = rv[counter];
                    counter++;
                }
            }
            return rm;
        }
        // =======================================================================
        //
        // =======================================================================
        // =======================================================================
        // =======================================================================
    }
}
