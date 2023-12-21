
using System.Diagnostics.SymbolStore;
using System.Runtime.Intrinsics.Arm;

double[,] datapoints = { { 3, 6 }, { 0, 3 }, { 2, 1 } };

double[,] datapoints1 = { { -1, 1.25 }, { 2, 3.5 } };
double[,] datapoints2 = { { 3, 6 }, { 3, 7 } };
double[,] datapoints3 = { { 3, 6 }, { 3, 6 }, { 2, 1 } };
Func<double, double> f = PolyInterpolate(datapoints);
Console.WriteLine(f(2));




Func<double, double> PolyInterpolate(double[,] DP)
{
    int n = DP.GetLength(0);
    int col = DP.GetLength(1);

    // Remove duplicates from the 2D array
    DP = RemoveDuplicates(DP);


    double[,] RemoveDuplicates(double[,] points)
    {
        int n = points.GetLength(0);
        double[,] result = new double[n,col];
        int count = 0;

        for (int i = 0; i < n; i++)
        {
            bool isDuplicate = false;
            for (int j = 0; j < i; j++)
            {
                // If a duplicate pair of points is found, set 'isDuplicate' to true.
                if (points[i, 0] == points[j, 0] && points[i, 1] == points[j, 1])
                {
                    isDuplicate = true;
                    break;
                }else if(points[i, 0] == points[j, 0])
                {
                    throw new Exception("Error: incorrect data points");
                }
            }

            // If no duplicate pair of points is found, add the unique pair of points to the 'result' array.
            if (!isDuplicate)
            {
                result[count, 0] = points[i, 0];
                result[count, 1] = points[i, 1];
                count++;
            }
        }
        return result;
    }

    double max = DP[0, 0];
    double min = DP[0, 0];

    // Finding the minimum and maximum x values in the given data points
    for (int i = 0; i < n; i++)
    {
        if (DP[i,0] > max) max = DP[i,0];
        if (DP[i,0] < min) min = DP[i,0];
    }

    // Generating the Vandermonde matrix
    double[,] vandermonde_matrix = VandermondeMatrix(DP);

    // Creating an augmented matrix by combining the Vandermonde matrix with the y-values from the data points
    double[,] augmentedMatrix = AugmentedMatrix(vandermonde_matrix, DP);

    // Solving the linear system to find the polynomial coefficients
    double[] PolyCoeffs = SystemSolve(augmentedMatrix);

    // Forming the Vandermonde matrix
    double[,] VandermondeMatrix(double[,] DP)
    {

        int degree = DP.GetLength(0) - 1;
        double[,] vMatrix = new double[n, n];

        // Filling in the Vandermonde matrix
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                vMatrix[i, j] = Math.Pow(DP[i, 0], degree - j);
            }
        }

        return vMatrix;
    }

    // Joining the Vandermonde matrix to the y-values from the DP 
    double[,] AugmentedMatrix(double[,] vm, double[,] DP)
    {
        int n = vm.GetLength(0);
        double[,] augmentedMatrix = new double[n, n + 1];

        // Copying the values from the Vandermonde matrix
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                augmentedMatrix[i, j] = vm[i, j];
            }
            // Appending the y-values from the data points to the augmented matrix
            augmentedMatrix[i, n] = DP[i, 1];
        }

        return augmentedMatrix;
    }

    // Solving the linear system to find the coefficients
    double[] SystemSolve(double[,] matrix)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);
        int n = Math.Min(rows, cols - 1); // the number of variables
        double[] solution = new double[n];

        // Perform the Gaussian elimination
        for (int k = 0; k < n; k++)
        {
            // Find a row that has a non-zero value in the k-th column
            int pivot = FindPivot(matrix, k);

            // If no such row exists, the system has no or infinitely many solutions
            if (pivot == -1)
            {

                return new double[0];

            }

            // Swap this row with the k-th row
            SwapRows(matrix, pivot, k);

            // Scale the k-th row, so the leading entry is equal to 1
            ScaleRow(matrix, k, 1 / matrix[k, k]);

            // Remove the k-th coefficient from the other rows
            for (int i = 0; i < rows; i++)
            {
                if (i != k)
                {
                    SubRow(matrix, i, k, matrix[i, k]);
                }
            }
        }

        // Check if the system is consistent
        for (int i = n; i < cols - 1; i++)
        {
            if (matrix[i, i] != 1)
            {
                // There is a row with non-zero constant term only, the system has no solution
                return new double[0];
            }
        }

        // function to find the first non-zero element in a column of a matrix
        int FindPivot(double[,] matrix, int col)
        {
            int rows = matrix.GetLength(0);
            for (int i = col; i < rows; i++)
            {
                if (matrix[i, col] != 0)
                {
                    return i;
                }
            }
            return -1;
        }

        // function to swap two rows of a matrix
        void SwapRows(double[,] matrix, int row1, int row2)
        {
            int cols = matrix.GetLength(1);
            for (int j = 0; j < cols; j++)
            {
                double temp = matrix[row1, j];
                matrix[row1, j] = matrix[row2, j];
                matrix[row2, j] = temp;
            }
        }

        // function to scale a row of a matrix by a factor
        void ScaleRow(double[,] matrix, int row, double factor)
        {
            int cols = matrix.GetLength(1);
            for (int j = 0; j < cols; j++)
            {
                matrix[row, j] *= factor;
            }
        }

        // function to subtract a multiple of one row to another row of a matrix
        void SubRow(double[,] matrix, int row1, int row2, double factor)
        {
            int cols = matrix.GetLength(1);
            for (int j = 0; j < cols; j++)
            {
                matrix[row1, j] -= matrix[row2, j] * factor;
            }
        }

        // Reading the solution from the last column

        for (int i = 0; i < n; i++)
        {
            solution[i] = matrix[i, cols - 1];
        }

        return solution;
    }

    // Return a lambda function that represents the interpolated polynomial
    return (x) => {
        double y = 0;
        /// create an if condition to check if the x augment falls within the range
        if(!(x <= max && x >= min))
        {
            throw new ArgumentException("Error: Out of bound argument");
        }

        // Evaluating the polynomial using the computed coefficients
        for (int i = 0; i < PolyCoeffs.Length; i++)
        {
            y = y * x + PolyCoeffs[i];            
        }
        return y;

    };
}


