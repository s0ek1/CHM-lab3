package processor;
import java.util.Arrays;
public class Solver implements LinearSystemSolver {
    @Override
    public void solveGauss(double[][] A, double[] b) {
        System.out.println("Решение методом Гаусса:\n" + "----------------------------------");
        int n = b.length;
        double[] x = new double[n];
        for (int k = 0; k < n - 1; k++) {
            for (int i = k + 1; i < n; i++) {
                double coef = A[i][k] / A[k][k];
                b[i] -= coef * b[k];
                for (int j = k; j < n; j++) {
                    A[i][j] -= coef * A[k][j];
                }
            }
        }
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        for (int i = 0; i < n; i++) System.out.println("x[" + i + "] = " + x[i]);
        System.out.println("Погрешность: ");
        estimateError(A, b, x);
    }
    @Override
    public void solveLU(double[][] A, double[] b) {
        System.out.println("Решение методом LU-разложения:\n" + "----------------------------------");
        int n = b.length;
        double[][] L = new double[n][n];
        double[][] U = new double[n][n];
        for (int i = 0; i < n; i++) {
            L[i][i] = 1.0;
        }
        for (int k = 0; k < n; k++) {
            U[k][k] = A[k][k];
            for (int i = k + 1; i < n; i++) {
                L[i][k] = A[i][k] / U[k][k];
                U[k][i] = A[k][i];
            }
            for (int i = k + 1; i < n; i++) {
                for (int j = k + 1; j < n; j++) {
                    A[i][j] = A[i][j] - L[i][k] * U[k][j];
                }
            }
        }
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            y[i] = b[i];
            for (int j = 0; j < i; j++) {
                y[i] = y[i] - L[i][j] * y[j];
            }
        }
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            x[i] = y[i];
            for (int j = i + 1; j < n; j++) {
                x[i] = x[i] - U[i][j] * x[j];
            }
            x[i] = x[i] / U[i][i];
        }
        for (int i = 0; i < n; i++) System.out.println("x[" + i + "] = " + x[i]);
        System.out.println("Погрешность: ");
        estimateError(A, b, x);
    }
    @Override
    public void solveSquareSqrt(double[][] A, double[] b) {
        System.out.println("Решение методом квадратного корня:\n" + "----------------------------------");
        loop:
        while (true) {
            for (int i = 0; i < A.length; i++) {
                for (int j = 0; j < A[i].length; j++) {
                    if (A[i][j] <= 0) { System.out.println("Матрица содержит отрицательные числа.");  break loop;}
                }
            }
        }
        double[][] At = new double[A.length][A.length];
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[i].length; j++) {
                At[j][i] = A[i][j];
            }
        }
        if (!Arrays.deepEquals(A, At))System.out.println("Матрица не является обратно симметричной.");
        int n = A.length;
        for (int i = 0; i < A.length; i++) {
            for (int j = i+1; j < A.length; j++) {
                A[j][i] = A[i][j];
            }
        }
        double[][] L = new double[A.length][A.length];
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j <= i; j++) {
                double sum = 0.0;
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * L[j][k];
                }
                if (i == j) {
                    L[i][i] = Math.sqrt(A[i][i] - sum);
                } else {
                    L[i][j] = (A[i][j] - sum) / L[j][j];
                }
            }
        }
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum += L[i][j] * y[j];
            }
            y[i] = (b[i] - sum) / L[i][i];
        }
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i+1; j < n; j++) {
                sum += L[j][i] * x[j];
            }
            x[i] = (y[i] - sum) / L[i][i];
        }
        for (int i = 0; i < n; i++) System.out.println("x[" + i + "] = " + x[i]);
        System.out.println("Погрешность: ");
        estimateError(A, b, x);
    }
    private void estimateError(double[][] A, double[] b, double[] x) {
        int n = b.length;
        double[] residual = new double[n];
        double[] Ax = new double[n];
        for (int i = 0; i < n; i++) {
            Ax[i] = 0.0;
            for (int j = 0; j < n; j++) {
                Ax[i] += A[i][j] * x[j];
            }
            residual[i] = Math.abs(b[i] - (Ax[i]+1e-15));
        }
        for (int i = 0; i < n; i++) {
            System.out.println(" для x[" + i + "]: " + residual[i]);
        }
        System.out.println();
    }
    @Override
    public void print(double[][] A) {
        for (double[] doubles : A) {for (int j = 0; j < A.length; j++) {
            System.out.print(String.format("%.2f", doubles[j])+"\t");
        } System.out.println(); }
        System.out.println();
    }
}