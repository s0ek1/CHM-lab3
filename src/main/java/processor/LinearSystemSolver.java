package processor;
public interface LinearSystemSolver {
    void solveGauss(double[][] A, double[] b);
    void solveLU(double[][] A, double[] b);
    void solveSquareSqrt(double[][] A, double[] b);
    void print(double[][] A);

}