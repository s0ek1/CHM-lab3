package main;
import processor.*;
import java.util.Arrays;
public class Main {
    public static void main(String[] args) {
        Main main = new Main();
        main.run();
    }
    private void run() {
        LinearSystemSolver solver = new Solver();
        double[][] A = {{1.73, 0.27, 0.56},               //18 вариант
                        {-0.83, 0.53, -0.48},
                        {1.82, -0.64, 1.95}};
        double[] b = {0.36, 1.23, -0.76};
//        double[][] A = {{3.08, -1.5, 2.72},              //30 вариант
//                        {-1.5, 1.72, -1.58},
//                        {2.72, -1.58, 1.28}};
//        double[] b = {4.9, 2.14, 1.08};
        solver.print(A);
        double[][] Acopy = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++) Acopy[i] = Arrays.copyOf(A[i], A[i].length);
        solver.solveGauss(Acopy, b);
        solver.solveLU(Acopy, b);
        solver.solveSquareSqrt(A, b);   //22 24 25 28 29 30
    }
}

