package ru.mrsu.solvers.jacobi.impl;

import ru.mrsu.cuncurrency.Executor;
import ru.mrsu.model.DecartCoords;
import ru.mrsu.output.VtkWriter;
import ru.mrsu.solvers.AbstractCylinderPoissonSolver;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;

public class CylinderPoissonJacobiParallel extends AbstractCylinderPoissonSolver {

    Executor executor = Executor.getExecutor();

    private double[][][] jacobi(int num, double r, double z0, double z1, double hr, double hfi, double hz, double eps) {
        long start = System.nanoTime();

        int nr = (int) (r / hr)+1;
        int nfi = (int) (2 * Math.PI / hfi)+1;
        int nz = (int) ((z1 - z0) / hz)+1;
        double[][][][] uin = {new double[nr][nfi][nz], new double[nr][nfi][nz]};
        double[][][] f = new double[nr][nfi][nz];

        // Нижняя и верхняя грань цилиндра
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nfi; j++) {
                // Нижняя
                uin[0][i][j][0] = u_r_fi_0(num, i * hr, j * hfi);
                uin[1][i][j][0] = uin[0][i][j][0];
                // Верхняя
                uin[0][i][j][nz - 1] = u_r_fi_H(num, i * hr, j * hfi, z1);
                uin[1][i][j][nz - 1] = uin[0][i][j][nz - 1];
            }
        }

        // Боковая грань
        for (int j = 0; j < nfi; j++) {
            for (int k = 1; k < nz - 1; k++) {
                uin[0][nr - 1][j][k] = u_R_fi_z(num, r, j * hfi, k * hz);
                uin[1][nr - 1][j][k] = uin[0][nr - 1][j][k];
            }
        }

        // Функции f
        for (int i = 1; i < nr-1; i++) {
            for (int j = 1; j < nfi-1; j++) {
                for (int k = 1; k < nz-1; k++) {
                    f[i][j][k] = f(num, i * hr, j * hfi, k * hz);
                }
            }
        }

        double[] tau = new double[nr];
        double[] rCoef2 = new double[nr];
        double[] fiCoef = new double[nr];
        for (int i = 0; i < nr; i++) {
            fiCoef[i] = 1. / (Math.pow((hr * i) * hfi, 2));
            tau[i] = 1. / (2. * ((1. / Math.pow(hr, 2)) + fiCoef[i] + (1. / Math.pow(hz, 2))));
            rCoef2[i] = 1. / (2 * (hr * i) * hr);
        }
        double rCoef1 = 1. / (Math.pow(hr, 2));
        double zCoef = 1. / (Math.pow(hz, 2));

        List<int[]> borderIndicesNotInclusive = splitArray(nfi, executor.getThreadsCount(), false);

        boolean stop = false;
        int iterations = 0;
        while (!stop){
            iterations++;
            stop = true;

            List<CompletableFuture<Boolean>> futures = new ArrayList<>();
            for (int threadIndex = 0; threadIndex < borderIndicesNotInclusive.size(); threadIndex++) {
                int finalThreadIndex = threadIndex;
                CompletableFuture<Boolean> future = executor.supplyAsync(() ->
                    calculatingCycle(uin[0], uin[1], borderIndicesNotInclusive.get(finalThreadIndex), nz, nr, nfi, f, tau, rCoef1, rCoef2, fiCoef, zCoef, eps)
                );
                futures.add(future);
            }
            for (CompletableFuture<Boolean> futureStop: futures) {
                if (!futureStop.join()) {
                    stop = false;
                }
            }

            // Заполнение границ матрицы значениями из соседних узлов при fi = 0 и fi = 2pi - hfi
            CompletableFuture<Void> fillPhiStartAndEndFuture = executor.runAsync(() ->
                fillPhiStartAndEnd(uin[1], nr, nfi, nz)
            );
            // Заполнение границ матрицы значениями из соседних узлов при r = 0
            CompletableFuture<Void> fillRStartFuture = executor.runAsync(() ->
                fillRStart(uin[1], nfi, nz)
            );

            fillPhiStartAndEndFuture.join();
            fillRStartFuture.join();

            double[][][] tmpMatrix = uin[0];
            uin[0] = uin[1];
            uin[1] = tmpMatrix;
        }

        long end = System.nanoTime();

        System.out.println("iterations: " + iterations);
        printMillis(start, end, "time: ");
        System.out.println("nodes: " + nr * nfi * nz);

        DecartCoords[][][] decartCoords3D = DecartCoords.convertToDecartCoords3D(uin[0], hfi, hr, hz);
        VtkWriter.writeVtk(decartCoords3D, nr, nfi, nz, "jacobiParallelRes");

        return uin[0];
    }

    private boolean calculatingCycle(double[][][] uin0, double[][][] uin1, int[] borderIndices, int nz, int nr, int nfi, double[][][] f, double[] tau, double rCoef1, double[] rCoef2, double[] fiCoef, double zCoef, double eps) {
        boolean stop = true;
        for (int k = 1; k < nz - 1; k++) {
            for (int i = 1; i < nr - 1; i++) {
                for (int j = borderIndices[0]; j <= borderIndices[1]; j++) {
                    uin1[i][j][k] = tau[i] * ((rCoef1 * (uin0[i+1][j][k] + uin0[i-1][j][k])) + (rCoef2[i] * (uin0[i+1][j][k] - uin0[i-1][j][k])) + (fiCoef[i] * (uin0[i][j+1][k] + uin0[i][j-1][k])) + (zCoef * (uin0[i][j][k+1] + uin0[i][j][k-1])) - f[i][j][k]);
                    if (stop && Math.abs(uin0[i][j][k] - uin1[i][j][k]) > eps)
                        stop = false;
                }
            }
        }
        return stop;
    }

    private void fillPhiStartAndEnd(double[][][] uin1, int nr, int nfi, int nz) {
        // Заполнение границ матрицы значениями из соседних узлов при fi = 0 и fi = 2pi - hfi
        for (int i = 0; i < nr - 1; i++) {
            for (int k = 1; k < nz - 1; k++) {
                uin1[i][0][k] = uin1[i][nfi - 2][k];
                uin1[i][nfi - 1][k] = uin1[i][1][k];
            }
        }
    }

    private void fillRStart(double[][][] uin1, int nfi, int nz) {
        // Заполнение границ матрицы значениями из соседних узлов при r = 0
        for (int j = 0; j < nfi; j++) {
            int oppositeFiIndex = (j + nfi / 2) % nfi;
            for (int k = 1; k < nz - 1; k++) {
                uin1[0][j][k] = uin1[1][oppositeFiIndex][k];
            }
        }
    }

    public double[][][] solve(int num, double r, double z0, double z1, double hr, double hfi, double hz, double eps){
        return jacobi(num, r, z0, z1, hr, hfi, hz, eps);
    }
}


