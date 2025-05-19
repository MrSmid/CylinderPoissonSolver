package ru.mrsu.solvers.jacobi.impl;

import ru.mrsu.model.DecartCoords;
import ru.mrsu.output.VtkWriter;
import ru.mrsu.solvers.AbstractCylinderPoissonSolver;

public class CylinderPoissonJacobiIterable extends AbstractCylinderPoissonSolver {

    private double[][][] jacobi(int num, double r, double z0, double z1, double hr, double hfi, double hz, double eps) {
        long start = System.nanoTime();

        int nr = (int) (r / hr) + 1;
        int nfi = (int) (2 * Math.PI / hfi) + 1;
        int nz = (int) ((z1 - z0) / hz) + 1;
        double[][][] uin0 = new double[nr][nfi][nz];
        double[][][] uin1 = new double[nr][nfi][nz];
        double[][][] f = new double[nr][nfi][nz];

        // Нижняя и верхняя грань цилиндра
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nfi; j++) {
                // Нижняя
                uin0[i][j][0] = u_r_fi_0(num, i * hr, j * hfi);
                uin1[i][j][0] = uin0[i][j][0];
                // Верхняя
                uin0[i][j][nz - 1] = u_r_fi_H(num, i * hr, j * hfi, z1);
                uin1[i][j][nz - 1] = uin0[i][j][nz - 1];
            }
        }

        // Боковая грань
        for (int j = 0; j < nfi; j++) {
            for (int k = 1; k < nz - 1; k++) {
                uin0[nr - 1][j][k] = u_R_fi_z(num, r, j * hfi,k * hz);
                uin1[nr - 1][j][k] = uin0[nr - 1][j][k];
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

        boolean stop = false;
        int iterations = 0;
        while (!stop){
            iterations++;
            stop = true;
            for (int k = 1; k < nz - 1; k++) {
                for (int i = 1; i < nr - 1; i++) {
                    for (int j = 1; j < nfi - 1; j++) {
                        uin1[i][j][k] = tau[i] * ((rCoef1 * (uin0[i+1][j][k] + uin0[i-1][j][k])) + (rCoef2[i] * (uin0[i+1][j][k] - uin0[i-1][j][k])) + (fiCoef[i] * (uin0[i][j+1][k] + uin0[i][j-1][k])) + (zCoef * (uin0[i][j][k+1] + uin0[i][j][k-1])) - f[i][j][k]);
                        if (stop && Math.abs(uin0[i][j][k] - uin1[i][j][k]) > eps)
                            stop = false;
                    }
                }
            }
            // Заполнение границ матрицы значениями из соседних узлов при fi = 0 и fi = 2pi - hfi
            for (int i = 0; i < nr - 1; i++) {
                for (int k = 1; k < nz - 1; k++) {
                    uin1[i][0][k] = uin1[i][1][k];
                    uin1[i][nfi - 1][k] = uin1[i][nfi - 2][k];
                }
            }
            // Заполнение границ матрицы значениями из соседних узлов при r = 0
            for (int j = 0; j < nfi; j++) {
                for (int k = 1; k < nz - 1; k++) {
                    uin1[0][j][k] = uin1[1][j][k];
                }
            }
            for (int k = 1; k < nz - 1; k++) {
                for (int i = 0; i < nr - 1; i++) {
                    for (int j = 0; j < nfi; j++) {
                        uin0[i][j][k] = uin1[i][j][k];
                    }
                }
            }
        }

        long end = System.nanoTime();

        System.out.println("iterations: " + iterations);
        printMillis(start, end, "time: ");
        System.out.println("nodes: " + nr * nfi * nz);

        DecartCoords[][][] decartCoords = DecartCoords.convertToDecartCoords3D(uin0, hr, hfi, hz);
        VtkWriter.writeVtk(decartCoords, nr, nfi, nz, "jacobiIterableRes");

        return uin0;
    }


    public double[][][] solve(int num, double r, double z0, double z1, double hr, double hfi, double hz, double eps){
        long start = System.nanoTime();
        double[][][] result = jacobi(num, r, z0, z1, hr, hfi, hz, eps);
        long end = System.nanoTime();
        writeCsvMetrics("jacobiIterableMetrics", result, start, end, hr, hfi, hz);
        return result;
    }
}
