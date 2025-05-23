package ru.mrsu.solvers.fftBSOR.impl;

import org.jtransforms.fft.DoubleFFT_1D;
import ru.mrsu.cuncurrency.Executor;
import ru.mrsu.model.DecartCoords;
import ru.mrsu.output.VtkWriter;
import ru.mrsu.solvers.fftBSOR.AbstractFftBsorCylinderPoissonSolver;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;

public class CylinderPoissonFftBsorSolverParallel extends AbstractFftBsorCylinderPoissonSolver {

    private final Executor executor = Executor.getExecutor();

    private double[][][] solveByFftBsor(int num, double r, double z0, double z1, double hr, double hfi, double hz, double eps) {
        int nr = (int) (r / hr)+1;
        int nfi = (int) (2 * Math.PI / hfi)+1;
        int nz = (int) ((z1 - z0) / hz)+1;

        double[][][] uin = new double[nr][nfi][nz];
        double[][][] f = new double[nr][nfi][nz];

        double[][][] uinRe = new double[nfi][nr][nz];
        double[][][] uinIm = new double[nfi][nr][nz];
        double[][][] fRe = new double[nfi][nr][nz];
        double[][][] fIm = new double[nfi][nr][nz];

        DoubleFFT_1D fftFi = new DoubleFFT_1D(nfi);

        // Нижняя и верхняя грань цилиндра
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nfi; j++) {
                // Нижняя
                uin[i][j][0] = u_r_fi_0(num, i * hr, j * hfi);
                // Верхняя
                uin[i][j][nz - 1] = u_r_fi_H(num, i * hr, j * hfi, z1);
            }
        }

        // Боковая грань
        for (int j = 0; j < nfi; j++) {
            for (int k = 1; k < nz - 1; k++) {
                uin[nr - 1][j][k] = u_R_fi_z(num, r, j * hfi,k * hz);
            }
        }

        // Функции f
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nfi; j++) {
                for (int k = 0; k < nz; k++) {
                    f[i][j][k] = f(num, i * hr, j * hfi, k * hz);
                }
            }
        }

        List<int[]> indicesR = splitArray(nr, executor.getThreadsCount(), true);

        // FFT
        long start1 = System.nanoTime();
        List<CompletableFuture<Void>> fftFutures = new ArrayList<>();
        for (int[] part: indicesR) {
            CompletableFuture<Void> fftFuture = executor.runAsync(() -> {
                for (int i = part[0]; i <= part[1]; i++) {
                    for (int k = 0; k < nz; k++) {
                        double[] phiSliceU = new double[nfi];
                        double[] phiSliceF = new double[nfi];
                        for (int j = 0; j < nfi; j++) {
                            phiSliceU[j] = uin[i][j][k];
                            phiSliceF[j] = f[i][j][k];
                        }

                        double[][] fftU = fft(fftFi, phiSliceU);
                        double[][] fftF = fft(fftFi, phiSliceF);
                        for (int m = 0; m < nfi; m++) {
                            uinRe[m][i][k] = fftU[0][m];
                            uinIm[m][i][k] = fftU[1][m];

                            fRe[m][i][k] = fftF[0][m];
                            fIm[m][i][k] = fftF[1][m];
                        }
                    }
                }
            });
            fftFutures.add(fftFuture);
        }
        for (CompletableFuture<Void> fftFuture : fftFutures) {
            fftFuture.join();
        }
        long end1 = System.nanoTime();
        printMillis(start1, end1, "FFT: ");


        // BSOR
        long start2 = System.nanoTime();

        double[] lambdas = getLambdas(nfi, hfi);
        double[] omegas = getOmegas(lambdas, nr, nz, hr, hz);

        List<CompletableFuture<Void>> resFutures = new ArrayList<>();
        for (int m = 0; m < nfi; m++) {
            int finalM = m;
            CompletableFuture<Void> resFuture = executor.runAsync(() -> {
                double[][][] resReIm = new double[2][nr][nz];
                blockSOR(lambdas[finalM], finalM, uinRe[finalM], uinIm[finalM], fRe[finalM], fIm[finalM], hr, hz, omegas[finalM], eps, resReIm);
                uinRe[finalM] = resReIm[0];
                uinIm[finalM] = resReIm[1];
            });
            resFutures.add(resFuture);
        }
        for (int m = 0; m < nfi; m++) {
            resFutures.get(m).join();
        }
        long end2 = System.nanoTime();
        printMillis(start2, end2, "all cores: ");


        // IFFT
        double[][][] uinRes = new double[nr][nfi][nz];
        long start3 = System.nanoTime();
        List<CompletableFuture<Void>> ifftFutures = new ArrayList<>();
        for (int[] part: indicesR) {
            CompletableFuture<Void> ifftFuture = executor.runAsync(() -> {
                for (int i = part[0]; i <= part[1]; i++) {
                    for (int k = 0; k < nz; k++) {
                        double[] phiSliceURe = new double[nfi];
                        double[] phiSliceUIm = new double[nfi];
                        for (int m = 0; m < nfi; m++) {
                            phiSliceURe[m] = uinRe[m][i][k];
                            phiSliceUIm[m] = uinIm[m][i][k];
                        }

                        double[] fftU = ifft(fftFi, phiSliceURe, phiSliceUIm);
                        for (int j = 0; j < nfi; j++) {
                            uinRes[i][j][k] = fftU[j];
                        }
                    }
                }
            });
            ifftFutures.add(ifftFuture);
        }
        for (CompletableFuture<Void> ifftFuture : ifftFutures) {
            ifftFuture.join();
        }
        long end3 = System.nanoTime();
        printMillis(start3, end3, "IFFT: ");

        printMillis(start1, end3, "total: ");

        DecartCoords[][][] decartCoords = DecartCoords.convertToDecartCoords3D(uinRes, hfi, hr, hz);
        VtkWriter.writeVtk(decartCoords, nr, nfi, nz, "fftBsorParallelRes");

        return uinRes;
    }

    public double[][][] solve(int num, double r, double z0, double z1, double hr, double hfi, double hz, double eps) {
        long start = System.nanoTime();
        double[][][] result = solveByFftBsor(num, r, z0, z1, hr, hfi, hz, eps);
        long end = System.nanoTime();
        writeCsvMetrics("fftBsorParallelMetrics", result, start, end, hr, hfi, hz);
        return result;
    }
}