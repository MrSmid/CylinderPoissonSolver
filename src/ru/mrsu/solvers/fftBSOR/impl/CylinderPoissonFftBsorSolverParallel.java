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

    private class ModeInfo {
        final int m;
        final double lambda;
        final double weight;

        double omega;

        ModeInfo(int m, double lambda) {
            this.m = m;
            this.lambda = lambda;
            if (Math.abs(lambda) < 1E-20) {
                this.weight = 100.;
            } else {
                this.weight = (1 / lambda) + (1. / m);
            }
        }
    }

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
        List<ModeInfo> modeInfoList = new ArrayList<>();
        for (int m = 0; m < nfi; m++) {
            int mVal = m <= nfi / 2 ? m : m - nfi;
            double lambda = (4 / (hfi*hfi)) * Math.pow(Math.sin(mVal * hfi / 2), 2);
            modeInfoList.add(new ModeInfo(m, lambda));
        }
        List<List<ModeInfo>> coreTasks = getCoreTasksLoad(modeInfoList);
        List<CompletableFuture<Void>> coreTasksFutures = new ArrayList<>();
        for (int i = 0; i < coreTasks.size(); i++) {
            int finalI = i;
            CompletableFuture<Void> resFuture = executor.runAsync(() -> {
                long coreStart = System.nanoTime();
                for (ModeInfo modeInfo : coreTasks.get(finalI)) {
                    double[][][] resReIm = new double[2][nr][nz];
                    blockSOR(modeInfo.lambda, modeInfo.m, uinRe[modeInfo.m], uinIm[modeInfo.m], fRe[modeInfo.m], fIm[modeInfo.m], hr, hz, modeInfo.omega, eps, resReIm);
                    uinRe[modeInfo.m] = resReIm[0];
                    uinIm[modeInfo.m] = resReIm[1];
                }
                long coreEnd = System.nanoTime();
                printMillis(coreStart, coreEnd, "core " + finalI + ", tasks: " + coreTasks.get(finalI).size() + ", time" + ": ");
            });
            coreTasksFutures.add(resFuture);
        }
        for (CompletableFuture<Void> coreTasksFuture : coreTasksFutures) {
            coreTasksFuture.join();
        }
        long end2 = System.nanoTime();
        printMillis(start2, end2, "BSOR: ");


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

    private List<List<ModeInfo>> getCoreTasksLoad(List<ModeInfo> modeInfoList) {
        modeInfoList.sort((a, b) -> Double.compare(b.weight, a.weight));
        double totalSumWeight = 0.0;
        for (ModeInfo mi : modeInfoList) {
            totalSumWeight += mi.weight;
        }
        List<List<ModeInfo>> coreTasksLoad = new ArrayList<>();
        int index = 0;
        while (index < modeInfoList.size()) {
            double avgLoad = totalSumWeight / (executor.getThreadsCount() - coreTasksLoad.size());
            List<ModeInfo> oneCoreTask = new ArrayList<>();
            double oneCoreTaskWeight = 0;
            while (oneCoreTaskWeight < avgLoad) {
                ModeInfo modeInfo = modeInfoList.get(index);
                oneCoreTask.add(modeInfo);
                index++;
                if (index == modeInfoList.size())
                    break;
                oneCoreTaskWeight += modeInfo.weight;
            }
            coreTasksLoad.add(oneCoreTask);
            for (ModeInfo modeTask : oneCoreTask) {
                totalSumWeight -= modeTask.weight;
                modeTask.omega = 1.95 - 0.95 * (coreTasksLoad.size() - 1) / executor.getThreadsCount();
            }
        }
        return coreTasksLoad;
    }

    public double[][][] solve(int num, double r, double z0, double z1, double hr, double hfi, double hz, double eps) {
        return solveByFftBsor(num, r, z0,z1, hr, hfi, hz, eps);
    }
}


