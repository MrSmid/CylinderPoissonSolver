package ru.mrsu.solvers.fftBSOR;

import org.jtransforms.fft.DoubleFFT_1D;
import ru.mrsu.solvers.AbstractCylinderPoissonSolver;

public abstract class AbstractFftBsorCylinderPoissonSolver extends AbstractCylinderPoissonSolver {

    protected void blockSOR(double lambda, int m, double[][] uRe, double[][] uIm, double[][] fRe, double[][] fIm,
                            double hr, double hz, double omega, double eps, double[][][] resReIm) {
        int nr = uRe.length;
        int nz = uRe[0].length;

        double[][] uSorRe = new double[nr][nz];
        double[][] uSorIm = new double[nr][nz];
        for (int i = 0; i < nr; i++) {
            System.arraycopy(uRe[i], 0, uSorRe[i], 0, nz);
            System.arraycopy(uIm[i], 0, uSorIm[i], 0, nz);
        }

        double invHr2 = 1.0 / (hr * hr);
        double invHz2 = 1.0 / (hz * hz);
        double a = invHz2;
        double c = a;
        double twoInvHr2 = 2.0 * invHr2;
        double twoInvHz2 = 2.0 * invHz2;

        int iters = 0;
        boolean stop = false;
        while (!stop) {
            stop = true;
            // Фиксация по i и решение системы прогонкой
            for (int i = 1; i < nr-1; i++) {
                double ri = i * hr;
                double inv2RiHr = 1.0 / (2.0 * ri * hr);

                double b = -1 * (twoInvHr2 + (lambda / (ri * ri)) + twoInvHz2);

                double[] dRe = new double[nz-2];
                double[] dIm = new double[nz-2];
                for (int k = 1; k < nz-1; k++) {
                    double diffURe = uSorRe[i+1][k] - uSorRe[i-1][k];
                    double diffUIm = uSorIm[i+1][k] - uSorIm[i-1][k];
                    double summURe = uSorRe[i+1][k] + uSorRe[i-1][k];
                    double summUIm = uSorIm[i+1][k] + uSorIm[i-1][k];
                    dRe[k-1] = fRe[i][k] - (summURe * invHr2) - (diffURe * inv2RiHr);
                    dIm[k-1] = fIm[i][k] - (summUIm * invHr2) - (diffUIm * inv2RiHr);
                }
                dRe[0] -= a * uSorRe[i][0];
                dRe[nz-3] -= a * uSorRe[i][nz-1];
                dIm[0] -= a * uSorIm[i][0];
                dIm[nz-3] -= a * uSorIm[i][nz-1];

                // Прогонка
                double[] resRe = tridiagonalMatrixSolve(a, b, c, dRe);
                double[] resIm = tridiagonalMatrixSolve(a, b, c, dIm);

                // Релаксация + проверка условия выхода
                for (int k = 1; k < nz-1; k++) {
                    double uOldSorRe = uSorRe[i][k];
                    double uOldSorIm = uSorIm[i][k];
                    double uNewSorRe = omega * resRe[k-1] + (1 - omega) * uOldSorRe;
                    double uNewSorIm = omega * resIm[k-1] + (1 - omega) * uOldSorIm;

                    if (stop && (Math.abs(uNewSorRe - uOldSorRe) > eps || Math.abs(uNewSorIm - uOldSorIm) > eps)) {
                        stop = false;
                    }

                    uSorRe[i][k] = uNewSorRe;
                    uSorIm[i][k] = uNewSorIm;
                }
            }

            // Заполнение границ матрицы значениями из соседних узлов при r = 0 для первой моды
            if (m == 0) {
                for (int k = 1; k < nz - 1; k++) {
                    uSorRe[0][k] = uSorRe[1][k];
                    uSorIm[0][k] = uSorIm[1][k];
                }
            }

            iters++;
        }

        //System.out.println("mode: " + m + ", iters: " + iters);
        resReIm[0] = uSorRe;
        resReIm[1] = uSorIm;
    }

    protected double[] tridiagonalMatrixSolve(double a, double b, double c, double[] d) {
        int n = d.length;

        // Прямой ход прогонки
        double[] alpha = new double[n - 1];
        double[] beta = new double[n];

        alpha[0] = -c / b;
        beta[0] = d[0] / b;

        for (int i = 1; i < n - 1; i++) {
            double denominator = b + a * alpha[i - 1];
            alpha[i] = -c / denominator;
            beta[i] = (d[i] - a * beta[i - 1]) / denominator;
        }

        beta[n - 1] = (d[n - 1] - a * beta[n - 2]) /
            (b + a * alpha[n - 2]);

        // Обратный ход прогонки
        double[] x = new double[n];
        x[n - 1] = beta[n - 1];

        for (int i = n - 2; i >= 0; i--) {
            x[i] = alpha[i] * x[i + 1] + beta[i];
        }

        return x;
    }

    protected double[][] fft(DoubleFFT_1D fft, double[] vec) {
        double[] vecRI1D = new double[vec.length * 2];
        for (int i = 0; i < vec.length; i++) {
            vecRI1D[2*i] = vec[i];
            vecRI1D[2*i+1] = 0.;
        }
        fft.complexForward(vecRI1D);
        double[][] dataRI = new double[2][vec.length];
        for (int i = 0; i < vec.length; i++) {
            dataRI[0][i] = vecRI1D[i*2];
            dataRI[1][i] = vecRI1D[i*2+1];
        }
        return dataRI;
    }

    protected double[] ifft(DoubleFFT_1D fft, double[] dataRe, double[] dataIm) {
        double[] vecRI1D = new double[dataRe.length * 2];
        for (int i = 0; i < dataRe.length; i++) {
            vecRI1D[2*i] = dataRe[i];
            vecRI1D[2*i+1] = dataIm[i];
        }
        fft.complexInverse(vecRI1D, true);
        double[] vecR = new double[dataRe.length];
        for (int i = 0; i < vecR.length; i++) {
            vecR[i] = vecRI1D[i*2];
        }
        return vecR;
    }

    protected double[] getOmegas(double[] lambdas, double nr, double nz, double hr, double hz) {
        double[] omegas = new double[lambdas.length];
        double l_0 = 4 * ((Math.pow(Math.sin(Math.PI / (2 * (nr + 1))) / (hr * hr), 2)) + ((Math.pow(Math.sin(Math.PI / (2 * (nz + 1))) / (hz * hz), 2))));

        omegas[0] = 1.95;

        for (int i = 1; i < lambdas.length; i++) {
            double l_max_i = l_0 + Math.abs(lambdas[i]) / hr;
            double l_min_i = l_0 + Math.abs(lambdas[i]) / (nr * hr);
            omegas[i] = 2 / (1 + Math.sqrt(1 - Math.pow(l_min_i / l_max_i, 2)));
        }

        return omegas;
    }

    protected double[] getLambdas(int nfi, double hfi) {
        double[] lambdas = new double[nfi];
        for (int m = 0; m < nfi; m++) {
            int mVal = m <= nfi / 2 ? m : m - nfi;
            lambdas[m] = (4 / (hfi*hfi)) * Math.pow(Math.sin(mVal * hfi / 2), 2);
        }
        return lambdas;
    }
}