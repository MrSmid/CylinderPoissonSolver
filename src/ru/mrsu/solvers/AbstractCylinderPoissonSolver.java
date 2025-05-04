package ru.mrsu.solvers;

import ru.mrsu.model.DecartCoords;
import ru.mrsu.output.VtkWriter;

import java.util.ArrayList;
import java.util.List;

public abstract class AbstractCylinderPoissonSolver {

    // Нижняя грань цилиндра
    protected double u_r_fi_0(int num, double r, double fi){
        return switch (num) {
            case 1 -> r * (1-r) * Math.cos(fi);
            default -> 0.;
        };
    }

    // Верхняя грань цилиндра
    protected double u_r_fi_H(int num, double r, double fi, double H){
        return switch (num) {
            case 1 -> 1 + r * (1-r) * Math.cos(fi);
            default -> 0.;
        };
    }

    // Боковая грань
    protected double u_R_fi_z(int num, double R, double fi, double z){
        return switch (num) {
            case 1 -> z;
            default -> 0.;
        };
    }

    // Правая часть уравнения Пуассона
    protected double f(int num, double r, double fi, double z){
        return switch (num) {
            case 1 -> -3 * Math.cos(fi);
            default -> 0;
        };
    }

    // Точное решение
    protected double u_tochn(int num, double r, double fi, double z) {
        return switch (num) {
            case 1 -> r * (1-r) * Math.cos(fi) + z;
            default -> 0;
        };
    }

    protected double[][][] getUTochn(int num, double r, double z0, double z1, double hr, double hfi, double hz) {
        int nr = (int) (r / hr)+1;
        int nfi = (int) (2 * Math.PI / hfi) + 1;
        int nz = (int) ((z1 - z0) / hz)+1;

        double[][][] uin = new double[nr][nfi][nz];
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nfi; j++) {
                for (int k = 0; k < nz; k++) {
                    uin[i][j][k] = u_tochn(num, i*hr, j*hfi, k*hz);
                }
            }
        }

        DecartCoords[][][] decartCoords = DecartCoords.convertToDecartCoords3D(uin, hfi, hr, hz);
        VtkWriter.writeVtk(decartCoords, nr, nfi, nz, "fftResTochn");
        return uin;
    }

    protected List<int[]> splitArray(int arrayLength, int n, boolean borderInclusive) {
        List<int[]> segments = new ArrayList<>();
        int partSize = arrayLength / n;
        int remainder = arrayLength % n;

        int startIndex = 0;

        for (int i = 0; i < n; i++) {
            int endIndex = startIndex + partSize + (i < remainder ? 1 : 0) - 1;

            if (i == 0 && !borderInclusive) {
                segments.add(new int[]{startIndex + 1, endIndex});
            } else if (i == n - 1 && !borderInclusive) {
                segments.add(new int[]{startIndex, endIndex - 1});
            } else {
                segments.add(new int[]{startIndex, endIndex});
            }
            startIndex = endIndex + 1;
        }

        return segments;
    }

    protected void printMillis(long start, long end) {
        double diff = (end - start) / 1000000.;
        System.out.println(diff);
    }

    protected void printMillis(long start, long end, String prefix) {
        double diff = (end - start) / 1000000.;
        System.out.println(prefix + diff);
    }

    public abstract double[][][] solve(int num, double r, double z0, double z1, double hr, double hfi, double hz, double eps);
}