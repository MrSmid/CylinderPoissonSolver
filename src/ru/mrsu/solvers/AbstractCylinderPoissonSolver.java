package ru.mrsu.solvers;

import ru.mrsu.model.DecartCoords;
import ru.mrsu.output.CsvWriter;
import ru.mrsu.output.VtkWriter;
import ru.mrsu.output.support.CsvHeaders;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

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
        VtkWriter.writeVtk(decartCoords, nr, nfi, nz, "resTochn");
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

    protected void writeCsvMetrics(String filename, double[][][] result, long start, long end, double hr, double hfi, double hz) {
        double diff = (end - start);
        double millis = diff / 1000000.;
        double seconds = millis / 1000.;
        double minutes = seconds / 60.;
        int nr = result.length;
        int nfi = result[0].length;
        int nz = result[0][0].length;
        int nodesTotal = nr * nfi * nz;

        Map<CsvHeaders, String> data = Map.of(
            CsvHeaders.COMMON_TIME_MS, String.valueOf(millis),
            CsvHeaders.COMMON_TIME_S, String.valueOf(seconds),
            CsvHeaders.COMMON_TIME_M, String.valueOf(minutes),
            CsvHeaders.H_R, String.valueOf(hr),
            CsvHeaders.H_FI, String.valueOf(hfi),
            CsvHeaders.H_Z, String.valueOf(hz),
            CsvHeaders.N_R, String.valueOf(nr),
            CsvHeaders.N_FI, String.valueOf(nfi),
            CsvHeaders.N_Z, String.valueOf(nz),
            CsvHeaders.NODES_TOTAL, String.valueOf(nodesTotal)
        );
        CsvWriter.appendToCsv(filename, data);
    }

    public abstract double[][][] solve(int num, double r, double z0, double z1, double hr, double hfi, double hz, double eps);
}