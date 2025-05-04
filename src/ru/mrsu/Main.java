package ru.mrsu;

import ru.mrsu.solvers.SolversProvider;

public class Main {

    public static void main(String[] args) {
        // Номер уравнения и граничных условий
        int equationNum = 1;
        // Радиус
        double r = 1;
        // Нижняя координата z
        double z0 = 0;
        // Верхняя координата z
        double z1 = 1;
        // Шаг сетки
        double h = 0.03;
        // Точность
        double eps = Math.pow(h, 3)/10;

        SolversProvider.JACOBI_ITERABLE.solve(equationNum, r, z0, z1, h, h, h, eps);
        SolversProvider.JACOBI_PARALLEL.solve(equationNum, r, z0, z1, h, h, h, eps);
        SolversProvider.FFT_BSOR_ITERABLE.solve(equationNum, r, z0, z1, h, h, h, eps);
        SolversProvider.FFT_BSOR_PARALLEL.solve(equationNum, r, z0, z1, h, h, h, eps);
    }
}