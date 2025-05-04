package ru.mrsu.model;

public class DecartCoords {
    public final double[] coords = new double[3];
    public final Double value;

    DecartCoords(double[] roFiZ, double value) {
        coords[0] = roFiZ[0] * Math.cos(roFiZ[1]);
        coords[1] = roFiZ[0] * Math.sin(roFiZ[1]);
        coords[2] = roFiZ[2];
        this.value = value;
    }

    public static DecartCoords[][][] convertToDecartCoords3D(double[][][] uij, double hfi, double hr, double hz) {
        DecartCoords[][][] decartCoords = new DecartCoords[uij.length][uij[0].length][uij[0][0].length];
        for (int i = 0; i < uij.length; i++) {
            for (int j = 0; j < uij[i].length; j++) {
                for (int k = 0; k < uij[i][j].length; k++) {
                    double[] coords = new double[3];
                    coords[0] = hr * i;
                    coords[1] = hfi * j;
                    coords[2] = hz * k;
                    decartCoords[i][j][k] = new DecartCoords(coords, uij[i][j][k]);
                }
            }
        }
        return decartCoords;
    }
}