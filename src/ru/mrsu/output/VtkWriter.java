package ru.mrsu.output;

import ru.mrsu.model.DecartCoords;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class VtkWriter {

    public static void writeVtk(DecartCoords[][][] u, int nr, int nfi, int nz, String fileName) {
        int roundCriteria = 1000000;
        int nfiNew = nfi+1;
        DecartCoords[][][] uNew = new DecartCoords[nr][nfiNew][nz];
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nfiNew; j++) {
                for (int k = 0; k < nz; k++) {
                    if (j == nfi) {
                        uNew[i][j][k] = u[i][0][k];
                    } else {
                        uNew[i][j][k] = u[i][j][k];
                    }
                }
            }
        }

        String fileNameFull = fileName + ".vtk";
        String folderPath = "src/ru/mrsu/results";

        File folder = new File(folderPath);

        if (!folder.exists()) {
            boolean created = folder.mkdir();
            if (!created) {
                System.err.println("Не удалось создать папку " + folderPath);
                return;
            }
        }

        File file = new File(folder, fileNameFull);

        try (FileWriter writer = new FileWriter(file)) {
            // Записываем заголовок файла
            writer.write("# vtk DataFile Version 3.0\n");
            writer.write("3D Structured Grid\n");
            writer.write("ASCII\n");
            writer.write("DATASET STRUCTURED_GRID\n");

            // Записываем размеры сетки
            writer.write(String.format("DIMENSIONS %d %d %d\n", nr, nfiNew, nz));

            // Записываем количество точек
            int numPoints = nr * nfiNew * nz;
            writer.write(String.format("POINTS %d float\n", numPoints));

            // Записываем координаты точек
            for (int k = 0; k < nz; k++) {
                for (int j = 0; j < nfiNew; j++) {
                    for (int i = 0; i < nr; i++) {
                        double[] coords = uNew[i][j][k].coords;
                        double x = round(coords[0], roundCriteria);
                        double y = round(coords[1], roundCriteria);
                        double z = round(coords[2], roundCriteria);
                        writer.write(x + " " + y + " " + z + "\n");
                    }
                }
            }

            // Записываем данные (значения функции в узлах)
            writer.write("POINT_DATA " + numPoints + "\n");
            writer.write("SCALARS value float 1\n");
            writer.write("LOOKUP_TABLE default\n");

            for (int k = 0; k < nz; k++) {
                for (int j = 0; j < nfiNew; j++) {
                    for (int i = 0; i < nr; i++) {
                        double value = round(uNew[i][j][k].value, roundCriteria);
                        writer.write(value + "\n");
                    }
                }
            }

            System.out.println("VTK файл успешно создан: " + fileName + ".vtk");

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static double round(double value, int roundCriteria) {
        double multValue = value * roundCriteria;
        long roundedValue = Math.round(multValue);
        return roundedValue / (double) roundCriteria;
    }
}
