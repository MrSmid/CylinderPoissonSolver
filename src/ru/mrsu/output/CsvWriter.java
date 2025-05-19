package ru.mrsu.output;

import ru.mrsu.output.support.CsvHeaders;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

public class CsvWriter {

    private final static String folderPath = "src/ru/mrsu/metrics";
    private static final List<CsvHeaders> headersOrdered = CsvHeaders.getHeadersOrdered();

    public static void appendToCsv(String filename, Map<CsvHeaders, String> data) {
        String fileNameFull = filename + ".csv";

        File folder = new File(folderPath);

        if (!folder.exists()) {
            boolean created = folder.mkdir();
            if (!created) {
                System.err.println("Не удалось создать папку " + folderPath);
                return;
            }
        }

        File file = new File(folder, fileNameFull);
        boolean isInitialWrite = !file.exists() || file.length() == 0;

        try (FileWriter writer = new FileWriter(file, true)) {
            if (isInitialWrite) {
                String headers = headersOrdered.stream()
                    .map(CsvHeaders::getValue)
                    .collect(Collectors.joining(","));
                writer.write(String.join(",", headers));
                writer.write("\n");
            }

            String joinedData = headersOrdered.stream()
                .map(it -> Optional.ofNullable(data.get(it)).orElse(""))
                .collect(Collectors.joining(","));

            writer.write(joinedData);
            writer.write("\n");
        } catch (IOException e) {
            System.err.println("Ошибка при записи в CSV: " + e.getMessage());
        }
    }
}
