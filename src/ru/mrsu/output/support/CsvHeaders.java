package ru.mrsu.output.support;

import java.util.List;

public enum CsvHeaders {
    COMMON_TIME_MS("time_ms"),
    COMMON_TIME_S("time_s"),
    COMMON_TIME_M("time_m"),
    H_R("hr"),
    H_FI("hfi"),
    H_Z("hz"),
    N_R("nr"),
    N_FI("nfi"),
    N_Z("nz"),
    NODES_TOTAL("nodes_total");

    private final String value;

    CsvHeaders(String value) {
        this.value = value;
    }

    public String getValue() {
        return value;
    }

    public static List<CsvHeaders> getHeadersOrdered() {
        return List.of(
            COMMON_TIME_MS,
            COMMON_TIME_S,
            COMMON_TIME_M,
            H_R,
            H_FI,
            H_Z,
            N_R,
            N_FI,
            N_Z,
            NODES_TOTAL
        );
    }
}
