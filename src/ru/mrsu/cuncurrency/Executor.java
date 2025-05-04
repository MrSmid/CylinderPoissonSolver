package ru.mrsu.cuncurrency;

import java.util.concurrent.CompletableFuture;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.function.Supplier;

public class Executor {

    private static Executor instance;
    private final ThreadPoolExecutor executor;
    private final int threadsCount;

    private Executor() {
        int threadCount = calculateThreadCount();
        this.threadsCount = threadCount;
        executor = new ThreadPoolExecutor(threadCount, threadCount,
            0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
    }

    public static Executor getExecutor() {
        if (instance == null) {
            instance = new Executor();
        }
        return instance;
    }

    public int getThreadsCount() {
        return threadsCount;
    }

    public <U> CompletableFuture<U> supplyAsync(Supplier<U> supplier) {
        return CompletableFuture.supplyAsync(supplier, executor);
    }

    public CompletableFuture<Void> runAsync(Runnable runnable) {
        return CompletableFuture.runAsync(runnable, executor);
    }

    public void shutdown() {
        executor.shutdown();
    }

    private int calculateThreadCount() {
        return Runtime.getRuntime().availableProcessors();
    }
}
