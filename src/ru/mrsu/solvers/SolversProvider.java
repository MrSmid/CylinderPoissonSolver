package ru.mrsu.solvers;

import ru.mrsu.solvers.fftBSOR.impl.CylinderPoissonFftBsorSolverIterable;
import ru.mrsu.solvers.fftBSOR.impl.CylinderPoissonFftBsorSolverParallel;
import ru.mrsu.solvers.jacobi.impl.CylinderPoissonJacobiIterable;
import ru.mrsu.solvers.jacobi.impl.CylinderPoissonJacobiParallel;

public enum SolversProvider {
    JACOBI_ITERABLE(new CylinderPoissonJacobiIterable()),
    JACOBI_PARALLEL(new CylinderPoissonJacobiParallel()),
    FFT_BSOR_ITERABLE(new CylinderPoissonFftBsorSolverIterable()),
    FFT_BSOR_PARALLEL(new CylinderPoissonFftBsorSolverParallel());

    private final AbstractCylinderPoissonSolver solver;

    SolversProvider(AbstractCylinderPoissonSolver solver) {
        this.solver = solver;
    }

    public double[][][] solve(int num, double r, double z0, double z1, double hr, double hfi, double hz, double eps) {
        return solver.solve(num, r, z0, z1, hr, hfi, hz, eps);
    }

    public AbstractCylinderPoissonSolver getInstance() {
        return solver;
    }
}
