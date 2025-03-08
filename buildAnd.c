#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <slepceps.h>

PetscScalar uniformRnd(){
    return (double)rand() / RAND_MAX;
}

PetscScalar gaussRnd(PetscScalar mu, PetscScalar sig){

    double u1, u2;
    const long double pi = acosl(-1.0L);

    u1 = uniformRnd();
    u2 = uniformRnd();

    return (sig*cos(2*pi*u2)*sqrt(-2.*log(u1)) + mu);
}

// Helper function to convert (x, y, z) to a 1D index
int arr_index(int x, int y, int z, int Lx, int Ly) {
    return x + y * Lx + z * Lx * Ly;
}

// Function to compute Manhattan distance
int manhattan_distance(int x1, int y1, int z1, int x2, int y2, int z2, int L) {
    // return abs(x1 - x2) + abs(y1 - y2) + abs(z1 - z2);
    int dx = abs((x1 - x2 + L) % L);
    if (dx > L / 2) dx = L - dx;  // PBC adjustment for x
    int dy = abs((y1 - y2 + L) % L);
    if (dy > L / 2) dy = L - dy;  // PBC adjustment for y
    int dz = abs((z1 - z2 + L) % L);
    if (dz > L / 2) dz = L - dz;  // PBC adjustment for z
    return dx + dy + dz;
}

int main(int argc, char **argv){
    EPS eps;
    EPS eps_si;
    Vec xr;
    ST st;
    KSP ksp;
    PC pc;

    PetscInt L=2, k=1, N, nconv, itr=1, nev=1, cnt = 0;
    Mat A;
    PetscScalar t = 0, w = 1.0, w1 = 1.0, val, sig = 0, emax, emin, kr;
    PetscErrorCode ierr;
    char fvecname[100], fvalname[100];
    
    ierr = SlepcInitialize(&argc, &argv, NULL, NULL); if (ierr) return ierr;
    // Get the process rank and size
    PetscMPIInt petscRank, petscSize;
    MPI_Comm_rank(PETSC_COMM_WORLD, &petscRank);
    MPI_Comm_size(PETSC_COMM_WORLD, &petscSize);

    // dimension L
    PetscOptionsGetInt(NULL,NULL,"-l",&L,NULL);
    // # eigenpairs at energy e
    PetscOptionsGetInt(NULL,NULL,"-nev",&nev,NULL);
    PetscOptionsGetInt(NULL,NULL,"-itr",&itr,NULL);
    PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL);
    PetscOptionsGetReal(NULL,NULL,"-w",&w,NULL);
    PetscOptionsGetReal(NULL,NULL,"-e",&sig,NULL);
    PetscOptionsGetReal(NULL,NULL,"-t",&t,NULL);
    PetscOptionsGetReal(NULL,NULL,"-w1",&w1,NULL);

    N = L * L * L; // Total number of lattice sites
    
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N);
    MatSetFromOptions(A);
    MatSetUp(A);

    srand(time(0));

    // Populate Matrix
    for (int z1 = 0; z1 < L; z1++) {
        for (int y1 = 0; y1 < L; y1++) {
            for (int x1 = 0; x1 < L; x1++) {
                int i = arr_index(x1, y1, z1, L, L);

                MatSetValue(A, i, i, gaussRnd(0., w), INSERT_VALUES);
                // Loop through potential neighbors
                for (int z2 = 0; z2 < L; z2++) {
                    for (int y2 = 0; y2 < L; y2++) {
                        for (int x2 = 0; x2 < L; x2++) {
                            int j = arr_index(x2, y2, z2, L, L);

                            // Check Manhattan distance and ensure j > i
                            if (j > i && manhattan_distance(x1, y1, z1, x2, y2, z2, L) <= k) {
                                val = gaussRnd(t, w1);
                                MatSetValue(A, i, j, val, INSERT_VALUES);
                                MatSetValue(A, j, i, val, INSERT_VALUES);
				cnt++;
                            }
                        }
                    }
                }

            }
        }
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    // MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    // PetscPrintf(PETSC_COMM_WORLD, "# off-diagonal = %d\n", 2*cnt);

    MatCreateVecs(A,NULL,&xr);  

    EPSCreate(PETSC_COMM_WORLD, &eps);
    EPSSetOperators(eps, A, NULL);
    EPSSetProblemType(eps, EPS_HEP);
    EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
    EPSSetFromOptions(eps);
    EPSSolve(eps);
    EPSGetConverged(eps, &nconv);

    EPSGetEigenvalue(eps, 0, &kr, NULL);
    emin = (double)(kr);
    PetscPrintf(PETSC_COMM_WORLD, "min = %f\n", emin);

    EPSSetWhichEigenpairs(eps, EPS_LARGEST_REAL);
    EPSSetFromOptions(eps);
    EPSSolve(eps);
    EPSGetConverged(eps, &nconv);

    EPSGetEigenvalue(eps, 0, &kr, NULL);
    emax = (double)(kr);
    PetscPrintf(PETSC_COMM_WORLD, "max = %f\n", emax);

    // PetscPrintf(PETSC_COMM_WORLD, "sigma = %f\n", sig);
    EPSDestroy(&eps);

    // sig = sig*(emax-emin)+emin;

    EPSCreate(PETSC_COMM_WORLD, &eps_si);
    EPSSetOperators(eps_si, A, NULL);
    EPSSetProblemType(eps_si, EPS_HEP);
    EPSSetFromOptions(eps_si);
    EPSSetWhichEigenpairs(eps_si, EPS_TARGET_REAL);
    EPSSetTarget(eps_si, sig);

    EPSSetDimensions(eps_si, nev, PETSC_DEFAULT, PETSC_DEFAULT);
    EPSGetST(eps_si, &st);
    STSetType(st, STSINVERT);
    STGetKSP(st,&ksp);
    KSPSetType(ksp, KSPPREONLY);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCLU);
    // PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
    PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);

    EPSSolve(eps_si);
    EPSGetConverged(eps_si, &nconv);
    PetscPrintf(PETSC_COMM_WORLD, "#converged = %d\n", nconv);

    // snprintf(fvalname, 100, "./dataBinAE_L=%d/eigval_L=%d_w=%0.1f_t=%0.1f_k=%d_itr_%d_e=%0.1f.dat", L, L, w, t, k, itr, sig);
    // snprintf(fvecname, 100, "./dataBinAE_L=%d/eigvec_L=%d_w=%0.1f_t=%0.1f_k=%d_itr_%d_e=%0.1f.dat", L, L, w, t, k, itr, sig);

    snprintf(fvalname, 100, "./dataBinAE_PBC_L=%d/eigval_L=%d_w=%0.1f_w1=%0.1f_t=%0.1f_k=%d_itr_%d_e=%0.1f.dat", L, L, w, w1, t, k, itr, sig);
    snprintf(fvecname, 100, "./dataBinAE_PBC_L=%d/eigvec_L=%d_w=%0.1f_w1=%0.1f_t=%0.1f_k=%d_itr_%d_e=%0.1f.dat", L, L, w, w1, t, k, itr, sig);



    FILE* file1 = fopen(fvecname, "wb");
    FILE* file2 = fopen(fvalname, "wb");

    for (int i=0; i<nconv; i++){

        EPSGetEigenpair(eps_si,i,&kr,PETSC_NULLPTR,xr,PETSC_NULLPTR);
        // PetscPrintf(PETSC_COMM_WORLD, "e = %f\n", kr);

        Vec Vec_local;
        if (petscSize==1) { VecCreateSeq(PETSC_COMM_SELF,N,&Vec_local); VecCopy(xr,Vec_local); }
        else
            {
                VecScatter ctx; VecScatterCreateToZero(xr,&ctx,&Vec_local);
                VecScatterBegin(ctx,xr,Vec_local,INSERT_VALUES,SCATTER_FORWARD);
                VecScatterEnd(ctx,xr,Vec_local,INSERT_VALUES,SCATTER_FORWARD);
                VecScatterDestroy(&ctx);
            }

        if (petscRank==0)
            {
                fwrite(&kr, sizeof(double), 1, file2);
                PetscScalar * state; VecGetArray(Vec_local, &state);

                // fwrite(&state, sizeof(double), n, file);
                for (int i = 0; i < N; i++) {
                fwrite(&state[i], sizeof(double), 1, file1);
                }

            }
        VecDestroy(&Vec_local);
    }

    fclose(file1);
    fclose(file2);

    EPSDestroy(&eps_si);
    VecDestroy(&xr);

    PetscPrintf(MPI_COMM_WORLD, "itr = %d\n", itr);

    MatDestroy(&A);
    SlepcFinalize();


    return 0;
}
