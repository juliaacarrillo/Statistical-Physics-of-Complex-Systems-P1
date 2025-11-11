
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define N 1000               // Number of hard disks
#define sigma 1.0            // Particle diameter (units of sigma) 
#define delta 0.1
double phi;                  // Packing fraction
#define mass 1.0
#define kT 1.0
#define N_MC_STEPS 10000     // Number of Monte Carlo steps
#define N_SIM 5
#define random01() ((double)rand()/(double)RAND_MAX) // Random number 0..1

static double Lx,Ly,L,g_field; // Box dimensions and gravity

// Function declarations
void initial_conditions(double *x, double *y);
void initial_conditions_closed_box(double *x, double *y);
double dist_pbc(double dx);
double distance_pbc(double x1, double y1, double x2, double y2);
int trial_move(double *x, double *y, double *xu, double *yu, int i);
void save_positions_to_file(double *x, double *y, const char *filename);
void triangular_lattice(double *x,double *y, int Npart,double L);
int trial_move_gravity(double *x, double *y, double *xu, double *yu, int i);
int trial_move_gravity_closed_box(double *x, double *y, int i);
void compute_box_size(void);
void triangular_conditions_safe(double *x, double *y, double L);
double compute_energy_per_particle(double *y);
double wrap(double coord, double L);

int main (){

    // ================== CODE BLOCK 1 ==================
    /*
    srand((unsigned)time(NULL));

    double pi = 3.14159265358979323846;
    double area_particle = pi * ((sigma/2.0)*(sigma/2.0));

    L = sqrt((double)N * area_particle / phi);
    printf("Packing fraction phi = %lf, computed box size L = %lf\n", phi, L);

    static double x[N], y[N], xo[N], yo[N], xu[N], yu[N]; // positions
    double msd, dx, dy;

    double deltas[] = {0.001, 0.003, 0.01, 0.03, 0.1, 0.3};
    int n_deltas = sizeof(deltas) / sizeof(deltas[0]);

    int n, i, d, n_mc;

    // Initialize all arrays to zero
    for(i = 0; i < N; i++){
        x[i] = y[i] = xo[i] = yo[i] = xu[i] = yu[i] = 0.0;
    }

    for (d = 0; d < n_deltas; d++){
        delta = deltas[d];
        char filename[64];
        sprintf(filename, "msd_delta_%0.3f.txt", delta);
        FILE *g = fopen(filename, "w");
        if(!g){ perror("error opening"); exit(1); }

        int n_samples = N_MC_STEPS / 10;
        double msd_avg[n_samples];
        for (int k = 0; k < n_samples; k++) msd_avg[k] = 0.0;

        printf("Running delta = %g\n", delta);

        for (int sim = 0; sim < N_SIM; sim++){
            // Set initial positions
            initial_conditions(x, y);
            for (i = 0; i < N; i++){
                xu[i] = x[i];
                yu[i] = y[i];
                xo[i] = xu[i];
                yo[i] = yu[i];
            }

            // Equilibrate system
            int N_EQ = 200;
            for (int k = 0; k < N_EQ; k++){
                for (n = 0; n < N; n++){
                    int idx = rand() % N;
                    trial_move(x, y, xu, yu, idx);
                }
            }

            // Reset reference after equilibration
            for (i = 0; i < N; i++){ xo[i] = xu[i]; yo[i] = yu[i]; }

            int step_index = 0;

            // Monte Carlo steps
            for (n_mc = 0; n_mc < N_MC_STEPS; n_mc++){
                for (n = 0; n < N; n++){
                    int idx = rand() % N;
                    trial_move(x, y, xu, yu, idx);
                }

                // Compute MSD every 10 steps
                if (n_mc % 10 == 0){
                    msd = 0.0;
                    for (i = 0; i < N; i++){
                        dx = xu[i] - xo[i];
                        dy = yu[i] - yo[i];
                        msd += dx*dx + dy*dy;
                    }
                    msd /= (double)N;
                    if (step_index < n_samples) msd_avg[step_index] += msd;
                    step_index++;
                }
            } // end MC steps
        } // end N_SIM

        // Write average MSD to file
        for (int k = 0; k < n_samples; k++){
            double avg = msd_avg[k] / (double)N_SIM;
            fprintf(g, "%lf\t%d\t%le\n", delta, k*10, avg);
        }
        fclose(g);
    } // end deltas

    printf("Simulation finished.\n");
    */

    // ================== CODE BLOCK 2 ==================
    /*
    srand(time(NULL));
    double phis[] = {0.05, 0.2, 0.5};
    int n_phis = 3;
    double pi = 3.14159265358979323846;
    double area_particle = pi * (sigma/2.0) * (sigma/2.0);

    static double x[N], y[N], xo[N], yo[N], xu[N], yu[N];

    for (int p = 0; p < n_phis; p++) {
        phi = phis[p];
        L = sqrt((double)N * area_particle / phi);
        printf("\n=== Running for phi = %.2f, box L = %.3f ===\n", phi, L);

        // Place particles in triangular lattice
        triangular_conditions(x, y);
        for (int i = 0; i < N; i++) { xu[i] = x[i]; yu[i] = y[i]; xo[i] = x[i]; yo[i] = y[i]; }

        // --- Compute MSD and save to file ---
        FILE *fmsd;
        char fname[64];
        sprintf(fname, "msd_phi_%.2f.txt", phi);
        fmsd = fopen(fname, "w");

        for (int step = 0; step < N_MC_STEPS; step++) {

            for (int n = 0; n < N; n++) {
                int idx = rand() % N;
                trial_move(x, y, xu, yu, idx);
            }

            double msd = 0.0;
            for (int i = 0; i < N; i++) {
                double dx = xu[i] - xo[i];
                double dy = yu[i] - yo[i];
                msd += dx*dx + dy*dy;
            }
            msd /= (double)N;
            fprintf(fmsd, "%d\t%lf\n", step+1, msd);
        }
        fclose(fmsd);

        // --- Save final snapshot ---
        char snapname[64];
        sprintf(snapname, "snapshot_phi_%.2f.txt", phi);
        save_positions_to_file(x, y, snapname);

        printf("Snapshot and MSD written for phi = %.2f\n", phi);
    }

    printf("\nSimulation complete.\n");
    return 0;
    */

    // ================== CODE BLOCK 4 ==================
    /*
    srand(time(NULL));

    double g_values[]={0.0,0.01,0.10,1.0,10.0};
    int n_g = sizeof(g_values) / sizeof(g_values[0]);

    static double x[N], y[N],xu[N],yu[N];

    compute_box_size();
    printf("N=%d sigma=%.3f phi=%.4f -> Lx=%.6f Ly=%.6f\n", N, sigma, phi, Lx, Ly);

    // --- File to store energy ---
    char fname_E[128];
    sprintf(fname_E, "energy_g_%.3f.txt", g_field);
    FILE *fE = fopen(fname_E, "w");
    if (!fE) { perror("error opening energy file"); exit(1); }

    for (int gidx = 0; gidx < n_g; gidx++) {
        g_field = g_values[gidx];
        printf("\n=== Running for g = %.3f ===\n", g_field);

        // Initialize positions randomly in closed box
        initial_conditions_closed_box(x, y);
        for (int i = 0; i < N; i++) {
            xu[i] = x[i];
            yu[i] = y[i];
        }

        int total_accepts = 0, total_trials= 0; 

        // Equilibrate system
        for (int step = 0; step < N_MC_STEPS; step++) {
            for (int n = 0; n < N; n++) {
                int idx = rand() % N;
                total_trials++; 
                total_accepts+=trial_move_gravity_closed_box(x, y,idx);
            }

            // Compute energy every 1000 steps
            if (step % 1000 == 0) {
                double EperN = compute_energy_per_particle(y);
                fprintf(fE, "%d\t%.8f\n", step, EperN);
            }

            // Print acceptance rate occasionally
            if (step % (N_MC_STEPS / 10) == 0 && step > 0) {
                double acc = 100.0 * total_accepts / (double)total_trials;
                printf("Step %d: acceptance = %.2f%%\n", step, acc);
                total_accepts = 0;
                total_trials = 0;
            }
        }

        fclose(fE);

        // Save final snapshot
        char snapname[128];
        sprintf(snapname, "snapshot_g_%.3f_N=1000000.txt", g_field);
        save_positions_to_file(x, y, snapname);

        printf("Snapshot saved for g = %.3f\n", g_field);
    }

    printf("\nSimulation complete.\n");
    return 0;
    */

    // ================== CODE FOR SNAPSHOTS ==================
    srand(time(NULL));
    double phis_snap[] = {0.05, 0.2, 0.5};
    int n_phis_snap = 3;
    double pi_snap = 3.14159265358979323846;
    double area_particle_snap = pi_snap * (sigma/2.0) * (sigma/2.0);

    static double x_snap[N], y_snap[N], xo_snap[N], yo_snap[N], xu_snap[N], yu_snap[N];

    for (int p = 0; p < n_phis_snap; p++) {
        phi = phis_snap[p];
        L = sqrt((double)N * area_particle_snap / phi);
        printf("\n=== Running for phi = %.2f, box L = %.3f ===\n", phi, L);

        // Initialize triangular lattice
        triangular_lattice(x_snap, y_snap,1000,L);
        for (int i = 0; i < N; i++) { 
            xu_snap[i] = x_snap[i]; yu_snap[i] = y_snap[i]; 
            xo_snap[i] = x_snap[i]; yo_snap[i] = y_snap[i]; 
        }

        // --- Save initial snapshot ---
        char snapname[64];
        sprintf(snapname, "phi_%.2f_step_%d.txt", phi, 0);
        save_positions_to_file(x_snap, y_snap, snapname);

        // --- Monte Carlo dynamics ---
        for (int step = 0; step < N_MC_STEPS; step++) {
            for (int n = 0; n < N; n++) {
                int idx = rand() % N;
                trial_move(x_snap, y_snap, xu_snap, yu_snap, idx);
            }

            // --- Save mid and final snapshots ---
            if (step == N_MC_STEPS/2 - 1 || step == N_MC_STEPS - 1) {
                sprintf(snapname, "phi_%.2f_step_%d.txt", phi, step+1);
                save_positions_to_file(x_snap, y_snap, snapname);
            }
        }

        printf("Snapshots written for phi = %.2f\n", phi);
    }

    return 0;
}

// ---------------- BOX SIZE COMPUTATION ----------------
void compute_box_size(void) {
    double pi = acos(-1.0);
    double area_particle = pi * (sigma/2.0) * (sigma/2.0);
    double area_box = ((double)N * area_particle) / (double)phi;
    Lx = sqrt(area_box / 10.0);   // Ly = 10*Lx
    Ly = 10.0 * Lx;
}

// ---------------- INITIAL CONDITIONS ----------------
void initial_conditions(double *x, double *y) {
    int grid_size = (int)ceil(sqrt((double)N));
    double spacing = L / grid_size;
    int count = 0;

    for (int i = 0; i < grid_size && count < N; i++) {
        for (int j = 0; j < grid_size && count < N; j++) {
            x[count] = -L/2 + (i + 0.5) * spacing;
            y[count] = -L/2 + (j + 0.5) * spacing;
            count++;
        }
    }

    // Small random perturbation and wrapping into box
    for (int i = 0; i < N; i++) {
        x[i] += (random01() - 0.5) * spacing * 0.1;
        y[i] += (random01() - 0.5) * spacing * 0.1;
        if (x[i] > L/2) x[i] -= L;
        else if (x[i] < -L/2) x[i] += L;
        if (y[i] > L/2) y[i] -= L;
        else if (y[i] < -L/2) y[i] += L;
    }
}

void initial_conditions_closed_box(double *x, double *y) {
    int cols = (int)floor(Lx / sigma);
    if (cols < 1) cols = 1;
    int rows = (int)ceil((double)N / cols);

    double spacing_x = Lx / (double)cols;
    double spacing_y = Ly / (double)rows;

    int count = 0;
    for (int j = 0; j < rows && count < N; j++) {
        for (int i = 0; i < cols && count < N; i++) {
            // Centered positions
            x[count] = -Lx/2.0 + (i + 0.5) * spacing_x;
            y[count] = -Ly/2.0 + (j + 0.5) * spacing_y;
            count++;
        }
    }

    // Small perturbation without leaving box
    for (int i = 0; i < count; i++) {
        double dx = (random01() - 0.5) * (0.1 * spacing_x);
        double dy = (random01() - 0.5) * (0.1 * spacing_y);
        x[i] += dx; y[i] += dy;

        // Clamp to valid region
        if (x[i] < -Lx/2.0 + sigma/2.0) x[i] = -Lx/2.0 + sigma/2.0;
        if (x[i] >  Lx/2.0 - sigma/2.0) x[i] =  Lx/2.0 - sigma/2.0;
        if (y[i] < -Ly/2.0 + sigma/2.0) y[i] = -Ly/2.0 + sigma/2.0;
        if (y[i] >  Ly/2.0 - sigma/2.0) y[i] =  Ly/2.0 - sigma/2.0;
    }
}

// ---------------- TRIANGULAR LATTICE ----------------
void triangular_lattice(double *x, double *y, int Npart, double L){
    double a = L*sqrt(2.0/sqrt(3.0))/sqrt(N);
    double ax = a;
    double ay = a*sqrt(3.0)/2.0;
    int ny = (int)(L/ay)+1;
    int nx = (int)(L/ax)+1;

    int part = 0; 
    for(int i = 0; i<ny && part<N ; i++){
        double yi = -L/2.0+ay*i;
        for(int j = 0; j<nx && part<N ; j++){
            double xi = -L/2.0 + ax*j + (i%2)*ax/2.0;
            x[part] = xi;
            y[part] = yi;
            part++;
        }
    }
}

// ---------------- PBC FUNCTIONS ----------------
double wrap(double x, double L) {
    if (L <= 0.0) return x;
    x = fmod(x + L / 2.0, L);
    if (x < 0.0) x += L;
    return x - L / 2.0;
}

double dist_pbc(double dx) {
    if (dx > L/2) dx -= L;
    else if (dx < -L/2) dx += L;
    return dx;
}

double distance_pbc(double x1, double y1, double x2, double y2){
    double dx = dist_pbc(x1 - x2);
    double dy = dist_pbc(y1 - y2);
    return sqrt(dx*dx + dy*dy);
}

// ---------------- TRIAL MOVE ----------------
int trial_move(double *x, double *y, double *xu, double *yu, int i) {
    double old_xu = xu[i];
    double old_yu = yu[i];

    double dx_prop = (random01() - 0.5) * delta;
    double dy_prop = (random01() - 0.5) * delta;
    double new_xu = old_xu + dx_prop;
    double new_yu = old_yu + dy_prop;

    double new_x = wrap(new_xu,L);
    double new_y = wrap(new_yu,L);

    double sigma2 = sigma * sigma;
    for (int j = 0; j < N; j++) {
        if (i == j) continue;
        double dx = new_x - x[j];
        if (dx > L/2) dx -= L;
        else if (dx < -L/2) dx += L;
        double dy = new_y - y[j];
        if (dy > L/2) dy -= L;
        else if (dy < -L/2) dy += L;

        if (dx*dx + dy*dy <= sigma2) return 0; // Reject if overlapping
    }

    xu[i] = new_xu;
    yu[i] = new_yu;
    x[i] = new_x;
    y[i] = new_y;
    return 1;
}

// ---------------- TRIAL MOVE WITH GRAVITY ----------------
int trial_move_gravity_closed_box(double *x, double *y, int i) {
    double old_x = x[i];
    double old_y = y[i];

    double dx = (random01() - 0.5) * delta;
    double dy = (random01() - 0.5) * delta;
    double new_x = old_x + dx;
    double new_y = old_y + dy;

    // Check wall collisions
    if (new_x - sigma/2.0 < -Lx/2.0) return 0;
    if (new_x + sigma/2.0 >  Lx/2.0) return 0;
    if (new_y - sigma/2.0 < -Ly/2.0) return 0;
    if (new_y + sigma/2.0 >  Ly/2.0) return 0;

    // Check overlaps
    double sigma2 = sigma * sigma;
    for (int j = 0; j < N; j++) {
        if (i == j) continue;
        double dxij = new_x - x[j];
        double dyij = new_y - y[j];
        if (dxij*dxij + dyij*dyij < sigma2) return 0;
    }

    double dE = mass * g_field * (new_y - old_y);
    if (dE <= 0.0) {
        x[i] = new_x; y[i] = new_y;
        return 1;
    } else {
        if (random01() < exp(-dE / kT)) {
            x[i] = new_x; y[i] = new_y;
            return 1;
        } else return 0;
    }
}

// ---------------- SAVE POSITIONS ----------------
void save_positions_to_file(double *x, double *y, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) { perror("Error opening"); return; }
    for (int i = 0; i < N; i++) fprintf(fp, "%lf %lf\n", x[i], y[i]);
    fclose(fp);
}

// ---------------- ENERGY PER PARTICLE ----------------
double compute_energy_per_particle(double *y) {
    double Etot = 0.0;
    for (int i = 0; i < N; i++) Etot += mass * g_field * y[i];
    return Etot / (double)N;
}