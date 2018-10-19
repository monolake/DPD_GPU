#include "sim_run.h"

void run() {
}

void compute_force(double** r, double** v, double** f, double rand_num, int N) {

    double m = 1.0;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
                if(i != j) {
                    double delx, dely, delz;
                    delx = r[i][0] - r[j][0];
                    dely = r[i][1] - r[j][1];
                    delz = r[i][2] - r[j][2];
                    double r;
                    r = sqrt(delx * delx + dely * dely + delz * delz);
                    if(r < rc) {
                        double fpair;
                        double wr;
                        wr = 1 - r / rc;
                        fpair = force_a0 * wr;

                        double delvx, delvy, delvz;
                        delvx = v[i][0] - v[j][0];
                        delvy = v[i][1] - v[j][1];
                        delvz = v[i][2] - v[j][2];

                        double dot;
                        dot = (delx * delvx + dely * delvy + delz * delvz) / r;
                        fpair -= force_gamma * wr * wr * dot;

                        //randum - gaussian random number with zero mean and unit variance
                        fpair += force_sigma * wr * rand_num * 1 / sqrt(dt);

                        f[i][0] += delx / r * fpair;
                        f[i][1] += dely / r * fpair;
                        f[i][2] += delz / r * fpair;
                    }

                }

        }
    }
}

