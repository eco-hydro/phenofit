// #ifndef WHIT_H
// #define WHIT_H

/* Smoothing and interpolation with second differences.
   Input:  weights (w), data (y): vector from 1 to m.
   Input:  smoothing parameter (lambda), length (m).
   Output: smoothed vector (z): vector from 1 to m. 
   */
void smooth2(double * w, double * y, double * z, double * lamb, int * mm,
    double * d, double * c, double * e)
{
    int i, i1, i2, m;
    double lambda;
    lambda = * lamb;
    m = * mm - 1;
    d[0] = w[0] + lambda;
    c[0] = -2 * lambda / d[0];
    e[0] = lambda / d[0];
    z[0] = w[0] * y[0];
    d[1] = w[1] + 5 * lambda - d[0] * c[0] * c[0];
    c[1] = (-4 * lambda - d[0] * c[0] * e[0]) / d[1];
    e[1] = lambda / d[1];
    z[1] = w[1] * y[1] - c[0] * z[0];
    for (i = 2; i < m - 1; i++) {
        i1 = i - 1;
        i2 = i - 2;
        d[i] = w[i] + 6 * lambda - c[i1] * c[i1] * d[i1] - e[i2] * e[i2] * d[i2];
        c[i] = (-4 * lambda - d[i1] * c[i1] * e[i1]) / d[i];
        e[i] = lambda / d[i];
        z[i] = w[i] * y[i] - c[i1] * z[i1] - e[i2] * z[i2];
    };
    i1 = m - 2;
    i2 = m - 3;
    d[m - 1] = w[m - 1] + 5 * lambda - c[i1] * c[i1] * d[i1] - e[i2] * e[i2] * d[i2];
    c[m - 1] = (-2 * lambda - d[i1] * c[i1] * e[i1]) / d[m - 1];
    z[m - 1] = w[m - 1] * y[m - 1] - c[i1] * z[i1] - e[i2] * z[i2];
    i1 = m - 1;
    i2 = m - 2;
    d[m] = w[m] + lambda - c[i1] * c[i1] * d[i1] - e[i2] * e[i2] * d[i2];
    z[m] = (w[m] * y[m] - c[i1] * z[i1] - e[i2] * z[i2]) / d[m];
    z[m - 1] = z[m - 1] / d[m - 1] - c[m - 1] * z[m];
    for (i = m - 2; 0 <= i; i--)
        z[i] = z[i] / d[i] - c[i] * z[i + 1] - e[i] * z[i + 2];
}

// #endif
