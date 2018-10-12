#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct Comp {
    /* comp of the form: a + bi */
    double a, b;
} Comp;

Comp comp_create(double a, double b) {
    Comp res;
    res.a = a;
    res.b = b;
    return res;
}

void comp_print(Comp comp) {
    printf("%.6f + %.6f i\n", comp.a, comp.b);
}


#define PI 3.141592653589793


/* Calculate e^(ix) */
Comp comp_exp(double x) {
    Comp res;
    res.a = cos(x);
    res.b = sin(x);
    return res;
}

#define comp_mul_self(c, c2) \
do { \
    double ca = c->a; \
    c->a = ca * c2->a - c->b * c2->b; \
    c->b = c->b * c2->a + ca * c2->b; \
} while (0)

void fft(const Comp *sig, Comp *f, int s, int n, int inv) {
    int i, hn = n >> 1;
    Comp ep = comp_exp((inv ? PI : -PI) / (double)hn), ei;
    Comp *pi = &ei, *pp = &ep;
    if (!hn) *f = *sig;
    else
    {
        fft(sig, f, s << 1, hn, inv);
        fft(sig + s, f + hn, s << 1, hn, inv);
        pi->a = 1;
        pi->b = 0;
        for (i = 0; i < hn; i++)
        {
            Comp even = f[i], *pe = f + i, *po = pe + hn;
            comp_mul_self(po, pi);
            pe->a += po->a;
            pe->b += po->b;
            po->a = even.a - po->a;
            po->b = even.b - po->b;
            comp_mul_self(pi, pp);
        }
    }
}
void test_fft(const Comp *sig, Comp *f, int n) {
    int i;
    puts("## ---------------- FFT ---------------- ##");
    fft(sig, f, 1, n, 0);
    for (i = 0; i < n; i++)
        comp_print(f[i]);
}

int main() {
    int n, i, k;
    Comp *sig, *f;
    printf("provide the value of k for n=2^k\n");
    scanf("%d", &k);
    n = 1 << k;
    sig = (Comp *)malloc(sizeof(Comp) * (size_t)n);
    f = (Comp *)malloc(sizeof(Comp) * (size_t)n);
    for (i = 0; i < n; i++)
    {
        sig[i].a = ((double)rand()/(double)(RAND_MAX)) * 10.0;
        sig[i].b = ((double)rand()/(double)(RAND_MAX)) * 10.0;
    }
    puts("## ------------ Original Signal ------------ ##");
    for (i = 0; i < n; i++)
        comp_print(sig[i]);
    puts("---------------------------------------------");
    
    test_fft(sig, f, n);

    return 0;
}
