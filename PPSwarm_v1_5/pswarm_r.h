#include <stdio.h>
#include <math.h>
#ifdef linux
#include <strings.h>
#endif
#include <signal.h>
#include <setjmp.h>
#include <R.h>
#include <Rdefines.h>

#include "pswarm.h"

SEXP createRRealVector(int n, double* x);
SEXP createRIntScalar(int x);
SEXP pswarm_r(SEXP r_problem_ptr, SEXP r_options_ptr, SEXP environment_ptr);
void r_objfun(int n, int m, double *x, double *lb, double *ub, double *fxs);
double r_outfcn(int n, int s, int iter, int gbest, struct swarm *pop);

SEXP r_environment, r_objf, r_outf;

extern struct Options opt;
extern struct Stats stats;
extern int PSwarm(int n, void (*objf)(), double *lb, double *ub, int lincons, double *A,
          double *b, double **sol, double *f, double *x);

static int saved_options=0;
static struct Options opt_backup;


