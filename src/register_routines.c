#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/
// https://stackoverflow.com/questions/42313373/r-cmd-check-note-found-no-calls-to-r-registerroutines-r-usedynamicsymbols

/* .C calls */
extern void smooth2(void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _phenofit_fix_dt(SEXP);
extern SEXP _phenofit_sgfitw_rcpp(SEXP, SEXP, SEXP);
extern SEXP _phenofit_sgolayB(SEXP, SEXP);
extern SEXP _phenofit_wTSM_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"smooth2", (DL_FUNC) &smooth2, 8},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_phenofit_fix_dt",      (DL_FUNC) &_phenofit_fix_dt,      1},
    {"_phenofit_sgfitw_rcpp", (DL_FUNC) &_phenofit_sgfitw_rcpp, 3},
    {"_phenofit_sgolayB",     (DL_FUNC) &_phenofit_sgolayB,     2},
    {"_phenofit_wTSM_cpp",    (DL_FUNC) &_phenofit_wTSM_cpp,    6},
    {NULL, NULL, 0}
};

void R_init_phenofit(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
