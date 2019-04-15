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
extern SEXP _phenofit_sgmat_S(SEXP, SEXP);
extern SEXP _phenofit_sgmat_B(SEXP);
extern SEXP _phenofit_sgmat_wB(SEXP, SEXP);
extern SEXP _phenofit_smooth_SG(SEXP, SEXP, SEXP);
extern SEXP _phenofit_smooth_wSG(SEXP, SEXP, SEXP, SEXP);
extern SEXP _phenofit_wTSM_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _phenofit_movmean(SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"smooth2", (DL_FUNC) &smooth2, 8},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_phenofit_fix_dt"    , (DL_FUNC) &_phenofit_fix_dt, 1},
    {"_phenofit_sgmat_S"   , (DL_FUNC) &_phenofit_sgmat_S, 2},
    {"_phenofit_sgmat_B"   , (DL_FUNC) &_phenofit_sgmat_B, 1},
    {"_phenofit_sgmat_wB"  , (DL_FUNC) &_phenofit_sgmat_wB, 2},
    {"_phenofit_smooth_wSG", (DL_FUNC) &_phenofit_smooth_wSG, 4},
    {"_phenofit_smooth_SG" , (DL_FUNC) &_phenofit_smooth_SG, 3},
    {"_phenofit_wTSM_cpp"  , (DL_FUNC) &_phenofit_wTSM_cpp, 6},
    {"_phenofit_movmean", (DL_FUNC) &_phenofit_movmean, 4},
    {NULL, NULL, 0}
};

void R_init_phenofit(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
