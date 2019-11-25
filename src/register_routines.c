#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/
// https://stackoverflow.com/questions/42313373/r-cmd-check-note-found-no-calls-to-r-registerroutines-r-usedynamicsymbols

/* .C calls */
extern void smooth2(void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _phenofit_check_season(SEXP, SEXP, SEXP, SEXP);
extern SEXP _phenofit_sgmat_S(SEXP, SEXP);
extern SEXP _phenofit_sgmat_B(SEXP);
extern SEXP _phenofit_sgmat_wB(SEXP, SEXP);
extern SEXP _phenofit_rcpp_SG(SEXP, SEXP, SEXP);
extern SEXP _phenofit_rcpp_wSG(SEXP, SEXP, SEXP, SEXP);
extern SEXP _phenofit_rcpp_wTSM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _phenofit_movmean(SEXP, SEXP, SEXP, SEXP);

extern SEXP _phenofit_f_goal_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _phenofit_clogistic(SEXP, SEXP, SEXP);
extern SEXP _phenofit_cdoubleLog_Zhang(SEXP, SEXP, SEXP);
extern SEXP _phenofit_cdoubleLog_AG(SEXP, SEXP, SEXP);
extern SEXP _phenofit_cdoubleLog_Beck(SEXP, SEXP, SEXP);
extern SEXP _phenofit_cdoubleLog_Elmore(SEXP, SEXP, SEXP);
extern SEXP _phenofit_cdoubleLog_Gu(SEXP, SEXP, SEXP);
extern SEXP _phenofit_cdoubleLog_Klos(SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    CALLDEF(smooth2, 8),
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(_phenofit_check_season, 4),
    CALLDEF(_phenofit_sgmat_S   , 2),
    CALLDEF(_phenofit_sgmat_B   , 1),
    CALLDEF(_phenofit_sgmat_wB  , 2),
    CALLDEF(_phenofit_rcpp_wSG, 4),
    CALLDEF(_phenofit_rcpp_SG , 3),
    CALLDEF(_phenofit_rcpp_wTSM  , 6),
    CALLDEF(_phenofit_movmean   , 4),

    CALLDEF(_phenofit_f_goal_cpp, 7),
    CALLDEF(_phenofit_clogistic        , 3),
    CALLDEF(_phenofit_cdoubleLog_Zhang , 3),
    CALLDEF(_phenofit_cdoubleLog_AG    , 3),
    CALLDEF(_phenofit_cdoubleLog_Beck  , 3),
    CALLDEF(_phenofit_cdoubleLog_Elmore, 3),
    CALLDEF(_phenofit_cdoubleLog_Gu    , 3),
    CALLDEF(_phenofit_cdoubleLog_Klos  , 3),
    {NULL, NULL, 0}
};

void R_init_phenofit(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
