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
extern SEXP _phenofit_fix_dt(SEXP);
extern SEXP _phenofit_sgmat_S(SEXP, SEXP);
extern SEXP _phenofit_sgmat_B(SEXP);
extern SEXP _phenofit_sgmat_wB(SEXP, SEXP);
extern SEXP _phenofit_smooth_SG(SEXP, SEXP, SEXP);
extern SEXP _phenofit_smooth_wSG(SEXP, SEXP, SEXP, SEXP);
extern SEXP _phenofit_wTSM_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _phenofit_movmean(SEXP, SEXP, SEXP, SEXP);

extern SEXP _phenofit_f_goal_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _phenofit_logistic(SEXP, SEXP);
extern SEXP _phenofit_doubleLog_Zhang(SEXP, SEXP);
extern SEXP _phenofit_doubleLog_AG(SEXP, SEXP);
extern SEXP _phenofit_doubleLog_Beck(SEXP, SEXP);
extern SEXP _phenofit_doubleLog_Elmore(SEXP, SEXP);
extern SEXP _phenofit_doubleLog_Gu(SEXP, SEXP);
extern SEXP _phenofit_doubleLog_Klos(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    CALLDEF(smooth2, 8),
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(_phenofit_fix_dt    , 1),
    CALLDEF(_phenofit_sgmat_S   , 2),
    CALLDEF(_phenofit_sgmat_B   , 1),
    CALLDEF(_phenofit_sgmat_wB  , 2),
    CALLDEF(_phenofit_smooth_wSG, 4),
    CALLDEF(_phenofit_smooth_SG , 3),
    CALLDEF(_phenofit_wTSM_cpp  , 6),
    CALLDEF(_phenofit_movmean   , 4),

    CALLDEF(_phenofit_f_goal_cpp, 6),
    CALLDEF(_phenofit_logistic        , 2),
    CALLDEF(_phenofit_doubleLog_Zhang , 2),
    CALLDEF(_phenofit_doubleLog_AG    , 2),
    CALLDEF(_phenofit_doubleLog_Beck  , 2),
    CALLDEF(_phenofit_doubleLog_Elmore, 2),
    CALLDEF(_phenofit_doubleLog_Gu    , 2),
    CALLDEF(_phenofit_doubleLog_Klos  , 2),
    {NULL, NULL, 0}
};

void R_init_phenofit(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
