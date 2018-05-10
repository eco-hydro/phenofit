// #include "whit.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP _phenofit_sgolayB(SEXP SSEXP, SEXP wSEXP);
SEXP _phenofit_sgfitw_rcpp(SEXP ySEXP, SEXP wSEXP, SEXP SSEXP);
SEXP _phenofit_wTSM_cpp(SEXP ySEXP, SEXP yfitSEXP, SEXP wSEXP, SEXP itersSEXP, SEXP nptperyearSEXP, SEXP wfactSEXP);
void smooth2(double * w, double * y, double * z, double * lamb, int * mm,
    double * d, double * c, double * e);

static const R_CallMethodDef CallEntries[] = {
    {"_phenofit_sgolayB", (DL_FUNC) &_phenofit_sgolayB, 2},
    {"_phenofit_sgfitw_rcpp", (DL_FUNC) &_phenofit_sgfitw_rcpp, 3},
    {"_phenofit_wTSM_cpp", (DL_FUNC) &_phenofit_wTSM_cpp, 6},
    {NULL, NULL, 0}
};

static const R_CMethodDef cEntries[] = {
   {"smooth2", (DL_FUNC) &smooth2, 8},
   {NULL, NULL, 0}
};

void R_init_phenofit(DllInfo *dll) {
    // R_RegisterCCallable("phenofit", "smooth2", (DL_FUNC) &smooth2),
    R_registerRoutines(dll, cEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
