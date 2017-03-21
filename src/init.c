#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP TauStar_eigenForDiscreteProbs(SEXP);
extern SEXP TauStar_HoeffIndCdfRCPP(SEXP, SEXP);
extern SEXP TauStar_HoeffIndDiscreteCdfRCPP(SEXP, SEXP, SEXP, SEXP);
extern SEXP TauStar_HoeffIndDiscretePdfRCPP(SEXP, SEXP, SEXP, SEXP);
extern SEXP TauStar_HoeffIndMixedCdfRCPP(SEXP, SEXP, SEXP);
extern SEXP TauStar_HoeffIndMixedPdfRCPP(SEXP, SEXP, SEXP);
extern SEXP TauStar_HoeffIndPdfRCPP(SEXP, SEXP);
extern SEXP TauStar_RcppExport_registerCCallable();
extern SEXP TauStar_TStarFastResampleRCPP(SEXP, SEXP, SEXP, SEXP);
extern SEXP TauStar_TStarHellerAndHellerRCPP(SEXP, SEXP);
extern SEXP TauStar_TStarNaiveRCPP(SEXP, SEXP, SEXP);
extern SEXP TauStar_TStarWeihsEtAlRCPP(SEXP, SEXP);
extern SEXP TauStar_VTStarHellerAndHellerRCPP(SEXP, SEXP);
extern SEXP TauStar_VTStarWeihsEtAlRCPP(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"TauStar_eigenForDiscreteProbs",        (DL_FUNC) &TauStar_eigenForDiscreteProbs,        1},
    {"TauStar_HoeffIndCdfRCPP",              (DL_FUNC) &TauStar_HoeffIndCdfRCPP,              2},
    {"TauStar_HoeffIndDiscreteCdfRCPP",      (DL_FUNC) &TauStar_HoeffIndDiscreteCdfRCPP,      4},
    {"TauStar_HoeffIndDiscretePdfRCPP",      (DL_FUNC) &TauStar_HoeffIndDiscretePdfRCPP,      4},
    {"TauStar_HoeffIndMixedCdfRCPP",         (DL_FUNC) &TauStar_HoeffIndMixedCdfRCPP,         3},
    {"TauStar_HoeffIndMixedPdfRCPP",         (DL_FUNC) &TauStar_HoeffIndMixedPdfRCPP,         3},
    {"TauStar_HoeffIndPdfRCPP",              (DL_FUNC) &TauStar_HoeffIndPdfRCPP,              2},
    {"TauStar_RcppExport_registerCCallable", (DL_FUNC) &TauStar_RcppExport_registerCCallable, 0},
    {"TauStar_TStarFastResampleRCPP",        (DL_FUNC) &TauStar_TStarFastResampleRCPP,        4},
    {"TauStar_TStarHellerAndHellerRCPP",     (DL_FUNC) &TauStar_TStarHellerAndHellerRCPP,     2},
    {"TauStar_TStarNaiveRCPP",               (DL_FUNC) &TauStar_TStarNaiveRCPP,               3},
    {"TauStar_TStarWeihsEtAlRCPP",           (DL_FUNC) &TauStar_TStarWeihsEtAlRCPP,           2},
    {"TauStar_VTStarHellerAndHellerRCPP",    (DL_FUNC) &TauStar_VTStarHellerAndHellerRCPP,    2},
    {"TauStar_VTStarWeihsEtAlRCPP",          (DL_FUNC) &TauStar_VTStarWeihsEtAlRCPP,          2},
    {NULL, NULL, 0}
};

void R_init_TauStar(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
