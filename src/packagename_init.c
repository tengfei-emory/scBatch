#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _scBatch_derif(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _scBatch_scBatchCpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_scBatch_derif",      (DL_FUNC) &_scBatch_derif,      5},
  {"_scBatch_scBatchCpp", (DL_FUNC) &_scBatch_scBatchCpp, 9},
  {NULL, NULL, 0}
};

void R_init_scBatch(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
