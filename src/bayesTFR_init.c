#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void doDLcurve(void *, void *, void *, void *, void *, void *);
extern void log_cond_Triangle_c4_trans(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"doDLcurve",                  (DL_FUNC) &doDLcurve,                  6},
  {"log_cond_Triangle_c4_trans", (DL_FUNC) &log_cond_Triangle_c4_trans, 8},
  {NULL, NULL, 0}
};

void R_init_bayesTFR(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
