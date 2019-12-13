#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "marqLevAlg.h"

static R_FortranMethodDef FortRout[] = {
  {"dsinv", (DL_FUNC) &F77_SUB(dsinv), 5},
  {"dchole", (DL_FUNC) &F77_SUB(dchole), 4},
  {NULL, NULL, 0}
};

void R_init_marqLevAlg(DllInfo * dll)
{
  R_registerRoutines(dll, NULL, NULL, FortRout, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
