#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void C_isomean(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"C_isomean", (DL_FUNC) &C_isomean, 4},
    {NULL, NULL, 0}
};

void R_init_fdrtool(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
