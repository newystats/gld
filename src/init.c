#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void fit_fkml(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gl_fm5_distfunc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gl_fmkl_distfunc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gl_rs_distfunc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gl_vsk_distfunc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"fit_fkml",         (DL_FUNC) &fit_fkml,         15},
    {"gl_fm5_distfunc",  (DL_FUNC) &gl_fm5_distfunc,  12},
    {"gl_fmkl_distfunc", (DL_FUNC) &gl_fmkl_distfunc, 11},
    {"gl_rs_distfunc",   (DL_FUNC) &gl_rs_distfunc,   11},
    {"gl_vsk_distfunc",  (DL_FUNC) &gl_vsk_distfunc,  11},
    {NULL, NULL, 0}
};

void R_init_gld(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
