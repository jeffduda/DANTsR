#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP antsrMesh(SEXP, SEXP, SEXP, SEXP);
extern SEXP antsrMesh_AddPoint(SEXP, SEXP, SEXP, SEXP);
extern SEXP antsrMesh_AddPolyline(SEXP, SEXP, SEXP, SEXP);
extern SEXP antsrMesh_GetCell(SEXP, SEXP, SEXP);
extern SEXP antsrMesh_GetCellPoints(SEXP, SEXP, SEXP);
extern SEXP antsrMesh_GetNumberOfCells(SEXP, SEXP);
extern SEXP antsrMesh_GetNumberOfPoints(SEXP, SEXP);
extern SEXP antsrMesh_GetPoint(SEXP, SEXP, SEXP);
extern SEXP antsrMesh_GetPoints(SEXP, SEXP, SEXP);
extern SEXP antsrMesh_ReadCamino(SEXP, SEXP, SEXP);
extern SEXP antsrMesh_ReadVTK(SEXP, SEXP, SEXP, SEXP);
extern SEXP antsrMesh_SetPoint(SEXP, SEXP, SEXP, SEXP);
extern SEXP antsrMesh_TransformMesh(SEXP, SEXP, SEXP);
extern SEXP antsrMesh_WriteCamino(SEXP, SEXP, SEXP, SEXP);
extern SEXP antsrMesh_WriteVTK(SEXP, SEXP, SEXP, SEXP);
extern SEXP antsrRegistrationRun(SEXP, SEXP, SEXP, SEXP);
extern SEXP deterministicTracking(SEXP, SEXP, SEXP, SEXP);
extern SEXP dtiFilters(SEXP, SEXP, SEXP);
extern SEXP dtiReconstruction(SEXP, SEXP, SEXP, SEXP);
extern SEXP interpolateImageValues(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP isInImage(SEXP, SEXP, SEXP);
extern SEXP labelsToPoints(SEXP, SEXP, SEXP);
extern SEXP pointCountImage(SEXP, SEXP, SEXP);
extern SEXP polylineCountImage(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"antsrMesh",                   (DL_FUNC) &antsrMesh,                   4},
    {"antsrMesh_AddPoint",          (DL_FUNC) &antsrMesh_AddPoint,          4},
    {"antsrMesh_AddPolyline",       (DL_FUNC) &antsrMesh_AddPolyline,       4},
    {"antsrMesh_GetCell",           (DL_FUNC) &antsrMesh_GetCell,           3},
    {"antsrMesh_GetCellPoints",     (DL_FUNC) &antsrMesh_GetCellPoints,     3},
    {"antsrMesh_GetNumberOfCells",  (DL_FUNC) &antsrMesh_GetNumberOfCells,  2},
    {"antsrMesh_GetNumberOfPoints", (DL_FUNC) &antsrMesh_GetNumberOfPoints, 2},
    {"antsrMesh_GetPoint",          (DL_FUNC) &antsrMesh_GetPoint,          3},
    {"antsrMesh_GetPoints",         (DL_FUNC) &antsrMesh_GetPoints,         3},
    {"antsrMesh_ReadCamino",        (DL_FUNC) &antsrMesh_ReadCamino,        3},
    {"antsrMesh_ReadVTK",           (DL_FUNC) &antsrMesh_ReadVTK,           4},
    {"antsrMesh_SetPoint",          (DL_FUNC) &antsrMesh_SetPoint,          4},
    {"antsrMesh_TransformMesh",     (DL_FUNC) &antsrMesh_TransformMesh,     3},
    {"antsrMesh_WriteCamino",       (DL_FUNC) &antsrMesh_WriteCamino,       4},
    {"antsrMesh_WriteVTK",          (DL_FUNC) &antsrMesh_WriteVTK,          4},
    {"antsrRegistrationRun",        (DL_FUNC) &antsrRegistrationRun,        4},
    {"deterministicTracking",       (DL_FUNC) &deterministicTracking,       4},
    {"dtiFilters",                  (DL_FUNC) &dtiFilters,                  3},
    {"dtiReconstruction",           (DL_FUNC) &dtiReconstruction,           4},
    {"interpolateImageValues",      (DL_FUNC) &interpolateImageValues,      5},
    {"isInImage",                   (DL_FUNC) &isInImage,                   4},
    {"labelsToPoints",              (DL_FUNC) &labelsToPoints,              3},
    {"polylineCountImage",          (DL_FUNC) &pointCountImage,             3},
    {"pointCountImage",             (DL_FUNC) &pointCountImage,             3},
    {NULL, NULL, 0}
};

void R_init_DANTsR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
