#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP antsrMesh(SEXP, SEXP, SEXP, SEXP);
extern SEXP antsrMesh_AddPoint(SEXP, SEXP, SEXP);
extern SEXP antsrMesh_AddPolyline(SEXP, SEXP, SEXP);
extern SEXP antsrMesh_GetCell(SEXP, SEXP);
extern SEXP antsrMesh_GetCellPoints(SEXP, SEXP);
extern SEXP antsrMesh_GetNumberOfCells(SEXP);
extern SEXP antsrMesh_GetNumberOfPoints(SEXP);
extern SEXP antsrMesh_GetPoint(SEXP, SEXP);
extern SEXP antsrMesh_GetPoints(SEXP, SEXP);
extern SEXP antsrMesh_ReadCamino(SEXP, SEXP);
extern SEXP antsrMesh_ReadTrk(SEXP,SEXP);
extern SEXP antsrMesh_ReadVTK(SEXP, SEXP, SEXP);
extern SEXP antsrMesh_SetPoint(SEXP, SEXP, SEXP);
extern SEXP antsrMesh_TransformMesh(SEXP, SEXP, SEXP);
extern SEXP antsrMesh_WriteCamino(SEXP, SEXP, SEXP);
extern SEXP antsrMesh_WriteTrk(SEXP, SEXP, SEXP);
extern SEXP antsrMesh_WriteVTK(SEXP, SEXP, SEXP);
extern SEXP antsrRegistrationRun(SEXP, SEXP, SEXP);
extern SEXP deterministicTracking(SEXP, SEXP, SEXP);
extern SEXP dtiFilters(SEXP, SEXP);
extern SEXP dtiReconstruction(SEXP, SEXP, SEXP);
extern SEXP interpolateImageValues(SEXP, SEXP, SEXP, SEXP);
extern SEXP isInImage(SEXP, SEXP, SEXP);
extern SEXP isInMesh(SEXP, SEXP, SEXP);
extern SEXP labelsToPoints(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"antsrMesh",                   (DL_FUNC) &antsrMesh,                   4},
    {"antsrMesh_AddPoint",          (DL_FUNC) &antsrMesh_AddPoint,          3},
    {"antsrMesh_AddPolyline",       (DL_FUNC) &antsrMesh_AddPolyline,       3},
    {"antsrMesh_GetCell",           (DL_FUNC) &antsrMesh_GetCell,           2},
    {"antsrMesh_GetCellPoints",     (DL_FUNC) &antsrMesh_GetCellPoints,     2},
    {"antsrMesh_GetNumberOfCells",  (DL_FUNC) &antsrMesh_GetNumberOfCells,  1},
    {"antsrMesh_GetNumberOfPoints", (DL_FUNC) &antsrMesh_GetNumberOfPoints, 1},
    {"antsrMesh_GetPoint",          (DL_FUNC) &antsrMesh_GetPoint,          2},
    {"antsrMesh_GetPoints",         (DL_FUNC) &antsrMesh_GetPoints,         2},
    {"antsrMesh_ReadCamino",        (DL_FUNC) &antsrMesh_ReadCamino,        2},
    {"antsrMesh_ReadTrk",           (DL_FUNC) &antsrMesh_ReadTrk,           2},
    {"antsrMesh_ReadVTK",           (DL_FUNC) &antsrMesh_ReadVTK,           3},
    {"antsrMesh_SetPoint",          (DL_FUNC) &antsrMesh_SetPoint,          3},
    {"antsrMesh_TransformMesh",     (DL_FUNC) &antsrMesh_TransformMesh,     3},
    {"antsrMesh_WriteCamino",       (DL_FUNC) &antsrMesh_WriteCamino,       3},
    {"antsrMesh_WriteTrk",          (DL_FUNC) &antsrMesh_WriteTrk,          3},
    {"antsrMesh_WriteVTK",          (DL_FUNC) &antsrMesh_WriteVTK,          3},
    {"antsrRegistrationRun",        (DL_FUNC) &antsrRegistrationRun,        3},
    {"deterministicTracking",       (DL_FUNC) &deterministicTracking,       3},
    {"dtiFilters",                  (DL_FUNC) &dtiFilters,                  2},
    {"dtiReconstruction",           (DL_FUNC) &dtiReconstruction,           3},
    {"interpolateImageValues",      (DL_FUNC) &interpolateImageValues,      4},
    {"isInImage",                   (DL_FUNC) &isInImage,                   3},
    {"isInMesh",                    (DL_FUNC) &isInMesh,                    3},
    {"labelsToPoints",              (DL_FUNC) &labelsToPoints,              2},
    {NULL, NULL, 0}
};

void R_init_DANTsR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
