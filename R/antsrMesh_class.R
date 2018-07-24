# this file defines the S4 classes related to 'antsrMesh and the associated
# methods



#' @rdname antsrMesh
#' @title antsrMesh
#'
#' @description class for point sets and meshes
#' @param object input object to convert
#' @param .Object input object to convert
#' @param precision string e.g. "float" or "double"
#' @param dimension dimensionality of the transform (2,3,or 4)
#' @slot dimension usually 2 or 3 but can be 4
#' @slot precision math precision is float or double'
#' @slot pointer to the memory location of the itk object
setClass(Class = "antsrMesh",
         representation(precision= "character", dimension = "integer",
         pointer = "externalptr"))

#' @rdname antsrMesh
#' @aliases show,antsrMesh-method
setMethod(f = "show", "antsrMesh", function(object){
    cat("antsrMesh\n")
    cat("  Dimension :", object@dimension, "\n")
    cat("  Precision :", object@precision, "\n")
    cat("  Points    :", antsrMeshGetNumberOfPoints(object), "\n")
    cat("  Cells     :", antsrMeshGetNumberOfCells(object), "\n")
    cat("\n")
})

#' @rdname antsrMesh
#' @aliases initialize,antsrTransform-method
#' @param reserve number of points to allocate on creation
setMethod(f = "initialize", signature(.Object = "antsrMesh"), definition = function(.Object,
  dimension = 3, precision = "float", reserve=0) {
  mesh = .Call("antsrMesh", precision, dimension, reserve, matrix(0), PACKAGE = "DANTsR")
  return( mesh )
})


#' @title antsrMeshCreate
#' @description create a mesh
#' @param dimension number of dimensions
#' @param precision use 'float' or 'double' for values
#' @param reserve number of points to allocate on creation
#' @param points matrix of points in mesh
#' @return antsrMesh
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=128 )
#' @export
antsrMeshCreate <- function(dimension=3, precision="float", reserve=0, points=NULL)
{

  # Check for valid dimension
  if ( (dimension < 2) | (dimension > 4) )
  {
    stop(paste("Unsupported dimension:", dimension))
  }

  if ( (precision != "float") & (precision != "double"))
  {
    stop(paste("Unsupported precision:", precision))
  }

  if ( reserve < 0 ) {
    stop(paste("Unsupported reserve number:", reserve))
  }

  if ( !is.null(points) )  {
    if ( reserve == 0 ) {
      reserve = dim(points)[1]
    }
    else if ( reserve < dim(points)[1] ) {
      stop( "reserve must be >= number of points passed")
    }

    if ( dim(points)[2] != dimension ) {
      stop( "dimension and point size matrix don't match")
    }
  }
  else {
    points = matrix(0)
  }

  mesh = .Call("antsrMesh", precision, dimension, reserve, points, PACKAGE = "DANTsR")

  return(mesh)
  }

#' @title antsrMeshGetNumberOfPoints
#' @description get number of points in an antsrMesh
#' @param mesh an 'antsrMesh'
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=128 )
#' antsrMeshGetNumberOfPoints(x)
#' @export
  antsrMeshGetNumberOfPoints = function( mesh) {
    .Call("antsrMesh_GetNumberOfPoints", mesh, package="DANTsR")
  }

#' @title antsrMeshGetNumberOfCells
#' @description get number of cells in an antsrMesh
#' @param mesh an 'antsrMesh'
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=128 )
#' antsrMeshGetNumberOfCells(x)
#' @export
  antsrMeshGetNumberOfCells = function( mesh) {
    .Call("antsrMesh_GetNumberOfCells", mesh, package="DANTsR")
  }


#' @title antsrMeshAddPoint
#' @description add point to mesh
#' @param mesh an 'antsrMesh'
#' @param point spatial point to add to mesh
#' @param identifier index of point to add
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=128 )
#' antsrMeshAddPoint( x, c(0,0,0) )
#' @export
  antsrMeshAddPoint = function( mesh, point, identifier=NA ) {
    if ( is.na(identifier) )  {
      identifier = antsrMeshGetNumberOfPoints(mesh)+1
    }
    invisible(.Call("antsrMesh_AddPoint", mesh, identifier, point, package="DANTsR"))
  }

#' @title antsrMeshSetPoint
#' @description set a given point in mesh
#' @param mesh an 'antsrMesh'
#' @param identifier identified for point to set
#' @param point spatial point to set in mesh
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=128 )
#' antsrMeshSetPoint( x, c(0,0,0), 1 )
#' @export
  antsrMeshSetPoint = function( mesh, point, identifier ) {
    invisible(.Call("antsrMesh_SetPoint", mesh, identifier, point, package="DANTsR"))
  }

#' @title antsrMeshGetPoint
#' @description get a given point in mesh
#' @param mesh an 'antsrMesh'
#' @param identifier identifier of point to get
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=128 )
#' antsrMeshSetPoint( x, c(0,0,0), 0 )
#' pt = antsrMeshGetPoint(x, 0)
#' @export
  antsrMeshGetPoint = function( mesh, identifier ) {
    .Call("antsrMesh_GetPoint", mesh, identifier, package="DANTsR")
  }

#' @title antsrMeshGetCell
#' @description get a given cell in mesh
#' @param mesh an 'antsrMesh'
#' @param identifier identifier of point to get
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=128 )
#' antsrMeshAddPoint( x, c(0,0,0) )
#'
#' @export
  antsrMeshGetCell = function( mesh, identifier ) {
    .Call("antsrMesh_GetCell", mesh, identifier, package="DANTsR")
  }

#' @title antsrMeshGetCellPoints
#' @description get points for a given cell in mesh
#' @param mesh an 'antsrMesh'
#' @param identifier identifier of point to get
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=128 )
#' antsrMeshAddPoint( x, c(0,0,0) )
#'
#' @export
  antsrMeshGetCellPoints = function( mesh, identifier ) {
    .Call("antsrMesh_GetCellPoints", mesh, identifier, package="DANTsR")
  }

#' @title antsrMeshGetPoints
#' @description get all points in mesh
#' @param mesh an 'antsrMesh'
#' @param identifiers identifiers of points to get
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=128 )
#' antsrMeshAddPoint( x, c(0,0,0) )
#' pt = antsrMeshGetPoints(x, 1)
#' @export
antsrMeshGetPoints = function( mesh, identifiers=NULL) {
  if ( is.null(identifiers) ) {
    identifiers = numeric(0)
  }
  .Call("antsrMesh_GetPoints", mesh, identifiers, package="DANTsR")
}

#' @title antsrMeshAddPolyline
#' @description add a polyline cell to the mesh
#' @param mesh an 'antsrMesh'
#' @param points array of point indices
#' @param identifier index of cell to add
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=128 )
#' antsrMeshAddPoint( x, c(0,0,0), 1 )
#' antsrMeshAddPoint( x, c(1,0,0), 2 )
#' antsrMeshAddPoint( x, c(1,1,0), 3 )
#' antsrMeshAddPolyline( x, c(0,1,2), 1)
#' @export
antsrMeshAddPolyline = function( mesh, points, identifier=NA ) {
  if ( is.na(identifier) )  {
    identifier = antsrMeshGetNumberOfCells(mesh)
  }

  invisible(.Call("antsrMesh_AddPolyline", mesh, identifier, points, package="DANTsR"))
}

#' @title antsrMeshAddCell
#' @description add a cell to the mesh
#' @param mesh an 'antsrMesh'
#' @param points array of point indices
#' @param type of cell to add
#' @param identifier index of cell to add
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=128 )
#' antsrMeshAddPoint( x, c(0,0,0), 1 )
#' antsrMeshAddPoint( x, c(1,0,0), 2 )
#' antsrMeshAddPoint( x, c(1,1,0), 3 )
#' antsrMeshAddCell( x, c(0,1,2), "polyline", 1)
#' @export
antsrMeshAddCell = function( mesh, points, type, identifier=NA ) {
  if ( is.na(identifier) )  {
    identifier = antsrMeshGetNumberOfCells(mesh)
  }

  if ( type=="polyline" ) {
    invisible(.Call("antsrMesh_AddPolyline", mesh, identifier, points, package="DANTsR"))
  }
  else {
    stop( "Unsupported cell type")
  }

}

#' @title applyAntsrTransformToMesh
#' @description Apply transform/s to an antsrMesh
#' @param transform antsrTransform
#' @param mesh antsrMesh to transform
#' @param in.place if true modify the input mesh, if false return a new mesh
#' @return antsrMesh
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=1 )
#' antsrMeshAddPoint( x, c(1,2,3) )
#' tx = new("antsrTransform")
#' params = getAntsrTransformParameters(tx)
#' setAntsrTransformParameters(tx, params*2)
#' x2 = applyAntsrTransformToMesh(tx, x)
#' @export
applyAntsrTransformToMesh <- function(transform, mesh, in.place=FALSE) {
  if ( typeof(transform) == "list")
  {
    transform <- composeAntsrTransforms(transform)
  }
  return(.Call("antsrMesh_TransformMesh", transform, mesh, in.place, PACKAGE = "DANTsR"))
  return(NA)
}

#' @title read.antsrMesh
#' @description read a file into an antsrMesh
#' @param filename name of the file to read
#' @param dimension dimension of point data in mesh
#' @param pixeltype float or double
#' @return antsrMesh
#' @export
read.antsrMesh = function( filename, dimension=3, pixeltype="float" ) {
  mesh = NA
  if ( grepl(".vtk", filename ) ) {
    mesh = .Call("antsrMesh_ReadVTK", filename, dimension, pixeltype, package="DANTsR")
  }
  else if ( grepl(".Bfloat", filename ) ) {
    mesh = .Call("antsrMesh_ReadCamino", filename, pixeltype="float", package="DANTsR")
  }
  return(mesh)
}

#' @title write.antsrMesh
#' @description write an antsrMesh to a file
#' @param mesh antsrMesh to write
#' @param filename name of the file to read
#' @param seeds seed indices (for Camino files)
#' @param cell.as what VTK-type of cell should cell data be written as ("NA", "polygon", "line")
#' @export
write.antsrMesh = function( mesh, filename, seeds=NULL, cells.as="NA" ) {
  if ( grepl(".vtk", filename ) ) {
    #print("Writing VTK mesh")
    invisible(.Call("antsrMesh_WriteVTK", mesh, filename, cells.as, package="DANTsR"))
  }
  else if ( grepl(".Bfloat", filename ) ) {
    if (is.null(seeds) ) {
      stop("Camino file needs seeds")
    }
    else {
      #print(typeof(mesh))
      #print(typeof(filename))
      #print(typeof(seeds))
      invisible(.Call("antsrMesh_WriteCamino", mesh, filename, seeds, package="DANTsR"))
    }
  }
  return(0)
}
