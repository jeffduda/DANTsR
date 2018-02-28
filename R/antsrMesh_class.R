# this file defines the S4 classes related to 'antsrTransform' and the associated
# methods



#' @rdname antsrMesh
#' @title antsrMesh
#'
#' @description class for point sets and meshes
#'
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
setMethod(f = "initialize", signature(.Object = "antsrMesh"), definition = function(.Object,
  dimension = 3, precision = "float", reserve=0) {
  mesh = .Call("antsrMesh", precision, dimension, reserve, PACKAGE = "DANTsR")
  return( mesh )
})


#' @title antsrMeshCreate
#' @description create a mesh
#' @param dimension number of dimensions
#' @param percision use 'float' or 'double' for values
#' @param reserve number of points to allocate on creation
#' @return antsrMesh
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=128 )
#' @export
antsrMeshCreate <- function(dimension=3, precision="float", reserve=0)
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

  mesh = .Call("antsrMesh", precision, dimension, reserve, PACKAGE = "DANTsR")

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
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=128 )
#' antsrMeshAddPoint( x, c(0,0,0) )
#' @export
  antsrMeshAddPoint = function( mesh, point, identifier=NA ) {
    if ( is.na(identifier) )  {
      identifer = antsrMeshGetNumberOfPoints(mesh)
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
#' antsrMeshSetPoint( x, c(0,0,0), 0 )
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
#' @param identifier identifier of point to get
#' @examples
#' x =  antsrMeshCreate( 3, "float", reserve=128 )
#' antsrMeshAddPoint( x, c(0,0,0) )
#' pt = antsrMeshGetPoints(x, 0)
#' @export
  antsrMeshGetPoints = function( mesh, identifiers=NULL) {
    if ( is.null(identifiers) ) {
      identifiers = numeric(0)
    }
    .Call("antsrMesh_GetPoints", mesh, identifiers, package="DANTsR")
  }

#' @title applyAntsrTransformToMesh
#' @description Apply transform/s to an antsrMesh
#' @param transform antsrTransform
#' @param image antsrMesh to transform
#' @param interpolation type of interpolator to use
#' @return antsImage
#' @export
applyAntsrTransformToMesh <- function(transform, image, reference, interpolation="linear") {
  if ( typeof(transform) == "list")
  {
    transform <- composeAntsrTransforms(transform)
  }
  #return(.Call("antsrTransform_TransformImage", transform, image, reference, tolower(interpolation), PACKAGE = "ANTsRCore"))
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
#' @export
write.antsrMesh = function( mesh, filename, seeds=NULL ) {
  if ( grepl(".vtk", filename ) ) {
    #mesh = .Call("antsrMesh_ReadVTK", filename, dimension, pixeltype, package="DANTsR")
  }
  else if ( grepl(".Bfloat", filename ) ) {
    if (is.null(seeds) ) {
      stop("Camino file needs seeds")
    }
    else {
      print(typeof(mesh))
      print(typeof(filename))
      print(typeof(seeds))
      .Call("antsrMesh_WriteCamino", mesh, filename, seeds, package="DANTsR")
    }
  }
  return(0)
}
