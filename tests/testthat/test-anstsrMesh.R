context("antsrMesh")

test_that("mesh can be created with new()", {
  x = new("antsrMesh")
  expect_true( class(x)=="antsrMesh")
})

test_that("mesh can be created with new()", {
  x = new("antsrMesh", reserve=1)
  expect_true( ( class(x)=="antsrMesh" ) )
})

# Test all combos for the C++ template parameters
for ( d in c(2,3,4) ) {
  for ( p in c("float", "double") ) {

    test_that("mesh can be created with new()", {
      x = new("antsrMesh", dimension=d, precision=p)
      expect_true( ( class(x)=="antsrMesh" ) &
                   ( x@dimension == d ) &
                   ( x@precision == p ) )
    })
  }
}

test_that("mesh can be created with antsrMeshCreate()", {
  x = antsrMeshCreate()
  expect_true( class(x)=="antsrMesh")
})

test_that("mesh can be created with antsrMeshCreate() using reserve", {
  x = antsrMeshCreate(reserve=1)
  expect_true( ( class(x)=="antsrMesh" )  )
})

test_that("mesh can be created with antsrMeshCreate() using point matrix", {
  pts = matrix(0, 10, 3)
  x = antsrMeshCreate(points=pts)
  expect_true( ( class(x)=="antsrMesh" )  )
})

test_that("mesh reports correct number of points", {
  x = new("antsrMesh", reserve=10)
  expect_true( antsrMeshGetNumberOfPoints(x)==10 )
})

test_that("mesh reports correct number of cells", {
  x = new("antsrMesh", reserve=10)
  expect_true( antsrMeshGetNumberOfCells(x)==0 )
})

# Test all combos for the C++ template parameters
for ( d in c(2,3,4) ) {
  for ( p in c("float", "double") ) {

    test_that("mesh can be created with antsrMeshCreate()", {
      x = antsrMeshCreate(dimension=d, precision=p)
      expect_true( ( class(x)=="antsrMesh" ) &
                   ( x@dimension == d ) &
                   ( x@precision == p ) )
    })
  }
}

test_that( "point can be added to a mesh with no index", {
  x = antsrMeshCreate()
  antsrMeshAddPoint(x, c(1,2,3) )
  antsrMeshAddPoint(x, c(1,2,3) )
  expect_true( antsrMeshGetNumberOfPoints(x) == 2)
})

test_that( "point can be retrieved from mesh", {
  x = antsrMeshCreate()
  antsrMeshAddPoint(x, c(1,2,3) )
  expect_true( sum(antsrMeshGetPoint(x,1) == c(1,2,3))==3 )
})

test_that( "point can be added to a mesh with an index", {
  x = antsrMeshCreate()
  antsrMeshAddPoint(x, c(1,2,3), 5 )
  expect_true( antsrMeshGetNumberOfPoints(x) == 5)
})

test_that( "mesh points can be retrieved as matrix", {
  x = antsrMeshCreate()
  for (i in 1:4) {
    antsrMeshAddPoint(x, i*c(1,2,3), i )
  }
  pts = antsrMeshGetPoints(x)
  expect_true( ( class(pts) == "matrix" ) &
               ( sum( dim(pts) == c(4,3) ) == 2 ) &
               ( sum(pts) == 60 ) )
})

test_that("mesh points can be transformed", {
  x =  antsrMeshCreate( 3, "float" )
  antsrMeshAddPoint( x, c(1,2,3) )
  tx = new("antsrTransform")
  params = getAntsrTransformParameters(tx)
  setAntsrTransformParameters(tx, params*2)
  x2 = applyAntsrTransformToMesh(tx, x)
  pts2  = antsrMeshGetPoints(x2)
  expect_true( sum(pts2)==12 )
})

test_that("mesh points can be transformed in place", {
  x =  antsrMeshCreate( 3, "float" )
  antsrMeshAddPoint( x, c(1,2,3) )
  tx = new("antsrTransform")
  params = getAntsrTransformParameters(tx)
  setAntsrTransformParameters(tx, params*2)
  invisible(applyAntsrTransformToMesh(tx, x, in.place=TRUE))
  pts2  = antsrMeshGetPoints(x)
  expect_true( sum(pts2)==12 )
})

test_that( "polyline can be added to mesh directly", {
  x =  antsrMeshCreate( 3, "float", reserve=3 )
  antsrMeshAddPoint( x, c(0,0,0), 1 )
  antsrMeshAddPoint( x, c(1,0,0), 2 )
  antsrMeshAddPoint( x, c(1,1,0), 3 )
  antsrMeshAddPolyline( x, c(0,1,2), 1)
  expect_true( antsrMeshGetNumberOfCells(x)==1 )
})

test_that( "polyline can be added to mesh by name", {
  x =  antsrMeshCreate( 3, "float", reserve=3 )
  antsrMeshAddPoint( x, c(0,0,0), 1 )
  antsrMeshAddPoint( x, c(1,0,0), 2 )
  antsrMeshAddPoint( x, c(1,1,0), 3 )
  antsrMeshAddCell( x, c(0,1,2), "polyline", 1)
  expect_true( antsrMeshGetNumberOfCells(x)==1 )
})

test_that( "cell can be retrieved from mesh", {
  x =  antsrMeshCreate( 3, "float", reserve=3 )
  antsrMeshAddPoint( x, c(0,0,0), 1 )
  antsrMeshAddPoint( x, c(1,0,0), 2 )
  antsrMeshAddPoint( x, c(1,1,0), 3 )
  antsrMeshAddPolyline( x, c(0,1,2), 1)
  c = antsrMeshGetCell(x, 1)
  expect_true( length(c)==3 )
})

test_that( "cell points can be retrieved from mesh", {
  x =  antsrMeshCreate( 3, "float", reserve=3 )
  antsrMeshAddPoint( x, c(0,0,0), 1 )
  antsrMeshAddPoint( x, c(1,0,0), 2 )
  antsrMeshAddPoint( x, c(1,1,0), 3 )
  antsrMeshAddPolyline( x, c(1,2,3), 1)
  pts = antsrMeshGetCellPoints(x, 1)
  expect_true( sum(pts)==3 )
})

test_that( "read vtk polydata", {
  tmp = tempfile()
  tmp = paste(sep="", tmp, ".vtk")
  cat( getExamplePolydata(), file=tmp)
  x = read.antsrMesh(tmp)
  expect_true( class(x$Mesh)=="antsrMesh" )
})

test_that( "write and read vtk polydata binary", {
  tmp = tempfile()
  tmp = paste(sep="", tmp, ".vtk")
  x = rnorm(300)
  dim(x) = c(100,3)
  m = antsrMeshCreate(points=x)
  for ( i in 1:10 ) {
    antsrMeshAddPolyline( m, (1:10)+(i-1), i)
  }

  write.antsrMesh( m, tmp, cells.as="lines" )
  m2 = read.antsrMesh(tmp)

  expect_equal( 0, sum( x-antsrMeshGetPoints( m2$Mesh) ), tolerance=0.00001 )
})

test_that( "write and read vtk polydata ascii (greater tolerance)", {
  tmp = tempfile()
  tmp = paste(sep="", tmp, ".vtk")
  x = rnorm(300)
  dim(x) = c(100,3)
  m = antsrMeshCreate(points=x)
  for ( i in 1:10 ) {
    antsrMeshAddPolyline( m, (1:10)+(i-1), i)
  }

  write.antsrMesh( m, tmp, cells.as="lines", binary=FALSE )
  m2 = read.antsrMesh(tmp)

  expect_equal( 0, sum( x-antsrMeshGetPoints( m2$Mesh) ), tolerance=0.0001 )
})

test_that( "write and read camino file of floats", {
  tmp = tempfile()
  tmp = paste(sep="", tmp, ".Bfloat")
  x = rnorm(300)
  dim(x) = c(100,3)
  m = antsrMeshCreate(points=x)
  for ( i in 1:10 ) {
    antsrMeshAddPolyline( m, (1:10)+10*(i-1), i)
  }
  seeds = 1:10

  write.antsrMesh( m, tmp, seeds=seeds )
  m2 = read.antsrMesh(tmp)

  expect_equal( 0, sum( x-antsrMeshGetPoints( m2$Mesh) ), tolerance=0.00001 )
})

test_that( "write and read camino file of doubles", {
  tmp = tempfile()
  tmp = "file"
  tmp = paste(sep="", tmp, ".Bdouble")
  x = rnorm(300)
  dim(x) = c(100,3)
  m = antsrMeshCreate(points=x)
  for ( i in 1:10 ) {
    antsrMeshAddPolyline( m, (1:10)+10*(i-1), i)
  }
  seeds = 1:10

  write.antsrMesh( m, tmp, seeds=seeds )
  m2 = read.antsrMesh(tmp)

  expect_equal( 0, sum( x-antsrMeshGetPoints( m2$Mesh) ), tolerance=0.00001 )
})

test_that( "write and read TrackVis .trk file", {
  tmp = tempfile()
  tmp = paste(sep="", tmp, ".trk")
  x = rnorm(300)
  dim(x) = c(100,3)
  m = antsrMeshCreate(points=x)
  for ( i in 1:10 ) {
    antsrMeshAddPolyline( m, (1:10)+10*(i-1), i)
  }
  seeds = 1:10

  img = makeImage( c(10,10,10) )
  write.antsrMesh( m, tmp, seeds=seeds, image=img )
  m2 = read.antsrMesh(tmp)

  expect_equal( 0, sum( x-antsrMeshGetPoints( m2$Mesh) ), tolerance=0.00001 )
})

test_that( "write and read MRTrix .tck file", {
  tmp = tempfile()
  tmp = paste(sep="", tmp, ".tck")
  x = rnorm(300)
  dim(x) = c(100,3)
  m = antsrMeshCreate(points=x)
  for ( i in 1:10 ) {
    antsrMeshAddPolyline( m, (1:10)+10*(i-1), i)
  }
  seeds = 1:10

  img = makeImage( c(10,10,10) )
  write.antsrMesh( m, tmp, seeds=seeds, image=img )
  m2 = read.antsrMesh(tmp)

  expect_equal( 0, sum( x-antsrMeshGetPoints( m2$Mesh) ), tolerance=0.00001 )
})

#test_that("matrix of points can be passed to mesh", {
#  x = antsrMeshCreate()
#  m = matrix(1:12, 4, 3)
#  antsrMeshSetPoints(x, m)
#}
