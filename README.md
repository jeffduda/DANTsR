![dantsr killers](http://www.picslyrics.net/images/255741-the-killers-are-we-human-or-are-we-dancer.jpg)

 [![Travis-CI Build Status] (https://travis-ci.org/jeffduda/DANTsR.svg?branch=master)] (https://travis-ci.org/jeffduda/DANTsR)

# DANTsR
Diffusion imaging extension to ANTsR.

img = antsImageRead("yourfile.nii.gz")

fa = dtiAnisotrpyImage(img, "FractionalAnisotropy")

ra = dtiAnisotropyImage(img, "RelativeAnisotropy")

eig = dtiEigenSystem(img)

dec = dtiColorMap(img)

seeds = labelsToPoints( fa > 0.1 )

fibers = deterministicTractograpy( v1, seeds )
