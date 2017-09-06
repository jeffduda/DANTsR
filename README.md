![dantsr killers](http://www.picslyrics.net/images/255741-the-killers-are-we-human-or-are-we-dancer.jpg)

# DANTsR
Diffusion imaging extension to ANTsR.

img = antsImageRead("yourfile.nii.gz")

fa = dtiAnisotrpyImage(img, "FractionalAnisotropy")

ra = dtiAnisotropyImage(img, "RelativeAnisotropy")

eig = dtiEigenSystem(img)

dec = dtiColorMap(img)
