# DANTsR
Diffusion imaging extension to ANTsR.

img = antsImageRead("yourfile.nii.gz")

fa = dtiAnisotrpyImage(img, "FractionalAnisotropy")

ra = dtiAnisotropyImage(img, "RelativeAnisotropy")

eig = dtiEigenSystem(img)

dec = dtiColorMap(img)
