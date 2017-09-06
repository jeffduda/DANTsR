# DANTsR
Diffusion imaging for ANTsR.

img = antsImageRead("yourfile.nii.gz")
fa = dtiAnisotrpyImage(img, "FractionalAnisotropy")
ra = dtiAnisotropyImage(img, "RelativeAnisotropy")
eig = dtiEigenSystem(img)
dec = dtiColorMap(img)
