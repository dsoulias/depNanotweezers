# Imaging parameters

pixelSize = 6.5 # camera default pixel size in um (Prime=6.5, Emc2=8, ORCA-ER=6.45)
magObjective = 60 # magnification of objective lens used during recording
magLens = 1 # lens between the camera and the microscope objective (none = 1)
cMount = 1 # C-mount lens between camera and microscope's optical bath (none = 1)

# Input the following:

binning = int(input('Binning: ')) # enter binning mode set during recording (file metadata)
pxImgWidth = int(input('Image width in pixels: ')) # enter image width in pixels
pxImgHeight = int(input('Image height in pixels: ')) # enter image height in pixels

# Equation for converging pixels into um

factor = pixelSize * (binning / (magObjective * magLens * cMount))
print("The factor conveging pixel to \u03BCm is:", round(factor,3))

# Image widht and height in um

umImgWidth = pxImgWidth * pixelSize * (binning / (magObjective * magLens * cMount))
print("Image width in \u03BCm =", round(umImgWidth,2))

umImgHeight = pxImgHeight * pixelSize * (binning / (magObjective * magLens * cMount))
print("Image height in \u03BCm =", round(umImgHeight,2))