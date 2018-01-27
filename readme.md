# False discovery rate control of cryo-EM maps

Confidence maps are complementary maps generated from cryo-EM maps by means of statistical hypothesis testing and subsequent FDR control. They allow thresholding of EM maps based on the expected amount of background noise visible at the respective threshold and thus allow rigorous error assesment of visible features in the density. 
Additional postprocessing like local filtering or local amplitude scaling (LocScale) can be incorporated in the framework in order to increase the power.


## Getting Started


### Prerequisites

The software is written in Python and itself dependent on EMAN2, NumPy and SciPy libraries. 

For incorporation of the local ampltidue scaling (LocScale) procedue, additional LocScale libraries are needed, as described here:  https://git.embl.de/jakobi/LocScale/wikis/home


### Installing

The software consists of four scripts, FDRcontrol.py, mapUtil.py, FDRUtil.py and locscaleUtil.py, that just have to be copied to your computer.

Once you have a working EMAN2 installation ( http://blake.bcm.edu/emanwiki/EMAN2 ), the software can be simply run by using your Python version that comes with EMAN2:   

For example if /programs/x86_64-linux/eman2/2.2/bin/python is the respective EMAN2 python installation: 
```
/programs/x86_64-linux/eman2/2.2/bin/python FDRcontrol.py -em yourMap.mrc -p thePixelSize
```

If you want to use the LocScale feature, you should have a running LocScale version, together with MPI,  as described here:  https://git.embl.de/jakobi/LocScale/wikis/home.  



## How to use

The simplest, and probably most important case, is the generation of a confidence map of a simple cryoEM density without incorporation of local resolution or atomic model information.

```
/programs/x86_64-linux/eman2/2.2/bin/python FDRcontrol.py -em yourMap.mrc -p thePixelSize
```

The output will be the corresponding confidence map ( yourMap_confidenceMap.mrc ) and a diagnostic image (diag_image.pdf), that shows three slices thorugh the map together with the regions used for noise estimation. In order to get good estimates of the background noise distribution, you should make sure that the region contains just noise and no particle.

Size of the noise estimation region can be adjusted with -w sizeOfRegion. If the default regions for noise estimation fall into noise, you can specify the center of region of your choice with -noiseBox x y z .
For example:

```
/programs/x86_64-linux/eman2/2.2/bin/python FDRcontrol.py -em yourMap.mrc -p 5.6 -w 20 -boxNoise 50 50 120
```

### Incorporation of local resolution information

A corresponding local resolution map can be supplied to the program by -locResMap yourLocalResolutionMap.mrc .

The ouput will the be a locally filtered map together with the diagnostic image and the confidence map. Adjustment of the noise estimation region can be done accordingly.

For example:
```
/programs/x86_64-linux/eman2/2.2/bin/python FDRcontrol.py -em yourMap.mrc -p thePixelSize -locResMap yourLocalResolutionMap
```

### Incorporation of local amplitude scaling

If 
