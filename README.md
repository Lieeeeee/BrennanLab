Brennan Lab ImageJ Plugins
==========================
These are ImageJ plugins developed for the Brennan Lab's use in analyzing optical intrinisc imaging experiments of spreading depression in mouse cortex. To install, download the jar file from releases and place it into plugins/ in your ImageJ directory. The jar file will place an entry under Plugins/Brennan Lab.

External dependencies
---------------------
* Cern Colt 1.2.0
* UJMP 0.3.0
* JChart2D 3.3.2
* Apache commons math 3.2
* ImageJ >=1.49

Plugins
-------
* Bright CSD Tracker from differences - CSD wavefront segmentation based on inter-frame difference image
* Utils
  * Differences - compute interframe differences
  * Filtered Differences - compute interframe differences with some saturation of outliers
  * Interactive Z Axis Profile - Z Axis profile sensitive to changes in the ROI selection
  * Open a substack of a Tiff - for selecting a particular range (and loading into memory) of a multipage tiff
  * Open multiple subsequences from a Tiff
  * Changes from baseline - Various forms of baseline subtraction
  * Tresholded changes from baseline - Take the changes, smooth, threshold
  * Blue and threshold versus baseline
  * Sign of differences - Sign of interframe differences
  * Cumulative sum - pixelwise cumulative sum function
* CSD
  * Peak Isolator - Isolate pixelwise peak in the image intensity
  * Montage CSDs from file - Make a montage based on the Peak Isolator
  * Levelset Method
    * Reconstruct Wave from Speed - Soft eikonal equation with a given intial condition
    * Signed distance to ROI - Get 0-level set

CSD Tracker
===========

![Segmented moving CSD wavefront](http://joshchang.github.io/_static/images/csd.gif)

The CSD tracker is based on regularized Chan-Vese segmentation (of a bright expanding foreground region) where inference is performed using graph cuts. The regularization comes from a running prediction of future front positions based on previous front positions by estimating a spatially-correlated wavespeed field. The paper in IEEE describes the overall method (and is completely reproducible as I fortunately discovered recently when rewriting this plugin!):

```
@article{chang2012tracking,
  title={Tracking monotonically advancing boundaries in image sequences using graph cuts and recursive kernel shape priors},
  author={Chang, Joshua C and Brennan, Kevin C and Chou, Tom},
  journal={Medical Imaging, IEEE Transactions on},
  volume={31},
  number={5},
  pages={1008--1020},
  year={2012},
  publisher={IEEE}
}
```

The shape priors are also normalized now using the method from my paper in JMIV

```
@article{chang2014iterative,
  title={Iterative graph cuts for image segmentation with a nonlinear statistical shape prior},
  author={Chang, Joshua C and Chou, Tom},
  journal={Journal of mathematical imaging and vision},
  volume={49},
  number={1},
  pages={87--97},
  year={2014},
  publisher={Springer}
}
```

If you find this plugin useful, please cite these papers in your work. Thanks!
