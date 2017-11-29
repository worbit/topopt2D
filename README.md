# topopt2D
## Topology optimization using images to control load, support and passive elements

Adaptation of python code taken from here: http://www.topopt.dtu.dk/?q=node/881.
The number of elements is not hardcoded but taken from the size of the images. The images should be the same size, the code does not test for that though.
In line 14 and 15 of the main python file (`topopt2D.py`), the name of the folder containing the load, support (and if necessary passive) maps can be specified.

### Loads
Loads and supports can now easily be drawn in any image processing software. The red channel defines vertical loads, the green channel horizontal loads, the blue channel is not taken into account. Color values are scaled to range from -1 to +1 `vl=(R-128)/128.0`, so all gray (128,128,128) represents no load. The folder `beam` contains an example with a continuous vertical load decreasing in intensity from left to right.

![load direction](https://github.com/worbit/topopt2D/blob/master/map.png)

A vertical load (down in y-direction) is an orange pixel `(255,128,128)`, a load from left to right is greenish `(128,255,128)`, one to the top left is dark blue `(0,0,128)` etc. Below is a low resolution sample for the half MBB beam example:

![load example](https://github.com/worbit/topopt2D/blob/master/load.png)

### Support
For y-support, the red channel must be 255, for x-support the green channel must be 255. For both x and y, the pixel must therefore be yellow. Below is a low resolution sample for the half MBB beam example:

![load example](https://github.com/worbit/topopt2D/blob/master/support.png)

### Passive
Nodes that need to be void must be red, solid ones must be green. The blue channel doesn't have any effect whatsoever.
