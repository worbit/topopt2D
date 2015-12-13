# topopt2D
##Topology optimization using images to control load, support and passive elements

Adaptation of python code taken from here: http://www.topopt.dtu.dk/?q=node/881.

###Loads
Loads and supports can now easily be drawn in any image processing software. The number of elements is not hardcoded but taken from the size of the images.

![load direction](https://github.com/worbit/topopt2D/blob/master/map.png)

A vertical load (down in y-direction) is an orange pixel, a load from left to right is greenish etc. Below is a low resolution sample for the half MBB beam example:

![load example](https://github.com/worbit/topopt2D/blob/master/load.png)

###Support
For y-support, the red channel must be 255, for x-support the green channel must be 255. For both x and y, the pixel must therefore be yellow. Below is a low resolution sample for the half MBB beam example:

![load example](https://github.com/worbit/topopt2D/blob/master/support.png)

###Passive
Nodes that need to be void must be red, solid ones must be green.
