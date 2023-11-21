# gencolormap

this is a forge of gencolormap program by Martin Lambers 
* [gencolormap home page](https://marlam.de/gencolormap/)
* [flathub](https://flathub.org/pl/apps/de.marlam.gencolormap)

These tools generate color maps for visualization.
A variety of methods for sequential, diverging, and qualitative maps is available.

The color map generation is contained in just two C++ files (`colormap.hpp` and
`colormap.cpp`) and requires no additional libraries. You can simply copy these
two files to your own project.

Two frontends are included: 
* a GUI for interactive use, the GUI requires Qt. 
* and a command line tool for scripting. The command line tool requires no libraries, 

The following papers are implemented:

- M. Lambers.
  [Interactive Creation of Perceptually Uniform Color Maps](https://dx.doi.org/10.2312/evs.20201048).
  Proc. EuroVis Short Papers, May 2020.
- M. Wijffelaars, R. Vliegen, J.J. van Wijk, E.-J. van der Linden.
  [Generating color palettes using intuitive parameters](https://dx.doi.org/10.1111/j.1467-8659.2008.01203.x).
  Computer Graphics Forum 27(3), May 2008.
- K. Moreland.
  [Diverging Color Maps for Scientific Visualization](https://dx.doi.org/10.1007/978-3-642-10520-3_9).
  Proc. Int. Symp. Visual Computing, December 2009.
- D. A. Green.
  [A colour scheme for the display of astronomical intensity](https://ui.adsabs.harvard.edu/abs/2011BASI...39..289G/abstract).
  Bulletin of the Astronomical Society of India 39(2), June 2011.
- J. McNames.
  [An Effective Color Scale for Simultaneous Color and Gray-Scale Publications](https://dx.doi.org/10.1109/MSP.2006.1593340).
  IEEE Signal Processing Magazine 23(1), January 2006.

![GUI screen shot](https://marlam.de/gencolormap/gencolormap-screenshot.png)



Notes about the color spaces used internally:
* We use D65 white everywhere
* RGB means linear RGB; we also have sRGB
* RGB and sRGB values are in [0,1]
* XYZ, LUV, and similar values are in the original range (not normalized);
* often this is [0,100]
* All angles (for hue) are measured in radian
 
 
 
Generate color maps for scientific visualization purposes.

Usage:
* Decide which type of color map you need and how many colors the map should contain.
* Allocate memory for you color map (3 * unsigned char for each color entry).
* Call the function that generates your color map.
* the return value is always the number of colors that had to be clipped to fit into sRGB; you want to keep that number low by adjusting parameters.
* All colors are represented as unsigned char sRGB triplets, with each value in [0,255] range



# src code

I have converted library and CLI program to [the one file program](./src/cmdline/g.cpp)






# compile

for all original versions

```
sudo apt-get install qt6-base-dev
cd gencolormap
mkdir build
cd build
cmake ..
make
```

My one file program

```
g++ g.cpp -Wall -Wextra -lm
```


# usage

./a.out --help
./a.out -H


```txt
Usage: ./a.out [option...]
Generates a color map and prints it to standard output.
Prints the number of colors that had to be clipped to standard error.
Common options:
  [-f|--format=csv|json|ppm]          Set output format
  [-n|--n=N]                          Set number of colors in the map
Brewer-like color maps:
  [-t|--type=brewer-sequential]       Generate a sequential color map
  [-t|--type=brewer-diverging]        Generate a diverging color map
  [-t|--type=brewer-qualitative]      Generate a qualitative color map
  [-h|--hue=H]                        Set default hue in [0,360] degrees
  [-c|--contrast=C]                   Set contrast in [0,1]
  [-s|--saturation=S]                 Set saturation in [0,1]
  [-b|--brightness=B]                 Set brightness in [0,1]
  [-w|--warmth=W]                     Set warmth in [0,1] for seq. and div. maps
  [-d|--divergence=D]                 Set diverg. in deg for div. and qual. maps
Perceptually uniform color maps:
  [-t|--type=pusequential-lightness]  Sequential map, varying lightness
  [-t|--type=pusequential-saturation] Sequential map, varying saturation
  [-t|--type=pusequential-rainbow]    Sequential map, varying hue (rainbow)
  [-t|--type=pusequential-blackbody]  Sequential map, varying hue (black body)
  [-t|--type=pusequential-multihue]   Sequential map, varying hue (custom)
  [-t|--type=pudiverging-lightness]   Diverging map, varying lightness
  [-t|--type=pudiverging-saturation]  Diverging map, varying saturation
  [-t|--type=puqualitative-hue]       Qualitative map, evenly distributed hue
  [-l|--lightness=L]                  Set lightness in [0,1]
  [-L|--lightness-range=LR]           Set lightness range in [0.7,1]
  [-s|--saturation=S]                 Set saturation in [0,1]
  [-S|--saturation-range=SR]          Set saturation range in [0.7,1]
  [-h|--hue=H]                        Set default hue in [0,360] degrees
  [-d|--divergence=D]                 Set diverg. in deg for div. and qual. maps
  [-r|--rotations=R]                  Set number of rotations for rainbow maps
  [-T|--temperature=T]                Set start temp. in K for black body maps
  [-R|--temperature-range=TR]         Set range for temperature in K
  [-V|--hue-values=H0,H1,...]         Set hue values in [0,360] for multi-hue maps
  [-P|--hue-positions=P0,P1,...]      Set hue positions in [0,1] for multi-hue maps
CubeHelix color maps:
  [-t|--type=cubehelix]               Generate a CubeHelix color map
  [-h|--hue=H]                        Set start hue in [0,180] degrees
  [-r|--rotations=R]                  Set number of rotations, in (-infty,infty)
  [-s|--saturation=S]                 Set saturation, in [0,1]
  [-g|--gamma=G]                      Set gamma correction, in (0,infty)
Moreland diverging color maps:
  [-t|--type=moreland]                Generate a Moreland diverging color map
  [-A|--color0=sr,sg,sb]              Set the first color as sRGB in [0,255]
  [-O|--color1=sr,sg,sb]              Set the last color as sRGB in [0,255]
McNames sequential color maps:
  [-t|--type=mcnames]                 Generate a McNames sequential color map
  [-p|--periods=P]                    Set the number of periods in (0, infty)
```  
  
Defaults: format=csv, n=256, type=brewer-sequential



https://marlam.de/gencolormap


# version


./a.out -v

```txt
gencolormap version 2.3
https://marlam.de/gencolormap
Copyright (C) 2022 Computer Graphics Group, University of Siegen.
Written by Martin Lambers <martin.lambers@uni-siegen.de>.
This is free software under the terms of the MIT/Expat License.
There is NO WARRANTY, to the extent permitted by law.
```
# Similar projects
* [1D-RGB-color-gradient](https://github.com/adammaj1/1D-RGB-color-gradient)
* [1D RGB gray gradients or wave forms](https://github.com/adammaj1/Waveform)
* [hsluv-color-gradient](https://github.com/adammaj1/hsluv-color-gradient)
* [1D-pastel-color-gradient-HSV](https://github.com/adammaj1/1D-pastel-color-gradient-HSV)
* [2D-color-gradient-or-Procedural-texture](https://github.com/adammaj1/2D-color-gradient-or-Procedural-texture)
* [Mandelbrot-set-with-blended-gradients](https://github.com/adammaj1/Mandelbrot-set-with-blended-gradients)


# Git

Git original repos ( the web frontend)
* [github](https://github.com/marlam/gencolormap-mirror)
* [git](https://git.marlam.de/gitweb/?p=gencolormap.git)


```
git clone https://git.marlam.de/git/gencolormap.git
```


```
git clone git@github.com:marlam/gencolormap-mirror.git
```

Init new repo and first commit

```
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin git@github.com:adammaj1/gencolormap.git
git push -u origin main
```



