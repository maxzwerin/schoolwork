# **SCHOOLWORK**
_an archive of some cool (and hopefully useful) work I have done_

## how to run programs (on linux)
compile with: _cc FILENAME.c_\
add _-lm_ to link with math library\
add _-lX11_ to link with X11 graphical library

this will create an executable _a.out_\
run the program with _./a.out_

in summary, the terminal should look something like this:
```
cc FILENAME.c -lm -lX11\
./a.out
```

## **layout**
### TOOLS
- FPToolkit : a simple set of graphical tools
- matrix_tools : some useful matrix manipulation tools for 2D and 3D usecases

### GRAPHICS
- all graphics labs (poorly documented)
- XYZ files

NOTE: most of the labs require XYZ file(s) to be passed in when running the executable:
```
./a.out XYZ/torus.xyz XYZ/sphere.xyz
```

### NUMERICAL
- all numerical analysis labs
- test files for labs 01 and 02

## NOTE
most of these files simply call _#include "FPToolkit.c"_\
This means the path will need to be modified accordingly to\
_#include "../TOOLS/FPToolkit.c"_ for each file.

ALTERNATIVELY you could copy/paste the contents of _/TOOLS/_\
into each folder (which is much easier). enjoy!
