In this tutorial we show you how to generate STM images under the
Tersoff-Hamann approximation, and visualise them in VESTA.  Here 
we use a simple model of the Si(001) surface.

The first step is to do a ground state Conquest calculation using
the files provided.  This will write out the wavefunction files
(Process000000?WF.dat) needed to simulate STM images.  The only
thing to note is that we have specified that wavefunctions between
-0.15Ha and 0.15Ha (relative to the Fermi level) are written out, 
which limits the range of bias voltages to between -4V and 4V.

The post-processing calculation generates two Cube files (which can
be read by VESTA).  STM.cube contains the current values on a 3D
grid while STMHeight.cube gives the height of each grid point 
(used for colouring the isosurface).

Opening STM.cube should give an image similar to 1JustSTM.png which
shows an isosurface of current.  While this is correct, it is not 
close to standard STM images, so we need to make further adjustments.
First, we colour the isosurface by height.  In VESTA, this is done
using Edit->Edit Data->Volumetric Data and then under Surface Coloring
selecting the STMHeight.cube file.  The image should now look like
2AddHeightSTM.png.

We now adjust the isosurface values using Properties->Isosurfaces.
Set Opacity 1 to 255 (from 127) and the isosurface value to 0.0008
(in this case - it will require different values for different biases
and different systems).  This should give an image like 3IsosurfaceSTM.png.

Finally we adjust the colouring, using Properties->Sections and
selecting Grayscale and then Properties->Isosurfaces and adjusting
the Saturation Level under Surface Coloring, setting the Minimum
level to two or three units less than the Maximum level (in this case
I chose 15).  The image should now resemble 4FinalSTM.png.
