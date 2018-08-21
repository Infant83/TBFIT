  POSCAR-TB-ribbon
  1 1 1 .5 .0 .5                                # nx,ny,nz
  1  Bi   1.61   0.2 0.5 0.5 1.6           # atom number,atom species,bond length,color
  0.03           0.5 0.5 0.5            # bond width,color
  0                                     # bounding box, 0: off, 1: x, 2:y, 3:z, 4: all 
  15 30                                 # bounding box, for 1~3: single line "min" to "max"
  -2 4                                  # bounding box, for 4 : 3 lines mix_i max_i (i=x,y,z)
  -2 4

################################################################################
Created,modyfied by Hyun-Jung Kim. Dept of Physics. Hanyang Univ. 2008.1.21 Mon. 
################################################################################
  2  H    0.50   0.7 0.0 0.7
1st line contains input file name. 
You can change into what you want to view. 
There may be three types you can use as input file, POSCAR CONTCAR CHGCAR.

After read the input file name, the next three parameters would be read.
It means unit cell periodicity. nx,ny,nz.

3rd line and 4th line contains atom species, bond length, atomic color.

last line means bond properties. bond width, color.

You have to change atom name as the POSCAR file contained. 
If there are three different atoms, then you have to specify by noting
into each lines. For example, let's assume there exists pyridine molecule.
pyridine's chemical formular : C5H5N 
The POSCAR file, maybe you would have wrote it through following order
------------------------------------------------------------------------
pyridine          ==> system 
some real_number  ==> scale factor
a b c
d e f             ==> unit vector matrix
g h j        
5 5 1             ==> five carbon, five hydrogen, one nitrogen
selective dynamics==> If you assumed constaraint into some atoms...
Direct            ==> way to describe atomic position
a a a T T T
b b b F T T 
.......
......            ==> atomic positions and constraint
....             
...
..
--------------------------------------------------------------------------
Look sixth line, ther are three integers. 
The leading term means Carbon, 2nd term : Hydrogen, last term : Nitrogen

And you must write the three lines  according to its order, like following sentence.

 5   C   1.0    0.5 0.5 0.5
 2   H   0.50   0.7 0.0 0.7
 1   N   1.1    0.7 0.1 0.3

___________________________________________________________________________

The OUTPUT file name : poscar.bs
The INPUT  file name : POSCAR or CONTCAR (selective)

