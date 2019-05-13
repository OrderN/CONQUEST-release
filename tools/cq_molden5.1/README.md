# How to use [Molden5.1](http://cheminf.cmbi.ru.nl/molden/) for editing [Conquest](http://www.order-n.org) coord file
This tool is to help Conquest users which want to:
- generate from scratch molecular stuctures
- edit xyz, pdb... molecular/crystal structures 
- convert/check the format without pain 

Typical application is for editing molecules/DNA strand/proteins... and save it directly into `Conquest_coord` format.

Remarks:
- this is not intended to handle millions of atoms
- this is more to help you for input easy edit (without scripts ;-)
- help for fast vizual inspection that the structure provided to the CQ code is correct!

Move to and edit the makefile in `src/`, then:

`$ make`  
`$ molden -l` or `$ gmolden` (activate OpenGL)

## Examples

### From pdb file
with the X-Ray crystal structure of HMGB1 domain A bound to a cisplatin, from [*Nature* 399, 708-712 (1999)](https://www.rcsb.org/structure/1CKT)

- in `examples/1.from_pdb/` directory, run:
`$ gmolden 1ckt.pdb` 
- in *Draw Mode* click on buttons `Solid` and `StickColor` to have a
  nice overview
- click on the unit cell button, then `Molecule` and `Mol+Cell` to
vizualise the initial guess for the unit cell
- you can setup cell dimensions with `Edit Cell Parameters`
- then `Write` and `write Conquest`

You should find the corresponding `Conquest_coord` file in the current
directory. You can check with:

`$ gmolden Conquest_coord` 

Remark (need `wget` package installed):
you can also load the `pdb` file directly in `molden`by:
- clicking on `Read`
- entering the `pdb` code; in this example, the code is `1CKT`
- clicking on `Get PDB`

### From Conquest_coord file

- in `examples/2.from_conquest/` directory, run:
`$ gmolden Conquest_coord` 

**WARNING**: the `Conquest_coord` file must contain the list of species
following the numbering of the 4th column 

Example with CH<sub>3</sub>COOH in a cube of length 10 Ang.

`    10.000000000000       0.000000000000       0.000000000000`
`      0.000000000000     10.000000000000       0.000000000000`
`      0.000000000000       0.000000000000      10.000000000000`
`           8     O  C  H`
`    0.5090523420     0.5000002846     0.4969778870        2   T  T  T`
`    0.3834785868     0.5000002846     0.4244778461        2   T  T  T`
`    0.5090523420     0.5000002846     0.6369779668        1   T  T  T`
`    0.6147075016     0.5000002846     0.4359778525        1   T  T  T`
`    0.5986192518     0.5000002846     0.6686446516        3   T  T  T`
`    0.4033778282     0.5000002846     0.3174113084        3   T  T  T`
`    0.3263738563     0.4110837574     0.4507860993        3   T  T  T`
`    0.3263738563     0.5889168118     0.4507860993        3   T  T  T`









