# How to use [Molden5.1](http://cheminf.cmbi.ru.nl/molden/) for
editing [Conquest](http://www.order-n.org) coord file
This tool is to help Conquest users which want to:
- generate from scratch molecular stuctures
- edit xyz, pdb, cif... molecular/crystal structures 
- convert/check without pain the format

Typical application is for building DNA sequence (in a water box) and save it directly into Conquest input format.

Remarks:
- this is not intended to handle millions of atoms
- this is more to help you for input easy edit (without scripts ;-)
- help for fast vizual inspection that the structure provided to the CQ code is correct!

Move an edit the makefile in `src/`, then:

`$ make`

`$ molden -l` or `$ gmolden` (activate OpenGL)

## Examples

### From pdb file
with the X-Ray crystal structure of HMGB1 domain A bound to a cisplatin, from [*Nature* 399, 708-712 (1999)](https://www.rcsb.org/structure/1CKT)

- In `examples/1.from_pdb/` directory, run:
`$ gmolden 1ckt.pdb` 
- In *Draw Mode* click on buttons `Solid` and `StickColor` to have a
  nice overview
- Click on the unit cell button, then `Molecule` and `Mol+Cell` to
vizualise the initial guess for the unit cell
- You can setup cell dimensions with `Edit Cell Parameters`
- Then `Write` and `write Conquest`

You should find the corresponding `Conquest_coord` file in the current
directory.

You can check with:
`$ gmolden Conquest_coord` 

Remark: (need `wget` package installed)

You can also load the `pdb` file directly in `molden`by:
- clicking on `Read`
- entering the `pdb` code, here `1CKT`
- and clicking on `Get PDB`

### From Conquest coord file
**WARNING**: the `Conquest_coord` file must contain the list of species
following the order of appearance in the input

`Example with H<sub>2</sub>O:

`Conquest_input` read:



