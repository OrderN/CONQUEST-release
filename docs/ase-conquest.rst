.. _ase-conquest:

==========================================
Managing Conquest with ASE
==========================================

Below we give an introduction how to setup the ASE environment with respect
to CONQUEST repository along with a few examples of ASE/Conquest capabilities.
We assume that a python script or jupyter-notebook is used. 

Environment variables
---------------------

The script will need to set environmental variables specifying the
locations of the CONQUEST executable ``Conquest``, and if required, the basis
set generation executable ``MakeIonFiles`` and pseudopotential database.
These variables are:

* ``ASE_ESPRESSO_COMMAND``: the Conquest executable command including MPI/openMPI prefix.
* ``CQ_PP_PATH``: the PAO path directory to where are located the the ``.ion`` files.
* (optional) ``CQ_GEN_BASIS_CMD`` : the PAO generation executable ``MakeIonFiles``.

..
  Given the Conquest root directory ``CQ_ROOT`` which contains

..
    ::
    
    bin/
    docs/
    pseudo-and-pao/
    src/
    testsuite/
    tools/
    
Given the Conquest root directory ``CQ_ROOT``, initialisation might look to something like
::

    import os

    CQ_ROOT = 'PATH_TO_CONQUEST_DIRECTORY'
    
    os.environ['ASE_CONQUEST_COMMAND'] = 'mpirun -np 4 '+CQ_ROOT+'/bin/Conquest'
    os.environ["CQ_GEN_BASIS_CMD"] = CQ_ROOT+'/bin/MakeIonFiles"
    os.environ['CQ_PP_PATH'] = CQ_ROOT+'/pseudo-and-pao/'


Pseudopotential/PAO files
-------------------------

Conquest atomic pseudotential and basis functions are store in the ``.ion`` 
files which will ne referred as to PAO files. Provided the pseudopotential files ``.pot`` available in ``CQ_PP_PATH``, 
automatic generation of numerical PAOs is possible using the program :ref:`MakeIonFiles <generating_paos>` 
available from the Conquest package.

Provided the PAO files, the basis set is specified through a python dictionary, 
for example:
::

    basis = {'O' : {'file': 'O_SZP.ion'},
             'H' : {'file': 'H_SZP.ion'},
             'C' : {'file': 'C_SZP.ion'}}
             
In this case they are all assumed to be obtained from Hamann pseudopotentials,
which are the default. Knowing the the exchange and correlation functional <XC>
from the Conquest input (*vide infra*) and the chemical symbol <X>, the *Calculator*
will search the ``.ion`` file in different places:
::

    CQ_PP_PATH
    CQ_PP_PATH/lib/
    CQ_PP_PATH/<XC>/<X>/
    
including the current directory and the ASE working directory (*vide infra*). If
your PAO file is located in a different place you can include the path in the
basis dictionary:
::

    basis = {'O' : {'file': 'O_SZP.ion',
                    'directory': '<PATH_TO_FILE>'},
             'H' : {'file' : 'H_SZP.ion'},
             'C' : {'file' : 'C_SZP.ion'}}
 

For generating the PAO files, the keyword ``gen_basis`` should be set to ``True``
(default is ``False``) and the size be provided (default is ``medium``). 
For instance:
::

    basis = {'O' : {'gen_basis' : True,
                    'basis_size': 'small'},
             'C' : {'gen_basis' : True,
                    'basis_size': 'medium'},                             
             'H' : {'gen_basis' : True,
                    'basis_size': 'large'}}

will create the ``O.ion``, ``C.ion`` and ``H.ion`` files where ``small``, ``medium`` and ``large`` 
are :ref:`default size basis set <pao_gen_default>`. You are allowed to choose the functional and 
to add options for basis set generation:
::

    basis = {'H' : {'gen_basis' : True,
                    'basis_size': 'small',
                    'xc'        : 'LDA',
                    'Atom.Perturbative_Polarised': False}}
                    
.. note::
    Only Hamann pseudopotentials for LDA, PBE and PBEsol are available within
    the CONQUEST distribution. For using other functionals see
    :ref:`Generating new pseudopotentials <pao_gen_pseudo>`.

Go to :ref:`top <ext-tools>`.

.. _et_ase_input:

Keywords for generating the Conquest_input file
+++++++++++++++++++++++++++++++++++++++++++++++

The *Calculator* object contains a dictionray containing a small number of
mandatory keywords, listed below:

::

    default_parameters = {
        'grid_cutoff'   : 100,     # DFT defaults
        'kpts'          : None,
        'xc'            : 'PBE',
        'scf_tolerance' : 1.0e-6,
        'nspin'         : 1,
        'general.pseudopotentialtype' : 'Hamann', # CONQUEST defaults
        'basis.basisset'              : 'PAOs',
        'io.iprint'                   : 2,
        'io.fractionalatomiccoords'   : True,
        'mine.selfconsistent'         : True,
        'sc.maxiters'                 : 50,
        'atommove.typeofrun'          : 'static',
        'dm.solutionmethod'           : 'diagon'}

The first five key/value pairs are special DFT parameters, the grid cutoff, the
k-point mesh, the exchange-correlation functional, the SCF tolerance and the
number of spins respectively. The rest are CONQUEST-specific input flags.

The atomic species blocks are handled slightly differently, with a dictionary of
their own. If the ``.ion`` files are present in the calculation directory, they
can be specified as follows:

::

  basis = {"H": {"valence_charge": 1.0,
                 "number_of_supports": 1,
                 "support_fn_range": 6.9},
           "O": {"valence_charge": 6.0,
                 "number_of_supports": 4,
                 "support_fn_range": 6.9}}

If the basis set ``.ion`` files are present in the directory containing the ASE
script are pressent and are named ``element.ion``, then the relevant parameters
will be parsed from the ``.ion`` files and included when the input file is
written and this dictionary can be omitted. It is more important when, for
example, setting up a multisite calculation, when the number of contracted
support functions is different from the number in the ``.ion`` file.

ASE can also invoke the CONQUEST basis set generation tool, although care should
be taken when generating basis sets:

::

  basis = {"H": {"basis_size": "minimal",
                 "pseudopotential_type": hamann",
                 "gen_basis": True},
           "O": {"basis_size": "minimal",
                 "pseudopotential_type": hamann",
                 "gen_basis": True}}

Finally, non-mandatory input flags can be defined in a new dictionary, and
passed as an expanded set of keyword arguments.

::

  conquest_flags = {'IO.Iprint'         : 1,         # CONQUEST keywords
                    'DM.SolutionMethod' : 'ordern',
                    'DM.L_range'        : 8.0,
                    'minE.LTolerance'   : 1.0e-6}

Here is an example, combining the above. We set up a cubic diamond cell
containing 8 atoms, and perform a single point energy calculation using the
order(N) method (the default is diagonalisation, so we must specify all of the
order(N) flags). We don't define a basis set, instead providing keywords that
specify that a minimal basis set should be constructed using the MakeIonFiles
basis generation tool.

::

  from ase.build import bulk
  from ase.calculators.conquest import Conquest

  os.environ["ASE_CONQUEST_COMMAND"] = "mpirun -np 4 Conquest_master"
  os.environ["CQ_PP_PATH"] = "/Users/zamaan/Conquest/PPDB/"
  os.environ["CQ_GEN_BASIS_CMD"] = "MakeIonFiles"

  diamond = bulk('C', 'diamond', a=3.6, cubic=True)  # The atoms object
  conquest_flags = {'IO.Iprint'         : 1,         # Conquest keywords
                    'DM.SolutionMethod' : 'ordern',
                    'DM.L_range'        : 8.0,
                    'minE.LTolerance'   : 1.0e-6}
  basis = {'C': {"basis_size"           : 'minimal', # Generate a minimal basis
                "gen_basis"             : True,
                "pseudopotential_type"  : "hamann"}}

  calc = Conquest(grid_cutoff = 80,    # Set the calculator keywords
                  xc="LDA",
                  self_consistent=True,
                  basis=basis,
                  nspin=1,
                  **conquest_flags)
  diamond.set_calculator(calc)             # attach the calculator to the atoms object
  energy = diamond.get_potential_energy()  # calculate the potential energy

Go to :ref:`top <ext-tools>`.

.. _et_ase_mssf:

Multisite support functions
+++++++++++++++++++++++++++

Multisite support functions require a few additional keywords in the atomic
species block, which can be specified as follows:

::

  basis = {'C': {"basis_size": 'medium',
                 "gen_basis": True,
                 "pseudopotential_type": "hamann",
                 "Atom.NumberofSupports": 4,
                 "Atom.MultisiteRange": 7.0,
                 "Atom.LFDRange": 7.0}}

Note that we are constructing a DZP basis set (size medium) with 13 primitive
support functions using ``MakeIonFiles``, and contracting it to multisite basis
of 4 support functions. The calculation requires a few more input flags, which
are specified in the ``other_keywords`` dictionary:

::

  other_keywords = {"Basis.MultisiteSF": True,
                    "Multisite.LFD": True,
                    "Multisite.LFD.Min.ThreshE": 1.0e-7,
                    "Multisite.LFD.Min.ThreshD": 1.0e-7,
                    "Multisite.LFD.Min.MaxIteration": 150,
                    }

Go to :ref:`top <ext-tools>`.

.. _et_ase_load_dm:

Loading the K/L matrix
++++++++++++++++++++++
   
Most calculation that involve incrementally moving atoms (molecular dynamics,
geometry optimisation, equations of state, nudged elastic band etc.) can be made
faster by using the K or L matrix from a previous calculation as the initial
guess for a subsequent calculation in which that atoms have been moved slightly.
This can be achieved by first performing a single point calculation to generate
the first K/L matrix, then adding the following keywords to the calculator:

::

  other_keywords = {"General.LoadL": True,
                    "SC.MakeInitialChargeFromK": True}

These keywords respectively cause the K or L matrix to be loaded from file(s)
``Kmatrix.i**.p*****``, and the initial charge density to be constructed from
this matrix. In all subsequent calculations, the K or L matrix will be written
at the end of the calculation and used as the initial guess for the subsequent
ionic step.

Go to :ref:`top <ext-tools>`.

.. _et_eos:

Equation of state
+++++++++++++++++

The following code computes the equation of state of diamond by doing single
point calculations on a uniform grid of the ``a`` lattice parameter. It then
interpolates the equation of state and uses ``matplotlib`` to generate a plot.

::

  import scipy as sp
  from ase.build import bulk
  from ase.io.trajectory import Trajectory
  from ase.calculators.conquest import Conquest


  # Construct a unit cell
  diamond = bulk('C', 'diamond', a=3.6, cubic=True)

  basis = {'C': {"basis_size": 'minimal', 
                 "gen_basis": True,
                 "pseudopotential_type": "hamann"}}
  calc = Conquest(grid_cutoff = 50,
                  xc = "LDA",
                  basis = basis,
                  kpts = [4,4,4]}
  diamond.set_calculator(calc)

  cell = diamond.get_cell()
  traj = Trajectory('diamond.traj', 'w') # save all results to trajectory

  for x in sp.linspace(0.95, 1.05, 5):   # grid for equation of state
    diamond.set_cell(cell*x, scale_atoms=True)
    diamond.get_potential_energy()
    traj.write(diamond)

  from ase.io import read
  from ase.eos import EquationOfState

  configs = read('diamond.traj@0:5')
  volumes = [diamond.get_volume() for diamond in configs]
  energies = [diamond.get_potential_energy() for diamond in configs]
  eos = EquationOfState(volumes, energies)
  v0, e0, B = eos.fit()

  import matplotlib
  eos.plot('diamond-eos.pdf')    # Plot the equation of state

Go to :ref:`top <ase-conquest>`.


