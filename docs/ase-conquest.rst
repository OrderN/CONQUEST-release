.. _ase-conquest:

==========================================
Managing Conquest with ASE
==========================================

Below we give an introduction how to setup the ASE environment with respect 
to CONQUEST repository along with a few examples of ASE/Conquest capabilities.
We assume that a python script or jupyter-notebook is used. 

.. _ase_setup:

Setup
-----

.. _ase_env:

Environment variables
+++++++++++++++++++++

The script will need to set environmental variables specifying the
locations of the CONQUEST executable ``Conquest``, and if required, the basis
set generation executable ``MakeIonFiles`` and pseudopotential database.
These variables are:

* ``ASE_CONQUEST_COMMAND``: the Conquest executable command including MPI/openMPI prefix.
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


Go to :ref:`top <ase-conquest>`.

.. _ase_pao:

Pseudopotential/PAO files
+++++++++++++++++++++++++

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

.. warning::
    Generating polarised PAOs for some atoms can be problematic (mainly group I 
    and II). Please review carefuly the ``MakeIonFiles`` input files
    named ``Conquest_ion_input`` which are collected in ``CQ_PP_PATH/<XC>/<X>/``
    if you are not sure about what you are doing, and check your PAOs.


Go to :ref:`top <ase-conquest>`.

.. _ase_calculator:

CONQUEST Calculator
--------------------

The CONQUEST *Calculator* class can be invoked from the ase Calculator set as described
in the example below:
::

    from ase.calculators.conquest import Conquest

A minimal example is given below for setting the CONQUEST *Calculator* (named ``calc``) 
of the `ASE Atoms object <https://wiki.fysik.dtu.dk/ase/ase/atoms.html#module-ase.atoms>`_ 
named `struct`:
::

    from ase.calculators.conquest import Conquest
    from ase.build import bulk

    struct = bulk('NaCl', crystalstructure='rocksalt', a=5.71, cubic=True)    
    basis  = {'Cl' : {'file' : 'Cl.ion'}, 'Na' : {'file' : 'Na.ion'}}
        
    calc = Conquest(basis=basis,atoms=struct)

or, equivalently,
::

    from ase.calculators.conquest import Conquest
    from ase.build import bulk

    struct = bulk('NaCl', crystalstructure='rocksalt', a=5.71, cubic=True)    
    basis  = {'Cl' : {'file' : 'Cl.ion'}, 'Na' : {'file' : 'Na.ion'}}
        
    struct.calc = Conquest(basis=basis)


In basic calculate mode (compute ``energy``), the *Calculator* comes with 3 *methods*:    

* ``write_input()``: 
    this function will setup the :ref:`input files <setting_up>`. For CONQUEST, the PAO basis will be generated/copied with respect to the dictionary key/value pairs, and ``Conquest_input`` file including the calculation parameters will be written, a long with the ``coordinate`` file, containing the lattice vectors (in Bohr Unit) and atomic positions (in fractional coordinates).

* ``execute()``: 
    this function execute the calculation. For CONQUEST, it will launch the ``ASE_CONQUEST_COMMAND`` setup in :ref:`the environment varaibles <ase_env>`.

* ``read_results()``: 
    this function post-process the the output file. For CONQUEST, the ``energy``, ``forces``, ``stress`` and ``eigenvalues`` will be extracted from the ``Conquest_out_ase`` output file.

.. note::
    The funtion ``read_results()`` operate on the  ``Conquest_out_ase`` file. 
    This output file is not created by default by CONQUEST. If you want to post-process
    a calculation with an input generated by hand you must add ``IO.WriteOutToASEFile True``
    in ``conquest_input``.
    

The **indirect** way for managing CONQUEST calculation with ASE is:
::
    
    struct.calc = Conquest(basis=basis)    

    struct.calc.write_input(struct)    
    struct.calc.execute()    
    struct.calc.read_results(struct)    


where ``struct.calc.execute()`` can be ignored when, for instance, the calculation
is performed on a supercomputer and the output file is then copied back to the current 
directory for post-processing.


The **direct** way is simply:
::
    
    struct.calc = Conquest(basis=basis)    

    struct.calc.calculate(struct)

or, equivalently,
::
    
    struct.calc = Conquest(basis=basis)    

    struct.get_potential_energy() 

Go to :ref:`top <ase-conquest>`.

.. _ase_input:

Keywords for generating the Conquest_input file
-----------------------------------------------

In principle all the `Conquest input parameters <https://conquest.readthedocs.io/en/latest/input_tags.html>`_
can be added to ``Conquest_out_ase`` using key/value pairs in a dictionary. There are 3 class of parameters:

* mandatory : they are parsed to the *Calculator* and have no defaults ; there are mandatory.
* important : they are parsed to the *Calculator* they can be freely modified. Some of them are pure ASE keywords.
* defaults : they are set as defaults ; some of them must not be modified. They are read by the *Calculator* through a dictionay ``conquest_flags``.


Mandatory keywords
++++++++++++++++++


=================== =====================  ===============  ================================
keyword             type                   default value    description
=================== =====================  ===============  ================================
``atoms``           ``atoms``              None             an atoms object constructed either via ASE or read from an input
``basis``           ``dict``               None             a dictionary specifying the pseudopotential/basis files
=================== =====================  ===============  ================================

Important keywords
++++++++++++++++++


=================== ==========================  =====================  ===============  ================================
keyword             CONQUEST equivalence        type                   default value    description
=================== ==========================  =====================  ===============  ================================
``directory``       None                        ``str``                None             directory used for storing input/output and calculation files
``label``           None                        ``str``                None             basename for working files (only used by ASE, eg. NEB)
``kpts``            None                        ``list`` or ``tuple``  None             k-points grid ; converted to CONQUEST Monkhorst-Pack grid 
``grid_cutoff``     ``Grid.GridCutoff``         ``float``              100              integration grid in Ha
``xc``              ``General.FunctionalType``  ``str``                'PBE'            exchange and correlation functional
``self_consistent`` ``minE.SelfConsistent``     ``bool``               True             choose either SCF or non-SCF
``scf_tolerance``   ``minE.SCTolerance``        ``float``              1e-6             Self-consistent-field convergence tolerance in Ha
``nspin``           ``Spin.SpinPolarised``      ``int``                 1               spin polarisation: 1 for unpolarized or 2 for polarised
``conquest_flags``  None                        ``dict``               None             other CONQUET keyword arguments
=================== ==========================  =====================  ===============  ================================


Defaults keywords
+++++++++++++++++

===============================  =========  ===============  ================================
keyword                          type       default value    description
===============================  =========  ===============  ================================
``IO.WriteOutToASEFile``         ``bool``   True             write ASE output file ; **must always be True when using ASE for post-processing**
``IO.Iprint``                    ``int``    1                verbose for the output ; **must always be 1 when using ASE for post-processing**
``DM.SolutionMethod``            ``str``    'diagon'         'diagon' stands for diagonalisation other is 'ordern' (base on density matrix)
``General.PseudopotentialType``  ``str``    'Hamann'         kind of pseudopotential other type are 'siesta' and 'abinit'
``SC.MaxIters``                  ``int``    50               maximum number SCF cycles
``AtomMove.TypeOfRun``           ``str``    'static'         'static' stands for single (non)SCF other are 'md' or optimisation algorithms.
``Diag.SmearingType``            ``int``    1                1 for Methfessel-Paxton ; 0 for Fermi-Dirac
``Diag.kT``                      ``float``  0.001            smearing temperature in Ha
===============================  =========  ===============  ================================

..
  ``io.fractionalatomiccoords``    ``bool``   True             atomic coordinates format for the structure file (fractional or cartesian)
  ``basis.basisset``               ``str``    'PAOs'           type of basis set ; always 'PAOs' with ASE 

Some examples
-------------

An example of more advanced Calculator setup is given below for a SCF calculation on BCC-Na
where for a PBE calculation using a k-point grid of :math:`6\times 6\times6` using the Fermi-Dirac
distribution for the occupation with a smearing of 0.005 Ha::

    struct = bulk('Na', crystalstructure='bcc', a=4.17, cubic=True)    
    basis  = {'Na' : {'file' : 'NaCQ.ion'}}
    
    conquest_flags = {'Diag.SmearingType': 0,
                      'Diag.kT'          : 0.005}
                      
    struct.calc = Conquest(directory      = 'Na_bcc_example',
                           grid_cutoff    = 90.0,
                           self_consistent= True,
                           xc    = 'PBE',
                           basis = basis,
                           kpts  = [6,6,6],
                           nspin = 1,
                           **conquest_flags)
    
    struct.get_potential_energy()

..
    A dictionary containing a small number of mandatory keywords, listed below:

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

..
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

Finally, defaults and other input flags can be defined in a new dictionary, and
passed as an expanded set of keyword arguments.
::

  conquest_flags = {'DM.SolutionMethod' : 'ordern',
                    'DM.L_range'        : 12.0,
                    'minE.LTolerance'   : 1.0e-6}

Here is an example, combining the above. We set up a cubic diamond cell
containing 8 atoms, and perform a single point energy calculation using the
order(*N*) method (the default is diagonalisation, so we must specify all of the
order(*N*) flags). We don't define a basis set, instead providing keywords that
specify that a minimal basis set should be constructed using the MakeIonFiles
basis generation tool.

::

  import os
  from ase.build import bulk
  from ase.calculators.conquest import Conquest

  CQ_ROOT = 'PATH_TO_CONQUEST_DIRECTORY'
    
  os.environ['ASE_CONQUEST_COMMAND'] = 'mpirun -np 4 '+CQ_ROOT+'/bin/Conquest'
  os.environ["CQ_GEN_BASIS_CMD"] = CQ_ROOT+'/bin/MakeIonFiles"
  os.environ['CQ_PP_PATH'] = CQ_ROOT+'/pseudo-and-pao/'

  diamond = bulk('C', 'diamond', a=3.6, cubic=True)  # The atoms object
  conquest_flags = {'DM.SolutionMethod' : 'ordern',  # Conquest keywords
                    'DM.L_range'        : 12.0,
                    'minE.LTolerance'   : 1.0e-6}
                    
  basis = {'C': {'basis_size' : 'minimal', # Generate a minimal basis
                 'gen_basis'  : True}

  calc = Conquest(grid_cutoff = 80,    # Set the calculator keywords
                  xc = 'PBE',
                  self_consistent=True,
                  basis = basis,
                  nspin = 1,
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
                 "gen_basis": True}}
                 
  calc = Conquest(grid_cutoff = 50,
                  xc = "PBE",
                  basis = basis,
                  kpts = [4,4,4])
                  
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


