      program ambfor
c
c  AMBER General Force Field for organic mol. (June, 2003)
c  The AMBFOR GAFF implementation:
c  
c   1 X   Placeholder for wildcard
c   2 c   Sp2 C carbonyl group 
c   3 c1  Sp C
c   4 c2  Sp2 C  
c   5 c3  Sp3 C
c   6 ca  Sp2 C in pure aromatic systems
c   7 cp  Head Sp2 C that connect two rings in biphenyl sys. 
c   8 cq  Head Sp2 C that connect two rings in biphenyl sys. identical to cp 
c   9 cc  Sp2 carbons in non-pure aromatic systems
c  10 cd  Sp2 carbons in non-pure aromatic systems, identical to cc
c  11 ce  Inner Sp2 carbons in conjugated systems
c  12 cf  Inner Sp2 carbons in conjugated systems, identical to ce
c  13 cg  Inner Sp carbons in conjugated systems
c  14 ch  Inner Sp carbons in conjugated systems, identical to cg
c  15 cx  Sp3 carbons in triangle systems
c  16 cy  Sp3 carbons in square systems
c  17 cu  Sp2 carbons in triangle systems
c  18 cv  Sp2 carbons in square systems
c  19 h1  H bonded to aliphatic carbon with 1 electrwd. group  
c  20 h2  H bonded to aliphatic carbon with 2 electrwd. group 
c  21 h3  H bonded to aliphatic carbon with 3 electrwd. group 
c  22 h4  H bonded to non-sp3 carbon with 1 electrwd. group 
c  23 h5  H bonded to non-sp3 carbon with 2 electrwd. group 
c  24 ha  H bonded to aromatic carbon  
c  25 hc  H bonded to aliphatic carbon without electrwd. group 
c  26 hn  H bonded to nitrogen atoms
c  27 ho  Hydroxyl group
c  28 hp  H bonded to phosphate 
c  29 hs  Hydrogen bonded to sulphur 
c  30 hw  Hydrogen in water 
c  31 hx  H bonded to C next to positively charged group  
c  32 f   Fluorine
c  33 cl  Chlorine 
c  34 br  Bromine 
c  35 i   Iodine 
c  36 n   Sp2 nitrogen in amide groups
c  37 n1  Sp N  
c  38 n2  aliphatic Sp2 N with two connected atoms 
c  39 n3  Sp3 N with three connected atoms
c  40 n4  Sp3 N with four connected atoms 
c  41 na  Sp2 N with three connected atoms 
c  42 nb  Sp2 N in pure aromatic systems 
c  43 nc  Sp2 N in non-pure aromatic systems
c  44 nd  Sp2 N in non-pure aromatic systems, identical to nc
c  45 ne  Inner Sp2 N in conjugated systems
c  46 nf  Inner Sp2 N in conjugated systems, identical to ne
c  47 nh  Amine N connected one or more aromatic rings 
c  48 no  Nitro N  
c  49 o   Oxygen with one connected atom
c  50 oh  Oxygen in hydroxyl group
c  51 os  Ether and ester oxygen
c  52 ow  Oxygen in water 
c  53 p2  Phosphate with two connected atoms 
c  54 p3  Phosphate with three connected atoms, such as PH3
c  55 p4  Phosphate with three connected atoms, such as O=P(CH3)2
c  56 p5  Phosphate with four connected atoms, such as O=P(OH)3
c  57 pb  Sp2 P in pure aromatic systems 
c  58 pc  Sp2 P in non-pure aromatic systems
c  59 pd  Sp2 P in non-pure aromatic systems, identical to pc
c  60 pe  Inner Sp2 P in conjugated systems
c  61 pf  Inner Sp2 P in conjugated systems, identical to pe
c  62 px  Special p4 in conjugated systems
c  63 py  Special p5 in conjugated systems
c  64 s   S with one connected atom 
c  65 s2  S with two connected atom, involved at least one double bond  
c  66 s4  S with three connected atoms 
c  67 s6  S with four connected atoms 
c  68 sh  Sp3 S connected with hydrogen 
c  69 ss  Sp3 S in thio-ester and thio-ether
c  70 sx  Special s4 in conjugated systems
c  71 sy  Special s6 in conjugated systems
c  72 cz  Sp2 carbon in guadinium
c
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxion=2000)
      character line*132,fniun*132
      character argstr*75
      logical opfil,gargpl,osingl,dolbfgs,oqscal
cmpi      include 'mpif.h'
      common /fnam/ fniun, lenf
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      common /opts/ idebug,imon,iarc,ilog,iout,osingl,dolbfgs,oqscal
      logical usecut,usesw
      common /limits/ cutoff, cutof2,cuton, cuton2,conof3,usecut,usesw
      logical box,cell,addbox,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /passe/ emin
      common /passq/ scalq
      common /athlp/  iatoms, mxnat
      common /prot/ ireswr
      logical doint
      common /intr/ ecold,evold,ewold,iprog,doint
      common /convg/ icyc
      common /debug/ ideb
      common /h2oer/ numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      integer taskid
      common /mpih/ nproc,taskid
      character*3 ambstr
      character*2 gffstr
      common /ffstr/ ambstr(mxamb), gffstr(mxgff)
      data (ambstr(i),i=1,100) /
     & 'N  ','CT ','C  ','H  ','O  ','H1 ','N  ','CT ','C  ','H  ',
     & 'O  ','H1 ','CT ','HC ','N  ','CT ','C  ','H  ','O  ','H1 ',
     & 'CT ','HC ','CT ','HC ','CT ','HC ','N  ','CT ','C  ','H  ',
     & 'O  ','H1 ','CT ','HC ','CT ','HC ','CT ','HC ','CT ','HC ',
     & 'N  ','CT ','C  ','H  ','O  ','H1 ','CT ','HC ','CT ','HC ',
     & 'CT ','HC ','CT ','HC ','N  ','CT ','C  ','H  ','O  ','H1 ',
     & 'CT ','H1 ','OH ','HO ','N  ','CT ','C  ','H  ','O  ','H1 ',
     & 'CT ','H1 ','OH ','HO ','CT ','HC ','N  ','CT ','C  ','H  ',
     & 'O  ','H1 ','CT ','H1 ','SH ','HS ','N  ','CT ','C  ','H  ',
     & 'O  ','H1 ','CT ','H1 ','S  ','N  ','CT ','C  ','O  ','H1 '/
      data (ambstr(i),i=101,200) /
     & 'CT ','HC ','CT ','HC ','CT ','H1 ','N  ','CT ','C  ','H  ',
     & 'O  ','H1 ','CT ','HC ','CA ','CA ','HA ','CA ','HA ','CA ',
     & 'HA ','N  ','CT ','C  ','H  ','O  ','H1 ','CT ','HC ','CA ',
     & 'CA ','HA ','CA ','HA ','C  ','OH ','HO ','N  ','CT ','C  ',
     & 'H  ','O  ','H1 ','CT ','HC ','C* ','CW ','H4 ','CB ','NA ',
     & 'H  ','CN ','CA ','HA ','CA ','HA ','CA ','HA ','CA ','HA ',
     & 'N  ','CT ','C  ','H  ','O  ','H1 ','CT ','HC ','CC ','NA ',
     & 'H  ','CW ','H4 ','CR ','H5 ','NA ','H  ','N  ','CT ','C  ',
     & 'H  ','O  ','H1 ','CT ','HC ','CC ','NA ','H  ','CV ','H4 ',
     & 'CR ','H5 ','NB ','N  ','CT ','C  ','H  ','O  ','H1 ','CT '/
      data (ambstr(i),i=201,300) /
     & 'HC ','CC ','NB ','CW ','H4 ','CR ','H5 ','NA ','H  ','N  ',
     & 'CT ','C  ','H  ','O  ','H1 ','CT ','HC ','C  ','O2 ','N  ',
     & 'CT ','C  ','H  ','O  ','H1 ','CT ','HC ','C  ','O  ','N  ',
     & 'H  ','N  ','CT ','C  ','H  ','O  ','H1 ','CT ','HC ','CT ',
     & 'HC ','C  ','O2 ','N  ','CT ','C  ','H  ','O  ','H1 ','CT ',
     & 'HC ','CT ','HC ','C  ','O  ','N  ','H  ','N  ','CT ','C  ',
     & 'H  ','O  ','H1 ','CT ','HC ','CT ','H1 ','S  ','CT ','H1 ',
     & 'N  ','CT ','C  ','H  ','O  ','H1 ','CT ','HC ','CT ','HC ',
     & 'CT ','HC ','CT ','HP ','N3 ','H  ','N  ','CT ','C  ','H  ',
     & 'O  ','H1 ','CT ','HC ','CT ','HC ','CT ','H1 ','N2 ','H  '/
      data (ambstr(i),i=301,400) /
     & 'CA ','N2 ','H  ','N  ','CT ','C  ','H  ','O  ','H1 ','CT ',
     & 'HC ','CT ','HC ','CT ','HP ','N3 ','H  ','N  ','CT ','C  ',
     & 'H  ','O  ','CT ','HC ','N  ','CT ','C  ','H  ','O  ','H1 ',
     & 'CT ','HC ','CT ','HC ','C  ','O  ','C  ','H  ','O  ','CT ',
     & 'HC ','C  ','O  ','N  ','H  ','N  ','H  ','CT ','H1 ','N3 ',
     & 'CT ','C  ','H  ','O  ','H1 ','N3 ','CT ','C  ','H  ','O  ',
     & 'H1 ','N3 ','CT ','C  ','H  ','O  ','H1 ','N3 ','CT ','C  ',
     & 'H  ','O  ','H1 ','N3 ','CT ','C  ','H  ','O  ','H1 ','N3 ',
     & 'CT ','C  ','H  ','O  ','H1 ','N3 ','CT ','C  ','H  ','O  ',
     & 'H1 ','N3 ','CT ','C  ','H  ','O  ','H1 ','N3 ','CT ','C  '/
      data (ambstr(i),i=401,500) /
     & 'H  ','O  ','H1 ','N3 ','CT ','C  ','H  ','O  ','H1 ','CT ',
     & 'HP ','N3 ','CT ','C  ','H  ','O  ','H1 ','N3 ','CT ','C  ',
     & 'H  ','O  ','H1 ','N3 ','CT ','C  ','H  ','O  ','H1 ','N3 ',
     & 'CT ','C  ','H  ','O  ','H1 ','N3 ','CT ','C  ','H  ','O  ',
     & 'H1 ','N3 ','CT ','C  ','H  ','O  ','H1 ','N3 ','CT ','C  ',
     & 'H  ','O  ','H1 ','N3 ','CT ','C  ','H  ','O  ','H1 ','N3 ',
     & 'CT ','C  ','H  ','O  ','H1 ','N3 ','CT ','C  ','H  ','O  ',
     & 'H1 ','N3 ','CT ','C  ','H  ','O  ','H1 ','N3 ','CT ','C  ',
     & 'H  ','O  ','H1 ','N3 ','CT ','C  ','H  ','O  ','H1 ','N3 ',
     & 'CT ','C  ','H  ','O  ','H1 ','N3 ','CT ','C  ','H  ','O  '/
      data (ambstr(i),i=501,600) /
     & 'N  ','CT ','C  ','H  ','O2 ','H1 ','N  ','CT ','C  ','H  ',
     & 'O2 ','H1 ','N  ','CT ','C  ','H  ','O2 ','H1 ','N  ','CT ',
     & 'C  ','H  ','O2 ','H1 ','N  ','CT ','C  ','H  ','O2 ','H1 ',
     & 'N  ','CT ','C  ','H  ','O2 ','H1 ','N  ','CT ','C  ','H  ',
     & 'O2 ','H1 ','N  ','CT ','C  ','H  ','O2 ','H1 ','N  ','CT ',
     & 'C  ','H  ','O2 ','H1 ','N  ','CT ','C  ','O2 ','H1 ','N  ',
     & 'CT ','C  ','H  ','O2 ','H1 ','N  ','CT ','C  ','H  ','O2 ',
     & 'H1 ','N  ','CT ','C  ','H  ','O2 ','H1 ','N  ','CT ','C  ',
     & 'H  ','O2 ','H1 ','N  ','CT ','C  ','H  ','O2 ','H1 ','N  ',
     & 'CT ','C  ','H  ','O2 ','H1 ','N  ','CT ','C  ','H  ','O2 '/
      data (ambstr(i),i=601,659) /
     & 'H1 ','N  ','CT ','C  ','H  ','O2 ','H1 ','N  ','CT ','C  ',
     & 'H  ','O2 ','H1 ','N  ','CT ','C  ','H  ','O2 ','H1 ','N  ',
     & 'CT ','C  ','H  ','O2 ','H1 ','N  ','CT ','C  ','H  ','O2 ',
     & 'H1 ','N  ','CT ','C  ','H  ','O2 ','H1 ','N  ','CT ','C  ',
     & 'H  ','O2 ','H1 ','N  ','CT ','C  ','H  ','O2 ','OW ','HW ',
     & 'Li+','Na+','K+ ','Rb+','Cs+','Mg+','Ca+','Zn+','Cl-'/

c     some neutral residues

      data (ambstr(i),i=660,671) /
     & 'N  ','CT ','C  ','H  ','O  ','H1 ','CT ','HC ',
     & 'C  ','O  ','OH ','HO '/ 
      data (ambstr(i),i=672,685) /
     & 'N  ','CT ','C  ','H  ','O  ','H1 ','CT ','HC ',
     & 'CT ','HC ','C  ','O  ','OH ','HO '/ 
      data (ambstr(i),i=686,701) /
     & 'N  ','CT ','C  ','H  ','O  ','H1 ','CT ','HC ',
     & 'CT ','HC ','CT ','HC ','CT ','HP ','N3 ','H  '/ 

      data (ambstr(i),i=1001,1100) /
     & 'OS ','CT ','H1 ','H1 ','CT ','H1 ','OS ','CT ','H2 ','CT ',
     & 'H1 ','CT ','H1 ','OH ','HO ','OS ','N* ','CB ','CB ','NB ',
     & 'CK ','NC ','CQ ','NC ','CA ','H5 ','N2 ','H  ','H  ','H5 ',
     & 'OS ','CT ','H1 ','H1 ','CT ','H1 ','OS ','CT ','H2 ','CT ',
     & 'H1 ','CT ','H1 ','OH ','HO ','OS ','N* ','CB ','CB ','NB ',
     & 'CK ','NC ','CA ','NA ','C  ','H  ','N2 ','H  ','H  ','O  ',
     & 'H5 ','OS ','CT ','H1 ','H1 ','CT ','H1 ','OS ','CT ','H2 ',
     & 'CT ','H1 ','CT ','H1 ','OH ','HO ','OS ','N* ','C  ','NC ',
     & 'CA ','CM ','CM ','O  ','N2 ','H  ','H  ','HA ','H4 ','OS ',
     & 'CT ','H1 ','H1 ','CT ','H1 ','OS ','CT ','H2 ','CT ','H1 '/

      data (ambstr(i),i=1101,1200) /
     & 'CT ','H1 ','OH ','HO ','OS ','N* ','C  ','NA ','C  ','CM ',
     & 'CM ','O  ','H  ','O  ','HA ','H4 ','OS ','CT ','H1 ','H1 ',
     & 'CT ','H1 ','OS ','CT ','H2 ','CT ','H1 ','CT ','HC ','HC ',
     & 'OS ','N* ','CB ','CB ','NB ','CK ','NC ','CQ ','NC ','CA ',
     & 'H5 ','N2 ','H  ','H  ','H5 ','OS ','CT ','H1 ','H1 ','CT ',
     & 'H1 ','OS ','CT ','H2 ','CT ','H1 ','CT ','HC ','HC ','OS ',
     & 'N* ','CB ','CB ','NB ','CK ','NC ','CA ','NA ','C  ','H  ',
     & 'N2 ','H  ','H  ','O  ','H5 ','OS ','CT ','H1 ','H1 ','CT ',
     & 'H1 ','OS ','CT ','H2 ','CT ','H1 ','CT ','HC ','HC ','OS ',
     & 'N* ','C  ','NC ','CA ','CM ','CM ','O  ','N2 ','H  ','H  '/

      data (ambstr(i),i=1201,1253) /
     & 'HA ','H4 ','OS ','CT ','H1 ','H1 ','CT ','H1 ','OS ','CT ',
     & 'H2 ','CT ','H1 ','CT ','HC ','HC ','OS ','N* ','C  ','NA ',
     & 'C  ','CM ','CM ','O  ','H  ','O  ','CT ','HC ','H4 ','P  ',
     & 'O2 ','OH ','HO ','OS ','P  ','O2 ','OH ','HO ','OS ','P  ',
     & 'O2 ','P  ','O2 ','OH ','HO ','OS ','P  ','O2 ','OH ','HO ',
     & 'OS ','P  ','O2 '/

c modified nucl

      data (ambstr(i),i=1254,1284) /
     & 'OS','CT','H1','H1','CT','H1','OS','CT','H2','CT','H1',
     & 'CT','H1','OH','HO','OS','N*','CB','CB','NB','CK','NC',
     & 'CQ','N*','C ','H1','N ','H ','H5','CT','H1'/

      data (ambstr(i),i=1285,1311) /
     & 'OS', 'CM', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'N*', 'C ', 'NC',
     & 'CA', 'CT', 'CM', 'H4', 'O ', 'N2', 'H ', 'CT', 'HC'/

      data (ambstr(i),i=1312,1338) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OS', 'OS', 'N*', 'C ', 'O ', 'NC',
     & 'CA', 'N2', 'H ', 'CM', 'HA', 'CM', 'H4', 'CT', 'H1'/
 
      data (ambstr(i),i=1339,1369) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'N*', 'CB', 'CB',
     & 'NB', 'CK', 'H5', 'NC', 'CA', 'NA', 'H ', 'C ', 'O ',
     & 'N ', 'H ', 'CT', 'H1'/ 

      data (ambstr(i),i=1370,1399) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'N*', 'CB', 'CB',
     & 'NB', 'CK', 'NC', 'CA', 'NA', 'C ', 'O ', 'H ', 'H5',
     & 'N2', 'CT', 'H1'/

      data (ambstr(i),i=1400,1430) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'N*', 'CB', 'CB',
     & 'NA', 'CK', 'NC', 'CA', 'NA', 'C ', 'O ', 'H ', 'H5',
     & 'N2', 'H ', 'CT', 'H1'/

      data (ambstr(i),i=1431,1460) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OS', 'OS', 'N*', 'CB', 'CB', 'NB',
     & 'CK', 'NC', 'CA', 'NA', 'C ', 'O ', 'H ', 'H5', 'N2',
     & 'H ', 'CT', 'H1'/ 

      data (ambstr(i),i=1461,1511) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H1', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'N*', 'CB', 'CB',
     & 'NB', 'CK', 'N*', 'CB', 'N*', 'C ', 'O ', 'NB', 'CC',
     & 'CC', 'H5', 'CT', 'H1', 'CT', 'HC', 'CT', 'HC', 'CT',
     & 'HC', 'CT', 'H1', 'C ', 'O ', 'OS', 'CT', 'H1', 'N ',
     & 'H ', 'C ', 'O ', 'OS', 'CT', 'H1'/ 

      data (ambstr(i),i=1512,1537) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'N ', 'C ', 'N ',
     & 'C ', 'CT', 'CT', 'O ', 'O ', 'H ', 'HC', 'H1'/ 

      data (ambstr(i),i=1538,1564) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'N*', 'C ', 'NA',
     & 'C ', 'CM', 'CM', 'O ', 'O ', 'H ', 'H4', 'CT', 'HC'/

      data (ambstr(i),i=1565,1590) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H1', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'CM', 'CM', 'NA',
     & 'C ', 'NA', 'C ', 'O ', 'O ', 'H4', 'H ', 'H '/

c GAFF
      data gffstr /
     &  'x ', 'c ', 'c1', 'c2', 'c3', 'ca', 'cp', 'cq', 'cc', 'cd',
     &  'ce', 'cf', 'cg', 'ch', 'cx', 'cy', 'cu', 'cv', 'h1', 'h2',
     &  'h3', 'h4', 'h5', 'ha', 'hc', 'hn', 'ho', 'hp', 'hs', 'hw',
     &  'hx', 'f ', 'cl', 'br', 'i ', 'n ', 'n1', 'n2', 'n3', 'n4',
     &  'na', 'nb', 'nc', 'nd', 'ne', 'nf', 'nh', 'no', 'o ', 'oh',
     &  'os', 'ow', 'p2', 'p3', 'p4', 'p5', 'pb', 'pc', 'pd', 'pe',
     &  'pf', 'px', 'py', 's ', 's2', 's4', 's6', 'sh', 'ss', 'sx',
     &  'sy', 'cz'/

      idebug  = 0
      osingl  = .false.
      oqscal  = .false.
      scalq   = 1.0e0
      dolbfgs   = .true.
      box     = .false.
      addbox  = .false.
      usecut  = .false.
      usesw   = .false.
      fast = .true.
      doint = .true.
      iprog = 1
      imon = 0
      iarc = 0
      iout = 1
      ilog = 0
      gtol = 0.03e0
      nsd  = 100000
      ireswr = 0
      iff99sb = 1
      iff99bc = 0
      ionpl = 1
      nproc = 1
      ideb = 0
  
      n = iargc()

      if (idebug.eq.1) then
         if (.not.opfil(61,"log",3,1,0,0)) then
            stop
         endif
      endif

      icnt = 1
      lenf = -1
      do while (icnt.le.n)
         call getarg(icnt,line)
         if (idebug.eq.1) write(61,*) "=",line,"="
         if (line(1:1).eq.'-') then

            if (line(1:2).eq.'-v') then
               idebug = 1
            elseif (line(1:2).eq.'-m') then
               imon = 1
            elseif (line(1:2).eq.'-M') then
               imon = 2
            elseif (line(1:2).eq.'-b') then
               addbox = .true.
            elseif (line(1:2).eq.'-n') then
               ilog = 1
            elseif (line(1:2).eq.'-i') then
               fast = .false.
            elseif (line(1:2).eq.'-N') then
               usecut = .true.
            elseif (line(1:2).eq.'-a') then
               iarc = 1
            elseif (line(1:2).eq.'-p') then
               dolbfgs = .false.
            elseif (line(1:2).eq.'-s') then
               osingl = .true.
            elseif (line(1:2).eq.'-S') then
               oqscal = .true.
               if (gargpl('-S',icnt,line,argstr)) then
                  scalq = reada(argstr,1,len(argstr))
               endif
            elseif (line(1:2).eq.'-x') then
               iff99sb = 0
            elseif (line(1:2).eq.'-y') then
               iff99bc = 1
            elseif (line(1:2).eq.'-z') then
               ionpl = 0
            elseif (line(1:2).eq.'-g') then
               if (gargpl('-g',icnt,line,argstr)) then
                  gtol = reada(argstr,1,len(argstr))
               endif
            elseif (line(1:2).eq.'-c') then
               if (gargpl('-c',icnt,line,argstr)) then
                  read(argstr,*) nsd
               endif
            elseif (line(1:2).eq.'-C') then
               usecut = .true.
               usesw  = .true.
            elseif (line(1:2).eq.'-h') then
               call prthlp
               stop
            else 
               if (idebug.eq.1) then
                 write(61,*) 'Unknown commandline option '//line(1:2)
               else
                 print*, 'Unknown commandline option '//line(1:2)
                 call prthlp
               endif
               stop
            endif

         else

            lenf = index(line,' ')-1
            if (linlen(line).gt.lenf) lenf = linlen(line)
            call fndchr(line(1:lenf),lenf,'.',idot)
            call fndchr(line(1:lenf),lenf,'/',isl)
            if (idot.gt.0.and.idot.gt.isl) then
               fniun = line(1:idot-1)
            else
               fniun = line(1:linlen(line))
            endif
            lenf = linlen(fniun)
            if (idebug.eq.1) 
     &         write(61,*) fniun(1:lenf),lenf
            call opfiles
            icnt = n

         endif

         icnt = icnt + 1

      end do

c     no input file found

      if (lenf.eq.-1) then
         print*, 'Missing input file '
         call prthlp
         stop
      endif

      if (idebug.eq.1) then
         write(61,*) "end of argument parsing"
         write(61,*) "gtol=",gtol," nsd=",nsd
         close(61)
      endif

cmpi      call mpi_init(ierr)
cmpi      call mpi_comm_size(mpi_comm_world, nproc, ierr)

      call param

      call parptr(6,fdum,fdum,nproc)
      call parptr(7,fdum,fdum,taskid)

      if (iff99sb.eq.1) call ff99sb
      if (iff99bc.eq.1) call ff99bc

      call getinp(istat)

      if (addbox.and.box) addbox = .false.

      call conn34(istat)

      call bndarr(istat)
      call angarr(istat)
      call itrarr(istat)
      call torarr(istat)
      call asschg(idebug)

      if (idebug.eq.1) then
         call prttor(istat)
         call pritor(istat)
      endif

      if (addbox) then
         box = .true.
         call makbox
         call allbox
         if (nion.ne.0) call cntwat
         call filbox
      endif

      if (cell) then
         call initxp
         dolbfgs   = .false.
      endif

      if (usesw) then
         cutoff = 8.0e0
         call inisw
         call bldlst
      endif

      if (iatoms.gt.333.and.dolbfgs) then
c      if (dolbfgs) then
         call lbfgs(emin,gtol,nsd)
      else
         dolbfgs = .false.
         call optimise(emin,gtol,nsd)
      endif

      if (box) call appbnx()
      if (imon.ne.0) then
         icyc = icyc + 1
         call wrmon(icyc,emin)
      endif

      call wrtout(iun3,emin)
      close(iun3)

      if (iarc.eq.1) close(iun4)
      if (ilog.eq.1) close(iun5)

cmpi      call mpi_finalize(ierr)

      end

      subroutine prthlp
      implicit real (a-h,o-z)

      print*, 
     &           'Usage: ambfor [commandline options] ambfor_file[.xyz]'
      print*, ' '
      print*, '  Commandline options:'
      print*, '    -m     - used in conjuction with molden'
      print*, '             produces a file per iteration:'
      print*, '             ambfor_file.001 etc.'
      print*, '    -M     - same as -m, only now binary'
      print*, '    -n     - write optimisation details to '
      print*, '             ambfor_file.log'
      print*, '    -N     - use cutoff'
      print*, '    -C     - use neighbourlist and switch functions'
      print*, '    -a     - concatenate intermediate structures'
      print*, '             into an archive file : ambfor_file.arc'
      print*, '    -p     - use Powell-Beale conjugate gradients'
      print*, '             (default is LBFGS (atoms > 333)'
      print*, '    -s     - calculate single point energy'
      print*, '    -c 10000 - maximum iterations'
      print*, '    -g 0.01  - gradient convergence criterium'
      print*, '    -b       - generate water filled box'
      print*, '    -x       - use amber ff99 force field '
      print*, '               (default amber ff99sb)'
      print*, '    -y       - use amber ff99bsc0 DNA/RNA force field '
      print*, '    -z       - random placement ions (with -b)'
      print*, '    -h       - this help message'
      print*, ' '

      return
      end

      subroutine prtene(nfun,f,gnorm)
      implicit real (a-h,o-z)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      integer taskid
      common /mpih/ nproc,taskid

      if (nproc.eq.1.or.(nproc.gt.1.and.taskid.eq.0)) then
          write(iun5,'(i6,a,f12.3,a,f12.3)') 
     &      nfun,' Emin=',f,' gnorm=',gnorm
      endif

      return
      end

