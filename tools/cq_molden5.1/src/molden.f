      program molden
      implicit double precision (a-h,o-z)
c
c     ##  COPYRIGHT (C) 1991 by G. Schaftenaar  ##
c
c     MOLecular DENsity Program
c     By G. Schaftenaar
c     CMBI, University of Nijmegen, The Netherlands
c     (formerly CAOS/CAMM Center)
c     1991
c
c
c     NUMATM         MAXIMUM NUMBER OF ATOMS ALLOWED.                   
c     MAXORB         MAXIMUM NUMBER OF ORBITALS ALLOWED.                
c
c
      parameter (numat1=100000)
      parameter (mxcon=10)
      parameter (numatm=100000)
      parameter (maxpt=1000)
      parameter (MAXPNT=2000)
      parameter (maxorb=256)
      parameter (mxcntr=100)
      parameter (maxat=1000)
      parameter (mxel=100)
      parameter (mxres=42)
      parameter (nacc=6)
      parameter (ndon=4)
      parameter (mxplev=5)
      parameter (maxdm=20)
      parameter (mxrss=20)
      parameter (mxatt=10)
      parameter (mxath=8)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxchtp=136)
      parameter (mxvalc=10)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (maxm3d=10)
      parameter (max3d2=max3d*max3d)
      parameter (mxppmf=16)
      parameter (mxlpmf=26)
      parameter (mxpmfl=5)
      parameter (mxmol2=41)
      parameter (mxmm3=164)
      parameter (mxmsf=235)
      parameter (mxamb=1590)
      parameter (mxambc=49)
      parameter (mxgff=72)
      parameter (mxamo=201)
      parameter (mxgmx=53)
      parameter (mxgmx2=57)
      parameter (mxg43=49)
      parameter (mxrsa=23)
      parameter (mxata=9)
      parameter (mxatha=9)
      parameter (mxrsn=19)
      parameter (mxnucl=39)
      parameter (mxhnuc=27)
      parameter (mxrso=24)
      parameter (mxato=10)
      parameter (mxatho=11)
      parameter (mxsg=238)
      parameter (ncmx=32)
      parameter (mxonh=100)
      parameter (mxmmul=100)
      parameter (sincr=1.5d0)
      parameter (aincr=5.0d0)
      parameter (mexcl=10)

***********************************************************************
      logical  bonds, fine, euclid ,ostep, yes, 
     &         molpot,ovrlap,befo,valenc,
     &         domol,statio,forces,odeflt,
     &         docnt,molonl,atomic,doori,ctoz
      logical  denok,bontem,atotem,ovrtem,negfrq,
     &         ssbond,oscal,do3dx,first,dconpt,odos,
     &         elpot,domax,zmat,obin,zsev,wxyz,wmolf

      logical oclose,osshad,ofill,ofrst,oeerst,
     &        oempty,onofil,gargpl,wmovie
      logical keyi,keyr,keyrv,keyirv,
     &        dozme,espchg,rdhb,forpdb,colpsp,monops,
     &        doxmol,opfil,dolap,lapdbl,dosrf2

      logical rdmult,domul,hesend,dovdw,dmachg,
     &        fortnk,doiso,chpot,chkmap,mapxyz

      integer space,do3d
      common /spacom/ grey,iscol,oclose,osshad,ofill,ofrst
      common /pntsta/ iproj,idrcol,nespt
      common /strt/   oeerst
      common /plot/   iplot,iplwin,icolps
      common /plspec/ step,stepf,scale,scalf,cut,pmin,pmax,one,
     &                oscal,ostep,fine
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /atopla/ xsym(numatm),ysym(numatm),zsym(numatm),
     &                isym(numatm)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      common /densty/ p(maxorb*maxorb), paa(maxorb)
      real stoalfa,stobnorm
      common /ADF/    istos(5,maxorb),naorbs,stoalfa(maxorb),
     &                stobnorm(maxorb)
      common /sphere/ averag(maxorb)
      common /dynhlp/ qd(maxorb*maxorb),pd(maxorb),gd(3,maxorb),
     &                hd(6,maxorb)
      common /orihlp/ ori(numatm),iuser(numatm),oalpha(numatm),
     &                obeta(numatm),norien,ibal
      common /spa3d/  denn(max3d*max3d*max3d),pmnn(max3d)
      common /grdhlp/ mx3d,mx3d2
      common /cnthlp/ rng1,rng2,vlcnt,vlcnt2
      common /map3d/  fmap(maxm3d*maxm3d*maxm3d),mxm3d
      common /mapcls/ valcol(5),colmap(3,5)
      character*80 mapfil
      character*80 grdfil
      common /maphlp/ mapfil,grdfil
      common /wrtcom/ zsev
      common /animo/  movie

      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)
      common /keywrd/ keywrd,keyori
      common /coord / xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /pseudo/ ipseud,ivale(numatm)
      common /oniomh/ ioni,nion,ionih(mxonh),natonh(mxonh),
     &                iomap(numatm),icntat(mxonh),fct(mxonh),
     &                xyzi(3,mxonh)
      common /typoni/ ioniad
      common /orbhlp/ mxorb,iuhf,ispd
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /totchg/ itot
      common /chgbck/ qat(numat1)
      common /extchg/ exchg(3,3),iexchg(3),nexchg
      common /qmchar/ qch(numat1),ihasesp
      common /cntval/ cntval, euclid, yes, icnt, iflag
      common /esprad/ vander(mxel)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /sclcom/ scal,fscal,scali,smag,iscupd
      common /pnthlp/ ipoints,ipnt
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /athlp/ iatoms, mxnat
      integer*2 ityp,ipdbt
      common /types/ iff
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /chmtyp/ ncc(mxrss),icc(2,mxatt,mxrss),
     &                nhh(mxrss),ihh(2,mxath,mxrss)
      integer amboff,ntca,ctca
      common /ambtyp/ amboff(mxrss*3),ntca(mxrss),ctca(mxrss),
     &                ncca(mxrsa),icca(2,mxata,mxrsa),
     &                nhha(mxrsa),ihha(2,mxatha,mxrsa),
     &                ncco(mxrso),icco(2,mxato,mxrso),
     &                nhho(mxrso),ihho(2,mxatho,mxrso),
     &                nnuc(mxrsn),nhnuc(mxrsn),irna(mxres-23),
     &                inuc(2,mxnucl,mxrsn),ihnuc(2,mxhnuc,mxrsn)

      character*2 ppmf, lpmf, gffstr
      character*5 mol2
      character*19 mm3
      character*20 chmtnk
      character*4 chmsf
      character*20 ambstr
      character*20 amostr
      character*28 gffext(mxgff)
      common /ftypes/ihasl(11),mol2(mxmol2),mm3(mxmm3),chmtnk(mxchtp),
     &               chmsf(mxmsf),ambstr(mxamb),amostr(mxamo),
     &               ppmf(mxppmf),lpmf(mxlpmf)
      character*3 pdbsym,hsym,chtnk,ambtnk
      character*2 amotnk
      common /symbol/ pdbsym(mxsym),hsym(mxhsym),chtnk(mxchtp),
     &                ambtnk(mxamb),amotnk(mxamo),gffstr(mxgff)
      character*5 gro43,grogmx,grog2x,gro43l,grogmxl,grog2xl
      character*35 gro43s, grogms, grog2s
      common /symgro/ gro43(mxg43),grogmx(mxgmx),grog2x(mxgmx2),
     &                gro43l(mxg43),grogmxl(mxgmx),grog2xl(mxgmx2),
     &                gro43s(mxg43),grogms(mxgmx),grog2s(mxgmx2)
      integer ambvdt
      common /fcharg/ ambchg(mxamb),ambvw1(mxambc),ambvw2(mxambc),
     &                gfvdw(2,mxgff),ambvdt(mxamb),cysneg(9)
      common /zmatin/ bl(maxat),alph(maxat),bet(maxat),
     &                ibl(maxat),ialph(maxat),ibet(maxat),
     &                ianzz(maxat),iz(4,maxat),imap(maxat)
      common /zmatst/ nwrit,nvar,isimpl
      common /zmsbck/ lwrit(numat1)
      common /zmfrst/ ihaszm, nz, mxzat
      common /zmtyp/  igztyp
      common /zmpart/ ipart,imn,imx,idcur
      common /tstoc/  cstoc(3,maxat),czstoc(3,maxat),
     &                astoc(maxat),bstoc(maxat),ianstc(maxat)
      common /hring/  lring(numat1)
      common /selatm/ jring(numat1)
      common /pltmp/  inat(numat1)
      common /gammus/ iold
      common /mopver/ mopopt
      common /cllab/  iclon,iclpnt(4)
      character*2 clstr
      common /clchr/  clstr(4)
      integer*2 ir,it
      common /cell/   natc,ndum(6),ichx,icrtp,nopr,
     &                ir(3,3,192),it(3,192),rdum(12)

      common /pers/xv,yv,zv,c0,pincr
      real eiga,eigb
      common /eigval/ eiga(maxorb),eigb(maxorb)
      common /nmr/    shlnuc(numatm),ihsnmr
      common /uvspec/ ihasex

      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      integer reson
      character*1 achain
      common /calf/ rphi(numcal),rpsi(numcal),
     &              icalf(6,numcal),ncalf,ianf(mxchai),islu(mxchai),
     &              nchain,iamino(numcal),ihet(mxheta),
     &              iclhet(mxheta),reson(numcal),issdon,
     &              icxp(numcal),icyp(numcal),
     &              iams(numcal),ihets(mxheta),ibck(4),
     &              isal(numcal),irsnr(numcal),ihashb,ishoh,
     &              achain(numcal)
      character*3 hetz
      common /clfstr/ ihashz,ihetq(mxheta),ihqset(mxheta),ihhadd(mxheta)
     &                ,hetz(mxheta)
      common /clfhlp/ isndcl(4),iamicl(mxres),ichcol(mxchai)
      character*3 aminos
      common /amino/aminos(mxres)
      parameter (mxaln=20)
      common /palign/ nalign,istst(mxaln),istch(mxaln),istres(mxaln),
     &                istcol(mxaln),istptr(mxaln)
      common /helpar/ helrod, ihtype
      logical finded
      integer srfmap, srfloc
      common /hlpsrf/ npts1,npts2,npts3,srfmap,srfloc,ifogl,itsrf
      common /psa/ psa, tsa, exs, pol, pol2, epmin, epmax, ipsa,
     &             icpsa, idtpsa

      character*2 elemnt
      common /elem/elemnt(mxel)
 
      common /savcor/ angs(3,maxpt),idon(maxpt)
      common /eulx/   ca,cb,sa,sb,cc,sc

c iftyp      file type
c
c 1          mopac
c 2          gamess
c 3          gamus
c 4          gauss
c 5          adfin
c 6          chemx
c 7          cpmd
c 8          qchem
c 9          orca
c 10         xyz
c 11         tinker + small
c 12         tinker protien
c 13         mol file
c 15         nwchem

      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character fniun*256,filecpmd*256,vrmlfil*256
      common /fnunit/ fniun
      common /vrunit/ vrmlfil
 
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common /lapcom/ lapbox,lapdbl
      common /surf/   natorg,nosncd
      common /srfbck/ isurf(numat1)
      real frmul
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifreq,normc
      integer dolabs,fancy,persp,shade,atcol,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo
      common /potlev/grdw, plevel(mxplev),ipcol(mxplev+1),nplev
      common /srfhlp/edge,ctval(mxvalc),pxyz(3),nvalc,nspts,istyp
      common /coars/ icoars
      common /backg/ ibgcol,ibgclo,ibgmod
      logical hon,notitle,extest,dompi
      real hamin,hamax
      common /hcon/  iacc(nacc),idn(ndon),hdmin,hdmax,hamin,hamax,hon
      common /vropt/ ivtwo,ihand,ivadd
      common /xyzopt/ ixyz,ipdbwh,iambch,nwramb
      character*80 glin1,glin2,gtitl,rungam
      character*15 jname,qname
      common /gauopt/ ito,imo,ibo,itotc,imult,ibatch,ihess,itime,
     &                iwxyz,ichh,ichm,ichl,imh,imm,iml,
     &                glin1,glin2,gtitl,jname,qname,rungam
      character*80 tnknm
      integer tnkbg,tnkit,tnkarc,tnkarf,tnkprg
      common /tnkopt/ rmsgrd,tnkbg,tnkit,tnkarc,tnkarf,tnkprg,tnknm,icst
      common /tnkpro/ iresrd

      common /getpnt/ irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &                iconv,ircus,dozme
      common /distmn/ rdm(maxdm),idmon(2,maxdm),ndm
      common /cllmat/ rrr(18),itz
      common /gracom/ uscl,colscd,colscpd,ivdwpl
      common /gtmfil/ igtfil
      common /strips/ qnormo(3),crpnto(3,ncmx),crnrmo(3,ncmx),
     &                numcir,nquad
      common /align/  vecs(3,3),nscnd,iscst,ialtyp,iocnt
      character*80 spcdfil
      character*7 spnam
      common /spgrnm/ spnam(mxsg)
      common /vrmlhl/ spcdfil
      common /pmflvl/ ipmfm,ipmfh,pmflev(mxpmfl),levcol(mxpmfl)
      common /comsrf/ vo(3), rt(3),v1t(3),v2t(3),v3t(3),wo(3),
     &                sl(3),isl
      common /doh/ idoh
      common /metexc/ qexcl(mexcl),ianexc(mexcl)
c     with some versions of gfortran it helps to uncomment this line:
c      external iargc
c GCC 4.0 or higher: 
c      external gfortran_iargc
      dimension fcnt(mxcntr),tmp(3)
      character title*80, keywrd*320, keyhlp*80, keyori*320
      character tmpstr*80,movfil*80,povfil*80,oglfil*80
      character esc, etx, lf, cr, eot, gs, ff, fs, us
      character*2 gstr
      character*5 hstr
      character line*80,liris*256
      character plfile*75, parfil*75, argstr*75, basfile*75
      character hom*500
      character lstr*137
      character str*100
      common /rotmat/rx(3),ry(3),rzz(3),trns(3),td(3)
      common /setogl/ idirogl
      common /denfir/ first
      common /multim/ imulm, nmulm,ihasqm(mxmmul)
      common /pbc/ abc(3),ibox,icell,igfmap
      real xx, yy,coloff

      data isndcl/ 3,9,10,8 /
      data iamicl /3,3,4,6,4,3,3,6,1,9,3,7,1,9,3,7,7,3,4,3,
     &             7,10,11,13,5,7,10,10,13,5,5,7,7,7,7,7,7,1,10,
     &             10,10,10/
      data ichcol/12,2,7,5,6,10,11,1,3,4,8,9,13,14,36*15/

      data aminos/'GLY','ALA','SER','CYS','THR','ILE','VAL','MET',
     &            'ASP','ASN','LEU','LYS','GLU','GLN','PRO','ARG',
     &            'HIS','PHE','TYR','TRP','ASX','GLX','HYP','  A',
     &            '  C','  G','  T','  U','1MA','5MC','OMC','1MG',
     &            '2MG','M2G','7MG','OMG',' YG','  I',' +U','H2U',
     &            '5MU','PSU'/

      data pdbsym/
     &   'N  ','CA ','C  ','O  ','CB ','CG ','CG1','CG2','CD ','CD1',
     &   'CD2','CE ','CE1','CE2','CE3','CH2','CZ ','CZ2','CZ3','ND1',
     &   'ND2','NE ','NE1','NE2','NH1','NH2','NZ ','OD ','OD1','OD2',
     &   'OG ','OG1','OH ','OE1','OE2','SD ','SG ','OXT','AD1','AD2',
     &   'AE1','AE2','P  ','O1P','O2P','O5*','C5*','C4*','O4*','C3*',
     &   'O3*','C2*','O2*','C1*','C2 ','C4 ','C5 ','C6 ','C8 ','N1 ',
     &   'N2 ','N3 ','N4 ','N6 ','N7 ','N9 ','O2 ','O4 ','O6 ','C5M',
     &   'C1 ','C5A','C2A','C2B','C7 ','O3P','O5T','O3T','CM1','CM2',
     &   'CM3','CM4','CM5','CM6','CM7','CM8','CM9','C3 ','C10','C11',
     &   'C12','C13','C14','C15','C16','O17','O18','C19','N20','C21',
     &   'O22','O23','C24'/
      data hsym/
     &   'H  ','HA ','HB ','HG ','HG1','HG2','HD ','HD1','HD2','HE ',
     &   'HE1','HE2','HE3','HZ ','HZ1','HZ2','HZ3','HH ','HH1','HH2',
     &   'H5*','H4*','H3*','H2*','H1*','HO2','H1 ','H2 ','H21','H22',
     &   'H3 ','H41','H42','H5 ','H6 ','H61','H62','H8 ','HO3','H5T',
     &   'HN1','HN2','HN3','HN4','HN5','HN6','HN7','HN8','HN9','HM1',
     &   'HM2','HM3','HM4','HM5','HM6','HM7','HM8','HM9','H10','H13',
     &   'H14','H15','H19','H24'/

      data elemnt/' H','He',
     2 'Li','Be',' B',' C',' N',' O',' F','Ne',
     3 'Na','Mg','Al','Si',' P',' S','Cl','Ar',
     4 ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu',
     4 'Zn','Ga','Ge','As','Se','Br','Kr',
     5 'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag',
     5 'Cd','In','Sn','Sb','Te',' I','Xe',
     6 'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',
     6 'Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re','Os','Ir','Pt',
     6 'Au','Hg','Tl','Pb','Bi','Po','At','Rn',
     7 'Fr','Ra','Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf',
     8 'XX','  '/

      data vrad/0.200d0,0.286d0,0.340d0,0.589d0,0.415d0,0.400d0,0.400d0,
     1          0.400d0,0.320d0,0.423d0,0.485d0,0.550d0,0.675d0,0.600d0,
     2          0.525d0,0.510d0,0.495d0,0.508d0,0.665d0,0.495d0,0.720d0,
     3          0.735d0,0.665d0,0.675d0,0.675d0,0.670d0,0.615d0,0.750d0,
     4          0.760d0,0.725d0,0.610d0,0.585d0,0.605d0,0.610d0,0.605d0,
     5          0.524d0,0.735d0,0.560d0,0.890d0,0.780d0,0.740d0,0.735d0,
     6          0.675d0,0.700d0,0.725d0,0.750d0,0.795d0,0.845d0,0.815d0,
     7          0.730d0,0.730d0,0.735d0,0.700d0,0.577d0,0.835d0,0.670d0,
     8          0.935d0,0.915d0,0.910d0,0.905d0,0.900d0,0.900d0,0.995d0,
     9          0.895d0,0.880d0,0.875d0,0.870d0,0.865d0,0.860d0,0.970d0,
     &          0.860d0,0.785d0,0.715d0,0.685d0,0.675d0,0.685d0,0.660d0,
     &          0.750d0,0.750d0,0.850d0,0.775d0,0.770d0,0.770d0,0.840d0,
     &          1.000d0,1.000d0,1.000d0,0.950d0,0.940d0,0.895d0,0.805d0,
     &          0.790d0,0.775d0,1.000d0,0.755d0,1.000d0,1.000d0,0.765d0,
     &          0.100d0,0.900d0/

c de laatste waarde van invloed op helix solid

c CSD atomic radii + .20

      data vdwr/0.430d0,0.741d0,0.880d0,0.550d0,1.030d0,0.900d0,0.880d0,
     1          0.880d0,0.840d0,0.815d0,1.170d0,1.300d0,1.550d0,1.400d0,
     2          1.250d0,1.220d0,1.190d0,0.995d0,1.530d0,1.190d0,1.640d0,
     3          1.670d0,1.530d0,1.550d0,1.550d0,1.540d0,1.530d0,1.700d0,
     4          1.720d0,1.650d0,1.420d0,1.370d0,1.410d0,1.420d0,1.410d0,
     5          1.069d0,1.670d0,1.320d0,1.980d0,1.760d0,1.680d0,1.670d0,
     6          1.550d0,1.600d0,1.650d0,1.700d0,1.790d0,1.890d0,1.830d0,
     7          1.660d0,1.660d0,1.670d0,1.600d0,1.750d0,1.870d0,1.540d0,
     8          2.070d0,2.030d0,2.020d0,2.010d0,2.000d0,2.000d0,2.190d0,
     9          1.990d0,1.960d0,1.950d0,1.940d0,1.930d0,1.920d0,2.140d0,
     &          1.920d0,1.770d0,1.630d0,1.570d0,1.550d0,1.570d0,1.520d0,
     &          1.700d0,1.700d0,1.900d0,1.750d0,1.740d0,1.740d0,1.880d0,
     &          0.200d0,0.200d0,0.200d0,2.100d0,2.080d0,1.990d0,1.810d0,
     &          1.780d0,1.750d0,0.200d0,1.710d0,0.200d0,0.200d0,1.730d0,
     &          0.100d0,0.200d0/

      data vander/1.20d0,1.20d0,1.37d0,1.45d0,1.45d0,1.50d0,1.50d0,
     &            1.40d0,1.35d0,1.30d0,1.57d0,1.36d0,1.24d0,1.17d0,
     &            1.80d0,1.75d0,1.70d0,17*2.5d0,2.3d0,65*2.5d0/

c de laatste waarde van invloed op helix solid
      data icol /15,4,
     1        8,8,9,14,7,1,2,12,
     2        8,8,7,11,2,6,3,12,
     3        2,1,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,12,
     4        2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,
     5        2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
     3        2,2,2,2,2,2,2,2,2,2,2,2,2,4,2,2,2,2,2,2,2,2,2,
     4        2,2,2,1,10/

      data clstr/' O',' A',' B',' C'/
      data iedlog/max3d2*1/
      data iacc /7,8,9,17,35,53/
      data idn /7,8,9,17/

      data ncc /1,1,2,2,3,4,3,4,4,4,4,5,5,5,5,7,6,7,8,10/
      data nhh /0,0,1,1,5,8,6,3,0,2,6,7,0,2,2,7,4,5,5,6/
c GLY
      data ((icc(i,j,1),i=1,2),j=1,mxatt) /
     & 2,28,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,1),i=1,2),j=1,mxath) /
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c ALA
      data ((icc(i,j,2),i=1,2),j=1,mxatt) /
     & 5,27,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,2),i=1,2),j=1,mxath) /
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c SER
      data ((icc(i,j,3),i=1,2),j=1,mxatt) /
     & 5,35,31,76,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,3),i=1,2),j=1,mxath) /
     & 10,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c CYS
      data ((icc(i,j,4),i=1,2),j=1,mxatt) /
     & 5,37,37,80,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,4),i=1,2),j=1,mxath) /
     & 10,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c THR
      data ((icc(i,j,5),i=1,2),j=1,mxatt) /
     & 5,36,8,27,32,76,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,5),i=1,2),j=1,mxath) /
     & 10,8,13,8,16,1,17,1,18,1,0,0,0,0,0,0/
c ILE
      data ((icc(i,j,6),i=1,2),j=1,mxatt) /
     & 5,25,7,26,8,27,10,27,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,6),i=1,2),j=1,mxath) /
     & 13,1,14,1,16,1,17,1,18,1,22,1,23,1,24,1/
c VAL
      data ((icc(i,j,7),i=1,2),j=1,mxatt) /
     & 5,25,7,27,8,27,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,7),i=1,2),j=1,mxath) /
     & 13,1,14,1,15,1,16,1,17,1,18,1,0,0,0,0/
c MET
      data ((icc(i,j,8),i=1,2),j=1,mxatt) /
     & 5,26,6,56,12,57,36,81,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,8),i=1,2),j=1,mxath) /
     & 28,1,29,1,30,1,0,0,0,0,0,0,0,0,0,0/
c ASP
      data ((icc(i,j,9),i=1,2),j=1,mxatt) /
     & 5,54,6,55,29,78,30,78,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,9),i=1,2),j=1,mxath) /
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c ASN
      data ((icc(i,j,10),i=1,2),j=1,mxatt) /
     & 5,26,6,49,29,75,21,64,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,10),i=1,2),j=1,mxath) /
     & 25,3,26,3,0,0,0,0,0,0,0,0,0,0,0,0/
c LEU
      data ((icc(i,j,11),i=1,2),j=1,mxatt) /
     & 5,26,6,25,10,27,11,27,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,11),i=1,2),j=1,mxath) /
     & 22,1,23,1,24,1,25,1,26,1,27,1,0,0,0,0/
c LYS
      data ((icc(i,j,12),i=1,2),j=1,mxatt) /
     & 5,26,6,26,9,26,12,58,27,65,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,12),i=1,2),j=1,mxath) /
     & 19,1,20,1,28,17,29,17,40,6,41,6,42,6,0,0/
c GLU
      data ((icc(i,j,13),i=1,2),j=1,mxatt) /
     & 5,26,6,54,9,55,34,78,35,78,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,13),i=1,2),j=1,mxath) /
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c GLN
      data ((icc(i,j,14),i=1,2),j=1,mxatt) /
     & 5,26,6,26,9,49,34,75,24,64,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,14),i=1,2),j=1,mxath) /
     & 34,3,35,3,0,0,0,0,0,0,0,0,0,0,0,0/
c PRO
      data ((icc(i,j,15),i=1,2),j=1,mxatt) /
     & 1,66,2,30,5,31,6,31,9,32,0,0,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,15),i=1,2),j=1,mxath) /
     & 19,1,20,1,0,0,0,0,0,0,0,0,0,0,0,0/
c ARG
      data ((icc(i,j,16),i=1,2),j=1,mxatt) /
     & 5,26,6,26,9,59,17,60,22,72,25,73,26,73,0,0,0,0,0,0/
      data ((ihh(i,j,16),i=1,2),j=1,mxath) /
     & 19,1,20,1,28,18,55,19,56,19,58,19,59,19,0,0/
c HIS
      data ((icc(i,j,17),i=1,2),j=1,mxatt) /
     & 5,39,6,42,11,42,13,45,20,70,24,70,0,0,0,0,0,0,0,0/
      data ((ihh(i,j,17),i=1,2),j=1,mxath) /
     & 22,10,25,12,31,13,34,10,0,0,0,0,0,0,0,0/
c PHE
      data ((icc(i,j,18),i=1,2),j=1,mxatt) /
     & 5,26,6,47,10,21,11,21,13,21,14,21,17,21,0,0,0,0,0,0/
      data ((ihh(i,j,18),i=1,2),j=1,mxath) /
     & 22,2,25,2,31,2,34,2,40,2,0,0,0,0,0,0/
c TYR
      data ((icc(i,j,19),i=1,2),j=1,mxatt) /
     & 5,26,6,47,10,21,11,21,13,21,14,21,17,48,33,77,0,0,0,0/
      data ((ihh(i,j,19),i=1,2),j=1,mxath) /
     & 22,2,25,2,31,2,34,2,52,8,0,0,0,0,0,0/
c TRP
      data ((icc(i,j,20),i=1,2),j=1,mxatt) /
     & 5,26,6,50,10,51,11,52,14,53,15,21,16,21,18,21,19,21,23,71/
      data ((ihh(i,j,20),i=1,2),j=1,mxath) /
     & 22,2,31,9,37,2,46,2,49,2,58,2,0,0,0,0/

      data chtnk /
     & 'HA ','HP ','H  ','HB ','HB ','HC ','HC ','H  ','H  ','H  ',
     & 'H  ','HR1','HR2','HR3','HR3','HS ','HA ','HC ','HC ','C  ',
     & 'CA ','CT1','CT1','CC ','CT1','CT2','CT3','CT2','CT2','CP1',
     & 'CP2','CP3','CP1','CP3','CT2','CT1','CT2','CT2','CT2','CT2',
     & 'CT2','CH1','CH1','CH1','CH2','CH2','CA ','CA ','CC ','CY ',
     & 'CA ','CPT','CPT','CT2','CC ','CT2','CT3','CT2','CT2','C  ',
     & 'CT3','CT ','NH1','NH2','NH3','N  ','NP ','NR1','NR2','NR3',
     & 'NY ','NC2','NC2','O  ','O  ','OH1','OH1','OC ','OC ','S  ',
     & 'S  ','SM ','HA ','HA ','HA ','HA1','HA2','HT ','H  ','CD ',
     & 'CT2','CPA','CPB','CPM','C  ','CM ','CS ','CE1','CE2','NPH',
     & 'OT ','OH1','OB ','OS ','OM ','SS ','FE ','HE ','NE ','SOD',
     & 'MG ','CLA','POT','CAL','ZN ','HL ','HL1','HL2','HL3','HEL',
     & 'CL ','CL1','CL2','CL2','CL2','CL3','CL2','CL2','CL5','CEL',
     & 'NTL','OSL','OBL','OSL','O2L','PL '/

      data amboff /1,7,55,77,65,41,15,258,210,220,
     &             27,271,232,244,96,287,161,107,122,138,
CNF N-Term
     &           350,356,380,392,386,374,362,472,448,454,
     &           368,478,460,466,404,484,430,412,418,424,
CNF C-Term
     &           501,507,531,543,537,525,513,620,596,602,
     &           519,626,608,614,555,632,578,560,566,572/
c cystine -SS- 83, HID 174, HIE 190, ORN 300, mALA 314, pyro 321 
c RES N-TERM 352, RES C-TERM 394, GLY N-TERM 346, GLY C-TERM 388
c PRO N-TERM 365, PRO C-TERM 407 (No CB,HB)
c PRO N-TERM CD 371, PRO N-TERM HD 372, PRO C-TERM CD 412,
c PRO C-TERM HD 413

c N CA C O CB (1,2,3,4,5) (gly, NTERM,CTERM; geen 5)
c HN HA HB (1-3,4-5,7-9) (gly, NTERM,CTERM; geen 7-9)
c PRO heeft een geen HN, dus de rest verschuifd

      data ntca /347,353,361,363,362,360,358,383,379,380,359,384,
     &           381,382,366,385,376,373,374,375/

c cystine -SS- 364, HID 377, HIE 378, ORN 386, AIB 387

      data ctca /389,395,403,405,404,402,400,424,420,421,401,425,
     &           422,423,408,426,417,414,415,416/
c cystine -SS- 406, HID 418, HIE 419, ORN 427, AIB 428

      data ncca /0,0,1,1,2,3,2,3,3,3,3,4,4,4,2,6,5,6,7,9,1,5,5/
      data nhha /0,0,1,1,5,8,6,5,0,2,7,9,2,4,8,9,4,5,5,6,0,3,3/

c GLY
      data ((icca(i,j,1),i=1,2),j=1,mxata) /
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,1),i=1,2),j=1,mxatha) /
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c ALA
      data ((icca(i,j,2),i=1,2),j=1,mxata) /
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,2),i=1,2),j=1,mxatha) /
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c SER
      data ((icca(i,j,3),i=1,2),j=1,mxata) /
     & 31,63,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,3),i=1,2),j=1,mxatha) /
     & 10,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c CYS -SH
      data ((icca(i,j,4),i=1,2),j=1,mxata) /
     & 37,85,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,4),i=1,2),j=1,mxatha) /
     & 10,86,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c THR
      data ((icca(i,j,5),i=1,2),j=1,mxata) /
     & 8,75,32,73,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,5),i=1,2),j=1,mxatha) /
     & 10,74,13,74,16,76,17,76,18,76,0,0,0,0,0,0,0,0/
c ILE
      data ((icca(i,j,6),i=1,2),j=1,mxata) /
     & 7,49,8,51,10,53,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,6),i=1,2),j=1,mxatha) /
     & 13,50,14,50,16,52,17,52,18,52,22,54,23,54,24,54,0,0/
c VAL
      data ((icca(i,j,7),i=1,2),j=1,mxata) /
     & 7,23,8,25,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,7),i=1,2),j=1,mxatha) /
     & 13,24,14,24,15,24,16,26,17,26,18,26,0,0,0,0,0,0/
c MET
      data ((icca(i,j,8),i=1,2),j=1,mxata) /
     & 6,266,12,269,36,268,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,8),i=1,2),j=1,mxatha) /
     & 10,267,11,267,28,270,29,270,30,270,0,0,0,0,0,0,0,0/
c ASP
      data ((icca(i,j,9),i=1,2),j=1,mxata) /
     & 6,218,29,219,30,219,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,9),i=1,2),j=1,mxatha) /
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c ASN
      data ((icca(i,j,10),i=1,2),j=1,mxata) /
     & 6,228,29,229,21,230,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,10),i=1,2),j=1,mxatha) /
     & 25,231,26,231,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c LEU
      data ((icca(i,j,11),i=1,2),j=1,mxata) /
     & 6,35,10,37,11,39,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,11),i=1,2),j=1,mxatha) /
     & 10,36,22,38,23,38,24,38,25,40,26,40,27,40,0,0,0,0/
c LYS
      data ((icca(i,j,12),i=1,2),j=1,mxata) /
     & 6,279,9,281,12,283,27,285,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,12),i=1,2),j=1,mxatha) /
     & 10,280,11,280,19,282,20,282,28,284,29,284,
     & 40,286,41,286,42,286/
c GLU
      data ((icca(i,j,13),i=1,2),j=1,mxata) /
     & 6,240,9,242,34,243,35,243,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,13),i=1,2),j=1,mxatha) /
     & 10,241,11,241,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c GLN
      data ((icca(i,j,14),i=1,2),j=1,mxata) /
     & 6,252,9,254,34,255,24,256,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,14),i=1,2),j=1,mxatha) /
     & 10,253,11,253,34,257,35,257,0,0,0,0,0,0,0,0,0,0/
c PRO
      data ((icca(i,j,15),i=1,2),j=1,mxata) /
     & 6,103,9,105,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,15),i=1,2),j=1,mxatha) /
     & 10,104,11,104,13,104,16,104,19,106,20,106,22,106,25,106,0,0/
c ARG
      data ((icca(i,j,16),i=1,2),j=1,mxata) /
     & 6,295,9,297,17,301,22,299,25,302,26,302,0,0,0,0,0,0/
      data ((ihha(i,j,16),i=1,2),j=1,mxatha) /
     & 10,296,11,296,19,298,20,298,28,300,
     & 55,303,56,303,58,303,59,303/
c HIS (+, HIP)
      data ((icca(i,j,17),i=1,2),j=1,mxata) /
     & 6,169,11,172,13,174,20,170,24,176,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,17),i=1,2),j=1,mxatha) /
     & 22,171,25,173,31,175,34,177,0,0,0,0,0,0,0,0,0,0/
c PHE
      data ((icca(i,j,18),i=1,2),j=1,mxata) /
     & 6,115,10,116,11,116,13,118,14,118,17,120,0,0,0,0,0,0/
      data ((ihha(i,j,18),i=1,2),j=1,mxatha) /
     & 22,117,25,117,31,119,34,119,40,121,0,0,0,0,0,0,0,0/
c TYR
      data ((icca(i,j,19),i=1,2),j=1,mxata) /
     & 6,130,10,131,11,131,13,133,14,133,17,135,33,136,0,0,0,0/
      data ((ihha(i,j,19),i=1,2),j=1,mxatha) /
     & 22,132,25,132,31,134,34,134,52,137,0,0,0,0,0,0,0,0/
c TRP
      data ((icca(i,j,20),i=1,2),j=1,mxata) /
     & 6,146,10,147,11,149,14,152,15,153,16,159,
     & 18,155,19,157,23,150/
      data ((ihha(i,j,20),i=1,2),j=1,mxatha) /
     & 22,148,31,151,37,154,46,156,49,158,58,160,0,0,0,0,0,0/
c CYS -SS-
      data ((icca(i,j,21),i=1,2),j=1,mxata) /
     & 37,95,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,21),i=1,2),j=1,mxatha) /
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c HISD
      data ((icca(i,j,22),i=1,2),j=1,mxata) /
     & 6,186,11,189,13,191,20,187,24,193,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,22),i=1,2),j=1,mxatha) /
     & 22,188,25,190,31,192,0,0,0,0,0,0,0,0,0,0,0,0/
c HISE
      data ((icca(i,j,23),i=1,2),j=1,mxata) /
     & 6,202,11,204,13,206,20,203,24,208,0,0,0,0,0,0,0,0/
      data ((ihha(i,j,23),i=1,2),j=1,mxatha) /
     & 25,205,31,207,34,209,0,0,0,0,0,0,0,0,0,0,0,0/


      data (ambtnk(i),i=1,100) /
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
      data (ambtnk(i),i=101,200) /
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
      data (ambtnk(i),i=201,300) /
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
      data (ambtnk(i),i=301,400) /
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
      data (ambtnk(i),i=401,500) /
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
      data (ambtnk(i),i=501,600) /
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
      data (ambtnk(i),i=601,659) /
     & 'H1 ','N  ','CT ','C  ','H  ','O2 ','H1 ','N  ','CT ','C  ',
     & 'H  ','O2 ','H1 ','N  ','CT ','C  ','H  ','O2 ','H1 ','N  ',
     & 'CT ','C  ','H  ','O2 ','H1 ','N  ','CT ','C  ','H  ','O2 ',
     & 'H1 ','N  ','CT ','C  ','H  ','O2 ','H1 ','N  ','CT ','C  ',
     & 'H  ','O2 ','H1 ','N  ','CT ','C  ','H  ','O2 ','OW ','HW ',
     & 'Li+','Na+','K+ ','Rb+','Cs+','Mg+','Ca+','Zn+','Cl-'/

c neutral asp,glu,lys
      data (ambtnk(i),i=660,671) /
     & 'N  ','CT ','C  ','H  ','O  ','H1 ','CT ','HC ',
     & 'C  ','O  ','OH ','HO '/ 
      data (ambtnk(i),i=672,685) /
     & 'N  ','CT ','C  ','H  ','O  ','H1 ','CT ','HC ',
     & 'CT ','HC ','C  ','O  ','OH ','HO '/ 
      data (ambtnk(i),i=686,701) /
     & 'N  ','CT ','C  ','H  ','O  ','H1 ','CT ','HC ',
     & 'CT ','HC ','CT ','HC ','CT ','HP ','N3 ','H  '/ 

      data (ambtnk(i),i=702,1000) / 299*"   "/

      data (ambtnk(i),i=1001,1100) /
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

      data (ambtnk(i),i=1101,1200) /
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

      data (ambtnk(i),i=1201,1253) /
     & 'HA ','H4 ','OS ','CT ','H1 ','H1 ','CT ','H1 ','OS ','CT ',
     & 'H2 ','CT ','H1 ','CT ','HC ','HC ','OS ','N* ','C  ','NA ',
     & 'C  ','CM ','CM ','O  ','H  ','O  ','CT ','HC ','H4 ','P  ',
     & 'O2 ','OH ','HO ','OS ','P  ','O2 ','OH ','HO ','OS ','P  ',
     & 'O2 ','P  ','O2 ','OH ','HO ','OS ','P  ','O2 ','OH ','HO ',
     & 'OS ','P  ','O2 '/

      data (ambchg(i),i=1,60) /
     &  -0.4157 , -0.0252 ,  0.5973 ,  0.2719 , -0.5679 ,  0.0698
     &, -0.4157 ,  0.0337 ,  0.5973 ,  0.2719 , -0.5679 ,  0.0823
     &, -0.1825 ,  0.0603 , -0.4157 , -0.0875 ,  0.5973 ,  0.2719
     &, -0.5679 ,  0.0969 ,  0.2985 , -0.0297 , -0.3192 ,  0.0791
     &, -0.3192 ,  0.0791 , -0.4157 , -0.0518 ,  0.5973 ,  0.2719
     &, -0.5679 ,  0.0922 , -0.1102 ,  0.0457 ,  0.3531 , -0.0361
     &, -0.4121 ,  0.1000 , -0.4121 ,  0.1000 , -0.4157 , -0.0597
     &,  0.5973 ,  0.2719 , -0.5679 ,  0.0869 ,  0.1303 ,  0.0187
     &, -0.0430 ,  0.0236 , -0.3204 ,  0.0882 , -0.0660 ,  0.0186
     &, -0.4157 , -0.0249 ,  0.5973 ,  0.2719 , -0.5679 ,  0.0843/

      data (ambchg(i),i=61,120) /
     &   0.2117 ,  0.0352 , -0.6546 ,  0.4275 , -0.4157 , -0.0389
     &,  0.5973 ,  0.2719 , -0.5679 ,  0.1007 ,  0.3654 ,  0.0043
     &, -0.6761 ,  0.4102 , -0.2438 ,  0.0642 , -0.4157 ,  0.0213
     &,  0.5973 ,  0.2719 , -0.5679 ,  0.1124 , -0.1231 ,  0.1112
     &, -0.3119 ,  0.1933 , -0.4157 ,  0.0429 ,  0.5973 ,  0.2719
     &, -0.5679 ,  0.0766 , -0.0790 ,  0.0910 , -0.1081 , -0.2548
     &, -0.0266 ,  0.5896 , -0.5748 ,  0.0641 , -0.0070 ,  0.0253
     &,  0.0189 ,  0.0213 ,  0.0192 ,  0.0391 , -0.4157 , -0.0024
     &,  0.5973 ,  0.2719 , -0.5679 ,  0.0978 , -0.0343 ,  0.0295
     &,  0.0118 , -0.1256 ,  0.1330 , -0.1704 ,  0.1430 , -0.1072/

      data (ambchg(i),i=121,180) /
     &   0.1297 , -0.4157 , -0.0014 ,  0.5973 ,  0.2719 , -0.5679
     &,  0.0876 , -0.0152 ,  0.0295 , -0.0011 , -0.1906 ,  0.1699
     &, -0.2341 ,  0.1656 ,  0.3226 , -0.5579 ,  0.3992 , -0.4157
     &, -0.0275 ,  0.5973 ,  0.2719 , -0.5679 ,  0.1123 , -0.0050
     &,  0.0339 , -0.1415 , -0.1638 ,  0.2062 ,  0.1243 , -0.3418
     &,  0.3412 ,  0.1380 , -0.2387 ,  0.1700 , -0.2601 ,  0.1572
     &, -0.1972 ,  0.1447 , -0.1134 ,  0.1417 , -0.3479 , -0.1354
     &,  0.7341 ,  0.2747 , -0.5894 ,  0.1212 , -0.0414 ,  0.0810
     &, -0.0012 , -0.1513 ,  0.3866 , -0.1141 ,  0.2317 , -0.0170
     &,  0.2681 , -0.1718 ,  0.3911 , -0.4157 ,  0.0188 ,  0.5973/

      data (ambchg(i),i=181,240) /
     &   0.2719 , -0.5679 ,  0.0881 , -0.0462 ,  0.0402 , -0.0266
     &, -0.3811 ,  0.3649 ,  0.1292 ,  0.1147 ,  0.2057 ,  0.1392
     &, -0.5727 , -0.4157 , -0.0581 ,  0.5973 ,  0.2719 , -0.5679
     &,  0.1360 , -0.0074 ,  0.0367 ,  0.1868 , -0.5432 , -0.2207
     &,  0.1862 ,  0.1635 ,  0.1435 , -0.2795 ,  0.3339 , -0.5163
     &,  0.0381 ,  0.5366 ,  0.2936 , -0.5819 ,  0.0880 , -0.0303
     &, -0.0122 ,  0.7994 , -0.8014 , -0.4157 ,  0.0143 ,  0.5973
     &,  0.2719 , -0.5679 ,  0.1048 , -0.2041 ,  0.0797 ,  0.7130
     &, -0.5931 , -0.9191 ,  0.4196 , -0.5163 ,  0.0397 ,  0.5366
     &,  0.2936 , -0.5819 ,  0.1105 ,  0.0560 , -0.0173 ,  0.0136/

      data (ambchg(i),i=241,300) /
     &  -0.0425 ,  0.8054 , -0.8188 , -0.4157 , -0.0031 ,  0.5973
     &,  0.2719 , -0.5679 ,  0.0850 , -0.0036 ,  0.0171 , -0.0645
     &,  0.0352 ,  0.6951 , -0.6086 , -0.9407 ,  0.4251 , -0.4157
     &, -0.0237 ,  0.5973 ,  0.2719 , -0.5679 ,  0.0880 ,  0.0342
     &,  0.0241 ,  0.0018 ,  0.0440 , -0.2737 , -0.0536 ,  0.0684
     &, -0.3479 , -0.2400 ,  0.7341 ,  0.2747 , -0.5894 ,  0.1426
     &, -0.0094 ,  0.0362 ,  0.0187 ,  0.0103 , -0.0479 ,  0.0621
     &, -0.0143 ,  0.1135 , -0.3854 ,  0.3400 , -0.3479 , -0.2637
     &,  0.7341 ,  0.2747 , -0.5894 ,  0.1560 , -0.0007 ,  0.0327
     &,  0.0390 ,  0.0285 ,  0.0486 ,  0.0687 , -0.5295 ,  0.3456/

      data (ambchg(i),i=301,360) /
     &   0.8076 , -0.8627 ,  0.4478 , -0.3479 , -0.2400 ,  0.7341
     &,  0.2747 , -0.5894 ,  0.1426 ,  0.0099 ,  0.0362 , -0.0279
     &,  0.0621 , -0.0143 ,  0.1135 , -0.3854 ,  0.3400 , -0.4157
     &,  0.1032 ,  0.5973 ,  0.2719 , -0.5679 , -0.2503 ,  0.0853
     &, -0.4157 , -0.0371 ,  0.5973 ,  0.2719 , -0.5679 ,  0.0850
     &, -0.0036 ,  0.0171 , -0.0645 ,  0.0352 ,  0.6000 , -0.5700
     &,  0.5500 ,  0.0000 , -0.5500 , -0.3662 ,  0.1123 ,  0.5972
     &, -0.5679 , -0.8400 ,  0.4200 , -0.4157 ,  0.2719 , -0.1490
     &,  0.0976 ,  0.2943 , -0.0100 ,  0.6163 ,  0.1642 , -0.5722
     &,  0.0895 ,  0.1414 ,  0.1281 ,  0.6163 ,  0.1997 , -0.5722/

      data (ambchg(i),i=361,420) /
     &   0.0889 ,  0.0577 ,  0.0023 ,  0.6163 ,  0.2272 , -0.5722
     &,  0.1093 ,  0.1010 ,  0.0343 ,  0.6123 ,  0.2148 , -0.5713
     &,  0.1053 ,  0.0311 ,  0.0389 ,  0.6123 ,  0.2329 , -0.5713
     &,  0.1031 ,  0.1849 ,  0.0684 ,  0.6163 ,  0.1898 , -0.5722
     &,  0.0782 ,  0.1812 ,  0.0332 ,  0.6163 ,  0.1934 , -0.5722
     &,  0.1087 ,  0.1325 ,  0.0978 ,  0.6123 ,  0.2023 , -0.5713
     &,  0.1411 ,  0.2096 ,  0.1138 ,  0.6163 ,  0.1815 , -0.5713
     &,  0.0922 , -0.2020 ,  0.1589 ,  0.5260 ,  0.3120 , -0.5000
     &,  0.1000 , -0.0120 ,  0.1000 ,  0.1737 ,  0.0859 ,  0.6123
     &,  0.1921 , -0.5713 ,  0.1041 ,  0.1940 ,  0.0766 ,  0.6123/

      data (ambchg(i),i=421,480) /
     &   0.1873 , -0.5713 ,  0.0983 ,  0.1913 ,  0.0555 ,  0.6123
     &,  0.1888 , -0.5713 ,  0.1162 ,  0.2560 ,  0.0653 ,  0.7214
     &,  0.1704 , -0.6013 ,  0.1047 ,  0.1542 ,  0.1126 ,  0.6123
     &,  0.1963 , -0.5713 ,  0.0958 ,  0.1472 ,  0.0325 ,  0.6123
     &,  0.2016 , -0.5713 ,  0.1380 ,  0.0782 ,  0.0326 ,  0.5621
     &,  0.2200 , -0.5889 ,  0.1141 ,  0.1801 ,  0.0811 ,  0.6163
     &,  0.1921 , -0.5722 ,  0.1231 ,  0.0017 ,  0.0698 ,  0.5621
     &,  0.2391 , -0.5889 ,  0.1202 ,  0.1493 ,  0.0769 ,  0.6123
     &,  0.1996 , -0.5713 ,  0.1015 ,  0.1592 ,  0.0429 ,  0.6123
     &,  0.1984 , -0.5713 ,  0.1116 ,  0.0966 , -0.0101 ,  0.7214/

      data (ambchg(i),i=481,540) /
     &   0.2165 , -0.6013 ,  0.1180 ,  0.1305 , -0.0359 ,  0.7214
     &,  0.2083 , -0.6013 ,  0.1242 ,  0.0966 , -0.0101 ,  0.7214
     &,  0.2165 , -0.6013 ,  0.1180 ,  0.1000 ,  0.2447 ,  0.6163
     &,  0.2000 , -0.5722 , -0.3821 , -0.2493 ,  0.7231 ,  0.2681
     &, -0.7855 ,  0.1056 , -0.3821 , -0.1532 ,  0.7731 ,  0.2681
     &, -0.8055 ,  0.1067 , -0.3821 , -0.3352 ,  0.8350 ,  0.2681
     &, -0.8173 ,  0.1438 , -0.3821 , -0.2874 ,  0.8326 ,  0.2681
     &, -0.8199 ,  0.1346 , -0.3821 , -0.3070 ,  0.8343 ,  0.2681
     &, -0.8190 ,  0.1375 , -0.3821 , -0.2563 ,  0.8113 ,  0.2681
     &, -0.8132 ,  0.1304 , -0.3821 , -0.2315 ,  0.7810 ,  0.2681/

      data (ambchg(i),i=541,600) /
     &  -0.8044 ,  0.1207 , -0.3821 , -0.1598 ,  0.7497 ,  0.2681
     &, -0.7981 ,  0.1396 , -0.3821 , -0.1283 ,  0.7618 ,  0.2681
     &, -0.8041 ,  0.0938 , -0.2802 , -0.1236 ,  0.6631 , -0.7697
     &,  0.0776 , -0.3821 , -0.1756 ,  0.7660 ,  0.2681 , -0.8026
     &,  0.1098 , -0.3821 , -0.1911 ,  0.7817 ,  0.2681 , -0.8070
     &,  0.1092 , -0.3821 , -0.2064 ,  0.7658 ,  0.2681 , -0.8011
     &,  0.1272 , -0.3481 , -0.1503 ,  0.8032 ,  0.2764 , -0.8177
     &,  0.1115 , -0.3821 , -0.1618 ,  0.7615 ,  0.2681 , -0.8016
     &,  0.1100 , -0.3821 , -0.2661 ,  0.7916 ,  0.2681 , -0.8065
     &,  0.1650 , -0.5192 , -0.1810 ,  0.7256 ,  0.3055 , -0.7887/

      data (ambchg(i),i=601,659) /
     &   0.1046 , -0.3821 , -0.1927 ,  0.8050 ,  0.2681 , -0.8147
     &,  0.1358 , -0.5192 , -0.2000 ,  0.7420 ,  0.3055 , -0.7930
     &,  0.1399 , -0.3821 , -0.2108 ,  0.7775 ,  0.2681 , -0.8042
     &,  0.1232 , -0.3821 , -0.2441 ,  0.8013 ,  0.2681 , -0.8105
     &,  0.1277 , -0.3481 , -0.2964 ,  0.8488 ,  0.2764 , -0.8252
     &,  0.1438 , -0.3481 , -0.3117 ,  0.8557 ,  0.2764 , -0.8266
     &,  0.1447 , -0.3481 , -0.2964 ,  0.8488 ,  0.2764 , -0.8252
     &,  0.1438 , -0.3821 , -0.0572 ,  0.8000 ,  0.2681 , -0.8200
     &, -0.8340 ,  0.4170 ,  1.0000 ,  1.0000 ,  1.0000 ,  1.0000
     &,  1.0000 ,  2.0000 ,  2.0000 ,  2.0000 , -1.0000/

c ASP neutral
      data (ambchg(i),i=660,671) /
     &  -0.4157 ,  0.0341 ,  0.5973 ,  0.2719 , -0.5679 ,  0.0864
     &, -0.0316 ,  0.0488 ,  0.6462 , -0.5554 , -0.6376 ,  0.4747/

c GLU neutral
      data (ambchg(i),i=672,685) /
     &  -0.4157 ,  0.0145 ,  0.5973 ,  0.2719 , -0.5679 ,  0.0779
     &, -0.0071 ,  0.0256 , -0.0174 ,  0.0430 ,  0.6801 , -0.5838
     &, -0.6511 ,  0.4641/

c LYS neutral
      data (ambchg(i),i=686,701) /
     &  -0.4157 , -0.07206,  0.5973 ,  0.2719 , -0.5679 ,  0.0994
     &, -0.04845,  0.0340 , 0.06612 , 0.01041 ,-0.03768 , 0.01155
     &,  0.32604, -0.03358, -1.03581, 0.38604/

c CYS Negative

      data (cysneg(i),i=1,9) /
     & -0.4937, -0.2870, 0.6731, 0.3018, -0.5854, 0.1203, 0.0935,
     &  0.0125, -0.8476/

      data (ambchg(i),i=702,1000) / 299*0.0/

      data (ambchg(i),i=1001,1030) /
     & -0.4989,  0.0558,  0.0679,  0.0679,  0.1065,  0.1174,
     & -0.3548,  0.0394,  0.2007,  0.2022,  0.0615,  0.0670,
     &  0.0972, -0.6139,  0.4186, -0.5246, -0.0251,  0.3053,
     &  0.0515, -0.6073,  0.2006, -0.6997,  0.5875, -0.7615,
     &  0.7009,  0.0473, -0.9019,  0.4115,  0.4115,  0.1553/

      data (ambchg(i),i=1031,1061) /
     & -0.4989,  0.0558,  0.0679,  0.0679,  0.1065,  0.1174,
     & -0.3548,  0.0191,  0.2006,  0.2022,  0.0615,  0.0670,
     &  0.0972, -0.6139,  0.4186, -0.5246,  0.0492,  0.1222,
     &  0.1744, -0.5709,  0.1374, -0.6323,  0.7657, -0.4787,
     &  0.4770,  0.3424, -0.9672,  0.4364,  0.4364, -0.5597,
     &  0.1640/

      data (ambchg(i),i=1062,1089) /
     & -0.4989,  0.0558,  0.0679,  0.0679,  0.1065,  0.1174,
     & -0.3548,  0.0066,  0.2029,  0.2022,  0.0615,  0.0670,
     &  0.0972, -0.6139,  0.4186, -0.5246, -0.0484,  0.7538,
     & -0.7584,  0.8185, -0.5215,  0.0053, -0.6252, -0.9530,
     &  0.4234,  0.4234,  0.1928,  0.1958/

      data (ambchg(i),i=1090,1116) /
     & -0.4989,  0.0558,  0.0679,  0.0679,  0.1065,  0.1174,
     & -0.3548,  0.0674,  0.1824,  0.2022,  0.0615,  0.0670,
     &  0.0972, -0.6139,  0.4186, -0.5246,  0.0418,  0.4687,
     & -0.3549,  0.5952, -0.3635, -0.1126, -0.5477,  0.3154,
     & -0.5761,  0.1811,  0.2188/

      data (ambchg(i),i=1117,1145) /
     & -0.4954, -0.0069,  0.0754,  0.0754,  0.1629,  0.1176,
     & -0.3691,  0.0431,  0.1838,  0.0713,  0.0985, -0.0854,
     &  0.0718,  0.0718, -0.5232, -0.0268,  0.3800,  0.0725,
     & -0.6175,  0.1607, -0.7417,  0.5716, -0.7624,  0.6897,
     &  0.0598, -0.9123,  0.4167,  0.4167,  0.1877/

      data (ambchg(i),i=1146,1175) /
     & -0.4954, -0.0069,  0.0754,  0.0754,  0.1629,  0.1176,
     & -0.3691,  0.0358,  0.1746,  0.0713,  0.0985, -0.0854,
     &  0.0718,  0.0718, -0.5232,  0.0577,  0.1814,  0.1991,
     & -0.5725,  0.0736, -0.6636,  0.7432, -0.5053,  0.4918,
     &  0.3520, -0.9230,  0.4235,  0.4235, -0.5699,  0.1997/

      data (ambchg(i),i=1176,1202) /
     & -0.4954, -0.0069,  0.0754,  0.0754,  0.1629,  0.1176,
     & -0.3691, -0.0116,  0.1963,  0.0713,  0.0985, -0.0854,
     &  0.0718,  0.0718, -0.5232, -0.0339,  0.7959, -0.7748,
     &  0.8439, -0.5222, -0.0183, -0.6548, -0.9773,  0.4314,
     &  0.4314,  0.1863,  0.2293/

      data (ambchg(i),i=1203,1229) /
     & -0.4954, -0.0069,  0.0754,  0.0754,  0.1629,  0.1176,
     & -0.3691,  0.0680,  0.1804,  0.0713,  0.0985, -0.0854,
     &  0.0718,  0.0718, -0.5232, -0.0239,  0.5677, -0.4340,
     &  0.5194,  0.0025, -0.2209, -0.5881,  0.3420, -0.5563,
     & -0.2269,  0.0770,  0.2607/

      data (ambchg(i),i=1230,1241) /
     &  1.1662, -0.7760, -0.6223,  0.4295, -0.6328,  0.9000,
     & -0.8200, -0.6541,  0.4376, -0.6565,  0.9000, -0.8200/

      data (ambchg(i),i=1242,1253) /
     &  1.1659, -0.7761, -0.6318,  0.4422, -0.6296,  0.9000,
     & -0.8200, -0.6549,  0.4396, -0.6553,  0.9000, -0.8200/

      data (ambvw1(i),i=1,49) /
     & 1.9080, 1.9080, 1.9080, 1.9080, 1.9080, 1.9080, 1.9080,
     & 1.9080, 1.9080, 1.9080, 1.9080, 1.9080, 1.9080, 1.8240,
     & 1.8240, 1.8240, 1.8240, 1.8240, 1.8240, 1.8750, 1.7683,
     & 1.7210, 1.6837, 1.6612, 1.6612, 2.0000, 2.0000, 2.1000,
     & 0.6000, 0.0001, 0.0001, 0.6000, 1.4590, 1.4870, 1.3870,
     & 1.2870, 1.1870, 1.1000, 1.4090, 1.3590, 1.1370, 1.8680,
     & 2.6580, 2.9560, 3.3950, 0.7926, 1.7131, 1.1000, 2.4700/

      data (ambvw2(i),i=1,49) /
     & 0.1094, 0.0860, 0.0860, 0.0860, 0.0860, 0.0860, 0.0860,
     & 0.0860, 0.0860, 0.0860, 0.0860, 0.0860, 0.0860, 0.1700,
     & 0.1700, 0.1700, 0.1700, 0.1700, 0.1700, 0.1700, 0.1520,
     & 0.2104, 0.1700, 0.2100, 0.2100, 0.2500, 0.2500, 0.2000,
     & 0.0157, 0.0000, 0.0000, 0.0157, 0.0150, 0.0157, 0.0157,
     & 0.0157, 0.0157, 0.0157, 0.0150, 0.0150, 0.0183, 0.00277,
     & 0.000328, 0.00017, 0.0000806, 0.8947, 0.459789, 
     & 0.0125, 0.1000/

      data (ambvdt(i),i=1,60) /
     & 14, 1, 2,29,24,35,14, 1, 2,29,
     & 24,35, 1,34,14, 1, 2,29,24,35,
     &  1,34, 1,34, 1,34,14, 1, 2,29,
     & 24,35, 1,34, 1,34, 1,34, 1,34,
     & 14, 1, 2,29,24,35, 1,34, 1,34,
     &  1,34, 1,34,14, 1, 2,29,24,35/

      data (ambvdt(i),i=61,120) /
     &  1,35,22,31,14, 1, 2,29,24,35,
     &  1,35,22,31, 1,34,14, 1, 2,29,
     & 24,35, 1,35,27,32,14, 1, 2,29,
     & 24,35, 1,35,26,14, 1, 2,24,35,
     &  1,34, 1,34, 1,35,14, 1, 2,29,
     & 24,35, 1,34, 3, 3,33, 3,33, 3/

      data (ambvdt(i),i=121,180) /
     & 33,14, 1, 2,29,24,35, 1,34, 3,
     &  3,33, 3,33, 3,22,31,14, 1, 2,
     & 29,24,35, 1,34,10, 7,39, 9,15,
     & 29,11, 3,33, 3,33, 3,33, 3,33,
     & 14, 1, 2,29,24,35, 1,34, 5,15,
     & 29, 7,39, 8,40,15,29,14, 1, 2/

      data (ambvdt(i),i=181,240) /
     & 29,24,35, 1,34, 5,15,29, 6,39,
     &  8,40,16,14, 1, 2,29,24,35, 1,
     & 34, 5,16, 7,39, 8,40,15,29,14,
     &  1, 2,29,24,35, 1,34, 2,25,14,
     &  1, 2,29,24,35, 1,34, 2,24,14,
     & 29,14, 1, 2,29,24,35, 1,34, 1/

      data (ambvdt(i),i=241,300) /
     & 34, 2,25,14, 1, 2,29,24,35, 1,
     & 34, 1,34, 2,24,14,29,14, 1, 2,
     & 29,24,35, 1,34, 1,35,26, 1,35,
     & 14, 1, 2,29,24,35, 1,34, 1,34,
     &  1,34, 1,38,20,29,14, 1, 2,29,
     & 24,35, 1,34, 1,34, 1,35,19,29/

      data (ambvdt(i),i=301,360) /
     &  3,19,29,14, 1, 2,29,24,35, 1,
     & 34, 1,34, 1,38,20,29,14, 1, 2,
     & 29,24, 1,34,14, 1, 2,29,24,35,
     &  1,34, 1,34, 2,24, 2, 0,24, 1,
     & 34, 2,24,14,29,14,29, 1,35,20,
     &  1, 2,29,24,35,20, 1, 2,29,24/

      data (ambvdt(i),i=361,420) /
     & 35,20, 1, 2,29,24,35,20, 1, 2,
     & 29,24,35,20, 1, 2,29,24,35,20,
     &  1, 2,29,24,35,20, 1, 2,29,24,
     & 35,20, 1, 2,29,24,35,20, 1, 2,
     & 29,24,35,20, 1, 2,29,24,35, 1,
     & 38,20, 1, 2,29,24,35,20, 1, 2/

      data (ambvdt(i),i=421,480) /
     & 29,24,35,20, 1, 2,29,24,35,20,
     &  1, 2,29,24,35,20, 1, 2,29,24,
     & 35,20, 1, 2,29,24,35,20, 1, 2,
     & 29,24,35,20, 1, 2,29,24,35,20,
     &  1, 2,29,24,35,20, 1, 2,29,24,
     & 35,20, 1, 2,29,24,35,20, 1, 2/

      data (ambvdt(i),i=481,540) /
     & 29,24,35,20, 1, 2,29,24,35,20,
     &  1, 2,29,24,35,20, 1, 2,29,24,
     & 14, 1, 2,29,25,35,14, 1, 2,29,
     & 25,35,14, 1, 2,29,25,35,14, 1,
     &  2,29,25,35,14, 1, 2,29,25,35,
     & 14, 1, 2,29,25,35,14, 1, 2,29/

      data (ambvdt(i),i=541,600) /
     & 25,35,14, 1, 2,29,25,35,14, 1,
     &  2,29,25,35,14, 1, 2,25,35,14,
     &  1, 2,29,25,35,14, 1, 2,29,25,
     & 35,14, 1, 2,29,25,35,14, 1, 2,
     & 29,25,35,14, 1, 2,29,25,35,14,
     &  1, 2,29,25,35,14, 1, 2,29,25/

      data (ambvdt(i),i=601,648) /
     & 35,14, 1, 2,29,25,35,14, 1, 2,
     & 29,25,35,14, 1, 2,29,25,35,14,
     &  1, 2,29,25,35,14, 1, 2,29,25,
     & 35,14, 1, 2,29,25,35,14, 1, 2,
     & 29,25,35,14, 1, 2,29,25/

      data (ambvdt(i),i=649,659) /
     & 32,30,41,42,43,44,45,46,47,48,
     & 49/

c neutral asp, glu, lys
      data (ambvdt(i),i=660,671) /
     & 14, 1, 2,29,24,35, 1,34,
     &  2,24,22,31/
      data (ambvdt(i),i=672,685) /
     & 14, 1, 2,29,24,35, 1,34,
     &  1,34, 2,24,22,31/
      data (ambvdt(i),i=686,701) /
     & 14, 1, 2,29,24,35, 1,34,
     &  1,34, 1,34, 1,38,20,29/

      data (ambvdt(i),i=702,1000) / 299*0/

      data (ambvdt(i),i=1001,1060) /
     & 23,  1, 35, 35,  1, 35, 23,  1, 36,  1,
     & 35,  1, 35, 22, 31, 23, 18,  9,  9, 16,
     & 12, 17, 13, 17,  3, 40, 19, 29, 29, 40,
     & 23,  1, 35, 35,  1, 35, 23,  1, 36,  1,
     & 35,  1, 35, 22, 31, 23, 18,  9,  9, 16,
     & 12, 17,  3, 15,  2, 29, 19, 29, 29, 24/
       
      data (ambvdt(i),i=1061,1120) /
     & 40, 23,  1, 35, 35,  1, 35, 23,  1, 36,
     &  1, 35,  1, 35, 22, 31, 23, 18,  2, 17,
     &  3,  4,  4, 24, 19, 29, 29, 33, 39, 23,
     &  1, 35, 35,  1, 35, 23,  1, 36,  1, 35,
     &  1, 35, 22, 31, 23, 18,  2, 15,  2,  4,
     &  4, 24, 29, 24, 33, 39, 23,  1, 35, 35/

      data (ambvdt(i),i=1121,1180) /
     &  1, 35, 23,  1, 36,  1, 35,  1, 34, 34,
     & 23, 18,  9,  9, 16, 12, 17, 13, 17,  3,
     & 40, 19, 29, 29, 40, 23,  1, 35, 35,  1,
     & 35, 23,  1, 36,  1, 35,  1, 34, 34, 23,
     & 18,  9,  9, 16, 12, 17,  3, 15,  2, 29,
     & 19, 29, 29, 24, 40, 23,  1, 35, 35,  1/
       
      data (ambvdt(i),i=1181,1253) /
     & 35, 23,  1, 36,  1, 35,  1, 34, 34, 23,
     & 18,  2, 17,  3,  4,  4, 24, 19, 29, 29,
     & 33, 39, 23,  1, 35, 35,  1, 35, 23,  1,
     & 36,  1, 35,  1, 34, 34, 23, 18,  2, 15,
     &  2,  4,  4, 24, 29, 24,  1, 34, 39, 28,
     & 25, 22, 31, 23, 28, 25, 22, 31, 23, 28,
     & 25, 28, 25, 22, 31, 23, 28, 25, 22, 31,
     & 23, 28, 25/

c GAFF parameters

      data gfvdw /
     &   0.0000, 0.0000,
     &   1.9080, 0.0860,
     &   1.9080, 0.0860,
     &   1.9080, 0.0860,
     &   1.9080, 0.1094,
     &   1.9080, 0.0860,
     &   1.9080, 0.0860,
     &   1.9080, 0.0860,
     &   1.9080, 0.0860,
     &   1.9080, 0.0860,
     &   1.9080, 0.0860,
     &   1.9080, 0.0860,
     &   1.9080, 0.0860,
     &   1.9080, 0.0860,
     &   1.9080, 0.0860,
     &   1.9080, 0.0860,
     &   1.9080, 0.0860,
     &   1.9080, 0.0860,
     &   1.3870, 0.0157,
     &   1.2870, 0.0157,
     &   1.1870, 0.0157,
     &   1.4090, 0.0150,
     &   1.3590, 0.0150,
     &   1.4590, 0.0150,
     &   1.4870, 0.0157,
     &   0.6000, 0.0157,
     &   0.0000, 0.0000,
     &   0.6000, 0.0157,
     &   0.6000, 0.0157,
     &   0.0000, 0.0000,
     &   1.1000, 0.0157,
     &   1.75  , 0.061 ,
     &   1.948 , 0.265 ,
     &   2.22  , 0.320 ,
     &   2.35  , 0.40  ,
     &   1.8240, 0.1700,
     &   1.8240, 0.1700,
     &   1.8240, 0.1700,
     &   1.8240, 0.1700,
     &   1.8240, 0.1700,
     &   1.8240, 0.1700,
     &   1.8240, 0.1700,
     &   1.8240, 0.1700,
     &   1.8240, 0.1700,
     &   1.8240, 0.1700,
     &   1.8240, 0.1700,
     &   1.8240, 0.1700,
     &   1.8240, 0.1700,
     &   1.6612, 0.2100,
     &   1.7210, 0.2104,
     &   1.6837, 0.1700,
     &   1.7683, 0.1520,
     &   2.1000, 0.2000,
     &   2.1000, 0.2000,
     &   2.1000, 0.2000,
     &   2.1000, 0.2000,
     &   2.1000, 0.2000,
     &   2.1000, 0.2000,
     &   2.1000, 0.2000,
     &   2.1000, 0.2000,
     &   2.1000, 0.2000,
     &   2.1000, 0.2000,
     &   2.1000, 0.2000,
     &   2.0000, 0.2500,
     &   2.0000, 0.2500,
     &   2.0000, 0.2500,
     &   2.0000, 0.2500,
     &   2.0000, 0.2500,
     &   2.0000, 0.2500,
     &   2.0000, 0.2500,
     &   2.0000, 0.2500,
     &   1.9080, 0.0860/

      data ihasl /-1,19,3,5,3,5,3,2,5,5,5/

      data gro43 /
     &     '    O','   OM','   OA','   OW','    N','   NT','   NL',
     &     '   NR','   NZ','   NE','    C','  CH1','  CH2','  CH3',
     &     '  CH4','  CR1','   HC','    H','  DUM','    S',' CU1+',
     &     ' CU2+','   FE',' ZN2+',' MG2+',' CA2+','    P','   AR',
     &     '    F','   CL','   BR',' CMET',' OMET','  NA+','  CL-',
     &     ' CCHL','CLCHL',' HCHL','SDMSO','CDMSO','ODMSO',' CCL4',
     &     'CLCL4','   SI',' MNH3','   MW','   IW',' OWT3',' OWT4'/
      data gro43l /
     &     'O    ','OM   ','OA   ','OW   ','N    ','NT   ','NL   ',
     &     'NR   ','NZ   ','NE   ','C    ','CH1  ','CH2  ','CH3  ',
     &     'CH4  ','CR1  ','HC   ','H    ','DUM  ','S    ','CU1+ ',
     &     'CU2+ ','FE   ','ZN2+ ','MG2+ ','CA2+ ','P    ','AR   ',
     &     'F    ','CL   ','BR   ','CMET ','OMET ','NA+  ','CL-  ',
     &     'CCHL ','CLCHL','HCHL ','SDMSO','CDMSO','ODMSO','CCL4 ',
     &     'CLCL4','SI   ','MNH3 ','MW   ','IW   ','OWT3 ','OWT4 '/
      data gro43s /
     &  'carbonyl oxygen (C=O)              ',
     &  'carboxyl oxygen (CO-)              ',
     &  'hydroxyl, sugar or ester oxygen    ',
     &  'water oxygen                       ',
     &  'peptide nitrogen (N or NH)         ',
     &  'terminal nitrogen (NH2)            ',
     &  'terminal nitrogen (NH3)            ',
     &  'aromatic nitrogen                  ',
     &  'Arg NH (NH2)                       ',
     &  'Arg NE (NH)                        ',
     &  'bare carbon                        ',
     &  'aliphatic or sugar CH-group        ',
     &  'aliphatic or sugar CH2-group       ',
     &  'aliphatic CH3-group                ',
     &  'methane                            ',
     &  'aromatic CH-group                  ',
     &  'hydrogen bound to carbon           ',
     &  'hydrogen not bound to carbon       ',
     &  'dummy atom                         ',
     &  'sulfur                             ',
     &  'copper (charge 1+)                 ',
     &  'copper (charge 2+)                 ',
     &  'iron (heme)                        ',
     &  'zinc (charge 2+)                   ',
     &  'magnesium (charge 2+)              ',
     &  'calcium (charge 2+)                ',
     &  'phosphor                           ',
     &  'argon                              ',
     &  'fluor (non-ionic)                  ',
     &  'chlorine (non-ionic)               ',
     &  'bromine (non-ionic)                ',
     &  'CH3-group in methanol (solvent)    ',
     &  'oxygen in methanol (solvent)       ',
     &  'sodium (charge 1+)                 ',
     &  'chlorine (charge 1-)               ',
     &  'carbon in chloroform (solvent)     ',
     &  'chloride in chloroform (solvent)   ',
     &  'hydrogen in chloroform (solvent)   ',
     &  'DMSO Sulphur (solvent)             ',
     &  'DMSO Carbon (solvent)              ',
     &  'DMSO Oxygen (solvent)              ',
     &  'carbon in CCl4 (solvent)           ',
     &  'chloride in CCl4 (solvent)         ',
     &  'silicon                            ',
     &  'Dummy mass in rig. tetraed. NH3 grp',
     &  'Dummy mass in rig. tyrosine rings  ',
     &  'Dummy particle in TIP4P etc.       ',
     &  'TIP3P WATER OXYGEN                 ',
     &  'TIP4P WATER OXYGEN                 '/

      data grogmx /
     &  '    O','   OM','   OA','   OW','    N','   NT','   NL',
     &  '  NR5',' NR5*','   NP','    C','  CH1','  CH2','  CH3',
     &  ' CR51',' CR61','   CB','    H','   HO','   HW','   HS',
     &  '    S','   FE','   ZN','   NZ','   NE','    P','   OS',
     &  '  CS1','  NR6',' NR6*','  CS2','   SI','   NA','    K',
     &  '   CL','   CA','   MG','    F','  CP2','  CP3','  CR5',
     &  '  CR6','  HCR',' OWT3',' OWT4','   SD','   OD','   CD',
     &  '  CHE',' MNH3','   MW','   IW'/
      data grogmxl /
     &  'O    ','OM   ','OA   ','OW   ','N    ','NT   ','NL   ',
     &  'NR5  ','NR5* ','NP   ','C    ','CH1  ','CH2  ','CH3  ',
     &  'CR51 ','CR61 ','CB   ','H    ','HO   ','HW   ','HS   ',
     &  'S    ','FE   ','ZN   ','NZ   ','NE   ','P    ','OS   ',
     &  'CS1  ','NR6  ','NR6* ','CS2  ','SI   ','NA   ','K    ',
     &  'CL   ','CA   ','MG   ','F    ','CP2  ','CP3  ','CR5  ',
     &  'CR6  ','HCR  ','OWT3 ','OWT4 ','SD   ','OD   ','CD   ',
     &  'CHE  ','MNH3 ','MW   ','IW   '/
      data grogms /
     & 'CARBONYL OXYGEN (C=O)              ',
     & 'CARBOXYL OXYGEN (CO-)              ',
     & 'HYDROXYL OXYGEN (OH)               ',
     & 'WATER OXYGEN                       ',
     & 'PEPTIDE NITROGEN (N OR NH)         ',
     & 'TERMINAL NITROGEN (NH2)            ',
     & 'TERMINAL NITROGEN (NH3)            ',
     & 'AROMATIC N (5-RING,2 BONDS)        ',
     & 'AROMATIC N (5-RING,3 BONDS)        ',
     & 'PORPHYRIN NITROGEN                 ',
     & 'BARE CARBON (PEPTIDE,C=O,C-N)      ',
     & 'ALIPHATIC CH-GROUP                 ',
     & 'ALIPHATIC CH2-GROUP                ',
     & 'ALIPHATIC CH3-GROUP                ',
     & 'AROMATIC CH-GROUP (5-RING), united ',
     & 'AROMATIC CH-GROUP (6-RING), united ',
     & 'BARE CARBON (5-,6-RING)            ',
     & 'HYDROGEN BONDED TO NITROGEN        ',
     & 'HYDROXYL HYDROGEN                  ',
     & 'WATER HYDROGEN                     ',
     & 'HYDROGEN BONDED TO SULFUR          ',
     & 'SULFUR                             ',
     & 'IRON                               ',
     & 'ZINC                               ',
     & 'ARG NH (NH2)                       ',
     & 'ARG NE (NH)                        ',
     & 'PHOSPHOR                           ',
     & 'SUGAR OR ESTER OXYGEN              ',
     & 'SUGAR CH-GROUP                     ',
     & 'AROMATIC N (6-RING,2 BONDS)        ',
     & 'AROMATIC N (6-RING,3 BONDS)        ',
     & 'SUGAR CH2-GROUP                    ',
     & 'SILICON                            ',
     & 'SODIUM (1+)                        ',
     & 'POTASSIUM (1+) (Thesis Straatsma)  ',
     & 'CHLORINE (1-)                      ',
     & 'CALCIUM (2+)                       ',
     & 'MAGNESIUM (2+)                     ',
     & 'FLUORINE (COV. BOUND)              ',
     & 'ALIPHATIC CH2-GROUP USING Ryckaert-',
     & 'ALIPHATIC CH3-GROUP USING Ryckaert-',
     & 'AROMATIC CH-GROUP (5-RING)+H       ',
     & 'AROMATIC C- bonded to H (6-RING)+H ',
     & 'H attached to aromatic C (5 or 6 ri',
     & 'TIP3P WATER OXYGEN                 ',
     & 'TIP4P WATER OXYGEN                 ',
     & 'DMSO Sulphur                       ',
     & 'DMSO Oxygen                        ',
     & 'DMSO Carbon                        ',
     & 'HEME RING CARBON                   ',
     & 'Dummy mass in rig. tetraed. NH3 grp',
     & 'Dummy mass in rig. tyrosine rings  ',
     & 'Dummy particle in TIP4P etc.       '/
      data grog2x /
     &  '    O','   OM','   OA','   OW','    N','   NT','   NL',
     &  '  NR5',' NR5*','   NP','    C','  CH1','  CH2','  CH3',
     &  ' CR51',' CR61','   CB','    H','   HO','   HW','   HS',
     &  '    S','   FE','   ZN','   NZ','   NE','    P','   OS',
     &  '  CS1','  NR6',' NR6*','  CS2','   SI','   NA','   CL',
     &  '   CA','   MG','    F','     ','  CP2','  CP3','  CR5',
     &  '  CR6','  HCR',' OWT3',' OWT4','   SD','   OD','   CD',
     &  '   HC',' MNH3',' MCH3','   MW','  CHE','  CH3','  CH2',
     &  '   IW'/
      data grog2xl /
     &  'O    ','OM   ','OA   ','OW   ','N    ','NT   ','NL   ',
     &  'NR5  ','NR5* ','NP   ','C    ','CH1  ','CH2  ','CH3  ',
     &  'CR51 ','CR61 ','CB   ','H    ','HO   ','HW   ','HS   ',
     &  'S    ','FE   ','ZN   ','NZ   ','NE   ','P    ','OS   ',
     &  'CS1  ','NR6  ','NR6* ','CS2  ','SI   ','NA   ','CL   ',
     &  'CA   ','MG   ','F    ','     ','CP2  ','CP3  ','CR5  ',
     &  'CR6  ','HCR  ','OWT3 ','OWT4 ','SD   ','OD   ','CD   ',
     &  'HC   ','MNH3 ','MCH3 ','MW   ','CHE  ','CH3  ','CH2  ',
     &  'IW   '/
      data grog2s /
     & 'CARBONYL OXYGEN (C=O)              ',
     & 'CARBOXYL OXYGEN (CO-)              ',
     & 'HYDROXYL OXYGEN (OH)               ',
     & 'WATER OXYGEN                       ',
     & 'PEPTIDE NITROGEN (N OR NH)         ',
     & 'TERMINAL NITROGEN (NH2)            ',
     & 'TERMINAL NITROGEN (NH3)            ',
     & 'AROMATIC N (5-RING,2 BONDS)        ',
     & 'AROMATIC N (5-RING,3 BONDS)        ',
     & 'PORPHYRIN NITROGEN                 ',
     & 'BARE CARBON (PEPTIDE,C=O,C-N)      ',
     & 'ALIPHATIC CH-GROUP                 ',
     & 'ALIPHATIC CH2-GROUP                ',
     & 'ALIPHATIC CH3-GROUP                ',
     & 'AROMATIC CH-GROUP (5-RING), united ',
     & 'AROMATIC CH-GROUP (6-RING), united ',
     & 'BARE CARBON (5-,6-RING)            ',
     & 'HYDROGEN BONDED TO NITROGEN        ',
     & 'HYDROXYL HYDROGEN                  ',
     & 'WATER HYDROGEN                     ',
     & 'HYDROGEN BONDED TO SULFUR          ',
     & 'SULFUR                             ',
     & 'IRON                               ',
     & 'ZINC                               ',
     & 'ARG NH (NH2)                       ',
     & 'ARG NE (NH)                        ',
     & 'PHOSPHOR                           ',
     & 'SUGAR OR ESTER OXYGEN              ',
     & 'SUGAR CH-GROUP                     ',
     & 'AROMATIC N (6-RING,2 BONDS)        ',
     & 'AROMATIC N (6-RING,3 BONDS)        ',
     & 'SUGAR CH2-GROUP                    ',
     & 'SILICON                            ',
     & 'SODIUM (1+)                        ',
     & 'CHLORINE (1-)                      ',
     & 'CALCIUM (2+)                       ',
     & 'MAGNESIUM (2+)                     ',
     & 'FLUORINE (COV. BOUND)              ',
     & '                                   ',
     & 'ALIPHATIC CH2-GROUP USING Ryckaert-',
     & 'ALIPHATIC CH3-GROUP USING Ryckaert-',
     & 'AROMATIC CH-GROUP (5-RING)+H       ',
     & 'AROMATIC C- bonded to H (6-RING)+H ',
     & 'H attached to aromatic C (5 or 6 ri',
     & 'TIP3P WATER OXYGEN                 ',
     & 'TIP4P WATER OXYGEN                 ',
     & 'DMSO Sulphur                       ',
     & 'DMSO Oxygen                        ',
     & 'DMSO Carbon                        ',
     & 'H attached to aliphatic carbon     ',
     & 'Dummy mass in rig. tetraed. NH3 grp',
     & 'Dummy mass in rig. tetraed. CH3 grp',
     & 'Dummy mass in rig. tyrosine rings  ',
     & 'HEME RING CARBON                   ',
     & 'HEME RING CARBON                   ',
     & 'HEME RING CARBON                   ',
     & 'Dummy particle in TIP4P etc.       '/

c Protein pmf types

      data ppmf /'CF','CP','cF','cP','CO','CN','NC','ND',
     &           'NR','OC','OA','OD','OW','SA','SD','HH'/
c Ligand pmf types
c see inipmf for types with no data
      data lpmf /'CF','CP','cF','cP','C3','CW','CO','CN',
     &           'NC','NP','NA','ND','NR','N0','NS','OC',
     &           'OA','OE','OR','OS','OD','P ','SA','SD',
     &           'HL','F '/

c sybyl mol2 types

      data mol2 /'Any','Hal','Het','Hev',
     & 'C.3  ','C.2  ','C.1  ','C.ar ','C.cat',
     & 'N.3  ','N.2  ','N.1  ','N.ar ','N.am ','N.pl3','N.4  ',
     & 'O.3  ','O.2  ','O.co2','O.spc','O.t3p',
     & 'S.3  ','S.2  ','S.O  ','S.O2 ','P.3  ',
     & 'H    ','H.spc','H.t3p',
     & 'F    ','Cl   ','Br   ','I    ','Si   ','LP   ','Du   ',
     & 'Na   ','K    ','Ca   ','Li   ','Al   '/

      data (mm3(i),i=1,58)/
     &'CSP3 ALKANE        ','CSP2 ALKENE        ','CSP2 CARBONYL      ',
     &'CSP  ALKYNE        ','H EXCEPT ON N,O,S  ','O: C-O-H, C-O-C    ',
     &'=O   CARBONYL      ','NSP3               ','NSP2               ',
     &'NSP                ','FLUORIDE           ','CHLORIDE           ',
     &'BROMIDE            ','IODIDE             ','-S-  SULFIDE       ',
     &'>S+  S-FONIUM      ','>S=O S-OXIDE       ','>SO2 SULFONE       ',
     &'SILANE             ',
     &'NOT USED           ','H: -OH  ALCOHOL    ','C: CYCLOPROPANE    ',
     &'H: NH  AMINE/IMINE ','H: COOH CARBOXYL   ','>P-  PHOSPHINE     ',
     &'TRIGONAL B         ','TETRAHEDRAL B      ','H: C=COH/AMIDE     ',
     &'C RADICAL          ','CARBONIUM ION      ','GERMANIUM          ',
     &'TIN                ','LEAD (IV)          ','SELENIUM           ',
     &'TELLURIUM          ','DEUTERIUM          ','-N=N-/AZO,-N=C-/PYR',
     &'CSP2 CYCPROPENE    ','NSP3 (AMMONIUM)    ','NSP2 (PYRROLE)     ',
     &'OSP2 (FURAN)       ','SSP2 (THIOPHENE)   ','-N=N-O AZOXY N     ',
     &'SH   THIOL         ','AZIDE (CENTER N)   ','N NITRO            ',
     &'O CARBOXYLATE      ','H AMMONIUM         ','O EPOXY            ',
     &'C BENZENE (LOC.)   ','HELIUM             ','NEON               ',
     &'ARGON              ','KRYPTON            ','XENON              ',
     &'CSP3 (CYCLOBUTANE) ','CSP2 (CYCLOBUTENE) ','C=O  (CYCLOBUTANE) '/
      data (mm3(i),i=59,118)/
     &'MAGNESIUM          ','P (V)              ','IRON (II)          ',
     &'IRON (III)         ','NICKEL (II)        ','NICKEL (III)       ',
     &'COBALT (II)        ','COBALT (III)       ','C=O CYCLOPROPANONE ',
     &'CUMULENE CSP =C=   ','AMINE OXIDE OXYGEN ','KETONIUM OXYGEN    ',
     &'KETONIUM CARBON    ','=N- IMINE (LOCALZD)','H-O ENOL/PHENOL    ',
     &'S=C                ','O-H, O-C (CARBOXYL)','O=C-C=O            ',
     &'O=C-O-H (ACID)     ','O=C-O-C (ESTER)    ','O=C-N< (AMIDE)     ',
     &'O=C-X (HALIDE)     ','O=C-C=C<           ','O=C-O-C=O          ',
     &'O=C(OH)(OH)        ','O=C(OH)(OC)        ','O=C(OH)(N<)        ',
     &'O=C(OH)(X)         ','O=C(OH)(C=C<)      ','O=C(OH)(O-C=O)     ',
     &'O=C(O-C)(OC)       ','O=C(O-C)(N<)       ','O=C(O-C)(X)        ',
     &'O=C(O-C)(C=C<)     ','O=C(O-C)(O-C=O)    ','O=C(N<)(N<)        ',
     &'O=C(N<)(X)         ','O=C(N<)(C=C<)      ','O=C(N<)(O-C=O)     ',
     &'O=C(X)(X)          ','O=C(X)(C=C<)       ','O=C(X)(O-C=O)      ',
     &'O=C(C=C)(C=C)      ','O=C(C=C)(O-C=O)    ','O=C(O-C=O)(O-C=O)  ',
     &'-S-S- DISULFIDE    ','-S- POLYSULFIDE    ','=C=O KETENE        ',
     &'-N=N- AZO (LOCAL)  ','=N-OH OXIME        ','-N= AZOXY (LOCAL)  ',
     &'-N(+)= IMMINIUM    ','-N(+)= PYRIDINIUM  ','NOT USED           ',
     &'C FERROCENE-H      ','C FERROCENE-C      ','O=C(C=O)(C=O)      ',
     &'O=C(C=O)(OH)       ','O=C(C=O)(OC)       ','O=C(C=O)(N<)       '/
      data (mm3(i),i=119,164)/
     &'O=C(C=O)(X)        ','O=C(C=O)(C=C<)     ','O=C(C=O)(O-C=O)    ',
     &'NOT USED           ','NOT USED           ','H-C ACETYLENE      ',
     &'CALCIUM            ','STRONTIUM          ','BARIUM             ',
     &'LANTHANUM          ','CERIUM             ','PRAESEODYMIUM      ',
     &'NEODYMIUM          ','PROMETHIUM         ','SAMARIUM           ',
     &'EUROPIUM           ','GADOLINIUM         ','TERBIUM            ',
     &'DYSPROSIUM         ','HOLMIUM            ','ERBIUM             ',
     &'THULIUM            ','YTTERBIUM          ','LUTETIUM           ',
     &'=N-O AXOXY (DELOC) ','-N= AZOXY (DELOC)  ','>N-OH HYDROXYAMINE ',
     &'>N-OH HYDROXYAMINE ','NOT USED           ','-O- ANHYDRIDE (LOC)',
     &'-O- ANHYDRIDE(DELO)','NSP3 HYDRAZINE     ','NSP2 AMIDE (DELOC) ',
     &'NOT USED           ','>P=O PHOSPHATE     ','>SO2 SULFONAMIDE   ',
     &'NSP3 SULFONAMIDE   ','NOT USED           ','NOT USED           ',
     &'NOT USED           ','O-P=O PHOSPHATE    ','C=C NUCLEIC ACID   ',
     &'C=C NUCLEIC ACID   ','C=O NUCLEIC ACID   ','LITHIUM            ',
     &'NSP3 LI-AMIDE      '/
      data (chmtnk(i),i=1,40)/
     & 'Nonpolar Hydrogen   ', 'Aromatic Hydrogen   ',
     & 'Peptide Amide HN    ', 'Peptide HCA         ',
     & 'N-Terminal HCA      ', 'N-Terminal Hydrogen ',
     & 'N-Terminal PRO HN   ', 'Hydroxyl Hydrogen   ',
     & 'TRP Indole HE1      ', 'HIS+ Ring NH        ',
     & 'HISDE Ring NH       ', 'HIS+ HD2/HISDE HE1  ',
     & 'HIS+ HE1            ', 'HISD HD2            ',
     & 'HISE HD2            ', 'Thiol Hydrogen      ',
     & 'LYS HCE/ORN HCD     ', 'ARG HE Hydrogen     ',
     & 'ARG HZ Hydrogen     ', 'Peptide Carbonyl    ',
     & 'Aromatic Carbon     ', 'C-Term Carboxylate  ',
     & 'Peptide Alpha Carbon', 'N-Term Alpha Carbon ',
     & 'Methine Carbon      ', 'Methylene Carbon    ',
     & 'Methyl Carbon       ', 'GLY Alpha Carbon    ',
     & 'N-Terminal GLY CA   ', 'PRO CA Carbon       ',
     & 'PRO CB and CG       ', 'PRO CD Carbon       ',
     & 'N-Terminal PRO CA   ', 'N-Terminal PRO CD   ',
     & 'SER CB Carbon       ', 'THR CB Carbon       ',
     & 'Cysteine CB Carbon  ', 'Cystine CB Carbon   ',
     & 'HIS+ CB Carbon      ', 'HISD CB Carbon      '/
      data (chmtnk(i),i=41,80)/
     & 'HISE CB Carbon      ', 'HIS+ CG and CD2     ',
     & 'HISD CG/HISE CD2    ', 'HISE CG/HISD CD2    ',
     & 'HIS+ CE1 Carbon     ', 'HISDE CE1 Carbon    ',
     & 'PHE/TYR CG Carbon   ', 'TYR CZ Carbon       ',
     & 'ASN/GLN Carbonyl    ', 'TRP CG Carbon       ',
     & 'TRP CD1 Carbon      ', 'TRP CD2 Carbon      ',
     & 'TRP CE2 Carbon      ', 'ASP CB/GLU CG       ',
     & 'ASP/GLU Carboxylate ', 'MET CG Carbon       ',
     & 'MET CE Carbon       ', 'LYS CE/ORN CD       ',
     & 'ARG CD Carbon       ', 'ARG CZ Carbon       ',
     & 'N-Methyl Amide CH3  ', 'AIB Alpha Carbon    ',
     & 'Peptide Nitrogen    ', 'Amide Nitrogen      ',
     & 'Ammonium Nitrogen   ', 'PRO Nitrogen        ',
     & 'N-Terminal PRO N    ', 'HIS Ring N          ',
     & 'HIS Ring NH         ', 'HIS+ Ring NH        ',
     & 'TRP Pyrrole NE1     ', 'ARG NE Nitrogen     ',
     & 'ARG NZ Nitrogen     ', 'Peptide Oxygen      ',
     & 'ASN/GLN Carbonyl    ', 'Hydroxyl Oxygen     ',
     & 'Phenol Oxygen       ', 'ASP/GLU Carboxylate ',
     & 'C-Term Carboxylate  ', 'Thiol Sulfur        '/
      data (chmtnk(i),i=81,120)/
     & 'Sulfide Sulfur      ', 'Disulfide Sulfur    ',
     & 'Heme Meso Hydrogen  ', 'Heme H2C= Hydrogen  ',
     & 'Heme RHC= Hydrogen  ', 'Alkene RHC= Hydrogen',
     & 'Alkene H2C= Hydrogen', 'TIP3P Hydrogen      ',
     & 'COOH Hydrogen       ', 'ASP+/GLU+ Carboxyl  ',
     & 'ASP+ CB/GLU+ CG     ', 'Heme Alpha-Carbon   ',
     & 'Heme Beta-Carbon    ', 'Heme Meso-Carbon    ',
     & 'Heme Alkene Carbon  ', 'Heme CO Carbon      ',
     & 'Thiolate Carbon     ', 'Alkene RHC= Carbon  ',
     & 'Alkene H2C= Carbon  ', 'Heme Pyrrole N      ',
     & 'TIP3P Water Oxygen  ', 'COOH Hydroxyl Oxygen',
     & 'COOH Carbonyl Oxygen', 'Ester Oxygen        ',
     & 'Heme CO/O2 Oxygen   ', 'Thiolate Sulfur     ',
     & 'Heme Group Iron     ', 'Helium Atom         ',
     & 'Neon Atom           ', 'Sodium Ion          ',
     & 'Magnesium Ion       ', 'Chloride Ion        ',
     & 'Potassium Ion       ', 'Calcium Ion         ',
     & 'Zinc Ion            ', '>CH-N+ Hydrogen     ',
     & 'Methine Hydrogen    ', 'Methylene Hydrogen  ',
     & 'Methyl Hydrogen     ', '-CH=CR2 Hydrogen    '/
      data (chmtnk(i),i=121,136)/
     & 'Ester Carbonyl      ', 'Ester RCOO-CH<      ',
     & 'Ester RCOO-CH2-     ', 'Ester -CH2-COOR     ',
     & 'Methylene Carbon    ', 'Methyl Carbon       ',
     & 'Phosphate -CH2-PO4- ', 'R4N+ Methylene      ',
     & 'R4N+ Methyl         ', '-CH=CR2 Carbon      ',
     & 'R4N+ Nitrogen       ', 'Ester -O-           ',
     & 'Ester Carbonyl      ', 'Phosphate -O-       ',
     & 'Phosphate =O        ', 'Phosphate >P<       '/
      data (ambstr(i),i=1,6) /
     & 'Glycine N           ', 'Glycine CA          ',
     & 'Glycine C           ', 'Glycine HN          ',
     & 'Glycine O           ', 'Glycine HA          '/
      data (ambstr(i),i=7,14) /
     & 'Alanine N           ', 'Alanine CA          ',
     & 'Alanine C           ', 'Alanine HN          ',
     & 'Alanine O           ', 'Alanine HA          ',
     & 'Alanine CB          ', 'Alanine HB          '/
      data (ambstr(i),i=15,26) /
     & 'Valine N            ', 'Valine CA           ',
     & 'Valine C            ', 'Valine HN           ',
     & 'Valine O            ', 'Valine HA           ',
     & 'Valine CB           ', 'Valine HB           ',
     & 'Valine CG1          ', 'Valine HG1          ',
     & 'Valine CG2          ', 'Valine HG2          '/
      data (ambstr(i),i=27,40) /
     & 'Leucine N           ', 'Leucine CA          ',
     & 'Leucine C           ', 'Leucine HN          ',
     & 'Leucine O           ', 'Leucine HA          ',
     & 'Leucine CB          ', 'Leucine HB          ',
     & 'Leucine CG          ', 'Leucine HG          ',
     & 'Leucine CD1         ', 'Leucine HD1         ',
     & 'Leucine CD2         ', 'Leucine HD2         '/
      data (ambstr(i),i=41,54) /
     & 'Isoleucine N        ', 'Isoleucine CA       ',
     & 'Isoleucine C        ', 'Isoleucine HN       ',
     & 'Isoleucine O        ', 'Isoleucine HA       ',
     & 'Isoleucine CB       ', 'Isoleucine HB       ',
     & 'Isoleucine CG1      ', 'Isoleucine HG1      ',
     & 'Isoleucine CG2      ', 'Isoleucine HG2      ',
     & 'Isoleucine CD       ', 'Isoleucine HD       '/
      data (ambstr(i),i=55,64) /
     & 'Serine N            ', 'Serine CA           ',
     & 'Serine C            ', 'Serine HN           ',
     & 'Serine O            ', 'Serine HA           ',
     & 'Serine CB           ', 'Serine HB           ',
     & 'Serine OG           ', 'Serine HG           '/
      data (ambstr(i),i=65,76) /
     & 'Threonine N         ', 'Threonine CA        ',
     & 'Threonine C         ', 'Threonine HN        ',
     & 'Threonine O         ', 'Threonine HA        ',
     & 'Threonine CB        ', 'Threonine HB        ',
     & 'Threonine OG1       ', 'Threonine HG1       ',
     & 'Threonine CG2       ', 'Threonine HG2       '/
      data (ambstr(i),i=77,86) /
     & 'Cysteine (-SH) N    ', 'Cysteine (-SH) CA   ',
     & 'Cysteine (-SH) C    ', 'Cysteine (-SH) HN   ',
     & 'Cysteine (-SH) O    ', 'Cysteine (-SH) HA   ',
     & 'Cysteine (-SH) CB   ', 'Cysteine (-SH) HB   ',
     & 'Cysteine (-SH) SG   ', 'Cysteine (-SH) HG   '/
      data (ambstr(i),i=87,95) /
     & 'Cystine (-SS-) N    ', 'Cystine (-SS-) CA   ',
     & 'Cystine (-SS-) C    ', 'Cystine (-SS-) HN   ',
     & 'Cystine (-SS-) O    ', 'Cystine (-SS-) HA   ',
     & 'Cystine (-SS-) CB   ', 'Cystine (-SS-) HB   ',
     & 'Cystine (-SS-) SG   '/
      data (ambstr(i),i=96,106) /
     & 'Proline N           ', 'Proline CA          ',
     & 'Proline C           ', 'Proline O           ',
     & 'Proline HA          ', 'Proline CB          ',
     & 'Proline HB          ', 'Proline CG          ',
     & 'Proline HG          ', 'Proline CD          ',
     & 'Proline HD          '/
      data (ambstr(i),i=107,121) /
     & 'Phenylalanine N     ', 'Phenylalanine CA    ',
     & 'Phenylalanine C     ', 'Phenylalanine HN    ',
     & 'Phenylalanine O     ', 'Phenylalanine HA    ',
     & 'Phenylalanine CB    ', 'Phenylalanine HB    ',
     & 'Phenylalanine CG    ', 'Phenylalanine CD    ',
     & 'Phenylalanine HD    ', 'Phenylalanine CE    ',
     & 'Phenylalanine HE    ', 'Phenylalanine CZ    ',
     & 'Phenylalanine HZ    '/
      data (ambstr(i),i=122,137) /
     & 'Tyrosine N          ', 'Tyrosine CA         ',
     & 'Tyrosine C          ', 'Tyrosine HN         ',
     & 'Tyrosine O          ', 'Tyrosine HA         ',
     & 'Tyrosine CB         ', 'Tyrosine HB         ',
     & 'Tyrosine CG         ', 'Tyrosine CD         ',
     & 'Tyrosine HD         ', 'Tyrosine CE         ',
     & 'Tyrosine HE         ', 'Tyrosine CZ         ',
     & 'Tyrosine OH         ', 'Tyrosine HH         '/
      data (ambstr(i),i=138,160) /
     & 'Tryptophan N        ', 'Tryptophan CA       ',
     & 'Tryptophan C        ', 'Tryptophan HN       ',
     & 'Tryptophan O        ', 'Tryptophan HA       ',
     & 'Tryptophan CB       ', 'Tryptophan HB       ',
     & 'Tryptophan CG       ', 'Tryptophan CD1      ',
     & 'Tryptophan HD1      ', 'Tryptophan CD2      ',
     & 'Tryptophan NE1      ', 'Tryptophan HE1      ',
     & 'Tryptophan CE2      ', 'Tryptophan CE3      ',
     & 'Tryptophan HE3      ', 'Tryptophan CZ2      ',
     & 'Tryptophan HZ2      ', 'Tryptophan CZ3      ',
     & 'Tryptophan HZ3      ', 'Tryptophan CH2      ',
     & 'Tryptophan HH2      '/
      data (ambstr(i),i=161,177) /
     & 'Histidine (+) N     ', 'Histidine (+) CA    ',
     & 'Histidine (+) C     ', 'Histidine (+) HN    ',
     & 'Histidine (+) O     ', 'Histidine (+) HA    ',
     & 'Histidine (+) CB    ', 'Histidine (+) HB    ',
     & 'Histidine (+) CG    ', 'Histidine (+) ND1   ',
     & 'Histidine (+) HD1   ', 'Histidine (+) CD2   ',
     & 'Histidine (+) HD2   ', 'Histidine (+) CE1   ',
     & 'Histidine (+) HE1   ', 'Histidine (+) NE2   ',
     & 'Histidine (+) HE2   '/
      data (ambstr(i),i=178,193) /
     & 'Histidine (HD) N    ', 'Histidine (HD) CA   ',
     & 'Histidine (HD) C    ', 'Histidine (HD) HN   ',
     & 'Histidine (HD) O    ', 'Histidine (HD) HA   ',
     & 'Histidine (HD) CB   ', 'Histidine (HD) HB   ',
     & 'Histidine (HD) CG   ', 'Histidine (HD) ND1  ',
     & 'Histidine (HD) HD1  ', 'Histidine (HD) CD2  ',
     & 'Histidine (HD) HD2  ', 'Histidine (HD) CE1  ',
     & 'Histidine (HD) HE1  ', 'Histidine (HD) NE2  '/
      data (ambstr(i),i=194,209) /
     & 'Histidine (HE) N    ', 'Histidine (HE) CA   ',
     & 'Histidine (HE) C    ', 'Histidine (HE) HN   ',
     & 'Histidine (HE) O    ', 'Histidine (HE) HA   ',
     & 'Histidine (HE) CB   ', 'Histidine (HE) HB   ',
     & 'Histidine (HE) CG   ', 'Histidine (HE) ND1  ',
     & 'Histidine (HE) CD2  ', 'Histidine (HE) HD2  ',
     & 'Histidine (HE) CE1  ', 'Histidine (HE) HE1  ',
     & 'Histidine (HE) NE2  ', 'Histidine (HE) HE2  '/
      data (ambstr(i),i=210,219) /
     & 'Aspartic Acid N     ', 'Aspartic Acid CA    ',
     & 'Aspartic Acid C     ', 'Aspartic Acid HN    ',
     & 'Aspartic Acid O     ', 'Aspartic Acid HA    ',
     & 'Aspartic Acid CB    ', 'Aspartic Acid HB    ',
     & 'Aspartic Acid CG    ', 'Aspartic Acid OD    '/
      data (ambstr(i),i=220,231) /
     & 'Asparagine N        ', 'Asparagine CA       ',
     & 'Asparagine C        ', 'Asparagine HN       ',
     & 'Asparagine O        ', 'Asparagine HA       ',
     & 'Asparagine CB       ', 'Asparagine HB       ',
     & 'Asparagine CG       ', 'Asparagine OD1      ',
     & 'Asparagine ND2      ', 'Asparagine HD2      '/
      data (ambstr(i),i=232,243) /
     & 'Glutamic Acid N     ', 'Glutamic Acid CA    ',
     & 'Glutamic Acid C     ', 'Glutamic Acid HN    ',
     & 'Glutamic Acid O     ', 'Glutamic Acid HA    ',
     & 'Glutamic Acid CB    ', 'Glutamic Acid HB    ',
     & 'Glutamic Acid CG    ', 'Glutamic Acid HG    ',
     & 'Glutamic Acid CD    ', 'Glutamic Acid OE    '/
      data (ambstr(i),i=244,257) /
     & 'Glutamine N         ', 'Glutamine CA        ',
     & 'Glutamine C         ', 'Glutamine HN        ',
     & 'Glutamine O         ', 'Glutamine HA        ',
     & 'Glutamine CB        ', 'Glutamine HB        ',
     & 'Glutamine CG        ', 'Glutamine HG        ',
     & 'Glutamine CD        ', 'Glutamine OE1       ',
     & 'Glutamine NE2       ', 'Glutamine HE2       '/
      data (ambstr(i),i=258,270) /
     & 'Methionine N        ', 'Methionine CA       ',
     & 'Methionine C        ', 'Methionine HN       ',
     & 'Methionine O        ', 'Methionine HA       ',
     & 'Methionine CB       ', 'Methionine HB       ',
     & 'Methionine CG       ', 'Methionine HG       ',
     & 'Methionine SD       ', 'Methionine CE       ',
     & 'Methionine HE       '/
      data (ambstr(i),i=271,286) /
     & 'Lysine N            ', 'Lysine CA           ',
     & 'Lysine C            ', 'Lysine HN           ',
     & 'Lysine O            ', 'Lysine HA           ',
     & 'Lysine CB           ', 'Lysine HB           ',
     & 'Lysine CG           ', 'Lysine HG           ',
     & 'Lysine CD           ', 'Lysine HD           ',
     & 'Lysine CE           ', 'Lysine HE           ',
     & 'Lysine NZ           ', 'Lysine HZ           '/
      data (ambstr(i),i=287,303) /
     & 'Arginine N          ', 'Arginine CA         ',
     & 'Arginine C          ', 'Arginine HN         ',
     & 'Arginine O          ', 'Arginine HA         ',
     & 'Arginine CB         ', 'Arginine HB         ',
     & 'Arginine CG         ', 'Arginine HG         ',
     & 'Arginine CD         ', 'Arginine HD         ',
     & 'Arginine NE         ', 'Arginine HE         ',
     & 'Arginine CZ         ', 'Arginine NH         ',
     & 'Arginine HH         '/
      data (ambstr(i),i=304,317) /
     & 'Ornithine N         ', 'Ornithine CA        ',
     & 'Ornithine C         ', 'Ornithine HN        ',
     & 'Ornithine O         ', 'Ornithine HA        ',
     & 'Ornithine CB        ', 'Ornithine HB        ',
     & 'Ornithine CG        ', 'Ornithine HG        ',
     & 'Ornithine CD        ', 'Ornithine HD        ',
     & 'Ornithine NE        ', 'Ornithine HE        '/
      data (ambstr(i),i=318,324) /
     & 'MethylAlanine N     ', 'MethylAlanine CA    ',
     & 'MethylAlanine C     ', 'MethylAlanine HN    ',
     & 'MethylAlanine O     ', 'MethylAlanine CB    ',
     & 'MethylAlanine HB    '/
      data (ambstr(i),i=325,336) /
     & 'Pyroglutamate N     ', 'Pyroglutamate CA    ',
     & 'Pyroglutamate C     ', 'Pyroglutamate HN    ',
     & 'Pyroglutamate O     ', 'Pyroglutamate HA    ',
     & 'Pyroglutamate CB    ', 'Pyroglutamate HB    ',
     & 'Pyroglutamate CG    ', 'Pyroglutamate HG    ',
     & 'Pyroglutamate CD    ', 'Pyroglutamate OE    '/
      data (ambstr(i),i=337,349) /
     & 'Formyl C            ', 'Formyl H            ',
     & 'Formyl O            ', 'Acetyl CA           ',
     & 'Acetyl HA           ', 'Acetyl C            ',
     & 'Acetyl O            ', 'C-Term Amide N      ',
     & 'C-Term Amide HN     ', 'N-MeAmide N         ',
     & 'N-MeAmide HN        ', 'N-MeAmide C         ',
     & 'N-MeAmide HC        '/
      data (ambstr(i),i=350,355) /
     & 'N-Term GLY N        ', 'N-Term GLY CA       ',
     & 'N-Term GLY C        ', 'N-Term GLY HN       ',
     & 'N-Term GLY O        ', 'N-Term GLY HA       '/
      data (ambstr(i),i=356,361) /
     & 'N-Term ALA N        ', 'N-Term ALA CA       ',
     & 'N-Term ALA C        ', 'N-Term ALA HN       ',
     & 'N-Term ALA O        ', 'N-Term ALA HA       '/
      data (ambstr(i),i=362,367) /
     & 'N-Term VAL N        ', 'N-Term VAL CA       ',
     & 'N-Term VAL C        ', 'N-Term VAL HN       ',
     & 'N-Term VAL O        ', 'N-Term VAL HA       '/
      data (ambstr(i),i=368,373) /
     & 'N-Term LEU N        ', 'N-Term LEU CA       ',
     & 'N-Term LEU C        ', 'N-Term LEU HN       ',
     & 'N-Term LEU O        ', 'N-Term LEU HA       '/
      data (ambstr(i),i=374,379) /
     & 'N-Term ILE N        ', 'N-Term ILE CA       ',
     & 'N-Term ILE C        ', 'N-Term ILE HN       ',
     & 'N-Term ILE O        ', 'N-Term ILE HA       '/
      data (ambstr(i),i=380,385) /
     & 'N-Term SER N        ', 'N-Term SER CA       ',
     & 'N-Term SER C        ', 'N-Term SER HN       ',
     & 'N-Term SER O        ', 'N-Term SER HA       '/
      data (ambstr(i),i=386,391) /
     & 'N-Term THR N        ', 'N-Term THR CA       ',
     & 'N-Term THR C        ', 'N-Term THR HN       ',
     & 'N-Term THR O        ', 'N-Term THR HA       '/
      data (ambstr(i),i=392,397) /
     & 'N-Term CYS (-SH) N  ', 'N-Term CYS (-SH) CA ',
     & 'N-Term CYS (-SH) C  ', 'N-Term CYS (-SH) HN ',
     & 'N-Term CYS (-SH) O  ', 'N-Term CYS (-SH) HA '/
      data (ambstr(i),i=398,403) /
     & 'N-Term CYS (-SS-) N ', 'N-Term CYS (-SS-) CA',
     & 'N-Term CYS (-SS-) C ', 'N-Term CYS (-SS-) HN',
     & 'N-Term CYS (-SS-) O ', 'N-Term CYS (-SS-) HA'/
      data (ambstr(i),i=404,411) /
     & 'N-Term PRO N        ', 'N-Term PRO CA       ',
     & 'N-Term PRO C        ', 'N-Term PRO HN       ',
     & 'N-Term PRO O        ', 'N-Term PRO HA       ',
     & 'N-Term PRO CD       ', 'N-Term PRO HD       '/
      data (ambstr(i),i=412,417) /
     & 'N-Term PHE N        ', 'N-Term PHE CA       ',
     & 'N-Term PHE C        ', 'N-Term PHE HN       ',
     & 'N-Term PHE O        ', 'N-Term PHE HA       '/
      data (ambstr(i),i=418,423) /
     & 'N-Term TYR N        ', 'N-Term TYR CA       ',
     & 'N-Term TYR C        ', 'N-Term TYR HN       ',
     & 'N-Term TYR O        ', 'N-Term TYR HA       '/
      data (ambstr(i),i=424,429) /
     & 'N-Term TRP N        ', 'N-Term TRP CA       ',
     & 'N-Term TRP C        ', 'N-Term TRP HN       ',
     & 'N-Term TRP O        ', 'N-Term TRP HA       '/
      data (ambstr(i),i=430,435) /
     & 'N-Term HIS (+) N    ', 'N-Term HIS (+) CA   ',
     & 'N-Term HIS (+) C    ', 'N-Term HIS (+) HN   ',
     & 'N-Term HIS (+) O    ', 'N-Term HIS (+) HA   '/
      data (ambstr(i),i=436,441) /
     & 'N-Term HIS (HD) N   ', 'N-Term HIS (HD) CA  ',
     & 'N-Term HIS (HD) C   ', 'N-Term HIS (HD) HN  ',
     & 'N-Term HIS (HD) O   ', 'N-Term HIS (HD) HA  '/
      data (ambstr(i),i=442,447) /
     & 'N-Term HIS (HE) N   ', 'N-Term HIS (HE) CA  ',
     & 'N-Term HIS (HE) C   ', 'N-Term HIS (HE) HN  ',
     & 'N-Term HIS (HE) O   ', 'N-Term HIS (HE) HA  '/
      data (ambstr(i),i=448,453) /
     & 'N-Term ASP N        ', 'N-Term ASP CA       ',
     & 'N-Term ASP C        ', 'N-Term ASP HN       ',
     & 'N-Term ASP O        ', 'N-Term ASP HA       '/
      data (ambstr(i),i=454,459) /
     & 'N-Term ASN N        ', 'N-Term ASN CA       ',
     & 'N-Term ASN C        ', 'N-Term ASN HN       ',
     & 'N-Term ASN O        ', 'N-Term ASN HA       '/
      data (ambstr(i),i=460,465) /
     & 'N-Term GLU N        ', 'N-Term GLU CA       ',
     & 'N-Term GLU C        ', 'N-Term GLU HN       ',
     & 'N-Term GLU O        ', 'N-Term GLU HA       '/
      data (ambstr(i),i=466,471) /
     & 'N-Term GLN N        ', 'N-Term GLN CA       ',
     & 'N-Term GLN C        ', 'N-Term GLN HN       ',
     & 'N-Term GLN O        ', 'N-Term GLN HA       '/
      data (ambstr(i),i=472,477) /
     & 'N-Term MET N        ', 'N-Term MET CA       ',
     & 'N-Term MET C        ', 'N-Term MET HN       ',
     & 'N-Term MET O        ', 'N-Term MET HA       '/
      data (ambstr(i),i=478,483) /
     & 'N-Term LYS N        ', 'N-Term LYS CA       ',
     & 'N-Term LYS C        ', 'N-Term LYS HN       ',
     & 'N-Term LYS O        ', 'N-Term LYS HA       '/
      data (ambstr(i),i=484,489) /
     & 'N-Term ARG N        ', 'N-Term ARG CA       ',
     & 'N-Term ARG C        ', 'N-Term ARG HN       ',
     & 'N-Term ARG O        ', 'N-Term ARG HA       '/
      data (ambstr(i),i=490,495) /
     & 'N-Term ORN N        ', 'N-Term ORN CA       ',
     & 'N-Term ORN C        ', 'N-Term ORN HN       ',
     & 'N-Term ORN O        ', 'N-Term ORN HA       '/
      data (ambstr(i),i=496,500) /
     & 'N-Term AIB N        ', 'N-Term AIB CA       ',
     & 'N-Term AIB C        ', 'N-Term AIB HN       ',
     & 'N-Term AIB O        '/
      data (ambstr(i),i=501,506) /
     & 'C-Term GLY N        ', 'C-Term GLY CA       ',
     & 'C-Term GLY C        ', 'C-Term GLY HN       ',
     & 'C-Term GLY OXT      ', 'C-Term GLY HA       '/
      data (ambstr(i),i=507,512) /
     & 'C-Term ALA N        ', 'C-Term ALA CA       ',
     & 'C-Term ALA C        ', 'C-Term ALA HN       ',
     & 'C-Term ALA OXT      ', 'C-Term ALA HA       '/
      data (ambstr(i),i=513,518) /
     & 'C-Term VAL N        ', 'C-Term VAL CA       ',
     & 'C-Term VAL C        ', 'C-Term VAL HN       ',
     & 'C-Term VAL OXT      ', 'C-Term VAL HA       '/
      data (ambstr(i),i=519,524) /
     & 'C-Term LEU N        ', 'C-Term LEU CA       ',
     & 'C-Term LEU C        ', 'C-Term LEU HN       ',
     & 'C-Term LEU OXT      ', 'C-Term LEU HA       '/
      data (ambstr(i),i=525,530) /
     & 'C-Term ILE N        ', 'C-Term ILE CA       ',
     & 'C-Term ILE C        ', 'C-Term ILE HN       ',
     & 'C-Term ILE OXT      ', 'C-Term ILE HA       '/
      data (ambstr(i),i=531,536) /
     & 'C-Term SER N        ', 'C-Term SER CA       ',
     & 'C-Term SER C        ', 'C-Term SER HN       ',
     & 'C-Term SER OXT      ', 'C-Term SER HA       '/
      data (ambstr(i),i=537,542) /
     & 'C-Term THR N        ', 'C-Term THR CA       ',
     & 'C-Term THR C        ', 'C-Term THR HN       ',
     & 'C-Term THR OXT      ', 'C-Term THR HA       '/
      data (ambstr(i),i=543,548) /
     & 'C-Term CYS (-SH) N  ', 'C-Term CYS (-SH) CA ',
     & 'C-Term CYS (-SH) C  ', 'C-Term CYS (-SH) HN ',
     & 'C-Term CYS (-SH) OXT', 'C-Term CYS (-SH) HA '/
      data (ambstr(i),i=549,554) /
     & 'C-Term CYS (-SS-) N ', 'C-Term CYS (-SS-) CA',
     & 'C-Term CYS (-SS-) C ', 'C-Term CYS (-SS-) HN',
     & 'C-Term CYS (-SS-)OXT', 'C-Term CYS (-SS-) HA'/
       data (ambstr(i),i=555,559) /
     & 'C-Term PRO N        ', 'C-Term PRO CA       ',
     & 'C-Term PRO C        ', 'C-Term PRO OXT      ',
     & 'C-Term PRO HA       '/
      data (ambstr(i),i=560,565) /
     & 'C-Term PHE N        ', 'C-Term PHE CA       ',
     & 'C-Term PHE C        ', 'C-Term PHE HN       ',
     & 'C-Term PHE OXT      ', 'C-Term PHE HA       '/
      data (ambstr(i),i=566,571) /
     & 'C-Term TYR N        ', 'C-Term TYR CA       ',
     & 'C-Term TYR C        ', 'C-Term TYR HN       ',
     & 'C-Term TYR OXT      ', 'C-Term TYR HA       '/
      data (ambstr(i),i=572,577) /
     & 'C-Term TRP N        ', 'C-Term TRP CA       ',
     & 'C-Term TRP C        ', 'C-Term TRP HN       ',
     & 'C-Term TRP OXT      ', 'C-Term TRP HA       '/
      data (ambstr(i),i=578,583) /
     & 'C-Term HIS (+) N    ', 'C-Term (+) HIS CA   ',
     & 'C-Term HIS (+) C    ', 'C-Term (+) HIS HN   ',
     & 'C-Term HIS (+) OXT  ', 'C-Term (+) HIS HA   '/
      data (ambstr(i),i=584,589) /
     & 'C-Term HIS (HD) N   ', 'C-Term HIS (HD) CA  ',
     & 'C-Term HIS (HD) C   ', 'C-Term HIS (HD) HN  ',
     & 'C-Term HIS (HD) OXT ', 'C-Term HIS (HD) HA  '/
      data (ambstr(i),i=590,595) /
     & 'C-Term HIS (HE) N   ', 'C-Term HIS (HE) CA  ',
     & 'C-Term HIS (HE) C   ', 'C-Term HIS (HE) HN  ',
     & 'C-Term HIS (HE) OXT ', 'C-Term HIS (HE) HA  '/
      data (ambstr(i),i=596,601) /
     & 'C-Term ASP N        ', 'C-Term ASP CA       ',
     & 'C-Term ASP C        ', 'C-Term ASP HN       ',
     & 'C-Term ASP OXT      ', 'C-Term ASP HA       '/
      data (ambstr(i),i=602,607) /
     & 'C-Term ASN N        ', 'C-Term ASN CA       ',
     & 'C-Term ASN C        ', 'C-Term ASN HN       ',
     & 'C-Term ASN OXT      ', 'C-Term ASN HA       '/
      data (ambstr(i),i=608,613) /
     & 'C-Term GLU N        ', 'C-Term GLU CA       ',
     & 'C-Term GLU C        ', 'C-Term GLU HN       ',
     & 'C-Term GLU OXT      ', 'C-Term GLU HA       '/
      data (ambstr(i),i=614,619) /
     & 'C-Term GLN N        ', 'C-Term GLN CA       ',
     & 'C-Term GLN C        ', 'C-Term GLN HN       ',
     & 'C-Term GLN OXT      ', 'C-Term GLN HA       '/
      data (ambstr(i),i=620,625) /
     & 'C-Term MET N        ', 'C-Term MET CA       ',
     & 'C-Term MET C        ', 'C-Term MET HN       ',
     & 'C-Term MET OXT      ', 'C-Term MET HA       '/
      data (ambstr(i),i=626,631) /
     & 'C-Term LYS N        ', 'C-Term LYS CA       ',
     & 'C-Term LYS C        ', 'C-Term LYS HN       ',
     & 'C-Term LYS OXT      ', 'C-Term LYS HA       '/
      data (ambstr(i),i=632,637) /
     & 'C-Term ARG N        ', 'C-Term ARG CA       ',
     & 'C-Term ARG C        ', 'C-Term ARG HN       ',
     & 'C-Term ARG OXT      ', 'C-Term ARG HA       '/
      data (ambstr(i),i=638,643) /
     & 'C-Term ORN N        ', 'C-Term ORN CA       ',
     & 'C-Term ORN C        ', 'C-Term ORN HN       ',
     & 'C-Term ORN OXT      ', 'C-Term ORN HA       '/
      data (ambstr(i),i=644,648) /
     & 'C-Term AIB N        ', 'C-Term AIB CA       ',
     & 'C-Term AIB C        ', 'C-Term AIB HN       ',
     & 'C-Term AIB OXT      '/
      data (ambstr(i),i=649,659) /
     & 'TIP3P Oxygen        ','TIP3P Hydrogen      ',
     & 'Li+ Lithium Ion     ','Na+ Sodium Ion      ',
     & 'K+ Potassium Ion    ','Rb+ Rubidium Ion    ',
     & 'Cs+ Cesium Ion      ','Mg+2 Magnesium Ion  ',
     & 'Ca+2 Calcium Ion    ','Zn+2 Zinc Ion       ',
     & 'Cl-  Chlorine Ion   '/
c some neutral residues
      data (ambstr(i),i=660,671) /
     & 'ASP Neutral N	    ', 'ASP Neutral CA	    ',
     & 'ASP Neutral C       ', 'ASP Neutral HN      ',
     & 'ASP Neutral O       ', 'ASP Neutral HA      ',
     & 'ASP Neutral CB	    ', 'ASP Neutral HB      ',
     & 'ASP Neutral CG	    ', 'ASP Neutral OD1     ',
     & 'ASP Neutral OD2	    ', 'ASP Neutral HD2     '/
      data (ambstr(i),i=672,685) /
     & 'GLU Neutral N       ', 'GLU Neutral CA      ',
     & 'GLU Neutral C       ', 'GLU Neutral HN      ',
     & 'GLU Neutral O       ', 'GLU Neutral HA      ',
     & 'GLU Neutral CB      ', 'GLU Neutral HB      ',
     & 'GLU Neutral CG      ', 'GLU Neutral HG      ',
     & 'GLU Neutral CD      ', 'GLU Neutral OE1     ',
     & 'GLU Neutral OE2     ', 'GLU Neutral HE2     '/
      data (ambstr(i),i=686,701) /
     & 'LYS Neutral N       ', 'LYS Neutral CA      ',
     & 'LYS Neutral C       ', 'LYS Neutral HN      ',
     & 'LYS Neutral O       ', 'LYS Neutral HA      ',
     & 'LYS Neutral CB      ', 'LYS Neutral HB      ',
     & 'LYS Neutral CG      ', 'LYS Neutral HG      ',
     & 'LYS Neutral CD      ', 'LYS Neutral HD      ',
     & 'LYS Neutral CE      ', 'LYS Neutral HE      ',	
     & 'LYS Neutral NZ      ', 'LYS Neutral HZ      '/

      data (ambstr(i),i=702,1000) / 299*"                    "/

      data (ambstr(i),i=1001,1030) /
     & "R-Adenosine O5'     ", "R-Adenosine C5'     ",
     & "R-Adenosine H5'1    ", "R-Adenosine H5'2    ",
     & "R-Adenosine C4'     ", "R-Adenosine H4'     ",
     & "R-Adenosine O4'     ", "R-Adenosine C1'     ",
     & "R-Adenosine H1'     ", "R-Adenosine C3'     ",
     & "R-Adenosine H3'     ", "R-Adenosine C2'     ",
     & "R-Adenosine H2'1    ", "R-Adenosine O2'     ",
     & "R-Adenosine HO'2    ", "R-Adenosine O3'     ",
     & "R-Adenosine N9      ", "R-Adenosine C4      ",
     & "R-Adenosine C5      ", "R-Adenosine N7      ",
     & "R-Adenosine C8      ", "R-Adenosine N3      ",
     & "R-Adenosine C2      ", "R-Adenosine N1      ",
     & "R-Adenosine C6      ", "R-Adenosine H2      ",
     & "R-Adenosine N6      ", "R-Adenosine H61     ",
     & "R-Adenosine H62     ", "R-Adenosine H8      "/

      data (ambstr(i),i=1031,1061) /
     & "R-Guanosine O5'     ", "R-Guanosine C5'     ",
     & "R-Guanosine H5'1    ", "R-Guanosine H5'2    ",
     & "R-Guanosine C4'     ", "R-Guanosine H4'     ",
     & "R-Guanosine O4'     ", "R-Guanosine C1'     ",
     & "R-Guanosine H1'     ", "R-Guanosine C3'     ",
     & "R-Guanosine H3'     ", "R-Guanosine C2'     ",
     & "R-Guanosine H2'1    ", "R-Guanosine O2'     ",
     & "R-Guanosine HO'2    ", "R-Guanosine O3'     ",
     & "R-Guanosine N9      ", "R-Guanosine C4      ",
     & "R-Guanosine C5      ", "R-Guanosine N7      ",
     & "R-Guanosine C8      ", "R-Guanosine N3      ",
     & "R-Guanosine C2      ", "R-Guanosine N1      ",
     & "R-Guanosine C6      ", "R-Guanosine H1      ",
     & "R-Guanosine N2      ", "R-Guanosine H21     ",
     & "R-Guanosine H22     ", "R-Guanosine O6      ",
     & "R-Guanosine H8      "/

      data (ambstr(i),i=1062,1089) /
     & "R-Cytosine O5'      ", "R-Cytosine C5'      ",
     & "R-Cytosine H5'1     ", "R-Cytosine H5'2     ",
     & "R-Cytosine C4'      ", "R-Cytosine H4'      ",
     & "R-Cytosine O4'      ", "R-Cytosine C1'      ",
     & "R-Cytosine H1'      ", "R-Cytosine C3'      ",
     & "R-Cytosine H3'      ", "R-Cytosine C2'      ",
     & "R-Cytosine H2'1     ", "R-Cytosine O2'      ",
     & "R-Cytosine HO'2     ", "R-Cytosine O3'      ",
     & "R-Cytosine N1       ", "R-Cytosine C2       ",
     & "R-Cytosine N3       ", "R-Cytosine C4       ",
     & "R-Cytosine C5       ", "R-Cytosine C6       ",
     & "R-Cytosine O2       ", "R-Cytosine N4       ",
     & "R-Cytosine H41      ", "R-Cytosine H42      ",
     & "R-Cytosine H5       ", "R-Cytosine H6       "/

      data (ambstr(i),i=1090,1116) /
     & "R-Uracil O5'        ", "R-Uracil C5'        ",
     & "R-Uracil H5'1       ", "R-Uracil H5'2       ",
     & "R-Uracil C4'        ", "R-Uracil H4'        ",
     & "R-Uracil O4'        ", "R-Uracil C1'        ",
     & "R-Uracil H1'        ", "R-Uracil C3'        ",
     & "R-Uracil H3'        ", "R-Uracil C2'        ",
     & "R-Uracil H2'1       ", "R-Uracil O2'        ",
     & "R-Uracil HO'2       ", "R-Uracil O3'        ",
     & "R-Uracil N1         ", "R-Uracil C2         ",
     & "R-Uracil N3         ", "R-Uracil C4         ",
     & "R-Uracil C5         ", "R-Uracil C6         ",
     & "R-Uracil O2         ", "R-Uracil H3         ",
     & "R-Uracil O4         ", "R-Uracil H5         ",
     & "R-Uracil H6         "/

      data (ambstr(i),i=1117,1145) /
     & "D-Adenosine O5'     ", "D-Adenosine C5'     ",
     & "D-Adenosine H5'1    ", "D-Adenosine H5'2    ",
     & "D-Adenosine C4'     ", "D-Adenosine H4'     ",
     & "D-Adenosine O4'     ", "D-Adenosine C1'     ",
     & "D-Adenosine H1'     ", "D-Adenosine C3'     ",
     & "D-Adenosine H3'     ", "D-Adenosine C2'     ",
     & "D-Adenosine H2'1    ", "D-Adenosine H2'2    ",
     & "D-Adenosine O3'     ", "D-Adenosine N9      ",
     & "D-Adenosine C4      ", "D-Adenosine C5      ",
     & "D-Adenosine N7      ", "D-Adenosine C8      ",
     & "D-Adenosine N3      ", "D-Adenosine C2      ",
     & "D-Adenosine N1      ", "D-Adenosine C6      ",
     & "D-Adenosine H2      ", "D-Adenosine N6      ",
     & "D-Adenosine H61     ", "D-Adenosine H62     ",
     & "D-Adenosine H8      "/

      data (ambstr(i),i=1146,1175) /
     & "D-Guanosine O5'     ", "D-Guanosine C5'     ",
     & "D-Guanosine H5'1    ", "D-Guanosine H5'2    ",
     & "D-Guanosine C4'     ", "D-Guanosine H4'     ",
     & "D-Guanosine O4'     ", "D-Guanosine C1'     ",
     & "D-Guanosine H1'     ", "D-Guanosine C3'     ",
     & "D-Guanosine H3'     ", "D-Guanosine C2'     ",
     & "D-Guanosine H2'1    ", "D-Guanosine H2'2    ",
     & "D-Guanosine O3'     ", "D-Guanosine N9      ",
     & "D-Guanosine C4      ", "D-Guanosine C5      ",
     & "D-Guanosine N7      ", "D-Guanosine C8      ",
     & "D-Guanosine N3      ", "D-Guanosine C2      ",
     & "D-Guanosine N1      ", "D-Guanosine C6      ",
     & "D-Guanosine H1      ", "D-Guanosine N2      ",
     & "D-Guanosine H21     ", "D-Guanosine H22     ",
     & "D-Guanosine O6      ", "D-Guanosine H8      "/

      data (ambstr(i),i=1176,1202) /
     & "D-Cytosine O5'      ", "D-Cytosine C5'      ",
     & "D-Cytosine H5'1     ", "D-Cytosine H5'2     ",
     & "D-Cytosine C4'      ", "D-Cytosine H4'      ",
     & "D-Cytosine O4'      ", "D-Cytosine C1'      ",
     & "D-Cytosine H1'      ", "D-Cytosine C3'      ",
     & "D-Cytosine H3'      ", "D-Cytosine C2'      ",
     & "D-Cytosine H2'1     ", "D-Cytosine H2'2     ",
     & "D-Cytosine O3'      ", "D-Cytosine N1       ",
     & "D-Cytosine C2       ", "D-Cytosine N3       ",
     & "D-Cytosine C4       ", "D-Cytosine C5       ",
     & "D-Cytosine C6       ", "D-Cytosine O2       ",
     & "D-Cytosine N4       ", "D-Cytosine H41      ",
     & "D-Cytosine H42      ", "D-Cytosine H5       ",
     & "D-Cytosine H6       "/

      data (ambstr(i),i=1203,1229) /
     & "D-Thymine O5'       ", "D-Thymine C5'       ",
     & "D-Thymine H5'1      ", "D-Thymine H5'2      ",
     & "D-Thymine C4'       ", "D-Thymine H4'       ",
     & "D-Thymine O4'       ", "D-Thymine C1'       ",
     & "D-Thymine H1'       ", "D-Thymine C3'       ",
     & "D-Thymine H3'       ", "D-Thymine C2'       ",
     & "D-Thymine H2'1      ", "D-Thymine H2'2      ",
     & "D-Thymine O3'       ", "D-Thymine N1        ",
     & "D-Thymine C2        ", "D-Thymine N3        ",
     & "D-Thymine C4        ", "D-Thymine C5        ",
     & "D-Thymine C6        ", "D-Thymine O2        ",
     & "D-Thymine H3        ", "D-Thymine O4        ",
     & "D-Thymine C7        ", "D-Thymine H7        ",
     & "D-Thymine H6        "/

      data (ambstr(i),i=1230,1241) /
     & "R-Phosphodiester P  ", "R-Phosphodiester OP ",
     & "R-5'-Hydroxyl O5'   ", "R-5'-Hydroxyl H5T   ",
     & "R-5'-Phosphate O5'  ", "R-5'-Phosphate P    ",
     & "R-5'-Phosphate OP   ", "R-3'-Hydroxyl O3'   ",
     & "R-3'-Hydroxyl H3T   ", "R-3'-Phosphate O3'  ",
     & "R-3'-Phosphate P    ", "R-3'-Phosphate OP   "/

      data (ambstr(i),i=1242,1253) /
     & "D-Phosphodiester P  ", "D-Phosphodiester OP ",
     & "D-5'-Hydroxyl O5'   ", "D-5'-Hydroxyl H5T   ",
     & "D-5'-Phosphate O5'  ", "D-5'-Phosphate P    ",
     & "D-5'-Phosphate OP   ", "D-3'-Hydroxyl O3'   ",
     & "D-3'-Hydroxyl H3T   ", "D-3'-Phosphate O3'  ",
     & "D-3'-Phosphate P    ", "D-3'-Phosphate OP   "/

      data (ambstr(i),i=1254,1284) /
     & "1MA         O5'     ", "1MA         C5'     ",
     & "1MA         H5'1    ", "1MA         H5'2    ",
     & "1MA         C4'     ", "1MA         H4'     ",
     & "1MA         O4'     ", "1MA         C1'     ",
     & "1MA         H1'     ", "1MA         C3'     ",
     & "1MA         H3'     ", "1MA         C2'     ",
     & "1MA         H2'1    ", "1MA         O2'     ",
     & "1MA         HO'2    ", "1MA         O3'     ",
     & "1MA         N9      ", "1MA         C4      ",
     & "1MA         C5      ", "1MA         N7      ",
     & "1MA         C8      ", "1MA         N3      ",
     & "1MA         C2      ", "1MA         N1      ",
     & "1MA         C6      ", "1MA         H2      ",
     & "1MA         N6      ", "1MA         HN6     ",
     & "1MA         H8      ", "1MA         CM1     ", 
     & "1MA         HM1     "/

      data (ambchg(i),i=1254,1284) /
     & -0.4989, 0.0558, 0.0679, 0.0679, 0.1065, 0.1174,-0.3548,
     &  0.0499, 0.1524, 0.2022, 0.0615, 0.0670, 0.0972,-0.6139,
     &  0.4186,-0.5246, 0.0097, 0.3796,-0.0003,-0.5030, 0.0832,
     & -0.6579, 0.3242,-0.2368, 0.6002, 0.1277,-0.8601, 0.3794,
     &  0.1608,-0.2077, 0.1049/

      data (ambtnk(i),i=1254,1284) /
     & 'OS','CT','H1','H1','CT','H1','OS','CT','H2','CT','H1',
     & 'CT','H1','OH','HO','OS','N*','CB','CB','NB','CK','NC',
     & 'CQ','N*','C ','H1','N ','H ','H5','CT','H1'/

      data (ambvdt(i),i=1254,1284) /
     & 23,1,35,35,1,35,23,1,36,1,35,
     & 1,35,22,31,23,18,9,9,16,12,17,
     & 13,18,2,35,14,29,40,1,35/
  
      data (ambstr(i),i=1285,1311) /
     & "5MC        O5'      ", "5MC        C5'      ",
     & "5MC        H5'      ", "5MC        C4'      ", 
     & "5MC        H4'      ", "5MC        O4'      ", 
     & "5MC        C1'      ", "5MC        H1'      ", 
     & "5MC        C3'      ", "5MC        H3'      ", 
     & "5MC        C2'      ", "5MC        H2'1     ", 
     & "5MC        O2'      ", "5MC        HO'2     ", 
     & "5MC        O3'      ", "5MC        N1       ", 
     & "5MC        C2       ", "5MC        N3       ", 
     & "5MC        C4       ", "5MC        C5       ", 
     & "5MC        C6       ", "5MC        H6       ", 
     & "5MC        O2       ", "5MC        N4       ", 
     & "5MC        HN4      ", "5MC        CM5      ", 
     & "5MC        HM5      "/

      data (ambchg(i),i=1285,1311) /
     & -0.4989,  0.0558, 0.0679, 0.1065, 0.1174, -0.3548,
     &  0.0361,  0.1824, 0.2022, 0.0615, 0.0670,  0.0972,
     & -0.6139,  0.4186,-0.5246,-0.1305, 0.8758, -0.8272,
     &  0.8361, -0.2290,-0.0694, 0.1986,-0.6431, -1.0421,
     &  0.4509, -0.2096, 0.0787/

      data (ambtnk(i),i=1285,1311) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'N*', 'C ', 'NC',
     & 'CA', 'CM', 'CM', 'H4', 'O ', 'N2', 'H ', 'CT', 'HC'/

      data (ambvdt(i),i=1285,1311) /
     & 23,1,35,1,35,23,1,36,1,
     & 35,1,35,22,31,23,18,2,17,
     & 3,4,4,39,24,19,29,1,34/

      data (ambstr(i),i=1312,1338) /
     & "OMC        O5'      ", "OMC        C5'      ",
     & "OMC        H5'      ", "OMC        C4'      ", 
     & "OMC        H4'      ", "OMC        O4'      ", 
     & "OMC        C1'      ", "OMC        H1'      ", 
     & "OMC        C3'      ", "OMC        H3'      ", 
     & "OMC        C2'      ", "OMC        H2'1     ", 
     & "OMC        O2'      ", "OMC        O3'      ", 
     & "OMC        N1       ", "OMC        C2       ", 
     & "OMC        O2       ", "OMC        N3       ", 
     & "OMC        C4       ", "OMC        N4       ", 
     & "OMC        HN4      ", "OMC        C5       ", 
     & "OMC        H5       ", "OMC        C6       ", 
     & "OMC        H6       ", "OMC        CM2      ", 
     & "OMC        HM2      "/

      data (ambchg(i),i=1312,1338) /
     & -0.4989, 0.0558, 0.0679, 0.1065, 0.1174, -0.3548,
     &  0.0026, 0.2029, 0.2022, 0.0615, 0.0847,  0.0490,
     & -0.3668,-0.5246,-0.0484, 0.7538,-0.6252, -0.7584,
     &  0.8185,-0.9530, 0.4234,-0.5215, 0.1928,  0.0053,
     &  0.1958,-0.0823, 0.0961/

      data (ambtnk(i),i=1312,1338) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OS', 'OS', 'N*', 'C ', 'O ', 'NC',
     & 'CA', 'N2', 'H ', 'CM', 'HA', 'CM', 'H4', 'CT', 'H1'/
 
      data (ambvdt(i),i=1312,1338) /
     & 23,1,35,1,35,23,1,36,1,
     & 35,1,35,23,23,18,2,24,17,
     & 3,19,29,4,33,4,39,1,35/


c H5, two protons, HM1 doubles for 3 H

      data (ambstr(i),i=1339,1369) /
     & "2MG         O5'     ", "2MG         C5'     ",
     & "2MG         H5'     ", "2MG         C4'     ", 
     & "2MG         H4'     ", "2MG         O4'     ", 
     & "2MG         C1'     ", "2MG         H1'     ", 
     & "2MG         C3'     ", "2MG         H3'     ", 
     & "2MG         C2'     ", "2MG         H2'1    ", 
     & "2MG         O2'     ", "2MG         HO'2    ", 
     & "2MG         O3'     ", "2MG         N9      ", 
     & "2MG         C4      ", "2MG         C5      ", 
     & "2MG         N7      ", "2MG         C8      ", 
     & "2MG         H8      ", "2MG         N3      ", 
     & "2MG         C2      ", "2MG         N1      ", 
     & "2MG         HN1     ", "2MG         C6      ",
     & "2MG         O6      ", "2MG         N2      ", 
     & "2MG         HN2     ", "2MG         CM2     ", 
     & "2MG         HM2     "/

      data (ambchg(i),i=1339,1369) /
     & -0.4989, 0.0558, 0.0679, 0.1065, 0.1174, -0.3548,
     &  0.0620, 0.1824, 0.2022, 0.0615, 0.0670,  0.0972,
     & -0.6139, 0.4186,-0.5246,-0.0330, 0.4152, -0.0128,
     & -0.5067, 0.0695, 0.1617,-0.7926, 0.6977, -0.5766,
     &  0.3624, 0.6344,-0.5746,-0.4566, 0.3423, -0.1416,
     &  0.0943/

      data (ambtnk(i),i=1339,1369) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'N*', 'CB', 'CB',
     & 'NB', 'CK', 'H5', 'NC', 'CA', 'NA', 'H ', 'C ', 'O ',
     & 'N2', 'H ', 'CT', 'H1'/ 

      data (ambvdt(i),i=1339,1369) /
     & 23,1,35,1,35,23,1,36,1,
     & 35,1,35,22,31,23,18,9,9,
     & 16,12,40,17,3,15,29,2,24,
     & 19,29,1,35/

c H5, two protons, CM1 doubles for CM2, HM1 doubles for 6 H

      data (ambstr(i),i=1370,1399) /
     & "M2G         O5'     ", "M2G         C5'     ",
     & "M2G         H5'     ", "M2G         C4'     ", 
     & "M2G         H4'     ", "M2G         O4'     ", 
     & "M2G         C1'     ", "M2G         H1'     ", 
     & "M2G         C3'     ", "M2G         H3'     ", 
     & "M2G         C2'     ", "M2G         H2'1    ", 
     & "M2G         O2'     ", "M2G         HO'2    ", 
     & "M2G         O3'     ", "M2G         N9      ", 
     & "M2G         C4      ", "M2G         C5      ", 
     & "M2G         N7      ", "M2G         C8      ", 
     & "M2G         N3      ", "M2G         C2      ", 
     & "M2G         N1      ", "M2G         C6      ", 
     & "M2G         O6      ", "M2G         HN1     ", 
     & "M2G         H8      ", "M2G         N2      ",
     & "M2G         CM1,2   ", "M2G         HM1     "/ 

      data (ambchg(i),i=1370,1399) /
     & -0.4989, 0.0558, 0.0679, 0.1065, 0.1174, -0.3548,
     &  0.0454, 0.1824, 0.2022, 0.0615, 0.0670,  0.0972,
     & -0.6139, 0.4186,-0.5246,-0.0152, 0.2683,  0.0452,
     & -0.5432, 0.1262,-0.6273, 0.5994,-0.6435,  0.6715,
     & -0.5861, 0.3676, 0.1463,-0.1106,-0.2214,  0.1054/

      data (ambtnk(i),i=1370,1399) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'N*', 'CB', 'CB',
     & 'NB', 'CK', 'NC', 'CA', 'NA', 'C ', 'O ', 'H ', 'H5',
     & 'N2', 'CT', 'H1'/

      data (ambvdt(i),i=1370,1399) /
     & 23,1,35,1,35,23,1,36,1,
     & 35,1,35,22,31,23,18,9,9,
     & 16,12,17,3,15,2,24,29,40,
     & 19,1,35/


c must change molden rdpdb addhd, 7mg, to only add 1 H8 and 3 HM7

      data (ambstr(i),i=1400,1430) /
     & "7MG         O5'     ", "7MG         C5'     ",
     & "7MG         H5'     ", "7MG         C4'     ", 
     & "7MG         H4'     ", "7MG         O4'     ", 
     & "7MG         C1'     ", "7MG         H1'     ", 
     & "7MG         C3'     ", "7MG         H3'     ", 
     & "7MG         C2'     ", "7MG         H2'1    ", 
     & "7MG         O2'     ", "7MG         HO'2    ", 
     & "7MG         O3'     ", "7MG         N9      ", 
     & "7MG         C4      ", "7MG         C5      ", 
     & "7MG         N7      ", "7MG         C8      ", 
     & "7MG         N3      ", "7MG         C2      ", 
     & "7MG         N1      ", "7MG         C6      ", 
     & "7MG         O6      ", "7MG         HN1     ", 
     & "7MG         H8      ", "7MG         N2      ", 
     & "7MG         HN2     ", "7MG         CM7     ",
     & "7MG         HM7     "/

      data (ambchg(i),i=1400,1430) /
     & -0.4989, 0.0558, 0.0679, 0.1065, 0.1174, -0.3548,
     &  0.1274, 0.1824, 0.2022, 0.0615, 0.0670,  0.0972,
     & -0.6139, 0.4186,-0.5246, 0.0532, 0.5114, -0.4118,
     &  0.1353,-0.0615,-0.8029, 1.1308,-0.8380,  0.8701,
     & -0.5650, 0.4351, 0.2571,-1.1211, 0.4985, -0.2704,
     &  0.1622/

      data (ambtnk(i),i=1400,1430) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'N*', 'CB', 'CB',
     & 'N*', 'CK', 'NC', 'CA', 'NA', 'C ', 'O ', 'H ', 'H5',
     & 'N2', 'H ', 'CT', 'H1'/

      data (ambvdt(i),i=1400,1430) /
     & 23,1,35,1,35,23,1,36,1,
     & 35,1,35,22,31,23,18,9,9,
     & 18,12,17,3,15,2,24,29,40,
     & 19,29,1,35/

      data (ambstr(i),i=1431,1460) /
     & "OMG         O5'     ", "OMG         C5'     ",
     & "OMG         H5'     ", "OMG         C4'     ", 
     & "OMG         H4'     ", "OMG         O4'     ", 
     & "OMG         C1'     ", "OMG         H1'     ", 
     & "OMG         C3'     ", "OMG         H3'     ", 
     & "OMG         C2'     ", "OMG         H2'1    ", 
     & "OMG         O2'     ", "OMG         O3'     ", 
     & "OMG         N9      ", "OMG         C4      ", 
     & "OMG         C5      ", "OMG         N7      ", 
     & "OMG         C8      ", "OMG         N3      ", 
     & "OMG         C2      ", "OMG         N1      ", 
     & "OMG         C6      ", "OMG         O6      ", 
     & "OMG         HN1     ", "OMG         H8      ", 
     & "OMG         N2      ", "OMG         HN2     ", 
     & "OMG         CM2     ", "OMG         HM2     "/

      data (ambchg(i),i=1431,1460) /
     & -0.4989, 0.0558, 0.0679, 0.1065, 0.1174, -0.3548,
     &  0.0151, 0.2006, 0.2022, 0.0615, 0.0847,  0.0490,
     & -0.3668,-0.5246, 0.0492, 0.1222, 0.1744, -0.5709,
     &  0.1374,-0.6323, 0.7657,-0.4787, 0.4770, -0.5597,
     &  0.3424, 0.1640,-0.9672, 0.4364,-0.0823,  0.0961/

      data (ambtnk(i),i=1431,1460) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OS', 'OS', 'N*', 'CB', 'CB', 'NB',
     & 'CK', 'NC', 'CA', 'NA', 'C ', 'O ', 'H ', 'H5', 'N2',
     & 'H ', 'CT', 'H1'/ 

      data (ambvdt(i),i=1431,1460) /
     & 23,1,35,1,35,23,1,36,1,
     & 35,1,35,23,23,18,9,9,16,
     & 12,17,3,15,2,24,29,40,19,
     & 29,1,35/

      data (ambstr(i),i=1461,1511) /
     & "YG          O5'     ", "YG          C5'     ",
     & "YG          H5'     ", "YG          C4'     ", 
     & "YG          H4'     ", "YG          O4'     ", 
     & "YG          C1'     ", "YG          H1'     ", 
     & "YG          C3'     ", "YG          H3'     ", 
     & "YG          C2'     ", "YG          H2'1    ", 
     & "YG          O2'     ", "YG          HO'2    ", 
     & "YG          O3'     ", "YG          N9      ", 
     & "YG          C4      ", "YG          C5      ", 
     & "YG          N7      ", "YG          C8      ", 
     & "YG          N3      ", "YG          C2      ", 
     & "YG          N1      ", "YG          C6      ", 
     & "YG          O6      ", "YG          N2      ", 
     & "YG          C11     ", "YG          C12     ", 
     & "YG          H8      ", "YG          C3      ",
     & "YG          H3      ", "YG          C10     ", 
     & "YG          H10     ", "YG          C13     ",
     & "YG          H13     ", "YG          C14     ",
     & "YG          H14     ", "YG          C15     ",
     & "YG          H15     ", "YG          C16     ",
     & "YG          O17     ", "YG          O18     ",
     & "YG          C19     ", "YG          H19     ",
     & "YG          N20     ", "YG          HN2     ",
     & "YG          C21     ", "YG          O22     ",
     & "YG          O23     ", "YG          C24     ",
     & "YG          H24     "/

      data (ambchg(i),i=1461,1511) /
     & -0.4989, 0.0558, 0.0679, 0.1065, 0.1174, -0.3548,
     &  0.0079, 0.1624, 0.2022, 0.0615, 0.0670,  0.0972,
     & -0.6139, 0.4186,-0.5246, 0.1012,-0.0231,  0.0387,
     & -0.5390, 0.1773,-0.0117, 0.2754,-0.0815,  0.5824,
     & -0.4941,-0.5478, 0.4200,-0.2749, 0.1373, -0.0739,
     &  0.0786,-0.4901, 0.1444,-0.0089, 0.0648,  0.0045,
     &  0.0515,-0.4169, 0.1872, 0.9330,-0.6264, -0.3510,
     & -0.1991, 0.1286,-0.3326, 0.2626, 0.8298, -0.6305,
     & -0.4476, 0.1348, 0.0396/

      data (ambtnk(i),i=1461,1511) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H1', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'N*', 'CB', 'CB',
     & 'NB', 'CK', 'N*', 'CB', 'N*', 'C ', 'O ', 'NB', 'CC',
     & 'CC', 'H5', 'CT', 'H1', 'CT', 'HC', 'CT', 'HC', 'CT',
     & 'HC', 'CT', 'H1', 'C ', 'O ', 'OS', 'CT', 'H1', 'N ',
     & 'H ', 'C ', 'O ', 'OS', 'CT', 'H1'/ 

      data (ambvdt(i),i=1461,1511) /
     & 23,1,35,1,35,23,1,35,1,
     & 35,1,35,22,31,23,18,9,9,
     & 16,12,18,9,18,2,24,16,5,
     & 5,40,1,35,1,34,1,34,1,
     & 34,1,35,2,24,23,1,35,14,
     & 29,2,24,23,1,35/

      data (ambstr(i),i=1512,1537) /
     & "H2U      O5'        ", "H2U      C5'        ",
     & "H2U      H5'        ", "H2U      C4'        ", 
     & "H2U      H4'        ", "H2U      O4'        ", 
     & "H2U      C1'        ", "H2U      H1'        ", 
     & "H2U      C3'        ", "H2U      H3'        ", 
     & "H2U      C2'        ", "H2U      H2'1       ", 
     & "H2U      O2'        ", "H2U      HO'2       ", 
     & "H2U      O3'        ", "H2U      N1         ", 
     & "H2U      C2         ", "H2U      N3         ", 
     & "H2U      C4         ", "H2U      C5         ", 
     & "H2U      C6         ", "H2U      O2         ", 
     & "H2U      O4         ", "H2U      HN3        ", 
     & "H2U      H5         ", "H2U      H6         "/

      data (ambchg(i),i=1512,1537) /
     & -0.4989, 0.0558, 0.0679, 0.1065, 0.1174, -0.3548,
     &  0.0445, 0.1824, 0.2022, 0.0615, 0.0670,  0.0972,
     & -0.6139, 0.4186,-0.5246,-0.2939, 0.8257, -0.7280,
     &  0.7942,-0.3335, 0.0910,-0.6148,-0.5941,  0.3861,
     &  0.1202, 0.0580/

      data (ambtnk(i),i=1512,1537) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'N*', 'C ', 'N ',
     & 'C ', 'CT', 'CT', 'O ', 'O ', 'H ', 'HC', 'H1'/ 

      data (ambvdt(i),i=1512,1537) /
     & 23,1,35,1,35,23,1,36,1,
     & 35,1,35,22,31,23,18,2,14,
     & 2,1,1,24,24,29,34,35/


      data (ambstr(i),i=1538,1564) /
     & "5MU      O5'        ", "5MU      C5'        ",
     & "5MU      H5'        ", "5MU      C4'        ", 
     & "5MU      H4'        ", "5MU      O4'        ", 
     & "5MU      C1'        ", "5MU      H1'        ", 
     & "5MU      C3'        ", "5MU      H3'        ", 
     & "5MU      C2'        ", "5MU      H2'1       ", 
     & "5MU      O2'        ", "5MU      HO'2       ", 
     & "5MU      O3'        ", "5MU      N1         ", 
     & "5MU      C2         ", "5MU      N3         ", 
     & "5MU      C4         ", "5MU      C5         ", 
     & "5MU      C6         ", "5MU      O2         ", 
     & "5MU      O4         ", "5MU      HN3        ", 
     & "5MU      H6         ", "5MU      CM5        ",
     & "5MU      HM5        "/

      data (ambchg(i),i=1538,1564) /
     & -0.4989, 0.0558, 0.0679, 0.1065, 0.1174, -0.3548,
     &  0.0370, 0.1824, 0.2022, 0.0615, 0.0670,  0.0972,
     & -0.6139, 0.4186,-0.5246,-0.0561, 0.7537, -0.7343,
     &  0.7714,-0.0700,-0.1825,-0.6073,-0.5888,  0.4033,
     &  0.2167,-0.3875, 0.1260/

      data (ambtnk(i),i=1538,1564) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H2', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'N*', 'C ', 'NA',
     & 'C ', 'CM', 'CM', 'O ', 'O ', 'H ', 'H4', 'CT', 'HC'/

      data (ambvdt(i),i=1538,1564) /
     & 23,1,35,1,35,23,1,36,1,
     & 35,1,35,22,31,23,18,2,15,
     & 2,4,4,24,24,29,39,1,34/

      data (ambstr(i),i=1565,1590) /
     & "PSU      O5'        ", "PSU      C5'        ",
     & "PSU      H5'        ", "PSU      C4'        ", 
     & "PSU      H4'        ", "PSU      O4'        ", 
     & "PSU      C1'        ", "PSU      H1'        ", 
     & "PSU      C3'        ", "PSU      H3'        ", 
     & "PSU      C2'        ", "PSU      H2'1       ", 
     & "PSU      O2'        ", "PSU      HO'2       ", 
     & "PSU      O3'        ", "PSU      C5         ", 
     & "PSU      C6         ", "PSU      N1         ", 
     & "PSU      C2         ", "PSU      N3         ", 
     & "PSU      C4         ", "PSU      O2         ", 
     & "PSU      O4         ", "PSU      H6         ",
     & "PSU      HN1        ", "PSU      HN3        "/

      data (ambchg(i),i=1565,1590) /
     & -0.4989, 0.0558, 0.0679, 0.1065, 0.1174, -0.3548,
     & -0.1081, 0.2405, 0.2022, 0.0615, 0.0670,  0.0972,
     & -0.6139, 0.4186,-0.5246,-0.5317, 0.3896, -0.6770,
     &  0.9290,-0.7255, 0.8922,-0.6096,-0.6230,  0.1226,
     &  0.4311, 0.3859/

      data (ambtnk(i),i=1565,1590) /
     & 'OS', 'CT', 'H1', 'CT', 'H1', 'OS', 'CT', 'H1', 'CT',
     & 'H1', 'CT', 'H1', 'OH', 'HO', 'OS', 'CM', 'CA', 'NA',
     & 'C ', 'NA', 'C ', 'O ', 'O ', 'H4', 'H ', 'H '/

      data (ambvdt(i),i=1565,1590) /
     & 23,1,35,1,35,23,1,35,1,
     & 35,1,35,22,31,23,4,3,15,
     & 2,15,2,24,24,39,29,29/

CNF  Nothing between 1254 and 2000
CNF  Finally some parameters for my personal use
CNF      data (ambstr(i),i=2012,2015)
CNF     & 'C sp2 Ret           ','C sp3 Ret           ',
CNF     & 'H Ret               ','Cl- Chloride Ion    '/
     
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

      data gffext /
     &   'placeholder wildcard        ',
     &   'C_Sp2 carbonyl group        ',
     &   'C_Sp                        ',
     &   'C_Sp2                       ',
     &   'C_Sp3                       ',
     &   'pure aromatic Sp2 C         ',
     &   'C_Sp2 conn. rings biphenyl A',
     &   'C_Sp2 conn. rings biphenyl B',
     &   'non-pure aromatic C_sp2    A',
     &   'non-pure aromatic C_sp2    B',
     &   'Inner conjugated C_Sp2     A',
     &   'Inner conjugated C_Sp2     B',
     &   'Inner conjugated C_Sp      A',
     &   'Inner conjugated C_Sp      B',
     &   'C_Sp3 in 3-membered ring    ',
     &   'C_Sp3 in 4-membered ring    ',
     &   'C_Sp2 in 3-membered ring    ',
     &   'C_Sp2 in 4-membered ring    ',
     &   'H -C_aliphatic-  elecwd.gr. ',
     &   'H -C_aliphatic- (el.wd.gr.)2',
     &   'H -C_aliphatic- (el.wd.gr.)3',
     &   'H -C_not_Sp3  -  elecwd.gr. ',
     &   'H -C_not_Sp3  - (el.wd.gr.)2',
     &   'H -C_ar                     ',
     &   'H -C_aliphatic- (el.wd.gr.)0',
     &   'H -N                        ',
     &   'H -O                        ',
     &   'H -P                        ',
     &   'H -S                        ',
     &   'H -O-H (water)              ',
     &   'H -C-(X+)                   ',
     &   'Fluorine                    ',
     &   'Chlorine                    ',
     &   'Bromine                     ',
     &   'Iodine                      ',
     &   'N_Sp2 amide group           ',
     &   'N_Sp                        ',
     &   'N_Sp2 aliphatic- (X)2       ',
     &   'N_Sp3 - (X)3                ',
     &   'N_Sp3 - (X)4                ',
     &   'N_Sp2 - (X)3                ',
     &   'pure aromatic N_Sp2         ',
     &   'non-pure aromatic N_Sp2    A',
     &   'non-pure aromatic N_Sp2    B',
     &   'Inner conjugated N_Sp2     A',
     &   'Inner conjugated N_Sp2     B',
     &   'N_Sp3 amine-arom.ring(s)    ',
     &   'Nitro N                     ',
     &   'O_Sp2                       ',
     &   'O_Sp3 OH                    ',
     &   'O_Sp3 Ether and ester       ',
     &   'O_Sp3 in water              ',
     &   'P-(X)2, like P=C            ',
     &   'P_Sp3-(X)3, such as PH3     ',
     &   'P_Sp2-(X)3, like O=P(CH3)2  ',
     &   'P_Sp2-(X)4, like O=P(OH)3   ',
     &   'pure aromatic P_Sp2         ',
     &   'non-pure aromatic P_Sp2    A',
     &   'non-pure aromatic P_Sp2    B',
     &   'Inner conjugated P_Sp2     A',
     &   'Inner conjugated P_Sp2     B',
     &   'P_Sp2-(X)3, conjugated     A',
     &   'P_Sp2-(X)3, conjugated     B',
     &   'S-X                         ',
     &   'S, R1-S=R2, R1=S=R2         ',
     &   'S-(X)3, 3 conn. atoms,      ',
     &   'S-(X)4, 4 conn. atoms.      ',
     &   'S_Sp3, -S-H                 ',
     &   'S_Sp3, thio-ester and -ether',
     &   'S-(X)3, conjugated          ',
     &   'S-(X)4, conjugated          ',
     &   'C Sp2 carbon in guadinium   '/

C Tinker AMOEBA PRO
      data (amotnk(i),i=1,100) / 
     & 'N ','CA','C ','HN','O ','H ','N ','CA','C ','HN',
     & 'O ','H ','C ','H ','C ','H ','C ','H ','C ','H ',
     & 'C ','H ','C ','H ','C ','H ','C ','H ','C ','H ',
     & 'C ','H ','CA','C ','H ','OH','HO','C ','H ','C ',
     & 'H ','OH','HO','CA','C ','H ','SH','HS','SS','N ',
     & 'CA','C ','O ','H ','C ','H ','C ','H ','C ','H ',
     & 'C ','H ','C ','C ','H ','C ','H ','C ','H ','C ',
     & 'H ','C ','C ','H ','C ','H ','C ','OH','HO','C ',
     & 'H ','C ','C ','H ','C ','N ','HN','C ','C ','H ',
     & 'C ','H ','C ','H ','C ','H ','C ','H ','C ','N '/

      data (amotnk(i),i=101,201) / 
     & 'HN','C ','H ','C ','H ','N ','HN','C ','H ','C ',
     & 'N ','HN','C ','H ','C ','H ','N ','C ','H ','C ',
     & 'N ','C ','H ','C ','H ','N ','HN','C ','H ','C ',
     & 'O ','C ','H ','C ','O ','OH','HO','C ','H ','C ',
     & 'O ','N ','HN','C ','H ','C ','H ','C ','O ','C ',
     & 'H ','C ','H ','C ','H ','S ','C ','H ','C ','H ',
     & 'C ','H ','C ','H ','C ','H ','N ','HN','C ','H ',
     & 'C ','H ','C ','H ','N ','HN','C ','N ','HN','C ',
     & 'H ','C ','O ','N ','HN','C ','H ','N ','HN','N ',
     & 'H ','C ','O ','N ','HN','CA','C ','O ','H ','C ',
     & 'H '/

      data (amostr(i),i=1,6) /
     & 'Glycine N           ', 'Glycine CA          ',
     & 'Glycine C           ', 'Glycine HN          ',
     & 'Glycine O           ', 'Glycine HA          '/
      data (amostr(i),i=7,14) /
     & 'Alanine N           ', 'Alanine CA          ',
     & 'Alanine C           ', 'Alanine HN          ',
     & 'Alanine O           ', 'Alanine HA          ',
     & 'Alanine CB          ', 'Alanine HB          '/
      data (amostr(i),i=15,18) /
     & 'Valine CB           ', 'Valine HB           ',
     & 'Valine CG           ', 'Valine HG           '/
      data (amostr(i),i=19,24) /
     & 'Leucine CB          ', 'Leucine HB          ',
     & 'Leucine CG          ', 'Leucine HG          ',
     & 'Leucine CD          ', 'Leucine HD          '/
      data (amostr(i),i=25,32) /
     & 'Isoleucine CB       ', 'Isoleucine HB       ',
     & 'Isoleucine CG       ', 'Isoleucine HG       ',
     & 'Isoleucine CG       ', 'Isoleucine HG       ',
     & 'Isoleucine CD       ', 'Isoleucine HD       '/
      data (amostr(i),i=33,37) /
     & 'Serine CA           ', 'Serine CB           ',
     & 'Serine HB           ', 'Serine OG           ',
     & 'Serine HG           '/
      data (amostr(i),i=38,43) /
     & 'Threonine CB        ', 'Threonine HB        ',
     & 'Threonine CG        ', 'Threonine HG        ',
     & 'Threonine O         ', 'Threonine H         '/
      data (amostr(i),i=44,49) /
     & 'Cysteine CA         ', 'Cysteine CB         ',
     & 'Cysteine HB         ', 'Cysteine SG         ',
     & 'Cysteine HG         ', 'Cystine SG          '/
      data (amostr(i),i=50,60) /
     & 'Proline N           ', 'Proline CA          ',
     & 'Proline C           ', 'Proline O           ',
     & 'Proline HA          ', 'Proline CB          ',
     & 'Proline HB          ', 'Proline CG          ',
     & 'Proline HG          ', 'Proline CD          ',
     & 'Proline HD          '/
      data (amostr(i),i=61,69) /
     & 'Phenylalanine CB    ', 'Phenylalanine HB    ',
     & 'Phenylalanine CG    ', 'Phenylalanine CD    ',
     & 'Phenylalanine HD    ', 'Phenylalanine CE    ',
     & 'Phenylalanine HE    ', 'Phenylalanine CZ    ',
     & 'Phenylalanine HZ    '/
      data (amostr(i),i=70,79) /
     & 'Tyrosine CB         ', 'Tyrosine HB         ',
     & 'Tyrosine CG         ', 'Tyrosine CD         ',
     & 'Tyrosine HD         ', 'Tyrosine CE         ',
     & 'Tyrosine HE         ', 'Tyrosine CZ         ',
     & 'Tyrosine OH         ', 'Tyrosine HH         '/
      data (amostr(i),i=80,96) /
     & 'Tryptophan CB       ', 'Tryptophan HB       ',
     & 'Tryptophan CG       ', 'Tryptophan CD1      ',
     & 'Tryptophan HD1      ', 'Tryptophan CD2      ',
     & 'Tryptophan NE1      ', 'Tryptophan HE1      ',
     & 'Tryptophan CE2      ', 'Tryptophan CE3      ',
     & 'Tryptophan HE3      ', 'Tryptophan CZ2      ',
     & 'Tryptophan HZ2      ', 'Tryptophan CZ3      ',
     & 'Tryptophan HZ3      ', 'Tryptophan CH2      ',
     & 'Tryptophan HH2      '/
      data (amostr(i),i=97,107) /
     & 'Histidine(+) CB     ', 'Histidine(+) HB     ',
     & 'Histidine(+) CG     ', 'Histidine(+) ND1    ',
     & 'Histidine(+) HD1    ', 'Histidine(+) CD2    ',
     & 'Histidine(+) HD2    ', 'Histidine(+) CE1    ',
     & 'Histidine(+) HE1    ', 'Histidine(+) NE2    ',
     & 'Histidine(+) HE2    '/
      data (amostr(i),i=108,117) /
     & 'Histidine/HID CB    ', 'Histidine/HID HB    ',
     & 'Histidine/HID CG    ', 'Histidine/HID ND1   ',
     & 'Histidine/HID HD1   ', 'Histidine/HID CD2   ',
     & 'Histidine/HID HD2   ', 'Histidine/HID CE1   ',
     & 'Histidine/HID HE1   ', 'Histidine/HID NE2   '/
      data (amostr(i),i=118,127) /
     & 'Histidine/HIE CB    ', 'Histidine/HIE HB    ',
     & 'Histidine/HIE CG    ', 'Histidine/HIE ND1   ',
     & 'Histidine/HIE CD2   ', 'Histidine/HIE HD2   ',
     & 'Histidine/HIE CE1   ', 'Histidine/HIE HE1   ',
     & 'Histidine/HIE NE2   ', 'Histidine/HIE HE2   '/
      data (amostr(i),i=128,131) /
     & 'Aspartate CB        ', 'Aspartate HB        ',
     & 'Aspartate CG        ', 'Aspartate OG        '/
      data (amostr(i),i=132,137) /
     & 'Aspartic Acid CB    ', 'Aspartic Acid HB    ',
     & 'Aspartic Acid CG    ', 'Aspartic Acid OG    ',
     & 'Aspartic Acid OH    ', 'Aspartic Acid HO    '/
      data (amostr(i),i=138,143) /
     & 'Asparagine CB       ', 'Asparagine HB       ',
     & 'Asparagine CG       ', 'Asparagine OD1      ',
     & 'Asparagine ND2      ', 'Asparagine HD2      '/
      data (amostr(i),i=144,149) /
     & 'Glutamate CB        ', 'Glutamate HB        ',
     & 'Glutamate CG        ', 'Glutamate HG        ',
     & 'Glutamate CD        ', 'Glutamate OD        '/
      data (amostr(i),i=150,151) /
     & 'Glutamine CB        ', 'Glutamine HB        '/
      data (amostr(i),i=152,158) /
     & 'Methionine CB       ', 'Methionine HB       ',
     & 'Methionine CG       ', 'Methionine HG       ',
     & 'Methionine SD       ', 'Methionine CE       ',
     & 'Methionine HE       '/
      data (amostr(i),i=159,168) /
     & 'Lysine CB           ', 'Lysine HB           ',
     & 'Lysine CG           ', 'Lysine HG           ',
     & 'Lysine CD           ', 'Lysine HD           ',
     & 'Lysine CE           ', 'Lysine HE           ',
     & 'Lysine NZ           ', 'Lysine HN           '/
      data (amostr(i),i=169,179) /
     & 'Arginine CB         ', 'Arginine HB         ',
     & 'Arginine CG         ', 'Arginine HG         ',
     & 'Arginine CD         ', 'Arginine HD         ',
     & 'Arginine NE         ', 'Arginine HE         ',
     & 'Arginine CZ         ', 'Arginine NH         ',
     & 'Arginine HH         '/
      data (amostr(i),i=180,183) /
     & 'Acetyl CH3          ', 'Acetyl H3C          ',
     & 'Acetyl C            ', 'Acetyl O            '/
      data (amostr(i),i=184,187) /
     & 'N-MeAmide N         ', 'N-MeAmide HN        ',
     & 'N-MeAmide CH3       ', 'N-MeAmide H3C       '/
      data (amostr(i),i=188,189) /
     & 'Amide NH2           ', 'Amide H2N           '/
      data (amostr(i),i=190,191) /
     & 'N-Terminal NH3+     ', 'N-Terminal H3N+     '/
      data (amostr(i),i=192,193) /
     & 'C-Terminal COO-     ', 'C-Terminal COO-     '/
      data (amostr(i),i=194,201) /
     & 'N-Terminal PRO NH2+ ', 'N-Terminal PRO H2N+ ',
     & 'N-Terminal PRO CA   ', 'N-Terminal PRO C    ',
     & 'N-Terminal PRO O    ', 'N-Terminal PRO HA   ',
     & 'N-Terminal PRO CD   ', 'N-Terminal PRO HD   '/

      data ncco /0,1,2,2,3,4,3,4,4,4,4,5,5,5,3,7,6,7,8,10,2,6,6,4/
      data nhho /0,3,3,3,6,9,7,7,2,4,9,11,6,6,6,11,6,7,7,8,2,5,5,2/

c GLY
      data ((icco(i,j,1),i=1,2),j=1,mxato) /
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,1),i=1,2),j=1,mxatho) /
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c ALA
      data ((icco(i,j,2),i=1,2),j=1,mxato) /
     & 5,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,2),i=1,2),j=1,mxatho) /
     & 7,14,8,14,9,14,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c SER
      data ((icco(i,j,3),i=1,2),j=1,mxato) /
     & 5,34,31,36,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,3),i=1,2),j=1,mxatho) /
     & 10,37,7,14,8,14,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c CYS -SH
      data ((icco(i,j,4),i=1,2),j=1,mxato) /
     & 5,45,37,47,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,4),i=1,2),j=1,mxatho) /
     & 10,48,7,46,8,46,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c THR
      data ((icco(i,j,5),i=1,2),j=1,mxato) /
     & 5,38,8,40,32,42,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,5),i=1,2),j=1,mxatho) /
     & 10,43,13,43,16,41,17,41,18,41,7,39,0,0,0,0,0,0,0,0,0,0/
c ILE
      data ((icco(i,j,6),i=1,2),j=1,mxato) /
     & 5,25,7,27,8,29,10,31,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,6),i=1,2),j=1,mxatho) /
     & 13,28,14,28,16,30,17,30,18,30,22,32,23,32,24,32,7,26,0,0,0,0/
c VAL
      data ((icco(i,j,7),i=1,2),j=1,mxato) /
     & 5,15,7,17,8,17,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,7),i=1,2),j=1,mxatho) /
     & 13,18,14,18,15,18,16,18,17,18,18,18,7,16,0,0,0,0,0,0,0,0/
c MET
      data ((icco(i,j,8),i=1,2),j=1,mxato) /
     & 5,152,6,154,12,157,36,156,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,8),i=1,2),j=1,mxatho) /
     & 10,155,11,155,28,158,29,158,30,158,7,153,8,153,0,0,0,0,0,0,0,0/
c ASP
      data ((icco(i,j,9),i=1,2),j=1,mxato) /
     & 5,128,6,130,29,131,30,131,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,9),i=1,2),j=1,mxatho) /
     & 7,129,8,129,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c ASN
      data ((icco(i,j,10),i=1,2),j=1,mxato) /
     & 5,138,6,140,29,141,21,142,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,10),i=1,2),j=1,mxatho) /
     & 25,143,26,143,7,139,8,139,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c LEU
      data ((icco(i,j,11),i=1,2),j=1,mxato) /
     & 5,19,6,21,10,23,11,23,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,11),i=1,2),j=1,mxatho) /
     & 10,22,22,24,23,24,24,24,25,24,26,24,27,24,7,20,8,20,0,0,0,0/
c LYS
      data ((icco(i,j,12),i=1,2),j=1,mxato) /
     & 5,159,6,161,9,163,12,165,27,167,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,12),i=1,2),j=1,mxatho) /
     & 10,162,11,162,19,164,20,164,28,166,29,166,
     & 40,168,41,168,42,168,7,160,8,160/
c GLU
      data ((icco(i,j,13),i=1,2),j=1,mxato) /
     & 5,144,6,146,9,148,34,149,35,149,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,13),i=1,2),j=1,mxatho) /
     & 10,147,11,147,7,145,8,145,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c GLN
      data ((icco(i,j,14),i=1,2),j=1,mxato) /
     & 5,150,6,138,9,140,34,141,24,142,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,14),i=1,2),j=1,mxatho) /
     & 10,139,11,139,34,143,35,143,7,151,8,151,0,0,0,0,0,0,0,0,0,0/
c PRO
      data ((icco(i,j,15),i=1,2),j=1,mxato) /
     & 5,55,6,57,9,59,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,15),i=1,2),j=1,mxatho) /
     & 10,58,11,58,19,60,20,60,7,56,8,56,0,0,0,0,0,0,0,0,0,0/
c ARG
      data ((icco(i,j,16),i=1,2),j=1,mxato) /
     & 5,169,6,171,9,173,17,177,22,175,25,178,26,178,0,0,0,0,0,0/
      data ((ihho(i,j,16),i=1,2),j=1,mxatho) /
     & 10,172,11,172,19,174,20,174,28,176,
     & 55,179,56,179,58,179,59,179,7,170,8,170/
c HIS (+, HIP)
      data ((icco(i,j,17),i=1,2),j=1,mxato) /
     & 5,97,6,99,11,102,13,104,20,100,24,106,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,17),i=1,2),j=1,mxatho) /
     & 22,101,25,103,31,105,34,107,7,98,8,98,0,0,0,0,0,0,0,0,0,0/
c PHE
      data ((icco(i,j,18),i=1,2),j=1,mxato) /
     & 5,61,6,63,10,64,11,64,13,66,14,66,17,68,0,0,0,0,0,0/
      data ((ihho(i,j,18),i=1,2),j=1,mxatho) /
     & 22,65,25,65,31,67,34,67,40,69,7,62,8,62,0,0,0,0,0,0,0,0/
c TYR
      data ((icco(i,j,19),i=1,2),j=1,mxato) /
     & 5,70,6,72,10,73,11,73,13,75,14,75,17,77,33,78,0,0,0,0/
      data ((ihho(i,j,19),i=1,2),j=1,mxatho) /
     & 22,74,25,74,31,76,34,76,52,79,7,71,8,71,0,0,0,0,0,0,0,0/
c TRP
      data ((icco(i,j,20),i=1,2),j=1,mxato) /
     & 5,80,6,82,10,83,11,85,14,88,15,89,16,95,
     & 18,91,19,93,23,86/
      data ((ihho(i,j,20),i=1,2),j=1,mxatho) /
     & 22,84,31,87,37,90,46,92,49,94,58,96,7,81,8,81,0,0,0,0,0,0/
c CYS -SS-
      data ((icco(i,j,21),i=1,2),j=1,mxato) /
     & 5,45,37,49,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,21),i=1,2),j=1,mxatho) /
     & 7,45,8,45,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c HISD
      data ((icco(i,j,22),i=1,2),j=1,mxato) /
     & 5,108,6,110,11,113,13,115,20,111,24,117,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,22),i=1,2),j=1,mxatho) /
     & 22,112,25,114,31,116,7,109,8,109,0,0,0,0,0,0,0,0,0,0,0,0/
c HISE
      data ((icco(i,j,23),i=1,2),j=1,mxato) /
     & 5,118,6,120,11,122,13,124,20,121,24,126,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,23),i=1,2),j=1,mxatho) /
     & 25,123,31,125,34,127,7,119,8,119,0,0,0,0,0,0,0,0,0,0,0,0/
c ASP-H
      data ((icco(i,j,24),i=1,2),j=1,mxato) /
     & 5,132,6,134,29,135,30,135,0,0,0,0,0,0,0,0,0,0,0,0/
      data ((ihho(i,j,24),i=1,2),j=1,mxatho) /
     & 7,133,8,133,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/



      data nnuc /22,23,20,20,21,22,19,23,23,21,21,24,25,24,24,39,
     &           20,21,20/
      data nhnuc /13,14,11,11,13,14,11,16,17,13,13,17,16,17,16,27,
     &            13,13,11/
      data irna /1,3,2,8,4,9,10,11,2,12,13,14,15,16,2,4,17,18,19/

c r-adenosine

      data ((inuc(i,j,1),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1001,47,1002,48,1005,49,1007,
     & 50,1010,51,1016,52,1012,53,1014,54,1008,55,1023,56,1018,
     & 57,1019,58,1025,59,1021,60,1024,62,1022,64,1027,65,1020,
     & 66,1017,34*0/

      data ((ihnuc(i,j,1),i=1,2),j=1,mxhnuc) /
     & 61,1003,62,1004,64,1006,67,1011,70,1013,73,1009,
     & 76,1015,82,1026,112,1030,136,1028,137,1029,
     & 103,1028,104,1029,28*0/

c r-guanosine

      data ((inuc(i,j,2),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1031,47,1032,48,1035,49,1037,
     & 50,1040,51,1046,52,1042,53,1044,54,1038,55,1053,56,1048,
     & 57,1049,58,1055,59,1051,60,1054,61,1057,62,1052,65,1050,
     & 66,1047,69,1060,32*0/

      data ((ihnuc(i,j,2),i=1,2),j=1,mxhnuc) /
     & 61,1033,62,1034,64,1036,67,1041,70,1043,73,1039,76,1045,
     & 112,1061,121,1056,124,1058,125,1059,
     & 79,1056,82,1058,83,1059,26*0/

c r-cytosine

      data ((inuc(i,j,3),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1062,47,1063,48,1066,49,1068,
     & 50,1071,51,1077,52,1073,53,1075,54,1069,55,1079,56,1081,
     & 57,1082,58,1083,60,1078,62,1080,63,1085,67,1084,38*0/

      data ((ihnuc(i,j,3),i=1,2),j=1,mxhnuc) /
     & 61,1064,62,1065,64,1067,67,1072,70,1074,73,1070,76,1076,
     & 100,1088,103,1089,130,1086,131,1087,32*0/

c r-uracil

      data ((inuc(i,j,4),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1090,47,1091,48,1094,49,1096,
     & 50,1099,51,1105,52,1101,53,1103,54,1097,55,1107,56,1109,
     & 57,1110,58,1111,60,1106,62,1108,67,1112,68,1114,38*0/
 
      data ((ihnuc(i,j,4),i=1,2),j=1,mxhnuc) /
     & 61,1092,62,1093,64,1095,67,1100,70,1102,73,1098,76,1104,
     & 100,1115,103,1116,127,1113,91,1113,32*0/

c d-adenosine

      data ((inuc(i,j,5),i=1,2),j=1,mxnucl) /
     & 43,1242,44,1243,45,1243,46,1117,47,1118,48,1121,49,1123,
     & 50,1126,51,1131,52,1128,54,1124,55,1138,56,1133,57,1134,
     & 58,1140,59,1136,60,1139,62,1137,64,1142,65,1135,66,1132,
     & 36*0/

      data ((ihnuc(i,j,5),i=1,2),j=1,mxhnuc) /
     & 61,1119,62,1120,64,1122,67,1127,70,1129,71,1130,73,1125,
     & 82,1141,112,1145,136,1143,137,1144,103,1143,104,1144,28*0/

c d-guanosine

      data ((inuc(i,j,6),i=1,2),j=1,mxnucl) /
     & 43,1242,44,1243,45,1243,46,1146,47,1147,48,1150,49,1152,
     & 50,1155,51,1160,52,1157,54,1153,55,1167,56,1162,57,1163,
     & 58,1169,59,1165,60,1168,61,1171,62,1166,65,1164,66,1161,
     & 69,1174,34*0/

      data ((ihnuc(i,j,6),i=1,2),j=1,mxhnuc) /
     & 61,1148,62,1149,64,1151,67,1156,70,1158,71,1159,73,1154,
     & 112,1175,121,1170,124,1172,125,1173,
     & 79,1170,82,1172,83,1173,26*0/

c d-cytosine

      data ((inuc(i,j,7),i=1,2),j=1,mxnucl) /
     & 43,1242,44,1243,45,1243,46,1176,47,1177,48,1180,49,1182,
     & 50,1185,51,1190,52,1187,54,1183,55,1192,56,1194,57,1195,
     & 58,1196,60,1191,62,1193,63,1198,67,1197,40*0/

      data ((ihnuc(i,j,7),i=1,2),j=1,mxhnuc) /
     & 61,1178,62,1179,64,1181,67,1186,70,1188,71,1189,73,1184,
     & 100,1201,103,1202,130,1199,131,1200,32*0/

c d-thymine

      data ((inuc(i,j,8),i=1,2),j=1,mxnucl) /
     & 43,1242,44,1243,45,1243,46,1203,47,1204,48,1207,49,1209,
     & 50,1212,51,1217,52,1214,54,1210,55,1219,56,1221,57,1222,
     & 58,1223,70,1227,75,1227,83,1227,85,1227,60,1218,62,1220,
     & 67,1224,68,1226,32*0/

      data ((ihnuc(i,j,8),i=1,2),j=1,mxhnuc) /
     & 61,1205,62,1206,64,1208,67,1213,70,1215,71,1216,73,1211,
     & 160,1228,161,1228,162,1228,103,1229,127,1225,91,1225,
     & 166,1228,167,1228,168,1228,22*0/

c 1MA

      data ((inuc(i,j,9),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1254,47,1255,48,1258,49,1260,
     & 50,1263,51,1269,52,1265,53,1267,54,1261,55,1276,56,1271,
     & 57,1272,58,1278,59,1274,60,1277,62,1275,64,1280,65,1273,
     & 66,1270,79,1283,32*0/

      data ((ihnuc(i,j,9),i=1,2),j=1,mxhnuc) /
     & 61,1256,62,1257,64,1259,67,1264,70,1266,73,1262,
     & 76,1268,82,1279,112,1282,136,1281,137,1281,
     & 148,1284,149,1284,150,1284,103,1281,104,1281,22*0/


c 5MC

      data ((inuc(i,j,10),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1285,47,1286,48,1288,49,1290,
     & 50,1293,51,1299,52,1295,53,1297,54,1291,55,1301,56,1303,
     & 57,1304,58,1305,60,1300,62,1302,63,1308,67,1307,83,1310,
     & 36*0/

      data ((ihnuc(i,j,10),i=1,2),j=1,mxhnuc) /
     & 61,1287,62,1287,64,1289,67,1294,70,1296,73,1292,76,1298,
     & 103,1306,130,1309,131,1309,160,1311,161,1311,162,1311,28*0/

c OMC

      data ((inuc(i,j,11),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1312,47,1313,48,1315,49,1317,
     & 50,1320,51,1325,52,1322,53,1324,54,1318,55,1327,56,1330,
     & 57,1333,58,1335,60,1326,62,1329,63,1331,67,1328,80,1337,
     & 36*0/

      data ((ihnuc(i,j,11),i=1,2),j=1,mxhnuc) /
     & 61,1314,62,1314,64,1316,67,1321,70,1323,73,1319,103,1336,
     & 100,1334,130,1332,131,1332,151,1338,152,1338,153,1338,28*0/

c 2MG

      data ((inuc(i,j,12),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1339,47,1340,48,1342,49,1344,
     & 50,1347,51,1353,52,1349,53,1351,54,1345,55,1361,56,1355,
     & 57,1356,58,1364,59,1358,60,1362,61,1366,62,1360,65,1357,
     & 66,1354,69,1365,80,1368,30*0/

      data ((ihnuc(i,j,12),i=1,2),j=1,mxhnuc) /
     & 61,1341,62,1341,64,1343,67,1348,70,1350,73,1346,76,1352,
     & 112,1359,121,1363,124,1367,125,1367,
     & 151,1369,152,1369,153,1369,79,1363,82,1367,83,1367,20*0/

c M2G

      data ((inuc(i,j,13),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1370,47,1371,48,1373,49,1375,
     & 50,1378,51,1384,52,1380,53,1382,54,1376,55,1391,56,1386,
     & 57,1387,58,1393,59,1389,60,1392,61,1397,62,1390,65,1388,
     & 66,1385,69,1394,79,1398,80,1398,28*0/


      data ((ihnuc(i,j,13),i=1,2),j=1,mxhnuc) /
     & 61,1372,62,1372,64,1374,67,1379,70,1381,73,1377,76,1383,
     & 112,1396,121,1395,148,1399,149,1399,
     & 150,1399,151,1399,152,1399,153,1399,
     & 79,1395,22*0/

c 7MG

      data ((inuc(i,j,14),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1400,47,1401,48,1403,49,1405,
     & 50,1408,51,1414,52,1410,53,1412,54,1406,55,1421,56,1416,
     & 57,1417,58,1423,59,1419,60,1422,61,1427,62,1420,65,1418,
     & 66,1415,69,1424,85,1429,30*0/

      data ((ihnuc(i,j,14),i=1,2),j=1,mxhnuc) /
     & 61,1402,62,1402,64,1404,67,1409,70,1411,73,1407,76,1413,
     & 112,1426,121,1425,124,1428,125,1428,
     & 166,1430,167,1430,168,1430,
     & 79,1425,82,1428,83,1428,20*0/

c OMG

      data ((inuc(i,j,15),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1431,47,1432,48,1434,49,1436,
     & 50,1439,51,1444,52,1441,53,1443,54,1437,55,1451,56,1446,
     & 57,1447,58,1453,59,1449,60,1452,61,1457,62,1450,65,1448,
     & 66,1445,69,1454,80,1459,30*0/

      data ((ihnuc(i,j,15),i=1,2),j=1,mxhnuc) /
     & 61,1433,62,1433,64,1435,67,1440,70,1442,73,1438,
     & 112,1456,121,1455,124,1458,125,1458,
     & 151,1460,152,1460,153,1460,
     & 79,1455,82,1458,83,1458,22*0/

c YG

      data ((inuc(i,j,16),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1461,47,1462,48,1464,49,1466,
     & 50,1469,51,1475,52,1471,53,1473,54,1467,55,1482,56,1477,
     & 57,1478,58,1484,59,1480,60,1483,61,1486,62,1481,65,1479,
     & 66,1476,69,1485,88,1490,89,1492,90,1487,91,1488,92,1494,
     & 93,1496,94,1498,95,1500,96,1501,97,1502,98,1503,99,1505,
     & 100,1507,101,1508,102,1509,103,1510/

      data ((ihnuc(i,j,16),i=1,2),j=1,mxhnuc) /
     & 61,1463,62,1463,64,1465,67,1470,70,1472,73,1468,76,1474,
     & 112,1489,91,1491,92,1491,93,1491,124,1506,175,1493,
     & 176,1493,177,1493,178,1495,179,1495,181,1497,182,1497,
     & 184,1499,187,1504,188,1504,189,1504,190,1511,191,1511,
     & 192,1511,82,1506/

c H2U

      data ((inuc(i,j,17),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1512,47,1513,48,1515,49,1517,
     & 50,1520,51,1526,52,1522,53,1524,54,1518,55,1528,56,1530,
     & 57,1531,58,1532,60,1527,62,1529,67,1533,68,1534,38*0/
 
      data ((ihnuc(i,j,17),i=1,2),j=1,mxhnuc) /
     & 61,1514,62,1514,64,1516,67,1521,70,1523,73,1519,76,1525,
     & 100,1536,101,1536,103,1537,104,1537,127,1535,91,1537,
     & 28*0/

c 5MU

      data ((inuc(i,j,18),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1538,47,1539,48,1541,49,1543,
     & 50,1546,51,1552,52,1548,53,1550,54,1544,55,1554,56,1556,
     & 57,1557,58,1558,60,1553,62,1555,67,1559,68,1560,83,1563,
     & 36*0/
 
      data ((ihnuc(i,j,18),i=1,2),j=1,mxhnuc) /
     & 61,1540,62,1540,64,1542,67,1547,70,1549,73,1545,76,1551,
     & 103,1562,127,1561,91,1561,160,1564,161,1564,162,1564,
     & 28*0/

c PSU

      data ((inuc(i,j,19),i=1,2),j=1,mxnucl) /
     & 43,1230,44,1231,45,1231,46,1565,47,1566,48,1568,49,1570,
     & 50,1573,51,1579,52,1575,53,1577,54,1571,55,1583,56,1585,
     & 57,1580,58,1581,60,1582,62,1584,67,1586,68,1587,38*0/
 
      data ((ihnuc(i,j,19),i=1,2),j=1,mxhnuc) /
     & 61,1567,62,1567,64,1569,67,1574,70,1576,73,1572,76,1578,
     & 103,1588,121,1589,127,1590,91,1590,32*0/

      data (chmsf(i),i=1,160) /
     & 'H   ', 'HC  ', 'HA  ', 'HT  ', 'LP  ', 'MBE ', 'B   ', 'HO  ',
     & 'HMU ', 'CT  ', 'CH1E', 'CH2E', 'CH3E', 'C   ', 'CM  ', 'CUA1',
     & 'CUA2', 'CUY1', 'CUY2', 'CUA3', 'C5R ', 'C6R ', 'C5RE', 'C6RE',
     & 'CR55', 'CR56', 'CR66', 'C5RP', 'C6RP', 'N5RP', 'N   ', 'NP  ',
     & 'NX  ', 'N5R ', 'N6R ', 'NT  ', 'NC  ', 'NO2 ', 'N6RP', 'O   ',
     & 'OA  ', 'OK  ', 'OC  ', 'DUM ', 'OT  ', 'OW  ', 'OH2 ', 'OM  ',
     & 'OS  ', 'OE  ', 'OAC ', 'O5R ', 'O6R ', 'DUM ', 'OSI ', 'O2M ',
     & 'OSH ', 'DUM ', 'DUM ', 'PT  ', 'PO3 ', 'PO4 ', 'PUA1', 'P6R ',
     & 'DUM ', 'DUM ', 'DUM ', 'DUM ', 'DUM ', 'ST  ', 'SH1E', 'S5R ',
     & 'S6R ', 'SE  ', 'SK  ', 'SO1 ', 'SO2 ', 'SO3 ', 'SO4 ', 'MLI ',
     & 'MNA ', 'MMG ', 'MK  ', 'MCA ', 'MMN ', 'MFE ', 'MZN ', 'MRB ',
     & 'MCS ', 'MSI ', 'MAL ', 'XF  ', 'XCL ', 'XBR ', 'XI  ', 'MCU ',
     & 'MV  ', 'MCR ', 'MCO ', 'MNI ', 'MAS ', 'MSE ', 'MSR ', 'MY  ',
     & 'MZR ', 'MNB ', 'MMO ', 'MRU ', 'MRH ', 'MPD ', 'MAG ', 'MCD ',
     & 'MSN ', 'MSB ', 'MBA ', 'MW  ', 'MOS ', 'MPT ', 'MAU ', 'MHG ',
     & 'MPB ', 'MBI ', 'MLA ', 'MCE ', 'MPR ', 'MAC ', 'MTH ', 'MU  ',
     & 'MTE ', 'MPO ', 'XAT ', 'MTC ', 'MSC ', 'MTI ', 'MGA ', 'MGE ',
     & 'MIN ', 'MHF ', 'MTA ', 'MIR ', 'MTL ', 'MFR ', 'MRA ', 'MND ',
     & 'MPM ', 'MSM ', 'MEU ', 'MGD ', 'MTB ', 'MDY ', 'MHO ', 'MER ',
     & 'MTM ', 'MYB ', 'MLU ', 'MPA ', 'MNP ', 'MPU ', 'MAM ', 'MCM '/
      data (chmsf(i),i=161,235)/
     & 'MBK ', 'MCF ', 'MES ', 'MFM ', 'MMD ', 'MNO ', 'MLR ', 'HE  ',
     & 'NE  ', 'AR  ', 'KR  ', 'XE  ', 'RN  ', 'MRF ', 'MHA ', 'MRE ',
     & 'MSIU', 'DUM ', 'DUM ', 'NR56', 'NR66', 'NR55', 'NR1 ', 'NR2 ',
     & 'NR3 ', 'NC2 ', 'DUM ', 'DUM ', 'DUM ', 'C3  ', 'CT3 ', 'C4  ',
     & 'CT4 ', 'C6RQ', 'CQ66', 'CS66', 'CPH1', 'CPH2', 'CP3 ', 'XCLU',
     & 'CU1E', 'CU2E', 'CW1E', 'CW2E', 'NP1E', 'N5RE', 'NT1E', 'NT2E',
     & 'NT3E', 'NC2E', 'OI1E', 'OT1E', 'OH2E', 'DUM ', 'DUM ', 'DUM ',
     & 'DUM ', 'DUM ', 'DUM ', 'N6RE', 'NP2E', 'NC1E', 'DUM ', 'DUM ',
     & 'DUM ', 'DUM ', 'DUM ', 'DUM ', 'DUM ', 'CF1 ', 'CF2 ', 'CF3 ',
     & 'PUY1', 'N3  ', 'DUM '/

      data spnam/
     &  'P1     ','P-1    ','P2     ','P21    ','C2     ','Pm     ',
     &  'Pc     ','Cm     ','Cc     ','P2/m   ','P21/m  ','C2/m   ',
     &  'P2/c   ','P21/c  ','C2/c   ','P222   ','P2221  ','P21212 ',
     &  'P212121','C2221  ','C222   ','F222   ','I222   ','I212121',
     &  'Pmm2   ','Pmc21  ','Pcc2   ','Pma2   ','Pca21  ','Pnc2   ',
     &  'Pmn21  ','Pba2   ','Pna21  ','Pnn2   ','Cmm2   ','Cmc21  ',
     &  'Ccc2   ','Amm2   ','Abm2   ','Ama2   ','Aba2   ','Fmm2   ',
     &  'Fdd2   ','Imm2   ','Iba2   ','Ima2   ','Pmmm   ','Pnnn   ',
     &  'Pccm   ','Pban   ','Pmma   ','Pnna   ','Pmna   ','Pcca   ',
     &  'Pbam   ','Pccn   ','Pbcm   ','Pnnm   ','Pmmn   ','Pbcn   ',
     &  'Pbca   ','Pnma   ','Cmcm   ','Cmca   ','Cmmm   ','Cccm   ',
     &  'Cmma   ','Ccca   ','Fmmm   ','Fddd   ','Immm   ','Ibam   ',
     &  'Ibca   ','Imma   ','P4     ','P41    ','P42    ','P43    ',
     &  'I4     ','I41    ','P-4    ','I-4    ','P4/m   ','P42/m  ',
     &  'P4/n   ','P42/n  ','I4/m   ','I41/a  ','P422   ','P4212  ',
     &  'P4122  ','P41212 ','P4222  ','P42212 ','P4322  ','P43212 ',
     &  'I422   ','I4122  ','P4mm   ','P4bm   ','P42cm  ','P42nm  ',
     &  'P4cc   ','P4nc   ','P42mc  ','P42bc  ','I4mm   ','I4cm   ',
     &  'I41md  ','I41cd  ','P-42m  ','P-42c  ','P-421m ','P-421c ',
     &  'P-4m2  ','P-4c2  ','P-4b2  ','P-4n2  ','I-4m2  ','I-4c2  ',
     &  'I-42m  ','I-42d  ','P4/mmm ','P4/mmc ','P4/nbm ','P4/nnc ',
     &  'P4/mbm ','P4/mnc ','P4/nmm ','P4/ncc ','P42/mmc','P42/mcm',
     &  'P42/nbc','P42/nnm','P42/mbc','P42/mnm','P42/nmc','P42/ncm',
     &  'I4/mmm ','I4/mcm ','I41/amd','I41/acd','P3     ','P31    ',
     &  'P32    ','R3     ','P-3    ','R-3    ','P312   ','P321   ',
     &  'P3112  ','P3121  ','P3212  ','P3221  ','R32    ','P3M1   ',
     &  'P31M   ','P3C1   ','P31C   ','R3M    ','R3C    ','P-31M  ',
     &  'P-31C  ','P-3M1  ','P-3C1  ','R-3M   ','R-3C   ','P6     ',
     &  'P61    ','P65    ','P62    ','P64    ','P63    ','P-6    ',
     &  'P6/M   ','P63/M  ','P622   ','P6122  ','P6522  ','P6222  ',
     &  'P6422  ','P6322  ','P6mm   ','P6cc   ','P63cm  ','P63mc  ',
     &  'P-6m2  ','P-6c2  ','P-62m  ','P-62c  ','P6/mmm ','P6/mcc ',
     &  'P63/mcm','P63/mmc','P23    ','F23    ','I23    ','P213   ',
     &  'I213   ','PM3    ','Pn3    ','Fm3    ','Fd3    ','Im3    ',
     &  'Pa3    ','Ia3    ','P432   ','P      ','F432   ','F4132  ',
     &  'I432   ','P4332  ','P4132  ','I4132  ','P-43M  ','F-43M  ',
     &  'I43M   ','P-43N  ','F-43C  ','I-43D  ','Pm3m   ','Pn3n   ',
     &  'Pm3n   ','Pn3m   ','Fm3m   ','Fm3c   ','Fd3m   ','Fd3c   ',
     &  'Im3m   ','Ia3d   ','P21/n  ','R3     ','R-3    ','R32    ',
     &  'R3M    ','R3C    ','R-3M   ','R-3C   '/

      data valcol / -0.1d0, -0.05d0, 0.0d0, 0.05d0, 0.1d0 /
      data ((colmap(i,j),i=1,3),j=1,5) /
     & 1.0d0,0.0d0,0.0d0, 1.0d0,1.0d0,0.0d0, 0.0d0,1.0d0,0.0d0,
     & 0.0d0,1.0d0,1.0d0, 0.0d0,0.0d0,1.0d0 /

c                 Li+   Na+  Mg+  K+    Ca+  Mn2+,Fe+   Zn+  Rb+   Cs+
      data ianexc / 3,  11,  12,  19,   20,  25,  26,   30,   37,  55/
      data qexcl /1.0, 1.0, 2.0, 1.0,  2.0,  2.0, 2.0, 2.0,  1.0, 1.0/


      esc    = char(27)
      etx    = char(3)
      lf     = char(10)
      cr     = char(13)
      eot    = char(4)
      fs     = char(28)
      gs     = char(29)
      ff     = char(12)
      us     = char(31)

      mxorb = maxorb
      mx3d = max3d
      mx3d2 = max3d2
      mxm3d = maxm3d
      mxnat = numat1
      mapopt = 0
      muls   = 1
      irdogl = 0
      ivdwpl = 0
      iesp = 0
      iproj = 1
      istyp = 0
      ibgcol = 0
      ibgclo = 13
      ibgmod = 1
      vander(30) = 1.00d0
      isimpl = 0
      istaro = 0
      movfil = 'movie.dat'
      vrmlfil = 'molden.wrl'
      povfil = 'molden.pov'
      oglfil = 'molden.ogl'
      call parogf(oglfil,linlen(oglfil))
      mapfil = '3dgridfile'
      grdfil = '3dgridfile'//char(0)
      jname = 'gauss'
      qname = 'long'
      rungam = 'edtest'
      ivtwo = 0
      noshad = 0
      ihand = 0
      ixyz = 0
      iixyz = 0
      idoi = 1
      iwropt = 0
      tokcal = 627.5095d0
      toang  = 0.52917706d0
      do i=1,mxel
         vrad(i) = vrad(i) / toang
      end do
      helrod = 1.0d0
      ihtype = 1
      vrad(100) = .0d0
      moddma = 0
      npts = mx3d
      npts1 = npts
      npts2 = npts
      npts3 = npts
      space = 0
      mxzat = maxat
      idelx = 1
      uscl = 1.0d0
      thetm = 1.0d-10
      phim = 1.0d-10
      theta = thetm
      phi = phim

      grdw = 0.6d0
      grdwdf = 0.6d0
      ctval(1) = 0.0001d0
      ctval(2) = 0.00001d0
      ctval(3) = 0.000001d0
      nvalc = 3
      nspts = max3d
      edge = 3.0d0
      hon = .true.
      hdmin = 1.50d0/toang
      hdmax = 3.15d0/toang
      hamin = 145.0e0
      hamax = 215.0e0

      npix = 100
      ivadd = 0
      icalcd = 0
      icalc = icalcd
      icubtp = 0
      iqon = 0
      iqopt = 0
      ibatch = 1
      logo   = 1
      ibell = 1
      iwxyz = 0
      atcol  = 1
      fancy  = 0
      dolabs = 0
      shade  = 1
      persp = 0
      oempty = .false.
      iforscd = 1
      istdbd = 0
      iuseab = 0
      idoest = 0
      ctoz = .true.
      rdhb = .true.
      forpdb = .false.
      fortnk = .false.
      iaddh  = 1
      colpsp = .false.
      monops = .false.
      lapdbl = .false.
      lapbox = 0
      doxmol = .false.
      iold = 0
      igztyp = 0
      ioxyz = 0
      frmul = 0.1e0
      frthr = 0.0d0
      wmovie = .false.
      iwopt = 4
      rdmult = .false.
      hesend = .false.
      dovdw = .false.
      befo = .false.
      idebug = 0
      natoms = 0
      idofil = 1
      iparfl = 0
      iball = 1
      ideltm = 0

      tnkprg = 0
      rmsgrd = 0.01d0
      tnkbg = 0
      tnkarc = 0
      tnkarf = 1
      tnkit = 999
      tnknm = 'molin'

      call setmf(0)
      call gcres
      impas = 1
      ialfa = 0
      icoars = 1
      c0 = 13.d0
      smag   = 1.0d0
      coloff = -1.0
      notitle = .false.
      numcir = 8
      nquad = 6
      ipmfm = 0
      ipmfh = 1
      ifogl = 0
      ifdogl = 0
      ishad = 1

      nxtmf = 0
      ioadd = 0
      namols = 1
      imulm = 0
      nmulm = 0
      ioatms = 0
      nstrt = 1
      iplwin = 0
      plfile = 'plot'
      iplot = 3
      iun1 = 5
      iun2 = 30
      iun3 = 6
      iun4 = 15
      iun5 = 77
      idoh = 0
      pol = 0.115d0
      pol2 = -0.028d0
      ifullg = 1
      ifntcl = 15
      zsev = .true.
      dompi = .false.

      call chktmp
      call pargeo

      call parptr(5,fdum,fdum,ipoints)
      call parptr(6,bl,fdum,idum)
      call parptr(7,coo,fdum,idum)
      call parptr(8,vdwr,fdum,idum)
      call parptr(9,fdum,fdum,nz)
      call parptr(51,rphi,fdum,idum)
      call parptr(11,fdum,fdum,iuhf)
      call parptr(15,fdum,fdum,iatoms)
      call parptr(17,fdum,fdum,iftyp)
      call parptr(19,fdum,fdum,iesp)
      call parptr(21,fdum,fdum,irtval)
      call parptr(22,fdum,fdum,ifntcl)
      call parptr(24,fdum,fdum,icolps)
      call parptr(25,fdum,fdum,ivtwo)
      call parptr(26,fdum,fdum,ixyz)
      call parptr(27,fdum,fdum,iixyz)
      call parptr(28,fdum,fdum,iwropt)
      call parptr(29,fdum,fdum,ito)
      call parptr(37,vander,fdum,idum)
      call parptr(38,fdum,fdum,ispd)
      call parptr(39,fdum,fdum,ibgcol)
      call parptr(40,fdum,fdum,natc)
      call parptr(41,fdum,fdum,jring)
      call parptr(42,fdum,fdum,isimpl)
      call parptr(43,edge,fdum,idum)
      call parptr(44,fdum,fdum,ivdwpl)
      call parptr(45,fdum,fdum,natoms)
      call parptr(46,fdum,fdum,space)
      call parptr(47,qat,fdum,idum)
      call parptr(48,fdum,fdum,iff)
      call parptr(49,fdum,fdum,ityp)
      call parptr(50,fdum,fdum,ipdbt)
      call parptr(53,fdum,fdum,fancy)
      call parptr(52,fdum,fdum,istaro)
      call parptr(54,fdum,fdum,atcol)
      call parptr(55,fdum,fdum,shade)
      call parptr(56,fdum,fdum,persp)
      call parptr(57,fdum,fdum,tnkprg)
      call parptr(58,rmsgrd,fdum,idum)
      call parptr(59,fdum,fdum,tnkbg)
      call parptr(60,fdum,fdum,tnkarc)
      call parptr(61,fdum,fdum,tnkarf)
      call parptr(62,fdum,fdum,tnkit)
      call parptr(63,fdum,fdum,iqopt)
      call parptr(64,fdum,fdum,icalc)
      call parptr(65,fdum,fdum,ipdbon)
      call parptr(66,fdum,fdum,icubtp)
      call parptr(67,fdum,fdum,normc)
      call parptr(68,frmul,fdum,idum)
      call parptr(69,fdum,fdum,ibell)
      call parptr(70,fdum,fdum,natorg)
      call parptr(71,fdum,fdum,ihaszm)
      call parptr(72,fdum,fdum,mxzat)
      call parptr(73,cstoc,fdum,idum)
      call parptr(74,qd,fdum,idum)
      call parptr(2,eiga,focc,ncols)
      call parptr(10,eigb,focb,ncolb)
      call parptr(75,vectrs,vectrb,idum)
      call parptr(76,fdum,fdum,nocc)
      call parptr(77,fdum,fdum,nocb)
      call parptr(78,p,paa,idum)
      call parptr(79,fdum,fdum,istos)
      call parptr(80,stoalfa,stobnorm,naorbs)
      call parptr(81,averag,fdum,idum)
      call parptr(82,fdum,fdum,mxorb)
      call parptr(83,denn,pmnn,mx3d)
      call parptr(84,dens,denst,mx3d2)
      call parptr(85,edx,edy,iedlog)
      call parptr(86,rz,bucket,ix)
      call parptr(87,fdum,fdum,iy)
      call parptr(88,fdum,fdum,mx3d)
      call parptr(89,fdum,fdum,mx3d2)
      call parptr(90,fdum,fdum,irtcel)
      call parptr(91,fdum,fdum,idelx)
      call parptr(92,uscl,fdum,idum)
      call parptr(93,fdum,fdum,iball)
      call parptr(94,fdum,fdum,npix)
      call parptr(95,fdum,fdum,iclon)
      call parptr(96,fdum,fdum,igztyp)
      call parptr(97,fdum,fdum,igtfil)
      call parptr(98,fmap,fdum,mxm3d)
      call parptr(99,fdum,fdum,mapit)
      call parptr(100,fdum,fdum,ipsi)
      call parptr(101,valcol,fdum,idum)
      call parptr(102,fdum,fdum,mapopt)
      call parptr(110,fdum,fdum,ipdbwh)
      call parptr(113,fdum,fdum,ideltm)
      call parptr(121,fdum,fdum,ialtyp)
      call parptr(122,fdum,fdum,iscst)
      call parptr(126,fdum,fdum,ipart)
      call parptr(130,rx,fdum,idum)
      call parptr(131,xv,fdum,idum)
      call parptr(132,fdum,fdum,logo)
      call parptr(133,fdum,fdum,nscnd)
      call parptr(134,fdum,fdum,dolabs)
      call parptr(136,fdum,fdum,fyesno)
      call parptr(137,fdum,fdum,ifogl)
      call parptr(138,fdum,fdum,npts1)
      call parptr(139,fdum,fdum,backb)
      call parptr(142,xyz,fdum,idum)
      call parptr(146,fdum,fdum,ipsa)
c     parse common /b/ to parptr 143
      call parcom
      call parptr(144,vo,fdum,idum)
      call parptr(145,scale,fdum,idum)
      call parptr(148,fdum,fdum,ifdogl)
      call parptr(149,fdum,fdum,isurf)
      call parptr(150,colmap,valcol,idum)
      call parptr(156,fdum,fdum,iocnt)
      call parptr(157,rdm,fdum,idum)
      call parptr(158,fdum,fdum,mxnat)
      call parptr(159,fdum,fdum,lwrit)
      call parptr(160,fdum,fdum,lring)
      call parptr(161,fdum,fdum,inat)
      call parptr(162,fdum,fdum,ibgclo)
      call parptr(163,dipo,fdum,idum)
      call parptr(164,fdum,fdum,ibgmod)
      call parptr(165,rng1,rng2,idum)
      call parptr(166,vlcnt,vlcnt2,idum)
      call parptr(167,scal,fdum,idum)
      call parptr(168,fdum,fdum,isndcl)
      call parptr(169,fdum,fdum,idcur)
      call parptr(170,fdum,fdum,idoh)
      call parptr(171,hdmin,fdum,idum)
      call parptr(172,fdum,fdum,itot)
      call parptr(173,fdum,fdum,do3d)
      call parptr(174,fdum,fdum,iambch)
      call parptr(180,fdum,fdum,muls)
      call parptr(182,pol,fdum,idum)
      call parptr(183,pol2,fdum,idum)
      call parptr(184,fdum,fdum,icpsa)
      call parptr(185,fdum,fdum,idtpsa)
      call parptr(186,fdum,fdum,movie)
      call parptr(187,fdum,fdum,ifullg)
      call parptr(188,fdum,fdum,icst)
      call parptr(189,fdum,fdum,ibox)
      call parptr(190,fdum,fdum,igfmap)
      call parptr(191,fdum,fdum,ihashz)
      call parptr(192,shlnuc(1),fdum,idum)
      call parptr(193,fdum,fdum,iresrd)
      call parptr(194,fdum,fdum,nalign)
      call parptr(196,fdum,fdum,ishad)
      call parptr(197,fdum,fdum,ihasex)
      call parptr(198,fdum,fdum,noshad)
      call parsfn(vrmlfil,linlen(vrmlfil),2)
      call parstr(glin1,0)
      call parstr(glin2,1)
      call parstr(gtitl,2)
      call parstr(jname,3)
      call parstr(qname,4)
      call parstr(rungam,5)
      call parstr(tnknm,6)
      call parstr(vrmlfil,7)
      call parstr(mapfil,8)
      call parstr(oglfil,9)
      call parstr(povfil,10)
      call parstr(grdfil,12)

      do i=1,mxmm3
         call parsfn(mm3(i),19,4)
      end do
      do i=1,mxchtp
         call parsfn(chmtnk(i),20,5)
      end do
      do i=1,mxmol2
         call parsfn(mol2(i),5,6)
      end do
      do i=1,mxmsf
         call parsfn(chmsf(i),4,7)
      end do
      do i=1,mxamb
         call parsfn(ambstr(i),20,8)
      end do
      do i=1,mxsg
         call parsfn(spnam(i),7,10)
      end do
      do i=1,numcal
         call parsfn(achain(i),1,11)
      end do
      do i=1,mxgmx
         call parsfn(grogms(i),35,14)
      end do
      do i=1,mxgmx2
         call parsfn(grog2s(i),35,15)
      end do
      do i=1,mxg43
         call parsfn(gro43s(i),35,16)
      end do
      do i=1,mxamo
         call parsfn(amostr(i),20,17)
      end do
      do i=1,mxgff
         call parsfn(gffext(i),28,19)
      end do

      iredir = 0
      nargs = iargc()
      call getarg(0,liris)
      call parsfn(liris,linlen(liris),20)

      ntargs = 0
      n = 0
      if (nargs.gt.0) then
2013         n = n + 1
             call getarg(n,liris)
             if (liris(1:1).eq.'-') then
                if (liris(1:2).eq.'-5') ido5d = 1
                if (liris(1:2).eq.'-7') ido7f = 1
                if (liris(1:2).eq.'-a') ctoz = .false.
                if (liris(1:2).eq.'-b') befo = .true.
                if (liris(1:2).eq.'-d') idebug = 1
                if (liris(1:2).eq.'-e') idoest = 1
                if (liris(1:2).eq.'-v') idebug = 1
                if (liris(1:2).eq.'-f') iforscd = 0
                if (liris(1:5).eq.'-geom') then
                   if (gargpl('-geom',n,liris,argstr)) then
                       call parsfn(argstr,linlen(argstr),9)
                   endif
                else
                   if (liris(1:2).eq.'-g') rdhb = .false.
                endif
                if (liris(1:2).eq.'-h') then
                   if (gargpl('-hdmin',n,liris,argstr)) then
                       hdmin = reada(argstr,1,len(argstr))/toang
                   elseif (gargpl('-hdmax',n,liris,argstr)) then
                       hdmax = reada(argstr,1,len(argstr))/toang
                   elseif (gargpl('-hamin',n,liris,argstr)) then
                       hamin = reada(argstr,1,len(argstr))
                   elseif (gargpl('-hamax',n,liris,argstr)) then
                       hamax = reada(argstr,1,len(argstr))
                   elseif (liris(1:5).eq.'-hoff') then
                       hon = .false.
                   else
                       call prtflg
                       stop
                   endif
                endif
                if (liris(1:2).eq.'-q') icalcd = 1
                if (liris(1:2).eq.'-l') logo = 0
                if (liris(1:2).eq.'-n') iaddh = 0
                if (liris(1:2).eq.'-m') ibell = 0
                if (liris(1:2).eq.'-t') impas = 2
                if (liris(1:2).eq.'-z') then
                    numcir = 16
                    nquad = 12
                endif
                if (liris(1:2).eq.'-A') isimpl = 1
                if (liris(1:2).eq.'-B') isimpl = 2
                if (liris(1:2).eq.'-C') colpsp = .true.
                if (liris(1:2).eq.'-E') iuseab = 1
                if (gargpl('-G',n,liris,argstr)) 
     &              grdwdf = reada(argstr,1,len(argstr))
                if (liris(1:2).eq.'-F') ifullg = 0
                if (liris(1:2).eq.'-H') hesend = .true.
                if (liris(1:2).eq.'-I') ishad  = 0
                if (liris(1:2).eq.'-K') ivadd = 1
                if (liris(1:2).eq.'-M') then
                   if (liris(1:3).eq.'-MB') then
                      logo = 0
                      notitle = .true.
                   else
                      monops = .true.
                   endif
                endif
                if (liris(1:2).eq.'-N') dompi = .true.
                if (liris(1:2).eq.'-L') lapdbl = .true.
                if (liris(1:2).eq.'-1') lapbox = 1
                if (liris(1:2).eq.'-2') lapbox = 2
                if (liris(1:2).eq.'-P') forpdb = .true.

c             use -R int for resolution during density plots

                if (gargpl('-R',n,liris,argstr)) then
                   npts = reada(argstr,1,len(argstr))
                   if (npts.lt.mx3d) then
                      call inferr('changed nuber of pts',0)
                      npts1 = npts
                      npts2 = npts
                      npts3 = npts
                   else
                      call inferr('npts exceeds maximum',0)
                      npts = mx3d
                   endif
                endif

                if (liris(1:2).eq.'-S') shade = 0
                if (liris(1:2).eq.'-T') fortnk = .true.
                if (liris(1:2).eq.'-X') doxmol = .true.
                if (liris(1:2).eq.'-u') iold = 1
                if (liris(1:2).eq.'-=') igztyp = 1
                if (liris(1:2).eq.'-U') noshad = 1
                if (liris(1:2).eq.'-W') ivtwo = 1
                if (liris(1:2).eq.'-Z') isimpl = 3
                if (liris(1:2).eq.'-O') muls = 0
                if (liris(1:2).eq.'-Q') istaro = 1
                if (gargpl('-i',n,liris,argstr)) then
                   istdbd = reada(argstr,1,len(argstr))
                   if (istdbd.gt.2) istdbd=0
                endif
                if (gargpl('-j',n,liris,argstr)) then
                   npix = reada(argstr,1,len(argstr))
                endif
                if (gargpl('-k',n,liris,argstr)) then
                   ifntcl = reada(argstr,1,len(argstr))
                   if (ifntcl.gt.15) ifntcl=15
                endif
                if (gargpl('-p',n,liris,argstr)) 
     &             c0 = reada(argstr,1,len(argstr))
                if (gargpl('-s',n,liris,argstr)) then
                   frmul = frmul*reada(argstr,1,len(argstr))
                endif
                if (gargpl('-c',n,liris,argstr)) then
                   coloff = reada(argstr,1,len(argstr))
                   if (coloff.lt.0.0) coloff = 0.0
                   if (coloff.gt.1.0) coloff = 1.0
                endif
                if (gargpl('-o',n,liris,argstr)) plfile = argstr
                if (gargpl('-V',n,liris,argstr)) then
                   idvrml = 1
                   vrmlfil = argstr
                   call parsfn(vrmlfil,linlen(vrmlfil),2)
                endif
                if (gargpl('-r',n,liris,argstr)) then
                   parfil = argstr
                   call getpar(parfil,1,iball,ideltm)
                   iparfl = 1
                endif
                if (gargpl('-x',n,liris,argstr)) then
                   basfile = argstr
                   call rdexbas(basfile)
                endif
                if (gargpl('-y',n,liris,argstr)) then
                   frthr = reada(argstr,1,len(argstr))
                endif
                if (gargpl('-w',n,liris,argstr)) then
                   wmovie = .true.
                   itmp = reada(argstr,1,len(argstr))
                   if (itmp.eq.1) then
                      iwopt = 4
                   elseif (itmp.eq.2) then
                      iwopt = 3
                   elseif (itmp.eq.3) then
                      iwopt = 5
                      ivtwo = 1
                      movfil = 'movie.wrl'
                   elseif (itmp.eq.4) then
                      iwopt = 6
                   endif
                endif
                if (gargpl('-D',n,liris,argstr)) then
                   moddma = reada(argstr,1,len(argstr))
                endif
                goto 2012
             else
                ntargs = ntargs + 1
             endif

             
c             lenf = index(liris,' ') - 1
             lenf = linlen(liris)
             fniun = liris(1:lenf)

             call parfns(fniun,lenf)
             if (odos(fniun)) call tounx

             if (ntargs.eq.1) then
                if (index(fniun,'.ogl').ne.0) then
                   call parsfn(liris,lenf,0)
                   irdogl = 1
                   oempty = .true.
                   iftyp = 0
                else
                   if (opfil(48,fniun,lenf,1,1,0)) then
                      iun1 = 48
                      call parsfn(liris,lenf,0)
                   else
                      oempty = .true.
                      iftyp = 0
                   endif
                endif
             else
                if (iredir.eq.1) then
                   if (opfil(49,fniun,lenf,1,0,0)) iun3 = 49
                   goto 10
                endif
                if (liris(1:2).eq.'> ') then
                   iredir = 1
                else
                   call parsfn(liris,lenf,0)
                endif
             endif
2012         continue
          if (n.lt.nargs) goto 2013
          if (ntargs.eq.0) oempty = .true.
      else
          oempty = .true.
      endif

      if (iparfl.eq.0) then
c UNIX
          call getenv('HOME',hom)
          if (hom(1:1).ne.' ') then
              itlen = linlen(hom)
              call getpar(hom(1:itlen)//'/.moldenrc',0,iball,
     &             ideltm)
          endif
c PC
c          call getpar('.moldenrc',0,iball,ideltm)
      endif

      if (dompi) call chkmpi

c     set unit numbers for read and write

10    iun2 = 30

c---- set default plotter tektronix

      istpdb = iforscd
      if (ioadd.eq.0) then
         call inirot
         call setis(4,0)
         call setarr(10,idum,idum)
         call setarr(11,idum,idum)
         backb = 0
      endif

      iscupd = 1
      ihsdp = 0
      idipon = 0
      fyesno = 0
      idirogl = 0
      igtfil = 0
      grdw = grdwdf
      mapit = 0
      icubtp = 0
      icalc = icalcd
      iloop = 0
      ifreq = 0
      isrf = 0
      ispd = 0
      isgau = 1
      ido5d = 0
      ido7f = 0
      ido9g = 0
      ihasq = 0
      ihasesp = 0
      ihasd = 0
      ihasg = 0
      itotc = 0
      imult = 1
      iopt = 1
      iesp = 0
      nesp = 0
      isfdon = 0
      isftyp = -1
      mopopt = 0
      domul = .true.
      onofil = .false.
      oeerst = .true.
      domax = .false.
      first = .true.
      idisml= 0
      domol = .true.
      iuhf   = 0
      yes   = .false.
      iftyp = 0
      odeflt= .false.
      ipsprt = 0
      wmolf = .false.
      wxyz = .false.
      if (ioadd.eq.0) ipdbon = 0
      if (forpdb) ipdbon = 1
      ifd = 1
      iqopt = 0
      ofrst = .true.
      movie = 0
      molonl = .false.
      normc = 0
      forces = .false.
      dozme = .false.
      ircus = 0
      irtype = 0
      iconv = 1
      icst = 0
      ibox = 0
      icell = 0
      igfmap = 1
      ihashz = 0
      ihsnmr = 0
      iresrd = 0
      nalign = 0
      ihashb = 0
      nexchg = 0

      isl = 0
      do i=1,3
         sl(i) = 1.0d0
      end do

      do i=1,mxheta
         hetz(i) = '   '
         ihetq(i) = 0
         ihqset(i) = 0
         ihhadd(i) = 0
      end do

c     set ipoints and ngeoms in setarr

      call setarr(15,idum,idum)
      r(1)  = 1.0d0
      r(2)  = 1.0d0
      r(3)  = 1.0d0
      inact = 0
      if (ioadd.eq.0) then
         iatoms = 0
         natoms = 0
      endif
      norbs = 0
      iplat = 0
      ihaszm = 0
      natorg = 0
      istsurf = 0
      srfmap = 0
      srfloc = 1
      iclon = 0
      nz = 0
      idvrml = 0
      idocub = 0
      dosrf2 = .false.
      irtcel = 0
      itz = 0
      irtval = 0
      ndm = 0
      iff = 0
      isofst = 1
      edget = 0.0d0
      ctvalt = 0.0d0
      nsptst = 0
      ichh = 0
      ichm = 0
      ichl = 0
      imh = 1
      imm = 1
      iml = 1
      ipseud = 0
      ioni = 0
      ioniad = 0
      nscnd = 0
      iscst = 0
      ialtyp = 0
      ipart = 0
      idcur = 1
      imn = 0
      imx = 0
      if (nxtmf.eq.0) call setmf(0)
      ipdbgro = 0
      itsrf = 0
      ipsa = 0
      icpsa = 0
      idtpsa = 0
      itot = 0
c default dummy color
      icol(100) = 10
      iambch = 1

      call setpp(idoest)

      i = 0
      if (oempty) then
         isbin = 1
         oempty = .false.
         goto 2828
      endif

      if (doxmol.or.ipdbon.eq.1) goto 2828

c     use defaults

      read(iun1,'(a)',end=2828,err=2828) line
      isbin = 0
      if (obin(line)) isbin = 1
      if (isbin.eq.1) goto 2828
      if (line(1:6).eq.'HEADER') ipdbon = 1
      if (line(1:6).eq.'HETATM') ipdbon = 1
      if (line(1:6).eq.'REMARK') ipdbon = 1
      if (line(1:6).eq.'CRYST1') ipdbon = 1
      if (line(1:6).eq.'ATOM  ') ipdbon = 1
      if (line(1:6).eq.'TITLE ') ipdbon = 1
      if (line(1:6).eq.'COMPND') ipdbon = 1
      if (line(1:8).eq.'[AMBFOR]') fortnk = .true.
      if (line(2:9).eq.'[AMBFOR]') fortnk = .true.
      if (line(1:7).eq.'[AMBMD]'.or.line(2:8).eq.'[AMBMD]') then
          k = index(line,'steps')
          if (k.ne.0) then
             nstruc = int(reada(line,k+1,len(line)))
             if (nstruc.gt.MAXPNT) then
                call allgeo(nstruc,0)
             endif
          endif
          fortnk = .true.
      endif

      do i=1,14
         read(iun1,'(a)',end=2828,err=2828) line
      end do

2828  if ((i.le.1.or.i.ge.5).or.ipdbon.eq.1.or.isbin.eq.1
     &     .or.doxmol) then
c
c  use defaults (xwindows only) to make keywords line
c
         if (.not.(ipdbon.eq.1.or.isbin.eq.1.or.doxmol)) then
            if (i.eq.0) isbin = 1
         endif
         odeflt = .true.
         idisml = 1
         title  = 'defaults used'
         keyori = 'homo 3d contour xwindows'
         keywrd = 'HOMO 3D CONTOUR XWINDOWS'
         iplot  = 6
         iun2   = iun1
         if (isbin.eq.0) rewind iun2
      else
c
c  read title and keywords
c
         rewind iun1
         read(iun1,'(a)') title
         write(iun3,'(//1x,a)') title

         read(iun1,'(a)') keyhlp
         write(iun3,'(//1x,a)') keyhlp
         keyori(1:80) = keyhlp
         call tocap(keyhlp,80)
         keywrd(1:80) = keyhlp

         read(iun1,'(a)',end=4) keyhlp
         write(iun3,'(1x,a)') keyhlp
         keyori(81:160) = keyhlp
         call tocap(keyhlp,80)
         keywrd(81:160) = keyhlp

         read(iun1,'(a)',end=4) keyhlp
         write(iun3,'(1x,a)') keyhlp
         keyori(161:240) = keyhlp
         call tocap(keyhlp,80)
         keywrd(161:240) = keyhlp

         read(iun1,'(a)',end=4) keyhlp
         write(iun3,'(1x,a)') keyhlp
         keyori(241:320) = keyhlp
         call tocap(keyhlp,80)
         keywrd(241:320) = keyhlp

4        continue
c
c     determine graphics language

         call plotin
c
c     open file
c
         call files(onofil,idocub)
         if (idocub.eq.1) then
             call filmap(mapit)
             if (mapit.eq.1) then
                mapopt = 1
                n = keyr3v(keywrd,'MAPVAL',valcol,5)
             endif
         endif
         if (onofil) then
             isbin = 1
         else
             read(iun2,'(a)',end=2829,err=2829) line
2829         isbin = 0
             if (obin(line)) isbin = 1
             rewind iun2
         endif
      endif

      if (ialfa.eq.0) then
          if (iplot.eq.6) then
              call xwin(coloff,yy,1,str,nstr,icells,idum2)
          else
              open(unit=iun4,form='formatted',file=plfile,
     &             status='unknown')
          endif
      endif


      if (colpsp) then
         icolps = 1
      else
         icolps = 0
         if (ipdbon.eq.1.and..not.monops) icolps = 1
      endif

      if (idocub.eq.1) goto 1233

      if (doxmol.and.isbin.eq.0) then
         call getxyz(igetxy,heat,0)
         if (igetxy.eq.1) then
            inact = inact + 4
            domol = .true.
            iftyp = 10
            doxmol = .false.
            call inferr('XYZ arc or input file',0)
            goto 1234
         else
            iftyp = 0
            call inferr('Error in XYZ file',0)
            stop
         endif
      endif

      if (fortnk.and.isbin.eq.0) then
         call gettnk(igttnk,idebug,ipdbon,0,iheat,heat)
         if (igttnk.eq.1) then
            domol = .true.
            if (ipdbon.eq.1) then
               if (iresrd.eq.1) then
                  istpdb = 1
               else
                  istpdb = 0
               endif
               call calfa(istat,istpdb,iaddh,ioatms,nstrt,ioadd)
c hydrogen bonds to slow for massive waterbox
c               call dohcon(1)
               call setchn(1)
               if (istat.ne.1) ipdbon = 0
               inact = inact + 4
               iftyp = 12
            else
               iftyp = 11
               inact = inact + 4
            endif
            if (ichx.eq.1) iclon = 1
            call inferr('Tinker/Ambfor xyz file',0)
            fortnk = .false.
            goto 1234
         else
            iftyp = 0
            call inferr('Error in Tinker XYZ file',0)
            onofil = .true.
            isbin = 1
            fortnk = .false.
         endif
      endif

      if (onofil) then
         inact = 3
         domol = .true.
         goto 1234
      endif

1233  continue
      if (idebug.eq.0) then
         if (index(keywrd,'DEBU').ne.0) idebug = 1
      endif

      rdmult = (index(keywrd,'RDMULT').ne.0)
      molpot = (index(keywrd,'MOLP').ne.0)

      elpot  = (index(keywrd,'ELPO').ne.0)
      if (elpot.or.molpot) then
        if (keyrv(keywrd,'EXTPOSCHG',tmp(1),tmp(2),tmp(3))) then
          nexchg = nexchg + 1
          iexchg(nexchg) = 1
          do j=1,3
             exchg(j,nexchg) = tmp(j)
          end do
          print*,"Added extra positive charge at coordinates (bohr):"
          print*,(exchg(jj,nexchg),jj=1,3)
        endif
        if (keyrv(keywrd,'EXTNEGCHG',tmp(1),tmp(2),tmp(3))) then
          nexchg = nexchg + 1
          iexchg(nexchg) = -1
          do j=1,3
             exchg(j,nexchg) = tmp(j)
          end do
          print*,"Added extra negative charge at coordinates (bohr):"
          print*,(exchg(jj,nexchg),jj=1,3)
        endif
      endif

      chpot  = (index(keywrd,'CHPO').ne.0)
      espchg = (index(keywrd,'ESPCH').ne.0)
      dmachg = (index(keywrd,'DMACH').ne.0)
      doiso  = (index(keywrd,'ISODEN').ne.0)
      ovrlap = (index(keywrd,'OVER').ne.0)
      atomic = (index(keywrd,'ATOM').ne.0)
      valenc = (index(keywrd,' VAL').ne.0)
      bonds  = (index(keywrd,'BOND').ne.0)
      if (.not.befo) befo   = (index(keywrd,'BEFO').ne.0)
      molonl = (index(keywrd,'MOLONLY').ne.0)
      domax  = (index(keywrd,'MAXM').ne.0)
      if (index(keywrd,'VRML').ne.0) then
         idvrml = 1
      else
         idvrml = 0
      endif
      if (index(keywrd,'POVRAY').ne.0) then
         idvrml = 1
         ivtwo = 2
      endif
      if (index(keywrd,'MOGL').ne.0) then
         idvrml = 1
         ivtwo = 3
      endif
      if ((index(keywrd,'VRML2.0').ne.0)) ivtwo = 1
      if (mapit.eq.1.and.ivtwo.eq.0) 
     &   print*,'Warning MAPPING only works with VRML2.0 and OpenGL'
      dovdw  = (index(keywrd,'VDWS').ne.0)
      dolap  = (index(keywrd,'LAPL').ne.0)
      if (index(keywrd,'PERSPEC').ne.0) persp = 1
      zmat = keyi(keywrd,'WRZMAT',iopt)
      if (zmat) zsev = .false.
      wxyz = (index(keywrd,'WRXYZ').ne.0)
      if (wxyz) zsev = .false.
      wmolf = (index(keywrd,'WRMOLF').ne.0)
      if (wmolf) zsev = .false.
      if (index(keywrd,'ECTRUM').ne.0) zsev = .false.

      i1 = index(keywrd,'XANG')
      if (i1.ne.0) then
         xr = reada(keywrd,i1+4,len(keywrd))
         call xyzrot(-3,xr)
      endif
      i1 = index(keywrd,'YANG')
      if (i1.ne.0) then
         xr = reada(keywrd,i1+4,len(keywrd))
         call xyzrot(-2,xr)
      endif
      i1 = index(keywrd,'ZANG')
      if (i1.ne.0) then
         xr = reada(keywrd,i1+4,len(keywrd))
         call xyzrot(-1,xr)
      endif

      if (idocub.eq.1) goto 1235

      if (isbin.eq.1) then
         inact = 3
         if (iun2.eq.5) goto 1234
c     try for binary mopac .gpt file
         iftyp = 1
         close(iun2)
         if (opfil(iun2,fniun,lenf,0,0,1)) then
            rewind(iun2)
            call rdnorb(norbs,impas)
            call mopin(istat,isbin,impas)
            rewind(iun2)
         endif
1236     if (istat.eq.0) then
            iftyp = 0
            domol = .true.
            call rdmsf(idebug,istat)
            if (istat.ge.1) then
               iftyp = 6
               inact = 7
               if (istat.eq.2) then
                   call fdat(ifd,1,istdbd,iuseab,moddma,idebug)
               endif
               goto 1234
            else
               inact = 3
            endif
            goto 1234
         else
            goto 1235
         endif
      endif


c-------- atomic density data to be generated -----

      i = index(keywrd,'GENAT')
      if (i.ne.0) then
         call densmat(idebug)
         stop
      endif

      i = index(keywrd,'GENE')
      if (i.ne.0) then
         call densmto(idebug)
         stop
      endif

      if (ipdbon.eq.1) then
         if (iplot.eq.6) call inferr('PDB file',0)
         if (istpdb.eq.1) then
            call pdbstd(istat,rdhb,ioadd)
         else
            call rdpdb(istat)
         endif
         call xyzcoo(0,1,ioadd)
         if (istpdb.eq.0) call doconn
         call calfa(istat,istpdb,iaddh,ioatms,nstrt,ioadd)
         call dohcon(1)
         if (istat.ne.1) ipdbon = 0
         domol = .true.
         inact = 7
         iqon = 0
         goto 1234
      endif

      call rdchx(idebug,ifd,istdbd,iuseab,moddma,istat,0)
      if (istat.gt.0) then
          iftyp = 6
          inact = 7
          if ((icrtp.le.3.and.icrtp.ge.1).or.icrtp.eq.5) inact = 6
          if (iplot.eq.6) then
             if (icrtp.eq.1) call inferr('Chemx file',0)
             if (icrtp.eq.2) call inferr('FDAT file',0)
             if (icrtp.eq.3) call inferr('BIOSYM CAR/ARC file',0)
             if (icrtp.eq.4) call inferr('SHELX file/SPF file',0)
             if (icrtp.eq.5) call inferr('VASP file',0)
             if (icrtp.eq.6) call inferr('MSI file',0)
             if (icrtp.eq.7) call inferr('CIF file',0)
             if (icrtp.eq.8) call inferr('CONQUEST file',0)
          endif
          domol = .true.
          goto 1234
      else
         rewind iun2
      endif

c---------is it gamess, gaussian or mopac ?-----------

      call searcht(lstr,'g a m e s s','M.W.SCHMIDT','@<TRIPOS>',
     &             istat)
      if (istat.ne.0) then
         if (index(lstr,'@<TRIPOS>').ne.0) then
c try pdb .mol2 first
            if (nxtmf.eq.0) call prsgmf(1)
            call rewmf
            istat = 0
            call rdmol(idebug,ipdbon,ioadd,istat)

c istat = 0 error
c istat = -1 error need to allocate more space
c       = 1 ok
c       = 2 cell info

            if (istat.eq.3) then
               call rewmf

c in istat = -1: try non pdb .mol2 

               istat = -1
               call rdmol(idebug,ipdbon,ioadd,istat)
            endif

            if (istat.gt.0) then
               if (ipdbon.eq.1) then
                  istpdb = 1
                  call calfa(istat,istpdb,iaddh,ioatms,nstrt,ioadd)
                  call setchn(1)
                  inact = inact + 4
               else
                  iftyp = 6
                  inact = 7
               endif
               if (istat.eq.2) then
                   call fdat(ifd,1,istdbd,iuseab,moddma,idebug)
               endif
               domol = .true.
               goto 1234
            endif
         else
            iftyp = 2
            if (index(lstr,'M.W.SCHMIDT').ne.0) iftyp = 3
         endif
      else
         rewind iun2
         call seartu(line,'Gaussian System',
     &       'part of the Gaussian','Gaussian(R)',istat)
         if (istat.ne.0) then
            iftyp = 4
         else
            call searchq(line,'PROGRAM CPMD','Welcome to Q-Chem',
     &                        '* O   R   C   A *',
     &                   'Northwest Computational Chemistry',istat)
            if (istat.ne.0) then
               if (icdex(line,'PROGRAM CPMD').ne.0) then
                  iftyp = 7
               elseif (icdex(line,'O   R').ne.0) then
                  iftyp = 9
               elseif (icdex(line,'Q-Chem').ne.0) then
                  iftyp = 8
               else
                  iftyp = 15
               endif
            else
               call searchu(lstr,'[Molden Format]',istat)
               if (istat.ne.0) then
                  iftyp = 5
               else
                  rewind iun2
                  call searchu(lstr,'[Atoms]',istat)
                  if (istat.ne.0) then
                     iftyp = 5
                  else
                     iftyp = 1
                  endif
               endif
            endif
         endif
      endif

      if ((iftyp.ge.1.and.iftyp.le.5).or.iftyp.eq.7) then
         call setarr(6,idum,ioatms)
      endif

      if (iftyp.le.1) then
         rewind iun2
         call getzm(iatoms,0,0,istat)
         if (istat.eq.1) then
            iftyp = 0
            iconv = 0
            inact = inact + 4
            call inferr('Gamess/Gaussian input file',0)
            domol = .true.
            goto 1234
         endif
      endif
      rewind iun2

c------- read gamess gaussin or mopac outputfile or mopac .gpt

      if (rdmult) then
          irtype = 0
          call rdmolg(istat)
          call getmul
          if (idebug.eq.1) call mulprt
          molpot = .true.
          idisml = 1
      else
          if (iftyp.eq.1) then
            inact = 3
            mopopt = 0
            if (iplot.eq.6.or.index(keywrd,'SPECTRUM').ne.0.or.
     &           index(keywrd,'PLECTRUM').ne.0.or.
     &           index(keywrd,'WRZMAT').ne.0) then
c              call searchd(line,'MOPAC:  VERSION','VAMP',istat)
              call searchd(line,'MOPAC','VAMP',istat)
              if (istat.eq.0) then
                 rewind iun2
                 call search(line,'AMPAC Version ',istat)
                 if (istat.eq.1) mopopt = 3
              else
                 mopopt = 1
                 if (index(line,'VAMP').ne.0) mopopt = 2
c Mopac2007 aux
                 if (index(line,'START').ne.0) mopopt = 4
              endif
              if (mopopt.eq.0) then
                rewind iun2
                call getmop(iatoms,heat,0,1,istat)
                if (istat.eq.1) then
                   iconv = 0
                   inact = inact + 3
                   call inferr('Mopac arc or input file',0)
                   domol = .true.
                   goto 1234
                else
                   rewind iun2
                endif
              elseif (mopopt.ne.4) then
                inact = inact + 4
                call inferr('Mopac outfile',0)
                call mopaco(istat,mopopt,irtype)
                if (istat.eq.0) then
                   call getmop(iatoms,heat,0,1,istat)
                   if (istat.eq.1) then
                      iconv = 0
                      inact = inact - 1
                      call inferr('Mopac arc or input file',0)
                      domol = .true.
                      goto 1234
                   else
                      call inferr('No coordinates on File!',1)
                   endif
                endif
                domol = .true.
                goto 1234
              endif
            endif
            call rdnorb(norbs,impas)
            if (mopopt.eq.4) then
               call rdmaux(idebug,statio,istat)
               if (istat.eq.-1) call rdmaux(idebug,statio,istat)
               if (istat.eq.1) inact = 0
            else
               call mopin(istat,isbin,impas)
            endif
            if (istat.eq.0) then
               iftyp = 0
               rewind iun2
               call getxyz(igetxy,heat,ioadd)
               if (igetxy.eq.1) then
                  inact = 4
                  domol = .true.
                  iftyp = 10
                  call inferr('XYZ arc or input file',0)
                  goto 1234
               endif
               rewind iun2
               ircus = 1
               call getxyz(igetxy,heat,0)
               if (igetxy.eq.1) then
                  inact = 4
                  domol = .true.
                  iftyp = 10
                  call inferr('XYZ arc or input file',0)
                  goto 1234
               endif
               ircus = 0

               call gettnk(igttnk,idebug,ipdbon,0,iheat,heat)
               if (igttnk.eq.1) then
                  domol = .true.
                  if (ipdbon.eq.1) then
                     if (iresrd.eq.1) then
                        istpdb = 1
                     else
                        istpdb = 0
                     endif
                     call calfa(istat,istpdb,iaddh,ioatms,nstrt,ioadd)
c hydrogen bonds to slow for massive waterbox
c                     call dohcon(1)
                     call setchn(1)
                     if (istat.ne.1) ipdbon = 0
                     inact = 4
                     iftyp = 12
                  else
                     iftyp = 11
                     inact = 4
                  endif
                  if (ichx.eq.1) iclon = 1
                  call inferr('Tinker/Ambfor xyz file',0)
                  goto 1234
               endif
               rewind iun2
               call getmol(igetmo)
               if (igetmo.eq.1) then
                  domol = .true.
                  if (nxtmf.eq.0) call prsgmf(3)
                  call rewmf
                  isbin = 0
                  iftyp = 13
                  inact = 4
                  call inferr('mol file',0)
                  goto 1234
               endif
               ipdbgro = 0
               call rdgro(idebug,istat)
               if (istat.eq.1) then
                  if (iplot.eq.6) call inferr('GRO file',0)
                  ipdbon = 1
                  ipdbgro = 1
                  call xyzcoo(0,1,ioadd)
                  if (istpdb.eq.0) call doconn
                  call calfa(istat,istpdb,0,ioatms,nstrt,ioadd)
                  call dohcon(1)
                  if (istat.ne.1) ipdbon = 0
                  domol = .true.
                  inact = 8
                  iqon = 0
                  goto 1234
               endif
               goto 1234
            endif
          else
c========= read from gamess or gaussian outputfile ======
            if (iftyp.eq.2.or.iftyp.eq.3) then
                if (iplot.eq.6) call inferr('Gamess file',0)
                call rdnorb(norbs,impas)
                if (iftyp.eq.3) then
                   call rdgamu(idebug,befo,statio,irtype,
     &                         hesend,istat)
                else
                   call rdgam(idebug,befo,statio,ioxyz,irtype,istat)
                endif
                inact = 0
                if ((molonl.and.iplot.eq.6).or.istat.eq.0) then
                    inact = inact + 4
                    domol = .true.
                    domul = .false.
                    goto 1234
                elseif (molonl) then
                    stop 'MOLONLY is only valid with XWINDOWS'
                endif
                ihasd = 1
            endif
            if (iftyp.eq.4) then
                if (iplot.eq.6) call inferr('Gaussian file',0)
                call rdnorb(norbs,impas)
                call rdgaus(idebug,befo,statio,irtype,istat)
                if (istat.eq.0) then
                   if (iplot.eq.6) then
                       inact = inact + 4
                       domol = .true.
                       domul = .false.
                       goto 1234
                   else
                       if (.not.zsev) then
                          goto 1234
                       else
                          stop
                       endif
                   endif
                else
                   ihasd = 1
                endif
            endif

            if (iftyp.eq.7) then
               if (iplot.eq.6) 
     &              call inferr('Car-Parrinello MD (CPMD) file',0)
               write(*,'(A)')'Reading CPMD Output file.'
               write(*,'(2A)')'Interface Developed at NEST - INFM',
     &              ' laboratories.'
C developer Teodoro Laino. GPL Licence 2002...
               call rdcpmd(idebug,befo,statio,ioxyz,irtype,
     &              hesend,istat)
c read cell information...
               call cpmdcell
               
               if (istat.eq.0) then
                  if (iplot.eq.6) then
                     inact = inact + 4
                     domol = .true.
                     domul = .false.
                     goto 1234
                  else
                     stop
                  endif
               else
                  ihasd = 1
               endif
            endif

            if (iftyp.eq.5) then

               if (nxtmf.eq.0) call parsmf

               isrf = -1
               call rewmf
               call rdmolf(idebug,statio,irtype,isrf,istat)
               if (naorbs.gt.mxorb) then
                   call rdnorb(naorbs,impas)
                   call rewmf
                   call rdmolf(idebug,statio,irtype,isrf,istat)
               endif
               if (isrf.ne.-1) then
                   iesp = isrf
                   isrf = 1
               else
                   isrf = 0
               endif
               iconv = 0

               if (iplot.eq.6) call inferr('Molden Format',0)

               inact = 0

               if (ihasd.eq.0) then
                   if (iplot.eq.6) then
                       inact = inact + 4
                       domol = .true.
                       domul = .false.
                       goto 1234
                   else
                       stop
                   endif
               endif

            endif
            if (iftyp.eq.8) then

               if (nxtmf.eq.0) call prsqmf

               call rewmf
               call rdqchm(idebug,irtype,istat)
               if (naorbs.gt.mxorb) then
               endif

               if (iplot.eq.6) call inferr('Qchem file',0)

               inact = 0

               if (istat.eq.0) then
                   if (iplot.eq.6) then
                       inact = inact + 4
                       domol = .true.
                       domul = .false.
                       goto 1234
                   else
                       stop
                   endif
               else
                   ihasd = 1
               endif
            endif

            if (iftyp.eq.9) then

               if (nxtmf.eq.0) call prsgmf(2)

               call rewmf
               call rdorca(idebug,irtype,istat)

               if (iplot.eq.6) call inferr('Orca file',0)

               inact = 0

               if (istat.eq.0) then
                   if (iplot.eq.6) then
                       inact = inact + 4
                       domol = .true.
                       domul = .false.
                       goto 1234
                   else
                       stop
                   endif
               else
                   ihasd = 1
               endif
            endif

            if (iftyp.eq.15) then


               if (nxtmf.eq.0) call prsgmf(4)
               if (index(keywrd,'PLECTRUM').ne.0) then
                  call setmf(2)
               endif

               call rewmf
               call rdnwch(idebug,irtype,istat)

               if (iplot.eq.6) call inferr('NWChem file',0)

               inact = 0

               if (istat.eq.0) then
                   if (iplot.eq.6) then
                       inact = inact + 4
                       domol = .true.
                       domul = .false.
                       goto 1234
                   else
                       stop
                   endif
               else
                   ihasd = 1
               endif
            endif

          endif
       endif

1235  continue

      do ii=1,3
      if (keyirv(keywrd,'EXTPOSCOO',iatm1,iatm2,rcoo)) then
          nexchg = nexchg + 1
          iexchg(nexchg) = 1
          do i=1,3
             tmp(i) = xyz(i,iatm1) - xyz(i,iatm2)
          end do
          rc = vlen(tmp)
          do i=1,3
             tmp(i) = (tmp(i)*rcoo)/rc
          end do
          do i=1,3
             exchg(i,nexchg) = xyz(i,iatm1) + tmp(i)
          end do
          print*,"Added extra positive charge at coordinates (bohr):"
          print*,(exchg(jj,nexchg),jj=1,3)
          ij = index(keywrd,'EXTPOSCOO')
          keywrd(ij:ij+9) = '         '
      endif
      if (keyirv(keywrd,'EXTNEGNH4',iatm1,iatm2,rcoo)) then
          nexchg = nexchg + 1
          iexchg(nexchg) = -1
          do i=1,3
             tmp(i) = xyz(i,iatm1) - xyz(i,iatm2)
          end do
          rc = vlen(tmp)
          do i=1,3
             tmp(i) = (tmp(i)*rcoo)/rc
          end do
          do i=1,3
             exchg(i,nexchg) = xyz(i,iatm1) + tmp(i)
          end do
          print*,'H ',(xyz(j,iatm1),j=1,3)
          print*,'N ',(xyz(j,iatm2),j=1,3)
          print*,"Added extra positive charge at coordinates (bohr):"
          print*,(exchg(jj,nexchg),jj=1,3)
          ij = index(keywrd,'EXTNEGNH4')
          keywrd(ij:ij+9) = '         '
      endif
      end do
c
c process keywords CENTER,LINE,EDGE,HOMO,LUMO,PSI,OCCU,PLANE,CUT,
c                  ROT,LIFT
      cut=1.0d0
      if ( (bonds.or.atomic) .and.(iftyp.ne.1)) cut = 0.1
      call defrad(.true.)
      call datin(npts1,npts2,npts3)


      write(iun3,'(//25x,a)')
     &            'COORDINATES' 
      write(iun3,'(20x,a,/)')
     &            'used for orbitals/density' 
      write(iun3,50)
  50  format(20X,'X', 9X,'Y',9X,'Z',/)
      write(iun3,101)(i,elemnt(nat(i)),(xyz(j,i),j=1,3),i=1,natoms)
  101 format(i14,a2,3f12.6)

*========================================================
* geometrical preparations
*
      fine   = (index(keywrd,'FINE').ne.0)

      one=1.d0
      if (.not.keyr(keywrd,'PHASE',one)) then
         if (index(keywrd,'PHASE').ne.0) one = -1.d0
      endif

c     set up view direction for the 3-D plots

      if(index(keywrd,'3D').ne.0.or.index(keywrd,'GRID').ne.0
     &   .or.index(keywrd,'SPACE').ne.0) then
          tx=0.5d0
          ty=0.5d0*0.33333d0
          tz=0.5d0
      else
          tx=0.d0
          ty=0.d0
          tz=1.d0
      endif

      if (keyr(keywrd,'AXIS',tz)) then
          tx = 1.d0 - tz
          ty = tx*0.33333d0
      endif
      call setang(tx,ty,tz,theta,phi)
c
c plot is euclidean if tz .gt. 0.999
c
      euclid=(tz.gt.0.999d0)

c     parse orientation of the plot plane to common /eulx/

      call pareul

      if (r(1).eq.0.d0) then
         r(1) = 4.d0
         r(2) = r(1)
         r(3) = r(1)
      endif

      if (index(keywrd,'ADDMOL').ne.0) idisml = 1

      oscal = .false.
      if (keyr(keywrd,'MULT',scalf)) then
         scalf = -1.0*scalf
         oscal = .true.
      endif

      ostep = .false.
      if (keyr(keywrd,'STEP',stepf)) ostep = .true.

      do3d  = 0
      do3dx = .false.
      docnt = .false.
      if(index(keywrd,'3D').ne.0.or.index(keywrd,'GRID').ne.0)
     &              do3d = 1
      if((index(keywrd,'3D').eq.0.and.index(keywrd,'GRID').eq.0)
     &   .or.index(keywrd,'CONT').ne.0) docnt = .true.

      space = 0
      ofrst = .true.
      idofil = 0
      if (iplot.eq.6.or.iplot.eq.4) idofil = 1
      valcnt = 0.05d0
      if (keyr(keywrd,'SPACE',valcnt)) then
         idisml = 1
         space = 1
         do3d  = 0
         docnt = .false.
      endif
*========================================================
      if (rdmult) goto 22

      if (idebug.eq.1) call prtvec

      if (espchg.or.dmachg) then

         if (iftyp.eq.2.or.iftyp.eq.4.or.
     &       (iftyp.eq.5.and.isgau.eq.1)) then

            ipsi = 0
            bonds = .false.
            ovrlap = .false.
            atomic = .false.
            dolap = .false.
            call denmak(denok)
            if (dmachg) call muldma(vdwr,moddma,.true.,.true.)
            ntmp = keyr3v(keywrd,'ISODEN',ctval,3)
            if (ntmp.gt.0) nvalc = ntmp
            call espchrg(ctval,nvalc,doiso,dmachg)
            stop 'end of esp charges calculation'
         else
            stop 'esp charges only for gamess and gaussian output'
         endif
      endif

c
c----- process keyword ORIENT ------
c
      doori = .false.
      i     = index(keywrd,'ORIENT')
      if (i.ne.0) doori = .true.
      if (doori) then
         do j=i+6,320
            if (keywrd(j:j).ne.' ') goto 1199
         end do
1199     if (keywrd(j:j).eq.'=') then
            call oriin(j)
         else
            call parori
         endif
      endif
22    continue

c     make projection of atoms onto the plot plane

      call proato


1234  continue

      
     
c
c--------- plot initialisation ----------------
c
      if (ialfa.eq.0) then
         if (iplot.ne.6) call plini(iplot,.true.,icolps)
         ialfa = 1
      endif

      adjus = 1.0d0
      if ((iftyp.ge.2.and.iftyp.le.5).or.iftyp.eq.8) adjus=toang
      call parptr(154,adjus,fdum,idum)

c all workings in the molecule mode are au
c all workings in the density mode are au except mopac,
c adjus only has effect on the density mode

      if (index(keywrd,'SPECTRUM').ne.0.or.
     &    index(keywrd,'PLECTRUM').ne.0) then
         call getpoi(-1,ifd,iff,1,ioatms,ioadd)
         if (irtype.eq.4) then
            if (iftyp.eq.1) then
                call getfrm(istat)
            elseif (iftyp.eq.4) then
                call getfr(istat)
            elseif (iftyp.eq.5) then
                call getfra(istat)
                call getint(istat)
            elseif (iftyp.eq.8) then
                call getqfr(istat)
            elseif (iftyp.eq.9) then
                call getofr(istat)
            elseif (iftyp.eq.15) then
                call getnfr(istat)
            endif
         endif
         if (index(keywrd,'SPECTRUM').ne.0) then
            call pltspec(0)
         else
            call pltspec(1)
            tmpstr = "gs -sDEVICE=jpeggray -g1000x500 -dBATCH -dNOPAUSE"  
     &               // " -sOutputFile=spec.jpg spec.ps"
            call sysstr(tmpstr,linlen(tmpstr))
            if (iftyp.eq.1) then
               frmul = 0.1e0
            endif
            ifr = index(keywrd,'FRMUL=')
            if (ifr.ne.0) then
               frmul = real(reada(keywrd,ifr+6,linlen(keywrd)))
            endif
            call wrpnt("mol.xyz",7,4,0,ipoints,
     &                 fancy,atcol,dolabs,backb)
            ixyz = 9
            call wrpnt("mol.mol",7,4,0,ipoints,
     &                 fancy,atcol,dolabs,backb)
            ixyz = 0
            nfrq = nfrqs()
            do i=1,nfrq
               call setnrm(i,4,.true.,dconpt)
               idirct = 1
               iframe = nframe - 1
               iloop = 0
               do j=1,20
                  call nxtpnt(4,0,1,0,0,.true.)
               end do
            end do
            call jdxwr
            call prthtm
         endif
         stop
      endif

      if (iplot.eq.6.or.zmat.or.wxyz.or.wmolf) then
         if (ioadd.eq.0) then
            call setarr(13,idum,idum)
            backb  = 0
            ssbond = .false.
            if (ipdbon.eq.1) then
               atcol = 0
            else 
               atcol = 1
            endif
            call setarr(11,idum,idum)
            call setarr(12,idum,idum)

	 endif

c
c establish the number of points on the gamess/gauss output
c
         call procnv
         call progeo(ipoints,iff,needm)
         if (needm.ge.0) then
             call allgeo(needm,0)
             call progeo(ipoints,iff,istat)
         endif

         if (iftyp.eq.3.and.irtype.eq.4) ipoints = 0
         if (iftyp.eq.7.and.irtype.eq.1) ipoints = 0

c Molecular dynamics .. show trajectory

         if (iftyp.eq.7.and.irtype.eq.3) then 

            call dyncpmd(ipoints)

c find directory 

            idir_now = 0
            finded = .false.
            do while (.not.finded)
               idir_old = idir_now
               idir_now = index(fniun(idir_now+1:),'/')
               if (idir_now.eq.0) finded = .true.
            end do
            filecpmd = fniun(1:idir_old)//'TRAJECTORY'
            indblank = index(filecpmd,' ')-1
            write(*,'(3A)')'Opening file "',filecpmd(1:indblank),
     &           '" to read trajectory data...'
            extest = .false.
            inquire(file=filecpmd,exist=extest)
            if (extest) then
               open(iun5,file=filecpmd,status='old',form='formatted')
            else
               call inferr(
     &              'Can''t find Trajectory file...',1)
            endif
         endif

         if (rdmult) ipoints = 0

         if ((iftyp.eq.2.or.iftyp.eq.4.or.iftyp.eq.5) .and. 
     &        ipoints.le.1 ) then
            inact = inact + 1
c
c For single points show mulliken charges
c
            if (.not.rdmult.and.domul) then
                ipsit = ipsi
                ipsi = 0
                bonds = .false.
                ovrlap = .false.
                atomic = .false.
                call denmak(denok)
                if (denok)
     &              call muldma(vdwr,moddma,.true.,.false.)
                ipsi = ipsit
            endif
         endif

         if ((icrtp.le.3.and.icrtp.ge.1).and.ipoints.le.1) 
     &      inact = inact + 1

         if (iftyp.eq.2.and.ipoints.gt.0.and.ioxyz.eq.0) then

            if (ipoints.gt.maxpt) then
               call inferr('Number of points exceeds maxpt!',1)
               call inferr('Skipping last points           !',0)
               ipoints = min0(ipoints,maxpt)
            endif

            do i=1,ipoints
               idon(i) = 0
            end do

            call rotfir(ioxyz)

         endif

         if (isrf.eq.1) then
            call doscal
            if (iesp.eq.1) call connlp(1.0d0,4,5)
            call docent
         else
            if (zmat) then
               call getpoi(-1,ifd,iff,1,ioatms,ioadd)
            else
c cpmd walks out here
               if (wxyz.or.wmolf) then
                  call getpoi(-3,ifd,iff,1,ioatms,ioadd)
               else
                  call getpoi(-1,ifd,iff,1,ioatms,ioadd)
               endif
            endif
         endif
c         call distchk


         if (ipdbon.eq.1) then
            call ribbs
            if (ioadd.eq.1) then
               if (backb.eq.1) then
                  do l=1,4
                     call acthlp(l,0,0)
                  end do
               else
                  do l=1,4
                     call acthel(1,l-1,0,0)
                  end do
               endif
            endif
         endif

         if ((iftyp.eq.2.or.iftyp.eq.3).and.irtype.eq.4) then
            if (iftyp.eq.3) then
                call ugetfr(istat)
            else
                call ggetfr(istat)
            endif
         endif

c read frequencies in CPMD....
         if (iftyp.eq.7.and.irtype.eq.4) then
            call cpmdgetfr(istat)
            if (istat.ne.0) then
c open MOLVIB file....
c find directory 
               idir_now = 0
               finded = .false.
               do while (.not.finded)
                  idir_old = idir_now
                  idir_now = index(fniun(idir_now+1:),'/')
                  if (idir_now.eq.0) finded = .true.
               end do
               filecpmd = fniun(1:idir_old)//'VIBEIGVEC'
               indblank = index(filecpmd,' ')-1
               write(*,'(3A)')'Opening file "',filecpmd(1:indblank),
     &              '" to read frequencies...'
               extest = .false.
               inquire(file=filecpmd,exist=extest)
               if (extest) then
                  open(iun5,file=filecpmd,status='old',form='formatted')
               else
                  call inferr(
     &                 'Can''t find Vibrational Eigenvectors .',1)
               endif
            endif
         endif

         if (iftyp.eq.4.and.irtype.eq.4.and.ipoints.gt.1) then
              write(iun3,*) ' '
              write(iun3,*)
     &          'WARNING: composite optimisation/frequency job'
              write(iun3,*)
     &          'The frequency job should be a separate file'
              write(iun3,*) 'With First line:'
              write(iun3,*) ' '
              write(iun3,*) ' Entering Gaussian System'
         endif

         if (irtype.eq.4) then
            if (iftyp.eq.1) then
                call getfrm(istat)
            elseif (iftyp.eq.4) then
                call getfr(istat)
            elseif (iftyp.eq.5) then
                call getfra(istat)
                call getint(istat)
            elseif (iftyp.eq.8) then
                call getqfr(istat)
            elseif (iftyp.eq.9) then
                call getofr(istat)
            elseif (iftyp.eq.15) then
                call getnfr(istat)
            endif
         endif

         if (ioadd.eq.0) call setarr(14,idum,idum)

         if (zmat) call wrpnt('zmat.out',8,iopt,0,ipoints,
     &                        fancy,atcol,dolabs,backb)
         if (zmat.and..not.iplot.eq.6) goto 3333

         if (wxyz) call wrpnt('mol.xyz',7,4,0,ipoints,
     &                        fancy,atcol,dolabs,backb)
         if (wxyz.and..not.iplot.eq.6) goto 3333

         if (wmolf) then
             ixyz = 8
             call wrpnt('mol.molf',8,4,0,ipoints,
     &                   fancy,atcol,dolabs,backb)
             if (.not.iplot.eq.6) goto 3333
         endif

         idum1 = ncols
         call xwin(xx,yy,16,str,nstr,idum1,idum2)
         xx=1.0
         call xwin(xx,yy,10,str,nstr,idum1,idum2)
         call xwin(xx,yy,8,str,nstr,idum1,idum2)
         xx= 12.0
         call xwin(xx,yy,99,str,nstr,idum1,idum2)

         fyesno = 0

         if (ioadd.eq.1) then
            namls = namols / 15
            iopval = namols - 15*namls
            call setarr(3,iopval,ioatms)
c            if (ipdbon.eq.1) call actcal(-1)
         else
            if (ipdbon.eq.1) then
               call actcal(0)
            else
               if (iftyp.ne.6.and..not.isrf.eq.1) then
                  call setarr(4,12,ioatms)
               else
                  if (iftyp.eq.6) then
                     call setarr(4,1,ioatms)
                  endif
               endif
            endif
         endif

         if (domol) then
             call xwin(xx,yy,6,str,nstr,idum1,idum2)
             call qupd
         endif
      endif
777   continue
999   if (iplot.eq.6.and.domol) then
          inct = inact
          if (ifav.eq.1) then
             incp = 0
          else
             incp = 1
          endif
          if (irtype.eq.4.or.ihsnmr.eq.2.or.ihasex.eq.1) incp = 2
          if (ipdbon.eq.1) incp = incp + 10
          if (irdogl.eq.1) then
              call prsogl
              irdogl = 0
          endif
          call xwin(xx,yy,15,str,nstr,inct,incp)
          if (inct.eq.15) goto 3333
          if (inct.lt.0.or.inct.eq.290.or.
     &         (inct.ge.420.and.inct.le.450)) goto 1000

          if (inct.eq.130)then
             if (incp.eq.0) then
                 forces = .true.
             else
                 forces = .false.
             endif
          endif

          if (inct.eq.400)then
             if (atcol.eq.1) then
                 atcol = 0
             else
                 atcol = 1
             endif
          endif

          if (inct.eq.410)then
             if (persp.eq.1) then
                 persp = 0
             else
                 persp = 1
             endif
          endif

          if (inct.eq.310) then
             if (backb.eq.0) ssbond = .false.
          endif

          if (inct.eq.490) then
             if (incp.eq.-1) then
                 call acthb(0,hbfilt)
             else
                 keywrd = str(1:nstr)
                 hbfilt = reada(keywrd,1,nstr)/toang
                 call acthlp(1,0,1)
                 call acthlp(2,0,1)
                 call acthb(1,hbfilt)
             endif
          endif

          if (inct.eq.320.and.backb.eq.1) then

             call clkbck(istsurf,incp,ifogl)

          endif

          if (inct.eq.330) then

             ihind = abs(incp)+1

             if (ihind.ge.1.and.ihind.le.4) then

                call acttog(ihind,nstr,0)

             else

                if (istsurf.eq.1) then
                      idosurf = 1
                else
                      idosurf = 0
                endif

                call acttag(incp,nstr,idosurf)

                if (idosurf.eq.1) then
                   if (ifogl.eq.1) then
                      call initsrf
                   else
                      call connlp(1.0d0,0,4)
                   endif
                endif

             endif
          endif

          if (inct.eq.340) then
             if (ssbond) then
                ssbond = .false.
                call actss(0)
             else
                ssbond = .true.
                call actss(1)
             endif
          endif

          if (inct.eq.345.and.backb.eq.1) then
             call inferr
     &       ('Click on the backbone to activate !',0)
          endif

          if (inct.eq.460.and.backb.eq.1) 
     &       call aacom(vrad(100),incp,str,nstr,istsurf)

          if (inct.eq.461) call proxic(incp,backb,0,1)

          if (inct.eq.470) then
             if (incp.ne.-1) then
                keywrd = str(1:nstr)
                call wrpnt(keywrd,nstr,iwropt,0,ipoints,
     &                     fancy,atcol,dolabs,backb)
             endif
          endif

          if (inct.eq.480.or.inct.eq.481) then

             keywrd = str(1:nstr)

             if (inct.eq.480.and.(nstr.eq.0.or.keywrd.eq.' ')) then
                call inferr('Invalid Filename !',0)
             else

                if (idebug.eq.1) write(iun3,*) 'reading file: ',keywrd
                if (inct.eq.480) then
                   if (odos(keywrd)) call tounx
                endif

                fniun = keywrd

                iof = 1
                if (inct.eq.480) then
                   nxtmf = 0
                else
                   nxtmf = 1
                   iof = 0
                endif

                call newfil(idebug,istat,incp,ioadd,ioatms,nstrt,namols,
     &                      nxtmf,ipdbon,namls,iof)
                if (istat.eq.0) goto 10

             endif

          endif

          if (inct.eq.140) then
             if ( incp.eq.1 ) then
                 fancy = 0
             else
                 fancy = 1
             endif
          endif

          if (inct.eq.150)then
              keywrd = str(1:nstr)
              if (idebug.eq.1) write(iun3,*) 'Opening file: ',keywrd
              if (opfil(iun4,keywrd,nstr,1,0,0)) then
                 if (iatoms.gt.50) then
                      call plini(4,.false.,icolps)
                 else
                      call plini(4,.true.,icolps)
                 endif
                 call plpost(backb,dolabs,icolps,fancy,atcol,persp,
     &                       shade,idelx)
                 if (ipoints.eq.0) then
                    keywrd = ' '
                    if (iftyp.ne.1) keywrd = 'single point'
                 elseif (ipnt.eq.1) then
                    keywrd = 'first point'
                 elseif (ipnt.eq.ipoints) then
                    keywrd = 'last point'
                 else
                    keywrd = 'point '//gstr(ipnt)
                 endif
                 if (ipdbon.eq.1) then
                    call plpend
                 else
                    call plend(title,notitle)
                 endif
                 close(iun4)
                 call inferr('PostScript Ready !',0)
                 call qupd
              else
                 call inferr('Error Opening File !',0)
              endif
          endif

          if (inct.eq.151) then
              if (iixyz.ne.10.and.iixyz.ne.11) then
                 if (chkmap(idum)) then
                    dozme = .true.
                    call getpoi(-1,ifd,iff,1,ioatms,ioadd)
                 endif
              endif
              if (iixyz.ge.4) then
                  call runjob(iixyz-3,iqopt,ihaszm)
              else
                 keywrd = str(1:nstr)
                 iform = 1
                 if (iixyz.eq.2) iform = 0
                 if (opfil(46,keywrd,nstr,iform,1,0)) then
                    if (mapxyz(46,iixyz,0,1)) then
                       dozme = .true.
                       ipdbon = 0
c check
                       iftyp = 0
                       movie = 0
                       call getpoi(-1,ifd,iff,1,ioatms,ioadd)
                    endif
                    close(46)
                 endif
              endif
          endif
          if (inct.eq.152.or.inct.eq.1520)then
              keywrd = str(1:nstr)
              if (opfil(46,keywrd,nstr,1,1,0)) then
                  iuntmp = iun2
                  iun2 = 46
                  if (inct.eq.152) then
                     call aln2ml(0,istat)
                  else
                     call aln2ml(1,istat)
                     call qupd
                  endif
                  close(46)
                  iun2 = iuntmp
                  if (istat.eq.1) then
                     irtcel = 2
                     irtval = 2
                     atcol = 0
                     call butset(1,21,0)
                  endif
              endif
          endif
          if (inct.eq.153) then
              if (movie.eq.1) then
                 movie = 0
              else
                 if (irtcel.eq.0) then
                    irtcel = irtval
                    if (ipdbon.eq.1.and.iff.ne.7.and.ialtyp.ne.0) then
                       call pmfass(0,0)
                       call totpmf(dum)
                       call upsco()
                    endif
                    call mtinv3
                    call rarbxi
                    call rarbyi
                    call rarbzi
                 else
                    irtcel = 0
                 endif
              endif
          endif
          if (inct.eq.154) then
              call fdat(ifd,0,istdbd,iuseab,moddma,idebug)
          endif
          if (inct.eq.540) then
              iuntmp = iun2
              iun2 = 48
              nztmp = nz
              keywrd = str(1:nstr)
              if (idebug.eq.1) write(iun3,*) 'Opening file: ',keywrd
              if (opfil(iun2,keywrd,nstr,1,1,0)) then
                 rewind iun2
                 call getzm(iatoms,1,0,istat)
                 if (istat.eq.1) then
                     close(iun2)
                 else
                     rewind iun2
                     call getmop(iatoms,heat,1,0,istat) 
                     if (ista.eq.1) close(iun2)
                 endif
                 call cnvfrg(nztmp)
              else
                 call zmterr('error opening fragment file !',1,0,1)
              endif
              iun2 = iuntmp
          endif
          if (inct.eq.550)then
              if (incp.eq.14) then
                 if (irtcel.eq.1) then
                    irtcel = 0
                 else
                    irtval = 1
                    irtcel = 1
                    itz = 0
                    call mtinv3
                    call rarbxi
                    call rarbyi
                    call rarbzi
                 endif
              elseif (incp.ne.15) then
                 if (incp.eq.22) then
                   if (opfil(46,'shelx.ins',9,1,0,0)) then
                      call wrshlx(46,0)
                      close(46)
                   endif
                 elseif (incp.eq.21) then
                   if (opfil(46,'crystal95.in',12,1,0,0)) then
                      call wrcrys(46)
                      close(46)
                   endif
                 elseif (incp.eq.23) then
                   if (opfil(46,'POSCAR',6,1,0,0)) then
                      call fdat(4,1,istdbd,iuseab,moddma,idebug)
                      call wrvasp(46)
                      call fdat(ifd,0,istdbd,iuseab,moddma,idebug)
                      close(46)
                   endif
                 elseif (incp.eq.24) then
                   if (opfil(46,'pluton.spf',9,1,0,0)) then
                      call wrshlx(46,1)
                      close(46)
                   endif
                 elseif (incp.eq.25) then
                   if (opfil(46,'mopac.dat',9,1,0,0)) then
                      call fdat(4,1,istdbd,iuseab,moddma,idebug)
                      call wrmopa(46)
                      call fdat(ifd,0,istdbd,iuseab,moddma,idebug)
                      close(46)
                   endif
                 elseif (incp.eq.26) then
                   if (opfil(46,'cpmd.inp',8,1,0,0)) then
                      call wrcpmd(46,.true.)
                      close(46)
                   endif
                 elseif (incp.eq.27) then
                   if (opfil(46,'mol.cif',7,1,0,0)) then
                      call wrcif(46)
                      close(46)
                   endif

c                =========================================
                 elseif (incp.eq.28) then
                   if (opfil(46,'conquest.coord',6,1,0,0)) then
                      call fdat(4,1,istdbd,iuseab,moddma,idebug)
                      call wrconquest(46)
                      call fdat(ifd,0,istdbd,iuseab,moddma,idebug)
                      close(46)
                   endif
c                =========================================

                 else
                   ifdt = ifd
                   ifd = incp
                   if (ifdt.ne.ifd) iscupd = 1
                   if (ifd.eq.20) 
     &                call fdat(15,0,istdbd,iuseab,moddma,idebug)
                   if (ifd.eq.19) then
                      tnknm = 'xxx'
                      tnkprg = 4
                      tnkbg = 1
                      tnkarc = 0
                      rmsgrd = 0.01d0
                      ifd = ifdt
                      call runjob(7,iqopt,ihaszm)
                   endif
                   if (ifd.eq.17) then
                      call fdat(1,0,istdbd,iuseab,moddma,idebug)
                      call doconn
                      call dohcon(0)
                      call setarr(5,idum,ioatms)
                      ifd = 2
                   endif
                   call fdat(ifd,0,istdbd,iuseab,moddma,idebug)
                 endif
              endif
          endif
          if (inct.eq.560) call setorg(incp)
          if (inct.eq.561) then
             call pmfinf(incp)
          endif
          if (inct.eq.562) then
             call pckrot(istat,incp)
             if (istat.eq.1) then
                 irtcel = 2
                 irtval = 2
                 atcol = 0
                 call butset(1,21,0)
             endif
          endif
          if (inct.eq.565) then
             if (ipdbon.ne.0) then
                  call setarr(1,incp,ioatms)
             else
                  call inferr('Not a protein !',0)
             endif
          endif

          if (inct.eq.576) then
             iesp = 0
             isfdon = 0
                 if (natorg.ne.0) iatoms = natorg
c                 call getpoi(-1,ifd,iff,1,ioatms,ioadd)
          endif
          if (inct.ge.620.and.inct.le.622) call disabh(inct-620)
          if (inct.eq.623) then
              call nohcon
              call dohcon(0)
          endif
          if (inct.eq.624) call nohcon
          if (inct.eq.626) then
              htmp = hdmax
              keywrd = str(1:nstr)
              hdmax = reada(keywrd,1,nstr)/toang
              call nohcon
              call dohcon(0)
              hdmax = htmp
          endif
          if (inct.eq.630) then
              if (isfdon.ne.0.and..not.iesp.eq.1) then
                 iatoms = natorg
                 natorg = 0
                 isfdon = 0
              endif
          endif
          if (inct.ge.570.and.inct.le.573) then
              if (( ((iftyp.ge.2.and.iftyp.le.4).or.
     &              (iftyp.eq.5.and.isgau.eq.1))
     &            .and.norbs.ne.0.and.inct.le.571)
     &            .or.((inct.eq.572.or.inct.eq.573).and.
     &            (ihasq.eq.1.or.(ihasq.eq.0.and.ifogl.eq.1)))) then

                  iesp = 1
                  ifreq = 0
                  if (inct.eq.572.or.inct.eq.573) then
                      if (natorg.eq.0) then
                          natorg = iatoms
                      else 
                          iatoms = natorg
                      endif
                  else
                      call xyzcoo(1,0,0)
                      natorg = natoms
                  endif
                  call haszm(.false.)
                  call doconn

                  inctt = inct-569
                  if (inct.eq.573) inctt = inctt - 1

                  if (inct-569.eq.isfdon.and.istyp.eq.isftyp.and.
     &               (istyp.eq.0.or.(istyp.eq.1.and.
     &                edge.eq.edget.and.ctvalt.eq.ctval(1)
     &                .and.nspts.eq.nsptst))) then
                     iatoms = nesp
                     call doscal
                  else
                     call curs(1)
                     if (inct.eq.570.or.inct.eq.571) then
                        ipsi = 0
                        dolap = .false.
                        bonds = .false.
                        ovrlap = .false.
                        atomic = .false.
                        call denmak(denok)
                     endif
                     if (inct.eq.571) then 
                        call muldma(vdwr,moddma,.true.,.true.)
                     endif
                     if ((inct.eq.570.or.inct.eq.571)
     &                  .and.istyp.eq.1) then
                        if (isofst.eq.1) isofst = 0
                        iatoms = natorg
                        call isoden(ctval,1,1.0d0,idum,0)
                        call connlp(1.0d0,1,5)
                        edget = edge
                        ctvalt = ctval(1)
                        nsptst = nspts
                     else
                        if (inct.eq.573) ipsa = 1
                        call allsrf(0,inctt,1,4)
                        ipsa = 0
                     endif
                     call doscal
                     call curs(0)
                     call inferr('Electrostatic Potential Done !',0)
                     space = 0
                     isfdon = inct - 569
                     isftyp = istyp
                     nesp = iatoms
                  endif
                  if (fancy.eq.1) then
                     atcol = 1
                     call butset(1,21,1)
                  else
                     atcol = 0
                     call butset(1,21,0)
                  endif
              else
                  call inferr('No Electrostatic Potential !',0)
              endif
          endif

          if (inct.eq.574) then
              keywrd = str(1:nstr)
              if (opfil(46,keywrd,nstr,1,0,0)) then
                 if (nesp.eq.0) then
                    call wrsrf(46,iatoms,0)
                 else
                    call wrsrf(46,nesp,iesp)
                 endif
                 close(46)
              endif
          endif

          if (inct.eq.575) then
              keywrd = str(1:nstr)
              if (opfil(46,keywrd,nstr,1,1,0)) then
                 call rdsrf(46,istat,iesp,0,idebug)
                 call doscal
                 if (istat.eq.1) then
                    nesp = iatoms
                    if (iesp.eq.1) call connlp(1.0d0,4,5)
                    if (fancy.eq.1) then
                       atcol = 1
                       call butset(1,21,1)
                    else
                       atcol = 0
                       call butset(1,21,0)
                    endif
                 endif
                 close(46)
              endif
          endif

          if (inct.ge.577.and.inct.le.579) then
              ifreq = 0
              inctt = inct - 577
              call calelc(inctt,isfdon,dolabs,space,moddma)
          endif

          if (inct.eq.580) then
              if (ihasesp.eq.1) then
                 ihasq = 1
                 jatoms = natoms
                 if (ioni.eq.1) jatoms = iatoms
                 call setarr(7,jatoms,ioatms)
              else
                 call inferr('NO esp charges available from output !',0)
              endif
          endif
c          if (inct.eq.581) then
c              ieemop = incp
c              call eem(ieemop,1,istat)
c          endif
c          if (inct.eq.582) then
c              call calgas
c          endif
          if (inct.eq.584) then
c omap file
              ostep = .false.
              scale = -1.0d0
              cut   = 1.0d0
              keywrd = str(1:nstr)
              if (.not.opfil(21,keywrd(1:nstr),nstr,1,1,0)) then
                 istat = 0
                 call inferr('Error Opening File !',0)
              else
                 istat = 0
                 ijag = 3
                 close(21)
                 call cubtst(21,ijag)
                 if (ijag.ne.-1) then
                    call rdomap(npts1,npts2,npts3,21,istat)
                    close(21)
                    idofil = 1
                    ipsi = 1
                    valcnt = 0.01d0
                    space = 1
                    if (keyr(keywrd,'SPACE',tmpr)) then
                        valcnt = tmpr
                    endif
                 endif
              endif
           endif

           
          if (inct.eq.586) then
c kont file
              ostep = .false.
              scale = -1.0d0
              cut   = 1.0d0
              keywrd = str(1:nstr)
              if (.not.opfil(21,keywrd(1:nstr),nstr,1,1,0)) then
                 istat = 0
                 call inferr('Error Opening File !',0)
              else
                 istat = 0
                 ijag = 2
                 call cubtst(21,ijag)
                 if (ijag.ne.-1) then
                    call rdgrd(npts1,npts2,npts3,21,istat)
                    close(21)
                    idofil = 1
                    ipsi = 1
                    valcnt = 0.01d0
                    space = 1
                    if (keyr(keywrd,'SPACE',tmpr)) then
                        valcnt = tmpr
                    endif
                 endif
              endif
          endif
          if (inct.eq.588) then
              if (((iftyp.ge.2.and.iftyp.le.4).or.
     &             (iftyp.eq.5.and.isgau.eq.1)
     &            ).and.norbs.ne.0) then
                  call xyzcoo(1,0,0)
                  call haszm(.false.)
                  call doconn
                  call curs(1)
                  ipsi = 0
                  ifreq = 0
                  dolap = .false.
                  bonds = .false.
                  ovrlap = .false.
                  atomic = .false.
                  dmachg = .false.
                  call denmak(denok)
                  if (incp.eq.1) then
                     dmachg = .true.
                     call muldma(vdwr,moddma,.true.,.true.)
                  endif
                  keywrd = 'MPFIT '
                  if (opfil(46,'keywrd',6,1,1,0)) then
                     read(46,'(a)') line
                     keywrd = keywrd(1:6)//line(1:linlen(line))
                     call tocap(keywrd,320)
                     close(46)
                  endif
                  call espchrg(ctval,nvalc,.false.,dmachg)
                  call curs(0)
              else
                  call inferr('Can not calculate El. Pot. !',0)
              endif
          endif
          if (inct.eq.589) then
              keywrd = str(1:nstr)
              iun = 46
              if (opfil(iun,keywrd,nstr,1,0,0)) then
                 call plvrml(iun,fancy,atcol,dolabs,1,backb,.true.)
                 call inferr('Wrote VRML file',0)
                 close(iun)
              endif
          endif
          if (inct.eq.640) then
              keywrd = str(1:nstr)
              iun = 46
              if (opfil(iun,keywrd,nstr,1,1,0)) then
                 iuntmp = iun2
                 iun2 = iun
                 call rdprot
                 close(iun)
                 iun2 = iuntmp
              else
                 call inferr('Can not open Protein file',0)
              endif
          endif
          if (inct.eq.590.or.inct.eq.591) then
              if (inct.eq.590) then
                 keywrd = str(1:nstr)
                 spacng = reada(keywrd,1,nstr)
                 if (spacng.eq.0.0d0) spacng = 0.05d0
                 call setcod(spacng)
                 if (itsrf.eq.1) then
                    call clrsrf
                    call allsrf(0,isfdon,0,4)
                 endif
              endif
              call clrcod(natorg,iatoms,1)
          endif
          if (inct.eq.610) then
              if (opfil(46,'esp.xyz',7,1,1,0)) then
                  iuntmp = iun2
                  iun2 = 46
                  call getchr(igetca)
                  if (igetca.eq.1) then
                     call inferr('Read esp.xyz',0)
                  else
                     call inferr('Error reading esp.xyz',0)
                  endif
                  close(46)
                  iun2 = iuntmp
              else
                  call inferr('Error reading esp.xyz',0)
              endif
          endif
          if (inct.eq.170) then
             if (incp.eq.1) then
                call inferr('Click on two atoms !',0)
             elseif (incp.eq.2) then
                call inferr('Click on three atoms !',0)
             elseif (incp.eq.3) then
                call inferr('Click on four atoms !',0)
             endif
          endif
          if (inct.eq.180) then
             if ( dolabs.ge.1 ) then
                 dolabs = 0
             else
                 dolabs = 1
                 if (ihasq.eq.1.and.iqon.eq.3) then
                    if (((iftyp.ge.2.and.iftyp.le.4).or.
     &                   (iftyp.eq.5.and.isgau.eq.1)
     &                  ).and.norbs.ne.0.and.ioni.ne.1) then
                        call xyzcoo(1,0,0)
                        call doconn
                    endif
                 endif
             endif
          endif
          if (inct.eq.500) then
             icoars = incp
             if (icoars.eq.2) then
               call inferr('45 degree rotation',0)
             elseif (icoars.eq.1) then
               call inferr('5 degree rotation',0)
             else
               call inferr('1 degree rotation',0)
             endif
          endif
510       continue
          if (inct.eq.510) then
             call setnrm(incp,iwopt,wmovie,dconpt)
          endif

          if (inct.eq.520) then
             if (ihasex.eq.0) call resfr
             call doconn
             normc = 0
          endif
          if (inct.eq.530) then
             ifreq = 0
             dozme = .true.
             movie = 0
             call getpoi(-1,ifd,iff,1,ioatms,ioadd)
          endif
          if (inct.eq.700) then
             nxtmf = 1
             goto 10
          endif
          if (inct.eq.100) then
             iesp = 0
             ifreq = 0
             movie = 0
             call getpoi(-1,ifd,iff,0,ioatms,ioadd)
          endif
          if (inct.eq.110) then
             movie = 0
             iesp = 0
             isfdon = 0
             ifreq = 0
             call getpoi(-2,ifd,iff,0,ioatms,ioadd)
          endif
          if (inct.eq.109) then
             movie = 0
             iesp = 0
             isfdon = 0
             ifreq = 0
             call getpoi(-4,ifd,iff,0,ioatms,ioadd)
          endif
          if (inct.eq.155) then
             if (movie.eq.1) then
                movie = 0
                call parptr(16,fdum,fdum,idum)
             else
                movie = 1
                iesp = 0
                isfdon = 0
                ifreq = 0
                inct = 0
                call parptr(129,fdum,fdum,idum)
             endif
          endif
          if (inct.eq.115) then
             movie = 1
             iesp = 0
             isfdon = 0
             ifreq = 0
             call getpoi(-1,ifd,iff,0,ioatms,ioadd)
             if (wmovie) call wrpnt(movfil,len(movfil),iwopt,1,ipoints,
     &                              fancy,atcol,dolabs,backb)
          endif
          if (inct.eq.116) then
             iwopt = 6
             if (normc.eq.1) then
                wmovie = .true.
                inct = 510
                goto 510
             endif
             iesp = 0
             isfdon = 0
             ifreq = 0
             if (ipoints.eq.0.or.ipoints.eq.1) then
                if (ivtwo.eq.2) then
                   call wrpnt(povfil,len(povfil),iwopt,0,ipoints,
     &                     fancy,atcol,dolabs,backb)
                   call inferr('wrote file: '//povfil,0)
                else
                   call wrpnt(oglfil,len(oglfil),iwopt,0,ipoints,
     &                     fancy,atcol,dolabs,backb)
                   call inferr('wrote file: '//oglfil,0)
                endif
                goto 116
             else
                call getpoi(-1,ifd,iff,0,ioatms,ioadd)
                call wrpnt(oglfil,len(oglfil),iwopt,1,ipoints,
     &                     fancy,atcol,dolabs,backb)
             endif
             do while (.true.)
                 call getpoi(-2,ifd,iff,0,ioatms,ioadd)
                 if (ipnt.lt.ipoints) then
                    call wrpnt(oglfil,len(oglfil),iwopt,3,ipoints,
     &                         fancy,atcol,dolabs,backb)
                 else
                    call wrpnt(oglfil,len(oglfil),iwopt,2,ipoints,
     &                         fancy,atcol,dolabs,backb)
                    call inferr('wrote file: '//oglfil,0)
                    goto 116
                 endif
             end do
116          continue
             if (ivtwo.eq.3) call viewer
          endif

          if (inct.eq.117) then
             movie = 0
             iesp = 0
             isfdon = 0
             if (incp.gt.0.and.incp.le.ipoints)
     &           call getpoi(incp,ifd,iff,1,ioatms,ioadd)
          endif
c          if (inct.eq.120) then
c             call getpoi(-3,ifd,iff,1,ioatms,ioadd)
c          endif
1000      continue
          fyesno = 0
          if ( forces .and. ifav.eq.1 ) fyesno = 1
          call xwin(xx,yy,6,str,nstr,idum1,idum2)
          if (inct.lt.0.or.inct.eq.290.or.
     &         (inct.ge.420.and.inct.le.450)) call mktrn(inct,incp)
          call qupd
          if (movie.eq.1.and.ipoints.gt.1) then
             call dlystr
             if (ipnt.eq.ipoints) movie = 0
             call getpoi(-2,ifd,iff,0,ioatms,ioadd)
             if (wmovie) then
                 if (ipnt.lt.ipoints) then
                    call wrpnt(movfil,len(movfil),iwopt,3,ipoints,
     &                         fancy,atcol,dolabs,backb) 
                 else
                    call wrpnt(movfil,len(movfil),iwopt,2,ipoints,
     &                         fancy,atcol,dolabs,backb)
                    call inferr('wrote file: '//movfil,0)
                 endif
             endif
             if (ipnt.lt.ipoints) goto 999
             call parptr(16,fdum,fdum,idum)
             movie = 0
             goto 1000
          endif
          if (irtype.eq.4.and.normc.eq.1) then
             call nxtpnt(iwopt,fancy,atcol,dolabs,backb,wmovie)
             if (dconpt) call doconn
             call dlystr
          endif
          if (inct.ne.75) goto 999

          if (irtype.eq.4.and.normc.eq.1) then
             call resfr
             call doconn
             normc = 0
          endif

          iplwin = 1
          if (molonl) stop

          domol = .false.
          if (euclid) then
             call parang(0.0d0,0.0d0)
          else
             call parang(theta,phi)
          endif
          xx=1.0
          call xwin(xx,yy,10,str,nstr,idum1,idum2)
          call xwin(xx,yy,8,str,nstr,idum1,idum2)
          call xwin(xx,yy,6,str,nstr,idum1,idum2)
      endif

      if (iplot.eq.6) then
        if (docnt)then
           xx= 15.0
        else
           xx= 2.0
        endif
        yy = 0.0
        call xwin(xx,yy,99,str,nstr,idum1,idum2)
        call xwin(xx,yy,8,str,nstr,idum1,idum2)
        if (ipsprt.eq.0) call xwin(xx,yy,6,str,nstr,idum1,idum2)
        ifreq = 0
      endif

      if (first) then
c
c     distribute electrons over mos
c
         if (dovdw) then
             if (idvrml.eq.1) then
                call epvrml(vdwr,moddma,natoms,norbs,0)
             else
                call epvrml(vdwr,moddma,natoms,norbs,0)
             endif
         endif

         if (molpot) then
             ipsi = 0
             dolap = .false.
         endif

         if (.not.rdmult) call denmak(denok)
         if (molpot.and..not.rdmult) then
             call muldma(vdwr,moddma,.true.,.true.)
         endif
c
c
c     calculate function values in the plot plane
c
         call grdcl(npts1,npts2,1,space)

         if (iplot.eq.5) then
            if (space.eq.1) then
             spcdfil = vrmlfil
             ivtwo = 3
             idvrml = 1
             call spaced(npts1,npts2,npts3,valcnt,idofil,
     &                adjus,ipsprt,idisml,idvrml,mapit)
            else
               call oginit(r,adjus,natoms,nat,icol,
     &                   xsym,ysym,zsym,vdwr,-scale,npts1,npts2,ipsi)
            endif
         endif

         if (.not.rdmult) then
            if (index(keywrd,'WRBAS').ne.0.and.space.eq.0)
     &      call wrinfo(npts1,npts2)
            if (index(keywrd,'RDBAS').ne.0.and.space.eq.0)
     &      call rdinfo(npts1,npts2,0,istat)

c
c     print density at atom centres
c
            call atmd(iftyp.eq.1,ipsi,iftyp.eq.5)

         endif
         first = .false.
      endif

      if (elpot.and.idebug.eq.1) then
         write(iun3,*) ' '
         write(iun3,*) 'electrostatic potential at the atoms'
         write(iun3,*) '(nuclear repulsion excluded)'
         write(iun3,*) ' '
         do l=1,natoms
             call espot(xyz(1,l),xyz(2,l),xyz(3,l),pot,0)
             write(iun3,'(4x,a2,3f10.5,2x,f13.5)')
     &            elemnt(nat(l)),(xyz(ig,l),ig=1,3),pot
         end do
         write(iun3,*) ' '
      endif

      if (space.eq.1.and.(natoms.gt.0.or.idocub.eq.1.or.
     &    (icubtp.eq.2.or.icubtp.eq.3))) then
          iesp = 0
          isfdon = 0
          if (idocub.eq.1) then

c istat = 0; Gaussian cube

             call cubtst(iun2,0)
             istat = 0
             call rdcube(npts1,npts2,npts3,iposng,ipsi,istat,
     &                   iun2,idebug)
             if (mapit.eq.1) then
                call rd3chk(npts1,npts2,npts3,mapopt,0,istat)
                if (istat.eq.1) then
                   call rd3chk(npts1,npts2,npts3,mapopt,1,istat)
                else
                   idvrml = 0
                endif
             endif
             if (idofil.eq.0) idofil = 1

             if (istat.eq.1) then
                if (iposng.eq.1) ipsi = 1
                adjus = toang
                idisml = 1
                valcnt = 0.01d0
                space = 1
                ispd = 1
                if (keyr(keywrd,'SPACE',tmpr)) then
                    valcnt = tmpr
                endif
             endif
          endif
          if (ivtwo.eq.-1) then
             if (ifdogl.eq.1) then
                ivtwo = 4
                idvrml = 1
             else
                idvrml = 0
             endif
          endif
          if (ivtwo.eq.4) idvrml = 1

          if (ivtwo.eq.2) then
             spcdfil = povfil
             call spaced(npts1,npts2,npts3,valcnt,idofil,
     &                adjus,ipsprt,idisml,idvrml,mapit)
          else
             spcdfil = vrmlfil
c             if (ofrst.and.elpot) then
c                 call denmak(denok)
c                 call espgrd(npts1,npts2,npts3,idebug)
c                 ofrst = .false.
c             endif
             call spaced(npts1,npts2,npts3,valcnt,idofil,
     &                adjus,ipsprt,idisml,idvrml,mapit)
          endif

          if (dosrf2) then
             if (adjus.eq.1.0d0) then
                 call xyzcoo(1,1,0)
             else
                 call xyzcoo(1,0,0)
             endif
             call doconn
             if (adjus.eq.1.0d0) call xyzcoo(1,0,0)
             natorg = natoms
             iatoms = natorg
             call setarr(2,idum,ioatms)
             call spasrf(npts1,npts2,npts3,valcnt)
             if (adjus.eq.1.0d0) call xyzcoo(0,1,0)
             call doscal
             call docent
             dosrf2 = .false.
          endif
          
          if (index(keywrd,'WRBAS').ne.0)
     &      call wr3inf(npts1,npts2,npts3,adjus)

          if (index(keywrd,'RDBAS').ne.0) then
            call tstrd3
            call rd3inf(npts1,npts2,npts3,1,adjus,istat)
            ispd = 1
            spcdfil = vrmlfil
            call spaced(npts1,npts2,npts3,valcnt,idofil,
     &                adjus,ipsprt,idisml,idvrml,mapit)
          endif

      else

         if (ifdogl.eq.1.and.(do3dx.or.do3d.eq.1)
     &      .and.inct.ne.651) then
            call bldlst
            goto 1699
         endif

         if (do3dx) then
            call plden(npts1,npts2,scale,icells,adjus,idisml)
            goto 1699
         endif

         if (do3d.eq.1) call den3d(npts1,npts2,scale)

         if (docnt.and.pmax.ne.0.0d0.and..not.do3dx) 
     &      call dencnt(npts1,npts2,fcnt)

      endif

      call selsol

      if (euclid.and.space.eq.0) call plbox

      if (iplot.eq.4) then
        write(iun4,'(''s'')')
        write(iun4,'(''   3 setlinewidth'')')
        write(iun4,'(''n'')')
      endif
      if (iplot.eq.3) then
         idum=1
         call plotgh(idum,0.0d0,0.0d0)
      endif
      if (iplot.eq.6) then
         idum=1
         call plotgh(idum,0.0d0,0.0d0)
      endif
 
      if ((iplot.eq.6.or.iplot.eq.4).and.
     &(.not.euclid).and.idisml.eq.1.and.space.eq.0) 
     & call pl3dm(adjus,.false.,rdum1)

      if (euclid.and.space.eq.0) call eucmol(vdwr,adjus)
      if (domax.and..not.do3dx.and.space.eq.0) 
     &    call maxmin(npts1,npts2,scale)

      call plhead(title)

      if (docnt.and.space.eq.0) call pltab(fcnt)


      if (ipsprt.eq.1) then
        ipsprt = 0
        keywrd = ' '
        if (ipsi.ne.0) keywrd = ' Psi = '//gstr(ipsi)//
     &                          keywrd
        if (molpot) keywrd = ' Molpot'//keywrd
        if (euclid) keywrd = ' Euclid'//keywrd
        if (do3d.eq.1  ) keywrd = ' 3D'//keywrd
        if (docnt ) keywrd = ' Contour'//keywrd
        if (bonds ) keywrd = ' Bonds'//keywrd
        if (doori ) keywrd = ' Orient'//keywrd
        if (ovrlap) keywrd = ' Overlap'//keywrd
        if (atomic) keywrd = ' Atomic'//keywrd
        if (space.eq.1) then
            tmpstr = ' Space ='
            write(tmpstr(9:16),'(f8.4)') valcnt
            keywrd = tmpstr(1:16)//keywrd
        endif
        tmpstr = ' Edge ='
        write(tmpstr(8:13),'(f6.2)') r(1)
        keywrd = tmpstr(1:13)//keywrd
        call plend(title,notitle)
        close(iun4)
        call inferr('PostScript Ready !',0)
        iplot = 6
      endif

1699  if (iplot.eq.6) then
          call qupd
1700      inct = 0
          if (iftyp.eq.1) inct = 1
          call xwin(xx,yy,0,str,nstr,inct,incp)
          iplwin = 1
          if (inct.ne.15) then
             if (inct.eq.30) then
                if (incp.eq.1) scale = scale*sincr
                if (incp.eq.-1) scale = scale/sincr
                scale = min(scale,-0.01d0)
                scale = max(scale,-100.0d0)
             elseif (inct.eq.45) then
                if (idisml.eq.1) then
                   idisml = 0
                else
                   idisml = 1
                endif
             elseif (inct.eq.80) then
                if (docnt) then
                   docnt = .false.
                else
                   docnt = .true.
                endif
             elseif (inct.eq.85) then
                euclid = .true.
                docnt = .true.
                do3d  = 0
                do3dx = .false.
                space = 0
                rdum1 = 0.0d0
                rdum2 = 0.0d0
                call parang(rdum1,rdum2)
             elseif (inct.eq.90) then
                do3d  = 1
                do3dx = .false.
                docnt = .true.
                euclid = .false.
                space = 0
                theta = 60.0d0
                phi = 45.0d0
                call parang(theta,phi)
             elseif (inct.eq.91) then
                do3d  = 0
                do3dx = .true.
                docnt = .false.
                euclid = .false.
                space = 0
                theta = 60.0d0
                phi = 45.0d0
                call parang(theta,phi)
             elseif (inct.eq.92.or.inct.eq.600.or.inct.eq.94) then
                do3d  = 0
                do3dx = .false.
                docnt = .false.
                euclid = .false.
                space = 1
                tx = 0.5d0
                ty = 0.5d0*0.33333d0
                tz = 0.5d0
                call setang(tx,ty,tz,theta,phi)
                if (inct.eq.94) then
                   dosrf2 = .true.
                else
                   keywrd = str(1:nstr)
                   valcnt = reada(keywrd,1,nstr)
                endif
                if (inct.eq.600) then
                   idvrml = 1
                   ipsa = 1
                endif
                if (inct.eq.92.and.ifdogl.ne.0) then
                   idvrml = 1
                   ivtwo = 4
                endif
                if (mapit.eq.1) then
                   mapfil = mapfil(1:linlen(mapfil))
                   if (ivtwo.eq.3) then
                      oglfil = oglfil(1:linlen(oglfil))
                      call parogf(oglfil,linlen(oglfil))
                   elseif (ivtwo.eq.1) then
                      vrmlfil = vrmlfil(1:linlen(vrmlfil))
                   endif
                   call rd3chk(npts1,npts2,npts3,mapopt,0,istat)
                   if (istat.eq.1) then
                      call rd3chk(npts1,npts2,npts3,mapopt,1,istat)
                   else
                      idvrml = 0
                   endif
                else
                   if (idvrml.eq.1) then
                      if (ivtwo.eq.0.or.ivtwo.eq.1) then
                         vrmlfil = vrmlfil(1:linlen(vrmlfil))
                      elseif (ivtwo.eq.2) then
                         povfil = povfil(1:linlen(povfil))
                      elseif (ivtwo.eq.3) then
                         oglfil = oglfil(1:linlen(oglfil))
                         call parogf(oglfil,linlen(oglfil))
                      endif
                   endif
                endif
             elseif (inct.eq.93) then
                if (idofil.eq.1) then
                   idofil = 0
                else
                   if (icubtp.eq.0) idofil = 1
                endif
             elseif (inct.eq.650) then
                keywrd = str(1:nstr)
                if (opfil(46,keywrd,nstr,1,0,0)) then
                    ivtwo = 1
                    call p3dv(46,scale,npts1,npts2,adjus)
                    close(46)
                endif
             elseif (inct.eq.651) then
                    ivtwo = 3
                    call oginit(r,adjus,natoms,nat,icol,
     &                   xsym,ysym,zsym,vdwr,-scale,npts1,npts2,ipsi)
             elseif (inct.eq.150) then
                keywrd = str(1:nstr)
                if (idebug.eq.1) write(iun3,*) 'reading file: ',keywrd
                if (opfil(iun4,keywrd,nstr,1,0,0)) then
                    call plini(4,.true.,icolps)
                    ipsprt = 1
                    iplot = 4
                else
                    call inferr('Error Opening File !',0)
                endif
             elseif (inct.eq.160) then
                if ((incp.gt.0.and.incp.le.ncols).or.
     &              (incp.lt.0.and.abs(incp).le.ncolb)) then
                   ostep = .false.
                   oscal = .false.
                   ofrst = .true.
                   cut   = 1.0d0
                   ipsi = incp
                   bonds = .false.
                   dolap = .false.
                   elpot = .false.
                   molpot = .false.
                   chpot = .false.
                   call denmak(denok)
                   call grdcl(npts1,npts2,1,space)
                   if (pmax.eq.0.0d0) then
                      call inferr('Nodal Plane ?!',0)
                   else
                      lstr = 'Orbital '//hstr(abs(incp))
                      call inferr(lstr,0)
                   endif
                else
                   if (incp.ne.0) then
                      call inferr('Invalid Orbital Number !',0)
                   endif
                endif
             elseif (inct.eq.170) then
                ostep = .false.
                oscal = .false.
                ofrst = .true.
                cut   = 1.0d0
                ipsi = 0
                bonds  = .false.
                atomic = .false.
                ovrlap = .false.
                elpot = .false.
                molpot = .false.
                chpot = .false.
                dolap = .false.
                call denmak(denok)
                if (denok) then
                   call grdcl(npts1,npts2,1,space)
                   if (ispd.eq.0) then
                      call inferr('Normal Density',0)
                   else
                      call inferr('Spin Density',0)
                   endif
                endif
             elseif (inct.eq.171) then
                if (iftyp.eq.1) then
                   call inferr('No Laplacian for Mopac',0)
                else
                   ostep = .false.
                   oscal = .false.
                   ofrst = .true.
                   cut   = 1.0d0
                   ipsi = 0
                   bonds  = .false.
                   atomic = .false.
                   ovrlap = .false.
                   elpot = .false.
                   molpot = .false.
                   chpot = .false.
                   dolap = .true.
                   call denmak(denok)
                   if (denok) then
                      call grdcl(npts1,npts2,1,space)
                      call inferr('Laplacian of the Density',0)
                   endif
                endif
             elseif (inct.ge.350.and.inct.le.352) then
                if (inct.eq.352.and.ihasq.eq.0) then
                     call inferr('Charge.Pot. with no charges',1)
                else
                   ostep = .false.
                   oscal = .false.
                   ofrst = .true.
                   cut   = 1.0d0
                   ipsi = 0
                   bonds  = .false.
                   dolap = .false.
                   atomic = .false.
                   ovrlap = .false.
                   if (inct.eq.350) then
                       elpot = .true.
                       molpot = .false.
                       chpot = .false.
                       call denmak(denok)
                   elseif (inct.eq.351) then
                       elpot = .false.
                       molpot = .true.
                       chpot = .false.
                       call denmak(denok)
                       call muldma(vdwr,moddma,.true.,.true.)
                   elseif (inct.eq.352) then
                       elpot = .false.
                       molpot = .false.
                       chpot = .true.
                       denok =.true.
                       call xyzcoo(1,0,0)
                   endif
                   if (denok) then
                      call grdcl(npts1,npts2,1,space)
                      call inferr('Electrostatic Potential',0)
                   endif
                endif
             elseif (inct.eq.180) then
                ipst = ipsi
                bontem = bonds
                atotem = atomic
                ovrtem = ovrlap
                ipsi = 0
                bonds  = .true.
                dolap = .false.
                atomic = .false.
                ovrlap = .false.
                elpot = .false.
                molpot = .false.
                chpot = .false.
                call denmak(denok)
                if (denok) then
                   ostep = .false.
                   oscal = .false.
                   if (iftyp.ne.1) cut = 0.1
                   call grdcl(npts1,npts2,1,space)
                   lstr = 'Molecular - spherical atomic density'
                   if (doori) then
                      lstr = 'Molecular - oriented atomic density'
                   endif
                   call inferr(lstr,0)
                   ofrst = .true.
                else
                   ipsi = ipst 
                   bonds = bontem 
                   atomic = atotem 
                   ovrlap = ovrtem 
                   call inferr('No Atomic Density available !',0)
                endif
             elseif (inct.eq.190) then
                ipst = ipsi
                bontem = bonds
                atotem = atomic
                ovrtem = ovrlap
                ipsi = 0
                bonds  = .true.
                dolap = .false.
                atomic = .true.
                ovrlap = .false.
                elpot = .false.
                molpot = .false.
                chpot = .false.
                call denmak(denok)
                if (denok) then
                   ostep = .false.
                   oscal = .false.
                   ofrst = .true.
                   if (iftyp.ne.1) cut = 0.1
                   call grdcl(npts1,npts2,1,space)
                   call inferr('Atomic Part Mol. - Atom. Density',0)
                else
                   ipsi = ipst 
                   bonds = bontem 
                   atomic = atotem 
                   ovrlap = ovrtem 
                   call inferr('No Atomic Density available !',0)
                endif
             elseif (inct.eq.200) then
                ipst = ipsi
                bontem = bonds
                atotem = atomic
                ovrtem = ovrlap
                ipsi = 0
                bonds  = .false.
                dolap = .false.
                atomic = .false.
                ovrlap = .true.
                elpot = .false.
                molpot = .false.
                chpot = .false.
                call denmak(denok)
                if (denok) then
                   ostep = .false.
                   oscal = .false.
                   ofrst = .true.
                   cut   = 1.0d0
                   call grdcl(npts1,npts2,1,space)
                   call inferr('Overlap Density',0)
                else
                   ipsi = ipst 
                   bonds = bontem 
                   atomic = atotem 
                   ovrlap = ovrtem 
                   call inferr('No Atomic Density available !',0)
                endif
             elseif (inct.eq.210) then
                if (doori) then
                   doori = .false.
                   norien = 0
                else
                   doori = .true.
                   call parori
                endif
                elpot = .false.
                molpot = .false.
                chpot = .false.
                call denmak(denok)
                if (denok) then
                   call grdcl(npts1,npts2,1,space)
                   ofrst = .true.
                endif
             elseif (inct.eq.220) then
                ostep = .false.
                keywrd = str(1:nstr)
                call tocap(keywrd,nstr)
                call planky(npts1,npts2,npts3,keywrd,.false.)
                call occup(istat)
                if (istat.eq.1) then
                    ipsi = 0
                    call denmak(denok)
                endif
                call pareul
                call proato
                call grdcl(npts1,npts2,1,space)
                ofrst = .true.
             elseif (inct.eq.230) then
                ostep = .false.
                oscal = .false.
                cut   = 1.0d0
                call homo(ipsi)
                bonds = .false.
                dolap = .false.
                ofrst = .true.
                elpot = .false.
                molpot = .false.
                chpot = .false.
                call denmak(denok)
                call grdcl(npts1,npts2,1,space)
                if (pmax.eq.0.0d0) then
                   call inferr('Nodal Plane ?!',0)
                else
                   call inferr('Highest Occupied MO',0)
                endif
             elseif (inct.eq.240) then
                ostep = .false.
                oscal = .false.
                cut   = 1.0d0
                call lumo(ipsi)
                bonds = .false.
                dolap = .false.
                ofrst = .true.
                elpot = .false.
                molpot = .false.
                chpot = .false.
                call denmak(denok)
                call grdcl(npts1,npts2,1,space)
                if (pmax.eq.0.0d0) then
                   call inferr('Nodal Plane ?!',0)
                else
                   call inferr('Lowest Unoccupied MO',0)
                endif
             elseif (inct.eq.250) then
                keywrd = str(1:nstr)
                ostep = .true.
                step = reada(keywrd,1,nstr)
             elseif (inct.eq.260) then
                keywrd = str(1:nstr)
                cut = reada(keywrd,1,nstr)
                cut = min(max(0.0d0,cut),1.0d0)
                call parstp
             elseif (inct.eq.270) then
                ostep = .false.
                scale = -1.0d0
                cut   = 1.0d0
                if (space.eq.1) then
                   call wr3inf(npts1,npts2,npts3,adjus)
                else
                   call wrinfo(npts1,npts2)
                endif
             elseif (inct.eq.280.or.inct.eq.281) then
                ostep = .false.
                scale = -1.0d0
                cut   = 1.0d0
                if (inct.eq.281) then
                   keywrd = str(1:nstr)
                   if (icubtp.eq.0.or.icubtp.eq.2) then
                      if (.not.opfil(21,keywrd(1:nstr),nstr,1,1,0)) 
     &                then
                          istat = 0
                          call inferr('Error Opening File !',0)
                      else

                          if (icubtp.eq.0) then
c Gaussian cube
                             istat = 0
                             ijag = 0
                          else
c Jaguar cube
                             istat = 1
                             ijag = 1
                          endif

                          call cubtst(21,ijag)
                          call rdcube(npts1,npts2,npts3,iposng,ipsi,
     &                                istat,21,idebug)
                          if (idofil.eq.0) idofil = 1
                          close(21)
                      endif
                   elseif (icubtp.eq.3) then
c GRID cube
                      istat = 0
                      ijag = 2
                      call cubtst(21,ijag)
                      if (ijag.ne.-1) then
                         call rdgrd(npts1,npts2,npts3,21,istat)
                         close(21)
                         idofil = 1
                         ipsi = 1
                         valcnt = 0.01d0
                         space = 1
                      endif
                   elseif (icubtp.eq.4) then
c omap file
                      istat = 0
                      ijag = 3
                      call cubtst(21,ijag)
                      if (ijag.ne.-1) then
                         call rdomap(npts1,npts2,npts3,21,istat)
                         close(21)
                         idofil = 1
                         ipsi = 1
                         valcnt = 0.01d0
                         space = 1
                      endif
                   else
c VASP cube
                      call rdvasp(npts1,npts2,npts3,iposng,istat,
     &                            nstr,1,idebug)
                      nptsmx = npts1
                      if (npts2.gt.nptsmx) nptsmx = npts2
                      if (npts3.gt.nptsmx) nptsmx = npts3

                      if (nptsmx.gt.mx3d) then
                          call allgrd(nptsmx)
                          call rdvasp(npts1,npts2,npts3,iposng,
     &                                istat,nstr,1,idebug)
                      endif

                      idofil = 0
                      ipsi = 0
                   endif
                   if (istat.eq.1) then
                      if (iposng.eq.1) ipsi = 1
                      adjus = toang
                      idisml = 1
                      valcnt = 0.01d0
                      space = 1
                      ispd = 1
                   endif
                else
                   if (space.eq.1) then
                      call tstrd3
                      call rd3inf(npts1,npts2,npts3,incp,adjus,istat)
                      ispd = 1
                   else
                      call rdinfo(npts1,npts2,incp,istat)
                      if (istat.eq.0.and.natoms.eq.0) then
                         call tstrd3
                         call rd3inf(npts1,npts2,npts3,0,adjus,istat)
                         if (istat.eq.1) then
                            valcnt = 0.01d0
                            space = 1
                            ispd = 1
                         endif
                      endif
                   endif
                endif
             elseif (inct.eq.282) then
                keywrd = str(1:nstr)
                if (.not.opfil(21,keywrd(1:nstr),nstr,1,0,0)) then
                    call inferr('Error Opening File !',0)
                else
                    call wrcube(npts1,npts2,npts3,ipsi)
                    close(21)
                endif
             elseif (inct.eq.360) then
                if (domax) then
                   domax = .false.
                else
                   domax = .true.
                endif
             elseif (inct.ne.75.and..not.
     &               (euclid.and.space.eq.0)) then
                theta = theta + inct*aincr
                phi   = phi + incp*aincr
                if (space.eq.0) then
                   if (theta.gt.90.0) then
                       theta = 90.0 - thetm
                       goto 1701
                   endif
                else
                   if (theta.gt.180.0) then
                       theta = 180.0
                       goto 1701
                   endif
                endif
                if (theta.lt.0.0) then
                    theta = thetm
                    goto 1701
                endif
                if (phi.gt.90.0) then
                    phi = 90.0 - phim
                    goto 1701
                endif
                if (phi.lt.0.0) then
                    phi = phim
                    goto 1701
                endif
1701            call parang(theta,phi)
             endif

             call resedl

             if (inct.eq.75) then
                iplwin = 0
                domol = .true.
                call xwin(xx,yy,6,str,nstr,idum1,idum2)
                call qupd
                goto 999
             endif
             goto 777
          endif
      endif
3333  continue

c--------- plot termination ---------------------------------

      if (iplot.eq.1) then
        write(iun4,'(a)')'PM2;RO0;EP;FP;SP;PG;'//esc//'.);'
      endif

      if (iplot.eq.4) call plend(title,notitle)
      if (iplot.eq.4) then
         write(iun4,'(''%%Trailer'',a)')
         write(iun4,'(''%%Pages: 1'',a)')
         write(iun4,'(a)')eot
      endif

      close(iun4)
      
      call vaxconv
      end

      subroutine prtflg
      implicit double precision (a-h,o-z)

      print*,'Molden5.1 (1996-2010) Dr. G.Schaftenaar, CMBI'
      print*,'URL: http://www.cmbi.ru.nl/molden/molden.html'
      print*,
     &'G. Schaftenaar,J.H. Noordik,J.Comp.-Aided.Mol.Des.,14 (2000) 123'
      print*,' '
      print*,'Usage: molden [ options ... ] [file1 file2 ... ]'
      print*,' '
      print*,'where options include:'
      print*,' '
      print*,'-a        no automatic cartesian -> zmat conversion'
      print*,'-b        use orbitals of first point on opt. runs'
      print*,'-c 0.5    change depth of shading, range [0.0-1.0]'
      print*,'-e        DMAREL INPUT: Use est set of parameters'
      print*,
     &'-f        PDB: build connectivity from cartesian coordinates'
      print*,
     &'-g        PDB: always calculate Helix/Sheet information'
      print*,'-geom XXXxYYY-xxx-yyy'
      print*,'          XXX and YYY are the size of the window'
      print*,'          xxx and yyy are the position of the window'
      print*,'-h        print commandline flags'
      print*,'-hoff     switch of hydrogen bonds'
      print*,'-hdmin x  mininum hydrogen bond distance (Ang)'
      print*,'-hdmax x  maximum hydrogen bond distance (Ang)'
      print*,'-hamin x  mininum hydrogen bond angle (Degrees)'
      print*,'-hamax x  maximum hydrogen bond angle (Degrees)'
      print*,'-i opt    fdat files: '
      print*,'          opt=1 standardise H-C, H-N'
      print*,'          opt=2 1 + standardise phenyl rings'
      print*,'-j num    maximum number of gifs to write'
      print*,'-k num    select color of labels (0-15)'
      print*,'          (gmolden only).'
      print*,'-l        dont display molden logo'
      print*,'-n        dont add hydrogens to PDB file'
      print*,'-m        turn off the beep sounds'
      print*,'-o fname  plotfilename (default=plot)'
      print*,'-p 2.0    change perspective value (def: 13.0)'
      print*,'-r fname  read file with per line; '
      print*,'          atom color(1-15) VandeWaalsRadius, (- = skip)'
      print*,'          background color(1-15)'
      print*,'          palette red #CF54FD ...   (14 colors)'
      print*,'-s 4.0    scale amplitude of normal vibrations'
      print*,'-t        read ascii MOPAC93 Chem3D style'
      print*,'-u        With GAMESS-US optimisation output, molden'
      print*,'          generates a z-matrix, (def: read from output)'
      print*,'-v        print verbose information'
      print*,'-w opt    write all points of a movie to a file:'
      print*,'          opt specifies format; xyz(=1) zmat(=2,mopac)'
      print*,'          VRML2.0(=3)'
      print*,'-x file   read in file with spherical atomic densities'
      print*,'-y 1.0    threshold for printing displacement vectors'
      print*,'          of normal modes to postscript file'
      print*,'-z        create high quality opengl coils'
      print*,'-A        Keep order of atoms when creating a Z-matrix'
      print*,'-C        Color postscript (default=mono, except PDB)'
      print*,'-D opt    DMA mode:'
      print*,'          0 = atomic sites only (default)'
      print*,'          1 = atomic+halfway-bond sites'
      print*,'          2 = no shift of overlap dens. of conn. atoms'
      print*,'-E        DMAREL input: use coordinates from multipoles'
      print*,'-F        gmolden: Use all opengl code, (line drawing)'
      print*,'-G 0.6    Grid width colour coded ESP potential map'
      print*,'-H        GAMESS-US: do normal modes when HSSEND=.TRUE.'
      print*,'-I        dont use shaders, if available'
      print*,'-L        display both neg. and pos. contour in space'
      print*,'          plot of the laplacian'
      print*,'-1        use only the lower half of the cubic grid'
      print*,'          used for the space type plot'
      print*,'-2        use only the upper half of the cubic grid'
      print*,'          used for the space type plot'
      print*,'-M        MonoChrome postscript'
      print*,'-N        Check for mpi, to run ambfor/ambmd '
      print*,'          in parallel '
      print*,'-O        switch off multiple structures handling'
      print*,'-P        PDB: treat all input files as PDB files'
      print*,'-Q        support for older StarNet xwin32 (ver. 6)'
      print*,'-R npts   adjust the gridsize in points'
      print*,'-S        start with shade off'
      print*,'-T        treat all input files as TINKER xyz files'
      print*,'-U        do not use opengl shaders'
      print*,'-X        use with XMOL cartesian format input'
      print*,'-V fname  VRML density filename'
      print*,'-W        Write VRML2.0 instead of VRML1.0'
      print*,'-Z        Map the Z-matrix file mapfile onto crystal'
      print*,'          mapfile contains Z-matrix followed by keyword'
      print*,'          MAP and per line an integer that maps a'
      print*,'          Z-matrix line onto a cartesian line'
      print*,'-=        Use gamess-us dialect of gaussian zmat writing'
      print*,' '

      return
      end

      subroutine setid(nres,istart,isurf,iresid,ipdbt)
      implicit double precision (a-h,o-z), integer ( i-n)
      common /athlp/ iatoms, mxnat
      integer*2 ipdbt
      dimension isurf(*),iresid(*),ipdbt(*)

      do i=istart+1,mxnat
         isurf(i) = 0
         ipdbt(i) = 0
         iresid(i) = -nres
      end do

      return
      end

      subroutine getpar(filenm,iverb,iball,ideltm)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (mxel=100)
      integer getlin
      character*(*) filenm
      character*137 str
      character*2 catom, catomt, tolowf
      character*2 elemnt
      common /elem/elemnt(mxel)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /backg/ ibgcol,ibgclo,ibgmod



      iuntmp = iun2
      open(unit=46,form='formatted',
     & file=filenm,status='unknown',err=300)
      iun2 = 46

      do while (getlin(1).eq.1)

c Get Atom label

          idly = 0
          ibg = 0
          ibg2 = 0
          ibgm = 0
          ipal = 0
          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.ne.1) goto 100
          if (nstr.ge.4) then
             if (icdex(str,'dela').ne.0) idly = 1
             if (icdex(str,'back').ne.0.and.icdex(str,'ogl').eq.0.and.
     &           icdex(str,'mode').eq.0)
     &           ibg = 1
             if (icdex(str,'oglback').ne.0) ibg2 = 1
             if (icdex(str,'mode').ne.0) ibgm = 1
             if (icdex(str,'pale').ne.0) ipal = 1
             if (icdex(str,'oldr').ne.0) then
                 iball = 0
                 goto 270
             endif
             if (icdex(str,'newr').ne.0) then
                 iball = 1
                 goto 270
             endif
          else
             iat = 0
             if (nstr.eq.1) then
                catomt(2:2) = str(1:1)
                catomt(1:1) = ' '
             else
                catomt = str(1:2)
             endif
             catom = tolowf(catomt)
             do j=1,100
                if (catom .eq. tolowf(elemnt(j))) iat = j
             end do
             if (catom.eq.'xx') iat = 99
             if (iat.le.0.or.iat.gt.mxel) then
                 print*,'unrecognized atom'
                 goto 100
             endif
             if (iverb.eq.1) then
               write(iun3,'(a,a2,a,i2,a,f6.3)') 'Defaults for Atom: ',
     &        catom,' color=',icol(iat),' VandeWaalsRadius=',vdwr(iat)
             endif
          endif

c Get Color Palette

          if (ipal.eq.1) then
             
             call parsfn('Black',5,3)
             do i=2,15
                 ktype = nxtwrd(str,nstr,itype,rtype)
                 if (ktype.eq.1) then
                     if (nstr.eq.1) then
                        if (str(1:1).eq.char(92)) then
                           if (getlin(1).eq.1) then
                              ktype = nxtwrd(str,nstr,itype,rtype)
                           endif
                        endif
                     endif
                     call parsfn(str,nstr,3)
                     ipal = ipal + 1
                 endif
             end do
             call parsfn('White',5,3)
             ipal = ipal + 1
             if (ipal.ne.16) print*, 'Not enough Colors, Need 14'
             goto 270
          endif

c Get Atom color

          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.ne.2.and.ktype.ne.1) goto 100
          if (ktype.eq.1) then
              if (nstr.eq.1.and.str(1:1).eq.'-') then
                 goto 200
              else
                 goto 100
              endif
          endif
          if (idly.eq.1) then
              ideltm = itype
              goto 270
          endif
          if (ibg.eq.1) then
              if (itype.le.15.and.itype.ge.0) then
                  ibgcol = itype
                  goto 270
              endif
          elseif (ibg2.eq.1) then
              if (itype.le.15.and.itype.ge.0) then
                  ibgclo = itype
                  goto 270
              endif
          elseif (ibgm.eq.1) then
              if (itype.le.1.and.itype.ge.0) then
                  ibgmod = itype
                  goto 270
              endif
          else
              if (itype.le.15.and.itype.gt.0) then
                  icol(iat) = itype
              else
                  print*,'Color value out of range [1-15]'
              endif
          endif
          
c Get Atom van de waals radius

200       ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.ne.3.and.ktype.ne.1) goto 100
          if (ktype.eq.1) then
              if (nstr.eq.1.and.str(1:1).eq.'-') then
                 goto 250
              else
                 goto 100
              endif
          endif
          if (rtype.gt.0.0d0) then
              vdwr(iat) = rtype
          else
              print*,'Negative Van de Waals Radius'
          endif
250       if (iverb.eq.1) then
              write(iun3,'(a,a2,a,i2,a,f6.3)') 'New values: ',catom,
     &               ' color=',icol(iat),' VandeWaalsRadius=',vdwr(iat)
          endif
270       continue
      end do

      close(46)
      iun2 = iuntmp

      return

300   print*,'Error opening file ',filenm
      iun2 = iuntmp
      return

100   print*,'Atom     Color     VandeWaalsRadius'
      print*,'-----------------------------------'
      print*,'Atom    0-15 or -       value or -'
      close(46)
      iun2 = iuntmp

      return
      end

      subroutine anid(coo,ianz)
c This is really anim
      implicit double precision (a-h,o-z)
      parameter (maxdm=20)
      integer dolabs,fancy,persp,shade,atcol,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo
      common /distmn/ rdm(maxdm),idmon(2,maxdm),ndm
      common /athlp/ iatoms, mxnat
      dimension coo(3,*),ianz(*)

      call dumzm(coo,ianz,iatoms)
      ndmtmp = ndm
      call clrmon
      ndm = ndmtmp
      call domcon(1,1)

      call qupd

      return
      end

      subroutine initim
      implicit double precision (a-h,o-z)
      common /times/ t1,t2,t3,t4,t5
    
      t1 = 0.0d0
      t2 = 0.0d0
      t3 = 0.0d0
      t4 = 0.0d0
      t5 = 0.0d0

c To call timimg routine
c 
c common /times/ t1,t2,t3,t4,t5
c dimension cpu(2)
C
c call dtime(cpu)
c t1 = t1 + cpu(1)
c
      return
      end

      subroutine prttim
      implicit double precision (a-h,o-z)
      common /times/ t1,t2,t3,t4,t5
    

      print*,'t1=',t1
      print*,'t2=',t2
      print*,'t3=',t3
      print*,'t4=',t4
      print*,'t5=',t5

      return
      end

      logical function opfil(iun,filenm,lenfn,iform,iold,isil)
      implicit double precision (a-h,o-z)
      character*(*) filenm
      character*7 stat

      opfil = .true.

      if (iold.eq.1) then
          stat = 'old'
      else
          stat = 'unknown'
      endif

      if (lenfn.eq.0.or.filenm.eq.' ') then
          call inferr('Invalid Filename !',0)
          opfil = .false.
      else
          close(iun)
          if (iform.eq.1) then
             open(unit=iun,form='formatted',file=filenm,
     &               status=stat,err=100)
          else
             open(unit=iun,form='unformatted',file=filenm,
     &               status=stat,err=100)
          endif
      endif

      return

100   if (isil.eq.0) then
         print*,'=',filenm,'='
         call inferr('Error Opening File ! ',0)
         call messg(12)
      endif

      opfil = .false.
      return
      end

      subroutine prtved(
     &                 vectrs,vectrb,eiga,eigb,ncols,ncolb)

c THIS IS REALLY prtvec

      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numatm=2000)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /orbhlp/ mxorb,iuhf,ispd
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      real eiga,eigb
      dimension vectrs(*),vectrb(*),eiga(*),eigb(*)


      write(iun3,*)' '
      write(iun3,*)'****** Vectors ***********'
      write(iun3,*)' '

      if (iuhf.eq.1) then
        write(iun3,*)'Alpha set'
        write(iun3,*)' '
      endif
      write(iun3,*) (eiga(i),i=1,ncols)
      write(iun3,*)' '
      call prev(vectrs,ncols,norbs,mxorb)

      if (iuhf.eq.1) then
        write(iun3,*)' '
        write(iun3,*)'Beta set'
        write(iun3,*)' '
        write(iun3,*) (eigb(i),i=1,ncolb)
        write(iun3,*)' '
        call prev(vectrb,ncolb,norbs,mxorb)
      endif

      return
      end

      subroutine clkbcd(istsurf,incp,ifogl,iresid,coo,
     &                  icalf,ncalf,reson)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (mxplev=5)
      integer reson
      common /potlev/grdw, plevel(mxplev),ipcol(mxplev+1),nplev
      dimension iresid(*),coo(3,*),icalf(6,*),reson(*)

      if (istsurf.eq.1) then
         idosurf = 1
         grdw = 1.3d0
      else
         idosurf = 0
      endif

      if (incp.ne.0) then

       if (iresid(incp).gt.0) then

          if (reson(iresid(incp)).eq.1) then
             call actami(iresid(incp),0,0,idosurf)
          else
             call actami(iresid(incp),0,1,idosurf)
          endif

       elseif (iresid(incp).gt.-4) then

          rdum = 10000.0d0
          jres = 0

          do i=1,ncalf

            dist = dist2(coo(1,incp),coo(1,icalf(1,i)))
            if (dist.lt.rdum) then
               rdum = dist
               jres = i
            endif

          end do

          if (jres.ne.0) then
             if (reson(jres).eq.1) then
                call actami(jres,0,0,idosurf)
             else
                call actami(jres,0,1,idosurf)
             endif
          endif

       endif

       if (idosurf.eq.1) then
           if (ifogl.eq.1) then
              call initsrf
           else
              call connlp(1.0d0,0,4)
           endif
       endif

      endif

      return
      end

      subroutine newfid(idebug,istats,incp,ioadd,ioatms,nstrt,namols,
     &                  nxtmf,ipdbon,namls,iof,iaton,iatclr,iresid,
     &                  ncalf)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (mxmmul=100)
      character fniun*256
      common /fnunit/ fniun
      common /multim/ imulm, nmulm,ihasqm(mxmmul)
      common /zmfrst/ ihaszm, nz, mxzat
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /athlp/  iatoms, mxnat
      logical opfil, filok
      dimension iaton(*),iatclr(*),iresid(*)

      istats = 1

      call addnew(incp,ioadd,ioatms,nstrt,namols,ncalf,iresid)

      if (iof.eq.1) then
         lenf = linlen(fniun)
         filok = opfil(48,fniun,lenf,1,1,0)
      else
         filok = .true.
      endif

      if (filok) then

         if (iof.eq.1) then
            iun1 = 48
            rewind iun1
         endif

         if (ioadd.eq.1) then

            ipdbont = ipdbon
            ihaszmt = ihaszm

            istat = 0
            if (nxtmf.eq.0) call prsgmf(1)
            call rewmf
            call rdmol(idebug,ipdbon,ioadd,istat)

c                     istat = 0 error
c                           = 1 ok
c                           = 2 cell info (not applicable here)

       
            if (istat.eq.3) then

               rewind iun1
               istat = -1
               call rdmol(idebug,ipdbon,ioadd,istat)
               goto 2001

            elseif (istat.eq.1) then

               ipdbon = ipdbont
               ihaszm = ihaszmt
               goto 2001

            endif

            rewind iun1

            call getxyz(igetxy,heat,1)
            if (igetxy.eq.1) then

               call numhet(nhmol)
               call setis(nhmol+1,ioatms)
               goto 2001

            else

               rewind iun1
               goto 10

            endif

         else
            goto 10
         endif

      else
         call inferr('Error Opening File !',0)
      endif

c found and read new mol (xyz,mol2), disable oadd
c for pdb files we dont get here but goto 10 and pass through 
c like any new file, but with oadd still set

2001  continue

      ioadd = 0
      namls = namols / 15
      do i=ioatms+1,iatoms
         iaton(i) = 1
         iatclr(i) = namols - 15*namls
      end do

      call ribbs

      do l=1,4
         call acthlp(l,0,0)
      end do

      return

10    istats = 0
      return
      end

      subroutine nwfil(str,nstr)
      implicit double precision (a-h,o-z), integer (i-n)
      character str*100
      character fniun*256
      common /fnunit/ fniun
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /pnthlp/ ipoints,ipnt
      logical ctoz,molpot,elpot,chpot,opfil
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      logical first
      common /denfir/ first

      first = .true.
c binary file ? zou eigenlijk bijgehouden moeten worden
c nxtmf ? misschien toevoegen aan common /multim/, zou communicatie
c met C makkelijker maken, check .trr etc files compatibility

      fniun = str(1:nstr)
      lenf = linlen(fniun)

      close(iun2,err=10)
10    continue

      if (opfil(48,fniun,lenf,1,1,0)) then

         iun2 = 48
         if (iftyp.eq.5.or.iftyp.eq.8.or.iftyp.eq.9) then
            call rewmf
         else
            rewind iun2
         endif

      else
         call inferr('Error Opening File !',0)
      endif

      return
      end

      subroutine addnew(iaddf,ioadd,ioatms,nstrt,namols,ncalf,iresid)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (mxmmul=100)
      common /multim/ imulm, nmulm,ihasqm(mxmmul)
      common /athlp/  iatoms, mxnat
      dimension iresid(*)

      if (iaddf.eq.1) then

c ADD NEW file to existing file

         ioadd = 1
         ioatms = iatoms

c   start adding new atoms at last real atom, disregard some fake atoms
c   like: 0 > iresid > -3 (secondary structure display help atoms)

         do i=1,iatoms
            if (iresid(i).le.0.and.iresid(i).ge.-3) then
               ioatms = i - 1
               goto 10
            endif
         end do

10       namols = namols + 1
         iatoms = ioatms
         nstrt = ncalf + 1

      else

c REPLACE with NEW file

         ioadd = 0
         ioatms = 0
         namols = 1
         nmulm = 1
         nstrt = 1

      endif

      return
      end 

      subroutine acthld(iop1,iop2,iop3,ihet)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension ihet(*)

      if (ihet(iop1).eq.1) call acthel(1,iop1-1,iop2,iop3)

      return
      end

      subroutine acttod(iop1,iop2,iop3,ihet)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension ihet(*)

c     Toggle visibility secundary structure elements

      if (ihet(iop1).eq.1) then
         ihet(iop1) = 0
         call acthel(0,iop1-1,iop2,iop3)
      else
         ihet(iop1) = 1
         call acthel(1,iop1-1,iop2,iop3)
      endif

      return
      end

      subroutine acttad(incp,kleur,idosurf,ihet)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension ihet(*)

c     Toggle visibility hetatm residue

      ihind = abs(incp)+1

      if (ihet(ihind).eq.1) then
         call actami(incp,kleur,0,idosurf)
      else
         call actami(incp,kleur,1,idosurf)
      endif

      return
      end

      subroutine setard(iopt,iopval,ioatms,
     &                  ianz,iaton,iatclr,iresid,iconn,qat,
     &                  ihet,iclhet,reson,iams,ihets,irsnr,
     &                  ncalf,issdon,
     &                  scal,scali,fscal,smag,
     &                  xv,yv,zv,pincr,natc,ichx,icrtp,
     &                  ipoints,ngeoms)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      parameter (mxel=100)
      parameter (numcal=50000)
      parameter (mxheta=150)
      common /athlp/  iatoms, mxnat
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /qmchar/ qch(numat1),ihasesp
      integer reson
      dimension ianz(*),iaton(*),iresid(*),iatclr(*),iconn(mxcon+1,*),
     &          qat(*),ihet(*),iclhet(*),reson(*),iams(*),ihets(*),
     &          irsnr(*)

      if (iopt.eq.1) then

          ir = iresid(iopval)
          do i = 1,iatoms
             if (iresid(i).eq.ir) iaton(i) = 2
          end do

      elseif (iopt.eq.2) then

          do i=1,iatoms
             iaton(i) = 1
             iatclr(i) = icol(ianz(i))
          end do

      elseif (iopt.eq.3) then

          do i=ioatms+1,iatoms
             ir = iresid(i)
             if (ir.gt.0.or.ir.lt.-3) then
                iatclr(i) = iopval
             endif
          end do

      elseif (iopt.eq.4) then

          do i=1,iatoms
             iatclr(i) = iopval
          end do

      elseif (iopt.eq.5) then

          nstor = mxnat - natc
          do i=1,natc
             do j=1,iconn(1,i)+1
                iconn(j,nstor+i) = iconn(j,i)
             end do
          end do

      elseif (iopt.eq.6) then

         do i=1,mxnat
            qat(i) = 0.0d0
         end do

      elseif (iopt.eq.7) then

         do i=1,iopval
            qat(i) = qch(i)
         end do

      elseif (iopt.eq.8) then

         do i=iopval+1,iatoms
            iaton(i) = 0
         end do

      elseif (iopt.eq.9) then

         do i=1,iatoms
             iaton(i) = 1
         end do

      elseif (iopt.eq.10) then

         issdon = 0
         ncalf = 0

         do i=1,numcal
            irsnr(i) = 0
            iams(i) = 0
         end do

      elseif (iopt.eq.11) then

         do i=1,mxheta
            ihet(i) = 1
            ihets(i) = 0
            iclhet(i) = 1
         end do

      elseif (iopt.eq.12) then

         do i=1,numcal
            reson(i) = 1
         end do

      elseif (iopt.eq.13) then

         smag   = 1.0d0
         xv     = 0.0d0
         yv     = 0.0d0
         fscal  = 1.0d0
         scal   = 1.0d0
         scali  = 1.0d0

      elseif (iopt.eq.14) then

         zv = scali
         pincr = 0.02d0*scali

      elseif (iopt.eq.15) then
 
         ichx = 0
         icrtp = 0
         ipoints = 0
         ngeoms = 0

      endif

      return
      end

      subroutine pckrod(istat,iatsel,iresid,iconn,ianf,islu,nchain)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      common /athlp/ iatoms, mxnat
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt
      dimension iresid(*),iconn(mxcon+1,*),icn(mxcon)
      dimension ipdb(mxsym),ihpdb(mxhsym*3),ianf(*),islu(*)

      istat = 1
      irsp = iresid(iatsel)

      if (irsp.lt.-3.or.irsp.gt.0) then

         iscst = 0
         ifrs = 1
         itmp = 0

         if (irsp.lt.-3) then

            do i=1,iatoms
               if (iresid(i).eq.irsp) then
                  if (ifrs.eq.1) then
                     iscst = i
                     ifrs = 0
                  endif
                  itmp = i
               endif
            end do

         else if (irsp.gt.0) then

              ich = 0
              do i=1,nchain
                 if (irsp.ge.ianf(i).and.irsp.le.islu(i)) then
                   ich = i
                 endif
              end do

              do i=ianf(ich),islu(ich)

                 call getpdb(i,ipdb,ihpdb)

                 do j=1,mxsym
                    if (ipdb(j).ne.0) then
                        if (ifrs.eq.1) then
                           iscst = ipdb(j)
                           ifrs = 0
                        endif
                    endif
                 end do

                 do j=1,mxhsym*3
                    if (ihpdb(j).ne.0) then
                        itmp = ihpdb(j)
                    endif
                 end do

              end do

         endif

         if (iscst.ne.0) iscst = iscst - 1

         if (itmp.ne.0) then
            nscnd = itmp - iscst
            do i=1,nscnd
               iat = iscst+i
               ncn = 0
               do j=1,iconn(1,iat)
                  if (iconn(1+j,iat).gt.0) then
                     if (irsp.lt.-3) then
                        if (iresid(iconn(1+j,iat)).eq.irsp) then
                           ncn = ncn + 1
                           icn(ncn) = iconn(1+j,iat)
                        endif
                     else
                        ir = iresid(iconn(1+j,iat))
                        if (ir.ge.ianf(ich).and.ir.le.islu(ich)) then
                           ncn = ncn + 1
                           icn(ncn) = iconn(1+j,iat)
                        endif
                     endif
                  endif
               end do
               iconn(1,iat) = ncn
               do j=1,ncn
                  iconn(1+j,iat) = icn(j)
               end do
            end do
            call aln2ml(2,istat)
         else
            print*,'error '
         endif
      else
         istat = 0
         print*,'only HETATM records are allowed'
         call messg(13)
      endif

      return
      end

      subroutine calelt
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg

      if (norbs.gt.0.and.isgau.eq.1) then
         itmp = 0
         call calelc(2,itmp,itmp,itmp,0)
      elseif (ihasq.gt.0) then
         call caldip
      endif

      return
      end

      subroutine putxyz
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg

c     restore wavefunction coordinates

      if (norbs.gt.0.and.isgau.eq.1) then
         if (ihsdp.ne.2) then
            call xyzcoo(1,0,0)
            call haszm(.false.)
            call doconn
         endif
      endif

      return
      end

      subroutine calelc(iopt,isfdon,dolabs,space,moddma)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numatm=2000)
      parameter (mxvalc=10)
      parameter (mxel=100)

      logical denok,dolap,bonds,ovrlap,atomic,
     &        doiso
      logical molpot,elpot,ctoz,chpot
      logical valenc,doori

      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /srfhlp/ edge,ctval(mxvalc),pxyz(3),nvalc,nspts,istyp
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      integer dolabs,space

      if (((iftyp.ge.2.and.iftyp.le.4).or.iftyp.eq.8.or.
     &     (iftyp.eq.5.and.isgau.eq.1)
     &    ).and.norbs.ne.0) then

          call xyzcoo(1,0,0)
          call haszm(.false.)
          call doconn

          call curs(1)

          ipsi = 0
          dolap = .false.
          bonds = .false.
          ovrlap = .false.
          atomic = .false.

          call denmak(denok)

          doiso = .false.
          if (istyp.eq.1) doiso = .true.

          if (iopt.eq.0) then
             call espchrg(ctval,nvalc,doiso,.false.)
          endif


          if (iopt.eq.1) then
             call muldma(vdwr,moddma,.true.,.true.)
             call espchrg(ctval,nvalc,doiso,.true.)
          endif

          if (iopt.eq.2) then
             call muldma(vdwr,moddma,.true.,.false.)
             call wrxyz(2)
          endif

          call curs(0)

          space = 0
          if (isfdon.eq.3) isfdon = 0

          iqon = 3
          dolabs = 1
          call butset(1,14,0)
          call inferr('Wrote file esp.xyz !',0)
      else
          call inferr('Can not calculate charges !',0)
      endif

      return
      end

      subroutine prntd(iaton)
      implicit double precision (a-h,o-z), integer (i-n)
      common /athlp/ iatoms, mxnat
      dimension iaton(*)

      print*,'mxnat=',mxnat,'iaton=',(iaton(i),i=1,mxnat)
      return
      end

      subroutine parcom
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (numatm=2000)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common /mopac/ nfirst(numatm),nlast(numatm),
     &               npq(numatm),pqn(numatm),emus(numatm),
     &               emup(numatm),emud(numatm),consts(numatm),
     &               constp(numatm),constd(numatm),npqref(54)
      common /pseudo / ipseud,ivale(numatm)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /eulx/   ca,cb,sa,sb,cc,sc
      common /atopla/ xsym(numatm),ysym(numatm),zsym(numatm),
     &                isym(numatm)
      real freq,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      parameter (maxmol=100)
      common /mfdata/ nmols,imol,iendmf,ielin,mollin(maxmol)
      parameter (mxplev=5)
      common /potlev/ grdw, plevel(mxplev),ipcol(mxplev+1),nplev
      common /condens/ icont(numat1),ncont
      common /pbc/ abc(3),ibox,icell,igfmap

      call parptr(1,freq,fdum,idum)
      call parptr(143,exx,fdum,ihasd)
      call parptr(112,frint,ramint,ihasi)
      call parptr(175,fdum,fdum,nfirst)
      call parptr(176,fdum,fdum,ipseud)
      call parptr(177,px,fdum,idum)
      call parptr(178,ca,fdum,idum)
      call parptr(179,xsym,fdum,idum)
      call parptr(181,fdum,fdum,nmols)
      call parptr(20,grdw,fdum,idum)
      call parptr(12,fdum,fdum,icont)
      call parptr(13,fdum,fdum,ncont)
      call parptr(195,abc,fdum,idum)

      return
      end

      subroutine parhet(i,str)
      implicit double precision (a-h,o-z)
      parameter (mxheta=150)
      character*137 str
      character*3 hetz
      common /clfstr/ ihashz,ihetq(mxheta),ihqset(mxheta),ihhadd(mxheta)
     &                ,hetz(mxheta)

      hetz(i) = str(1:3)

      return
      end

      subroutine prttyd(ianz,ipdbt,ityp)
      implicit double precision (a-h,o-z)
      parameter (mxamb=1590)
      integer*2 ityp,ipdbt
      common /athlp/ iatoms, mxnat
      dimension ianz(*),ipdbt(*),ityp(*),imap(mxamb)

      do i=1,mxamb
         imap(i) = 0
      end do

      do i=1,iatoms
         if (ityp(i).gt.0.and.ipdbt(i).gt.0) then
             if (ianz(i).eq.1) then
                im = imap(ityp(i))
                if (im.ne.0) then
                   if (ipdbt(i).lt.im) then
                      imap(ityp(i)) = ipdbt(i)
                   endif
                else
                   imap(ityp(i)) = ipdbt(i)
                endif
             else
                imap(ityp(i)) = ipdbt(i)
             endif
         endif
      end do

      do i=1,mxamb
         write(*,'(i4,1x,i3)') i,imap(i)
      end do

      return
      end

      subroutine prthtm()
      implicit double precision (a-h,o-z), integer (i-n)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /athlp/ iatoms, mxnat
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      logical opfil
      character*4 ani
      character*11 movfil
      character*28 st1
      character*44 st2
      character*31 st3
      character*16 st4
      character*49 st5
      character*37 st6
      character*24 st7
      character*2 st8
      character*1 qs,qd
      character*26 jav0  
      character*22 jav1  
      character*17 jav2  
      character*48 jav3  
      character*31 jav4  
      character*34 jav5  
      character*42 jav6  
      character*40 jav7  
      character*9  jav8  
      character*6  jav9  
      character*26 jav10  
      character*1  jav11  
      character*1  jav12  
      character*59 jav13  
      character*1  jav14  
      character*17 jav15  
      character*46 jav16  
      character*70 jav17  
      character*65 jav17b 
      character*1  jav18  
      character*1  jav19  
      character*29 jav20  
      character*24 jav21  
      character*24 jav22  
      character*53 jav23  
      character*36 jav24  
      character*54 jav26  
      character*18 jav27  
      character*19 jav28  
      character*11 jav29  
      character*52 jav30  
      character*51 jav31  
      character*51 jav32  
      character*50 jav33  
      character*4  jav34  
      character*11 jav35  
      character*1  jav36  
      character*26 jav37  
      character*26 jav38  
      character*25 jav39  
      character*1  jav40  
      character*26 jav41  
      character*27 jav42  
      character*26 jav43  
      character*53 jav44  
      character*26 jav45  
      character*12 jav46  
      character*43 jav47  
      character*28 jav48  
      character*28 jav49  
      character*45 jav50  
      character*37 jav51  
      character*7  jav52  
      character*39 jav53  
      character*18 jav54  
      character*37 jav55  
      character*16 jav56  
      character*7  jav57  
      character*4  jav58  
      character*19 jav59  
      character*28 jav60  
      character*18 jav61  
      character*4  jav62  
      character*1  jav63  
      character*1  jav64  
      character*22 jav65  
      character*30 jav66  
      character*43 jav67  
      character*45 jav68  
      character*41 jav69  
      character*25 jav69b 
      character*28 jav70  
      character*30 jav71  
      character*27 jav72  
      character*37 jav73  
      character*22 jav74  
      character*7  jav75  
      character*1  jav76  
      character*1  jav77  
      character*20 jav78  
      character*26 jav79  
      character*26 jav80  
      character*42 jav81  
      character*21 jav82  
      character*25 jav83  
      character*31 jav84  
      character*31 jav85  
      character*27 jav86  
      character*14 jav87  
      character*1  jav88  
      character*1  jav89  
      
      qs = char(39)
      qd = char(34)

      jav0 = "var canvas = new Object();"
      jav1 = "var img = new Image();"
      jav2 = "function draw() {"
      jav3 = "     canvas = document.getElementById("//
     &       qd//"canvas"//qd//");"
      jav4 = "     img.onload = function () {"
      jav5 = "        if (canvas.getContext) {"
      jav6 = "        var ctx = canvas.getContext("//qd//"2d"//qd//");"
      jav7 = "        ctx.drawImage(img,0,0,1000,500);"
      jav8 = "        }"
      jav9 = "     }"
      jav10 = "     img.src = "//qd//"spec.jpg"//qd//";"
      jav11 = "}"
      jav12 = " "
      jav15 = "function init() {"
      jav16 = "   canvas = document.getElementById("//qd
     &        //"canvas"//qd//");"
      jav17 = "   if (canvas.attachEvent) canvas.attachEvent("//qd//
     &        "onclick"//qd//", getPosition);"
      jav17b = "   else canvas.addEventListener("//qd//"mousedown"//
     &         qd//", getPosition, false);"
      jav18 = "}"
      jav19 = " "
      jav20 = "function getPosition(event) {"
      jav21 = "   var x = new Number();"
      jav22 = "   var y = new Number();"
      jav23 = "   var canvas = document.getElementById("//qd
     &        //"canvas"//qd//");"
      jav24 = "   var parent = canvas.offsetParent;"
      jav26 = "   if (event.x != undefined && event.y != undefined) {"
      jav27 = "      x = event.x;"
      jav28 = "      y = event.y;"
      jav29 = "   } else {"
      jav30 = "      x = event.clientX + document.body.scrollLeft +"
      jav31 = "               document.documentElement.scrollLeft;"
      jav32 = "      y = event.clientY + document.body.scrollTop +"
      jav33 = "               document.documentElement.scrollTop;"
      jav34 = "   }"
      jav35 = "   x -= 10;"
      jav36 = " "
      jav37 = "   x -= parent.offsetLeft;"
      jav38 = "   x -= canvas.offsetLeft;"
      jav39 = "   y -= canvas.offsetTop;"
      jav40 = " "
      jav41 = "   var fnd = new Number();"
      jav42 = "   var xfnd = new Number();"
      jav43 = "   var temp = new Array();"
      jav44 = "   var lists = document.getElementsByTagName("//qd
     &        //"area"//qd//");"
      jav45 = "   var crd = new String();"
      jav46 = "   fnd = -1;"
      jav47 = "   for (var i = 0; i < lists.length; i++) {"
      jav48 = "      crd = lists[i].coords;"
      jav49 = "      temp = crd.split("//qd//","//qd//");"
      jav50 = "      for (var j = 0; j < temp.length; j++) {"
      jav51 = "         temp[j] = parseInt(temp[j]);"
      jav52 = "      }"
      jav53 = "      if (x > temp[0] && x < temp[2]) {"
      jav54 = "          fnd = i;"
      jav55 = "          xfnd = (temp[0]+temp[2])/2;"
      jav56 = "          break;"
      jav57 = "      }"
      jav58 = "   }"
      jav59 = "   if (fnd != -1) {"
      jav60 = "      eval(lists[fnd].href);"
      jav61 = "      dline(xfnd);"
      jav62 = "   }"
      jav63 = "}"
      jav64 = " "
      jav65 = "function dline(xfnd) {"
      jav66 = "      if (canvas.getContext) {"
      jav67 = "         var ctx = canvas.getContext("//qd
     &                            //"2d"//qd//");"
      jav68 = "         ctx.canvas.width = ctx.canvas.width;"
      jav69 = "         ctx.drawImage(img,0,0,1000,500);"
      jav69b = "         ctx.beginPath();"
      jav70 = "         ctx.moveTo(xfnd,0);"
      jav71 = "         ctx.lineTo(xfnd,500);"
      jav72 = "         ctx.lineWidth = 2;"
      jav73 = "         ctx.strokeStyle = "//qd//"#ff0000"//qd//";"
      jav74 = "         ctx.stroke();"
      jav75 = "      }"
      jav76 = "}"
      jav77 = " "
      jav78 = "function dolin(id) {"
      jav79 = "   var crd = new String();"
      jav80 = "   var temp = new Array();"
      jav81 = "   var area = document.getElementById(id);"
      jav82 = "   crd = area.coords;"
      jav83 = "   temp = crd.split("//qd//","//qd//");"
      jav84 = "   var linl = parseInt(temp[0]);"
      jav85 = "   var linu = parseInt(temp[2]);"
      jav86 = "   var lin = (linl+linu)/2;"
      jav87 = "   dline(lin);"
      jav88 = "}"
      jav89 = " "

      st1 = "document.jmola.script("//qs//"load "
      st2 = ";wireframe 30;spacefill 100;select hydrogen;"
      st3 = "spacefill 70;zoom 100;select *;"
      st4 = "set specular on;"
      st5 = "anim mode loop 0 0;animation on;animation fps 50;"
      st6 = "set echo bottom left;echo "//qd//"frequency "
      st7 = ";set frank off;"//qs//");dolin("
      st8 = ");"


      if (opfil(48,"index.html",10,1,0,0)) then
         write(48,'(a)') "<HTML>"
         write(48,'(a)') "<HEAD>"
         write(48,'(a)') "<!--[if IE]><script type="//qd//
     &                    "text/javascript"//qd//" src="//qd//
     &                    "/excanvas.js"//qd//"></script><![endif]-->"
         write(48,'(a)') "<BODY onload="//qd//"init();draw();"//qd//">"
         write(48,'(a)') "<TITLE>IR frequency</TITLE>"
         write(48,'(a)') "<script language="//qd//"javascript"//qd//">"

         write(48,'(a)')  jav0  
         write(48,'(a)')  jav1  
         write(48,'(a)')  jav2  
         write(48,'(a)')  jav3  
         write(48,'(a)')  jav4  
         write(48,'(a)')  jav5  
         write(48,'(a)')  jav6  
         write(48,'(a)')  jav7  
         write(48,'(a)')  jav8  
         write(48,'(a)')  jav9  
         write(48,'(a)')  jav10  
         write(48,'(a)')  jav11  
         write(48,'(a)')  jav12  
         write(48,'(a)')  jav15  
         write(48,'(a)')  jav16  
         write(48,'(a)')  jav17  
         write(48,'(a)')  jav17b 
         write(48,'(a)')  jav18  
         write(48,'(a)')  jav19  
         write(48,'(a)')  jav20  
         write(48,'(a)')  jav21  
         write(48,'(a)')  jav22  
         write(48,'(a)')  jav23  
         write(48,'(a)')  jav24  
         write(48,'(a)')  jav26  
         write(48,'(a)')  jav27  
         write(48,'(a)')  jav28  
         write(48,'(a)')  jav29  
         write(48,'(a)')  jav30  
         write(48,'(a)')  jav31  
         write(48,'(a)')  jav32  
         write(48,'(a)')  jav33  
         write(48,'(a)')  jav34  
         write(48,'(a)')  jav35  
         write(48,'(a)')  jav36  
         write(48,'(a)')  jav37  
         write(48,'(a)')  jav38  
         write(48,'(a)')  jav39  
         write(48,'(a)')  jav40  
         write(48,'(a)')  jav41  
         write(48,'(a)')  jav42  
         write(48,'(a)')  jav43  
         write(48,'(a)')  jav44  
         write(48,'(a)')  jav45  
         write(48,'(a)')  jav46  
         write(48,'(a)')  jav47  
         write(48,'(a)')  jav48  
         write(48,'(a)')  jav49  
         write(48,'(a)')  jav50  
         write(48,'(a)')  jav51  
         write(48,'(a)')  jav52  
         write(48,'(a)')  jav53  
         write(48,'(a)')  jav54  
         write(48,'(a)')  jav55  
         write(48,'(a)')  jav56  
         write(48,'(a)')  jav57  
         write(48,'(a)')  jav58  
         write(48,'(a)')  jav59  
         write(48,'(a)')  jav60  
         write(48,'(a)')  jav61  
         write(48,'(a)')  jav62  
         write(48,'(a)')  jav63  
         write(48,'(a)')  jav64  
         write(48,'(a)')  jav65  
         write(48,'(a)')  jav66  
         write(48,'(a)')  jav67  
         write(48,'(a)')  jav68  
         write(48,'(a)')  jav69  
         write(48,'(a)')  jav69b 
         write(48,'(a)')  jav70  
         write(48,'(a)')  jav71  
         write(48,'(a)')  jav72  
         write(48,'(a)')  jav73  
         write(48,'(a)')  jav74  
         write(48,'(a)')  jav75  
         write(48,'(a)')  jav76  
         write(48,'(a)')  jav77  
         write(48,'(a)')  jav78  
         write(48,'(a)')  jav79  
         write(48,'(a)')  jav80  
         write(48,'(a)')  jav81  
         write(48,'(a)')  jav82  
         write(48,'(a)')  jav83  
         write(48,'(a)')  jav84  
         write(48,'(a)')  jav85  
         write(48,'(a)')  jav86  
         write(48,'(a)')  jav87  
         write(48,'(a)')  jav88  
         write(48,'(a)')  jav89  

         nfrq = nfrqs()
         do i=1,nfrq
            call zerstr(i,ani,4,0)
            call frqstr(i,movfil,mlen)
            write(48,'(a,a,a)') "function animate",ani,"() {"
            write(48,'(a,a,a,a,a,a,a,f8.1,a,i4,a)') 
     &          st1,movfil(1:mlen),st2,st3,st4,st5,st6,freq(i),
     &          qd//st7,i,st8
            write(48,'(a)') "}"
         end do
         write(48,'(a)') "</script>"
         write(48,'(a)') "</HEAD>"
         write(48,'(a)') "<TABLE border=5 cellspacing=5>"
         write(48,'(a)') "<TR>"
         write(48,'(a)') "<TD>"
         write(48,'(a,a,a)') 
     &    "<applet name="//qd//"jmola"//qd//" code="//qd//"JmolApplet"
     &    //qd," archive="//qd//"/jmolnew/JmolApplet0.jar"//qd//
     &    " width="//qd//"350"//qd//" height="//qd//"500"//qd,
     &    " mayScript="//qd//"true"//qd//">"
         write(48,'(a)') "<param name="//qd//"style"//qd//
     &    " value="//qd//"shaded"//qd//">"
         write(48,'(a)') "<param name="//qd//"progressbar"//qd//
     &    " value="//qd//"true"//qd//">"
         write(48,'(a)') "<param name="//qd//"frank"//qd//
     &    " value="//qd//"no"//qd//">"
         write(48,'(a)') "<param name="//qd//"bgcolor"//qd//
     &    " value="//qd//"orange"//qd//">"
         write(48,'(a)') "<param name="//qd//"script"//qd//
     &    " value="//qd//"load mol.xyz;"//qd//">"
         write(48,'(a)') "</applet>"
         write(48,'(a)') "<TD width=1000>"
         write(48,'(a)') "<canvas id="//qd//"canvas"//qd//
     &    " width=1000 height=500 mouseEnabled="//qd//"true"//qd//
     &    "></canvas>"
         write(48,'(a)') "<!--[if !IE]> <img src="//qd//
     &    " width=1000 height=500 usemap="//qd//"#sel"//qd,
     &    " name="//qd//"irspec"//qd//"> <![endif]-->"

         write(48,'(a)') "<map name="//qd//"sel"//qd//">"
         do i=1,nfrq
            if (frint(i).ne.0.0d0) then
               call zerstr(i,ani,4,0)
               call gtfcor(i,ix1,iy1,ix2,iy2)
               write(48,'(a,i4,a,i4,a,i4,a,i4,a,i5,a,a4,a)') 
     &         "<area shape="//qd//"rect"//qd//" coords="//qd,
     &         ix1,",",iy1,",",ix2,",",iy2,qd//
     &         " id=",i," href="//qd
     &         //"javascript:animate",ani,"()"//qd//">"
            endif
         end do
         write(48,'(a)') "</map>"
         write(48,'(a)') "</TABLE>"

         write(48,'(a)') "<TABLE border=5>"
         write(48,'(a)') "<TR>"
         write(48,'(a)') "<TD> "
         write(48,'(a)') "<TD>Frequency"
         write(48,'(a)') "<TD>Intensity"
         do i=1,nfrq
            if (frint(i).ge.0.0d0) then
            call zerstr(i,ani,4,0)
            write(48,'(a)') "<TR>"
            write(48,'(a)') 
     &      "<td nowrap> <input type=radio name=frq value="//ani//
     &      " onclick="//qd//"javascript:animate"//ani//"()"//qd//">"
            write(48,'(a,f8.1,a)') "<td> ",freq(i),"</td>"
            write(48,'(a,f8.3,a)') "<td> ",frint(i),"</td>"
            end if
         end do
         write(48,'(a)') "</TABLE>"

      endif

      return
      end

