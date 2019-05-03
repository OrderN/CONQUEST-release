$! 
$! If you don't have Xwindows remove comment on line:
$!
$! fortran/contin=90 dummys.f
$!
$! And add comment on line:
$!
$! cc/define=("VMS","DOBACK") xwin.c
$!
$! in routine fmt.f change tol = 10.0d-66 to tol = 10.0d-39
$!
$ fortran/contin=90 atomdens.f
$ fortran/contin=90 molden.f
$ fortran/contin=90 above.f
$ fortran/contin=90 actcal.f
$ fortran/contin=90 basprt.f
$ fortran/contin=90 calc.f
$ fortran/contin=90 caldis.f
$ fortran/contin=90 calfa.f
$ fortran/contin=90 cntour.f
$ fortran/contin=90 cnvgam.f
$ fortran/contin=90 cnvgau.f
$ fortran/contin=90 convzmat.f
$ fortran/contin=90 cross.f
$ fortran/contin=90 crprod.f
$ fortran/contin=90 datin.f
$ fortran/contin=90 defpc.f
$ fortran/contin=90 defrad.f
$ fortran/contin=90 del.f
$ fortran/contin=90 denmak.f
$ fortran/contin=90 densmat.f
$ fortran/contin=90 distot.f
$ fortran/contin=90 dmat.f
$ fortran/contin=90 docent.f
$ fortran/contin=90 draw.f
$ fortran/contin=90 euler.f
$ fortran/contin=90 eulerh.f
$ fortran/contin=90 files.f
$ fortran/contin=90 fndcal.f
$ fortran/contin=90 gampoi.f
$ fortran/contin=90 gaupoi.f
$ fortran/contin=90 gaussian.f
$ fortran/contin=90 geogam.f
$ fortran/contin=90 geogau.f
$ fortran/contin=90 getmul.f
$ fortran/contin=90 getpoi.f
$ fortran/contin=90 getreal.f
$ fortran/contin=90 getarg.f
$ fortran/contin=90 gmmcnv.f
$ fortran/contin=90 grdcal.f
$ fortran/contin=90 gstr.f
$ fortran/contin=90 hidedr.f
$ fortran/contin=90 impsc.f
$ fortran/contin=90 locatc.f
$ fortran/contin=90 maxmin.f
$ fortran/contin=90 mdout.f
$ fortran/contin=90 mmcnv.f
$ fortran/contin=90 mopaco.f
$ fortran/contin=90 mopin.f
$ fortran/contin=90 mulprt.f
$ fortran/contin=90 obin.f
$ fortran/contin=90 occin.f
$ fortran/contin=90 oriin.f
$ fortran/contin=90 parang.f
$ fortran/contin=90 pareul.f
$ fortran/contin=90 parfc.f
$ fortran/contin=90 parori.f
$ fortran/contin=90 parpla.f
$ fortran/contin=90 parstp.f
$ fortran/contin=90 planky.f
$ fortran/contin=90 plend.f
$ fortran/contin=90 plini.f
$ fortran/contin=90 plmol.f
$ fortran/contin=90 plmolp.f
$ fortran/contin=90 plotgh.f
$ fortran/contin=90 plotgr.f
$ fortran/contin=90 plotin.f
$ fortran/contin=90 plpost.f
$ fortran/contin=90 pred.f
$ fortran/contin=90 prev.f
$ fortran/contin=90 proato.f
$ fortran/contin=90 procnv.f
$ fortran/contin=90 progeo.f
$ fortran/contin=90 rdbas.f
$ fortran/contin=90 rdcor.f
$ fortran/contin=90 rdgam.f
$ fortran/contin=90 rdgaus.f
$ fortran/contin=90 rdinfo.f
$ fortran/contin=90 rdpdb.f
$ fortran/contin=90 rdvect.f
$ fortran/contin=90 reada.f
$ fortran/contin=90 readel.f
$ fortran/contin=90 readvv.f
$ fortran/contin=90 renorm.f
$ fortran/contin=90 rmomen.f
$ fortran/contin=90 rota.f
$ fortran/contin=90 rotatg.f
$ fortran/contin=90 rotb.f
$ fortran/contin=90 rotc.f
$ fortran/contin=90 rotcor.f
$ fortran/contin=90 rotd.f
$ fortran/contin=90 rotfir.f
$ fortran/contin=90 rotm.f
$ fortran/contin=90 rotmom.f
$ fortran/contin=90 rott.f
$ fortran/contin=90 scback.f
$ fortran/contin=90 search.f
$ fortran/contin=90 searchd.f
$ fortran/contin=90 setang.f
$ fortran/contin=90 setbas.f
$ fortran/contin=90 settc.f
$ fortran/contin=90 shsort.f
$ fortran/contin=90 sillydum.f
$ fortran/contin=90 site.f
$ fortran/contin=90 slater.f
$ fortran/contin=90 stoc.f
$ fortran/contin=90 tessa.f
$ fortran/contin=90 tk4014.f
$ fortran/contin=90 tocap.f
$ fortran/contin=90 tocapf.f
$ fortran/contin=90 tomold.f
$ fortran/contin=90 under.f
$ fortran/contin=90 vaxconv.f
$ fortran/contin=90 vclr.f
$ fortran/contin=90 vec.f
$ fortran/contin=90 vlen.f
$ fortran/contin=90 vsc1.f
$ fortran/contin=90 wrinfo.f
$ fortran/contin=90 samino.f
$ fortran/contin=90 prtcal.f
$ fortran/contin=90 actss.f
$ fortran/contin=90 actami.f
$ fortran/contin=90 zread.f
$ fortran/contin=90 plden.f
$ fortran/contin=90 heaps.f
$ fortran/contin=90 den3d.f
$ fortran/contin=90 dencnt.f
$ fortran/contin=90 plhead.f
$ fortran/contin=90 pltab.f
$ fortran/contin=90 eucmol.f
$ fortran/contin=90 pl3dm.f
$ fortran/contin=90 plbox.f
$ fortran/contin=90 selsol.f
$ fortran/contin=90 atmd.f
$ fortran/contin=90 dolift.f
$ fortran/contin=90 spaced.f
$ fortran/contin=90 snypnt.f
$ fortran/contin=90 eulstr.f
$ fortran/contin=90 calct.f
$ fortran/contin=90 coeffs.f
$ fortran/contin=90 epint.f
$ fortran/contin=90 espot.f
$ fortran/contin=90 fcij.f
$ fortran/contin=90 fmt.f
$ fortran/contin=90 genaos.f
$ fortran/contin=90 rys.f
$ fortran/contin=90 ryspol.f
$ fortran/contin=90 rysrot.f
$ fortran/contin=90 thrcen.f
$ fortran/contin=90 twocen.f
$ fortran/contin=90 ifblen.f
$ fortran/contin=90 wrzmat.f
$ fortran/contin=90 rdchx.f
$ fortran/contin=90 gargpl.f
$ fortran/contin=90 freqs.f
$ fortran/contin=90 getmop.f
$ fortran/contin=90 getzm.f
$ fortran/contin=90 brklin.f
$ fortran/contin=90 xyzcoo.f
$ fortran/contin=90 geomop.f
$ fortran/contin=90 dumzm.f
$ fortran/contin=90 getxyz.f
$ fortran/contin=90 inferr.f
$ fortran/contin=90 espchrg.f
$ fortran/contin=90 proxim.f
$ fortran/contin=90 rdgamu.f
$ fortran/contin=90 plvrml.f
$ fortran/contin=90 molsint.f
$ fortran/contin=90 rdmsf.f
$ fortran/contin=90 wrmsf.f
$ fortran/contin=90 rdmolf.f
$ fortran/contin=90 adf_fun.f
$ fortran/contin=90 mpdum.f
$ fortran/contin=90 rotpol.f
$ fortran/contin=90 extbas.f
$!
$! fortran/contin=90 dummys.f
$! 
$ if f$getsyi ("ARCH_NAME") .eqs. "Alpha"
$ then
$  OS = "ALPHA"
$ else
$  OS = "VAX"
$ endif
$ if OS .eqs. "ALPHA"
$ then
$    cc/define=("VMS","DOBACK")/stand=vaxc xwin.c
$    cc/stand=vaxc strcasecmp.c
$ else
$    cc/define=("VMS","DOBACK") xwin.c
$    cc strcasecmp.c
$ endif
$!
$ delete molden.olb;*
$ library/create/object molden.olb *.obj
$ 
$ define lnk$library sys$library:vaxcrtl
$!
$ link/notrace/exe=molden sys$input/opt
molden/include=molden/library
sys$share:decw$xlibshr/share
