c!!!!!! Start of the xwindow dummys !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Delete or put a "C" in the first colomn of the following lines
c     when the xwindow driver is to be installed 
c!!!!!!
      subroutine xwin(xx,yy,ind,str,nstr,inct,incp)
      implicit double precision (a-h,o-z)
      real xx,yy
      character str*100
      print*,'XWindow Driver not installed'
      stop
      end
      subroutine molstr(str,idum1,idum2)
      implicit double precision (a-h,o-z)
      character str*(*)
      stop
      end
      subroutine errzme(str,idum1,idum2,idum3)
      implicit double precision (a-h,o-z)
      character str*(*)
      stop
      end
      subroutine butset(idum1,idum2,idum3)
      implicit double precision (a-h,o-z)
      stop
      end
      subroutine drawseg(iseg,idum1,idum2)
      implicit double precision (a-h,o-z)
      integer*2 iseg(*)
      stop
      end
      subroutine drwseg(iseg,idum1,idum2)
      implicit double precision (a-h,o-z)
      integer*2 iseg(*)
      stop
      end
      subroutine dash(ion)
      implicit double precision (a-h,o-z)
      stop
      end
      subroutine setcol(ion)
      implicit double precision (a-h,o-z)
      stop
      end
      subroutine curs(ion)
      implicit double precision (a-h,o-z)
      return
      end
      subroutine messg(ion)
      implicit double precision (a-h,o-z)
      return
      end
      subroutine doexp
      stop
      end
      subroutine dlogo
      stop
      end
      subroutine unstip
      stop
      end
      subroutine ststip
      stop
      end
      subroutine drwgeo
      stop
      end
      subroutine upzme
      stop
      end
      subroutine qupd
      stop
      end
      subroutine chkmpi
      return
      end
      subroutine ogribb(i)
      integer i
      stop
      end
      subroutine stowat(i)
      integer i
      return
      end
      subroutine nmrcpl(idebug)
      integer i
      return
      end
      subroutine ogendd(i)
      integer i
      stop
      end
      subroutine crsco
      stop
      end
      subroutine upsco
      stop
      end
      subroutine ogreset
      return
      end
      subroutine bldlst
      return
      end
      subroutine dlystr
      return
      end
      subroutine cpmf
      stop
      end
      subroutine srfclr
      stop
      end
      subroutine initsrf
      return
      end
      subroutine asurf(i,j)
      integer i,j
      stop
      end
      subroutine sribcol(i)
      integer i
      stop
      end
      subroutine setcll
      return
      end
      subroutine ribpnt(i,j)
      integer i,j
      return
      end
      subroutine ogbegg(i,j,k,l,r,m,s)
      integer i,j,k,l
      double precision r
      character s*(*)
      stop
      end
      subroutine ogvrt(r1,r2,r3)
      double precision r1,r2,r3
      stop
      end
      subroutine ognrm(r1,r2,r3)
      double precision r1,r2,r3
      stop
      end
      subroutine ogcoll(r1,r2,r3)
      double precision r1,r2,r3
      stop
      end
      subroutine rqsrt(i1,r1,i2)
      double precision r1
      stop
      end
      subroutine getene(i,j)
      return
      end
      subroutine gettrr(i,j)
      return
      end
      subroutine gtfrm(i)
      return
      end
      subroutine gethei(ihigh)
      implicit double precision (a-h,o-z)
      stop
      end
      subroutine plsph(ixp,iyp,iatoms,ia)
      implicit double precision (a-h,o-z)
      stop
      end
      subroutine plsph3(ixp,iyp,iatoms,ia,i1,i2)
      implicit double precision (a-h,o-z)
      stop
      end
      subroutine plsel(ixp,iyp,ia)
      implicit double precision (a-h,o-z)
      stop
      end
      subroutine labnr(ia)
      implicit double precision (a-h,o-z)
      ia = 0
      return
      end
      subroutine cnvfrg(ia)
      implicit double precision (a-h,o-z)
      stop
      end
      subroutine plrodx(ixp,iyp,d1,id2,id3,irod,i1,ca,sa)
      implicit double precision (a-h,o-z)
      real ca,sa
      stop
      end
      subroutine plrod3(ixp,iyp,d1,id2,id3,irod,i1,i2,ca,sa)
      implicit double precision (a-h,o-z)
      real ca,sa
      stop
      end
      integer function islck(i)
      implicit double precision (a-h,o-z)
      islck = 0
      return
      end
      integer function mseed()
      implicit double precision (a-h,o-z)
      mseed = 0
      return
      end
      subroutine drwstr(ixp,iyp,str,nstr)
      implicit double precision (a-h,o-z)
      character*2 str
      stop
      end
      subroutine drwqstr(ixp,iyp,ian,q)
      implicit double precision (a-h,o-z)
      stop
      end
      subroutine drwgl(x1,y1,x2,y2,i1,i2,i3,i4)
      real x1,y1,x2,y2
      stop
      end
      subroutine drwpol(iarray,i1,i2,i3,i4)
      integer*2 iarray(8)
      stop
      end
      subroutine parsfn(str,idum1,idum2)
      implicit double precision (a-h,o-z)
      character str*(*)
      return
      end
      subroutine prslab(str,idum1,idum2)
      implicit double precision (a-h,o-z)
      character str*(*)
      return
      end
      subroutine parfns(str,idum1)
      implicit double precision (a-h,o-z)
      character str*(*)
      return
      end
      subroutine parpoi(
     &            nzm,nso,nio,nzo,ioropt,ifor,ixyz98,iopr,isymm,irc,
     &            imp2,icntp,msucc,ioni,ipart,mopopt,isbin,irtype,
     &            ipdbgro,ifav,ioxyz,iconv,ircus,nscnd,iscst,ialtyp)
      return
      end
      subroutine parptr(i1,x,y,i2)
      real x,y
      return
      end
      subroutine parstr(str,idum1)
      implicit double precision (a-h,o-z)
      character str*(*)
      return
      end
      subroutine confrm(iop,istat)
      stop
      end
      subroutine exstr(str,idum1,idum2)
      implicit double precision (a-h,o-z)
      character str*(*)
      return
      end
      subroutine sysstr(str,idum1)
      implicit double precision (a-h,o-z)
      character str*(*)
      return
      end
      subroutine readsq(iamin,angs,namin)
      implicit double precision (a-h,o-z)
      dimension iamin(*),angs(7,*)
      return
      end

      subroutine cwidth(idum1)
      implicit double precision (a-h,o-z)
      stop
      end

      subroutine sollin
      stop
      end

      subroutine tounx
      implicit double precision (a-h,o-z)
      idum = 0
      return
      end
c dummys for zmatrix mallocs

      subroutine chkmaz(istat)
      implicit double precision (a-h,o-z)
      parameter (maxat=1000)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /zmatin/ bl(maxat),alph(maxat),bet(maxat),
     &                ibl(maxat),ialph(maxat),ibet(maxat),
     &                ianzz(maxat),iz(4,maxat),imap(maxat)
      common /atcom/  coo(3,numat1),rzp(numat1),ianz(numat1),
     &                iaton(numat1),iatclr(numat1),iresid(numat1),
     &                ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /hring/  lring(numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /chgbck/ qat(numat1)

      call chkmzz(istat,qat,rzp,ianzz,imap,ianz,lring,ityp,ipdbt)

      return
      end

      subroutine ogmon()
      implicit double precision (a-h,o-z)
      print*,'XWindow Driver not installed'
      stop
      end

      subroutine fndmaz(ixyz,istat)
      implicit double precision (a-h,o-z)
      parameter (maxat=1000)
      common /zmatin/ bl(maxat),alph(maxat),bet(maxat),
     &                ibl(maxat),ialph(maxat),ibet(maxat),
     &                ianz(maxat),iz(4,maxat),imap(maxat)

      call fndmzz(ixyz,istat,imap)

      return
      end

      subroutine mapxzz(iun,imod,iff,izmtmp,istat)
      implicit double precision (a-h,o-z)
      parameter (maxat=1000)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /zmatin/ bl(maxat),alph(maxat),bet(maxat),
     &                ibl(maxat),ialph(maxat),ibet(maxat),
     &                iianz(maxat),iz(4,maxat),imap(maxat)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call mapzzz(iun,imod,iff,izmtmp,istat,bl,alph,bet,ibl,ialph,
     &            ibet,imap,iianz,iz,ianz)

      return
      end

      subroutine wrzmat(iun,iopt)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      parameter (maxat=1000)
      parameter (MAXPNT=2000)
      common /zmatin/ bl(maxat),alph(maxat),bet(maxat),
     &                ibl(maxat),ialph(maxat),ibet(maxat),
     &                ianzz(maxat),iz(4,maxat),imap(maxat)
      common /geocnv2/ formax(MAXPNT),forrms(MAXPNT),dismax(MAXPNT),
     &                 disrms(MAXPNT),epoints(MAXPNT),isav(MAXPNT)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call wrzmaz(iun,iopt,bl,alph,bet,ibl,ialph,
     &            ibet,imap,ianzz,iz,epoints,ianz)

      return
      end

      subroutine wline(iun,iopt,igamb)
      implicit double precision (a-h,o-z)
      parameter (maxat=1000)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /zmatin/ bl(maxat),alph(maxat),bet(maxat),
     &                ibl(maxat),ialph(maxat),ibet(maxat),
     &                ianzz(maxat),iz(4,maxat),imap(maxat)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /chgbck/ qat(numat1)

      call wlinz(iun,iopt,igamb,bl,alph,bet,ibl,ialph,
     &           ibet,imap,ianzz,iz,iconn,ianz,ityp,qat)

      return
      end

      subroutine dumlin(isel,blv,alphv,betv)
      implicit double precision (a-h,o-z)
      parameter (maxat=1000)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /zmatin/ bl(maxat),alph(maxat),bet(maxat),
     &                ibl(maxat),ialph(maxat),ibet(maxat),
     &                ianzz(maxat),iz(4,maxat),imap(maxat)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /zmsbck/ lwrit(numat1)
      dimension isel(*)

      call dumliz(isel,blv,alphv,betv,bl,alph,bet,ibl,ialph,
     &            ibet,imap,ianzz,iz,lwrit,ianz)

      return
      end

      subroutine plinzz(isel,istat)
      implicit double precision (a-h,o-z)
      parameter (maxat=1000)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /zmatin/ bl(maxat),alph(maxat),bet(maxat),
     &                ibl(maxat),ialph(maxat),ibet(maxat),
     &                ianzz(maxat),iz(4,maxat),imap(maxat)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /zmsbck/ lwrit(numat1)
      dimension isel(*)

      call plinz(isel,istat,bl,alph,bet,ibl,ialph,
     &            ibet,imap,ianzz,iz,lwrit,ianz)

      return
      end

      subroutine prtzm
      implicit double precision (a-h,o-z)
      parameter (maxat=1000)
      common /zmatin/ bl(maxat),alph(maxat),bet(maxat),
     &                ibl(maxat),ialph(maxat),ibet(maxat),
     &                ianz(maxat),iz(4,maxat),imap(maxat)

      call prtzz(bl,alph,bet,ibl,ialph,ibet,imap,ianz,iz)

      return
      end

      subroutine getzmz(istat)
      implicit double precision (a-h,o-z)
      parameter (maxat=1000)
      common /zmatin/ bl(maxat),alph(maxat),bet(maxat),
     &                ibl(maxat),ialph(maxat),ibet(maxat),
     &                ianz(maxat),iz(4,maxat),imap(maxat)

      call getzzz(istat,bl,alph,bet,ibl,ialph,ibet,imap,ianz,iz)

      return
      end

      subroutine dumzm(cc,ianc,nnatoms)
      implicit double precision (a-h,o-z)
      parameter (maxat=1000)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /zmatin/ bl(maxat),alph(maxat),bet(maxat),
     &                ibl(maxat),ialph(maxat),ibet(maxat),
     &                ianzz(maxat),iz(4,maxat),imap(maxat)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      dimension cc(3,*),ianc(*)
      dimension ian(maxat),c(3,maxat),cz(3,maxat),
     &          alpha(maxat),beta(maxat)

      call dumzz(cc,ianc,nnatoms,bl,alph,bet,ibl,ialph,ibet,
     &           imap,ianzz,iz,c,cz,alpha,beta,ian,coo,iresid,
     &           issdon)

      return
      end

      subroutine convzmat(cc,ianc,nnatoms,igo,ico,ido)
      implicit double precision (a-h,o-z)
      parameter (maxat=1000)
      common /zmatin/ bl(maxat),alph(maxat),bet(maxat),
     &                ibl(maxat),ialph(maxat),ibet(maxat),
     &                ianz(maxat),iz(4,maxat),imap(maxat)
      dimension cc(3,*),ianc(*)
      dimension ian(maxat),c(3,maxat),cz(3,maxat),
     &          alpha(maxat),beta(maxat)

      call convzmzz(cc,ianc,nnatoms,igo,ico,ido,
     &           bl,alph,bet,ibl,ialph,ibet,
     &           imap,ianz,iz,c,cz,alpha,beta,ian)

      return
      end

      subroutine getmop(nnatoms,heat,igo,ico,istat)
      implicit double precision (a-h,o-z)
      parameter (maxat=1000)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /zmatin/ bl(maxat),alph(maxat),bet(maxat),
     &                ibl(maxat),ialph(maxat),ibet(maxat),
     &                ianzz(maxat),iz(4,maxat),imap(maxat)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      dimension ian(maxat),c(3,maxat),cz(3,maxat),
     &          alpha(maxat),beta(maxat)

      call getmdp(coo,ianz,nnatoms,heat,igo,ico,istat,
     &           bl,alph,bet,ibl,ialph,ibet,
     &           imap,ianzz,iz,c,cz,alpha,beta,ian)

      return
      end

      subroutine getzm(nnatoms,igo,ico,istat)
      implicit double precision (a-h,o-z)
      parameter (maxat=1000)
      common /zmatin/ bl(maxat),alph(maxat),bet(maxat),
     &                ibl(maxat),ialph(maxat),ibet(maxat),
     &                ianzz(maxat),iz(4,maxat),imap(maxat)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      dimension ian(maxat),c(3,maxat),cz(3,maxat),
     &          alpha(maxat),beta(maxat)

      call getzd(coo,ianz,nnatoms,igo,ico,istat,
     &           bl,alph,bet,ibl,ialph,ibet,
     &           imap,ianzz,iz,c,cz,alpha,beta,ian)

      return
      end

      subroutine allorb(isiz,ifrst)
      implicit double precision (a-h,o-z), integer (i-n)

      idum = 0

      return
      end

      subroutine allgrd(isiz)
      implicit double precision (a-h,o-z), integer (i-n)

      idum = 0

      return
      end

      subroutine allgeo(isiz,ifrst)
      implicit double precision (a-h,o-z), integer (i-n)

      idum = 0

      return
      end

      subroutine allcoo(isiz,ifrst)
      implicit double precision (a-h,o-z), integer (i-n)

      idum = 0

      return
      end

      subroutine almgrd()
      implicit double precision (a-h,o-z), integer (i-n)

      idum = 0

      return
      end

      subroutine dozmt(istat)
      implicit double precision (a-h,o-z), integer (i-n)

      istat = 0

      return
      end

      subroutine adffun(x,y,z,psi)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      real stoalfa,stobnorm
      common /ADF/    istos(5,maxorb),naorbs,stoalfa(maxorb),
     &                stobnorm(maxorb)
      dimension psi(*)

      call adffud(x,y,z,psi,
     &            stoalfa,stobnorm,istos,naorbs)
      return
      end

      subroutine atdd(imo,ipsi,iao)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /densty/ p(maxorb*maxorb), paa(maxorb)
      common /dynhlp/ qd(maxorb*maxorb),pd(maxorb),gd(3,maxorb),
     &                hd(6,maxorb)

      call addd(imo,ipsi,iao,p,paa,pd)

      return
      end

      subroutine denmad(ido)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      common /densty/ p(maxorb*maxorb), paa(maxorb)
      common /sphere/ averag(maxorb)
      common /dynhlp/ qd(maxorb*maxorb),pd(maxorb),gd(3,maxorb),
     &                hd(6,maxorb)


      call denmdd(ido,vectrs,vectrb,focc,focb,p,paa,averag,qd)

      return
      end

      subroutine densmat(ido)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      common /densty/ p(maxorb*maxorb), paa(maxorb)

      call densmad(ido,vectrs,vectrb,focc,focb,p,nocc)

      return
      end

      subroutine densmto(ido)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /densty/ p(maxorb*maxorb), paa(maxorb)

      call densmtd(ido,vectrs,focc,p,nocc)

      return
      end

      subroutine espot(x,y,z,epot,idebug)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /densty/ p(maxorb*maxorb), paa(maxorb)

      call espod(x,y,z,epot,idebug,p)

      return
      end

      subroutine espgrd(npts1,npts2,npts3,idebug)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      parameter (max3d=61)
      common /spa3d/  denn(max3d*max3d*max3d),pmnn(max3d)
      common /densty/ p(maxorb*maxorb), paa(maxorb)

c arguments missing, xden,yden,zden !!!
c      call espgrdd(npts1,npts2,npts3,idebug,denn,p)

      return
      end

      subroutine boys(norb,nor,dmao,vecs)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension nor(*),dmao(*),vecs(*)

      print*,"Boys not available in noxwin"

      return
      end

      subroutine newdenmak
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /densty/ p(maxorb*maxorb), paa(maxorb)
      common /dynhlp/ qd(maxorb*maxorb),pd(maxorb),gd(3,maxorb),
     &                hd(6,maxorb)

      call newdenmad(p,qd)

      return
      end

      subroutine calhes(psi,grd,hess,den,g,hes)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /densty/ p(maxorb*maxorb), paa(maxorb)
      dimension psi(*),grd(3,*),hess(6,*), g(3),hes(6)

      call calhed(psi,grd,hess,den,g,hes,p)

      return
      end

      subroutine grdcal(dens,npts1,npts2,iprnt,space)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /densty/ p(maxorb*maxorb), paa(maxorb)
      common /dynhlp/ qd(maxorb*maxorb),pd(maxorb),gd(3,maxorb),
     &                hd(6,maxorb)
      integer space
      dimension dens(*)

      call grdcad(dens,npts1,npts2,iprnt,space,
     &            p,paa,pd,gd,hd)

      return
      end

      subroutine prtvec
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      real eiga,eigb
      common /eigval/ eiga(maxorb),eigb(maxorb)

      call prtved(vectrs,vectrb,eiga,eigb,ncols,ncolb)

      return
      end

      subroutine stint
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      parameter (numat1=20000)
      common /densty/ p(maxorb*maxorb), paa(maxorb)
      common /chgbck/ qat(numat1)

      call stind(p,qat)

      return
      end

      subroutine muldmd(vdwr,moddma,idm,idd)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /densty/ p(maxorb*maxorb), paa(maxorb)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      dimension vdwr(*)

      call mulddd(vdwr,moddma,idm,idd,p,qat,iconn)

      return
      end

      subroutine mopin(istat,ibin,impas)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /densty/ p(maxorb*maxorb), paa(maxorb)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      real eiga,eigb
      common /eigval/ eiga(maxorb),eigb(maxorb)
      common /dynhlp/ qd(maxorb*maxorb),pd(maxorb),gd(3,maxorb),
     &                hd(6,maxorb)
      common /sphere/ averag(maxorb)

      call mopdd(istat,ibin,impas,vectrs,averag,p,focc,eiga,qd,pd,
     &           nocc,ncols)

      return
      end

      subroutine rdgad(idebug,ibefo,istatio,ioxyz,irtype,istats)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      real eiga,eigb
      common /eigval/ eiga(maxorb),eigb(maxorb)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call rdgdd(idebug,ibefo,istatio,ioxyz,irtype,istats,
     &           vectrs,vectrb,focc,focb,eiga,eigb,nocc,nocb,ncols,
     &           ncolb,coo,ianz)

      return
      end

      subroutine rdgamd(idebug,ibefo,istatio,ioxyz,irtype,
     &                  ihsend,istats)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      real eiga,eigb
      common /eigval/ eiga(maxorb),eigb(maxorb)

      call rdgadd(idebug,ibefo,istatio,ioxyz,irtype,ihsend,istats,
     &            vectrs,vectrb,focc,focb,eiga,eigb,
     &            nocc,nocb,ncols,ncolb)

      return
      end

      subroutine rdcpmdd(idebug,ibefo,istatio,ioxyz,irtype,
     &                  ihsend,istats)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      real eiga,eigb
      common /eigval/ eiga(maxorb),eigb(maxorb)

      call rdcpmddd(idebug,ibefo,istatio,ioxyz,irtype,ihsend,istats,
     &     vectrs,vectrb,focc,focb,eiga,eigb,
     &     nocc,nocb,ncols,ncolb)
 
      return
      end

      subroutine rdgaud(idebug,ibefo,istatio,irtype,istats)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      call rdgdud(idebug,ibefo,istatio,irtype,istats,focc,focb,
     &            nocc,nocb,ncols,ncolb,coo,ianz)

      return
      end

      subroutine rdmaud(idebug,istatio,istats)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      parameter (MAXPNT=2000)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      real eiga,eigb
      common /eigval/ eiga(maxorb),eigb(maxorb)
      common /geocnv2/ formax(MAXPNT),forrms(MAXPNT),dismax(MAXPNT),
     &                 disrms(MAXPNT),epoints(MAXPNT),isav(MAXPNT)
      common /dynhlp/ qd(maxorb*maxorb),pd(maxorb),gd(3,maxorb),
     &                hd(6,maxorb)
      common /sphere/ averag(maxorb)

      call rdmadd(idebug,istatio,istats,
     &            vectrs,vectrb,focc,focb,eiga,eigb,
     &            averag,p,qd,pd,nocc,nocb,ncols,ncolb,
     &            formax,forrms,dismax,disrms,epoints,isav,qat)

      return
      end

      subroutine rdmold(idebug,istatio,irtype,iesp,istats)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      parameter (MAXPNT=2000)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      real eiga,eigb
      common /eigval/ eiga(maxorb),eigb(maxorb)
      real stoalfa,stobnorm
      common /ADF/    istos(5,maxorb),naorbs,stoalfa(maxorb),
     &                stobnorm(maxorb)
      common /geocnv2/ formax(MAXPNT),forrms(MAXPNT),dismax(MAXPNT),
     &                 disrms(MAXPNT),epoints(MAXPNT),isav(MAXPNT)

      call rdmodd(idebug,istatio,irtype,iesp,istats,
     &            vectrs,vectrb,focc,focb,eiga,eigb,nocc,nocb,
     &            ncols,ncolb,stoalfa,stobnorm,istos,naorbs,
     &            formax,forrms,dismax,disrms,epoints,isav)

      return
      end

      subroutine prtmolf(iun,ihaszm,ipoints)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      parameter (numat1=20000)
      parameter (mxcon=10)
      parameter (MAXPNT=2000)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      real eiga,eigb
      common /eigval/ eiga(maxorb),eigb(maxorb)
      real stoalfa,stobnorm
      common /ADF/    istos(5,maxorb),naorbs,stoalfa(maxorb),
     &                stobnorm(maxorb)
      common /geocnv2/ formax(MAXPNT),forrms(MAXPNT),dismax(MAXPNT),
     &                 disrms(MAXPNT),epoints(MAXPNT),isav(MAXPNT)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call prtmold(iun,ihaszm,ipoints,
     &            vectrs,vectrb,focc,focb,eiga,eigb,nocc,nocb,
     &            ncols,ncolb,stoalfa,stobnorm,istos,naorbs,
     &            formax,forrms,dismax,disrms,epoints,isav,ianz)

      return
      end

      subroutine rdvect(idebug,ig94,istats)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      real eiga,eigb
      common /eigval/ eiga(maxorb),eigb(maxorb)

      call rdvecd(idebug,ig94,istats,vectrs,vectrb,focc,eiga,eigb,
     &            ncols,ncolb)

      return
      end

      subroutine datin(npts1,npts2,npts3)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      real eiga,eigb
      common /eigval/ eiga(maxorb),eigb(maxorb)


      call datid(npts1,npts2,npts3,ncols,ncolb,focc,focb,eiga,eigb)

      return
      end

      subroutine homo(ipsi)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      real eiga,eigb
      common /eigval/ eiga(maxorb),eigb(maxorb)

      call homod(ipsi,focc,focb,eiga,eigb)

      return
      end

      subroutine lumo(ipsi)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      real eiga,eigb
      common /eigval/ eiga(maxorb),eigb(maxorb)

      call lumod(ipsi,focc,focb,eiga,eigb)

      return
      end

      subroutine occup(istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb

      call occupd(istat,focc,focb)

      return
      end

      subroutine cntour(a,mdim,imax,jmax,pz,value,r11,id)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)
      dimension a(mdim,*), id(mdim,*)

      call cntoud(a,mdim,imax,jmax,pz,value,r11,id,ix)

      return
      end

      subroutine den3d(npts1,npts2,scale)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)

      call dendd(npts1,npts2,scale,dens,edx,edy,iedlog)

      return
      end

      subroutine dencnt(npts1,npts2,fcnt)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)
      dimension fcnt(*)

      call dencnd(npts1,npts2,fcnt,dens,iedlog)

      return
      end

      subroutine grdcl(npts1,npts2,iprnt,ispace)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)

      call grdcd(npts1,npts2,iprnt,ispace,dens)

      return
      end

      subroutine rdqchm(idebug,irtype,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb

      call rdqchd(idebug,irtype,istat,focc,focb,nocc,nocb,ncols,ncolb)

      return
      end

      subroutine rdqvec(idebug,istats)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      real eiga,eigb
      common /eigval/ eiga(maxorb),eigb(maxorb)

      call rdqvcd(idebug,istats,
     &                  vectrs,vectrb,eiga,eigb,ncols,ncolb)

      return
      end

      subroutine rdorca(idebug,irtype,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb

      call rdorcd(idebug,irtype,istat,ianz,iatoms,
     &            focc,focb,nocc,nocb,ncols,ncolb)

      return
      end

      subroutine pltspec(inum)
      implicit double precision (a-h,o-z), integer (i-n)

      return
      end

      subroutine prsogl
      implicit double precision (a-h,o-z), integer (i-n)

      return
      end

      subroutine rdnwch(idebug,irtype,istats)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb

      call rdnwcd(idebug,irtype,istats,
     &                 ianz,iatoms,focc,focb,nocc,nocb,ncols,ncolb)

      return
      end

      subroutine nmcshl
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (maxcpl=100)
      dimension couplj(maxcpl*maxcpl)

      call nmcshd(couplj)

      return
      end

      subroutine resedl
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)

      call resedd(iedlog)

      return
      end

      subroutine maxmin(npts1,npts2,scale)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)

      call maxmid(npts1,npts2,scale,dens)

      return
      end

      subroutine plden(ndim1,ndim2,scale,icells,adjus,idisml)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)

      call plded(ndim1,ndim2,scale,icells,adjus,idisml,
     &           dens,ix,iy,rz)

      return
      end

      subroutine p3dv(iun,scale,ndimx,ndimz,adjus)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)

      call p3dd(iun,scale,ndimx,ndimz,adjus,dens)

      return
      end

      subroutine rdcube(npts1,npts2,npts3,iposng,ipsi,istat,iun,idebug)
      implicit double precision (a-h,o-z)
      parameter (max3d=61)
      common /spa3d/  denn(max3d*max3d*max3d),pmnn(max3d)

      call rdcubd(npts1,npts2,npts3,iposng,ipsi,istat,iun,idebug,
     &            denn,pmnn)

      return
      end

      subroutine wrcube(npts1,npts2,npts3,ipsi)
      implicit double precision (a-h,o-z)
      parameter (max3d=61)
      common /spa3d/  denn(max3d*max3d*max3d),pmnn(max3d)

      call wrcubd(npts1,npts2,npts3,ipsi,denn)

      return
      end

      subroutine rdgrd(npts1,npts2,npts3,iun,istat)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /spa3d/  denn(max3d*max3d*max3d),pmnn(max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma
      character*4 buff(max3d2)

      call rdgrdd(npts1,npts2,npts3,iun,istat,
     &            denn,dens,pmnn,buff,ichx)

      return
      end

      subroutine dpomap(iopt)
      implicit double precision (a-h,o-z)
      parameter (max3d=61)
      common /spa3d/  denn(max3d*max3d*max3d),pmnn(max3d)

      call dpomad(iopt,denn)

      return
      end

      subroutine rdomap(npts1,npts2,npts3,iun,istat)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /spa3d/  denn(max3d*max3d*max3d),pmnn(max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call rdomad(npts1,npts2,npts3,iun,istat,
     &            denn,dens,pmnn,ichx)

      return
      end

      subroutine rdvasp(npts1,npts2,npts3,iposng,istat,lenf,idocub,
     &                  idebug)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /spa3d/  denn(max3d*max3d*max3d),pmnn(max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call rdvasd(npts1,npts2,npts3,iposng,istat,lenf,idocub,
     &            idebug,denn,pmnn,bucket,coo,ianz,iatclr,iconn,
     &            nat,norg,icent,inorm,ncon,nspg,ichx,
     &            nopr,ir,it,
     &            xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      
      return
      end
      
      subroutine rdconquest(npts1,npts2,npts3,iposng,istat,lenf,idocub,
     &                  idebug)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /spa3d/  denn(max3d*max3d*max3d),pmnn(max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma
      
      call rdcqd(npts1,npts2,npts3,iposng,istat,lenf,idocub,
     &            idebug,denn,pmnn,bucket,coo,ianz,iatclr,iconn,
     &            nat,norg,icent,inorm,ncon,nspg,ichx,
     &            nopr,ir,it,
     &            xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)

      return
      end
      
      subroutine rdinfo(npts1,npts2,isubtr,istat)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)

      call rdinfd(npts1,npts2,isubtr,istat,adjus,dens)

      return
      end

      subroutine rd3inf(npts1,npts2,npts3,isubtr,adjus,istat)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /spa3d/  denn(max3d*max3d*max3d),pmnn(max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)

      call rd3ind(npts1,npts2,npts3,isubtr,adjus,istat,
     &            denn,pmnn,denst)

      return
      end

      subroutine rd3chk(npts1,npts2,npts3,igauss,impas,istat)
      implicit double precision (a-h,o-z)
      parameter (maxm3d=10)
      common /map3d/  fmap(maxm3d*maxm3d*maxm3d),mxm3d

      call rd3chd(npts1,npts2,npts3,igauss,impas,istat,fmap)

      return
      end

      subroutine wrinfo(npts1,npts2)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)

      call wrinfd(npts1,npts2,dens)

      return
      end

      subroutine wr3inf(npts1,npts2,npts3,adjus)
      implicit double precision (a-h,o-z)
      parameter (max3d=61)
      common /spa3d/  denn(max3d*max3d*max3d),pmnn(max3d)

      call wr3ind(npts1,npts2,npts3,adjus,denn)

      return
      end

      subroutine denfst(sum,psi)
      implicit double precision (a-h,o-z)
      parameter (maxorb=256)
      common /vectrs/ vectrs(maxorb*maxorb),focc(maxorb),nocc,ncols
      common /vectrb/ vectrb(maxorb*maxorb),focb(maxorb),nocb,ncolb
      dimension psi(*)

      call denfsd(sum,psi,vectrs,vectrb,focc,focb)

      return
      end

      subroutine spaced(npts1,npts2,npts3,valcnt,idofil,adjus,
     &                  ipsprt,idisml,idvrml,mapit)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      parameter (maxm3d=10)

      common /spa3d/  denn(max3d*max3d*max3d),pmnn(max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)
      common /map3d/  fmap(maxm3d*maxm3d*maxm3d),mxm3d

      call spacdd(npts1,npts2,npts3,valcnt,idofil,adjus,ipsprt,
     &            idisml,idvrml,mapit,denn,pmnn,iedlog,fmap)

      return
      end

      subroutine spasrf(npts1,npts2,npts3,valcnt)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /spa3d/  denn(max3d*max3d*max3d),pmnn(max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)

      call spasrd(npts1,npts2,npts3,valcnt,denn,pmnn,dens,iedlog)

      return
      end

      subroutine isoden(valc,nvalc,scincr,nespt,iwhere)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)

      call isoded(valc,nvalc,scincr,nespt,iwhere,dens,iedlog)

      return
      end

      subroutine wrcart(iun,dopdb,idogau,ipdbwh)
      parameter (MAXPNT=2000)
      implicit double precision ( a-h,o-z )
      common /geocnv2/ formax(MAXPNT),forrms(MAXPNT),dismax(MAXPNT),
     &                 disrms(MAXPNT),epoints(MAXPNT),isav(MAXPNT)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call wrcard(iun,dopdb,idogau,ipdbwh,epoints,
     &            ncalf,ianf,islu,nchain,iamino,reson,irsnr,achain)

      return
      end

      subroutine plvend(iun,iloop)
      parameter (MAXPNT=2000)
      implicit double precision ( a-h,o-z )
      common /geocnv2/ formax(MAXPNT),forrms(MAXPNT),dismax(MAXPNT),
     &                 disrms(MAXPNT),epoints(MAXPNT),isav(MAXPNT)
      common /sclcom/ scal,fscal,scali,smag,iscupd

      call plvedd(iun,iloop,epoints,scal)

      return
      end

      subroutine dyncpmd(ipoints)
      parameter (MAXPNT=2000)
      implicit double precision ( a-h,o-z )
      common /geocnv2/ formax(MAXPNT),forrms(MAXPNT),dismax(MAXPNT),
     &                 disrms(MAXPNT),epoints(MAXPNT),isav(MAXPNT)

      call dyncpdd(ipoints,isav,epoints)

      return
      end

      subroutine progeo(ipoints,iff,istat)
      parameter (MAXPNT=2000)
      implicit double precision ( a-h,o-z )
      common /geocnv2/ formax(MAXPNT),forrms(MAXPNT),dismax(MAXPNT),
     &                 disrms(MAXPNT),epoints(MAXPNT),isav(MAXPNT)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call proged(ipoints,iff,istat,
     &            formax,forrms,dismax,disrms,epoints,isav,ichx)

      return
      end

      subroutine oginit(r,adjus,natoms,nat,icol,xsym,ysym,zsym,
     &                  vdwr,cnst,npts1,npts2,iorb)
      implicit double precision (a-h,o-z)
      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (max3d2=max3d*max3d)
      common /grdall/ dens(max3d2),denst(max3d2),edx(max3d2),
     &                edy(max3d2),rz(max3d2),bucket(max3d+lnbuck),
     &                iedlog(max3d2),ix(max3d2),iy(max3d2)
      dimension r(*),xsym(*),ysym(*),zsym(*),vdwr(*)
      dimension nat(*),icol(*)

      call oginid(r,adjus,natoms,nat,icol,xsym,ysym,zsym,
     &            vdwr,cnst,npts1,npts2,iorb,dens)

      return
      end

      subroutine anim
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call anid(coo,ianz)

      return
      end

      subroutine actami(inum,ikleur,iopt,idosrf)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /srfbck/ isurf(numat1)

      call actamd(inum,ikleur,iopt,idosrf,iaton,iatclr,iresid,isurf,
     &            icalf,ncalf,ianf,islu,nchain,iamino,ihet,
     &            iclhet,reson,iams,isal)

      return
      end

      subroutine actexp(que,lque,ikleur,idosrf)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)

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
      call actexd(que,lque,ikleur,idosrf,
     &            ncalf,iamino)

      return
      end

      subroutine aacom(vrad,ires,str,nstr,istsrf)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /srfbck/ isurf(numat1)
      dimension vrad(*)

      call aacod(vrad,ires,str,nstr,istsrf,iresid,isurf,
     &           ncalf,ianf,islu,nchain,iamino,ihet,reson,
     &           isal,irsnr,achain);

      return
      end

      subroutine actcal(iopt)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      call actcad(iopt,iaton,iatclr,iresid,icalf,ianf,islu,nchain,
     &            iamino,ibck);

      return
      end

      subroutine acthel(iopt,iscnd,jcolsp,inclbb)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      call acthed(iopt,iscnd,jcolsp,inclbb,iaton,iatclr,iresid,iconn,
     &            icalf,ianf,islu,nchain,iamino,ihet,reson,
     &            isal,ibck)

      return
      end

      subroutine actss(iopt)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call actsd(iopt,ianz,iaton,iatclr,iresid,iconn)

      return
      end

      subroutine pmfass(iopt,dochg)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      integer dochg
      common /typbck/ ityp(numat1),ipdbt(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      call pmfasd(iopt,dochg,ianz,iresid,iconn,ityp,ncalf,iamino)

      return
      end

      subroutine dfiass
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call dfiasd(ityp,ncalf,iamino)

      return
      end

      subroutine updres
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
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

      call updred(coo,ncalf,icalf)

      return
      end

      subroutine dfires
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call dfired(coo,icalf,ncalf,iamino)

      return
      end

      subroutine twodfi(tdfi,ires1,ires2)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)

      call twodfd(tdfi,ires1,ires2,coo,ityp)

      return
      end

      subroutine twodfib(tdfi,ires1,ires2)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)

      call twodfid(tdfi,ires1,ires2,coo,ityp)

      return
      end

      subroutine onedfi(tdfi,ires1)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)

      call onedfd(tdfi,ires1,coo,ityp)

      return
      end

      subroutine totpmf(ttpmf)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)

      call totpmd(ttpmf,coo,ianz,iaton,iatclr,ityp)

      return
      end

      subroutine pmfinf(iatm)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)

      call pmfind(iatm,coo,ianz,ityp)

      return
      end

      subroutine ipmtyp(iptm,iat,ian,idochg)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)

      call ipmtyd(iptm,iat,ian,idochg,ianz,iconn,qat)

      return
      end

      subroutine clmon(car,pot)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      dimension car(*)

      call clmod(car,pot,coo,qat)

      return
      end

      subroutine clmons(car,pot,idoloc)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      common /condens/ icont(numat1),ncont
      dimension car(*)

      call clmond(car,pot,idoloc,coo,qat,icont,ncont)

      return
      end

      subroutine calfa(istat,istpdb,iaddh,ioatms,nstrt,ioadd)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call calfd(istat,istpdb,iaddh,ioatms,nstrt,ioadd,ianz,iconn,ityp,
     &           icalf,ncalf,ianf,islu,nchain,iamino,
     &           isal,irsnr,ihashb)

      return
      end

      subroutine docct
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call doccd(coo,ianz)

      return
      end

      subroutine docent
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /rotmat/rx(3),ry(3),rz(3),t(3),td(3)

      call docend(t,coo,ianz)

      return
      end

      subroutine valdis(var,ipntr,numat,iop,iasel,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (mxeat=300)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      dimension var(mxeat,2),ipntr(mxeat)

      call valdid(var,ipntr,numat,iop,iasel,istat,ianz,iresid,qat)

      return
      end

      subroutine eemcalc(var,ipntr,numat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (mxeat=300)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      dimension var(mxeat,2),ipntr(mxeat)

      call eemcald(var,ipntr,numat,ianz,coo,qat)

      return
      end

      subroutine distchk
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call distchd(coo)

      return
      end

      subroutine calgas(iasel,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)

      call calgad(iasel,istat,qat,ianz,iconn,iresid,ityp)

      return
      end

      subroutine mkscon(c,rad,grdw,iptr,npts)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      dimension c(*)

      call mkscod(c,rad,grdw,iptr,npts,coo,iconn)

      return
      end

      subroutine connlp(dens,idomap,isp)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /srfbck/ isurf(numat1)

      call connld(dens,idomap,isp,coo,ianz,iaton,iatclr,iresid,iconn,
     &            isurf)

      return
      end

      subroutine allsrf(idocol,idomap,idocal,isp)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /srfbck/ isurf(numat1)

      call allsrd(idocol,idomap,idocal,isp,isurf,iaton,ianz,iatclr)

      return
      end

      subroutine alasrf
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /srfbck/ isurf(numat1)
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

      call alasrd(iresid,isurf,iams,ihets)

      return
      end

      subroutine clrsrf
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      common /srfbck/ isurf(numat1)

      call clrsrd(isurf)

      return
      end

      subroutine propnt(qx,qy,qz,ipen)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call propnd(qx,qy,qz,ipen,
     &            coo,ianz,iaton,iatclr,iresid,iconn)

      return
      end

      subroutine getfr(istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call getfd(istat,coo)

      return
      end

      subroutine ggetfr(istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call ggetfd(istat,coo)

      return
      end

      subroutine getfrm(istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call getfrd(istat,coo,ianz)

      return
      end

      subroutine getfra(istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call getfrad(istat,ianz)

      return
      end

      subroutine nxtpnt(iopt,fancy,atcol,dolabs,backb,wpnt)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer dolabs,fancy,atcol,backb
      logical wpnt

      call nxtpnd(iopt,fancy,atcol,dolabs,backb,wpnt,coo)

      return
      end

      subroutine scalfr(iopt,fancy,atcol,dolabs,backb,wpnt)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer dolabs,fancy,shade,atcol,backb
      logical wpnt
      common /rotmat/rx(3),ry(3),rz(3),t(3),td(3)
      common /sclcom/ scal,fscal,scali,smag,iscupd

      call scalfd(iopt,fancy,atcol,dolabs,backb,wpnt,t,coo,
     &            scal,scali,smag)

      return
      end

      subroutine resfr
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call resfd(coo)

      return
      end

      subroutine tofcoo
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call tofcod(coo)

      return
      end

      subroutine iatnox
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call iatnod(ianz)

      return
      end

      subroutine gampoi(ipoint,istat,ioxyz)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call gampod(ipoint,istat,ioxyz,coo,ianz)

      return
      end

      subroutine gaupoi(istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call gaupod(istat,coo,ianz)

      return
      end

      subroutine getpoi(icomm,ifd,iff,idoscl,ioatms,ioadd)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /sclcom/ scal,fscal,scali,smag,iscupd
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call getpod(icomm,ifd,iff,idoscl,ioatms,ioadd,coo,ianz,iaton,
     &            iconn,scal,scali,smag,nat,ichx,icrtp)

      return
      end

      subroutine doconn
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call docond(coo,iconn,ianz)

      return
      end

      subroutine dohcon(ihflag)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /srfbck/ isurf(numat1)

      call dohcod(ihflag,coo,ianz,iconn,isurf)

      return
      end

      subroutine nohcon
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call nohcod(iconn)

      return
      end

      subroutine doscal
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /rotmat/rx(3),ry(3),rz(3),t(3),td(3)
      common /pers/  xv,yv,zv,c0,pincr
      common /sclcom/ scal,fscal,scali,smag,iscupd

      call doscad(t,coo,zv,pincr,scal,scali,smag,iscupd)

      return
      end

      subroutine redcon
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call redcod(iconn)

      return
      end

      subroutine domcon(ipt,idoall)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      monmod = 2

      call domcod(ipt,idoall,monmod,iconn)

      return
      end

      subroutine getxyz(igetxy,heat,iaddprv)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call getxyd(igetxy,heat,iaddprv,coo,ianz,iaton,iatclr,iconn,qat,
     &                  nat,norg,icent,nspg,ichx,nopr,ir,it,
     &                  xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      if (igetxy.eq.-1) then
          print*,'exceded maxnum of atoms!'
          igetxy = 0
      endif

      return
      end

      subroutine aln2ml(iopt,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call aln2md(iopt,istat,coo,ianz,iaton,iatclr,iresid)

      return
      end

      subroutine alnrot(vec,irot)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      dimension vec(*)

      call alnrod(vec,irot,coo)

      return
      end

      subroutine alnwrt
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call alnwrd(coo,ianz)

      return
      end

      subroutine alnsel(isel)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      dimension isel(4)

      call alnsed(isel,coo)

      return
      end

      subroutine getchr(igetca)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      common /chgbck/ qat(numat1)

      call getchd(igetca,qat)

      return
      end

      subroutine getmol(igetmo)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)

      call getmod(igetmo,coo,qat,ianz,iconn)

      return
      end

      subroutine outmol(iun)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call outmod(iun,coo,ianz,iconn)

      return
      end

      subroutine mopaco(istats,mopopt,irtype)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call mopacd(istats,mopopt,irtype,ianz,iatclr,iconn,coo,qat,
     &            nat,icent,nspg,ichx,nopr,ir,it,
     &            xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine parfc
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /sclcom/ scal,fscal,scali,smag,iscupd

      call parfd(coo,fscal)

      return
      end

      subroutine plmol
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      common /pltmp/  inat(numat1)
      common /pers/xv,yv,zv,c0,pincr
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      common /sclcom/ scal,fscal,scali,smag,iscupd

      call plmod(coo,qat,rzp,ixp,iyp,iconn,ianz,iaton,iatclr,iresid,
     &           inat,xv,yv,icalf,ncalf,icxp,icyp,
     &           scal,scali)

      return
      end

      subroutine plmolp
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      common /pltmp/  inat(numat1)
      common /pers/xv,yv,zv,c0,pincr
      common /sclcom/ scal,fscal,scali,smag,iscupd
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call plmold(coo,qat,rzp,ixp,iyp,iconn,ianz,iaton,iatclr,iresid,
     &           inat,xv,yv,zv,c0,
     &           icalf,ncalf,ianf,islu,nchain,iamino,reson,
     &           icxp,icyp,irsnr,scali)

      return
      end

      subroutine plfc(shade,ixx,k,ihigh,colsc,icltan,zvect)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /pers/xv,yv,zv,c0,pincr
      common /sclcom/ scal,fscal,scali,smag,iscupd
      integer shade
      dimension zvect(*)

      call plfd(shade,ixx,k,ihigh,colsc,icltan,zvect,
     &          coo,rzp,ixp,iyp,xv,yv,scal,scali)

      return
      end

      subroutine plfcp(shade,ixx,k,ihigh,colsc,icltan,zvect,scalp)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /pers/xv,yv,zv,c0,pincr
      common /sclcom/ scal,fscal,scali,smag,iscupd
      integer shade
      dimension zvect(*)

      call plfcd(shade,ixx,k,ihigh,colsc,icltan,zvect,scalp,
     &           coo,rzp,izp,iyp,xv,yv,zv,c0,scali)

      return
      end

      subroutine mstick(zvect,colsc,icltan,m,k,shade,idash)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /pers/xv,yv,zv,c0,pincr
      common /sclcom/ scal,fscal,scali,smag,iscupd
      integer shade
      dimension zvect(*)

      call msticd(zvect,colsc,icltan,m,k,shade,idash,
     &            coo,rzp,ixp,iyp,iatclr,zv,scali)

      return
      end

      subroutine mstck(zvect,colsc,icltan,m,k,shade,idash)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /sclcom/ scal,fscal,scali,smag,iscupd
      integer shade
      dimension zvect(*)

      call mstcd(zvect,colsc,icltan,m,k,shade,idash,
     &           coo,rzp,ixp,iyp,iatclr,scali)

      return
      end

      subroutine sstick
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /pers/xv,yv,zv,c0,pincr
      common /sclcom/ scal,fscal,scali,smag,iscupd

      call ssticd(coo,ixp,iyp,iaton,iconn,xv,yv,scal)

      return
      end

      subroutine astick(zvect,colsc,icltan,m,k,ihigh,shade,idash,imon)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /pers/xv,yv,zv,c0,pincr
      common /sclcom/ scal,fscal,scali,smag,iscupd
      integer shade
      dimension zvect(*)

      call asticd(zvect,colsc,icltan,m,k,ihigh,shade,idash,imon,
     &            coo,rzp,ixp,iyp,ianz,iatclr,xv,yv,scal,scali)

      return
      end

      subroutine astck(zvect,colsc,scalp,icltan,m,k,ihigh,shade,idash,
     &                 imon)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /pers/xv,yv,zv,c0,pincr
      common /sclcom/ scal,fscal,scali,smag,iscupd
      integer shade
      dimension zvect(*)

      call astcd(zvect,colsc,scalp,icltan,m,k,ihigh,shade,idash,imon,
     &           coo,rzp,ixp,iyp,ianz,iatclr,
     &           xv,yv,zv,c0,scali)

      return
      end

      subroutine snglat(zvect,scalp,colsc,icltan,k,ihigh,shade,atcol,
     &                  ipersp,ipost,icolps)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /pers/xv,yv,zv,c0,pincr
      common /sclcom/ scal,fscal,scali,smag,iscupd
      integer shade,atcol
      dimension zvect(*)

      call snglad(zvect,scalp,colsc,icltan,k,ihigh,shade,atcol,
     &            ipersp,ipost,icolps,coo,rzp,ixp,iyp,ianz,iatclr,
     &            xv,yv,zv,c0,scal,scali)

      return
      end

      subroutine fstick(roddef,rfac,zvect,colsc,icltan,m,k,ihigh,shade,
     &                  atcol,scnd)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /pers/xv,yv,zv,c0,pincr
      common /sclcom/ scal,fscal,scali,smag,iscupd
      integer shade,atcol
      logical scnd
      dimension zvect(*)

      call fsticd(roddef,rfac,zvect,colsc,icltan,m,k,ihigh,shade,
     &            atcol,scnd,coo,rzp,ianz,iatclr,xv,yv,scal,scali)

      return
      end

      subroutine fstck(roddef,rfac,scalp,zvect,colsc,icltan,m,k,ihigh,
     &                 shade,atcol,scnd)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /pers/xv,yv,zv,c0,pincr
      common /sclcom/ scal,fscal,scali,smag,iscupd
      integer shade,atcol
      logical scnd
      dimension zvect(*)

      call fstcd(roddef,rfac,scalp,zvect,colsc,icltan,m,k,ihigh,
     &           shade,atcol,scnd,coo,rzp,ianz,iatclr,
     &           xv,yv,zv,c0,scali)

      return
      end

      subroutine pldstg
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call pldstd(ixp,iyp)

      return
      end

      subroutine plalab(iqon)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      common /sclcom/ scal,fscal,scali,smag,iscupd

      call plalad(iqon,rzp,ixp,iyp,
     &            icalf,ianf,islu,nchain,iamino,reson,irsnr,achain,
     &            scali)

      return
      end

      subroutine plalb(iami)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call plald(iami,ixp,iyp,icalf,iamino,irsnr)

      return
      end

      subroutine pllab(ixp,iyp,ianz,k,qat,inr,iqon,ires,ipost)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call pllad(ixp,iyp,ianz,k,qat,inr,iqon,ires,ipost,ityp,ipbdt,
     &           iamino,reson)

      return
      end

      subroutine plpost(backb,dolabs,icolps,fancy,atcol,persp,shade,
     &                  idelx)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /pltmp/  inat(numat1)
      common /chgbck/ qat(numat1)
      common /pers/  xv,yv,zv,c0,pincr
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      common /sclcom/ scal,fscal,scali,smag,iscupd
      integer shade,atcol,persp,fancy,dolabs,backb

      call plposd(backb,dolabs,icolps,fancy,atcol,persp,shade,idelx,
     &            ixp,iyp,rzp,inat,qat,xv,yv,zv,
     &            icalf,ncalf,ianf,islu,nchain,iamino,reson,scal)

      return
      end

      subroutine plvrml(iun,fancy,atcol,dolabs,ihnd,backb,dohead)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      logical dohead
      integer atcol,fancy,dolabs,backb
      common /rotmat/rx(3),ry(3),rz(3),t(3),td(3)
      common /sclcom/ scal,fscal,scali,smag,iscupd
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call plvrmd(iun,fancy,atcol,dolabs,ihnd,backb,dohead,
     &            rz,t,coo,ianz,iaton,iatclr,iresid,iconn,
     &            icalf,ianf,islu,nchain,iamino,reson,scal,scali)

      return
      end

      subroutine plvst(iun,jcol,k,atcol,hndexl)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer atcol
      logical hndexl
      dimension jcol(3,16)

      call plvsd(iun,jcol,k,atcol,hndexl,coo,ianz,iaton,iatclr,iconn)

      return
      end

      subroutine hcoord(ioatms,nstrt)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      istat = 1

      call hcoodd(istat,ioatms,nstrt,ipdbt,coo,ianz,iaton,iresid,iconn,
     &            icalf,ncalf,ianf,islu,nchain,iamino)

      if (istat.eq.0) print*,'No room to add hydrogens'

      return
      end

      subroutine hang(idx1,idx2,hng)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call hand(idx1,idx2,hng,coo,ianz,iconn,icalf)

      return
      end

      subroutine hbond(ioatms,nstrt)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call hbodd(ioatms,nstrt,coo,ianz,iconn,
     &           icalf,ncalf,iamino)

      return
      end

      subroutine acthb(iopt,hbfilt)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call acthd(iopt,hbfilt,coo,iconn,ianz,icalf,ncalf)

      return
      end

      subroutine disabh(iopt)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call disabd(iopt,ianz,iaton,iconn)

      return
      end

      subroutine hbconn(iopt,iat1,iat2)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call hbcond(iopt,iat1,iat2,ianz,iatclr,iaton,iconn)

      return
      end

      subroutine getchn(ires,ichain)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call getcdh(ires,ichain,ianf,islu,nchain)

      return
      end

      subroutine ribbon(iscnd,idogl,nr,iungl,ipart,ist,ichain)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /rotmat/rx(3),ry(3),rz(3),t(3),td(3)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      istat = 1

      call ribbod(istat,iscnd,iidogl,nr,iungl,ipart,ist,ichain,t,
     &            coo,ianz,iaton,iatclr,iresid,iconn,
     &            icalf,ncalf,ianf,islu,nchain,iamino,isal,ihashb)

      if (istat.eq.0) then
         if (iscnd.eq.0) then
            print*,'No room to add helices'
         elseif (iscnd.eq.1) then
            print*,'No room to add beta sheets'
         elseif (iscnd.eq.2) then
            print*,'No room to add DNA backbone'
         endif
      endif

      return
      end

      subroutine proxim(itarg,thresh,ikleur,idosrf)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
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

      call proxid(itarg,thresh,ikleur,idosrf,coo,iresid,ishoh)

      return
      end

      subroutine proxic(itarg,backb,adds,idocom,thresh)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
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
      integer backb,adds

      call proxd(itarg,backb,adds,idocom,thresh,
     &           coo,iresid,iaton,iconn,iams,ishoh)

      return
      end

      subroutine rdchx(idebug,iop,istdbd,iuseab,moddma,istat,
     &                 icssr)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /chgbck/ qat(numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call rdchd(idebug,iop,istdbd,iuseab,moddma,istat,
     &           icrtp,icssr,coo,iconn,ianz,iatclr,ityp,qat,
     &           nat,icent,inorm,ncon,nspg,kz,ichx,icrtp,nopr,ir,it,
     &           xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)


      return
      end

      subroutine rfbio(idebug,ifrst,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call rfbid(idebug,ifrst,istat,coo,qat,ianz,iatclr,iconn,
     &                 nat,icent,inorm,ncon,nspg,kz,nopr,ir,it,
     &                 xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine getxdt(ipnt,nat,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call getxdd(ipnt,nat,istat,coo)

      return
      end

      subroutine rfdat(idebug,doest,istat,refcod)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma
      character*8 refcod
      logical doest

      call rfdad(idebug,doest,istat,refcod,coo,ianz,iatclr,iconn,
     &           nat,norg,icent,inorm,ncon,nspg,kz,nopr,ir,it,
     &           xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine fdat(iop,ifrst,istdbd,iuseab,moddma,idebug)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /zmsbck/ lwrit(numat1)
      common /hring/  lring(numat1)
      common /sclcom/ scal,fscal,scali,smag,iscupd
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call fdad(iop,ifrst,istdbd,iuseab,moddma,idebug,
     &     coo,qat,ianz,iaton,iatclr,iresid,iconn,lwrit,lring,ityp,
     &     scal,scali,smag,nat,norg,icent,inorm,ncon,nspg,nopr,ir,it,
     &     xa,ya,yb,za,zb,zc,a,b,c)


      return
      end

      subroutine delat(iat,iacn,iorg,mxdma)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)

      call delad(iat,iacn,iorg,mxdma,coo,ianz,iatclr,qat,ityp)

      return
      end

      subroutine polh(nat,ipolh,iatom)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call pold(nat,ipolh,iatom,ianz,iconn)

      return
      end

      subroutine ifnn(ifn,noff,ia1,ia2)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call ifnd(ifn,noff,ia1,ia2,iconn)

      return
      end

      subroutine rdgro(idebug,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /chgbck/ qat(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call rdgrod(idebug,istat,coo,qat,ianz,iatclr,iresid,iconn,
     &            ityp,ipdbt,
     &            icalf,ncalf,ianf,islu,nchain,iamino,reson,
     &            isal,irsnr,ishoh)

      if (istat.eq.-1) then
          print*,'exceeded maximum atoms !'
          istat = 0
      endif

      return
      end

      subroutine addbox(v1,v2,v3)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call addbod(v1,v2,v3,coo,ianz,iatclr,iconn)

      return
      end

      subroutine addtbx(v1,v2,v3)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call addtbd(v1,v2,v3,coo,ianz,iatclr,iconn,iaton,iresid)

      return
      end

      subroutine gropt(ipnt)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call gropd(ipnt,coo)

      return
      end

      subroutine rdmol(idebug,ipdbon,ioadd,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /chgbck/ qat(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      istat = 1

      call rdmod(idebug,ipdbon,ioadd,istat,coo,qat,
     &           ianz,iatclr,iresid,iconn,ityp,ipdbt,
     &           icalf,ncalf,ianf,islu,nchain,iamino,reson,
     &           isal,ishoh,ihashb,
     &           nat,norg,icent,ncon,nspg,kz,ichx,nopr,ir,it,
     &           xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)

      if (istat.eq.-1) then
          print*,'exceeded maximum atoms !'
          istat = 0
      endif

      return
      end

      subroutine wrmol(iun)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /chgbck/ qat(numat1)
      common /pltmp/  inat(numat1)
      common /hring/  lring(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call wrmod(iun,coo,qat,ianz,iaton,iatclr,iconn,
     &           iresid,lring,inat,ityp,ipdbt,
     &           icalf,ncalf,iamino,ishoh,
     &           nat,nspg,ichx,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine ispn(ispt,iat,irng,idochg,ifive)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      common /zmsbck/ lwrit(numat1)
      common /hring/  lring(numat1)

      call ispnd(ispt,iat,irng,idochg,ifive,
     &          qat,ianz,iaton,iconn,lwrit,lring)

      return
      end

      subroutine wrcrys(iun)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call wrcryd(iun,ianz,coo,
     &            nat,nspg,ichx,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine rdcif(idebug,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call rdcifd(idebug,istat,coo,ianz,iconn,iatclr,
     &            nat,icent,ncon,nspg,nopr,ir,it,
     &            xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine rdshlx(idebug,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call rdshld(idebug,istat,coo,ianz,iconn,iatclr,
     &            nat,icent,ncon,nspg,kz,nopr,ir,it,
     &            xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine wrshlx(iun,idospf)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call wrshld(iun,idospf,coo,ianz,
     &            nat,icent,nspg,kz,ichx,nopr,ir,it,
     &            a,b,c,alpha,beta,gamma)

      return
      end

      subroutine wrcif(iun)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call wrcifd(iun,coo,ianz,
     &            nat,nspg,ichx,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine cllrot(vec,irot,ifd)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma
      dimension vec(*)

      call cllrod(vec,irot,ifd,coo,ianz,
     &            nat,xa,ya,yb,za,zb,zc)

      return
      end

      subroutine zm2fr(cm,ctmp,imkeep)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma
      dimension cm(3,3),ctmp(3,*),imkeep(3)

      call zm2fd(cm,ctmp,imkeep,coo,ianz,iatclr,iconn,
     &           nat,ichx,xa,ya,yb,za,zb,zc)

      return
      end

      subroutine wrchx(iun)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /chgbck/ qat(numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call wrchd(iun,coo,ianz,iatclr,iconn,qat,ityp,
     &           nat,inorm,nspg,ichx,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine rdmsi(idebug,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call rdmsd(idebug,istat,coo,ianz,iatclr,iconn,qat,
     &           nat,icent,ncon,nspg,ichx,nopr,ir,it,
     &           xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine wrmsi(iun)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call wrmsd(iun,coo,ianz,iconn,qat,
     &           nat,nspg,ichx,xa,ya,yb,za,zb,zc)

      return
      end

      subroutine rdcpmolu(istats)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call rdcpmold(istats,coo,ianz)

      return
      end

      subroutine cpmdpt(ipnt,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call cpmdpd(ipnt,istat,coo,ianz)

      return
      end

      subroutine cpmdgetfr(istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call cpmdgetfd(istat,coo)

      return
      end

      subroutine cpmdcell
      implicit double precision (a-h,o-z), integer (i-n)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call cpmdceld(nspg,xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine cpmdptdyn(ipnt,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call cpmdptdyd(ipnt,istat,coo)

      return
      end

      subroutine wrcpmd(iun)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call wrcpdd(iun,ianz,coo,
     &            nat,ichx,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine rdmolu(sline,iemlin,idocoo,idobohr,istats)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call rdmodu(sline,iemlin,idocoo,idobohr,istats,coo,ianz)

      return
      end

      subroutine gamupt(ipnt,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call gamupd(ipnt,istat,coo,ianz)

      return
      end

      subroutine ugetfr(istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call ugetfd(istat,coo)

      return
      end

      subroutine wrvasp(iun)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call wrvasd(iun,coo,ianz,
     &            xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine wrconquest(iun)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call wrcqd(iun,coo,ianz,
     &            xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      
      return
      end

      subroutine wrmopa(iun)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call wrmopd(iun,coo,ianz,
     &            xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine wrfc(iun3)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call wrfd(iun3,coo)

      return
      end

      subroutine rdfc(ipnt,istats)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call rdfd(ipnt,istats,coo)

      return
      end

      subroutine rdmsf(idebug,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call rdmdf(idebug,istat,coo,ianz,iatclr,iconn,qat,ityp,
     &           nat,norg,icent,inorm,ncon,nspg,kz,ichx,nopr,ir,it,
     &           xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine rdpdb(istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
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

      call rdpdd(istat,coo,ianz,ihashb)

      return
      end

      subroutine pdbstd(istat,doscnd,ioadd)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      logical doscnd
   
      istat = 1

      call pdbsdd(istat,doscnd,ioadd,coo,ianz,iaton,iresid,iconn,ityp,
     &            ncalf,ianf,islu,nchain,iamino,ihet,
     &            isal,irsnr,achain,ihashb,ishoh)

      if (istat.eq.-1) then
          print*,'exceeded maximum atoms !'
          istat = 0
      endif

      return
      end

      subroutine conpdb
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call conpdd(ianz,iconn,iresid,ncalf,iamino)

      return
      end

      subroutine convpdb
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call convpdd(ianz,iconn,iresid,ncalf,iamino)

      return
      end

      subroutine conslv
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call consld(iconn,iresid)

      return
      end

      subroutine connij(idcon,i,j,idoconv)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call connid(idcon,i,j,idoconv,iconn,ianz,coo)

      return
      end

      subroutine setchg(iat,iopt)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /chgbck/ qat(numat1)

      call setchd(iat,iopt,qat,ityp)

      return
      end

      subroutine flagh(ihpdb,iat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      common /srfbck/ isurf(numat1)
      dimension ihpdb(*)

      call flagd(ihpdb,iat,isurf)

      return
      end

      subroutine conat(ipdb,iat1,iat2,iop)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      dimension ipdb(*)

      call conad(ipdb,iat1,iat2,iop,iconn)

      return
      end

      subroutine conath(ipdb,ihpdb,iat1,iat2)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      dimension ipdb(*),ihpdb(*)

      call conatd(ipdb,ihpdb,iat1,iat2,iconn)

      return
      end

      subroutine pdbtyp(ipdb,ihpdb,jres,ihashy)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      dimension ipdb(*),ihpdb(*)

      call pdbtyd(ipdb,ihpdb,jres,ihashy,ipdbt)

      return
      end

      subroutine typeit(ipdb,jres,ihpdb,ihashy)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      dimension ipdb(*),ihpdb(*)

      call typeid(ipdb,jres,ihpdb,ihashy,ianz,iconn,ityp)

      return
      end

      subroutine typamb(ipdb,jres,ihpdb,ihashy)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      dimension ipdb(*),ihpdb(*)

      call typamd(ipdb,jres,ihpdb,ihashy,ianz,iconn,ityp)

      return
      end

      subroutine typamo(ipdb,jres,ihpdb,ihashy)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      dimension ipdb(*),ihpdb(*)

      call typado(ipdb,jres,ihpdb,ihashy,ianz,iconn,ityp)

      return
      end

      subroutine mkback(ipdb,ihpdb,jres,icres,ihashy,idoconv)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      dimension ipdb(*),ihpdb(*)

      call mkbacd(ipdb,ihpdb,jres,icres,ihashy,idoconv,iconn,coo,
     &            icalf,ianf,islu,nchain,iamino)

      return
      end

      subroutine mknbck(ipdb,ihpdb,jres,icres,ihashy,idoconv)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      dimension ipdb(*),ihpdb(*)

      call mknbcd(ipdb,ihpdb,jres,icres,ihashy,idoconv,iconn,coo,
     &            icalf,ianf,islu,nchain,iamino)

      return
      end

      subroutine sftlab(i)
      implicit double precision (a-h,o-z), integer (i-n)

      idum = i

      return
      end

      subroutine addhs(ires,jres,ipdb,ihpdb,nterm)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /srfbck/ isurf(numat1)
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
      dimension ipdb(*),ihpdb(*)

      call addhd(ires,jres,ipdb,ihpdb,nterm,ianz,iaton,iatclr,iresid,
     &           iconn,isurf,ipdbt,ityp,ncalf,icalf,coo)

      return
      end

      subroutine numhet(num)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call numhed(num,iresid)

      return
      end

      subroutine qcxyz(idebug,nuclear,ipnt,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call qcxyd(idebug,nuclear,ipnt,istat,ianz,coo)

      return
      end

      subroutine orcxyz(idebug,ipnt,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call orcxyd(idebug,ipnt,istat,ianz,coo)

      return
      end

      subroutine nwxyz(idebug,ipnt,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call nwxyd(idebug,ipnt,istat,ianz,coo)

      return
      end

      subroutine fnwxyz(idebug,istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call fnwxyd(idebug,istat,ianz,coo)

      return
      end

      subroutine getqfr(istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call getqfd(istat,coo)

      return
      end

      subroutine getofr(istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call getofd(istat,coo)

      return
      end

      subroutine getnfr(istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call getnfd(istat,coo)

      return
      end

      subroutine rotcor(b)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      dimension b(3,*)

      call rotcod(b,coo)

      return
      end

      subroutine rotmom(ipoint,ifav)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call rotmod(ipoint,ifav,coo)

      return
      end

      subroutine rotfir(ioxyz)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call rotfid(ioxyz,ianz,coo)

      return
      end

      subroutine setorg(iatom)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /rotmat/rx(3),ry(3),rz(3),t(3), td(3)

      call setord(iatom,t,coo)

      return
      end

      subroutine samino(istat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call samind(istat,ianz,iresid,iconn,ipdbt,icalf,ncalf,iamino)

      return
      end

      subroutine getpdb(ires,ipdb,ihpdb)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      dimension ipdb(*),ihpdb(*)

      call getpdd(ires,ipdb,ihpdb,ianz,iresid,ipdbt)

      return
      end

      subroutine bckok(ibckok,ires,iop)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call bckod(ibckok,ires,iop,ianz,iresid,ipdbt,iamino)

      return
      end

      subroutine wrsrf(iun,nesp,iesp)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call wrsrd(iun,nesp,iesp,ianz,iatclr,iconn,coo)

      return
      end

      subroutine rdsrf(iun,istats,iesp,iaddprv,idebug)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call rdsrd(iun,istats,iesp,iaddprv,idebug,
     &           ianz,iaton,iatclr,iresid,iconn,coo)

      return
      end

      subroutine srfden(x,y,z,f)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /srfbck/ isurf(numat1)
      common /condens/ icont(numat1),ncont

      call srfded(x,y,z,f,ianz,isurf,coo,icont,ncont)

      return
      end

      subroutine defsrf
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /srfbck/ isurf(numat1)

      call defsrd(isurf,coo)

      return
      end

      subroutine clrmon
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call clrmod(iconn)

      return
      end

      subroutine intcor(intc,rout,isel,inum)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      real rout
      dimension isel(4)

      call intcod(intc,rout,isel,inum,coo)

      return
      end

      subroutine xyzcoo(idocopy,idoconv,ioadd)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call xyzcod(idocopy,idoconv,ioadd,ianz,coo)

      return
      end

      subroutine wrmsf(iun)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /chgbck/ qat(numat1)
      common /hring/  lring(numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call wrmsfd(iun,coo,qat,ianz,iaton,iatclr,iresid,
     &           iconn,lring,ityp,nat,nspg,ichx,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine wrtnk(iun)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /hring/  lring(numat1)

      call wrtnd(iun,ianz,iaton,iconn,lring,ityp,coo)

      return
      end

      subroutine wrgff(iun)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /chgbck/ qat(numat1)
      common /srfbck/ isurf(numat1)
      common /hring/  lring(numat1)
      common /zmsbck/ lwrit(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call wrgfd(iun,ianz,iaton,iconn,isurf,lring,lwrit,
     &           ncalf,ishoh,iresid,ityp,coo,qat,
     &           xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine wrogl(iun)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /hring/  lring(numat1)
      common /rotmat/rx(3),ry(3),rz(3),t(3), td(3)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call wrogd(iun,rx,ry,rz,lring,coo,rzp,reson,ianf,nchain,ncalf)

      return
      end

      subroutine rdbin(iun,heat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /hring/  lring(numat1)

      call rdbid(iun,heat,
     &            coo,ianz,iatclr,iconn,lring,
     &            ichx,nat,norg,
     &            xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine getnat(natoms)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /athlp/ iatoms, mxnat

      natoms = 0

      do i=1,iatoms
         if (ianz(i).lt.100.and.ianz(i).gt.0) natoms = natoms + 1
      end do

      return
      end

      subroutine allon(ioff,newat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /athlp/ iatoms, mxnat

      do i=ioff,newat
         iaton(i) = 1
         iatclr(i) = 1
      end do

      return
      end

      subroutine gettyp(ires,iat1,iat2)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call gettyd(ires,iat1,iat2,ityp,ipdbt,ianz,iresid,
     &            iamino,icalf,ncalf)

      return
      end

      subroutine chkbrk
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call chkbrd(iconn,icalf,ianf,islu,iamino,isal,reson,ncalf,nchain)

      return
      end

      subroutine gettnk(igttnk,idebug,ipdbon,iffset,iheat,heat)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /chgbck/ qat(numat1)
      common /srfbck/ isurf(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma
c      integer*2 ir,it
c      common /cell/   natc,ndum(6),ichx,icrtp,nopr,
c     &                ir(3,3,192),it(3,192),rdum(12)

      call gettnd(igttnk,idebug,ipdbon,iffset,iheat,heat,
     &            ianz,iconn,iatclr,ityp,coo,qat,isurf,
     &            issdon,iclon,ichx,ishoh,nspg,
     &            nat,norg,xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      if (igttnk.eq.-1) igttnk = 0

      return
      end

      subroutine tnkfst(igttnk,idebug)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /hring/  lring(numat1)
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma

      call tnkfsd(igttnk,idebug,coo,ianz,iatclr,iconn,lring,
     &            ichx,nat,norg,
     &            xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      if (igttnk.eq.-1) igttnk = 0

      return
      end

      subroutine dotyp(icel)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /zmsbck/ lwrit(numat1)
      common /hring/  lring(numat1)
      common /chgbck/ qat(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma
      common /condens/ icont(numat1),ncont

      call dotyd(icel,ianz,iaton,iatclr,iconn,iresid,
     &           lwrit,lring,ityp,coo,qat,icont,
     &           icalf,ncalf,ianf,islu,nchain,iamino,ishoh,
     &           nat,a,b,c,alpha,beta,gamma)

      return
      end

      subroutine espfit(idip,nesp,esp,connl,dx,dy,dz,iz,dmachg,ichadd)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      common /chgbck/ qat(numat1)
      logical dmachg
      dimension esp(*),connl(3,*)

      call espfid(idip,nesp,esp,connl,dx,dy,dz,iz,dmachg,ichadd,qat)

      return
      end

      subroutine allzmt(ipdbon)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call allzmd(ipdbon,ianz,iaton,coo)

      return
      end

      subroutine ligzmt
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call ligzmd(ianz,iaton)

      return
      end

      subroutine pdbzmt
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call pdbzmd(ianz,iaton,iresid,iconn)

      return
      end

      subroutine haswat(ino)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call haswad(ino,ianz,iaton,iresid,iconn)

      return
      end

      subroutine intzmt(ispdb)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /hring/  lring(numat1)
      common /srfbck/ isurf(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call intzmd(ispdb,ianz,iaton,iresid,iconn,lwrit,lring,
     &            icalf,ianf,islu,nchain,iamino)

      return
      end

      subroutine icrcon(icrcn,isel,idisc,ndisc,nanz,ispdb)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /zmsbck/ lwrit(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
      dimension isel(*),idisc(*)
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

      call icrcod(icrcn,isel,idisc,ndisc,nanz,ispdb,
     &            ianz,iaton,iconn,lwrit,coo,icalf,ncalf)

      return
      end

      subroutine calcx(ical,isel,nx)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      dimension isel(*)

      call calcd(ical,isel,nx,ianz,iaton,iconn,coo)

      return
      end

      subroutine prelea(iprel,ilead,isel,ispdb,ithree)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /zmsbck/ lwrit(numat1)
      dimension isel(*)

      call prelead(iprel,ilead,isel,ispdb,ithree,
     &             ianz,iaton,iresid,iconn,lwrit)

      return
      end

      subroutine preleh(iprel,ilead,isel,ispdb,ithree)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /zmsbck/ lwrit(numat1)
      dimension isel(*)

      call prelehd(iprel,ilead,isel,ispdb,ithree,
     &            ianz,iaton,iresid,iconn,lwrit)

      return
      end

      subroutine rdbas(idebug,idfree,istats)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numat1=20000)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)

      call rdbad(idebug,idfree,istats,ityp)

      return
      end

      subroutine setis(nres,istart)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /srfbck/ isurf(numat1)

      call setid(nres,istart,isurf,iresid,ipdbt)

      return
      end

      subroutine clkbck(istsurf,incp,ifogl)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call clkbcd(istsurf,incp,ifogl,iresid,coo,icalf,ncalf,reson)

      return
      end

      subroutine newfil(idebug,istat,inc,ioadd,ioatms,nstrt,namols,
     &                  nxtmf,ipdbon,namls,iof)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call newfid(idebug,istat,inc,ioadd,ioatms,nstrt,namols,
     &            nxtmf,ipdbon,namls,iof,iaton,iatclr,iresid,
     &            ncalf)

      return
      end

      subroutine setarr(iop,iopval,ioatms)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      common /chgbck/ qat(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      common /sclcom/ scal,fscal,scali,smag,iscupd
      common /pers/  xv,yv,zv,c0,pincr
      integer*2 ir,it
      common /cell/   nat,norg,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                nopr,ir(3,3,192),it(3,192),
     &                xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      common /pnthlp/ ipoints,ipnt

      call setard(iop,iopval,ioatms,
     &            ianz,iaton,iatclr,iresid,iconn,qat,
     &            ihet,iclhet,reson,iams,ihets,irsnr,
     &            ncalf,issdon,scal,scali,fscal,smag,
     &            xv,yv,xv,pincr,nat,ichx,icrtp,ipoints,ngeoms)

      return
      end

      subroutine epvrml(vdwr,moddma,natoms,norbs,idops)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call epvrmd(vdwr,moddma,natoms,norbs,idops,iaton)

      return
      end

      subroutine clrcod(natorg,natoms,idebug)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call clrcdd(natorg,natoms,idebug,iatclr)

      return
      end

      subroutine pckrot(istat,iatsel)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call pckrod(istat,iatsel,iresid,iconn)

      return
      end

      subroutine wrxyz(jmode)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      common /chgbck/ qat(numat1)

      call wrxyd(jmode,qat)

      return
      end

      subroutine caldip
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /chgbck/ qat(numat1)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call caldid(qat,coo)

      return
      end

      subroutine hetnum(idum)
      implicit double precision (a-h,o-z)

      idum = 0

      return
      end

      subroutine chkmol2(iok)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/ coo(3,numat1),rzp(numat1),ianz(numat1),
     &               iaton(numat1),iatclr(numat1),iresid(numat1),
     &               ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      
      call chkmld(iok,ianz,ityp)

      return
      end

      subroutine chktmp
      implicit double precision (a-h,o-z)

      idum = 0

      return
      end

      subroutine xyzrot(inct,theang)
      implicit double precision (a-h,o-z)
      common /rotmat/rx(3),ry(3),rz(3),t(3),td(3)

      call xyzrod(inct,theang,rx,ry,rz)

      return
      end

      subroutine rotts(x,y,z,xc,yc,zc,itran)
      implicit double precision (a-h,o-z)
      real x,y,z
      common /rotmat/rx(3),ry(3),rz(3),t(3),td(3)

      call rottd(x,y,z,xc,yc,zc,itran,
     &                 rx,ry,rz,t)

      return
      end

      subroutine rott(x,y,z,xc,yc,zc,itran)
      implicit double precision (a-h,o-z)
      common /rotmat/rx(3),ry(3),rz(3),t(3),td(3)

      call rotd(x,y,z,xc,yc,zc,itran,
     &                 rx,ry,rz,t)

      return
      end

      subroutine inirot
      implicit double precision (a-h,o-z)
      common /rotmat/rx(3),ry(3),rz(3),t(3),td(3)

      call inirod(rx,ry,rz,t)

      return
      end

      subroutine mtinv3
      implicit double precision (a-h,o-z)
      common /rotmat/rx(3),ry(3),rz(3),t(3),td(3)

      call mtind3(rx,ry,rz)

      return
      end

      subroutine mktrn(inct,incp)
      implicit double precision (a-h,o-z)
      common /pers/  xv,yv,zv,c0,pincr
      common /sclcom/ scal,fscal,scali,smag,iscupd

      call mktrd(inct,incp,xv,yv,zv,pincr,scal,scali,smag)

      return
      end

      subroutine setxyv
      implicit double precision (a-h,o-z)
      common /pers/  xv,yv,zv,c0,pincr

      call setxyd(xv,yv)

      return
      end

      subroutine chkbck(iupogl)
      implicit double precision (a-h,o-z)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call chkbcd(ncalf,ihet,reson)

      return
      end

      subroutine ribbs
      implicit double precision (a-h,o-z)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call ribbd(issdon)

      return
      end

      subroutine setchn(iopt)
      implicit double precision (a-h,o-z)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call setcdd(iopt,nchain,ianf,islu,irsnr,achain)

      return
      end

      subroutine acthlp(iop1,iop2,iop3)
      implicit double precision (a-h,o-z)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call acthld(iop1,iop2,iop3,ihet)

      return
      end

      subroutine acttog(iop1,iop2,iop3)
      implicit double precision (a-h,o-z)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call acttod(iop1,iop2,iop3,ihet)

      return
      end

      subroutine acttag(iop1,iop2,iop3)
      implicit double precision (a-h,o-z)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call acttad(iop1,iop2,iop3,ihet)

      return
      end

      subroutine calpsa(v,cval,pol,psa,tsa,exs,idepth)
      implicit double precision (a-h,o-z)

      idum = idepth

      return
      end

      subroutine tpsa
      implicit double precision (a-h,o-z)

      idum = i

      return
      end

      subroutine tstpsa(p1,p2,p3,i)
      implicit double precision (a-h,o-z)

      idum = i

      return
      end

      subroutine tomap(imap,iambfr)
      implicit double precision (a-h,o-z)

      imap = 0
      iambfr = 0

      return
      end

      subroutine upajob
      implicit double precision (a-h,o-z)

      idum = 0

      return
      end

      subroutine fndoh(itar,ang,ires,copt,istat)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/  coo(3,numat1),rzp(numat1),ianz(numat1),
     &                iaton(numat1),iatclr(numat1),iresid(numat1),
     &                ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /chgbck/ qat(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      dimension copt(3)

      call fndod(itar,ang,ires,copt,istat,
     &           coo,qat,iconn,ityp,iresid,irsnr,ncalf,icalf)


      return
      end

      subroutine evwat
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/  coo(3,numat1),rzp(numat1),ianz(numat1),
     &                iaton(numat1),iatclr(numat1),iresid(numat1),
     &                ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      common /chgbck/ qat(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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

      call evwad(coo,qat,iresid,irsnr,iatclr,iaton,iconn,ianz,
     &           ncalf,icalf,ityp,ipdbt)


      return
      end

      subroutine quwat(nwat)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/  coo(3,numat1),rzp(numat1),ianz(numat1),
     &                iaton(numat1),iatclr(numat1),iresid(numat1),
     &                ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call quwad(nwat,iresid,iconn,ianz)


      return
      end

      subroutine opthyd
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/  coo(3,numat1),rzp(numat1),ianz(numat1),
     &                iaton(numat1),iatclr(numat1),iresid(numat1),
     &                ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)
      integer*2 ityp,ipdbt
      common /typbck/ ityp(numat1),ipdbt(numat1)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
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
      common /chgbck/ qat(numat1)
      dimension iupres(500)

      call opthdd(iupres,nupres,irsnr,
     &            ncalf,icalf,coo,qat,iconn,ianz,iresid,iamino,
     &            ityp,ipdbt)

      return
      end

      subroutine gtfcor(i1,i2,i3,i4,i5)
      implicit double precision (a-h,o-z)

      idum = i1

      return
      end

      subroutine jdxwr
      implicit double precision (a-h,o-z)

      idum = i1

      return
      end

      subroutine chkcoo(kcoo,kcooh)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      common /atcom/  coo(3,numat1),rzp(numat1),ianz(numat1),
     &                iaton(numat1),iatclr(numat1),iresid(numat1),
     &                ixp(numat1),iyp(numat1),iconn(mxcon+1,numat1)

      call chkcod(kcoo,kcooh,ianz,iconn)

      return
      end

!!!!!! End of the xwindow dummys !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
