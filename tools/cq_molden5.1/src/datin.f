      subroutine datid(npts1,npts2,npts3,
     &                 ncols,ncolb,focc,focb,eiga,eigb)
c
c THIS IS REALLY datin
c
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      logical valenc,bonds,ovrlap,atomic,doori,ostep,
     &        fine,oscal,dolap
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /plspec/ step,stepf,scale,scalf,cut,pmin,pmax,one,
     &                oscal,ostep,fine
      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap
      common /keywrd/ keywrd,keyori
      common /orbhlp/ mxorb,iuhf,ispd
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      character*320 keywrd,keyori
      logical keyi,keyr
      real eiga,eigb
      dimension focc(*),focb(*),eiga(*),eigb(*)

      write(iun3,'(80(''*''))')
      write(iun3,'(''   DATA  PROVIDED BY INPUT FILE '')')

c----- Proces PlotPlane keywords -------------------------------

      call planky(npts1,npts2,npts3,keywrd,.true.)

      if (keyr(keywrd,'CUT',cut)) then
         if ((cut.gt.1.0d0).or.(cut.lt.0.0d0)) then
             cut = min(max(0.0d0,cut),1.0d0)
             call inferr('ERROR  0.0 <= CUT <= 1.0 ',1)
         endif
      endif

c---- use negative numbers to specify beta orbital with uhf

      imos = 0      
      ipsi = 0

      if (keyi(keywrd,'PSI',ipsi)) then
          if ((ipsi.gt.0.and.ipsi.gt.ncols).or.
     &        (ipsi.lt.0.and.iabs(ipsi).gt.ncolb)) then
             call inferr('Invalid Orbital Number!',1)
          else
             imos = imos + 1
          endif
      endif

      i=index(keywrd,'SPINDENS')
      if (i.ne.0) ispd = 1

      i=index(keywrd,'HOMO')
      if (i.ne.0) then
          imos = imos + 1
          call homo(ipsi)
      endif
      i=index(keywrd,'LUMO')
      if (i.ne.0) then
          imos = imos + 1
          call lumo(ipsi)
      endif
      if (imos.gt.1) then
*
*   MADE A MISTAKE
* 
          call inferr('INCOMPATIBLE KEY-WORDS !',1)
      endif
c
c----- the occupancy directive -----------------------

      call occup(istat)
      if (istat.eq.1) ipsi = 0

      write(iun3,'(80(''*''))')
      return
      end

      subroutine occupd(istat,focc,focb)
c
c----- the occupancy directive -----------------------
c----- example occu = (1-22/0,9/1.0,10/2 ) 
c----          sets occupancies of orbital 1 through 22 to zero
c----                                      9 to 1
c----                                     10 to 2
c----          in sequential order
c
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /keywrd/ keywrd,keyori
      character*320 keywrd,keyori
      dimension focc(*),focb(*)

      istat = 0

      i = index(keywrd,'OCCU')
      if (i.ne.0) then
          call occin(i,focc,norbs)
          istat = 1
      endif

      i = index(keywrd,'OCCA')
      if (i.ne.0) then
          call occin(i,focc,norbs)
          istat = 1
      endif

      i = index(keywrd,'OCCB')
      if (i.ne.0) then
          call occin(i,focb,norbs)
          istat = 1
      endif

      return
      end

      subroutine homod(ipsi,focc,focb,eiga,eigb)
      implicit double precision (a-h,o-z)
      common /orbhlp/ mxorb,iuhf,ispd
      real eiga,eigb
      dimension focc(*),focb(*),eiga(*),eigb(*)

      ipsi = 0

      if (iuhf.eq.1) then
         ipsa = 0
         do ii = 1,mxorb
            if (focc(ii).eq.0.0d0.and.ipsa.eq.0) ipsa = ii - 1
         end do
         ipsb = 0
         do ii = 1,mxorb
            if (focb(ii).eq.0.0d0.and.ipsb.eq.0) ipsb = ii - 1
         end do
         if (eiga(ipsa).lt.eigb(ipsb)) then
            ipsi = -ipsb
         else
            ipsi = ipsa
         endif
      else
         do ii = 1,mxorb
           if (focc(ii).eq.0.0d0.and.ipsi.eq.0) ipsi = ii - 1
         end do
      endif

      return
      end

      subroutine lumod(ipsi,focc,focb,eiga,eigb)
      implicit double precision (a-h,o-z)
      common /orbhlp/ mxorb,iuhf,ispd
      real eiga,eigb
      dimension focc(*),focb(*),eiga(*),eigb(*)

      ipsi = 0

      if (iuhf.eq.1) then
         ipsa = 0
         do ii = 1,mxorb
            if (focc(ii).eq.0.0d0.and.ipsa.eq.0) ipsa = ii
         end do
         ipsb = 0
         do ii = 1,mxorb
            if (focb(ii).eq.0.0d0.and.ipsb.eq.0) ipsb = ii
         end do
         if (eiga(ipsa).lt.eigb(ipsb)) then
            ipsi = ipsa
         else
            ipsi = -ipsb
         endif
      else
         do ii = 1,mxorb
           if (focc(ii).eq.0.0d0.and.ipsi.eq.0) ipsi = ii 
         end do
      endif

      return
      end
