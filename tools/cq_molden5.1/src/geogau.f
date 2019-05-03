      subroutine geogau(formax,forrms,dismax,disrms,epoints,isav)
      implicit double precision (a-h,o-z)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      common /gauori/ nzm,nso,nio,nzo,ioropt,ifor,
     &                ixyz98,iopr,isymm,irc,imp2,icntp,itd
      common /gauver/ ivers
      character*137 line,str
      character*25 sirc
      common /curlin/ line
      integer getlin
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*)

      sirc = 'Item               Value'
      nis = 24
      if (irc.eq.2) then
         sirc = '  Delta-x Convergence'
         nis = 22
      endif

      rewind iun2
      isemi = 0
      ioni = 0
      call searchd(line,'ONIOM: generating','calculation of energy',
     &             istat)
      if (istat.ne.0) then
         if (icdex(line,'ONIOM:').ne.0) then
            ioni = 1
         elseif (icdex(line,'External').ne.0) then
            isemi = 1
         elseif (icdex(line,'second derivatives').eq.0) then
            if (icdex(line,'AM1').ne.0) isemi = 1
            if (icdex(line,'MNDO').ne.0) isemi = 1
            if (icdex(line,'MINDO3').ne.0) isemi = 1
            if (icdex(line,'PM3').ne.0) isemi = 1
            if (icdex(line,'CNDO').ne.0) isemi = 1
            if (icdex(line,'INDO').ne.0) isemi = 1
            if (icdex(line,'AMBER').ne.0) isemi = 1
            if (icdex(line,'UFF').ne.0) isemi = 1
            if (icdex(line,'Dreiding').ne.0) isemi = 1
            if (icdex(line,'MM2').ne.0) isemi = 1
            if (icdex(line,'MM3').ne.0) isemi = 1

         endif
      endif
      rewind iun2
      
      igcvav = 1
      ifmxav = 1
      ifrmav = 1
      idmxav = 1
      idrmav = 1
      ieav = 1


      nepnts = 0
      ngeoms = 0
      igeoms = 0
      do i=1,mxpnt
         if (isemi.eq.1) then
            call srclit(line,' Energy=',istat)
         elseif (ioni.eq.1) then
            call search(line,'ONIOM: extrapolated energy =',istat)
         elseif (icntp.eq.1) then
            call search(line,'Counterpoise: corrected energy =',istat)
         elseif (imp2.eq.1) then
            call searchd(line,'EUMP2 =','E(RMP2)=',istat)
         elseif (itd.eq.1) then
            call search(line,'E(TD-HF/TD-KS)',istat)
         else
            call searchd(line,'SCF Done:','MCSCF converged',istat)
         endif
         if (istat.eq.0) goto 100

10       if (nepnts.eq.mxpnt) then
            call inferr('exceeded MAXPNT !',1)
            return
         endif

         etmp = 0.0d0
         if (irc.eq.0) nepnts = nepnts + 1

         imcscf = 0
         if (icdex(line,'MCSCF').ne.0) then
             imcscf = 1
c             call scback(line,'... DO AN EXTRA-',istat)
             call scback(line,'... Do an extra-',istat)
             backspace iun2
             if (ivers.lt.2009) then
                backspace iun2
             endif
             read(iun2,'(a)') line
             i1 = index(line,'E=')
             if (i1.gt.0) i1 = i1 + 1
         elseif (icdex(line,'MP2').ne.0) then
             it1 = icdex(line,'EUMP2 =')
             if (it1.ne.0) i1 = it1 + 6
             it1 = icdex(line,'E(RMP2)=')
             if (it1.ne.0) i1 = it1 + 9
         else
             i1 = index(line,'=')
         endif
         if (i1.gt.0) then
             if (irc.gt.0) then
                etmp = reada(line,i1+1,len(line))
             else
                epoints(nepnts) = reada(line,i1+1,len(line))
             endif
         endif

         if (imcscf.eq.1) then
             call search(line,'MCSCF converged',istat)
             call searchd(line,
     &           'MCSCF converged',sirc(1:nis),istat)
         else
             if (isemi.eq.1) then
                call searchd(line,
     &           'Energy=',sirc(1:nis),istat)
             elseif (ioni.eq.1) then
                call searchd(line,
     &           'ONIOM: extrapolated energy =',sirc(1:nis),istat)
             elseif (icntp.eq.1) then
                call searchd(line,'Counterpoise: corrected energy',
     &           sirc(1:nis),istat)
             elseif (imp2.eq.1) then
                call searcht(line,'EUMP2 =','E(RMP2)=',
     &           sirc(1:nis),istat)
             elseif (itd.eq.1) then
                call searchd(line,'E(TD-HF/TD-KS)',
     &           sirc(1:nis),istat)
             else
                call searcht(line,
     &           'SCF Done:',sirc(1:nis),
     &           'Corrected End Point Energy ',istat)
             endif
         endif
         if (istat.eq.0) goto 100

         if (icdex(line,'Corrected').ne.0) then

             i1 = index(line,'=')
             etmp = reada(line,i1+1,len(line))
             nepnts = nepnts + 1
             igeoms = igeoms + 1
             epoints(nepnts) = etmp
             isav(nepnts) = 0 
    
         else if (icdex(line,'Item').ne.0) then


c            read(iun2,'(26x,f8.6,5x,f8.6)',end=100,err=100) tmp1,fmaxt
            icv1 = 0
            icv2 = 0
            icv3 = 0
            icv4 = 0
            if (getlin(0).eq.1) then
                if (icdex(line,'YES').ne.0) icv1 = 1
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.1) goto 100
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.1) goto 100
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.1) tmp1 = 0.0d0
                if (ktype.eq.3) tmp1 = rtype
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.1) fmaxt = 0.0d0
                if (ktype.eq.3) fmaxt = rtype
            else
                goto 100
            endif

c            read(iun2,'(26x,f8.6,5x,f8.6)',end=100,err=100) tmp2,frmst
            if (getlin(0).eq.1) then
                if (icdex(line,'YES').ne.0) icv2 = 1
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.1) goto 100
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.1) goto 100
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.1) tmp2 = 0.0d0
                if (ktype.eq.3) tmp2 = rtype
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.1) frmst = 0.0d0
                if (ktype.eq.3) frmst = rtype
            else
                goto 100
            endif

c            read(iun2,'(26x,f8.6,5x,f8.6)',end=100,err=100) tmp3,dmaxt
            if (getlin(0).eq.1) then
                if (icdex(line,'YES').ne.0) icv3 = 1
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.1) goto 100
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.1) goto 100
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.1) tmp3 = 0.0d0
                if (ktype.eq.3) tmp3 = rtype
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.1) dmaxt = 0.0d0
                if (ktype.eq.3) dmaxt = rtype
            else
                goto 100
            endif

c            read(iun2,'(26x,f8.6,5x,f8.6)',end=100,err=100) tmp4,drmst
            if (getlin(0).eq.1) then
                if (icdex(line,'YES').ne.0) icv4 = 1
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.1) goto 100
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.1) goto 100
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.1) tmp4 = 0.0d0
                if (ktype.eq.3) tmp4 = rtype
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.1) drmst = 0.0d0
                if (ktype.eq.3) drmst = rtype
            else
                goto 100
            endif

            icv = icv1 + icv2 + icv3 + icv4
            if (irc.eq.1) then
                if (icv.ne.4) then
                   goto 10
                else
                   nepnts = nepnts + 1
                   epoints(nepnts) = etmp
                endif
            endif

            igeoms = igeoms + 1
            isav(nepnts) = 1
            formax(nepnts) = tmp1
            forrms(nepnts) = tmp2
            dismax(nepnts) = tmp3
            disrms(nepnts) = tmp4

         else if (icdex(line,'Delta-x').ne.0) then

            if (icdex(line,' NOT met').ne.0) goto 10
            nepnts = nepnts + 1
            igeoms = igeoms + 1
            epoints(nepnts) = etmp
            isav(nepnts) = 0 

         else

            isav(nepnts) = 0 
            goto 10

         endif

      end do

100   continue

      if (igeoms.eq.0.or.igeoms.eq.1) then
         ngeoms = igeoms
         igcvav = 0
         ifmxav = 0
         ifrmav = 0
         idmxav = 0
         idrmav = 0
         ieav = 0
      else
         ngeoms = nepnts
      endif

      return
      end
