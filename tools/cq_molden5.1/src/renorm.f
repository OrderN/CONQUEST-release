      subroutine renorm(exx,c,iang)
      implicit double precision (a-h,o-z)
      data pi32 /5.56832799683170d0/
      data pt187 /1.875d+00/
      data pt5,pt75,pt656 /0.5d0,0.75d0,6.5625d0/

      ee = exx + exx
      if (ee.ne.0.0d0) then
        facs = pi32 / (ee * dsqrt(ee))
        if (iang.eq.0) fac = facs
        if (iang.eq.1) fac = pt5   * facs / ee
        if (iang.eq.2) fac = pt75  * facs / (ee*ee)
        if (iang.eq.3) fac = pt187 * facs / (ee*ee*ee)
        if (iang.eq.4) fac = pt656 * facs / (ee*ee*ee*ee)
        c = c / dsqrt(fac)
      endif

      return
      end

      subroutine norml
      implicit double precision (a-h,o-z)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp


      do i = 1, nshell

c what type of shell: s, p, d, f, sp, spd

         iscod = 0
         if(shellt(i).eq.0) iscod = 0
         if(shellt(i).eq.1.and.shellc(i).eq.1) iscod = 2
         if(shellt(i).eq.1.and.shellc(i).ne.1) iscod = 1
         if(shellt(i).eq.2.and.shellc(i).eq.0) iscod = 4
         if(shellt(i).eq.2.and.shellc(i).ne.0) iscod = 3
         if(shellt(i).eq.3) iscod = 5
         if(shellt(i).eq.4) iscod = 6

         mm = shella(i)
         mmdf = shladf(i)
         n = shelln(i)

         if (iscod.eq.0) then
            call normp(exx(mm),c1(mm),n,0)
         elseif (iscod.eq.1) then
            call normp(exx(mm),c1(mm),n,0)
            call normp(exx(mm),c2(mm),n,1)
         elseif (iscod.eq.2) then
            call normp(exx(mm),c2(mm),n,1)
         elseif (iscod.eq.4) then
            call normp(exx(mm),c1(mm),n,0)
            call normp(exx(mm),c2(mm),n,1)
            call normp(exx(mm),c3(mmdf),n,2)
         elseif (iscod.eq.3) then
            call normp(exx(mm),c3(mmdf),n,2)
         elseif (iscod.eq.5) then
            call normp(exx(mm),c4(mmdf),n,3)
         elseif (iscod.eq.6) then
            call normp(exx(mm),c5(mmdf),n,4)
         endif

      end do

      return
      end

      subroutine normp(exx,c,n,iang)
      implicit double precision (a-h,o-z)
      dimension exx(n),c(n)

      small = 1.0d-15

      pow  = (3.0d0 + 2.0d0*dble(iang))/4.0d0
      sint = 0.0d0
      do i=1,n
         ei = exx(i)
         do j=1,n
            ej = exx(j)
            sint = sint + 
     &       c(i)*c(j)*(4.0d0*ei*ej/((ei+ej)**2.0d0))**pow
         end do
      end do

      if (sint.gt.small) then
         do i=1,n
            c(i) = c(i) / dsqrt(sint)
         end do
      endif

      do i=1,n
         call renorm(exx(i),c(i),iang)
      end do

      return
      end

