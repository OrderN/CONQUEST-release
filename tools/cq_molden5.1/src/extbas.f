      subroutine rdexbas(basfile)
      implicit double precision (a-h,o-z)
      parameter (maxbas=15)
      parameter (mxel=100)
      parameter (mxaorb=40)
      common /basext/ ibasin(mxel),eborbs(maxbas),
     &     pbas(mxaorb,mxaorb,maxbas)
      common /basopt/ extbas
      logical extbas
      character*75 basfile,line
      character*8 aname

      extbas = .false.

      open (unit=99,file=basfile,status='old',err=90)
      num = 1
      do while(.true.)
         
         read (99,'(A8)',end=10,err=100) aname
         read (99,*,err=100) line
         write(6,*) aname
         read (99,'(i4,i10)',err=100) icharge,norbs
         ibasin(icharge) = num
         eborbs(num) = norbs
         write(6,*) 'norbs=', eborbs(num),'icharge=', icharge
         do i=1,norbs
            read (99,'(6f15.10)',err=100) (pbas(i,j,num),j=1,norbs)
         end do
         num = num + 1
         read(99,*,err=100) line 
      end do
      close(99)

10    continue
      if (num.gt.1) extbas = .true.
      return

90    continue
      print*,'error opening file'
      return

100   continue
      return
      end

      subroutine newdenmad(
     &                     p,q)
c THIS IS REALLY newdenmat

c     returns density matrix in p common /densty/
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (maxbas=15)
      parameter (mxel=100)
      parameter (mxaorb=40)
      logical valenc,bonds,ovrlap,atomic,doori,orient,dolap
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap

      common /basext/ ibasin(mxel),eborbs(maxbas),
     &     pbas(mxaorb,mxaorb,maxbas)
      common /orbhlp/ mxorb,iuhf,ispd

      character*2 elemnt
      common /elem/elemnt(mxel)
      dimension p(*),q(*)

      do i=1,norbs
         do j=1,norbs
            q((j-1)*mxorb + i) = 0.0d0
         end do
      end do

      ix = 0
      do 700 k=1,natoms
         n = eborbs(ibasin(nat(k)))
         write(6,*) 'n=', n, 'K=',k
         if (ovrlap) then
c     
c============ do overlap densities only ========================
c
            do i=1,n
               do j=1,n
                     p((ix+j-1)*mxorb + (ix+i)) = 0.0d0
                  end do
               end do
            else
c
c============ molecule minus atoms =============================
c
               orient = .false.
C               do i=1,norien
C                  if (ori(i).eq.k) then
C                     orient = .true.
C                     ibal   = i
C                  endif
C               end do
               if (orient.and.doori) then
                  write (6,*) 'Orienting '
C                  
               else
c
c----- for all other atoms substract spherical density -----
c
                  do i=1,n
                     do j=1,n
c                        write (6,*) 'before:', ix,i,j,p((ix+j-1)*mxorb 
c     $                  + (ix+i)), pbas((ix+i),(ix+j),ibasin(nat(k))) 
                        p((ix+j-1)*mxorb + (ix+i)) = 
     $                   p((ix+j-1)*mxorb + (ix+i)) - pbas((ix+i),(ix+j)
     $                                 ,ibasin(nat(k)))
                        q((ix+j-1)*mxorb + (ix+i)) = 
     $                        pbas((ix+i),(ix+j),ibasin(nat(k)))
C                        write (6,*) 'after:' ,p(ix+i,ix+j)
                 end do
              end do
           endif
c
c============== clear out the overlap region of the density matrix
c               if ATOMIC is supplied
c
           if (atomic) then
              idavje = n 
              do i=ix+1,ix+idavje
                 if (k.ne.0) then
                    do j=1,ix
                       p((j-1)*mxorb + i) = 0.0d0
                    end do
                 endif
                 if (k.ne.natoms) then
                    do j=ix+idavje+1,norbs
                       p((j-1)*mxorb + i) = 0.0d0
                    end do
                 endif
              end do
           endif
        endif
        ix = ix + n 
 700  continue
      if (bonds.and.idebug.eq.1) then
         write(iun3,*)' '
         write(iun3,*)'***** Atomic Density Matrix *****'
         write(iun3,*)' '
         call prev(q,norbs,norbs,mxorb)
         write(iun3,*)' '
      endif
      
c############### end bonds , overlap , atomic ##################
      if (idebug.eq.1) then
          write(iun3,'(''   DENSITY MATRIX USED BY MAP'')')
          call prev(p,norbs,norbs,mxorb)
      endif

      return
      end
