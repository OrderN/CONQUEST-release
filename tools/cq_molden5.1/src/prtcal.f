      subroutine prtcal(rphi,rpsi,
     &                  icalf,ncalf,ianf,islu,nchain,iamino,isal)
      implicit double precision (a-h,o-z)
      parameter (numcal=50000)
      parameter (mxres=42)
      character*3 aminos
      character*1 scnd(4)
      common /amino/aminos(mxres)
      common /hbonds/ hbd(2,numcal),ihb(2,numcal)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      dimension icalf(6,*),ianf(*),islu(*),iamino(*),isal(*)
      dimension rphi(*),rpsi(*)
      data scnd/'H','B','R','C'/

      do i=1,nchain
          nlen = islu(i)-ianf(i)+1
          write(iun3,*) 'Chain ',i,' Length = ',nlen
          write(iun3,*) '       '
          do j=ianf(i),islu(i),13
            kl = min(islu(i)-j+1,13)
            write(iun3,'(13(a3,1x))')(aminos(iamino(k+j-1)),k=1,kl)
          end do
      end do

      toang = 0.52917706d0
      write(iun3,*) ' nr res    hbond                scnd phi    psi'
      do i=1,ncalf
        if (iamino(i).le.23) then
         p1 = phi(i,icalf,ncalf,ianf,nchain)
         p2 = psi(i,icalf,ncalf,islu,nchain)
         rphi(i) = p1
         rpsi(i) = p2
         h1 = hbd(1,i)*toang
         h2 = hbd(2,i)*toang
         if (ihb(1,i).eq.0) h1 = 0.0d0
         if (ihb(2,i).eq.0) h2 = 0.0d0
         write(iun3,100) i,' ',aminos(iamino(i)),' ',
     &   ihb(1,i),',',ihb(2,i),' ',h1,',',h2,
     &   ' ',scnd(isal(i)+1),' ',p1,' ',p2
        endif
      end do

100   format(i4,a,a3,a,i4,a,i4,a,f6.3,a,f6.3,a,a,a,f6.1,a,f6.1)
      return
      end
