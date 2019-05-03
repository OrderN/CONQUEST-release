      subroutine getarg(i,line)
      parameter (maxarg=80)
      integer i
      character*80 words
      character*(*) line
      common /gargs/ words(maxarg)

      if (i.le.maxarg) line = words(i)    

      return
      end

      integer function iargc()
      parameter (maxarg=80)
      integer j,k
      character*80 line,words
      common /gargs/ words(maxarg)

      iargc = 0
      call lib$get_foreign(line,,,)

      if (ifblen(line).eq.0) return

      if (index(line,' ').eq.0) then

         iargc = 1
         words(1) = line//' '
         return

      else

         do k=1,maxarg
            j = index(line,' ')
            if (j.eq.0) then
               words(k) = line(1:)//' '
               iargc = k
               return
            else
               words(k) = line(1:j-1)//' '
               iargc = k
               line = line(j:)
               if (ifblen(line).eq.0) return
               do i=1,80
                  if (line(1:1).eq.' ') then
                     line = line(2:)
                  else
                     goto 10
                  endif
               end do
10             continue
            endif
         end do

      endif
      
      iargc = maxarg

      return
      end
