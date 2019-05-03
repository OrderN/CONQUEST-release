      subroutine inferr(str,isev)
      implicit double precision (a-h,o-z)
c
c     print info/error
c
      character*(*) str
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /plot / iplot,iplwin,icolps
      logical zsev
      common /wrtcom/ zsev
cd     write(iun3,'(a)')'enter subroutine inferr'

      if (iplot.eq.6) then
          call molstr(str,len(str),iplwin)
          if (isev.eq.1) then
               write(iun3,*) ' '
               write(iun3,*) str
          endif
      else
          write(iun3,*) ' '
          write(iun3,*) str
          if (isev.eq.1.and.zsev) stop
      endif

cd     write(iun3,'(a)')'leave subroutine inferr'
      return
      end

      subroutine zmterr(str,icard,ivar,isev)
      implicit double precision (a-h,o-z)
c
c     print zmat parse error
c
      character*(*) str
      character*137 dumstr
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      common /plot / iplot,iplwin,icolps

      if (iplot.eq.6) then
          if (isev.eq.1) then
              dumstr = 'ZMAT NOT parsed: '//str
              call errzme(dumstr,len(str)+17,icard,ivar)
              write(iun3,*) str
              write(iun3,*) 'Center ',icard 
          else
              call errzme(str,len(str),icard,ivar)
          endif
      else
          write(iun3,*) str
c          if (isev.eq.1) stop
      endif

      return
      end
