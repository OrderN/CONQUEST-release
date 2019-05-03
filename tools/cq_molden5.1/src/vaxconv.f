      subroutine vaxconv
      open(unit=45,form='formatted',file='scratch.com',
     &status='unknown')
      write(45,'(a)')'$ convert/fdl=sys$input plot.dat scratch.dat'
      write(45,'(a)')'FILE'
      write(45,'(a)')'ORGANIZATION SEQUENTIAL'
      write(45,'(a)')'RECORD'
      write(45,'(a)')'CARRIAGE_CONTROL CARRIAGE_RETURN'
      write(45,'(a)')'FORMAT VARIABLE'
      write(45,'(a)')'$ delete plot.dat;'
      write(45,'(a)')'$ rename scratch.dat plot.dat'
      write(45,'(a)')'$ delete scratch.com;'
      close(45)
      call Lib$Do_Command('@[]scratch.com')
      return
      end

      subroutine runjob(idum,iqopt,ihaszm)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      idum = 1
      print*,'RUNJOB not implemented on VMS'
      return
      end

      subroutine tnkpnt(ipnt,iret)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)

      print*,'RUNJOB not implemented on VMS'

      return
      end

