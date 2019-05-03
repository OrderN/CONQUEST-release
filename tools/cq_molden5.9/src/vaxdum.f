      subroutine vaxconv
cvax      open(unit=45,form='formatted',file='scratch.com',
cvax     &status='unknown')
cvax      write(45,'(a)')'$ convert/fdl=sys$input plot.dat scratch.dat'
cvax      write(45,'(a)')'FILE'
cvax      write(45,'(a)')'ORGANIZATION SEQUENTIAL'
cvax      write(45,'(a)')'RECORD'
cvax      write(45,'(a)')'CARRIAGE_CONTROL CARRIAGE_RETURN'
cvax      write(45,'(a)')'FORMAT VARIABLE'
cvax      write(45,'(a)')'$ delete plot.dat;'
cvax      write(45,'(a)')'$ rename scratch.dat plot.dat'
cvax      write(45,'(a)')'$ delete scratch.com;'
cvax      close(45)
cvax      call Lib$Do_Command('@[]scratch.com')
      return
      end
