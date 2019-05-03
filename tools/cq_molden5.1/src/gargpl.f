       logical function gargpl(fl,n,strng1,strng2)
       integer n,fln
       character*(*) fl
       character*(*) strng1,strng2

       gargpl = .false.
       
       fln = len(fl)
       if (strng1(1:fln).ne.fl) return 
       if (strng1(fln+1:fln+1).eq.' ') then
           n = n + 1
           call getarg(n,strng2)
           if (strng2(1:1).eq.'-'.or.strng2(1:1).eq.' ') return
           strng2 = strng2(1:index(strng2,' ')-1)
       else
           strng2 = strng1(fln+1:index(strng1,' ')-1)
       endif
       gargpl = .true.

       return
       end
