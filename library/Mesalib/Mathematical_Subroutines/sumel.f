*deck sumel
c***begin prologue     sumel
c***date written       940322   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           sum
c***author             schneider, barry (nsf)
c***source             matha
c***purpose            sum a vector.
c
c***references         none
c
c***routines called
c***end prologue       sumel
      function sumel (v,n)
      real*8 v, sumel
      dimension v(n)
      sumel=0.d+00
      do 10 i=1,n
         sumel=sumel+v(i)
   10    continue
      return
      end
