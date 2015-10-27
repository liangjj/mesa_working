*deck @(#)modnum.f	1.1 9/8/91
c***begin prologue     modnum
c***date written       890605   (yymmdd)
c***revision date               (yymmdd)
c***keywords           modnum , m6008
c***author             schneider, barry (lanl)
c***source             m6008
c***purpose            modify bound-bound numerator matrix
c
c***references         none
c
c***routines called    modnum
c***end prologue      
      subroutine modnum(hamnum,diag,eng,matbb)
      implicit integer (a-z)
      real *8 hamnum, diag, eng
      dimension hamnum(matbb,matbb), diag(matbb)
c----------------------------------------------------------------------c
c              subtract energy and store diagonal                      c
c----------------------------------------------------------------------c
      do 10 i=1,matbb
         diag(i)=hamnum(i,i)
         hamnum(i,i)=hamnum(i,i)-eng
   10 continue
      return
      end
