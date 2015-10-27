*deck @(#)resnum.f	1.1 9/8/91
c***begin prologue     resnum
c***date written       890605   (yymmdd)
c***revision date               (yymmdd)
c***keywords           resnum , m6008
c***author             schneider, barry (lanl)
c***source             m6008
c***purpose            restore diagonal bound-bound numerator matrix
c
c***references         none
c
c***routines called    none
c***end prologue       resnum
      subroutine resnum(hamnum,diag,matbb)
      implicit integer (a-z)
      real *8 hamnum, diag
      dimension hamnum(matbb,matbb), diag(matbb)
      do 10 i=1,matbb
         hamnum(i,i)=diag(i)
   10 continue
      return
      end
