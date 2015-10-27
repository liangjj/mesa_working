*deck sqrtwt.f
c***begin prologue     sqrtwt
c***date written       930922   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           sqrtwt, link m6201, legendre integrals
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            calculate square root of integration weights.
c***references         none
c
c***routines called
c***end prologue       sqrtwt
      subroutine sqrtwt (wtthe,wtphi,nth,nph)
      implicit integer (a-z)
      real*8 wtthe, wtphi
      dimension wtthe(nth), wtphi(nph)
      do 10 n=1,nth
         wtthe(n)=sqrt(wtthe(n))          
   10 continue
      do 20 n=1,nph
         wtphi(n)=sqrt(wtphi(n))
   20 continue
      return
      end

