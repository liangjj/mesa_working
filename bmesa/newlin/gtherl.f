c $Header$
*deck gtherl.f
c***begin prologue     gtherl
c***date written       930119   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           gtherl, link 6201
c***author             schneider, barry (lanl)
c***source             m6203
c***purpose            pack functions needed for angular momentum
c***                   projections on a given channel in order with
c***                   no zeros.
c***references         
c***routines called
c***end prologue       gtherl
      subroutine gtherl (plm,plmc,wthet,nth,lval,nl,m)
      implicit integer(a-z)
      real*8 plm, plmc, wthet
      dimension plm(*), plmc(*), wthet(nth), lval(nl)
c     the plm are passed into this routine starting at l = abs ( m )
c     so to avoid a data error the minimum value of lval is examined
c     to make sure it is not smaller than m.
      minl=0
      do 10 l=1,nl
         minl=min(minl,lval(l))
   10 continue
      if (minl.lt.m) then
          call lnkerr('error in l values passed to gtherl')
      endif
      locpl=0
      do 20 l=1,nl
         locte=(lval(l)-m)*nth
         call copy(plm(locte+1),plmc(locpl+1),nth)
c        scale the function by the weights for numerical integration
         call vmul(plmc(locpl+1),plmc(locpl+1),wthet,nth)
         locpl=locpl+nth
   20 continue
      return
      end
