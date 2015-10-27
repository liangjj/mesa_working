c $Header$
*deck gtherr.f
c***begin prologue     gtherr
c***date written       930119   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           gtherr, link 6201
c***author             schneider, barry (nsf)
c***source             m6203
c***purpose            pack radial functions for each sector, channel
c***                   and l with no zeros,
c***references         
c***routines called
c***end prologue       gtherr
      subroutine gtherr (glft,grght,gr1,gr2,lval,nr,nl,ldim)
      implicit integer(a-z)
      real*8 gr1, gr2, glft, grght
      dimension glft(nr,nl), grght(nr,nl)
      dimension gr1(nr,0:ldim), gr2(nr,0:ldim), lval(nl)
      do 10 l=1,nl
         call copy(gr1(1,lval(l)),glft(1,l),nr)
         call copy(gr2(1,lval(l)),grght(1,l),nr)
   10 continue
      return
      end
