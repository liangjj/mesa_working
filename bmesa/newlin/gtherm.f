c $Header$
*deck gtherm.f
c***begin prologue     gtherm
c***date written       930119   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           gtherm, link 6201
c***author             schneider, barry (lanl)
c***source             m6203
c***purpose            get the phi(m) function corresponding to a given
c***                   m. the odd function always follows the even one
c***                   in storage.
c***references         
c***routines called
c***end prologue       gtherm
      subroutine gtherm (phifn,phich,wphi,nph,mval)
      implicit integer(a-z)
      real*8 phifn, phich, wphi
      dimension phifn(*), wphi(nph), phich(nph)
c        the first word address of phifn is passed as the location of
c        m = 0
      mabs=abs(mval)
c        locate the address of the first word of the m value needed
c        assuming that we want the even ( cosine ) function. if m is
c        zero there is nothing preceeding it in memory so location
c        cannot be smaller than zero.
      locte=max(0,(mabs+mabs-1)*nph)
c        if the odd ( sine ) function is wanted skip over nph words 
      if (mval(mu).lt.0) then
          locte=locte+nph
      endif
      call copy(phifn(locte+1),phich,nph)
c     scale the function by the weight for numerical integration
      call vmul(phich,phich,wphi,nph)
      return
      end
