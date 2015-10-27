*deck srtleg.f
c***begin prologue     srtleg
c***date written       930612   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           legend, link 6201, legendre functions
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            legendre functions
c***description        sort of p(l,m) functions
c***references         none
c
c***routines called
c***end prologue       legend
      subroutine srtleg(plmin,plmout,lval,npt,lmax,m,nl)
      implicit integer (a-z)
      real *8 plmin, plmout
      dimension plmin(npt,m:lmax), plmout(npt,nl), lval(nl)
c
      do 10 i=1,nl
         do 20 j=1,npt
            plmout(j,i)=plmin(j,lval(i))
   20    continue
   10 continue
      return
c
      end




