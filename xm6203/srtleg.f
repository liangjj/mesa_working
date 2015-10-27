*deck @(#)srtleg.f	1.2  10/27/94
c***begin prologue     srtleg
c***date written       930612   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           legend, legendre functions
c***author             schneider, barry (nsf)
c***source             
c***purpose            legendre functions
c***description        sort of p(l,m) functions
c***references         none
c
c***routines called
c***end prologue       srtleg
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




