*deck legend
c***begin prologue     legend
c***date written       880721   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           legend, link 2702, legendre functions
c***author             schneider, barry (lanl)
c***source             m2702
c***purpose            legendre functions
c***description        calculation of p(l,m) functions
c***references         none
c                      plm are the legendre functions l=m to l=lmax    
c                      x are the values of cos(theta)
c                      dfct are the factorials from 0 to maxfac
c                      ddfct are the double factorials from 0 to maxfac    
c***routines called
c***end prologue       legend
      subroutine legend (plm,x,dfct,ddfct,npt,lmax,m,maxfac)
      implicit integer (a-z)
      real *8 plm, x, dfct, ddfct, fm, facx, f1
      dimension plm(npt,m:lmax), x(npt)
      dimension  dfct(0:maxfac), ddfct(0:maxfac)
c----------------------------------------------------------------------c
c           start recursion with plm(m,m) and plm(m+1,m)               c
c                      and recur upward                                c
c----------------------------------------------------------------------c
      do 10 i=m,lmax
         do 10 j=1,npt
            plm(j,i)=0.d+00
   10 continue
      fm=.5d+00*m
      do 20 i=1,npt
         facx=1.d+00
         if (fm.ne.0.d+00) then
             facx=(1.d+00-x(i)*x(i))**fm
         endif
         plm(i,m)=ddfct(m)*facx
   20 continue
      if (lmax.ne.m) then
          mm=m+m+1
          mpls1=m+1
          do 30 i=1,npt
             plm(i,mpls1)=mm*x(i)*plm(i,m)
   30     continue
          if (lmax.ne.mpls1) then
              lind=m+2
              n1=2
              n2=m+m+3
              n3=n2-2
              do 50 i=lind,lmax
                 ii=i-1
                 jj=i-2
                 do 40 j=1,npt
                    f1=n2*x(j)*plm(j,ii)-n3*plm(j,jj)
                    f1=f1/n1
                    plm(j,i)=f1
   40            continue
                 n1=n1+1
                 n2=n2+2
                 n3=n3+1
   50       continue
          endif
      endif
      return
c
      end
