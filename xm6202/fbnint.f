*deck fbnint
c***begin prologue     fbnint
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            greens function quadrature
c***
c***description        performs the integration of the radial
c***                   greens function and an inhomogeneity using
c***                   newton-cotes quadrature               
c***references
c
c***routines called
c
c***end prologue       fbnint
c
      subroutine fbnint(flm,psilm,r,j,y,wt,scr,lmax,m,n,nshell,nr,
     1                  nints,prnt)
      implicit integer (a-z)
      real*8 flm, psilm, wt, j, y, scr, int0f, int0b, r
      logical prnt
      dimension flm(n,m:lmax), wt(*), psilm(nints,m:lmax), j(n,0:lmax)
      dimension y(n,0:lmax), nr(nshell), scr(nints,2), r(n)
      common /io/ inp, iout
c      do 75 i=1,n
c      j(i,0)=r(i)
c      j(i,1)=r(i)*j(i,0)
c      j(i,2)=j(i,1)*r(i)
c      y(i,0)=1.d0
c      y(i,1)=y(i,0)/r(i)
c      y(i,2)=y(i,1)/r(i)
c 75   continue
c----------------------------------------------------------------------c     
c                when we are done we have in psilm                     c
c                                                                      c
c            psilm(i) = y(i) * sum ( k=1..i ) j(k)*flm(k)*wt(k) +      c
c                       j(i) * sum ( k=i+1...n ) y(k)*flm(k)*wt(k)     c
c            which is the radial wavefunction at all of the quadrature c
c            points not including the origin.                          c
c            thus we have constructed a quadrature approximation to;   c
c                                                                      c
c            psilm(x) = y(x) * integral(0,x) [ j(x) * flm(x) ] +       c
c                       j(x) * integral(x,infinity) [ y(x)*flm(x0 ]    c
c                                                                      c
c----------------------------------------------------------------------c
      do 10 l=m,lmax
c----------------------------------------------------------------------c
c             initialize the integral at its first point               c
c----------------------------------------------------------------------c      
         int0f=0.d0
         locf=1
         locwt=1
         locint=1
         do 20 ns=1,nshell
            call forwrd(flm(locf,l),j(locf,l),y(locf,l),wt(locwt),
     1                  int0f,scr(locint,1),scr(locint,2),nr(ns))
            locf=locf+nr(ns)
            locwt=locwt+nr(ns)*(nr(ns)-1)
            locint=locint+nr(ns)-1
  20     continue
         if (prnt) then
             write(iout,1) m, l, (scr(i,2),i=1,nints)
         endif    
         call copy(scr(1,1),psilm(1,l),nints) 
c----------------------------------------------------------------------c
c             initialize the integral at the last point                c
c----------------------------------------------------------------------c
         int0b=0.d0
         do 30 ns=nshell,1,-1
            locf=locf-nr(ns)
            locwt=locwt-nr(ns)*(nr(ns)-1)
            locint=locint-nr(ns)+1
            call bkwrd(flm(locf,l),j(locf,l),y(locf,l),wt(locwt),
     1                 int0b,scr(locint,1),scr(locint,2),nr(ns))
   30    continue
         if (prnt) then
             write(iout,2) m, l, (scr(i,2),i=nints,1,-1)
         endif   
         do 40 int=1,nints-1
            psilm(int,l)=psilm(int,l)+scr(int+1,1)
   40    continue
   10 continue
      return
    1 format(/,'forward partial integral for m = ',i3,1x,'l = ',i3,
     1          (/,5(e15.8,1x)))
    2 format(/,'backward partial integral for m = ',i3,1x,'l = ',i3,
     1          (/,5(e15.8,1x)))
      end















