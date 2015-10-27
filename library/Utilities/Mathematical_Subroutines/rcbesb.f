*deck rcbesb
c***begin prologue     rcbesb
c***date written       xxxxxx   (yymmdd)
c***revision date      890422   (yymmdd)
c***keywords           m6004, link 6004, bessel, spline
c***author             schneider, barry (lanl)
c***source             m6004
c***purpose            generate bessel functions
c***description        bessel functions calculated backward
c***                   recursion. if greater accuracy needed change
c***                   parameter statement.
c***references       
c
c***routines called
c***end prologue       rcbesb
      subroutine rcbesb(x,xinv,j,jp,y,yp,norm,np,ptdim,lmax,ltop,prnt)
      implicit integer (a-z)
      parameter (acc=30)
      common /io/ inp, iout
      real *8 one, two, x, j, jp, y, yp, xinv, norm, rl
      logical prnt
      dimension x(ptdim), xinv(ptdim), j(ptdim,0:ltop), jp(ptdim,0:ltop)
      dimension y(ptdim,0:ltop), yp(ptdim,0:ltop), norm(ptdim)
      data one, two/ 1.d+00, 2.d+00 /
c----------------------------------------------------------------------c
c            estimate starting l                                       c
c----------------------------------------------------------------------c
      rl=lmax*acc
      strtl=lmax+sqrt(rl)
      strtl=max(strtl,ltop)
      if (strtl.gt.ltop) then
          call lnkerr('starting l bigger than ltop')
      endif
c----------------------------------------------------------------------c
c               make first two bessel functions                        c
c----------------------------------------------------------------------c
      do 10 i=1,np
         j(i,0)=sin(x(i))
         norm(i)=j(i,0)
         y(i,0)=-cos(x(i))
         jp(i,0)=-y(i,0)
         yp(i,0)=j(i,0)
   10 continue
      if (lmax.eq.0) then
          return
      endif
      do 20 i=1,np
         j(i,1)=j(i,0)*xinv(i)+y(i,0)
         y(i,1)=y(i,0)*xinv(i)-j(i,0)
         jp(i,1)=j(i,0)-one*j(i,1)*xinv(i)
         yp(i,1)=y(i,0)-one*y(i,1)*xinv(i)    
   20 continue
      if (lmax.eq.1) then
          return
      endif       
c----------------------------------------------------------------------c
c                   lmax is greater than one                           c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c              calculate y by upward recursion                         c
c----------------------------------------------------------------------c
       lfinal=min(ltop,lmax+1)     
       do 30 l=1,lfinal-1
          ll1=l+l+1
          lm1=l-1
          lp1=l+1
          do 40 i=1,np
             y(i,lp1)=ll1*y(i,l)*xinv(i)-y(i,lm1)
   40     continue
   30  continue
c----------------------------------------------------------------------c
c             calculate j by downward recursion                        c
c----------------------------------------------------------------------c
       call rzero(j(1,strtl),np)
       onelss=strtl-1
       do 50 i=1,np
          j(i,onelss)=one
   50  continue
       do 60 l=onelss-1,0,-1
          ll3=l+l+3
          lp1=l+1
          lp2=l+2
          do 70 i=1,np
             j(i,l)=ll3*j(i,lp1)*xinv(i)-j(i,lp2)
   70     continue
   60  continue
c----------------------------------------------------------------------c
c                normalize the j                                       c
c----------------------------------------------------------------------c
       do 80 i=1,np
          norm(i)=norm(i)/j(i,0)
  80   continue
       do 90 l=0,lfinal
          do 100 i=1,np
             j(i,l)=j(i,l)*norm(i)
  100     continue
   90  continue
c---------------------------------------------------------------------c
c             finish calculation by getting derivatives               c
c---------------------------------------------------------------------c
       do 200 l=0,lmax
          lp1=l+1
          do 300 i=1,np
             jp(i,l)=lp1*j(i,l)*xinv(i)-j(i,lp1)            
             yp(i,l)=lp1*y(i,l)*xinv(i)-y(i,lp1) 
  300    continue
  200 continue
      if (prnt) then
          do 800 l=0,lmax
             write (iout,900) l
             write (iout,1000) (j(i,l),i=1,np)
             write (iout,1010) l
             write (iout,1000) (jp(i,l),i=1,np)
             write (iout,1020) l
             write (iout,1000) (y(i,l),i=1,np)
             write (iout,1030) l
             write (iout,1000) (yp(i,l),i=1,np)
  800     continue
      endif
  900 format(/,5x,'regular function l = ',1x,i3)
 1000 format( (/,5x,5(e15.8,1x) ) )
 1010 format(/,5x,'derivative regular function l = ',1x,i3)
 1020 format(/,5x,'irregular function l = ',1x,i3)
 1030 format(/,5x,'derivative irregular function l = ',1x,i3)
      return
      end
