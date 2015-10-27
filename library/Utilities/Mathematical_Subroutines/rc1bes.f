*deck rc1bes
c***begin prologue     rc1bes
c***date written       xxxxxx   (yymmdd)
c***revision date      890422   (yymmdd)
c***keywords           m6004, link 6004, bessel, spline
c***author             schneider, barry (lanl)
c***source             m6004
c***purpose            generate bessel functions
c***description        bessel functions calculated using forward or
c***                   backward recursion. if greater accuracy needed
c***                   change parameter statement.
c***references       
c
c***routines called
c***end prologue       rc1bes
      subroutine rc1bes(x,j,jp,y,yp,lmax,ltop,dir,prnt)
      implicit integer (a-z)
      parameter (acc=30)
      common /io/ inp, iout
      real*8 one, two, x, scr, j, jp, y, yp, tmp
      logical prnt
      character*(*) dir
      dimension j(0:ltop), jp(0:ltop), y(0:ltop), yp(0:ltop)
      data one, two/ 1.d+00, 2.d+00 /
c----------------------------------------------------------------------c
c            estimate starting l                                       c
c----------------------------------------------------------------------c
      tmp=lmax*acc
      strtl=lmax+sqrt(tmp)
      strtl=max(strtl,lmax)
      if (strtl.gt.ltop) then
          call lnkerr('starting l bigger than ltop')
      endif
c----------------------------------------------------------------------c
c               make first two bessel functions                        c
c----------------------------------------------------------------------c
      j(0)=sin(x)
      y(0)=-cos(x)
      scr=j(0)
      if (lmax.ge.1) then
          j(1)=j(0)/x+y(0)
          y(1)=y(0)/x-j(0)
      endif
      if (lmax.ge.2) then
c----------------------------------------------------------------------c
c                   lmax is greater than one                           c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c              calculate y by upward recursion                         c
c                    its always stable                                 c
c----------------------------------------------------------------------c
          do 10 l=1,strtl-1
             ll1=l+l+1
             lm1=l-1
             lp1=l+1
             y(lp1)=ll1*y(l)/x-y(lm1)
   10     continue
c----------------------------------------------------------------------c
c             calculate j by upward or downward recursion              c
c             depending on the value of x                              c
c----------------------------------------------------------------------c
          j(strtl)=0.d0
          onelss=strtl-1
          j(onelss)=1.d-06             
          if (x.gt.dfloat(lmax)) then
c----------------------------------------------------------------------c
c               upward recursion                                       c
c----------------------------------------------------------------------c
              do 20 l=1,strtl-1
                 ll1=l+l+1
                 lm1=l-1
                 lp1=l+1
                 j(lp1)=ll1*j(l)/x-j(lm1)
   20         continue
          else
c----------------------------------------------------------------------c
c               downward recursion                                     c
c----------------------------------------------------------------------c
              do 30 l=onelss-1,0,-1
                 ll3=l+l+3
                 lp1=l+1
                 lp2=l+2
                 j(l)=ll3*j(lp1)/x-j(lp2)
   30         continue
c----------------------------------------------------------------------c
c                normalize the j                                       c
c----------------------------------------------------------------------c
                 scr=scr/j(0)
                 do 40 l=0,strtl
                    j(l)=j(l)*scr
   40            continue
          endif
      endif
c---------------------------------------------------------------------c
c             finish calculation by getting derivatives               c
c---------------------------------------------------------------------c
      if (dir.eq.'derivatives') then
          jp(0)=-y(0)
          yp(0)=j(0)
          do 50 l=1,lmax
             lm=l-1
             jp(l)=j(lm)-l*j(l)/x
             yp(l)=y(lm)-l*y(l)/x 
   50     continue
      endif
      if (prnt) then
          write (iout,60) lmax
          write (iout,100) (j(l),l=0,lmax)
          write (iout,70) lmax
          write (iout,100) (y(l),l=0,lmax)
          if (dir.eq.'derivatives') then
              write (iout,80) lmax
              write (iout,100) (jp(l),l=0,lmax)
              write (iout,90) lmax
              write (iout,100) (yp(l),l=0,lmax)
          endif 
      endif
   60 format(/,5x,'regular ricatti-bessel function l = 0 to ',i3)
   70 format(/,5x,'irregular ricatti-bessel function l = 0 to ',i3)
   80 format(/,5x,'derivative of regular ricatti-bessel function  l = 0
     1to ',i3)
   90 format(/,5x,'derivative of irregular ricatti-bessel function  l 
     1= 0 to ',i3)
  100 format( (/,5x,5(e15.8,1x) ) )
      return
      end
