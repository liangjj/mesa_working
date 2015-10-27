*deck @(#)rcbes.f	1.2  10/27/94
c***begin prologue     rcbes
c***date written       xxxxxx   (yymmdd)
c***revision date      890422   (yymmdd)
c***keywords           m6004, link 6004, bessel, spline
c***author             schneider, barry (lanl)
c***source             m6004
c***purpose            generate ricatti-bessel functions
c***description        ricatti-bessel functions calculated using forward or
c***                   backward recursion. if greater accuracy needed
c***                   change parameter statement.
c***references       
c
c***routines called
c***end prologue       rcbes
      subroutine rcbes(x,j,jp,y,yp,wron,scr,np,lmax,ltop,dir,prnt)
      implicit integer (a-z)
      parameter (acc=30)
      common /io/ inp, iout
      real*8 one, two, x, scr, j, jp, y, yp, tmp, wron
      logical prnt
      character*(*) dir
      dimension x(np), j(np,0:ltop), jp(np,0:ltop)
      dimension y(np,0:ltop), yp(np,0:ltop), scr(np,2), wron(0:ltop)
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
      do 10 i=1,np
         j(i,0)=sin(x(i))
         y(i,0)=-cos(x(i))
   10 continue
      call copy(j(1,0),scr(1,1),np)
      if (lmax.ge.1) then
          call vinv(scr(1,2),x,np)
          do 20 i=1,np
             j(i,1)=j(i,0)*scr(i,2)+y(i,0)
             y(i,1)=y(i,0)*scr(i,2)-j(i,0)
   20     continue
      endif
      if (lmax.ge.2) then
c----------------------------------------------------------------------c
c                   lmax is greater than one                           c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c              calculate y by upward recursion                         c
c                    its always stable                                 c
c----------------------------------------------------------------------c
          do 30 l=1,strtl-1
             ll1=l+l+1
             lm1=l-1
             lp1=l+1
             do 40 i=1,np
                y(i,lp1)=ll1*y(i,l)*scr(i,2)-y(i,lm1)
   40        continue
   30     continue
c----------------------------------------------------------------------c
c             calculate j by upward or downward recursion              c
c             depending on the value of x                              c
c----------------------------------------------------------------------c
          call rzero(j(1,strtl),np)
          onelss=strtl-1
          call vfill(j(1,onelss),1.d-06,np)             
          do 60 i=1,np
             if (x(i).gt.dfloat(lmax)) then
c----------------------------------------------------------------------c
c               upward recursion                                       c
c----------------------------------------------------------------------c
                 do 70 l=1,strtl-1
                    ll1=l+l+1
                    lm1=l-1
                    lp1=l+1
                    j(i,lp1)=ll1*j(i,l)*scr(i,2)-j(i,lm1)
   70            continue
             else
c----------------------------------------------------------------------c
c               downward recursion                                     c
c----------------------------------------------------------------------c
                 do 80 l=onelss-1,0,-1
                    ll3=l+l+3
                    lp1=l+1
                    lp2=l+2
                    j(i,l)=ll3*j(i,lp1)*scr(i,2)-j(i,lp2)
   80            continue
c----------------------------------------------------------------------c
c                normalize the j                                       c
c----------------------------------------------------------------------c
                 scr(i,1)=scr(i,1)/j(i,0)
                 do 90 l=0,strtl
                    j(i,l)=j(i,l)*scr(i,1)
   90            continue
             endif
   60     continue
      endif
c---------------------------------------------------------------------c
c             finish calculation by getting derivatives               c
c---------------------------------------------------------------------c
      if (dir.eq.'derivatives') then
          call vneg(jp(1,0),y(1,0),np)
          call copy(j(1,0),yp(1,0),np)
          wron(0)=j(1,0)*yp(1,0)-jp(1,0)*y(1,0)
          do 200 l=1,lmax
             lm=l-1
             do 300 i=1,np
                jp(i,l)=j(i,lm)-scr(i,2)*j(i,l)
                yp(i,l)=y(i,lm)-scr(i,2)*y(i,l) 
  300        continue
             wron(l)=j(1,l)*yp(1,l)-jp(1,l)*y(1,l)
  200     continue
      endif
      if (prnt) then
          do 800 l=0,lmax
             write (iout,900) l
             write (iout,1000) (j(i,l),i=1,np)
             write (iout,1020) l
             write (iout,1000) (y(i,l),i=1,np)
             if (dir.eq.'derivatives') then
                 write (iout,1010) l
                 write (iout,1000) (jp(i,l),i=1,np)
                 write (iout,1030) l
                 write (iout,1000) (yp(i,l),i=1,np)
                 write (iout,1040) l, wron(l)
             endif 
  800     continue
      endif
  900 format(/,5x,'regular function l = ',1x,i3)
 1000 format( (/,5x,5(e15.8,1x) ) )
 1010 format(/,5x,'derivative regular function l = ',1x,i3)
 1020 format(/,5x,'irregular function l = ',1x,i3)
 1030 format(/,5x,'derivative irregular function l = ',1x,i3)
 1040 format(/,5x,'wronskian for l = ',1x,i3,1x,e15.8)
      return
      end
