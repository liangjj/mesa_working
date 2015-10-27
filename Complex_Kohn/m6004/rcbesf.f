*deck rcbesf
c***begin prologue     rcbesf
c***date written       xxxxxx   (yymmdd)
c***revision date      890422   (yymmdd)
c***keywords           m6004, link 6004, bessel, spline
c***author             schneider, barry (lanl)
c***source             m6004
c***purpose            generate bessel functions
c***description        bessel functions calculated by forward
c***                   recursion. use only when lmax lt arg
c***references       
c
c***routines called
c***end prologue       rcbesf
      subroutine rcbesf(x,xinv,j,jp,y,yp,nr,ptdim,lmax,ltop,prnt) 
      implicit integer (a-z)
      common /io/ inp, iout
      logical prnt
      real *8 x, j, jp, y, yp, xinv
      dimension j(ptdim,0:lmax), jp(ptdim,0:lmax), x(*), xinv(*)
      dimension y(ptdim,0:lmax), yp(ptdim,0:lmax)
      do 10 i=1,nr
         j(i,0)=sin(x(i))
         y(i,0)=-cos(x(i))
         jp(i,0)=-y(i,0)
         yp(i,0)=j(i,0)
   10 continue
      if (lmax.eq.0) then
          return
      else
          do 20 i=1,nr
             y(i,1)=y(i,0)*xinv(i)-j(i,0)
             j(i,1)=j(i,0)*xinv(i)+y(i,0)
   20     continue
          if (lmax.ne.1) then
              do 30 l=1,lmax-1
                 lp1=l+1
                 ll1=l+l+1
                 lm1=l-1
                 do 40 i=1,nr
                    j(i,lp1)=ll1*j(i,l)*xinv(i)-j(i,lm1)
                    y(i,lp1)=ll1*y(i,l)*xinv(i)-y(i,lm1)
   40            continue
   30         continue
          endif
          do  50 l=1,lmax
              lm1=l-1
              do 60 i=1,nr
                 jp(i,l)=j(i,lm1)-l*j(i,l)*xinv(i)
                 yp(i,l)=y(i,lm1)-l*y(i,l)*xinv(i)
   60         continue
   50     continue
      endif
      if (prnt) then
          do 800 l=0,lmax
             write (iout,900) l
             write (iout,1000) (j(i,l),i=1,nr)
             write (iout,1010) l
             write (iout,1000) (jp(i,l),i=1,nr)
             write (iout,1020) l
             write (iout,1000) (y(i,l),i=1,nr)
             write (iout,1030) l
             write (iout,1000) (yp(i,l),i=1,nr)
  800     continue
      endif
  900 format(/,5x,'regular function l = ',1x,i3)
 1000 format( (/,5x,5(d15.8,1x) ) )
 1010 format(/,5x,'derivative regular function l = ',1x,i3)
 1020 format(/,5x,'irregular function l = ',1x,i3)
 1030 format(/,5x,'derivative irregular function l = ',1x,i3)
      return
      end
