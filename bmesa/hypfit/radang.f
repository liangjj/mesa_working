*deck radang.f
c***begin prologue     radang
c***date written       001230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            radial and angular hyperspherical variables
c***                   in cartesian coordinates.
c***                   
c***references         
c
c***routines called    
c***end prologue       radang
      subroutine radang(rho,ang,r1,r2,n1,n2)
      implicit integer (a-z)
      real*8 rho, ang, r1, r2
      real*8 pi2, tanalf, zero, atan, nrzero
      character*80 title
      dimension rho(n2,n1), ang(n2,n1), r1(n1), r2(n2)
      common/io/inp, iout
      data pi2 / 1.5707963267948966192313215d+00 /
      data zero / 0.d0 /
      data nrzero / 1.d-14 /
      do 10 i=1,n1
         do 20 j=1,n2
            rho(j,i) = sqrt( r1(i)*r1(i)+r2(j)*r2(j) )
            if(rho(j,i).gt.nrzero) then
               if(r2(j).ne.zero) then
                  tanalf=r1(i)/r2(j)
                  ang(j,i)=atan(tanalf)
               else
                  ang(j,i)=pi2
               endif
            else
               ang(j,i)=zero
            endif
 20      continue
 10   continue   
      return
      end       















