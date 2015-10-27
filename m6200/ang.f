*deck ang.f
c***begin prologue     ang
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           ang, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            fill angular quantities in lebedev quadrature
c***references         
c
c***routines called
c***end prologue       ang
      subroutine ang (pt,cthet,sthet,cphi,sphi,phpt,nleb)
      implicit integer (a-z)
      real*8 pt, cthet, sthet, cphi, sphi, phpt, twopi, rsq, r
      dimension pt(3,nleb), cthet(nleb), sthet(nleb), cphi(nleb)
      dimension sphi(nleb), phpt(nleb)
      common /io/ inp, iout
      data twopi /  6.283185307179586d0  /
      do 10 i=1,nleb
         rsq=pt(1,i)*pt(1,i)+pt(2,i)*pt(2,i)+pt(3,i)*pt(3,i)
         r=sqrt(rsq)
         cthet(i)=pt(3,i)/r
         sthet(i)=sqrt(1.d0-cthet(i)*cthet(i))
         if (sthet(i).ne.0.d0) then
             cphi(i)=pt(1,i)/(r*sthet(i))
             sphi(i)=pt(2,i)/(r*sthet(i))
         else
             cphi(i)=1.d0
             sphi(i)=0.d0
         endif         
         phpt(i)=atan2(sphi(i),cphi(i))
   10 continue 
      return
      end

