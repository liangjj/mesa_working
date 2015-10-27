*deck ovusym.f
c***begin prologue     ovusym
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            projection of perturbed symmetrized states
c***                   onto unperturbed states.
c***                   
c***references         
c
c***routines called    
c***end prologue       ovusym
      subroutine ovusym(ov,scr,vec,vecx0,vecy0,n1,n2,n,prn)
      implicit integer (a-z)
      real*8 ov, scr, vec, vecx0, vecy0
      character*80 title
      logical prn
      dimension ov(n,n), scr(n2,n1,n), vec(n2,n1,n)
      dimension vecx0(n1,n1), vecy0(n2,n2)
      common/io/inp, iout
c
c     project onto coordinate 1. 
c     scr(j,t,q) = Sum vec(j,i,q) * vecx0(i,t)
c
      do 10 i=1,n
         call ebc(scr(1,1,i),vec(1,1,i),vecx0,n2,n1,n1)
 10   continue   
      nprd=n1*n
c
c     finish projecting onto coordinate 2
c     ov(I,J) = Sum vecy0(j,t') * scr(j,t,q)
c
      call ebtc(ov,vecy0,scr,n2,n2,nprd)
      if(prn) then
         title='< psi0(t) | psi(q) >'
         call prntfm(title,ov,n,n,n,n,iout)
      endif          
      return      
      end       






