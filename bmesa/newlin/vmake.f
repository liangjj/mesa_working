c $Header$
*deck vmake.f
c***begin prologue     vmake
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           vmake, link 6201, wavefunction
c***author             schneider, barry (lanl)
c***source             m6203
c***purpose            product of local potential and wavefunction
c***                   in co-ordinate space.       
c***references         
c***routines called
c***end prologue       vmake
      subroutine vmake (vloc,psi,rhs,npts,nchan,ntri)
      implicit integer(a-z)
      real*8 psi, vloc, rhs
      dimension psi(npts,nchan), vloc(npts,ntri), rhs(npts,nchan)
      ij=0
      do 20 ic=1,nchan
         do 30 jc=1,ic 
            ij=ij+1
            do 40 pts=1,npts
               rhs(pts,ic)=rhs(pts,ic)+vloc(pts,ij)*psi(pts,jc)
               rhs(pts,jc)=rhs(pts,jc)+vloc(pts,ij)*psi(pts,ic)
   40       continue
   30    continue
   20 continue
      return
      end
