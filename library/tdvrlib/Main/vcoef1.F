*deck vcoef1.f
c***begin prologue     vcoef1
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***references         
c                        
c***routines called    
c***end prologue       vcoef1
      subroutine vcoef1(psi,u,vec,t,n,nt,nc)
      implicit integer (a-z)
      real*8 psi, u, vec, t
      character*80 title
      character*2 itoc
      dimension psi(n,nt,nc,2), vec(n), u(n,n), t(nt)
      common/io/inp, iout
      do 10 i=1,nt
         do 20 ic=1,nc
            write(iout,1) t(i), ic
            call ebtc(vec,u,psi(1,i,ic,1),n,n,1)
            write(iout,2) (vec(j),j=1,n)
            call ebtc(vec,u,psi(1,i,ic,2),n,n,1)
            write(iout,3) (vec(j),j=1,n)
 20      continue   
 10   continue   
      return
 1    format(/,'time = ',e15.8,1x,'channel = ',i3)
 2    format(/,1x,'projections for real part of psi',/,5(1x,e15.8))
 3    format(/,1x,'projections for imag part of psi',/,5(1x,e15.8))
      end       






