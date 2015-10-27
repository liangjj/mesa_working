*deck vcoef2.f
c***begin prologue     vcoef2
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***references         
c                        
c***routines called    
c***end prologue       vcoef2
      subroutine vcoef2(psi,u1,u2,vec,t,n1,n2,nt,nc)
      implicit integer (a-z)
      real*8 psi, u1, u2, vec, t
      character*80 title
      dimension psi(n2,n1,nt,nc,2), vec(n2,n1,2)
      dimension u1(n1,n1), u2(n2,n2), t(nt)
      common/io/inp, iout
      do 10 i=1,nt
         do 20 ic=1,nc
            write(iout,1) t(i), ic
            call ebtc(vec(1,1,2),u2,psi(1,1,i,ic,1),n2,n2,n1)
            call ebc(vec(1,1,1),vec(1,1,2),u1,n2,n1,n1)
            title='projections for real part of psi'
            call prntrm(title,vec(1,1,1),n2,n1,n2,n1,iout) 
            call ebtc(vec(1,1,2),u2,psi(1,1,i,ic,2),n2,n2,n1)
            call ebc(vec(1,1,1),vec(1,1,2),u1,n2,n1,n1)
            title='projections for imag part of psi'
            call prntrm(title,vec(1,1,1),n2,n1,n2,n1,iout) 
 20      continue   
 10   continue   
      return
 1    format(/,'time = ',e15.8,1x,'channel = ',i3)
      end       

