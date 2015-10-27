*deck vcoef3.f
c***begin prologue     vcoef3
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***references         
c                        
c***routines called    
c***end prologue       vcoef3
      subroutine vcoef3(psi,u1,u2,u3,vec,t,n1,n2,n3,nt,nc)
      implicit integer (a-z)
      real*8 psi, u1, u2, u3, vec, t
      character*80 title
      dimension psi(n3,n2,n1,nt,nc,2), vec(n3,n2,n1,3)
      dimension u1(n1,n1), u2(n2,n2), u3(n3,n3), t(nt)
      common/io/inp, iout
      do 10 i=1,nt
         do 20 ic=1,nc
            write(iout,1) t(i), ic
             call ebtc(vec(1,1,1,3),u3,psi(1,1,1,i,ic,1),n3,n3,n2*n1)
            do 30 j=1,n1
               call ebc(vec(1,1,1,2),vec(1,1,j,3),u2,n3,n2,n2)
 30         continue
            call ebc(vec(1,1,1,1),vec(1,1,1,2),u1,n3*n2,n1,n1)
            title='projections for real part of psi'
            call prntrm(title,vec(1,1,1,1),n2,n1,n2,n1,iout) 
            call ebtc(vec(1,1,1,3),u3,psi(1,1,1,i,ic,2),n3,n3,n2*n1)
            do 40 j=1,n1
               call ebc(vec(1,1,1,2),vec(1,1,j,3),u2,n3,n2,n2)
 40         continue
            call ebc(vec(1,1,1,1),vec(1,1,1,2),u1,n3*n2,n1,n1)
            title='projections for imag part of psi'
            call prntrm(title,vec(1,1,1,2),n2,n1,n2,n1,iout) 
 20      continue   
 10   continue   
      return
 1    format(/,'time = ',e15.8,1x,'channel = ',i3)
 2    format(/,1x,'projections for real part of psi',/,5(1x,e15.8))
 3    format(/,1x,'projections for imag part of psi',/,5(1x,e15.8))
      end       

