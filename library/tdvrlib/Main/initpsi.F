*deck initpsi.f
c***begin prologue     initpsi
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            calculate time independent initial wavepacket.
c***references         
c
c***routines called    
c***end prologue       initpsi
      subroutine initpsi(psi,psi0,pt,n,nt,nc,prn)
      implicit integer (a-z)
      real*8 psi, psi0, pt
      character*80 title
      logical prn
      dimension psi(n,nt,nc,2), psi0(n,nc,2)
      dimension pt(nt,nt)
      common/io/inp, iout
      do 10 i=1,n
         do 20 j=1,nc
            psi0(i,j,1)=psi(i,nt,j,1)
            psi0(i,j,2)=psi(i,nt,j,2)
 20      continue	    
 10   continue
      prn=.true.
      if(prn) then   
         title='psi0 written to disk'
         call prntrm(title,psi0,n*nc,2,n*nc,2,iout)
      endif	 
      call iosys ('write real solution to bec',n*nc*2,psi0,0,' ')
      return
      end       
