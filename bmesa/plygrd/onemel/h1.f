*deck h1.f
c***begin prologue     h1
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            find time-dependent wavefunction by expansion.
c***                   
c***description        simple time-dependent hamiltonian using polynomial
c***                   basis.  
c         
c                           d
c                      hbar _    
c                           dt
c***references         
c
c***routines called    
c***end prologue       h1
      subroutine h1(p,dp,q,wt,hamt,hbar,phrse,reg,n,prn)
      implicit integer (a-z)
      real*8 p, dp, q, wt, hamt, hbar
      real*8 y
      complex*16 eye
      character*(*) reg, phrse
      logical prn
      character*80 title
      dimension p(n,n), dp(n,n), q(n)
      dimension hamt(n,n), wt(n)
      common/io/inp, iout
      pointer(py,y(1))
      mat=1
      eigr=mat+n*n
      eigi=eigr+n
      vl=eigi+n
      ovlp=vl
      vr=vl+n*n
      work=vr+n*n
      lwork=12*n
      eigc=work
      vlc=eigc+lwork
      vrc=vlc+2*n*n
      need=wpadti(vrc+2*n*n)
      call memory(need,py,ngot,'y',0)
      call rzero(hamt,n*n)
      do 10 i=1,n
         do 20 j=1,n
            hamt(i,j) = wt(i)*p(i,i)*dp(i,j)
 20      continue
 10   continue
      call sscal(n*n,hbar,hamt,1)
      if (prn) then
          title='time-dependent hamiltonian'
          call prntrm(title,hamt,n,n,n,n,iout)
      endif     
      call copy(hamt,y(mat),n*n)
      call dgeev('v','v',n,y(mat),n,y(eigr),y(eigi),y(vl),n,
     1            y(vr),n,y(work),lwork,info)
      if(info.ne.0) then
         write(iout,*) info
         call lnkerr('error from direct diagonalization routine')
      endif
c
      call tocmplx(y(eigc),y(vlc),y(vrc),y(eigr),y(eigi),
     1             y(vl),y(vr),y(ovlp),n,prn)
c
c          test orthonormality
c
      call iosys('write real "dt mtrx-'//reg//' for '//
     1            phrse//'" to bec',n*n,hamt,0,' ')
      call iosys('write real "t-eigvals-'//reg//
     1           ' for '//phrse//'" to bec',2*n,y(eigc),0,' ')
      call iosys('write real "left t-eigvecs-'//reg//' for '//
     1            phrse//'" to bec',2*n*n,y(vlc),0,' ')
      call iosys('write real "right t-eigvecs-'//reg//' for '//
     1            phrse//'" to bec',2*n*n,y(vrc),0,' ')
      call memory(-ngot,py,idum,'y',idum)
      return
      end       









