*deck rpadiag.f
c***begin prologue     rpadiag
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           diagonalization
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for real non-symmetric diagonalization.
c***                   
c***references         
c
c***routines called    
c***end prologue       rpadiag
      subroutine rpadiag(a,er,ei,vl,vr,work,n,m,lwork)
      implicit integer (a-z)
      character*80 title
      real*8 a, er, ei, vl, vr, work, temp, sdot
      dimension a(n,*), er(*), ei(*), vl(n,*), vr(n,*), work(*)
      common/io/inp, iout
      call dgeev('v','v',m,a,n,er,ei,vl,n,vr,n,work,lwork,info)
c
c     check that eigenvalues are all real
c
      do 10 i=1,m
         if(abs(ei(i)).gt.1.d-10) then
            write(iout,1) i, ei(i)
            call lnkerr('rpa eigenvalues are not positive. quit')
         endif
 10   continue
c
c     order the real parts
c
      do 20 ii=2,m
         i=ii-1
         k=i
         temp=er(i)
         do 30 j=ii,m
            if(er(j).lt.temp) then
               k=j
               temp=er(j)
            endif   
 30      continue
         if(k.ne.i) then
            er(k) = er(i)
            er(i) = temp
            call sswap (m,vl(1,i),1,vl(1,k),1)
            call sswap (m,vr(1,i),1,vr(1,k),1)
         endif
 20   continue   
      do 40 i=1,m
         temp=1.d0/sqrt(sdot(m,vl(1,i),1,vr(1,i),1))
         call vscale(vl(1,i),vl(1,i),temp,m)
         call vscale(vr(1,i),vr(1,i),temp,m)
 40   continue   
c      title='left vectors'
c      call prntrm(title,vl,m,m,n,n,iout)
c      title='right vectors'
c      call prntrm(title,vr,m,m,n,n,iout)
      return
 1    format(/,5x,'the imaginary part of the ',i4,'d rpa eigenvalue = ',
     1            d15.8,/,5x,' it is above 1.d-10. quit')
      end       






