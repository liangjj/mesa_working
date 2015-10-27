*deck @(#)lowdin.f	1.1  11/30/90
      subroutine lowdin(s,uin,uout,eigval,t1,t2,tol,moin,moout,n,prn)
c***begin prologue     lowdin
c***date written       850601  (yymmdd)
c***revision date      860813  (yymmdd)
c
c***keywords           matrix, orthonormalize
c***author             schneider, barry (nsf)
c***source
c***purpose                                             
c***                   form an orthonormal set of vectors by diagonalization
c***                   of the overlap
c***description
c
c***references
c***routines called    
c***end prologue       lowdin
      implicit integer (a-z)
c
      real*8 s(n,n), uin(n,n), uout(n,n), eigval(*), t1(n,n), t2(n,n)
      real*8 fac, tol 
      logical prn(2)
      character*80 title
c
      common /io/ inp,iout
c
c     ----- print the overlap matrix if requested -----
c
      if (prn(2)) then
         title=' the a.o. overlap matrix '
         call prntrm (title,s,n,n,n,n,iout)
      end if
c
c     form overlap matrix of input orbitals
c
      call ebcxx(t1,s,uin,n,n,moin,n,n,n)
      call ebtcxx(t2,uin,t1,moin,n,moin,n,n,n)
c
c     ----- diagonalize s and get eigenvalues and vectors -----
c
      call dsyev('v','l',moin,t2,n,eigval,t1,5*moin,info)
      if(prn(1)) then
         title='eigenvalues'
         call prntrm(title,eigval,moin,1,moin,1,iout)
      endif
      moout=0
      do 10 i=1,moin
         if(abs(eigval(i)).gt.tol) then
             moout=moout+1
             eigval(moout)=eigval(i)          
             fac=1.d0/sqrt(eigval(moout))
             do 20 j=1,moin
                t2(j,moout)=t2(j,i)*fac
 20          continue   
         else
             write(iout,1) eigval(i), tol
         endif
 10   continue   
      write(iout,2) moin, moout
c
      call ebcxx(uout,uin,t2,n,moin,moout,n,n,n)
      if(prn(2)) then      
         title='eigenvectors'
         call prntrm(title,uout,n,moout,n,n,iout)
      endif
      return
 1    format(/,1x,'bad eigenvalue = ',e15.8,
     1       /,1x,'tol            = ',e15.8)
 2    format(/,1x,'number of vectors in  = ',i4,/,1x,
     1            'number of vectors out = ',i4)
      end
