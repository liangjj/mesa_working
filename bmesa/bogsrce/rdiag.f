*deck rdiag.f
c***begin prologue     rdiag
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           diagonalization
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for real diagonalization.
c***                   
c***references         
c
c***routines called    
c***end prologue       rdiag
      subroutine rdiag(v,cv,eig,ceig,vec,cvecl,cvecr,work,n,m,nonsym)
      implicit integer (a-z)
      logical nonsym
      real*8 v, eig, vec, work
      complex*16 cv, ceig, cvecl, cvecr
      character*80 title
      dimension v(n,*), eig(*), work(*), vec(n,*)
      dimension cv(n,*), ceig(*), cvecl(n,*), cvecr(n,*) 
      common/io/inp, iout
      if(nonsym) then 
        call cgeev(cv,n,m,ceig,cvecr,n,work,1,ierr)
        call cc2opy(cvecr,cvecl,n*m)
c         call cgevlrv(cv,n,m,ceig,cvecl,cvecr,n,work,
c     1                'left eigenvectors',ierr)
         call eigord(eig,ceig,vec,cvecl,cvecr,m,n,'complex')
c         title='eigenvalues'
c         call prntcm(title,ceig,m,1,m,1,iout)
c         title='right eigenvectors'
c         call prntcm(title,cvecr,n,m,n,m,iout)
c         title='left eigenvectors'
c         call prntcm(title,cvecl,n,m,n,m,iout)
      else
         call tred2(n,m,v,eig,work,vec)
         call tql2(n,m,eig,work,vec,ierr)
      endif
      if(ierr.ne.0) then
         call lnkerr('error from direct diagonalization routine')
      endif         
      return
      end       
