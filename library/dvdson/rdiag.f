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
***end prologue       rdiag
      subroutine rdiag(v,cv,eig,ceig,vec,cvecl,cvecr,work,
     1                 lwork,rwork,n,m,mattyp)
      implicit integer (a-z)
      character*(*) mattyp
      character*80 title
      real*8 v, eig, vec, rwork
      complex*16 cv, ceig, cvecl, cvecr, work
      dimension v(n,*), eig(*), work(*), vec(n,*)
      dimension cv(n,*), ceig(*), cvecl(n,*), cvecr(n,*), rwork(*) 
      common/io/inp, iout
      if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
        call zgeev('v','v',m,cv,n,ceig,cvecl,n,cvecr,n,work,lwork,
     1              rwork,info)
        call eigord(eig,eig,eig,ceig,vec,cvecl,cvecr,cvecl,cvecr,
     1              m,n,'complex')
        if(info.ne.0) then
           call lnkerr('error from direct diagonalization routine')
        endif
c        title='left vectors'
c        call prntcm(title,cvecl,n,m,n,m,iout)
c        title='right vectors'
c        call prntcm(title,cvecr,n,m,n,m,iout)
        call renorm(cvecl,cvecr,n,m)
      else
         call tred2(n,m,v,eig,work,vec)
         call tql2(n,m,eig,work,vec,ierr)
         if(ierr.ne.0) then
            call lnkerr('error from direct diagonalization routine')
         endif   
      endif
      return
      end       

