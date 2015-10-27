*deck hdiag.f
c***begin prologue     hdiag
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           finite difference, band, eigenvalues
c***author             schneider, barry(nsf)
c***source             
c***purpose            calculate eigenvalues/eigenvectors of tridiagonal
c***                   matrix.
c***
c***references
c
c***routines called
c***end prologue      hdiag
      subroutine hdiag(band,eig,vec,work,dum,s,ind,n,order,bw,evecs,
     1                 nroots,prnt)
      implicit integer (a-z)
      dimension band(n,-bw:bw), eig(n), vec(n,*) 
      dimension work(n,-bw:bw), dum(n,3), ind(n), s(nroots,nroots)
      real*8 band, vec, work, eps1, lb, ub, dum, eig, s
      logical evecs, prnt
      character*80 title
      common /io/ inp, iout
      if(evecs) then
c
c                 get eigenvalues as well as vectors
c
c        get all eigenvalues and eigenvectors-this is expensive for
c        big matrices
c       
         if(nroots.eq.n-2) then
            call figi2(n,n-2,band(2,-1),work(2,0),work(2,-1),vec(2,1),
     1                                                       ierr)
            call imtql2(n,n-2,work(2,0),work(2,-1),vec(2,1),ierr)
            nel=max(1,n-2)
            write(iout,1) (work(i,0),i=2,nel+1)
            call copy(work(2,0),eig(2),nel)
         else
c         
c        get subset of eigenvalues and eigenvectors
c
            call figi(n,n-2,band(2,-1),work(2,0),work(2,-1),
     1                                        work(2,1),ierr)
            eps1=0.d0
            call tridib(n-2,eps1,work(2,0),work(2,-1),work(2,1),lb,ub,
     1                  1,nroots,dum(2,1),ind(2),ierr,dum(1,2),
     2                  dum(1,3))
            nel=max(1,nroots)
            write(iout,1) (dum(i,1),i=2,nel+1)
            call copy(dum(2,1),eig(2),nel)
            call tinvit(n,n-2,work(2,0),work(2,-1),work(2,1),nroots,
     1                  dum(2,1),ind(2),vec(2,1),ierr,dum(1,2),dum(1,3),
     2                  dum(1,4),dum(1,5),dum(1,6))
            call bakvec(n,n-2,band(2,-1),work(2,-1),nroots,
     1                  vec(2,1),ierr)
            if (prnt) then
                call ebtcxx(s,vec(2,1),vec(2,1),nroots,n-2,nroots,
     1                      nroots,n,n)
                title='overlap in hdiag'
                call prntrm(title,s,nroots,nroots,nroots,nroots,iout)
                title='vectors in hdiag'
                call prntrm(title,vec,n,nroots,n,nroots,iout)
            endif
         endif                
      else 
         call figi(n,n-2,band(2,-1),work(2,0),work(2,-1),
     1                                        work(2,-1),ierr)
         call imtql1(n-2,work(2,0),work(2,-1),ierr)
         nel=max(1,n-2)
         write(iout,1) (work(i,0),i=2,nel+1)         
         call copy(work(2,0),eig(2),nel)
      endif
      return
 1    format(/,'the eigenvalues',/,(5e15.8))
      end

