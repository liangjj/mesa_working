*deck @(#)tomobb.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001, spline
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            transform bound functions to orthonormal basis 
c***
c***description        an orthonormal basis is found by diagonalization
c***                   of the overlap matrix. near linearly dependent
c***                   vectors are discarded and the functions transformed.
c***                 
c***                 
c***references       
c
c***routines called    tred2(clams) tql2(clams)
c***end prologue       m7001
      subroutine tomobb(fns,ddfns,mask,lb,mb,fmo,s,eig,work,tol,nsym,
     1                  n,nbf,nout,dimsym,prnt,orth)
      implicit integer (a-z)
      real *8 fns, ddfns, fmo, s, eig, work, tol
      logical prnt, orth
      dimension fns(n,nbf), ddfns(n,nbf), s(nbf,nbf), work(nbf)
      dimension eig(nbf), nsym(dimsym), fmo(n,nbf), lb(nbf), mb(nbf)
      dimension mask(nbf,nbf)
      common /io/ inp, iout
c**********************************************************************c
c                   orthonormalize the bound orbitals                  c
c                         maintain symmetry                            c
c**********************************************************************c
      if (orth) then
          call schmdt(fns,ddfns,fns,ddfns,mask,eig,eig,n,nbf,
     1                nout,'real')
          call iosys ('write real "bound molecular basis functions" '//
     1                'to atomci',n*nout,fns,0,' ')
          call iosys ('write real "kinetic energy of bound molecular '//
     1                'basis functions" to atomci',n*nout,ddfns,0,' ')
      else
          call iosys ('read integer "bound symmetry list" from atomci',
     1                 dimsym,nsym,0,' ')
          count=0
          do 10 i=1,dimsym
             if (nsym(i).ne.0) then
                 countp=count+1
                 call tred2(nbf,nsym(i),s(countp,countp),eig(countp),
     1                                  work,s(countp,countp))
                 call tql2(nbf,nsym(i),eig(countp),work,
     1                                 s(countp,countp),ierr) 
                 count=count+nsym(i)
             endif
   10     continue      
          if (prnt) then
              write(iout,20) (eig(i),i=1,nbf)
          endif
          nout=0
          do 30 i=1,nbf      
             if (abs(eig(i)).gt.tol) then
                 nout=nout+1
                 eig(i)=1.d0/sqrt(eig(i))
                 lb(nout)=lb(i)
                 mb(nout)=mb(i)
                 do 40 j=1,nbf
                    s(j,nout)=eig(i)*s(j,i)
   40            continue
             endif
   30     continue
          if (prnt) then
              write(iout,50) tol, nout
          endif
          call iosys('write real "bound transformation matrix" to '//
     1               'atomci',nbf*nout,s,0,' ')
c**********************************************************************c
c             transform orbitals and derivatives to orthogonal         c
c                            basis                                     c
c**********************************************************************c
          call ebc(fmo,fns,s,n,nbf,nout)
          call iosys ('write real "bound molecular basis functions" '//
     1                'to atomci',n*nout,fmo,0,' ')
          call copy(fmo,fns,n*nout)
          call ebc(fmo,ddfns,s,n,nbf,nout)
          call iosys ('write real "kinetic energy of bound molecular '//
     1                'basis functions" to atomci',n*nout,fmo,0,' ')
          call copy(fmo,ddfns,n*nout) 
      endif
      call iosys('write integer "number of orthogonal bound '//
     1           'functions" to atomci',1,nout,0,' ')
      return
   20 format(//,5x,'eigenvalues of bound overlap matrix',
     1              (/,5x,5(e15.8,1x)))
   50 format(//,5x,'number of eigenvectors above tolerance of',1x,
     1              e15.8,1x,'is',1x,i3)  
      end







