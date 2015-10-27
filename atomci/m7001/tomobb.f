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
      real *8 fns, ddfns, fmo, s, eig, work, tol, mask
      logical prnt, orth
      dimension fns(n,nbf), ddfns(n,nbf), s(nbf,nbf), work(nbf)
      dimension eig(nbf), nsym(dimsym), fmo(n,nbf), lb(nbf), mb(nbf)
      dimension mask(nbf,nbf)
      common /io/ inp, iout
c**********************************************************************c
c                   orthonormalize the bound orbitals                  c
c                         maintain symmetry                            c
c**********************************************************************c
      call iosys ('read integer "bound symmetry list" from atomci',
     1             dimsym,nsym,0,' ')
      count=1
      do 10 i=1,dimsym
         if (nsym(i).ne.0) then
c**********************************************************************c
c                use schmidt procedure                                 c
c**********************************************************************c
             if (orth) then
                 call schmdt(fns(1,count),ddfns(1,count),
     1                       fns(1,count),ddfns(1,count),
     2                       mask(count,count),eig,eig,work(count),
     3                       n,nsym(i),'real')
             else
c**********************************************************************c
c                     diagonalize the overlap                          c
c**********************************************************************c
                 call tred2(nbf,nsym(i),s(count,count),
     1                      eig(count),work,s(count,count))
                 call tql2(nbf,nsym(i),eig(count),work,
     1                                 s(count,count),ierr) 
             endif
             count=count+nsym(i)
         endif
   10 continue
      if (orth) then
          call vsarr(fns,fns,ddfns,ddfns,work,nsym,lb,mb,nout,nbf,
     1               n,dimsym,'real')
      else
          call vearr(s,s,eig,eig,tol,nsym,lb,mb,nout,nbf,
     1               dimsym,'real')
      endif 
      if (orth) then
          if (nout.ne.0) then
              call iosys ('write real "bound molecular basis '//
     1                    'functions" to atomci',n*nout,fns,0,' ')
              call iosys ('write real "kinetic energy of bound '//
     1                    'molecular basis functions" to atomci',
     2                     n*nout,ddfns,0,' ')
          endif
      else
          if (prnt) then
              write(iout,70) tol, nout
              write(iout,80) (eig(i),i=1,nbf)
          endif
          if (nout.ne.0) then
              call iosys('write real "bound transformation matrix" '//
     1                   'to atomci',nbf*nout,s,0,' ')
c**********************************************************************c
c             transform orbitals and derivatives to orthogonal         c
c                            basis                                     c
c**********************************************************************c
              call ebc(fmo,fns,s,n,nbf,nout)
              call iosys ('write real "bound molecular basis '//
     1                    'functions" to atomci',n*nout,fmo,0,' ')
              call copy(fmo,fns,n*nout)
              call ebc(fmo,ddfns,s,n,nbf,nout)
              call iosys ('write real "kinetic energy of bound '//
     1                    'molecular basis functions" to atomci',
     2                     n*nout,fmo,0,' ')
              call copy(fmo,ddfns,n*nout) 
          endif
      endif
      call iosys('write integer "number of orthogonal bound '//
     1           'functions" to atomci',1,nout,0,' ')
       call iosys('write integer "l values of orthogonal bound '//
     1            'orbitals" to atomci',nout,lb,0,' ')
       call iosys('write integer "m values of orthogonal bound '//
     1            'orbitals" to atomci',nout,mb,0,' ')
       call iosys ('write integer "bound orthogonal symmetry '//
     1             'list" to atomci',dimsym,nsym,0,' ')
      return
   70 format(//,5x,'number of eigenvectors above tolerance of',1x,
     1              e15.8,1x,'is',1x,i3)  
   80 format(//,5x,'eigenvalues of bound overlap matrix',
     1              (/,5x,5(e15.8,1x)))
      end
