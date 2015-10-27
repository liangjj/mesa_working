*deck @(#)tomocc.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001, spline
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            transform complex functions to orthonormal basis 
c***                   which is orthogonal to bound subspace.
c***
c***description        the complex functions are first schmidt orthogonalized
c***                   to the bound orbitals. then the remaining set is
c***                   orthonormalized by diagonalizing the overlap matix.
c***                   near linearly dependent vectors are discarded and 
c***                   the functions transformed.
c***                 
c***                 
c***references       
c
c***routines called    cgeev(clams)
c***end prologue       m7001
      subroutine tomocc(fns,ddfns,fnsc,ddfnsc,fmo,lb,mb,lc,mc,mask,
     1                  s1,s,eig,v,work,tol,nsym,n,nbfb,nbfc,nout,
     2                  dimsym,prnt,orth)
      implicit integer (a-z)
      real *8 fns, ddfns, mask, tol
      complex *16 fnsc, ddfnsc, fmo, s, s1, eig, work, v, tmp
      logical prnt, orth
      dimension fns(n,nbfb), ddfns(n,nbfb), fnsc(n,nbfc), s(nbfc,nbfc)
      dimension fmo(n,nbfc), lb(nbfb), mb(nbfb), lc(nbfc), mc(nbfc)
      dimension eig(nbfc), work(3*nbfc), v(nbfc,nbfc), s1(nbfb,nbfc)
      dimension  nsym(dimsym), ddfnsc(n,nbfc), mask(*)
      common /io/ inp, iout
c**********************************************************************c
c                schmidt the complex to the bound orbitals             c
c**********************************************************************c
      if (nbfb.ne.0) then
          call mskone(mask,lb,mb,nbfb,lc,mc,nbfc)
          call ebtcc(s1,fns,fnsc,nbfb,n,nbfc)
          call cvmul(s1,s1,mask,mask,nbfb*nbfc,'real')
          call ambcc(fnsc,fns,s1,n,nbfb,nbfc)
          call ambcc(ddfnsc,ddfns,s1,n,nbfb,nbfc)
      endif
c**********************************************************************c
c                calculate complex overlap                             c
c**********************************************************************c
      call mskone(mask,lc,mc,nbfc,lc,mc,nbfc)
      call cebtc(s,fnsc,fnsc,nbfc,n,nbfc)
      call cvmul(s,s,mask,mask,nbfc*nbfc,'real')  
      call iosys ('read integer "complex symmetry list" from atomci',
     1             dimsym,nsym,0,' ')
      if (orth) then
c**********************************************************************c
c               schmidt orthogonalize the complex functions            c
c**********************************************************************c
          call schmdt(fnsc,ddfnsc,fnsc,ddfnsc,mask,eig,eig,n,nbfc,
     1                nout,'complex')
          call iosys ('write real "complex molecular basis '//
     1                'functions" to atomci',2*n*nout,fnsc,0,' ')
          call iosys ('write real "kinetic energy of complex '//
     1                'molecular basis functions" to atomci',2*n*nout,
     2                 ddfnsc,0,' ')
      else
c**********************************************************************c
c             diagonalize the overlap of the complex functions         c
c**********************************************************************c
          count=0
          do 10 i=1,dimsym
             if (nsym(i).ne.0) then
                 countp=count+1
                 call cgeev(s(countp,countp),nbfc,nsym(i),eig(countp),
     1                      v(countp,countp),nbfc,work,1,info)
                 count=count+nsym(i)
             endif
   10     continue      
          if (prnt) then
              write(iout,20) (eig(i),i=1,nbfc)
          endif
          nout=0
          do 30 i=1,nbfc      
             if (abs(eig(i)).gt.tol) then
                 nout=nout+1
                 lc(nout)=lc(i)
                 mc(nout)=mc(i)
                 tmp=(0.d0,0.d0)
                 do 40 j=1,nbfc
                    tmp=tmp+v(j,i)*v(j,i)
   40            continue
                 tmp=1.d0/sqrt(tmp*eig(i))     
                 do 45 j=1,nbfc
                    v(j,nout)=tmp*v(j,i)
   45            continue
             endif
   30     continue    
          if (prnt) then
              write(iout,50) tol, nout
          endif
          call iosys('write real "complex transformation matrix" to '//
     1               'atomci',2*nbfc*nout,v,0,' ')
c**********************************************************************c
c             transform orbitals and derivatives to orthogonal         c
c                            basis                                     c
c**********************************************************************c
          call cebc(fmo,fnsc,v,n,nbfc,nout)
          call iosys ('write real "complex molecular basis '//
     1                'functions" to atomci',2*n*nout,fmo,0,' ')
          call copy(fmo,fnsc,2*n*nout) 
          call cebc(fmo,ddfnsc,v,n,nbfc,nout)
          call iosys ('write real "kinetic energy of complex '//
     1                'molecular basis functions" to atomci',2*n*nout,
     2                 fmo,0,' ')
          call copy(fmo,ddfnsc,2*n*nout) 
      endif
      call iosys('write integer "number of orthogonal complex '//
     1           'functions" to atomci',1,nout,0,' ')
      return
   20 format(//,5x,'eigenvalues of complex overlap matrix',
     1              (/,5x,'(',e15.8,',',e15.8,')',1x,
     2                    '(',e15.8,',',e15.8,')'))
   50 format(//,5x,'number of eigenvectors above tolerance of',1x,
     1              e15.8,1x,'is',1x,i3)  
      end







