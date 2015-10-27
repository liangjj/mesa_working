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
     1                  s1,s,eig,v,work,work1,tol,nsym,n,nbfb,nbfc,
     2                  nout,dimsym,prnt,orth)
      implicit integer (a-z)
      real *8 fns, ddfns, mask, tol, work1
      complex *16 fnsc, ddfnsc, fmo, s, s1, eig, work, v
      logical prnt, orth
      dimension fns(n,nbfb), ddfns(n,nbfb), fnsc(n,nbfc), s(nbfc,nbfc)
      dimension fmo(n,nbfc), lb(nbfb), mb(nbfb), lc(nbfc), mc(nbfc)
      dimension eig(nbfc), work(3*nbfc), v(nbfc,nbfc), s1(nbfb,nbfc)
      dimension  nsym(dimsym), ddfnsc(n,nbfc), mask(nbfc,nbfc)
      dimension work1(3*nbfc)
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
      count=1
      do 10 i=1,dimsym
         if (nsym(i).ne.0) then
             if (orth) then
c**********************************************************************c
c                 use schmidt procedure                                c
c**********************************************************************c
                 call schmdt(fnsc(1,count),ddfnsc(1,count),
     1                       fnsc(1,count),ddfnsc(1,count),
     2                       mask(count,count),eig,eig,work1(count),
     3                       n,nsym(i),'complex')
             else
                 call cgeev(s(count,count),nbfc,nsym(i),eig(count),
     1                      v(count,count),nbfc,work,1,info)
             endif
             count=count+nsym(i)
         endif
   10 continue      
      if (orth) then
          call vsarr(fnsc,fnsc,ddfnsc,ddfnsc,work1,nsym,lc,mc,nout,
     1               nbfc,n,dimsym,'complex')
      else
          call vearr(v,v,eig,eig,tol,nsym,lc,mc,nout,nbfc,
     1               dimsym,'complex')
      endif 
      if (orth) then 
          if (nout.ne.0) then
             call iosys ('write real "complex molecular basis '//
     1                   'functions" to atomci',2*n*nout,fnsc,0,' ')
             call iosys ('write real "kinetic energy of complex '//
     1                   'molecular basis functions" to atomci',
     2                    2*n*nout,ddfnsc,0,' ')
          endif 
      else
          if (prnt) then
              write(iout,50) tol, nout
          endif
          if (nout.ne.0) then
             call iosys('write real "complex transformation matrix" '//
     1                  'to atomci',2*nbfc*nout,v,0,' ')
c**********************************************************************c
c             transform orbitals and derivatives to orthogonal         c
c                            basis                                     c
c**********************************************************************c
             call cebc(fmo,fnsc,v,n,nbfc,nout)
             call iosys ('write real "complex molecular basis '//
     1                   'functions" to atomci',2*n*nout,fmo,0,' ')
             call copy(fmo,fnsc,2*n*nout) 
             call cebc(fmo,ddfnsc,v,n,nbfc,nout)
             call iosys ('write real "kinetic energy of complex '//
     1                   'molecular basis functions" to atomci',
     2                   2*n*nout,fmo,0,' ')
             call copy(fmo,ddfnsc,2*n*nout) 
          endif
      endif
      call iosys('write integer "number of orthogonal complex '//
     1           'functions" to atomci',1,nout,0,' ')
       call iosys('write integer "l values of orthogonal complex '//
     1            'orbitals" to atomci',nout,lc,0,' ')
       call iosys('write integer "m values of orthogonal complex '//
     1            'orbitals" to atomci',nout,mc,0,' ')
       call iosys ('write integer "complex orthogonal symmetry '//
     1             'list" to atomci',dimsym,nsym,0,' ')
      return
   20 format(//,25x,'symmetry block',1x,i3)
   25 format(//,5x,'eigenvalues of complex overlap matrix',
     1              (/,5x,:'(',e15.8,',',e15.8,')',1x,
     2                    :'(',e15.8,',',e15.8,')'))
   50 format(//,5x,'number of eigenvectors above tolerance of',1x,
     1              e15.8,1x,'is',1x,i3)  
      end







