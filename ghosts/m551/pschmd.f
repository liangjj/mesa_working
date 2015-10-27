*deck %W%  %G%
      subroutine pschmd(c,stri,s,t2,t3,nfunc,num,nnp,maxerr)
c***begin prologue     schmdt
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           gram-schmdt, orthogonalization
c***author             saxe, paul (lanl)
c***source             util
c***purpose            vectorized gram-schmidt orthogonalization.
c***description
c                      call schmdt(c,stri,s,t2,t3,nfunc,num,nnp,maxerr)
c                        c        input matrix of eigenvectors (num,nfunc).
c                                 on output, the vectors are orthonormal
c                                 on csc=i.
c                        stri     scratch (nnp).
c                        s        the overlap matrix (num,num).
c                        t2       scratch (num,num).
c                        t3       scratch (nfunc,nfunc).
c                        nfunc    the number of vectors.
c                        num      the dimension of the basis.
c                        nnp      num*(num+1)/2
c                        maxerr   the maximum acceptable deviation
c                                 from orthogonality. typically 10**(-6).
c
c***references
c***routines called    trtosq(math),saxpy(clams), sgemv(clams),sscal(clams),
c                      ebc(math),ebct(math),max,lnkerr(mdutil)
c***end prologue       schmdt
      implicit integer (a-z)
c
      real*8 c(num,nfunc),stri(nnp),s(num,num),t2(num,num)
      real*8 t3(nfunc,nfunc)
      real*8 sdot,maxerr,err
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
      common /io/     inp,iout
c
 1000 format(' schmdt finds largest deviation from orthogonality',e10.1,
     $     '      maximum acceptable is ',e10.1)
c
c
      call trtosq(s,stri,num,nnp)
c
c     ----- loop over orbitals -----
c
      do 10 vec=1,nfunc
c
c     ----- orthogonalize to previous vectors -----
c
         do 5 old=1,vec-1
            call saxpy(num,-sdot(num,t2(1,old),1,c(1,vec),1),c(1,old),
     $           1,c(1,vec),1)
    5    continue
c
c     ----- normalize -----
c
         call sgemv('n',num,num,one,s,num,c(1,vec),1,zero,t2(1,vec),1)
         call sscal(num,1.0/sqrt(sdot(num,t2(1,vec),1,c(1,vec),1)),
     $        c(1,vec),1,c(1,vec),1)
c
c     ----- put c(vec)s in t2 -----
c
         call sgemv('n',num,num,one,s,num,c(1,vec),1,zero,t2(1,vec),1)
 10   continue
c
c     ----- check that we've done it correctly -----
c
      call ebc(t2,s,c,num,num,nfunc)
      call ebtc(t3,c,t2,nfunc,num,nfunc)
      err=0.0d+00
      do 20 j=1,nfunc
         do 15 i=1,j-1
            err=max(err,abs(t3(i,j)))
 15      continue
         err=max(err,abs(t3(j,j)-1.0e00))
 20   continue
c
c
      if(err.gt.maxerr) then
         write(iout,1000) err,maxerr
         call lnkerr(' non orthogonal vectors.')
      endif
c
c
      return
      end
