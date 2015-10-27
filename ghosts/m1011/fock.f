*deck %W%  %G%
      subroutine fock(values,d,nnp,nbf,jmat,kmat,ncoul,
     #                  nexch,ntriang,nfock,isq,h,ndmat,
     #                  dsq)
c
c***module to form the two-electron coulomb and exchange matrices 
c   for general rhf using an ordered integral list.
c
c paul saxe                  6 september 1984                    lanl
c
      implicit integer (a-z)
c
      real*8 values(nnp,ntriang),jmat(nnp,ncoul)
      real*8 kmat(nnp,nexch),h(nnp),dsq(nbf,nbf,nexch)
      real*8 isq(nbf,nbf)
      real*8 d(nnp,ndmat)
      real*8 zero,one
c
      logical read
c
      parameter (zero=0.0d+00,one=1.0d+00)
c
      data call /0/
c
      save call
c
c
      call=call+1
      call iosys('rewind "sorted ao integrals" on ints',0,0,0,' ')
      call rzero(jmat,nnp*ncoul)
      call rzero(kmat,nnp*nexch)
c
c     ----- square up the density matrices -----
c
      do 60 i=1,nexch
         call trtosq(dsq(1,1,i),d(1,i),nbf,nnp)
   60 continue
c
c     ----- read through integrals if first iteration, or if
c           cannot hold all the integrals in core -----
c
      if (call.eq.1.or.ntriang.lt.nnp) then
         maxkl=0
         read=.true.
      else
         minkl=1
         maxkl=nnp
         read=.false.
      end if
c
      kl=0
      do 6 k=1,nbf
         k1=k*(k-1)/2+1
         do 5 l=1,k
            l1=l*(l-1)/2+1
            kl=kl+1
c
c     ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnp,maxkl+ntriang)
               lnread=(maxkl-minkl+1)*nnp
               call iosys('read real "sorted ao integrals" from ints '
     #              //'without rewinding',lnread,values,0,' ')
            end if
c
c    ----- form coulomb terms -----
c
            if (k.eq.l.and.read) then
               call sscal(nnp,0.5d+00,values(1,kl-minkl+1),1)
            end if
            if (kl.eq.maxkl) then
               do 1 type=1,ncoul
               call sgemv('n',nnp,maxkl-minkl+1,one,values,nnp,
     #                    d(minkl,type),1,one,jmat(1,type),1)
    1          continue
            end if
c
c    ----- form exchange terms -----
c
            call trtosq(isq,values(1,kl-minkl+1),nbf,nnp)
c
            do 2 type=1,nexch
               call sgemv('n',k,nbf,one,isq,nbf,dsq(1,l,type),1,
     $                    one,kmat(k1,type),1)
               call sgemv('n',l,nbf,one,isq,nbf,dsq(1,k,type),1,
     $                    one,kmat(l1,type),1)
    2       continue
    5    continue
    6 continue
c
c     ----- double coulomb terms to account for [ij;kl]
c                   and [ij;lk] ------
c
      do 7 type=1,ncoul
         call sscal(nnp,2.0d+00,jmat(1,type),1)
    7 continue
c
c
      return
      end
