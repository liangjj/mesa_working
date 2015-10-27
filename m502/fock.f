*deck  @(#)fock.f	5.1 11/6/94
      subroutine fock (values,d,nnp,num,jmat,kmat,ncoul,nexch,ntriang
     1                 ,isq,ndmat,dsq)
c
c***module to form the two-electron contribution to the fock matrix
c   for general rhf using an ordered integral list.
c
c paul saxe                  6 september 1984                    lanl
c
      implicit integer(a-z)
c
      real*8 values(nnp,ntriang), jmat(nnp,ncoul)
      real*8 kmat(nnp,nexch), dsq(num,num,nexch)
      real*8 isq(num,num), djl, djk, dkl
      real*8 d(nnp,ndmat)
      real*8 zero,one
      logical read
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
      data call /0/
c
      save  call
c
      call=call+1
      call iosys ('rewind "sorted ao integrals" on ints',0,0,0,' ')
      call rzero (jmat,nnp*ncoul)
      call rzero (kmat,nnp*nexch)
      minex=1
      maxex=ndmat
c
c     ----- square up the density matrices -----
c
      do 10 i=minex,maxex
         call trtosq (dsq(1,1,i-minex+1),d(1,i),num,nnp)
   10 continue
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
      endif
c
      kl=0
      do 50 k=1,num
         k1=k*(k-1)/2+1
         do 40 l=1,k
            l1=l*(l-1)/2+1
            kl=kl+1
c
c     ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnp,maxkl+ntriang)
               lnread=(maxkl-minkl+1)*nnp
               call iosys ('read real "sorted ao integrals"'
     1                   //' from ints without rewinding',
     2                     lnread,values,0,' ')
            endif
c
c    ----- form coulomb terms -----
c
            if (k.eq.l.and.read) then
               call sscal (nnp,0.5d+00,values(1,kl-minkl+1),1)
            endif
c           do 1 type=1,ncoul
c              dkl=d(kl,type)*2.0d+00
c              call saxpy(nnp,dkl,values(1,kl-minkl+1),1,
c    #c         c         c                jmat(1,type),1)
c   1       continue
            if (kl.eq.maxkl) then
               do 20 type=1,ncoul
                  call sgemv ('n',nnp,maxkl-minkl+1,one,values,nnp,
     $                        d(minkl,type),1,one,jmat(1,type),1)
   20          continue
            endif
c
c    ----- form exchange terms -----
c
            call trtosq (isq,values(1,kl-minkl+1),num,nnp)
c
cc            do 4 j=1,num
cc               junk=max(j,l)
cc               jl=junk*(junk-1)/2+min(j,l)
cc               junk=max(j,k)
cc               jk=junk*(junk-1)/2+min(j,k)
cc               do 2 type=minex,maxex
cc                  djl=d(jl,type)
cc                  djk=d(jk,type)
cc                  junk=type-minex+1
cc                  call saxpy(k,djl,isq(1,j),1,kmat(k1,junk),1)
cc                  call saxpy(l,djk,isq(1,j),1,kmat(l1,junk),1)
cc    2          continue
cc    4       continue
            do 30 type=1,nexch
               call sgemv ('n',k,num,one,isq,num,dsq(1,l,type),1,
     $                     one,kmat(k1,type),1)
               call sgemv ('n',l,num,one,isq,num,dsq(1,k,type),1,
     $                     one,kmat(l1,type),1)
   30       continue
   40    continue
   50 continue
c
c     ----- double coulomb terms to account for [ij;kl]
c                   and [ij;lk] ------
c
      do 60 type=1,ncoul
         call sscal (nnp,2.0d+00,jmat(1,type),1)
   60 continue
c
      return
      end
