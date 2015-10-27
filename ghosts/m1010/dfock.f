*deck @(#)dfock.f	1.1  11/30/90
      subroutine dfock(values,d,f,nnp,num,jmat,kmat,ncoul,
     #                  nexch,ntriang,scfnum,nfock,isq,h,ndmat,
     #                  dsq,der)
c
c***module to form the two-electron contribution to the fock matrix
c   for general rhf using an ordered integral list.
c
c paul saxe                  6 september 1984                    lanl
c
c   8 december 1986   pws at lanl
c      modifying to work with derivative integrals
c
      implicit integer (a-z)
c
      real*8 values(nnp,ntriang),f(nnp,nfock),jmat(nnp,ncoul)
      real*8 kmat(nnp,nexch),h(nnp),dsq(num,num,nexch)
      real*8 isq(num,num),djl,djk,dkl,dil,dik
      real*8 d(nnp,ndmat)
      real*8 sdot
      logical read
c
      data call /0/
c
      save start,call
c
c
      call=call+1
      call iosys('rewind "sorted ao derivative integrals" on ints',
     $     0,0,0,' ')
      offset=(der-1)*nnp**2
c
c     ----- make sure we read the derivative integrals -----
c
      call=1
c
      call rzero(jmat,nnp*ncoul)
      call rzero(kmat,nnp*nexch)
      if (scfnum.eq.0) then
         minex=1
         maxex=1
      else if (scfnum.eq.1) then
         minex=2
         maxex=3
      else if (scfnum.lt.0) then
         minex=1
         maxex=-scfnum
      end if
c
c     ----- square up the density matrices -----
c
      do 60 i=minex,maxex
         call trtosq(dsq(1,1,i-minex+1),d(1,i),num,nnp)
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
      do 6 k=1,num
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
               call iosys('read real "sorted ao derivative integrals" '
     $              //'from ints without rewinding',
     $              lnread,values,offset,' ')
               offset=0
            end if
c
c    ----- form coulomb terms -----
c
            if (k.eq.l.and.read) then
               call sscal(nnp,0.5d+00,values(1,kl-minkl+1),1)
            end if
c
            if (kl.eq.maxkl) then
               do 1 type=1,ncoul
               call sgemv(nnp,maxkl-minkl+1,values,nnp,
     #                    d(minkl,type),1,jmat(1,type),1,2)
    1          continue
            end if
c
c    ----- form exchange terms -----
c
            call trtosq(isq,values(1,kl-minkl+1),num,nnp)
c
            do 2 type=1,nexch
               call sgemv(k,num,isq,num,dsq(1,l,type),1,kmat(k1,type),
     #                                                            1,2)
               call sgemv(l,num,isq,num,dsq(1,k,type),1,kmat(l1,type),
     #                                                            1,2)
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
c     ----- combine coulomb and exchange to form fock-matrix -----
c
      if (scfnum.eq.0) then
c
c     ----- closed shell scf -----
c
         do 21 i=1,nnp
            f(i,1)=h(i)+2.0d+00*jmat(i,1)-kmat(i,1)
   21    continue
      else if (scfnum.eq.1) then
c
c     ----- high-spin open-shell scf -----
c
         do 22 i=1,nnp
            f(i,1)=h(i)+2.0d+00*jmat(i,1)-kmat(i,1)-0.5d+00*kmat(i,2)
            f(i,2)=0.5d+00*(h(i)+2.0*jmat(i,1)-kmat(i,1)-kmat(i,2))
   22    continue
c
      end if
c
c
      return
      end
