*deck @(#)lagrng.f	5.1  11/6/94
      subroutine lagrng(numshl,minshl,maxshl,f,alpha,beta,h,c,d,jmat,
     #                  kmat,values,nbf,nnp,nshell,ndmat,ncoul,
     #                  nexch,ntriang,t1,t2,lag,ops)
c
c***begin prologue     lagrng
c***date written       860819  (yymmdd)
c***revision date      861209  (yymmdd)
c
c   9 december 1986   pws at lanl
c       modified the normal lagrangian routine to calculate the
c       generalised lagrangians.
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)lagrng.f	5.1   11/6/94
c***purpose            formation of the generalised scf lagrangian
c***description
c
c***references
c***routines called    jandks     (m1011)
c***end prologue       lagrng
c
      implicit integer (a-z)
c
      character*(*) ops
      character*2 ians
      integer numshl(nshell),minshl(nshell),maxshl(nshell)
      logical logkey
      real*8 f(nshell),alpha(nshell,nshell)
      real*8 beta(nshell,nshell),h(nnp),c(nbf,nbf)
      real*8 d(nnp,ndmat),jmat(nnp,ncoul),kmat(nnp,ncoul)
      real*8 t1(nbf*nbf),t2(nbf,nbf),lag(nbf,nbf),values(nnp,ntriang)
c
      common /io/ inp,iout
c
c     ----- get the scf vector -----
c
      call iosys('read real "scf vector" from rwf',nbf**2,c,0,' ')
c
c     ----- transform the one-electron integrals to the mo basis -----
c
      call trtosq(t1,h,nbf,nnp)
      call ebc(t2,t1,c,nbf,nbf,nbf)
      call ebtc(t1,c,t2,nbf,nbf,nbf)
      call sqtotr(h,t1,nbf,nnp)
c
c     ----- form the density matrices for each orbital type -----
c
      do 1 dmat=1,ndmat
         call gdmat(d(1,dmat),c,nbf,nnp,minshl(dmat),maxshl(dmat))
    1 continue
c
c     ----- form the coulomb and exchange matrices -----
c
      call jandks(values,d,nnp,nbf,jmat,kmat,ncoul,nexch,ntriang,
     #            t1,ndmat,t2)
c
c     ----- transform to the mo basis and form the generalized lagrangian
c            as we go
c
c           l(ij,k)=f(k)h(ij) + sum(l=occ) {alpha(kl)[ij;ll]+beta(kl)[il;jl]}
c
      call iosys('does "generalised hf lagrangian" exist on rwf',
     #            0,0,0,ians)
      if(ians.eq.'no') then
         call iosys('create real "generalised hf lagrangian" on rwf',
     #               nbf**2*(nshell-1),0,0,' ')
      else
         call iosys('rewind "generalised hf lagrangian" on rwf',
     #               0,0,0,' ')
      end if
c
      do 1000 kshl=1,nshell-1
         call rzero(lag,nbf**2)
c
         do 100 ishl=1,nshell
            mini=minshl(ishl)
            nbfi=numshl(ishl)
c
c           ----- the one-electron term -----
c
            call trtosq(t2,h,nbf,nnp)
            do 5 i=mini,maxshl(ishl)
               do 4 j=1,nbf
                  lag(i,j)=lag(i,j)+f(kshl)*t2(i,j)
    4          continue
    5       continue
c
            do 80 lshl=1,nshell-1
c
c              ----- transform j(ll) and k(ll) over i and j orbital ranges
c
               call trtosq(t1,jmat(1,lshl),nbf,nnp)
               call ebc(t2,t1,c,nbf,nbf,nbf)
               call ebtc(t1,c(1,mini),t2,nbfi,nbf,nbf)
c
c              ----- put j in correct place -----
c
               do 10 i=1,nbfi
                  ii=i+mini-1
                  do 9 j=1,nbf
                     ij=(j-1)*nbfi+i
                     lag(ii,j)=lag(ii,j)+alpha(kshl,lshl)*t1(ij)
    9             continue
   10          continue
c
c              ----- exchange (k) part -----
c
               call trtosq(t1,kmat(1,lshl),nbf,nnp)
               call ebc(t2,t1,c,nbf,nbf,nbf)
               call ebtc(t1,c(1,mini),t2,nbfi,nbf,nbf)
c
c              ----- put k in correct place -----
c
               do 20 i=1,nbfi
                  ii=i+mini-1
                  do 19 j=1,nbf
                     ij=(j-1)*nbfi+i
                     lag(ii,j)=lag(ii,j)+beta(kshl,lshl)*t1(ij)
   19             continue
   20          continue
   80       continue
  100    continue
c
         if (logkey(ops,'print=generalised-lagrangians',.false.,' '))
     $        then
            write (iout,120) kshl
  120       format('1',//,t20,'mo generalised lagrangian',i5,/)
            call matout(lag,nbf,nbf,nbf,nbf,iout)
         end if
c
         call iosys('write real "generalised hf lagrangian" to rwf '//
     #              'without rewinding',nbf**2,lag,0,' ')
c
 1000 continue
c
c
      return
      end
