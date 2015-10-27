*deck @(#)lagrng.f	5.1  11/6/94
      subroutine lagrng(numshl,minshl,maxshl,f,alpha,beta,h,c,d,jmat,
     #                  kmat,values,num,nnp,nshell,ndmat,ncoul,
     #                  nexch,scfnum,ntriang,t1,t2,lag,ops,der,grad,
     #                  nder,clag)
c
c***begin prologue     lagrng
c***date written       860819  (yymmdd)
c***revision date      871115  (yymmdd)
c
c   15 november 1987   pws at lanl
c       forming the derivative core lagrangians also.
c
c    8 december 1986   pws at lanl
c        modifying for the calculation of 'derivative' lagrangians using
c        ao derivative integrals.
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)lagrng.f	5.1   11/6/94
c***purpose            formation of the general scf lagrangian
c***description
c
c***references
c***routines called    fock     (m1010)
c***end prologue       lagrng
c
      implicit integer (a-z)
c
      character*(*) ops
      integer numshl(nshell),minshl(nshell),maxshl(nshell)
      logical logkey
      real*8 f(nshell),alpha(nshell,nshell)
      real*8 beta(nshell,nshell),h(nnp),c(num,num)
      real*8 d(nnp,ndmat),jmat(nnp,ncoul),kmat(nnp,ncoul)
      real*8 t1(num*num),t2(num,num),lag(num,num),values(nnp,ntriang)
      real*8 clag(num,num)
      real*8 grad(nder)
c
      common /io/ inp,iout
c
c     ----- get the derivative one-electron integrals -----
c
      call iosys('read real "ao derivative one-electron integrals" '//
     $     'from rdints',nnp,h,(der-1)*nnp,' ')
c
c     ----- get the scf vector -----
c
      call iosys('read real "scf vector" from rwf',num**2,c,0,' ')
c
c     ----- transform the one-electron integrals to the mo basis -----
c
      call trtosq(t1,h,num,nnp)
      call ebc(t2,t1,c,num,num,num)
      call ebtc(t1,c,t2,num,num,num)
      call sqtotr(h,t1,num,nnp)
c
c     ----- get the density matrices for each orbital type -----
c
      call iosys('read real "hf density matrix" from rwf',nnp*ndmat,
     #            d,0,' ')
c
c     ----- form the coulomb and exchange matrices -----
c
      call derfoc(values,d,values,nnp,num,jmat,kmat,ncoul,nexch,ntriang,
     #          -nexch,1,t1,h,ndmat,t2,der)
c
c     ----- transform to the mo basis and form the lagrangian as we go
c       l(ij)=f(i)h(ij) + sum(l=occ) {alpha(il)[ij;ll]+beta(il)[il;jl]}
c      cl(ij)=h(ij) + sum(l=occ) f(i){2[ij;ll]-[il;jl]}
c
      call rzero(lag,num**2)
      call rzero(clag,num**2)
c
cps      do 100 ishl=1,nshell-1
      do 100 ishl=1,nshell
         mini=minshl(ishl)
         numi=numshl(ishl)
c..bhl
         if(numi.eq.0) go to 100
c..bhl
c
c        ----- the one-electron term -----
c
         call trtosq(t2,h,num,nnp)
         do 5 i=mini,maxshl(ishl)
            do 4 j=1,maxshl(nshell)
               lag(i,j)=lag(i,j)+f(ishl)*t2(i,j)
               clag(i,j)=clag(i,j)+t2(i,j)
    4       continue
c
c           ----- one-electron term of scf gradient -----
c
            if (ishl.lt.nshell) grad(der)=grad(der)+lag(i,i)
    5    continue
c
         do 80 lshl=1,nshell-1
c
c           ----- transform j(ll) and k(ll) over i and j orbital ranges
c
            call trtosq(t1,jmat(1,lshl),num,nnp)
            call ebc(t2,t1,c,num,num,num)
            call ebtc(t1,c(1,mini),t2,numi,num,num)
c
c           ----- put j in correct place -----
c
            do 10 i=1,numi
               ii=i+mini-1
               do 9 j=1,maxshl(nshell)
                  ij=(j-1)*numi+i
                  lag(ii,j)=lag(ii,j)+alpha(ishl,lshl)*t1(ij)
                  clag(ii,j)=clag(ii,j)+2*f(lshl)*t1(ij)
    9          continue
   10       continue
c
c           ----- exchange (k) part -----
c
            call trtosq(t1,kmat(1,lshl),num,nnp)
            call ebc(t2,t1,c,num,num,num)
            call ebtc(t1,c(1,mini),t2,numi,num,num)
c
c           ----- put k in correct place -----
c
            do 20 i=1,numi
               ii=i+mini-1
               do 19 j=1,maxshl(nshell)
                  ij=(j-1)*numi+i
                  lag(ii,j)=lag(ii,j)+beta(ishl,lshl)*t1(ij)
                  clag(ii,j)=clag(ii,j)-f(lshl)*t1(ij)
   19          continue
   20       continue
   80    continue
c
c        ----- lagrangian term of scf gradient -----
c
         if (ishl.lt.nshell) then
            do 85 i=mini,maxshl(ishl)
               grad(der)=grad(der)+lag(i,i)
 85         continue
         end if
  100 continue
c
      if (logkey(ops,'print=derivative=mo-lagrangian',.false.,' '))
     $     then
         write (iout,120) der
  120    format('1',//,t20,'mo derivative hartree-fock lagrangian',i4)
         call matout(lag,num,num,num,num,iout)
      end if
c
      call iosys('write real "mo derivative hf lagrangian" to rwf '//
     #           'without rewinding',num**2,lag,0,' ')
c
      if (logkey(ops,'print=derivative=mo-core-lagrangian',.false.,
     $     ' ')) then
         write (iout,121) der
 121     format('1',//,t20,'mo derivative core lagrangian',i4)
         call matout(clag,num,num,num,num,iout)
      end if
c
      call iosys('write real "mo derivative core lagrangian" to rwf '//
     #           'without rewinding',num**2,clag,0,' ')
c
c     ----- transform to the ao basis -----
c
      call ebct(t1,lag,c,num,num,num)
      call ebc(t2,c,t1,num,num,num)
c
      if (logkey(ops,'print=derivative=ao-lagrangian',.false.,' '))
     $     then
         write (iout,130) der
  130    format ('1',//,t20,'ao derivative hartree-fock lagrangian',i4)
         call matout(t2,num,num,num,num,iout)
      end if
c
      call iosys('write real "ao derivative hf lagrangian" to rwf '//
     #           'without rewinding',num**2,t1,0,' ')
c
c
      return
      end
