*deck %W%  %G%
      subroutine lagrng(shlnbf,shlmin,shlmax,f,alpha,beta,h,c,d,jmat,
     $                  kmat,values,nbf,nnp,nshell,ndmat,ncoul,
     $                  nexch,ntriang,t1,t2,lag,ops,nbfx,nnpx)
c
c***begin prologue     lagrng
c***date written       860819  (yymmdd)
c***revision date      861211  (yymmdd)
c
c   11 december 1986   pws at lanl
c      modifying to calculate virtual-virtual block of the lagrangian
c      but to only transform the occupied-occupied block to the ao basis
c      since the gradient terms only use that block.
c
c***keywords
c***author             saxe, paul (lanl)
c***source             %W%   %G%
c***purpose            formation of the general scf lagrangian
c***description
c
c***references
c***routines called    fock     (l601)
c***end prologue       lagrng
c
      implicit integer (a-z)
c
      character*(*) ops
      integer shlnbf(nshell),shlmin(nshell),shlmax(nshell)
      real*8 f(nshell),alpha(nshell,nshell)
      real*8 beta(nshell,nshell),h(nnp),c(nbf,nbf)
      real*8 d(nnp,ndmat),jmat(nnp,ncoul),kmat(nnp,ncoul)
      real*8 t1(nbf,nbf),t2(nbf,nbf),lag(nbf,nbf),values(*)
      real*8 c1,c2,haa,hbb,hab
      logical printj,printk
      logical logkey
c
      common /io/ inp,iout
c
      data printj,printk /.false.,.false./
      data c1,c2 /0.0d+00, 0.0d+00/
      save printj,printk,c1,c2
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
c     ----- form the density matrices for each orbital shell -----
c
      do 1 dmat=1,ndmat
         call gdmat(d(1,dmat),c,nbf,nnp,shlmin(dmat),shlmax(dmat))
    1 continue
c
      call iosys('write real "hf density matrix" to rwf',nnp*ndmat,
     $            d,0,' ')
c
c     ----- form the coulomb and exchange matrices -----
c
      xs=1
      x2=xs+nnp*nbfx
      t6=x2+nnpx 
      t7=t6+nbf*nbf*nbfx
      t8=t7+nbfx*nbfx
      t9=t8+nbf*nbf*nbfx
      t10=t9+nbf*nbf*nbfx
      t11=t10+nbf*nbf*ncoul
      top=t11+nbf*nbf*nexch
      call iosys('read real "expansion overlap integrals" from rwf',
     $           nnp*nbfx,values(xs),0,' ')
      call iosys('read real "auxiliary two-electron integrals"'
     $           //' from rwf',nnpx,values(x2),0,' ')
      call formjk(values(xs),values(x2),nbfx,nnpx,d,nnp,nbf,
     $            jmat,kmat,ncoul,nexch,ndmat,t2,values(t6),
     $            values(t7),values(t8),values(t9),values(t10),
     $            values(t11))
c      call fock(values,d,values,nnp,nbf,jmat,kmat,ncoul,nexch,ntriang,
c     $          -nexch,1,t1,h,ndmat,t2,printj,printk,c1,c2,haa,hbb,hab)
c
c     ----- transform to the mo basis and form the lagrangian as we go
c           l(ij)=f(i)h(ij) + sum(l=occ) {alpha(il)[ij;ll]+beta(il)[il;jl]}
c
      call rzero(lag,nbf**2)
c
      do 100 ishell=1,nshell-1
         mini=shlmin(ishell)
         nbfi=shlnbf(ishell)
c
c        ----- the one-electron term -----
c
         call trtosq(t2,h,nbf,nnp)
         do 5 i=mini,shlmax(ishell)
            do 4 j=1,nbf
               lag(i,j)=lag(i,j)+f(ishell)*t2(i,j)
    4       continue
    5    continue
c
         do 80 lshell=1,nshell-1
c
c           ----- transform j(ll) and k(ll) over i and j orbital ranges
c
            call trtosq(t1,jmat(1,lshell),nbf,nnp)
            call ebc(t2,t1,c,nbf,nbf,nbf)
            call ebtc(t1,c,t2,nbf,nbf,nbf)
c
c           ----- put j in correct place -----
c
            do 10 i=mini,shlmax(ishell)
               do 9 j=1,nbf
                  lag(i,j)=lag(i,j)+alpha(ishell,lshell)*t1(i,j)
    9          continue
   10       continue
c
c           ----- exchange (k) part -----
c
            call trtosq(t1,kmat(1,lshell),nbf,nnp)
            call ebc(t2,t1,c,nbf,nbf,nbf)
            call ebtc(t1,c,t2,nbf,nbf,nbf)
c
c           ----- put k in correct place -----
c
            do 20 i=mini,shlmax(ishell)
               do 19 j=1,nbf
                  lag(i,j)=lag(i,j)+beta(ishell,lshell)*t1(i,j)
   19          continue
   20       continue
   80    continue
  100 continue
c
      if (logkey(ops,'print=scf=lagrangian=mo',.false.,' ')) then
         write (iout,120)
  120    format('1',//,t20,'mo hartree-fock lagrangian',/)
         call matout(lag,nbf,nbf,nbf,nbf,iout)
      end if
c
      call iosys('write real "scf mo lagrangian" to rwf',nbf**2,lag,
     $     0,' ')
c
c     ----- transform to the ao basis -----
c
      call ebct(t1,lag,c,nbf,nbf,nbf)
      call ebc(t2,c,t1,nbf,nbf,nbf)
c..bhl
c..... symmetrize lagrangian
c..bhl
       call fixlag(nbf,t2,t1)
c..bhl
      call sqtotr(t1,t2,nbf,nnp)
c
      if (logkey(ops,'print=scf=lagrangian=ao',.false.,' ')) then
         write (iout,130)
  130    format ('1',//,t20,'ao hartree-fock lagrangian',/)
         call print(t1,nnp,nbf,iout)
      end if
c
      call iosys('write real "scf ao lagrangian" to rwf',nnp,t1,0,' ')
c
c
      return
      end
