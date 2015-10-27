*deck @(#)lagrng.f	5.1 11/6/94
      subroutine lagrng(shlmin,shlmax,f,alpha,beta,h,c,d,jmat,
     $                  kmat,values,nbf,nnp,nshell,ndmat,ncoul,
     $                  nexch,ntriang,t1,t2,lag,slater,becke,vwn,
     $                  lyp,natoms,ngrid,mxgrd,grdwts,calc,dmcut,
     $                  dencut,nzptrs,nzptrs2,tmpgwt,tmpgwt2,ops)
c
c***begin prologue     lagrng.f
c***date written       860819  
c***revision date      11/6/94      
c
c   17 june,     1993  rlm at lanl
c      modifying to calculate kohn-sham lagrangian.
c   11 december, 1986  pws at lanl
c      modifying to calculate virtual-virtual block of the lagrangian
c      but to only transform the occupied-occupied block to the ao basis
c      since the gradient terms only use that block.
c
c***keywords           
c***author             saxe, paul (lanl) 
c***source             @(#)lagrng.f	5.1   11/6/94
c***purpose            formation of the general scf lagrangian  
c***description       
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       lagrng.f
      implicit none
c     --- input variables -----
      integer nbf,nnp,nshell,ndmat,ncoul,nexch,ntriang
      integer natoms,mxgrd
      real*8 dencut,dmcut
      logical slater,becke,vwn,lyp
      character*(*) calc
c     --- input arrays (unmodified) ---
      integer shlmin(nshell),shlmax(nshell)
      integer ngrid(natoms)
      character*(*) ops
      real*8 f(nshell),alpha(nshell,nshell),beta(nshell,nshell)
      real*8 h(nnp),c(nbf,nbf),d(nnp,ndmat)
      real*8 grdwts(mxgrd,natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 lag(nbf,nbf)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 jmat(nnp,ncoul),kmat(nnp,ncoul)
      real*8 t1(nbf,nbf),t2(nbf,nbf),values(nnp,ntriang)
      integer nzptrs(mxgrd),nzptrs2(mxgrd)
      real*8 tmpgwt(mxgrd),tmpgwt2(mxgrd)
c     --- local variables ---
c
      character*8 calce
      logical printj,printk
      logical logkey
      integer inp,iout,ishell,lshell,i,j,dmat
      real*8 exc
c
      common /io/ inp,iout
c
      data printj,printk /.false.,.false./
      save printj,printk
c
c     --- get the scf vector ---
      call iosys('read real "scf vector" from rwf',nbf**2,c,0,' ')
c
c     --- form the density matrices for each orbital shell ---
      do 10 dmat=1,ndmat
         call gdmat(d(1,dmat),c,nbf,nnp,shlmin(dmat),shlmax(dmat))
   10 continue
      call iosys('write real "hf density matrix" to rwf',nnp*ndmat,
     $            d,0,' ')
c
c     --- form the coulomb and exchange matrices ---
      call jmatrix(values,d,nbf,nnp,jmat,ncoul,ntriang,ndmat)
      call kmatrix(values,d,nbf,nnp,kmat,nexch,grdwts,ndmat,natoms,
     $             mxgrd,ngrid,exc,slater,becke,lyp,vwn,calce,calc,
     $             dmcut,dencut,nzptrs,nzptrs2,tmpgwt,tmpgwt2)
c
c     --- transform to the mo basis ---
      call trtosq(lag,h,nbf,nnp)
      call tomo(nbf,nnp,nshell,ncoul,nexch,lag,h,jmat,kmat,
     $          c,t1,t2)
      call rzero(lag,nbf**2)
c
c     --- form the lagrangian ---
c         l(ij)=f(i)h(ij) + sum(l=occ) {alpha(il)[ij;ll]+beta(il)Vxc[ij;ll]}
      do 200 ishell=1,nshell-1
c
c        --- the one-electron term ---
         call trtosq(t2,h,nbf,nnp)
         do 50 i=shlmin(ishell),shlmax(ishell)
            do 40 j=1,nbf
               lag(i,j)=lag(i,j)+f(ishell)*t2(i,j)
   40       continue
   50    continue
c
         do 100 lshell=1,nshell-1
c
c           --- put j in correct place ---
            call trtosq(t1,jmat(1,lshell),nbf,nnp)
            do 70 i=shlmin(ishell),shlmax(ishell)
               do 60 j=1,nbf
                  lag(i,j)=lag(i,j)+alpha(ishell,lshell)*t1(i,j)
   60          continue
   70       continue
c
c           --- put k in correct place ---
            call trtosq(t1,kmat(1,lshell),nbf,nnp)
            do 90 i=shlmin(ishell),shlmax(ishell)
               do 80 j=1,nbf
                  lag(i,j)=lag(i,j)+beta(ishell,lshell)*t1(i,j)
   80          continue
   90       continue
  100    continue
  200 continue
c
      if (logkey(ops,'print=scf=lagrangian=mo',.false.,' ')) then
         write (iout,120)
  120    format('1',//,t20,'mo hartree-fock lagrangian',/)
         call matout(lag,nbf,nbf,nbf,nbf,iout)
      end if
c
      call iosys('write real "scf mo lagrangian" to rwf',nbf**2,lag,
     $            0,' ')
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
