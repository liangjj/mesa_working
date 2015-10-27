*deck %W%  %G%
      subroutine pm715(z,a)
c***begin prologue     pm715.f
c***date written       851113   (yymmdd)
c***revision date      11/6/94 
c    6 may     1994    rlm at lanl
c      modifying m712 to do just the coulomb derivatives for dft.
c***keywords           m715, link 711, derivative integrals
c***author             saxe, paul and martin, richard (lanl)
c***source             %W%   %G%
c***purpose            computes coulomb integral(J) derivatives
c                      general contraction scheme.
c***description
c
c***references
c
c***routines called
c***end prologue       pm715.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer a(*)
      real*8 z(*)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      character*4096 ops
      character*8 prtflg
      logical prnt,logkey
      integer inp,iout
      integer intkey
      integer nderiv
      integer nat,nbf,nnp
      integer grad,d2e,zan,c,nd1e,nd2e
      integer nprim,ncont,ntypes,nbtype,mxcont,maxl,lenxyz
      integer ptprim,noprim,ptcont,nocont,start,nocart,nobf
      integer minmom,maxmom,mintyp,nx,ny,nz,cont,ex,top
      integer npf,nnprim,pstart,prtoao
      integer ndmat,nshell,alpha,beta,dpr,need,rneed
      integer maxcor
      integer wpadti,iadtwp
      real*8 cutexp,fpkey
c
c 
      integer mynodeid, nodeid, mdtob, mitob, nproc, nnodes
      include 'msgtypesf.h'
      logical ispar
      common /tcgmesa/ ispar
c
      common /io/     inp,iout
c
      data prnt/.true./
      save prnt
c
 1000 format(1x,'m715: two-electron derivative integrals',
     $       /,5x,'gradient integral pre-exponential cutoff:',1pe8.1)
 1010 format (5x,'the coulomb contribution to the gradients')
 1020 format (5x,'the coulomb contribution to the force ',
     $       'constants')
c
c     
      mynodeid=nodeid()
      ispar=.true.
      if (mynodeid .eq. 0) then
c     
c     --- recover the options string ---
         call iosys('read character options from rwf',-1,0,0,ops)
         ispar=logkey(ops,'parallel',.false.,' ')
c     
c     --- set cutoffs for preexponential and final integral ---
         cutexp=fpkey(ops,'gradient=preexponential',1.0d-12,' ')
c
c     --- see if print has been turned off externally.
         call iosys('read character "print flag" from rwf',-1,0,0,
     $        prtflg)
         if(prtflg.eq.'minimum') prnt=.false.
c
c     --- what level of derivatives should we do.
         nderiv=intkey(ops,'nderiv',1,' ')
         if(logkey(ops,'force-constants',.false.,' ')) then
            if(logkey(ops,'force-constants=numerical',.false.,' '))
     $           then
c     numerical force constants, do nothing
            else
c     analytic force constants
               nderiv=2
            endif 
         endif
c
c     --- announce our presence.
         if(prnt) write (iout,1000) cutexp
      endif
      if (ispar) then
         call brdcst(1+MSGDBL,cutexp,mdtob(1),0)
         call brdcst(2+MSGINT,nderiv,mitob(1),0)
      endif
c
c     --- how much core have we
      call getscm(0,z,maxcor,' ',0)
      if (mynodeid.eq.0) then
c
c     --- get lengths, dimensions, etc. needed for core allocation --
         call iosys('read integer "number of atoms" from rwf',1,nat,0,
     $        ' ')
         call iosys(
     $        'read integer "number of basis functions" from rwf',
     $        1,nbf,0,' ')
      endif
      if (ispar) then
         call brdcst(3+MSGINT,nat,mitob(1),0)
         call brdcst(4+MSGINT,nbf,mitob(1),0)
      endif
      nnp=(nbf+1)*nbf/2
c
c     --- allocate core for gradients and hessians ----
      if (nderiv.eq.0) then
         grad=1
         d2e=1
         zan=1
      else if (nderiv.eq.1) then
         grad=1
         d2e=1
         zan=grad+3*nat
      else if (nderiv.eq.2) then
         nd1e=3*nat
         nd2e=nd1e*(nd1e+1)/2
         grad=1
         d2e=grad+nd1e
         zan=d2e+nd2e
      end if
      c=zan+nat
      top=wpadti(c+3*nat)
c
c     --- retrieve basis set information; returns pointers as well.
      call pbasis(nat,nbf,nprim,ncont,ntypes,nbtype,lenxyz,
     $           mxcont,maxl,ptprim,noprim,ptcont,nocont,start,
     $           nocart,nobf,minmom,maxmom,mintyp,nx,ny,nz,
     $           cont,ex,top,z,a)
c
c     --- get additional basis set information
      if (mynodeid.eq.0) then
         call iosys('read integer "number of primitive functions" '//
     $        'from rwf',1,npf,0,' ')
      endif
      if (ispar) call brdcst(5+MSGINT,npf,mitob(1),0)
      nnprim=(npf+1)*npf/2
c
      pstart=top
      prtoao=iadtwp(pstart+ntypes*nat)
      top=wpadti(prtoao+npf*nbf)
      if (mynodeid .eq. 0) then
         call iosys(
     $        'read integer "pointer to first primitive" from rwf',
     $        -1,a(pstart),0,' ')
         call iosys('read real t(prim,cont) from rwf',
     $        npf*nbf,z(prtoao),0,' ')
      endif
      if (ispar) then
         call brdcst(6+MSGINT,a(pstart),mitob(ntypes*nat),0)
         call brdcst(7+MSGDBL,z(prtoao),mdtob(npf*nbf),0)
      endif
c
c     --- get room for scf descriptors.
      if (mynodeid .eq. 0) then
         call iosys('read integer "number of hf density matrices" '//
     $        'from rwf',1,ndmat,0,' ')
         call iosys('read integer "number of shells" from rwf',
     $        1,nshell,0,' ')
      endif
      if (ispar) then
         call brdcst(8+MSGINT,ndmat,mitob(1),0)
         call brdcst(9+MSGINT,nshell,mitob(1),0)
      endif
      alpha=iadtwp(top)
      beta=alpha+nshell**2
      top=wpadti(beta+nshell**2)
c
c     --- room for the primitive density matrix.
      dpr=top
c
      need=wpadti(dpr+nnprim*ndmat)
      rneed=iadtwp(need)
c
      if (need.gt.maxcor) then
         call plnkerr('need more core for m715',100)
      end if
c
c     --- read in basis set information from the read-write file ---
c
      if (mynodeid .eq. 0) then
         call iosys('read real coordinates from rwf',-1,z(c),0,' ')
         call iosys('read real alpha from rwf',
     $        nshell**2,z(alpha),0,' ')
         call iosys('read real beta from rwf',
     $        nshell**2,z(beta),0,' ')
      endif
      if (ispar) then
         call brdcst(10+MSGDBL,z(c),mdtob(3*nat),0)
         call brdcst(11+MSGDBL,z(alpha),mdtob(nshell**2),0)
         call brdcst(12+MSGDBL,z(beta),mdtob(nshell**2),0)
      endif
c
c     --- calculate two electron coulomb integrals ---
      call driver(a(ptprim),a(noprim),nbtype,z(ex),z(c),
     $     a(nx),a(ny),a(nz),lenxyz,a(nobf),a(nocart),
     $     a(mintyp),a(minmom),a(maxmom),
     $     z(rneed),iadtwp(maxcor)-need+1,nat,npf,nnprim,nprim,
     $     ops,cutexp,a(need),z(grad),nderiv,
     $     prnt,ndmat,z(alpha),z(beta),nshell,a(pstart),z(dpr),
     $     z(d2e),nd2e,
     $     nbf,
     $     z(cont),ncont,a(ptcont),a(nocont),a(start))
      if (mynodeid .ne. 0) goto 999
c
c     --- print gradient contribution from coulomb integrals
      if (logkey(ops,'print=gradient=coulomb',.false.,' ')) then
         write (iout,1010)
         call matout(z(grad),3,nat,3,nat,iout)
         if (nderiv.ge.2) then
            write (iout,1020)
            call print(z(d2e),nd2e,nd1e,iout)
         end if
      end if
c
c     --- store them on the rwf.
      call iosys('write real "coulomb integral derivatives"'
     $           //' to rwf',3*nat,z(grad),0,' ')
      if (nderiv.ge.2) then
         call iosys('write real "coulomb integral force constants"'
     $              //' to rwf', nd2e,z(d2e),0,' ')
      end if
c
c     --- and exit gracefully ---
      call chainx(0)
c all nodes branch here except node 0
 999  continue
c
c
      return
      end