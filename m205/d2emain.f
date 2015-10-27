*deck @(#)d2emain.f	5.1  11/6/94
      subroutine d2emain(nvar,nvv,nz,natoms,maxpt,toang,ops,x,f,xx,ff,
     $                frcnst,fs,vname,lbl,ian,atmass,d2edone,abnrml,
     $                z,zsq,cmass,atmchg,c,ianz,iz,bl,alpha,beta,
     $                klbl,lalpha,lbeta,scr1,scr2,scr3,scr4,scr5,scr6,
     $                scr7,scr8,scr9,scr10,iscr11,xc,fc,xxc,ffc,zf,
     $                frcnsc,zsqc,nvvc)
c***begin prologue     d2emain.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)d2emain.f	5.1   11/6/94
c***purpose            
c***description
c
c     --- local arrays of note.
c     x      ... coordinate vector; real(nvar).
c     f      ... forces(-de/dx); real(nvar).
c     frcnst ... second derivative matrix if evaluated externally.
c                real(nvv).
c     vname  ... variable names vector; character*16(nvar).
c     fs     ... value of the function at earlier points; real(maxpt).
c     xx     ... previous coordinates; real(nvar,maxpt)
c     ff     ... previous forces; real(nvar,maxpt)
c     xxc    ... previous cartesian coordinates (if cartesian)
c     ffc    ... previous forces                  "    "
c     ian    ... atomic numbers
c     atmass ... atomic masses
c
c      common/d2einf/energy,rmax,rmin,rlim,d2ecycl,
c     $               prnt,chkpt,singpt,stpsize,cartfx
c     energy ... the value of the function at the current point.
c     rmax   ... unused
c     rmin   ... unused
c     rlim   ... unused
c     d2ecycl ... current cycle number.
c     prnt   ... print switch.
c     updrwf ... update the rwf?
c     chkpt  ... checkpoint the run?
c     freq   ... calculate vibrational frequencies ?
c     singpt ... single point finite difference formula if true
c                (default)
c     stpsize... step size for finite difference calculation
c
c***references
c
c***routines called
c
c***end prologue       d2emain.f
      implicit none
c     --- input variables -----
      integer nvar,nvv,nz,natoms,maxpt,nvvc
      logical d2edone,abnrml
      real*8 toang
c     --- input arrays (unmodified) ---
      character*(*) ops
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
      integer lbl(nz),lalpha(nz),lbeta(nz),klbl(nz)
      integer ianz(nz),iz(4,nz),ian(nz)
      integer iscr11(3*natoms)
      character*(*) vname(nvar)
      real*8 x(nvar),f(nvar),frcnst(nvv),atmass(natoms)
      real*8 frcnsc(nvvc),zsqc(3*natoms,3*natoms)
      real*8 cmass(3*natoms)
      real*8 xx(nvar,maxpt),ff(nvar,maxpt),fs(maxpt)
      real*8 xc(3*natoms),fc(3*natoms)
      real*8 xxc(natoms*3,maxpt),ffc(natoms*3,maxpt)
      real*8 z(9*nvar),zsq(2*nvar*nvar)
      real*8 bl(nz),alpha(nz),beta(nz),atmchg(nz)
      real*8 c(3*natoms)
      real*8 scr1(nz),scr2(nz),scr3(nz),scr4(nz),scr5(nz),scr6(nz)
      real*8 scr7(3*natoms*nvar),scr8(3*natoms),scr9(3*natoms)
      real*8 scr10(3*natoms)
      real*8 zf(*)
c     --- local variables ---
      integer d2ecycl
      integer inp,iout
      integer nat3,nnsq,nnp,ndegf,numtet
      logical exit,freq
      logical prnt,chkpt,singpt,cartfx
      logical logkey
      character*4 ians
      real*8 energy,rmax,rmin,rlim,stpsize
c
      common/io/inp,iout
      common/d2einf/energy,rmax,rmin,rlim,d2ecycl,
     $               prnt,chkpt,singpt,stpsize,cartfx
c
 1000 format(1x,'m205:numerical second derivatives')
 1010 format(5x,'initialization pass')
 1020 format(5x,'restart from: ',a8)
 1030 format(5x,'d2e cycle',i5)
 2000 format(5x,'calculation of vibrational frequencies'/)
 3000 format(/,5x,'warning: not all internal motions are',
     $      /5x,'being considered in the frequency calculation',
     $      /5x,'there are ',i2,' degrees of freedom',
     $      /5x,' and only ',i2,' variables ')
c
c     --- announce our presence
      write(iout,1000)
      exit=.false.
      call iosys('read integer zlbl from rwf',nz,lbl,0,' ')
c
c     --- determine the type of entry.  'd2e_cycle' has not
c         been written to the rwf is this is a restart or the first 
c         point of a new optimization.
      call iosys('does d2e_cycle exist on rwf',0,0,0,ians)
      if (ians.eq.'no') then
c        --- initial entry
         write(iout,1010)
         call d2einit(nvar,nvv,maxpt,ops,
     $                x,vname,toang,xc,natoms)
         call d2einfo('write',nvar,nvv,maxpt,vname,x,f,frcnst,xx,
     $                 ff,xc,fc,xxc,ffc,fs,natoms)
         d2edone=.false.
         return
      else
c        --- a continuation of the current run.
         call iosys('read integer d2e_cycle from rwf',1,d2ecycl,
     $              0,' ')
         write(iout,1030) d2ecycl
         call iosys('write character "print flag" to rwf',0,0,0,
     $              'minimum')
         call d2einfo('read',nvar,nvv,maxpt,vname,x,f,frcnst,xx,
     $                ff,xc,fc,xxc,ffc,fs,natoms)
      endif
c
      if(cartfx) then
         call d2ecstp(maxpt,xc,fc,xxc,ffc,d2edone,z,
     $                natoms)
      else
         call d2estep(nz,nvar,nvv,maxpt,toang,vname,x,f,xx,
     $                ff,d2edone,z)
      endif
c
      abnrml=exit.and.(.not.d2edone)
      if(d2edone) then
c        --- form the force constant matrix from the gradients
c            which are stored in ff
         if(cartfx) then
            nat3=natoms*3
            nnsq=nat3*nat3
            nvvc=nat3*(nat3+1)/2
            call d2ecfrm(nat3,nnsq,maxpt,xxc,ffc,nvvc,frcnsc,zsqc)
         else
            call d2eform(nvar,nvv,maxpt,frcnst,xx,ff,zsq,vname)
         endif
c
c        --- write the original z values and coordinates to rwf
         call iosys('write real zvalues to rwf',nvar,x,0,' ')
         call iosys('write real coordinates to rwf',nat3,xc,0,' ')
c
         freq=logkey(ops,'frequencies',.false.,' ')
         if(freq) then
            if(cartfx) then
               call iosys('read integer "atomic numbers" from rwf',
     $                    -1,ian,0,' ')
               call filmas(0,iout,ian,natoms,.false.,atmass,' ')
               nat3=natoms*3
               nnp=nat3*(nat3+1)/2
               call vibfrc(zsqc,natoms,ian,atmass,cmass,
     $                     nvar,xc,nat3,nnp,zf,vname)
            else
               write(iout,2000)
c              --- check to see if a full set of internal coordinates 
c                  is used in the difference procedure
               ndegf=3*natoms-6
               if(nvar.lt.ndegf) then
                  write(iout,3000) ndegf,nvar
                  call lnkerr(' hope you saved the check file ')
               endif
c              --- form the internal-to-cartesian transformation matrix
               call zctran(nz,nvar,ianz,iz,bl,
     $                 alpha,beta,.false.,numtet,natoms,
     $                 klbl,lalpha,lbeta,ian,atmchg,c,scr6,
     $                 scr1,scr2,scr3,scr4,scr5,scr7,xx,scr8,
     $                 scr9,scr10,iscr11,zf)
               call iosys('read integer "atomic numbers" from rwf',
     $                    -1,ian,0,' ')
               call filmas(0,iout,ian,natoms,.false.,atmass,' ')
               nat3=natoms*3
               nnp=nat3*(nat3+1)/2
               call vibfreq(frcnst,natoms,ian,atmass,cmass,xx,ff,
     $                      nvv,nvar,maxpt,scr7,scr8,nat3,nnp,
     $                      iscr11,zf,vname)
            endif
         endif
c
c        --- store the information for this step.
         call d2einfo('write',nvar,nvv,maxpt,vname,x,f,frcnst,xx,
     $               ff,xc,fc,xxc,ffc,fs,natoms)
      else
c        --- update the rwf and exit the link.
         if(cartfx) then
            call iosys('write real coordinates to rwf',natoms*3,xc,
     $                 0,' ')
         else
            call iosys('write real zvalues to rwf',nvar,x,0,' ')
         endif
         call d2einfo('write',nvar,nvv,maxpt,vname,x,f,frcnst,xx,
     $               ff,xc,fc,xxc,ffc,fs,natoms)
      endif
c
c
      return
      end
