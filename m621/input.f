*deck @(#)input.f	5.1 11/6/94
      subroutine input(c,ian,nat,x,y,z,q,qd,qq,xcnt,ycnt,zcnt,rad,
     $                 xobs,yobs,zobs,
     $                 epsilon,npoints,nsrcs,ncnt,nobs,prnt,
     $                 nfocus,toang,ops)
c***begin prologue     input.f
c***date written       yymmdd
c***revision date      11/6/94
c   january 25, 1994   rlm at lanl
c      converting to atomic units throughout
c***keywords
c***author             tawa, greg(lanl)
c***source             @(#)input.f	5.1   11/6/94
c***purpose
c***description
c
c
c
c***references
c
c***routines called
c
c***end prologue       input.f
      implicit none
c     --- input variables -----
      integer nat,npoints,nsrcs,ncnt,nobs,nfocus
      logical prnt
      real*8 epsilon,toang
c     --- input arrays (unmodified) ---
      character*(*) ops
      integer ian(nat)
      real*8 c(3,nat)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 x(nsrcs), y(nsrcs), z(nsrcs)
      real*8 q(nsrcs),qd(3*nsrcs), qq(6*nsrcs)
      real*8 xcnt(ncnt), ycnt(ncnt), zcnt(ncnt), rad(ncnt)
      real*8 xobs(nobs), yobs(nobs), zobs(nobs)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
      integer mxatm
      integer inp,iout
      parameter (mxatm=103)
      logical logkey
c
c     --- define centers and radii for the bubble ---
c         the atomic coordinates enter in c(3,nat),
c         atomic numbers in ian(nat)
      do 10 i=1,nat
         xcnt(i)=c(1,i)
         ycnt(i)=c(2,i)
         zcnt(i)=c(3,i)
   10 continue
      call iosys('read real "solvent radii" from rwf',
     $           -1,rad,0,' ')
c
c     --- for now assume the sources and charges are atomic
c         centers and default charges.
      if(nsrcs.ne.nat) then
         call lnkerr('do not know how to do more sources than atoms')
      else
         if(logkey(ops,'solvent=mulliken',.false.,' ')) then
            call iosys('read real "scf mulliken charges" from rwf',
     $                  nat,q,0,' ')
c           --- must remove nuclear contribution ---
            do 20 i=1,nat
               q(i)=-(q(i)-float(ian(i)))
   20       continue
         else 
c           default to cox-williams charges.
            call iosys('read real "cox-williams monopoles" from rwf',
     $                  nat,q,0,' ')
            call iosys('read real "cox-williams dipoles" from rwf',
     $                  3*nat,qd,0,' ')
            call iosys('read real "cox-williams quadrupoles" from rwf',
     $                  6*nat,qq,0,' ')
         endif
         call vmove(x,xcnt,nsrcs)
         call vmove(y,ycnt,nsrcs)
         call vmove(z,zcnt,nsrcs)
      endif
c
c     --- observation points; to calculate electrostatic
c         contribution to the free energy of solvation, the
c         coordinates of the observations points must be identical
c         to the coordinates of the source points ---
      if(nobs.ne.nsrcs) then
         call lnkerr('number of observation points .ne. number sources')
      else
         call vmove(xobs,x,nobs)
         call vmove(yobs,y,nobs)
         call vmove(zobs,z,nobs)
      endif

c
c
      return
      end
