*deck @(#)ftogrid0.f	5.1  11/28/95
      subroutine ftogrid0(ylm,ptlm,nang,vlmax,nrad,rpts, vlm,nlm,v)
c***begin prologue	ftogrid0.f
c***date written	951026	(yymmdd)
c***revision date	11/28/95
c
c***keywords		potential, spherical, harmonics
c***author		russo, thomas (lanl)
c***source		@(#)ftogrid0.f	5.1 11/28/95
c***purpose		forms the potential contributions from the spherical
c                       harmonic expansion around a single atom
c***description         vlm contains the radial solutions to the poisson
c                       equation as solved by rsolver.  Each spherical harmonic
c                       channel has a separate radial function.  The total
c                       potential at a point is therefore:
c
c v(r) = sum(atom) sum(l=1 to maxl,m=-l to l) Vlm(r,l,m,atom)*Ylm(omega)
c
c but we're only going to form one atom's contribution on its own grid, which
c will be integrated on that grid only.
c
c***references		we doan need no steenkin references
c***routines called
c***end prologue	ftogrid0.f
      implicit none
c     --- input variables (unchanged) ---
      integer nang, vlmax, nrad, nlm
c     --- input arrays (unchanged) ---
      real*8 ylm(nang,nlm),rpts(nrad),vlm(nrad,nlm)
      integer ptlm(0:vlmax,-vlmax:vlmax)
c     --- output variables ---
c     --- output arrays ---
      real*8 v(nang*nrad)
c     --- scratch space ---
      integer l,m,ir,omega,plm,ioff,lm,ierr,i
      integer stderr
      real*8 r,term
      real*8 pi,one,y00,four
      data one/1.0d0/
      data four/4.0d0/
c
      ierr=stderr()
      call rzero(v,nang*nrad)
c
c for the time being, we're gonna do this stupidly, with no attempt at
c cutoffs---we're already saving tons of time, we can save more when we're
c sure this is the Way To Go
c
c I'm also going to write out loops and not use vmul/vwxy/etc. because 
c I'll try to do it with no scratch space at all.  Then again, this may
c be ridiculously inefficient.  We'll see.
c
c for each lm channel, we need to multiply the vlm at each radius by
c all the angular values of ylm to form that channel's contribution at that
c radius, then stick it into v.
c
c first we do the l=0 term, which just needs y00=sqrt(1/(4pi))
c also remember that vlm is really ulm, which is r*vlm, so to get vlm need 
c to divide by r!
      pi=four*atan(one)
      y00=one/sqrt(four*pi)
      ioff=0
      do 10 ir=1,nrad
         r=rpts(ir)
         if (r.eq.0 .and. vlm(ir,1).ne.0) then
            write(ierr,*) 'division by zero imminent!'
            call plnkerr('div0 ftogrid0',1934)
         endif
         if (r.eq.0 .and. vlm(ir,1).eq.0) then
            go to 21
         else
            term=vlm(ir,1)*y00/r
         endif
         do 20 omega=1,nang
            v(ioff+omega)=v(ioff+omega)+term
 20      continue 
 21      continue 
         ioff=ioff+nang
 10   continue 
c      write(ierr,*)' after y00 term, v is:'
c      write(ierr,*) (v(i),i=1,nrad*nang)
c
c     now for the real spherical harmonicas
c
      lm=2
      do 30 l=1,vlmax
         do 40 m=-l,l
            plm=ptlm(l,m)
            ioff=0
            do 50 ir=1,nrad
               r=rpts(ir)
               if (r.eq.0.and.vlm(ir,lm).ne.0) then
                  write(ierr,*) 'division by zero imminent!'
                  call plnkerr('div0 ftogrid0',1935)
               endif
               if (r.eq.0) then
c                 vlm is zero on this radial shell, so skip it.
                  goto 61
               else
                  term=vlm(ir,lm)/r
               endif
               do 60 omega=1,nang
                  v(ioff+omega)=v(ioff+omega)+term*ylm(omega,plm)
 60            continue 
 61            continue 
               ioff=ioff+nang
 50         continue 
            lm=lm+1
 40      continue 
 30   continue 
c$$$      write(ierr,*)' v is:'
c$$$      write(ierr,*) (v(i),i=1,nrad*nang)

      if (lm.ne.nlm+1) then
         write(ierr,*)' error in nlm: ',lm-1,nlm
         call plnkerr('error in ftogrid0',1963)
      endif
      return
      end
