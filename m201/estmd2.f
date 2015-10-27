*deck @(#)estmd2.f	5.1  11/6/94
      subroutine estmd2(nvar,nvv,nz,toang,x,fc,ic,ianz,iz,lbl,lalpha,
     $                lbeta,xangst,isave)
c***begin prologue     estmd2.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)estmd2.f	5.1   11/6/94
c***purpose            
c***description
c     make guesses at the diagonal second derivatives.  bending force
c     constants are 1.0.  for stretches the value chosen depends upon
c     the internuclear distance and the rows of the periodic table
c     in which the two affected atoms reside.  dummies and atoms higher
c     than cl are considered the same as hydrogens.
c
c***references
c
c***routines called
c
c***end prologue       estmd2.f
      implicit none
c     --- input variables -----
      integer nvar,nvv,nz
      real*8 toang
c     --- input arrays (unmodified) ---
      integer ianz(nz),iz(4,nz),lbl(nz),lalpha(nz),lbeta(nz)
      integer ic(nvar),isave(nvar)
      real*8 x(nvar),xangst(nvar)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 fc(nvv)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer irow(18)
      integer i,ivar,iatno,ia,jatom,jatno,ib,ii,nrep,nsave
      real*8 aa(3,3)
      real*8 hartre,constr,conbnd,aaa,bbb,one
c
      parameter (aaa=4.0d+00,bbb=1.0d+00,one=1.0d+00)
c
      data aa
     $/-.129d0,.186d0,.349d0,.186d0,.574d0,.805d0,.349d0,.805d0,1.094d0/
      data irow/2*1,8*2,8*3/
      save aa,irow
c
c     --- get hartree conversion.
      call iosys('read real j/hartree from rwf',1,hartre,0,' ')
      hartre=hartre*(1.0d+18)
c
c
      if(nz.ge.2) then
c        --- convert bond lengths to angstroms.
         do 15 i=1,nvar
            xangst(i)=x(i)
            if(nrep(i,lbl,nz).ne.0) xangst(i)=xangst(i)*toang
   15    continue
c        --- guess diagonal second derivatives in mdyne units.
         do 20 i=2,nz
            ivar=abs(lbl(i))
            if(ivar.ne.0) then
               iatno = ianz(i)
               ia    = irow(iatno)
               if(iatno.lt.1.or.iatno.gt.18) ia = 1
               jatom = iz(1,i)
               jatno = ianz(jatom)
               ib    = irow(jatno)
               if(jatno.lt.1.or.jatno.gt.18) ib = 1
               ii=(ivar*(ivar+1))/2
               if(ic(ivar).eq.0)
     $                 fc(ii)=fc(ii) +aaa/((xangst(ivar)-aa(ia,ib))**3)
            endif
   20    continue
         if(nz.ge.3) then
            do 30 i=3,nz
               ivar=abs(lalpha(i))
               if(ivar.ne.0) then
                  ii=(ivar*(ivar+1))/2
                  if(ic(ivar).eq.0) fc(ii)=fc(ii)+bbb
               endif
   30       continue
            if(nz.ge.4) then
               do 40 i=4,nz
                  ivar=abs(lbeta(i))
                  if(ivar.ne.0) then
                     ii=(ivar*(ivar+1))/2
                     if(ic(ivar).eq.0) fc(ii)=fc(ii)+bbb
                  endif
   40          continue
            endif
         endif
      endif
c
c     --- convert to atomic units.
      constr=toang**2/hartre
      conbnd=one/hartre
      call izero(isave,nvar)
      nsave = 0
      if(nz.ge.2) then
         do 60 i=2,nz
            ivar=abs(lbl(i))
            if(ivar.ne.0 .and. nrep(ivar,isave,nsave).eq.0) then
               nsave=nsave+1
               isave(nsave)=ivar
               ii=(ivar*(ivar+1))/2
               if(ic(ivar).eq.0) fc(ii)=fc(ii)*constr
            endif
   60    continue
c
         if(nz.ge.3) then
            call izero(isave,nsave)
            nsave=0
            do 80 i=3,nz
               ivar=abs(lalpha(i))
               if(ivar.ne.0 .and. nrep(ivar,isave,nsave).eq.0) then
                  nsave=nsave+1
                  isave(nsave)=ivar
                  ii=(ivar*(ivar+1))/2
                  if(ic(ivar).eq.0) fc(ii)=fc(ii)*conbnd
               endif
   80       continue
c
            if(nz.ge.4) then
               do 100 i=4,nz
                  ivar=abs(lbeta(i))
                  if(ivar.ne.0 .and. nrep(ivar,isave,nsave).eq.0) then
                     nsave=nsave+1
                     isave(nsave)=ivar
                     ii=(ivar*(ivar+1))/2
                     if(ic(ivar).eq.0) fc(ii)=fc(ii)*conbnd
                  endif
  100          continue
            endif
         endif
      endif
c
c
      return
      end
