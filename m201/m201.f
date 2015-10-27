*deck @(#)m201.f	5.1  11/6/94
      program m201
c***begin prologue     pm201.f
c***date written       851001  
c***revision date      11/6/94      
c   may 30, 1991       rlm at lanl
c     fixing optimization status iosys write during fp optimization
c
c***keywords           
c   m201, link 201, geometry optimization, berny, fletcher-powell
c***author             martin, richard(lanl)
c***source             @(#)pm201.f	5.1   11/6/94
c***purpose            
c***description
c   contains fletcher-powell algorithm : if analytical gradients are unavailable
c            berny                     : if analytical gradients available
c            murtaugh-sargeant         : if analytical gradients available
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       pm201.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer a
      real*8 z
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer idum, ngot
      integer maxcor,mxpt
      integer nvar,nz
      integer lbl,pool0,top,pool1,delvar,yold,d1var,d2var,d1vold
      integer xi,h,lalpha,lbeta,intvec,fpvec,scr,x,frcnst,fc,nvv
      integer ic,fs,iadtwp,wpadti
      integer f,xx,ff,scrsq,hinv,iscr
      character*4096 ops
      character*4096 ref
      character*4 ians
      character*16 vname(500)
      logical convrg,abnrml,logkey
      real*8 toang
c
      data mxpt/49/, convrg/.false./, abnrml/.false./
c
      common/io/inp,iout
      pointer (p,z(1)),(p,a(1))
c
      call drum
      ref='opt=(murtagh-sargent,berny,fletcher-powell,force,'//
     $     'force-constants)'
c
c     --- get the link options.
      call iosys('read character options from rwf',-1,0,0,ops)
      call iosys('read integer "number of variable coordinates" '//
     $           'from rwf',1,nvar,0,' ')
      call iosys('read integer "number of z-matrix entries" from rwf',
     $           1,nz,0,' ')
      call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
      call iosys('does "optimization status" exist on rwf',0,0,0,ians)
      if (ians.eq.'no') then
         convrg=.false.
      else
         call iosys('read integer "optimization status" from rwf',
     $              1,convrg,0,' ')
      end if
      if(logkey(ops,'opt=fletcher-powell',.false.,ref)) then
c        --- fletcher-powell optimization.
c            should we kick out of the optimization loop?
         if(convrg) then
c           --- we converged the previous cycle. execute the route jump.
         else
c           --- allocate core.
            lbl=1
            lalpha=lbl+nz
            lbeta=lalpha+nz
            intvec=lbeta+nz
            pool0=iadtwp(intvec+nvar)
            pool1=pool0+nvar
            delvar=pool1+nvar
            yold=delvar+nvar
            d1var=yold+nvar
            d2var=d1var+nvar
            d1vold=d2var+nvar
            xi=d1vold+nvar
            h=xi+nvar
            fpvec=h+nvar*nvar
            scr=fpvec+nvar
            top=wpadti(scr+2*(nvar+nvar*nvar))
            call getmem(top,p,ngot,'m201',0)
c
c
            call fpmain(ops,vname,nvar,nz,convrg,abnrml,toang,z(pool0),
     $                  z(pool1),z(delvar),z(yold),z(d1var),z(d2var),
     $                  z(d1vold),z(xi),z(h),a(lbl),a(lalpha),a(lbeta),
     $                  a(intvec),z(fpvec),z(scr))
c
c           --- even if we've converged, we want to do one last iteration to
c               restore the wavefunction at the minimum, etc.  in order to
c               do this, we write the 'converged' status to the rwf where
c               it will be picked up the next time through here. 
            call iosys('write integer "optimization status" to rwf',
     $                 1,convrg,0,' ')
c
c           --- in order to beat the convergence test at the end of this
c               block, set convrg to false.
            convrg=.false.
            call getmem(-ngot,p,idum,'m201',idum)
         endif
      else if(logkey(ops,'opt=berny',.false.,ref)) then
c           --- berny optimization.
c               allocate core.
         nvv=nvar*(nvar+1)/2
         mxpt=49
         x=1
         f=x+nvar
         xx=f+nvar
         ff=xx+nvar*mxpt
         fc=ff+nvar*mxpt
         frcnst=fc+2*nvv
         fs=frcnst+nvv
         ic=wpadti(fs+mxpt)
         lbl=ic+mxpt
         scr=iadtwp(lbl+nz)
         scrsq=scr+9*nvar
         iscr=wpadti(scrsq+2*nvar*nvar)
         top=iscr+7*nz+2*nvar
         call getmem(top,p,ngot,'m201',0)
c
         call bsmain(nvar,nvv,nz,mxpt,toang,ops,z(x),z(f),z(xx),
     $        z(ff),z(fc),z(frcnst),z(fs),vname,
     $        a(ic),a(lbl),convrg,abnrml,z(scr),z(scrsq),
     $        a(iscr))
         call iosys('write integer "optimization status" to rwf',
     $              1,convrg,0,' ')
         call getmem(-ngot,p,idum,'m201',idum)
      else
c        --- murtaugh-sargent optimization.
c            allocate core.
         nvv=nvar*(nvar+1)/2
         mxpt=2
         x=1
         f=x+nvar
         xx=f+nvar
         ff=xx+nvar*mxpt
         fc=ff+nvar*mxpt
         frcnst=fc+2*nvv
         fs=frcnst+nvv
         ic=wpadti(fs+mxpt)
         hinv=iadtwp(ic+mxpt)
         lbl=wpadti(hinv+nvar*nvar)
         scr=iadtwp(lbl+nz)
         scrsq=scr+6*nvar
         iscr=wpadti(scrsq+nvar*nvar)
         top=iscr+7*nz+2*nvar
         call getscm(top,p,ngot,'m201',0)
c
c
         call msmain(nvar,nvv,nz,mxpt,toang,ops,z(x),z(f),z(xx),
     $        z(ff),z(fc),z(frcnst),z(fs),vname,
     $        a(ic),z(hinv),a(lbl),convrg,abnrml,z(scr),
     $        z(scrsq),a(iscr))
         call iosys('write integer "optimization status" to rwf',
     $               1,convrg,0,' ')
         call getscm(-ngot,p,idum,'m201',idum)
      endif
c
c     --- if we're not converged, override the jump in the route.
c         if we have converged, or have an exit flag
c         for exceeding the maximum number of iterations,etc.,
c         execute the jump.
      if(convrg.or.abnrml) then
         call chainx(0)
      else
         call chainx(1)
      endif
c
c
      stop
      end
