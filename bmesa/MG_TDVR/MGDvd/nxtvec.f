*deck nxtvec.f
c***begin prologue     nxtvec
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           iterative eigenvalue
c***author             schneider, barry (nsf)
c***source             
c***purpose            next davidson vector.
c***                   
c***description        the next vector in the sequence is computed
c***                   by taking the residual, dropping to the lowest
c***                   grid, solving exactly and then working back to
c***                   the top grid.
c***
c***                   
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       nxtvec
      subroutine nxtvec(vecout,resid,diag,basis,pm,array,ham,voft,pvt,
     1                  h0,h0ind,scale,cnverg,thresh,eps,nmax,dim,
     2                  spdim,msze,ntot,n3d,mgr,ntrial,maxit,maxvec,
     3                  mcycle,nx,ny,nz,nt,m,n,hx,hy,hz,ht,v,
     4                  grdi,prnt,mprnt)
      implicit integer (a-z)
      real*8 vecout, resid, diag, basis, array, voft, h0, soln, rhs
      real*8 scale, cnverg, thresh, eps, error, maxerr
      logical prnt, mprnt
      character*3 itoc
      character*80 title
      dimension nmax(4,mgr), msze(mgr), ntot(mgr), n3d(mgr)
      dimension vecout(*), resid(*), diag(*), basis(*)
      dimension array(*), voft(*)
      dimension pm(mgr,mgr,*), ham(mgr,*), pvt(mgr)
      dimension h0(*), h0ind(*)
      dimension maxit(mgr), maxvec(mgr), prnt(*), mprnt(*), ngot(2)
      common/io/inp, iout
      pointer(prhs,rhs(1))
      pointer(psoln,soln(1))
c
c
c     get space for the largest right hand side and two scratch vectors
c     to hold the solutions at previous and current grids.
c     
      need=wptoin(2*msze(grdi))
      call memory(need,prhs,ngot(1),'rhs',0)   
      call memory(need,psoln,ngot(2),'soln',0)   
c     
      do 10 cycle=1,mcycle
         write(iout,1) cycle
c
c        interpolate current value of residual on grdi to coarsest 
c        grid and store in rhs.
c
         if(mprnt(1)) then
            title='current grid residual on grid = '//itoc(grdi)
            call prntrm(title,resid,msze(grdi),1,msze(grdi),1,iout)
         endif 
         call gi2gj(resid,rhs,
     1              basis(pm(1,grdi,1)),basis(pm(1,grdi,2)),
     2              basis(pm(1,grdi,3)),basis(pm(1,grdi,4)),
     3              nmax(1,1),nmax(2,1),nmax(3,1),nmax(4,1),
     4              nmax(1,grdi),nmax(2,grdi),nmax(3,grdi),nmax(4,grdi),
     5              ntot(1),ntot(grdi),2,dim)       
         call vivo(rhs,rhs,
     1             basis(pm(1,1,1)),basis(pm(1,1,2)),
     2             basis(pm(1,1,3)),basis(pm(1,1,4)),
     3             nmax(1,1),nmax(2,1),nmax(3,1),nmax(4,1),ntot(1),
     4             2,dim)
         if(mprnt(2)) then
            title='residual interpolated to coarest grid'
            call prntrm(title,rhs,msze(1),1,msze(1),1,iout)
         endif 
c
c        solve exactly on coarest grid. put in soln
c
         call hamslv(h0,rhs,h0ind,msze(1),.false.)
         if(mprnt(1)) then
            title='coarest grid solution to residual equation'        
            call prntrm(title,rhs,msze(1),1,msze(1),1,iout)
         endif
c
c        we are now ready to go upward to the vector we want 
c        to approximate.
c
         if(grdi.gt.2) then
            lsol1 = 1
            lsol2 = lsol1 + msze(grdi)
            call copy(rhs,soln(lsol1),msze(1))
         else
            lsol2 = 1
            lsol1 = lsol2 + msze(grdi)
            call copy(rhs,soln(lsol2),msze(1))
         endif
         do 20 grdj=2,grdi-1
c
c           interpolate solution from from grdj-1 to grdj to provide
c           the initial guess 
c
            call gi2gj(soln(lsol1),soln(lsol2),
     1                 basis(pm(grdj,grdj-1,1)),
     2                 basis(pm(grdj,grdj-1,2)),
     3                 basis(pm(grdj,grdj-1,3)),
     4                 basis(pm(grdj,grdj-1,4)),
     5                 nmax(1,grdj),nmax(2,grdj),
     6                 nmax(3,grdj),nmax(4,grdj),
     7                 nmax(1,grdj-1),nmax(2,grdj-1),
     8                 nmax(3,grdj-1),nmax(4,grdj-1),
     9                 ntot(grdj),ntot(grdj-1),2,dim)
            call vivo(soln(lsol2),soln(lsol2),
     1                basis(pm(grdj,grdj,1)),
     2                basis(pm(grdj,grdj,2)),
     3                basis(pm(grdj,grdj,3)),
     4                basis(pm(grdj,grdj,4)),
     5                nmax(1,grdj),nmax(2,grdj),
     6                nmax(3,grdj),nmax(4,grdj),
     7                ntot(grdj),2,dim)
c
c           interpolate residual from grdi to this grid.
c
            call gi2gj(resid,rhs,
     1                 basis(pm(grdj,grdi,1)),basis(pm(grdj,grdi,2)),
     2                 basis(pm(grdj,grdi,3)),basis(pm(grdj,grdi,4)),
     3                 nmax(1,grdj),nmax(2,grdj),
     4                 nmax(3,grdj),nmax(4,grdj),
     5                 nmax(1,grdi),nmax(2,grdi),
     6                 nmax(3,grdi),nmax(4,grdi),
     7                 ntot(grdj),ntot(grdi),2,dim)       
            call vivo(rhs,rhs,
     1                 basis(pm(grdj,grdj,1)),basis(pm(grdj,grdj,2)),
     2                 basis(pm(grdj,grdj,3)),basis(pm(grdj,grdj,4)),
     3                 nmax(1,grdj),nmax(2,grdj),
     4                 nmax(3,grdj),nmax(4,grdj),
     5                 ntot(grdj),2,dim)
c
c           we now have a starting approximation for the grdj equation
c           and the residual.  call the davidson routine.
c
            title='solving residual equation for grid = '//itoc(grdj)
            call qcklin(soln(lsol2),rhs,basis,pm,array,ham,
     1                  voft,pvt,h0,h0ind,scale,cnverg,thresh,eps,
     2                  nmax,dim,spdim,msze,ntot,n3d,mgr,ntrial,
     3                  maxit,maxvec,mcycle,grdj,title,
     4                  prnt,mprnt)
            keep=lsol2
            lsol2 = lsol1
            lsol1 = keep 
 20      continue   
c
c        interpolate solution from grdi-1 to gridi
c
         call gi2gj(soln(lsol2),soln(lsol1),
     1              basis(pm(grdi,grdi-1,1)),
     2              basis(pm(grdi,grdi-1,2)),
     3              basis(pm(grdi,grdi-1,3)),
     4              basis(pm(grdi,grdi-1,4)),
     5              nmax(1,grdi),nmax(2,grdi),
     6              nmax(3,grdi),nmax(4,grdi),
     7              nmax(1,grdi-1),nmax(2,grdi-1),
     8              nmax(3,grdi-1),nmax(4,grdi-1),
     9              ntot(grdi),ntot(grdi-1),2,dim)
         call vivo(soln(lsol1),soln(lsol1),
     1             basis(pm(grdi,grdi,1)),
     2             basis(pm(grdi,grdi,2)),
     3             basis(pm(grdi,grdi,3)),
     4             basis(pm(grdi,grdi,4)),
     5             nmax(1,grdi),nmax(2,grdi),
     6             nmax(3,grdi),nmax(4,grdi),
     7             ntot(grdi),2,dim)
c
c        calculate the new residual
c
         call honv(array(ht),array(hx),array(hy),array(hz),
     1             voft(v),soln(lsol1),soln(lsol2),n,nt,nx,ny,nz,
     2             1,spdim)
         call vsub(resid,resid,soln(lsol2),m)         
 10   continue   
      call copy(soln(lsol1),vecout,m)
      if(mprnt(1)) then
         title='interpolated residual for next vector'
         call prntrm(title,vecout,m,1,m,1,iout)
      endif
      call memory(-ngot(1),pprhs,idum,'rhs',0)   
      call memory(-ngot(2),psoln,idum,'soln',0)   
      return
 1    format(/,5x,'multi-grid cycle = ',i3)
      end       

