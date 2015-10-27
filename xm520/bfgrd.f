*deck @(#)bfgrd.f	5.1 6/20/95
      subroutine bfgrd(c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     nat,nprim,ntypes,nbtype,nnp,ncont,
     $     start,nbf,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,xyzgrid,mxgrd,xyzpow,scr,rsq,
     $     s,r,t,phi,grad,hess,dograd,dohess,mxcont,mxgbsiz,
     $     ng,ioff,maxl,bfcut)
c***begin prologue     bfgrd.f
c***date written       900513   (yymmdd)
c***revision date      6/20/95
c   may 1,  1995       rlm at lanl
c      adding cutoffs to speed up exponential evaluation
c   may 24, 1994       russo at lanl
c      modifying to compute hessian of basis functions.
c   may 13, 1994       rlm at lanl
c      modifying to return grad shaped as (mxgbsiz,nbf,3)
c
c***keywords           basis, grid, amplitude
c***author             martin, richard (lanl), RUSSO, thomas (lanl)
c***source             @(#)bfgrd.f	5.1  6/20/95
c***description
c     evaluate all of the contracted basis functions on a grid given by
c     xyzgrid.
c
c***references
c
c***routines called
c
c***end prologue       bfgrd.f
      implicit none
c
c     --- input variables -----
      integer nat,nprim,nbtype,nnp,ncont,ntypes,mxcont,mxgrd,mxgbsiz,
     $     ng,ioff,nbf
      logical dograd,dohess
      real*8 bfcut
c     --- input arrays (unmodified) ---
      real*8 c(3,nat),ex(nprim),cont(ncont)
      real*8 xyzgrid(mxgrd,3)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(0:*)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nx(*),ny(*),nz(*)
c     --- input arrays (scratch) ---
      real*8 xyzpow(mxgbsiz,3,*),scr(mxgbsiz),s(mxgbsiz,mxcont),
     $     r(mxgbsiz,mxcont),t(mxgbsiz,mxcont),rsq(mxgbsiz)
      integer maxl(nat)
c     --- output arrays ---
      real*8 phi(mxgbsiz,nbf),grad(mxgbsiz,nbf,3),hess(mxgbsiz,nbf,6)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer grid,iatom,itype,jatom,jtype,bf,power,contr,angmom,cart,
     $     prim
      integer n,l,m,pptr,contptr
      integer inp,iout
      logical timeit
      real*8 foo,zero,one,two,four
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      real*8 timpow,timexp,tim40,tim70
      real*8 cut

      common/io/inp,iout

      data timpow,timexp,tim40,tim70/4*0.d0/
      save timpow,timexp,tim40,tim70

      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00,four=4.0d+00)
      parameter (timeit=.false.)
c
c     --- convert the basis function cutoff into an exponent
      cut=log(bfcut)
c     --- determine, for each atom, the maximum angular momentum
c         of any shell defined for it.
      call izero(maxl,nat)
      do 10 iatom=1,nat
         do 5 itype=1,nbtype
            if(noprim(iatom,itype).ne.0) then
               maxl(iatom)=max(maxl(iatom),maxmom(itype))
            endif
    5    continue
c        --- allow for first derivatives
         if(dograd) then
            maxl(iatom)=maxl(iatom)+1
         endif
c        --- and for second derivatives
         if(dohess) then
            if (.not.dograd) then
               maxl(iatom)=maxl(iatom)+2
            else
               maxl(iatom)=maxl(iatom)+1
            endif
         endif
   10 continue
c
c     ----- the main loop -----
c     the plan is: for each atom, determine the amplitude of all basis 
c     functions on the grid for that atom. 
      call rzero(phi,mxgbsiz*nbf)
      if(dograd) then
         call rzero(grad,mxgbsiz*nbf*3)
      endif
      if (dohess) then
         call rzero(hess,mxgbsiz*6*nbf)
      endif
c
c     --- loop over all basis functions.
      bf=0
      do 90 jatom=1,nat
c
c        determine some preliminary quantities: 
c        distances and powers of distances from atom j to all
c        grid points on atom i.
         if(timeit) then
            call timing(dum1,dum2,dum3)
         endif
         do 20 grid=1,ng
            xyzpow(grid,1,1)=xyzgrid(ioff+grid,1)-c(1,jatom)
            xyzpow(grid,2,1)=xyzgrid(ioff+grid,2)-c(2,jatom)
            xyzpow(grid,3,1)=xyzgrid(ioff+grid,3)-c(3,jatom)
 20      continue
         do 30 power=2,max(maxl(jatom),2)
            call vmul(xyzpow(1,1,power),xyzpow(1,1,1),
     $                xyzpow(1,1,power-1),ng)
            call vmul(xyzpow(1,2,power),xyzpow(1,2,1),
     $                xyzpow(1,2,power-1),ng)
            call vmul(xyzpow(1,3,power),xyzpow(1,3,1),
     $                xyzpow(1,3,power-1),ng)
 30      continue
         call vadd(rsq,xyzpow(1,1,2),xyzpow(1,2,2),ng)
         call vadd(rsq,rsq,xyzpow(1,3,2),ng)
         if(timeit) then
            call timing(dum4,dum5,dum6)
            timpow=timpow+dum4-dum1
         endif
c
c        --- loop over the shell types on this atom; s,p,d,sp,etc.
         do 80 jtype=1,nbtype
c           are there any primitive functions on this center of this type
            if (noprim(jatom,jtype).gt.0) then
               call rzero(s,mxgbsiz*mxcont)
               if(dograd) then
                  call rzero(r,mxgbsiz*mxcont)
               endif
               if (dohess) then
                  call rzero(t,mxgbsiz*mxcont)
               endif
c
c              --- loop over all the primitives in this shell.
               pptr=0
               do 50 prim=ptprim(jatom,jtype),
     $                    ptprim(jatom,jtype)+noprim(jatom,jtype)-1
                  if(timeit) then
                     call timing(dum1,dum2,dum3)
                  endif
                  call smul(scr,rsq,ex(prim),ng)
                  call vneg(scr,scr,ng)
c                 evaluate only the functions in the vicinity.
                  do 45 grid=1,ng
                     if(scr(grid).ge.cut) then
                        scr(grid)=exp(scr(grid))
                     else
                        scr(grid)=zero
                     endif
   45             continue
                  if(timeit) then
                     call timing(dum4,dum5,dum6)
                     timexp=timexp+dum4-dum1
                  endif
c
c                 --- loop over all the contractions in this shell.
                  if(timeit) then
                     call timing(dum1,dum2,dum3)
                  endif
                  contptr=0
                  do 40 contr=1,nocont(jatom,jtype)
                     call vwxs(s(1,contr),s(1,contr),scr,
     $                         cont(ptcont(jatom,jtype)+contptr+pptr),
     $                         +1,ng)
                     if(dograd) then
                        foo=two*
     $                     cont(ptcont(jatom,jtype)+contptr+pptr)
     $                     *ex(prim)
                        call vwxs(r(1,contr),r(1,contr),scr,foo,+1,ng)
                     endif
                     if (dohess) then
                        foo=four*cont(ptcont(jatom,jtype)+contptr+pptr)
     $                       *ex(prim)*ex(prim)
                        call vwxs(t(1,contr),t(1,contr),scr,foo,+1,ng)
                     endif
                     contptr=contptr+noprim(jatom,jtype)
 40               continue
                  pptr=pptr+1
                  if(timeit) then
                     call timing(dum4,dum5,dum6)
                     tim40=tim40+dum4-dum1
                  endif
 50            continue
c
               if(timeit) then
                  call timing(dum1,dum2,dum3)
               endif
               do 70 contr=1,nocont(jatom,jtype)
c                 loop over all the angular momentum components,
c                 e.g. s and p for an sp shell.
                  do 60 angmom=minmom(jtype),maxmom(jtype)
c                    loop over all the individual cartesian components
c                    of this type.
                     do 55 cart=1,nocart(angmom)
                        n=nx(mintyp(jtype)+cart-1)
                        l=ny(mintyp(jtype)+cart-1)
                        m=nz(mintyp(jtype)+cart-1)
                        bf=bf+1
c                       load the exponential portion.
                        call vmove(phi(1,bf),s(1,contr),ng)
c                       apply the polynomial part to the basis
c                       and form some of the gradient polynomials
                        if(n.ne.0) then
                           call vmul(phi(1,bf),phi(1,bf),
     $                               xyzpow(1,1,n),ng)
                        endif
                        if(l.ne.0) then
                           call vmul(phi(1,bf),phi(1,bf),
     $                              xyzpow(1,2,l),ng)
                        endif
                        if(m.ne.0) then
                           call vmul(phi(1,bf),phi(1,bf),
     $                               xyzpow(1,3,m),ng)
                        endif
c
c                       --- do the derivatives
                        if(dograd) then
                           call dgrad(n,l,m,bf,contr,grad,s,r,xyzpow,ng,
     $                                nbf,mxgbsiz,mxcont)
                        endif
c
c                       --- second derivatives
                        if(dohess) then
                           call dhess(n,l,m,bf,contr,hess,s,r,t,xyzpow,
     $                                scr,ng,nbf,mxgbsiz,mxcont)
                        endif
 55                  continue
 60               continue
 70            continue
               if(timeit) then
                  call timing(dum4,dum5,dum6)
                  tim70=tim70+dum4-dum1
               endif
            endif
 80      continue
 90   continue
c
c     --- check the number of basis functions processed in the loop.
      if(bf.ne.nbf) then
         write(iout,*)'nbf=',nbf,' bf = ',bf
         call  lnkerr(
     $        'bfgrd: calculated bf not same as passed nbf.  Bleah.')
      endif
c
      if(timeit) then
         write(iout,*) 'timpow,timexp,tim40,tim70',
     $                  timpow,timexp,tim40,tim70
      endif
c
c
      return
      end
*deck @(#)bfgrd.f	5.1 6/20/95
      subroutine dgrad(n,l,m,bf,contr,grad,s,r,xyzpow,ng,
     $                 nbf,mxgbsiz,mxcont)
c***begin prologue     bfgrd.f
c***date written       yymmdd  
c***revision date      6/20/95      
c
c***keywords           
c***author             
c***source             @(#)bfgrd.f	5.1   6/20/95
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       bfgrd.f
      implicit none
c     --- input variables -----
      integer n,l,m,bf,contr,ng
      integer nbf,mxgbsiz,mxcont
c     --- input arrays (unmodified) ---
      real*8 s(mxgbsiz,mxcont),r(mxgbsiz,mxcont)
      real*8 xyzpow(mxgbsiz,3,*)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 grad(mxgbsiz,nbf,3)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      real*8 foo
c
c     --- do the derivatives
      if(n.ne.0) then
         if (n .ne. 1) then
             call vmul(grad(1,bf,1),s(1,contr),xyzpow(1,1,n-1),ng)
             foo=float(n)
             call smul(grad(1,bf,1),grad(1,bf,1),foo,ng)
         else
             call vmove(grad(1,bf,1),s(1,contr),ng)
         endif
      endif
      call vwxy(grad(1,bf,1),grad(1,bf,1),xyzpow(1,1,n+1),
     $          r(1,contr),-1,ng)
c
      if(l.ne.0) then
         if (l .ne. 1) then
            call vmul(grad(1,bf,2),s(1,contr),xyzpow(1,2,l-1),ng)
            foo=float(l)
            call smul(grad(1,bf,2),grad(1,bf,2),foo,ng)
         else
            call vmove(grad(1,bf,2),s(1,contr),ng)
         endif
      endif
      call vwxy(grad(1,bf,2),grad(1,bf,2),xyzpow(1,2,l+1),
     $          r(1,contr),-1,ng)
c
      if(m.ne.0) then
          if (m .ne. 1) then
              call vmul(grad(1,bf,3),s(1,contr),xyzpow(1,3,m-1),ng)
              foo=float(m)
              call smul(grad(1,bf,3),grad(1,bf,3),foo,ng)
           else
              call vmove(grad(1,bf,3),s(1,contr),ng)
           endif
      endif
      call vwxy(grad(1,bf,3),grad(1,bf,3),xyzpow(1,3,m+1),
     $          r(1,contr),-1,ng)
c
c     ---finish off the gradients
      if (n .ne. 0) then
         call vmul(grad(1,bf,2),grad(1,bf,2),xyzpow(1,1,n),ng)
         call vmul(grad(1,bf,3),grad(1,bf,3),xyzpow(1,1,n),ng)
      endif
      if (l .ne. 0) then
         call vmul(grad(1,bf,1),grad(1,bf,1),xyzpow(1,2,l),ng)
         call vmul(grad(1,bf,3),grad(1,bf,3),xyzpow(1,2,l),ng)
      endif
      if (m .ne. 0) then
         call vmul(grad(1,bf,1),grad(1,bf,1),xyzpow(1,3,m),ng)
         call vmul(grad(1,bf,2),grad(1,bf,2),xyzpow(1,3,m),ng)
      endif
c
c
      return
      end
*deck @(#)bfgrd.f	5.1 6/20/95
      subroutine dhess(n,l,m,bf,contr,hess,s,r,t,xyzpow,scr,ng,
     $                 nbf,mxgbsiz,mxcont)
c***begin prologue     bfgrd.f
c***date written       yymmdd  
c***revision date      6/20/95      
c
c***keywords           
c***author             
c***source             @(#)bfgrd.f	5.1   6/20/95
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       bfgrd.f
      implicit none
c     --- input variables -----
      integer n,l,m,bf,contr,ng
      integer nbf,mxgbsiz,mxcont
c     --- input arrays (unmodified) ---
      real*8 s(mxgbsiz,mxcont),r(mxgbsiz,mxcont),t(mxgbsiz,mxcont)
      real*8 xyzpow(mxgbsiz,3,*)
c     --- input arrays (scratch) ---
      real*8 scr(ng)
c     --- output arrays ---
      real*8 hess(mxgbsiz,nbf,6)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer fribble,d1,d2,d3,p1,p2,p3,grid
      real*8 foo,two,three,five
c
      data two/2.0d0/,three/3.0d0/,five/5.0d0/
c
c     We'll lay out the hessian like this:
c        hess(.,.,i) has:
c        i=1: d2phi/dx2,
c        i=2: d2phi/dy2,
c        i=3:d2phi/dz2
c        i=4: d2phi/dxdy
c        i=5:d2phi/dxdz,
c        i=6:d2phi/dydz
c
      fribble=1
      do 50 d1=1,3
         if (d1.eq.1) then
            p1=n
            p2=l
            d2=2
            p3=m
            d3=3
         else if (d1.eq.2) then
            p1=l
            p2=m
            d2=3
            p3=n
            d3=1
         else
            p1=m
            p2=n
            d2=1
            p3=l
            d3=2
         endif
         call vmove(scr,t(1,contr),ng)
         call vmul(scr,scr,xyzpow(1,d1,p1+2),ng)
         if (p1.ne.0) then
            if (p1.ne.1) then
               if (p1.eq.2) then
c                 p1==2                               
                  do 10 grid=1,ng
                     scr(grid)=scr(grid)+two*s(grid,contr)
     $                        -five*r(grid,contr)*xyzpow(grid,d1,2)
   10             continue 
               else
c                 p1>2
                  do 20 grid=1,ng
                     scr(grid)=scr(grid)+p1*(p1-1)*
     $                         s(grid,contr)*xyzpow(grid,d1,p1-2)
     $                         -(2*p1+1)*r(grid,contr)*
     $                         xyzpow(grid,d1,p1)
   20             continue 
               endif
            else
c              p1==1
               do 30 grid=1,ng
                  scr(grid)=scr(grid)-three*r(grid,contr)
     $                                *xyzpow(grid,d1,1)
   30          continue 
            endif
         else
c           p1==0
            do 40 grid=1,ng
               scr(grid)=scr(grid)-r(grid,contr)
   40       continue 
         endif
         if (p2.ne.0) then
            call vmul(scr,scr,xyzpow(1,d2,p2),ng)
         endif
         if (p3.ne.0) then
            call vmul(scr,scr,xyzpow(1,d3,p3),ng)
         endif
         call vmove(hess(1,bf,fribble),scr,ng)
         fribble=fribble+1
  50  continue 
c
c     --- now the mixed partials:
c
c     Ok, so I chose an unreadable way to do it.  We'll generate
c     d2/dxdy, d2/dxdz and d2/dydz in that order.  p1, p2, and p3 are 
c     n,l,m appropriately permuted so we can form the terms we need without
c     actually coding term-specific special cases, which is what I started
c     to do but which quickly became ridiculously unwieldy.
c     Fribble is the offset into hess which will receive the values.
c
c     the general case is:
c     deriv w.r.t. d1&d2 = (omitting grid index, etc)
c     [ p1*p2*s*xyzpow(d1,p1-1)*xyzpow(d2,p2-1)
c     -(p1*xyzpow(d1,p1-1)*xyzpow(d2,p2+1)+p2*xyzpow(d1,p1+1)*xyzpow(d2,p2-1)*r
c     +t*xyzpow(d1,p1+1)*xyzpow(d2,p2+1) ]*xyzpow(d3,p3)
c 
c     obviously, some special cases must be taken care of to avoid accessing
c     xyzpow(...,-1) etc, and to avoid unnecessary flops.
c
      fribble=4
      do 100 d1=1,2
         p1=n
         if (d1.eq.2) p1=l
         do 90 d2=d1+1,3
            p2=l
            if (d2.eq.3) p2=m
            p3=m
            d3=3
            if (d1.eq.1 .and. d2.eq.3) then
               p3=l
               d3=2
            endif
            if (d1.eq.2 .and. d2.eq.3) then
                p3=n
                d3=1
            endif
c           third term above
            call vmul(scr,t(1,contr),xyzpow(1,d1,p1+1),ng)
            call vmul(hess(1,bf,fribble),scr,xyzpow(1,d2,p2+1),ng)
c           there are two "second terms",
c           and each one is nonzero only if its
c           respective coeff is nonzero
            call rzero(scr,ng)
            if (p1.ne.0) then
               if (p1.ne.1) then
                  call vmul(scr,xyzpow(1,d1,p1-1),xyzpow(1,d2,p2+1),ng)
                  foo=float(p1)
                  call smul(scr,scr,foo,ng)
               else
                  call vmove(scr,xyzpow(1,d2,p2+1),ng)
               endif
            endif
            if (p2 .ne.0) then
               if (p2 .ne.1) then
                  do 80 grid=1,ng
                     scr(grid)=scr(grid)
     $                        +p2*xyzpow(grid,d1,p1+1)
     $                           *xyzpow(grid,d2,p2-1)
   80             continue 
               else
                  call vadd(scr,scr,xyzpow(1,d1,p1+1),ng)
               endif
            endif
            call vmul(scr,scr,r(1,contr),ng)
            call vsub(hess(1,bf,fribble),hess(1,bf,fribble),scr,ng)
c
c           first term only nonzero in same case
            if (p1.ne.0 .and. p2.ne.0) then
               call vmove(scr,s(1,contr),ng)
               foo=p1*p2
               call smul(scr,scr,foo,ng)
               if (p1 .ne. 1) then
                  call vmul(scr,scr,xyzpow(1,d1,p1-1),ng)
               endif
               if (p2.ne.1) then
                  call vmul(scr,scr,xyzpow(1,d2,p2-1),ng)
               endif
               call vadd(hess(1,bf,fribble),hess(1,bf,fribble),scr,ng)
            endif
            if (p3 .ne. 0) then
                call vmul(hess(1,bf,fribble),hess(1,bf,fribble),
     $                    xyzpow(1,d3,p3),ng)
            endif
         fribble=fribble+1
   90    continue 
  100 continue 
c
c
      return
      end
