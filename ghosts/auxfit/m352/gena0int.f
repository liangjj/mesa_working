*deck %W% %G%
      subroutine gena0int(iatom,jatom,imax,jmax,nprimi,nprimj,
     $     i1,j1,ex,alpha,c,ainv,prctr,kappa,f,expon,wma,outarry,
     $     tmparry,tmpsiz,t,nv,nprim,nat,mmax,pi252,lval,mycart,yn0,
     $     amaxm,mindx)
c***begin prologue     %M%
c***date written       930129 (yymmdd)
c***revision date      %G%
c
c***keywords
c***author             russo, thomas (lanl)
c***source             %W% %G%
c***purpose
c     Generate two-center, two-electron integrals of form (A,0)
c     which will be used to build other types of two-electron integrals
c***description
c    uses basic obara-saika recurrences, taking advantage of the two-center
c    special case and using translational invariance to build the (a,c+1)
c    recurrence
c
c***references
c     obara and saika, j.chem.phys., 84,3963 (1986)
c     
c***routines called
c
c***end prologue     %M%
      implicit integer(a-z)
c
      logical debug
      real*8 half,one,two
      parameter (debug=.false.)
      parameter (half=0.5d0,one=1.0d0,two=2.0d0)
      common /io/ inp,iout
c
c
c
      integer toff,noff,loffm,loffp
      real*8 jx,ix
      real*8 pi252
c input arrays (unmodified)
      real*8 ex(nprim),c(3,nat)
cScratch arrays:
c
c a place to build the integrals.  we do 0:2 because at each step
c one of the places will hold the current generation, one will hold
c the last, and we'll build the next generation into the third.
c I don't like that 16, but that should be ok, that's k+k functions!
      real*8 tmparry(nv,tmpsiz,0:2),tmp(0:16),t(nv),f(nv,0:mmax),
     $     kappa(nv)
c
cOutput arrays:
c a place to put the (a,0)^(0) as they're made
c
      real*8 outarry(nv,*)
      real*8 alpha(nv,2),ainv(nv,5),prctr(nv,3),wma(nv,3)
      real*8 scalri,scalar
c
      integer lval(0:amaxm,mindx,3),yn0(0:amaxm),mycart(0:amaxm)
c------------------------------------------------------------------
c enough of that!
c Define a statement function to compute the offset of the (l,0)^(m)
c------------------------------------------------------------------
      intoff(l,m)=(m)*((l+1)*(l+2))/2
c
c Now, Rock and Roll!
c
      last=2
      this=0
      next=1
      putit=0
      if (debug) then
         write(iout,*) "****************"
         write(iout,*) "Doing (",imax,"|",jmax,") integrals"
      endif
c
c collect the pair data
c
      i2=i1+nprimi-1
      j2=j1+nprimj-1
      n=0
      do 10 jprim=j1,j2
         jx=ex(jprim)
         do 20 iprim=i1,i2
            ix=ex(iprim)
            n=n+1
            alpha(n,1)=ix
            alpha(n,2)=jx
 20      continue 
 10   continue 
      if(debug) then
         write(iout,*) 'alpha'
         write(iout,*) (alpha(i,1),i=1,n)
         write(iout,*) (alpha(i,2),i=1,n)
      endif
      if (n.ne.nv) then
         call lnkerr('gena0int: error in vector lengths')
      endif
c
c ainv(.,1)=1/alpha1
c ainv(.,2)=1/alpha2
c ainv(.,3)=1/(alpha1+alpha2)
c ainv(.,4)=alpha1/(alpha1+alpha2)
c ainv(.,5)=alpha2/(alpha1+alpha2)
c
      call vinv(ainv(1,1),alpha(1,1),n)
      call vinv(ainv(1,2),alpha(1,2),n)
      call vadd(ainv(1,3),alpha(1,1),alpha(1,2),n)
      call vinv(ainv(1,3),ainv(1,3),n)
      call vmul(ainv(1,4),alpha(1,1),ainv(1,3),n)
      call vmul(ainv(1,5),alpha(1,2),ainv(1,3),n)
      if(debug) then
         write(iout,*) 'ainv'
         write(iout,*) (ainv(ij,1),ij=1,n)
         write(iout,*) (ainv(ij,2),ij=1,n)
         write(iout,*) (ainv(ij,3),ij=1,n)
         write(iout,*) (ainv(ij,4),ij=1,n)
         write(iout,*) (ainv(ij,5),ij=1,n)
      endif
c
c     ----- form the pair center W, (ix*xi + jx*xj)/(ix+jx) -----
c
      do 30 coord=1,3
         call smul(prctr(1,coord),alpha(1,1),c(coord,iatom),n)
         call saxpy(n,c(coord,jatom),alpha(1,2),1,prctr(1,coord),1)
         call vmul(prctr(1,coord),prctr(1,coord),ainv(1,3),n)
   30 continue
      if(debug) then
         write(iout,*) 'natural center'
         write(iout,*) (prctr(ij,1),ij=1,n)
         write(iout,*) (prctr(ij,2),ij=1,n)
         write(iout,*) (prctr(ij,3),ij=1,n)
      endif
c
c     ----- form exponential parameter T -----
c
      call vmul(expon,alpha(1,1),alpha(1,2),n)
      call vmul(t,expon,ainv(1,3),n)
      scalar=((c(1,iatom)-c(1,jatom))**2+(c(2,iatom)-c(2,jatom))**2+
     $         (c(3,iatom)-c(3,jatom))**2)
      call smul(t,t,scalar,n)
      if(debug) then
         write(iout,*) 't:'
         call matout(t,nprimi,nprimj,nprimi,nprimj,iout)
      write(iout,*) 'ainv3',(ainv(ij,3),ij=1,n)
      endif
c
c     ----- form exponential prefactor -----
      call vinv(kappa,expon,n)
      call vsqrt(expon,ainv(1,3),n)
      call vmul(kappa,kappa,expon,n)
      call smul(kappa,kappa,pi252,n)
      if(debug) then
         write(iout,*) 'kappa:'
         write(iout,*) (kappa(ij),ij=1,n)
      endif
c
c     ----- form W-A -----
      do 450 coord=1,3
         call ssub(wma(1,coord),prctr(1,coord),c(coord,iatom),n)
 450  continue 
      if (debug) then
         write(iout,*)'wma:'
         write(iout,*) ((wma(j,coord),coord=1,3),j=1,n)
      endif
c
c     ----- get the f(m) table for this primitive block. -----
c
c     the dimensioning here is awkward.
      do 50 i=1,n
         call fmoft(tmp,mmax,t(i))
         do 40 m=0,mmax
            f(i,m)=tmp(m)
   40    continue
   50 continue
      if(debug) then
         write(iout,*) 'f(m,t):'
         do 60 m=0,mmax
            write(iout,*) 'm:',m
            call matout(f(1,m),nprimi,nprimj,nprimi,nprimj,iout)
   60    continue
      endif
c
c     ----- form necessary (s|s)(m) integrals to begin recursion
c
      do 70 m=0,mmax
         call vmul(t,kappa,f(1,m),n)
         call vmove(tmparry(1,m+1,this),t,n)
   70 continue
      if(debug) then
         write(iout,*) '(0,0)(m) integrals'
         write(iout,*) ((tmparry(ij,m+1,this),ij=1,n),m=0,mmax)
      endif


c*
c Form all (a,0)^(0) integrals, where a is a triple of integers.
c Start with (0,0)^(m) and build them up using
c 
c (a+1_i,0)^(m)=(W_i-A_i)(a,0)^(m+1)+1/(2\alpha)N_i(a) 
c        *[(a-1_i,0)^(m) - \rho/\alpha (a-1_i,0)^(m)]
c
c where N_i(a) is the i'th component of a, ie N_x(a_x,a_y,a_z)=a_x
c and \rho=\alpha\beta/(\alpha+\beta)
c
c la is the current generation's total angular momentum, it will be
c bumped up by one.
c
c oh yeah, and copy the (0,0)^(0) to the output array.  In general the
c (a,0)^(0) will go into outarry(.,putit) through outarry(.,putit+mycart(la)-1)
c
      putit=1
      call vmove(outarry(1,putit),tmparry(1,1,this),n)
      putit=2
      do 1000 la=0,mmax-1
         call rzero(tmparry(1,1,next),n*tmpsiz)
         do 1010 m=0,mmax-la-1
c
c figure out where to get the (la,0)^(m+1), etc.
c also, where to put the (la+1_x,0), (la+1_y,0) and (la+1_z,0)
c
            toff=intoff(la,m+1)
c
            loffm=intoff(la-1,m)
c
            loffp=intoff(la-1,m+1)
c
            noff=intoff(la+1,m)
c
c  add 1_x to the  value of all of the (a,0)'s, 
c  add 1_y to the values of the (a,0) with ax=0
c  add 1_z to the value of the (a,0) with ax=0,ay=0
            do 1020 coord=1,3
               if (coord.eq.1) then
                  minint=1
                  minm1=1
               else if (coord.eq.2) then
                  minint=yn0(la)+1
                  minm1=yn0(la-1)+1
               else
                  minint=mycart(la)
                  minm1=mycart(la-1)
               endif
               do 1030 theint=minint,mycart(la)
                  them1=theint-minint+minm1
                  asubi=lval(la,theint,coord)
                  if (debug) then
                     write(iout,*)"bumping coord ",coord,
     $                    "on  (",lval(la,theint,1),",",
     $                    lval(la,theint,2),",",lval(la,theint,3),")"
                     write(iout,*)"The (",la,")(",m+1,") are in ",
     $                    theint+toff
                     write(iout,*)(tmparry(i,theint+toff,this),
     $                    i=1,n)
                     write(iout,*) "putting in ",noff+1
                     if (asubi .ne.0) then
                     write(iout,*) "the (",la-1,")(",m,")'s are in ",
     $                    them1+loffm
                     write(iout,*) (tmparry(i,them1+loffm,last),
     $                    i=1,n)
                     write(iout,*) "the (",la-1,")(",m+1,")'s are in ",
     $                    them1+loffp
                     write(iout,*) (tmparry(i,them1+loffp,last),
     $                    i=1,n)
                     endif
                  endif
c
c form (w_i-a_i)(la,0)(m+1)+1/2alpha1 a_x[(a-1_i)(m)-rho/alpha1(a-1_i)(m+1)]
c also, rho/alpha1=alpha2/(alpha1+alpha2)=ainv(.,5)
c
                  call vwxy(tmparry(1,noff+1,next),
     $                 tmparry(1,noff+1,next),wma(1,coord),
     $                 tmparry(1,theint+toff,this),+1,n)
                  if (asubi .ne.0) then
                     scalri=float(asubi)*half
                     call smul(t,ainv(1,1),scalri,n)
                     call vwxy(tmparry(1,noff+1,next),
     $                    tmparry(1,noff+1,next),t,
     $                    tmparry(1,them1+loffm,last),+1,n)
                     call vmul(t,t,ainv(1,5),n)
                     call vwxy(tmparry(1,noff+1,next),
     $                    tmparry(1,noff+1,next),t,
     $                    tmparry(1,them1+loffp,last),-1,n)
                  endif
                  if (debug) then
                     ax=lval(la,theint,1)
                     ay=lval(la,theint,2)
                     az=lval(la,theint,3)
                     if (coord .eq.1)ax=ax+1
                     if (coord .eq.2)ay=ay+1
                     if (coord .eq.3)az=az+1
                     write(iout,*) "and the winner is (",
     $                    ax,",",ay,",",az,")(",m,
     $                    "):"
                     write(iout,*) (tmparry(i,noff+1,next),
     $                    i=1,n)
                  endif
                  noff=noff+1
 1030          continue 
 1020       continue 
c
c we now have done the (a,0)^(m).  If m's zero copy over to the
c output array
c
            if (m.eq.0) then
               call vmove(outarry(1,putit),tmparry(1,1,next),
     $              n*mycart(la+1))
               putit=putit+mycart(la+1)
            endif
 1010    continue 
c
c now we're done with the current generation.  Diddle with the pointers
c to change the meaning of next, last, this
c
         temp=last
         last=this
         this=next
         next=temp
c
c now go do the next la
c
 1000 continue 

      if (debug) then
         theoff=0
         do 666 la=0,mmax
            write(iout,*)"The (",la,",0) integrals are (first prim):"
            write(iout,*) (outarry(1,theoff+i),i=1,mycart(la))
            theoff=theoff+mycart(la)
            write(iout,*)"------------"
 666     continue 
      endif
      return
      end
