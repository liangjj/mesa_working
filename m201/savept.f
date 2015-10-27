*deck @(#)savept.f	5.1  11/6/94
      subroutine savept(nvar,nvv,mxpt,rises,ff,xx,f,x,fc,frcnst,ic,
     $                  fs,tempxx,tempff,srcd2e,np,cycle,neg,
     $                  ipsav1,ipsav2,energy,fncerr)
c***begin prologue     savept.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)savept.f	5.1   11/6/94
c***purpose            
c***description
c     
c     original version by h.b.schlegel, rewritten by r.l.martin.
c     push the new energy, position, forces, and possibly force
c     constants onto the stacks.  be careful to keep the best point
c     at the top of the stack.  rises is returned .true. if the
c     energy(neg=0) or gradient(neg=1) of the new point is not
c     the lowest so far.
c
c***references
c
c***routines called
c
c***end prologue       savept.f
      implicit none
c     --- input variables -----
      integer nvar,nvv,mxpt,np,cycle,neg,ipsav1,ipsav2
      character*(*) srcd2e
      real*8 energy,fncerr
c     --- input arrays (unmodified) ---
      real*8 fc(2*nvv),frcnst(nvv)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer ic(mxpt)
      real*8 fs(mxpt)
c     --- output variables ---
      logical rises
c     --- scratch arrays ---
      real*8 ff(nvar,mxpt),xx(nvar,mxpt),f(nvar),x(nvar)
      real*8 tempxx(nvar),tempff(nvar)
c     --- local variables ---
      integer inp,iout
      integer ip,i,iip,ifrom,ito,ictemp
      real*8 sdot,gni,gnnew,gnold,gnsave,temp
c
      common/io/inp,iout
c
c
      np=min(np+1,mxpt)
      if(np.ge.2) then
         do 10 iip=2, np
            ifrom=np+1-iip
            ito=ifrom+1
            fs(ito)=fs(ifrom)
            ic(ito)=ic(ifrom)
            call vmove(ff(1,ito),ff(1,ifrom),nvar)
            call vmove(xx(1,ito),xx(1,ifrom),nvar)
   10    continue
      endif
c
c
      fs(1)=energy
      ic(1)=cycle
      call vmove(xx(1,1),x,nvar)
      call vmove(ff(1,1),f,nvar)
c
c     --- make sure that the best point other than the current
c         one is second from the top of the stack.
      if(np.gt.2) then
         ip=2
         gnsave=sdot(nvar,ff(1,ip),1,ff(1,ip),1)
         do 20 i=3,np
            gni=sdot(nvar,ff(1,i),1,ff(1,i),1)
            if(neg.eq.0) then
               if(fs(i).lt.fs(ip)) ip=i
            else
               if(gni.lt.gnsave) ip=i
            endif
            gnsave=min(gni,gnsave)
  20     continue
         if(ip.ne.2) then
            temp=fs(ip)
            ictemp=ic(ip)
            call vmove(tempxx,xx(1,ip),nvar)
            call vmove(tempff,ff(1,ip),nvar)
            do 30 i=3,ip
               ito=ip+3-i
               ifrom=ito-1
               fs(ito)=fs(ifrom)
               ic(ito)=ic(ifrom)
               call vmove(xx(1,ito),xx(1,ifrom),nvar)
               call vmove(ff(1,ito),ff(1,ifrom),nvar)
   30       continue
            fs(2)=temp
            ic(2)=ictemp
            call vmove(xx(1,2),tempxx,nvar)
            call vmove(ff(1,2),tempff,nvar)
         endif
      endif
c
      if(ic(2).eq.ipsav1) then
         call vmove(fc(nvv+1),fc(1),nvv)
         ipsav2=ipsav1
      endif
      ipsav1=0
      if(ic(2).ne.ipsav2) ipsav2=0
      if(srcd2e.ne.' ') then
         call vmove(fc,frcnst,nvv)
         ipsav1=np
      endif
c
c
      rises=.false.
      if(np.gt.1) then
         gnnew=sdot(nvar,ff(1,1),1,ff(1,1),1)
         gnold=sdot(nvar,ff(1,2),1,ff(1,2),1)
         rises=(neg.eq.0.and.(fs(1)-fncerr).gt.(fs(2)+fncerr))
     $         .or.(neg.ne.0.and.(gnnew.gt.gnold))
      endif
c
c
      return
      end
