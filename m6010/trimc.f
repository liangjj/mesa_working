*deck @(#)trimc.f	1.1 9/7/91
c***begin prologue     trimc
c***date written       910304   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6010, link 6010, kohn data
c***author             schneider, barry (lanl)
c***source             m6010
c***purpose            re-arrange optical potential matrices
c***                   for kohn codes.
c
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       trimc
      subroutine trimc(hin,hout,scr,list,n,m,nchan,nl2)
      implicit integer (a-z)
      real *8 hin, hout, scr
      dimension hin(n,n), hout(m,m), scr(m,m), list(nchan)
      common /io/ inp, iout
      inew=0
      do 10 i=1,nchan
         oldi=nl2*(i-1)+1  
         if (list(i).ne.0) then
             inew=inew+1
             newi=nl2*(inew-1)+1
             jnew=0
             do 20 j=1,nchan
                oldj=nl2*(j-1)+1
                if (list(j).ne.0) then
                    jnew=jnew+1
                    newj=nl2*(jnew-1)+1
                    do 30 k=1,nl2
                       call copy(hin(oldi,oldj),scr(newi,newj),nl2)
                       oldj=oldj+1 
                       newj=newj+1
   30               continue                     
                endif
   20        continue
         endif
   10 continue
      call copy(scr,hout,m*m)
      return
      end    
