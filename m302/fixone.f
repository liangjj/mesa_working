*deck @(#)fixone.f	5.1  11/6/94
      subroutine fixone(z,num,pkindx)
c***begin prologue     fixone.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             schneider, barry (lanl)
c***source             @(#)fixone.f	5.1   11/6/94
c***purpose            
c***description
c      repack the one-electron integrals according to
c      the drop information in pkindx.
c
c***references
c
c***routines called
c
c***end prologue       fixone.f
      implicit none
c     --- input variables -----
      integer num
c     --- input arrays (unmodified) ---
      integer pkindx(num)
c     --- input arrays (scratch) ---
      real*8 z(*)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer old,new,ix,i,j
      logical debug
c
      parameter (debug=.false.)
c
      common/io/inp,iout
c
c
      ix=0
      new=0
      do 1 i=1,num
         if(pkindx(i).ne.0) then
            old=ix
            do 2 j=1,i
               old=old+1
               if(pkindx(j).ne.0) then
                  new=new+1
                  z(new)=z(old)
               end if
  2         continue
         endif
         ix=ix+i
  1   continue
c
      if (debug) then
         write(iout,*) 'after fixone'
         write(iout,*) (z(i),i=1,new)
      endif
c
c
      return
      end
