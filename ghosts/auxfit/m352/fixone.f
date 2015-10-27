*deck  @(#)fixone.f	2.1 4/28/92
      subroutine fixone(a,num,pkindx)
c
c***begin prologue     fixone.f
c***date written       910101  
c***revision date      4/28/92      
c
c***keywords           
c***author             schneider, barry(lanl) 
c***source             @(#)fixone.f	2.1   4/28/92
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
c
c
      implicit integer(a-z)
      integer pkindx(num),old
      real*8 a(*)
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
                  a(new)=a(old)
               end if
  2         continue
         endif
         ix=ix+i
  1   continue
c
      if (debug) then
         write(iout,*) 'after fixone'
         write(iout,*) (a(i),i=1,new)
      endif
c
c
      return
      end
