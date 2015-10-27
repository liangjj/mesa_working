*deck @(#)fixone.f	5.1  11/6/94
      subroutine fixone(a,num,pkindx)
c
c      repack the one-electron integrals according to
c      the drop information in pkindx.
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
