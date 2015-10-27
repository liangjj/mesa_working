*deck filprj.f
c***begin prologue     filprj
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           three-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            fill array with value of polynomial on surface
c***references         
c
c***routines called    
c***end prologue       filprj
      subroutine filprj(p,srfp,where,nply,npts)
      implicit integer (a-z)
      real*8 p, srfp
      character*1 where
      dimension p(npts,nply), srfp(nply)
      common/io/inp, iout 
      if(where.eq.'r') then
         point=npts
      else
         point=1
      endif
      do 10 i=1,nply
         srfp(i) = p(point,i)
 10   continue   
      return
      end       

