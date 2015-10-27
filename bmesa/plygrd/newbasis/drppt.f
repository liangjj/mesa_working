*deck drppt.f
c***begin prologue     drppt
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       drppt
      subroutine drppt(pptwt,fval,ngrid,nin,nout)
      implicit integer (a-z)
      real*8 pt
      dimension nin(ngrid), nout(ngrid)
      common/io/inp, iout
      pointer (pptwt,pt(1))
      count=1
      do 10 grd=1,ngrid
         q=count
         wt=q+nin(grd)
         count=wt+nin(grd)
         if(fval.eq.0) then
            nout(grd)=nout(grd)-1
            call copy(pt(q+1),pt(q),nout(grd))
            call copy(pt(wt+1),pt(wt),nout(grd))   
            write(iout,1)
         endif
 10   continue   
      return
 1    format(/,5x,'dropping function at left endpoint')
      end



