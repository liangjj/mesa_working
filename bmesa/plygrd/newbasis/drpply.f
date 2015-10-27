*deck drpply.f
c***begin prologue     drpply
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
c***end prologue       drpply
      subroutine drpply(ppoly,pptmp,fval,ngrid,nin,nout,which)
      implicit integer (a-z)
      real*8 poly, ptmp
      character*(*) which
      dimension nin(ngrid), nout(ngrid)
      common/io/inp, iout
      pointer (ppoly,poly(1))
      pointer (pptmp,ptmp(1))
      if(which.eq.'first') then
         count=1
         do 10 grdi=1,ngrid
            do 20 grdj=1,ngrid
               p=count
               dp=p+nin(grdi)*nin(grdj)
               ddp=dp+nin(grdi)*nin(grdj)
               count=ddp+nin(grdi)*nin(grdj)         
               if(fval.eq.0) then
c
c                 remove the first function and first component
c                 and copy back into original location
c    
                  call trmply(poly(p),ptmp,nin(grdi),nin(grdj),
     1                        nout(grdi),nout(grdj),which)  
                  call copy(ptmp,poly(p),nout(grdi)*nout(grdj))
                  call trmply(poly(dp),ptmp,nin(grdi),nin(grdj),  
     1                        nout(grdi),nout(grdj),which)  
                  call copy(ptmp,poly(dp),nout(grdi)*nout(grdj))
                  call trmply(poly(ddp),ptmp,nin(grdi),nin(grdj),
     1                        nout(grdi),nout(grdj),which)  
                  call copy(ptmp,poly(ddp),nout(grdi)*nout(grdj))
               endif
 30         continue   
 20      continue      
c
c        change nin array to reflect any decrease in number of functions.
c
         call icopy(nout,nin,ngrid)
      elseif(which.eq.'last') then
         do 40 grdi=1,ngrid
            if(fval.eq.0) then
               nout(grdi)=nout(grdi)-1
            endif
 40      continue
         count=1   
         do 50 grdi=1,ngrid
            do 60 grdj=1,ngrid
               p=count
               dp=p+nin(grdi)*nin(grdj)
               ddp=dp+nin(grdi)*nin(grdj)
               count=ddp+nin(grdi)*nin(grdj) 
               call trmply(poly(p),ptmp,nin(grdi),nin(grdj),
     1                     nout(grdi),nout(grdj),which)
               call copy(ptmp,poly(p),nout(grdi)*nout(grdj))
               call trmply(poly(dp),ptmp,nin(grdi),nin(grdj),
     1                     nout(grdi),nout(grdj),which)
               call copy(ptmp,poly(dp),nout(grdi)*nout(grdj))
               call trmply(poly(ddp),ptmp,nin(grdi),nin(grdj),
     1                     nout(grdi),nout(grdj),which)
               call copy(ptmp,poly(ddp),nout(grdi)*nout(grdj))
 60         continue   
 50      continue               
         call icopy(nout,nin,ngrid)
      endif
      return
 1    format(/,5x,'dropping function at left endpoint')
      end



