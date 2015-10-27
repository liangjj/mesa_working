*deck drpwt.f
c***begin prologue     drpwt
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
c***end prologue       drpwt
      subroutine drpwt(ppt,fl,nin,nout)
      implicit integer (a-z)
      real*8 pt
      common/io/inp, iout
      pointer (ppt,pt(1))
      nout=nin
      q=1
      wt=q+nin
      qm=wt+nin
      wtm=qm+nin  
      if(fl.ne.0) then
         call copy(pt(q),pt(qm),nin)
         call copy(pt(wt),pt(wtm),nin)   
      elseif 
         nout=nin-1
         call copy(pt(q+1),pt(qm),nout)
         call copy(pt(wt+1),pt(wtm),nout)   
         write(iout,1)
      endif
      return
 1    format(/,5x,'dropping function at left endpoint')
      end
