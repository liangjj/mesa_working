*deck @(#)sumvw.f	
c***begin prologue     sumvw
c***date written       920402   (yymmdd)
c***revision date      
c                      bis at nsf june 25, 1992
c                      definition changed to reflect new projection operator
c                      factor of two multiplying fourth vlamda and subtraction
c                      of new wlamda term.
c***keywords           sumvw, link 6050
c***authors            Schneider,B (NSF)
c***                   
c***source             m6050
c***purpose            sum the vlamdas and wlamda
c***references       
c***routines called    
c***end prologue       sumvw
      subroutine sumvw (vlamda,wlamda,npt,nlam,prnt)
      implicit real *8 (a-h,o-z)
      complex *16 vlamda, wlamda
      character *80 title
      logical prnt
      dimension vlamda(npt,nlam,5), wlamda(npt,nlam)
      common /io/ inp, iout 
      do 10 lam=1,nlam
         do 20 i=1,npt
            vlamda(i,lam,1) = vlamda(i,lam,2) - vlamda(i,lam,3) -
     1                        vlamda(i,lam,4) + 2.d0*vlamda(i,lam,5) -
     2                        wlamda(i,lam)
   20    continue
   10 continue     
      if (prnt) then
          title='total wlamda'
          call prntcmn(title,vlamda,npt,nlam,npt,nlam,iout,'e')
      endif
      return
      end














