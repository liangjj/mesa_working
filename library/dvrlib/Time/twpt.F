*deck twpt.f
c***begin prologue     twpt
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            lobatto points, weights, polynomials
c***                   and their first and second derivatives.
c***                   
c***references         
c
c***routines called    
c***end prologue       twpt
      subroutine twpt(q,wt,p,dp,ddp,qdtyp,edge,n,nq,prn)
      implicit integer (a-z)
      real*8 q, wt, p, dp, ddp, edge
      logical prn
      character*80 title
      character*(*) qdtyp
      dimension q(n), wt(n), p(n,n), dp(n,n), ddp(n,n), edge(2)
      dimension prn(*)
      common/io/inp, iout
      call drvply(q,wt,p,dp,ddp,edge,qdtyp,'none',0,n-1,nq,prn)
      if(prn(1)) then
         title='points'
         call prntrm(title,q,n,1,n,1,iout)
         title='weights'
         call prntrm(title,wt,n,1,n,1,iout)
      endif
      if(prn(2)) then
         title='polynomials'                  
         call prntrm(title,p,n,n,n,n,iout)  
         title='first derivative polynomials'
         call prntrm(title,dp,n,n,n,n,iout)  
         title='second derivative polynomials'        
         call prntrm(title,ddp,n,n,n,n,iout) 
      endif
      return
      end       


