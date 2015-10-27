*deck apdr1r2.f
c***begin prologue     apdr1r2
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           numerical derivatives wrt (r1,r2) from
c**                    DVR expansion coefficients
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***description          
c***references       
c
c***routines called
c***end prologue       apdr1r2
      subroutine apdr1r2(dfdr1,dfdr2,c,f1,f2,df1,df2,n1,n2)
      implicit integer (a-z)
      real*8 dfdr1, dfdr2, c, f1, f2, df1, df2
      dimension dfdr1(n2,n1), dfdr2(n2,n1), c(n2,n1)
      dimension f1(n1,n1), f2(n2,n2), df1(n1,n1), df2(n2,n2)
      character*80 title
      logical prnt, srf
      common/io/inp,iout
      call ebct(dfdr1,c,df1,n2,n1,n1)   
      do 10 i=1,n1
         do 20 j=1,n2
            dfdr1(j,i)=f2(j,j)*dfdr1(j,i)
 20      continue
 10   continue
      call ebc(dfdr2,df2,c,n2,n2,n1)
      do 30 i=1,n1
         do 40 j=1,n2
            dfdr2(j,i)=f1(i,i)*dfdr2(j,i)
 40      continue
 30   continue
      return
      end











