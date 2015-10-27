*deck @(#)mxma.f	5.1  11/6/94
      subroutine mxma(a,na,iad,b,nb,ibd,c,nc,icd,nar,nac,nbc)
c***begin prologue     mxma.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             lengsfield, byron iii (llnl)
c***source             @(#)mxma.f	5.1   11/6/94
c***purpose            matrix multiply 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       mxma.f
c
      implicit integer (a-z)
c
      real*8 a(*),b(*),c(*),t
c
      do 3 i=1,nar
         pc=(i-1)*nc+1
         do 2 j=1,nbc
            t=0.0d+00
            pa=(i-1)*na+1
            pb=(j-1)*ibd+1
            do 1 k=1,nac
               t=t+a(pa)*b(pb)
               pa=pa+iad
               pb=pb+nb
 1          continue
            c(pc)=t
            pc=pc+icd
 2       continue
 3    continue
c
c
      return
      end
