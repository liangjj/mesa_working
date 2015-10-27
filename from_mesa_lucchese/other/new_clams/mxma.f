*deck %W%  %G%
      subroutine mxma(a,na,iad,b,nb,ibd,c,nc,icd,nar,nac,nbc)
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
