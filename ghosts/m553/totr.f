      subroutine totr(triang,square,num,nnp)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
c
      implicit integer (a-z)
c
      real*8 triang(nnp),square(num,num)
c
      ij=0
      do 2 i=1,num
         do 1 j=1,i
            ij=ij+1
            triang(ij)=square(j,i)
    1    continue
    2 continue
c
      return
      end
