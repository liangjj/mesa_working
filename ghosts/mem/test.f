*deck %W% %G%
      program test
      implicit integer(a-z)
c
c
      real*8 z
      integer a
      parameter (maxsiz=10000)
      common//a(maxsiz)
      equivalence(a,z)
c
c
      call inimem(a,maxsiz)
      i1=getpt('real',10,'i1',.false.)
      i2=getpt('real',100,'i2',.false.)
      i3=getpt('integer',10,'i3',.false.)
      i4=getpt('integer',100,'i4',.false.)
      i5=getpt('integer',200,'i5',.false.)
      i2=relpt('i2')
      i4=relpt('i4')  
      i3=relpt('i3')
      i2=getpt('integer',200,'i2',.false.)
      i1=relpt('i1')
      i7=getpt('integer',3,'i7',.false.)
      i8=getpt('real',3,'i8',.false.)
      i9=getpt('real',3,'i9',.true.)
      i9=relpt('i9')
c
c
      stop
      end
      
