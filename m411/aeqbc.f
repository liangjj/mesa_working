*deck @(#)aeqbc.f	5.1  11/6/94
c***begin prologue     aeqbc
c***date written       890602   (yymmdd)
c***revision date               (yymmdd)
c***keywords           matrix muliply, arbitrary spacing
c***author             schneider, barry (lanl)
c***source             mylib
c***purpose            multiply two submatrices of a large matrix
c***                   with arbitrary spacing between the elements.
c***
c***description        a = bc
c***                   where a, b, and c are matrices or the transposes
c***                   of matrices. if the matrices are stored in the
c***                   normal fashion then it is sufficient to interchange
c***                   the row and column spacing parameters below to
c***                   get the transpose. 
c***references       
c
c***routines called    mxma (scilib)
c***end prologue      
      subroutine aeqbc (a,spcola,sprowa,b,spcolb,sprowb,
     1                  c,spcolc,sprowc,nb,nc,na)
      implicit integer (a-z)
      real*8 a, b, c
      dimension a(nb,na), b(nb,nc), c(nc,na)
      call mxma(b,spcolb,sprowb,c,spcolc,sprowc,a,spcola,sprowa,
     1          nb,nc,na)
      return
      end
