*deck @(#)rminv.f	5.1  11/6/94
      subroutine rminv(m,n,a,nb,b,det,dir)
c***begin prologue     rminv
c***date written       901991   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m935, link 935, linear systems, inversion
c***author             unknown
c***source             m935
c***purpose            driver for linpack rountines
c***description        calculates the solution of linear systems
c***                   or inverse od determinantst
c*** 
c
c
c      on input:               
c                      a: the matrix whose inverse is to be found.
c                      m: the 'row dimension' of the matrix a.
c                      n: the order of the matrix.
c                      nb:the number of columns in matrix b.
c                      b: the rhs of the linear system of equations
c
c      on output:
c                      a: the inverse of the orginal matrix a.
c                      b: the solution to the linear equations

c***references       
c
c***routines called    sgeco(linpack), sgedi(linpack), 
c***                   sgesl(linpack)
c***
c***end prologue       rminv
      implicit integer (a-z)
      parameter (maxdim=1000)
      integer ipvt(maxdim)
      real*8  a(m,n), b(m,nb), det(2), work(maxdim)
      real*8 one, zero, rcond
      character *(*) dir
      parameter (zero=0.0d+00,one=1.0d+00)
      common /io/inp,iout
c
c
      if (n.gt.maxdim) then
          write(iout,*) 'rminv: maxdim exceeded',n,'vs.',maxdim
          call lnkerr('m935: stop in rminv')
      endif
c
      call sgeco(a,m,n,ipvt,rcond,work)
      if ((rcond+one).eq.one) then
           call lnkerr ('matrix singular in sgeco')
      endif
      if(dir.eq.'all'.or.dir.eq.'solve') then
         job=0
         do 10 i=1,nb
 10         call sgesl(a,m,n,ipvt,b(1,i),job)
      endif
      if (dir.ne.'solve'.or.dir.ne.'all') then
          if (dir.eq.'inverse') then
              job=1
          elseif (dir.eq.'determinant') then
              job=10
          elseif (dir.eq.'inverse and determinant') then
              job=11
          elseif (dir.eq.'all') then
              job=11
          else
              call lnkerr('error in rminv call for directive')
          endif
          call sgedi(a,m,n,ipvt,det,work,job)
      endif
c
      return
      end
