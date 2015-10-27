*deck rprop1
c***begin prologue     rprop1
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           sector, r-matrix
c***author             schneider, barry (nsf)
c***source             
c***purpose            construct global r-matrix
c***description        global r-matrix is constructed from previous
c***                   global r-matrix and new sector r-matrices.
c***references       
c
c***routines called
c***end prologue       rprop1
      subroutine rprop1 (r4,dchan,tot,ipvt,nstep,n,space)
      implicit real *8 (a-h,o-z)
      common /io/ inp, iout
      dimension r4(n,n), tot(n,n), dchan(n,4), space(n,n,3)
      dimension ipvt(n)
c
c     space(1,1,1) holds z inverse
c     space(1,1,2) holds sector r(2)
c     space(1,1,3) holds sector r(3)
      call vinv(dchan(1,3),dchan(1,3),n)
      call vmul(dchan(1,4),dchan(1,4),dchan(1,3),n)
c
c        first step
c
      if (nstep.eq.1) return
c
c     rmatrix propagation after the first step
c
c     build sector r(2) and start building sector r(1)
c
      call mvmul(dchan(1,3),tot,space(1,1,2),n,n)
      call mvmul(dchan(1,4),tot,space(1,1,3),n,n)
c
c     finish sector r(1) and build matrix for inversion (z inverse)
c
      call apbct(r4,space(1,1,3),tot,n,n,n)
      call copy(r4,space(1,1,1),n*n)
c
c     save sector r(3) and invert z inverse
c
      call transp (space(1,1,2),space(1,1,3),n,n)
      call sgefa (space(1,1,1),n,n,ipvt,info)
      do 10 i=1,n
         call sgesl (space(1,1,1),n,n,ipvt,space(1,i,2),0)
   10 continue     
c
c     z now computed from z inverse. z*r(2) also computed in space(1,1,2)
c
c     finish recursion by left mulitplication with r(3) and add old r4
c
      call rzero(r4,n*n)
      call ambc(r4,space(1,1,3),space(1,1,2),n,n,n)
      do 20 i=1,n
         r4(i,i)=dchan(i,4)+r4(i,i)
   20 continue
      return
      end
