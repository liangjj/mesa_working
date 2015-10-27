*deck condit.f
c***begin prologue     condit
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           condition, matrix
c***author             schneider, barry(nsf)
c***source             
c***purpose            calculates condition of tridiagonal.
c***
c***references
c
c***routines called
c***end prologue      condition
      subroutine condit(e,rcond,band,r,v,ipiv,f,stp,n,order,bw,bc,type)
      implicit integer (a-z)
      dimension band(n,-bw:bw), r(n), v(n), f(n,3)
      dimension det(2), ipiv(n,2)
      real*8 e, band, r, v, f, stp, det, deter
      real*8 slangt, anorm, rcond
      character*(*) type, bc
      common /io/ inp, iout
      type='standard-numerov'
      do 10 i=1,n
         f(i,1)=v(i)-e
 10   continue     
      if (type.eq.'standard-finite-difference'
     1                   .or.
     2    type.eq.'standard-numerov') then
          call fdiff(band,f(1,1),stp,n,order,bw,type)
          call bndcnd(band,n,order,bw,bc)
      else
          call gfdiff(band,r,f(1,1),n,order,bw,type)
      endif
      info=10   
      anorm=slangt('1',n-2,band(3,-1),band(2,0),band(2,1))
      call sgttrf(n-2,band(3,-1),band(2,0),band(2,1),f(1,1),
     1            ipiv,det,info)
      call sgtcon('1',n-2,band(3,-1),band(2,0),band(2,1),f,
     1            ipiv(1,1),anorm,rcond,f(1,2),ipiv(1,2),info)
      deter=det(1)*10.d0**det(2)
      write(iout,1) e, deter, rcond
      rcond=deter  
      return
 1    format(/,'energy = ',e15.8,2x,'determinant = ',e15.8,
     1       /,'condition parameter = ',e15.8)
      end

