*deck exdr1r2.f
c***begin prologue     exdr1r2
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           chain rule 
c***author             schneider, barry (nsf)
c***source             
c***purpose            chain rule for df/dr1 and df/dr2 
c***                   from df/dr and df/dang
c***description          
c***references       
c
c***routines called
c***end prologue       exdr1r2
      subroutine exdr1r2(dfdr1,dfdr2,dfdang,dfdr,rho,ang,
     1                   type,n1,n2)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 dfdr1, dfdr2, dfdr, dfdang, rho, ang, sn, cn, xinv
      character*(*) type
      dimension dfdr1(n2,n1), dfdr2(n2,n1), dfdr(n2,n1), dfdang(n2,n1)
      dimension rho(n2,n1), ang(n2,n1)
      if(type.eq.'Angular') then
         do 10 i=1,n1
            do 20 j=1,n2
               sn = sin(ang(j,i)) 
               cn = cos(ang(j,i))
               xinv = 1.d0/rho(j,i)
               dfdr1(j,i) = dfdang(j,i)*cn*xinv
               dfdr2(j,i) = - dfdang(j,i)*sn*xinv
 20         continue   
 10      continue
      elseif(type.eq.'Radial') then
         do 30 i=1,n1
            do 40 j=1,n2
               sn = sin(ang(j,i)) 
               cn = cos(ang(j,i))
               dfdr1(j,i) = dfdr(j,i)*sn 
               dfdr2(j,i) = dfdr(j,i)*cn
 40         continue   
 30      continue
      elseif(type.eq.'General') then
         do 50 i=1,n1
            do 60 j=1,n2
               sn = sin(ang(j,i)) 
               cn = cos(ang(j,i))
               xinv = 1.d0/rho(j,i)
               dfdr1(j,i) = dfdr(j,i)*sn + dfdang(j,i)*cn*xinv
               dfdr2(j,i) = dfdr(j,i)*cn - dfdang(j,i)*sn*xinv
 60         continue   
 50      continue   
      else
         call lnkerr('error in type')
      endif
      return
      end

