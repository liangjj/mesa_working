*deck diff.f
c***begin prologue     diff
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           orthogonal polynomials
c***author             schneider, barry (nsf)
c***source             math
c***purpose            solve differential equation by expanding in
c***                   orthogonal polynomials.
c***                   
c***references         
c
c***routines called    
c***end prologue       diff
      subroutine diff(p,dp,ddp,x,wts,ham,eig,der,l,n,npts,pottyp,
     1                prange,coord,prnt)
      implicit integer (a-z)
      real*8 p, dp, ddp, x, ham, eig, wts, der, prange
      character*80 title
      character*(*) pottyp
      logical prnt, coord, typ
      dimension p(npts,0:n-1), dp(npts,0:n-1), ddp(npts,0:n-1)
      dimension x(npts), wts(npts), ham(n,n), eig(n), der(2)
      common/io/inp, iout 
      if(coord) then
         write(iout,1)
      endif
      call rzero(ham,n*n)
      do 10 i=1,n
         do 20 j=1,i
            do 30 k=2,npts
               ham(i,j) = ham(i,j) - .5d0*wts(k)*p(k,i-1)*ddp(k,j-1)
 30         continue
            ham(i,j) = ham(i,j) +.5d0 *
     1            ( p(npts,i-1)*( dp(npts,j-1) -der(2)*p(npts,j-1) ) 
     2                          - 
     3              p(1,i-1)*( dp(1,j-1) - der(1)*p(1,j-1) ) )
            ham(j,i)=ham(i,j)
 20      continue   
 10   continue
      if(coord) then
         do 40 i=1,n
            if(eig(i).le.prange) then
               ham(i,i) = ham(i,i) + .5d0*l*(l+1)/(eig(i)*eig(i))
            endif   
 40      continue
      else
         do 50 i=1,n
            do 60 j=1,i
               do 70 k=1,npts
                  if(x(k).le.prange) then
                     ham(i,j) = ham(i,j) + .5d0*wts(k)*p(k,i-1)*
     1                                      l*(l+1)*p(k,j-1)/
     2                                                  ( x(k)*x(k) )
                  endif
 70            continue
 60         continue
 50      continue
      endif        
      if (pottyp.eq.'exponential') then
          if(coord) then
             do 80 i=1,n
                if(eig(i).le.prange) then
                   ham(i,i) = ham(i,i) -exp(-eig(i))
                endif   
 80          continue
          else   
             do 90 i=1,n
                do 100 j=1,i
                   do 110 k=1,npts
                      if(x(k).le.prange) then
                         ham(i,j) = ham(i,j) - 
     1                              p(k,i-1)*wts(k)*exp(-x(k))*p(k,j-1)
                      endif
 110               continue
                   ham(j,i) = ham(i,j)
 100            continue
 90          continue   
          endif 
      elseif(pottyp.eq.'one') then
          if(coord) then
             do 200 i=1,n
                if(eig(i).le.prange) then
                   ham(i,i) = ham(i,i) - 1.d0
                endif   
 200         continue
          else   
             do 210 i=1,n
                do 220 j=1,i
                   do 230 k=1,npts
                      if(x(k).le.prange) then
                         ham(i,j) = ham(i,j) - p(k,i-1)*wts(k)*p(k,j-1)
                      endif   
 230               continue
                   ham(j,i) = ham(i,j)
 220            continue
 210          continue   
          endif
      elseif(pottyp.eq.'coulomb') then
          if(coord) then
             do 300 i=1,n
                if(eig(i).le.prange) then 
                   ham(i,i) = ham(i,i) - 1.d0/eig(i)
                endif   
 300         continue
          else   
             do 310 i=1,n
                do 320 j=1,i
                   do 330 k=2,npts
                      if(x(k).le.prange) then
                          ham(i,j) = ham(i,j) - 
     1                               p(k,i-1)*wts(k)*p(k,j-1)/x(k)
                      endif   
 330               continue
                   ham(j,i) = ham(i,j)
 320            continue
 310         continue   
          endif
      endif
      if (prnt) then
          title='hamiltonian'   
          call prntrm(title,ham,n,n,n,n,iout)
      endif
      return
 1    format(/,5x,'using co-ordinate representation of hamiltonian')
      end       

