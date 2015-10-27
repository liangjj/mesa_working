*deck potntl.f
c***begin prologue     potntl
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           potential
c***author             schneider, barry (nsf)
c***source             mesa
c***purpose            calculate potential
c***                   
c***references         
c
c***routines called    
c***end prologue       potntl
      subroutine potntl(ham,eigr,eigth,pottyp,lwr,upr,nr,nth,n)
      implicit integer (a-z)
      real*8 ham, eigr, eigth, denom, den1, den2
      character*(*) pottyp
      dimension ham(n,n), eigr(nr), eigth(nth)
      common/io/inp, iout
      count=0 
      if (pottyp.eq.'exponential') then
          do 10 i=lwr,upr
             do 20 j=1,nr
                count=count+1                
                ham(count,count) = ham(count,count) - exp(-eigr(j))
 20          continue
 10       continue
      elseif(pottyp.eq.'one') then
          do 30 i=lwr,upr
             do 40 j=1,nr
                count=count+1
                ham(count,count) = ham(count,count) - 1.d0
 40          continue
 30       continue
      elseif(pottyp.eq.'coulomb') then
          do 50 i=lwr,upr
             do 60 j=1,nr
                count=count+1
c                denom=sqrt(eigr(j)*eigr(j)-2.d0*eigr(j)*eigth(i)+1.d0)
c                ham(count,count) = ham(count,count) - 1.d0/denom
                 ham(count,count) = ham(count,count) - 1.d0/eigr(j)
 60          continue
 50       continue 
      elseif(pottyp.eq.'h2-plus') then
          do 70 i=lwr,upr
             do 80 j=1,nr
                count=count+1
                den1=sqrt(eigr(j)*eigr(j)-2.d0*eigr(j)*eigth(i)+1.d0)
                den2=sqrt(eigr(j)*eigr(j)+2.d0*eigr(j)*eigth(i)+1.d0)
                ham(count,count) = ham(count,count) - 1.d0/den1
     1                                              -1.d0/den2
 80          continue
 70       continue 
      endif
      return
 1    format(/,5x,'using co-ordinate representation of hamiltonian')
      end       

