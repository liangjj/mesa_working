*deck dagfil.f
c***begin prologue     dagfil
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           dagfililtonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            fill the diagonal of a matrix
c***                   
c***description                           
c***                   
c***references         
c
c***routines called    
c***end prologue       dagfil
      subroutine dagfil(fulmat,trimat,vec,n,type)
      implicit integer (a-z)
      real*8 fulmat, trimat, vec
      character*(*) type
      dimension fulmat(n,n), trimat(*), vec(*)
      common/io/inp, iout 
      if(type.eq.'full') then
         do 10 i=1,n
            fulmat(i,i)=vec(i)
 10      continue   
      else if(type.eq.'triangle') then
         ii=1
         trimat(1)=vec(1)
         do 20 i=2,n
            ii=ii+i
            trimat(ii)=vec(i)
 20      continue
      else
         call lnkerr('quit')
      endif   
      return
      end       


















