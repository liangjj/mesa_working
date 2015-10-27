*deck fromin.f
c***begin prologue     fromin
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           input matrix
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       fromin
      subroutine fromin(ham,cham,work,n,mattyp)
      implicit integer (a-z)
      complex*16 cham
      real*8 ham, work
      character*80 title
      character*1600 card
      character*1 itoc
      logical dollar
      character*(*) mattyp
      dimension ham(n,*), cham(n,*), work(n)
      common/io/inp, iout
      if ( dollar('$input',card,title,inp) ) then
          write(iout,1)
      endif                
      if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
         do 10 i=1,n
            call fparr(card,'column-'//itoc(i)//'-of-matrix',
     1                 work,n,' ')
            do 20 j=1,n
               cham(i,j)=work(j)  
 20         continue
 10      continue
         title='input matrix'
         call prntcm(title,cham,n,n,n,n,iout)
      else
         do 30 i=1,n
            call fparr(card,'column-'//itoc(i)//'-of-matrix',
     1                 ham(1,i),n,' ') 
 30      continue
         title='input matrix'
         call prntrm(title,ham,n,n,n,n,iout)
      endif
      return
 1    format(/,5x,'inputting matrix by hand')      
      end       






