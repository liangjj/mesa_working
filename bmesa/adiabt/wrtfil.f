
*deck wrtfil.f
c***begin prologue     wrtfil
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***description          
c***references       
c
c***routines called
c***end prologue       wrtfil
      subroutine wrtfil(eig,n,nr)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 eig, rad, dum
      character*80 title
      character*16 rowlab
      character*8 collab
      character*4 itoc
      dimension eig(nr,n+1), rowlab(nr), collab(n+1)
      write(90,*) 'Hyperradius'
      write(90,1) (eig(i,1),i=1,nr)
      do 10 i=2,n+1
         write(90,2) i
         write(90,1) (eig(j,i),j=1,nr)
 10   continue
      collab(1)='radius'
      do 30 i=2,n+1
         collab(i)='Psi-'//itoc(i-1)
 30   continue   
      call matpre(eig,nr,n+1,nr,n+1,0,1,rowlab,collab,0,dum,.false.)
      return
 1    format(5e15.8)
 2    format(/,'adiabatic eigenstate = ',i4)
      end
