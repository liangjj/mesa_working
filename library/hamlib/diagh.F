*deck diagh.f
c***begin prologue     diagh
c***date written       000714   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           diagonalize
c***author             schneider, barry (nsf)
c***source             
c***purpose            diagonalize packed hamiltonian
c***                   
c***description          
c***references         
c
c***routines called    
c***end prologue       diagh
      subroutine diagh(ham,eig,work,hbuf,buf,diag,len,
     1                 nonzro,incore,drctv,n)
      implicit integer (a-z)
      real*8 ham, eig, work, hbuf, diag
      character*80 title
      character*(*) drctv
      logical incore
      dimension ham(n,n), eig(n), work(n), hbuf(len), buf(len,2)
      dimension diag(n)
      common/io/inp, iout      
      if(incore) then
         call rzero(ham,n*n)
         do 10 nel=1,nonzro
            i=buf(nel,1)
            j=buf(nel,2)
            ham(i,j)=hbuf(nel)
            ham(j,i)=ham(i,j)
 10      continue   
      else
         call iosys('read integer "number of non-zero '//
     1              'matrix elements" from ham',1,nonzro,0,' ')
         trips=nonzro/len
         left=nonzro-trips*len
         do 20 trp=1,trips
            call iosys('read integer "hamiltonian buffers" from ham '//
     1                 'without rewinding',2*len,buf,0,' ')
            call iosys('read integer "hamiltonian buffers" from ham '//
     1                 'without rewinding',wptoin(len),hbuf,0,' ')
            do 30 nel=1,len
               i=buf(nel,1)
               j=buf(nel,2)
               ham(i,j)=hbuf(nel)
               ham(j,i)=ham(i,j)
 30         continue   
 20      continue   
         if(left.ne.0) then
            call iosys('read integer "hamiltonian buffers" from ham '//
     1                 'without rewinding',2*left,buf,0,' ')
            call iosys('read integer "hamiltonian buffers" from ham '//
     1                 'without rewinding',wptoin(left),hbuf,0,' ')
            do 40 nel=1,left
               i=buf(nel,1)
               j=buf(nel,2)
               ham(i,j)=hbuf(nel)
               ham(j,i)=ham(i,j)
 40         continue   
         endif
      endif
      do 50 i=1,n
         ham(i,i)=diag(i)
 50   continue   
      if(drctv.eq.'diagonalize') then
         call dsyev('v','l',n,ham,n,eig,work,5*n,info)
         title='eigenvalues'
         call prntrm(title,eig,n,1,n,1,iout)
      endif
c      title='eigenvectors'
c      call prntrm(title,ham,n,n,n,n,iout)   
      return
      end       
