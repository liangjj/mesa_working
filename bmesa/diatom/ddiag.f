*deck ddiag.f
c***begin prologue     ddiag
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            direct in-core diagonaliztion of hamiltonian.
c***references         
c
c***routines called    
c***end prologue       ddiag
      subroutine ddiag(tr,tth,eigr,eigth,ham,eig,work,dum,lwr,upr,
     1                 nr,nth,n,pottyp,typ,prn)
      implicit integer (a-z)
      real*8 tr, tth, eigr, eigth, ham, eig, work, dum
      logical prn
      character*(*) typ, pottyp
      character*3 itoc, li, lj
      character*80 title
      dimension tr(nr,nr), tth(nth,nth), eigr(nr), eigth(nth)
      dimension eig(n), ham(n,n), work(*), dum(*)
      common/io/inp, iout 
      call rzero(ham,n*n)
c     matrix will be computed as if blocked by angle
      icnt=0
      do 10 i=lwr,upr
         li=itoc(i)
         call pakstr(li,ii)
         jcnt=0
         do 20 j=lwr,i
            lj=itoc(j)
            call pakstr(lj,jj)
            call filham(ham(icnt+1,jcnt+1),tr,tth(i,j),eigr,i,j,nr,n)
            if (prn) then
                title='h('//li(1:ii)//','//lj(1:jj)//') block of '//
     1                'matrix'
                call prntrm(title,ham(icnt+1,jcnt+1),nr,nr,n,n,iout)
            endif
            jcnt=jcnt+nr 
 20      continue
         icnt=icnt+nr
 10   continue
      call potntl(ham,eigr,eigth,pottyp,lwr,upr,nr,nth,n)
      do 30 i=1,n
         do 40 j=1,i
            ham(j,i)=ham(i,j)
 40      continue
 30   continue   
      if(typ.eq.'eigenvalues-and-eigenvectors') then
         write(iout,1) n
         call tred2(n,n,ham,eig,work,ham)
         call tql2(n,n,eig,work,ham,ierr)
         title='eigenvalues of hamiltonian'
         call prntrm(title,eig,n,1,n,1,iout)
      elseif(typ.eq.'eigenvalues') then
         write(iout,2) n
         call tred1(n,n,ham,eig,work,dum)
         call tql1(n,eig,work,ierr)
         write(iout,*) 'ierr = ',ierr
         title='eigenvalues of hamiltonian'
         call prntrm(title,eig,n,1,n,1,iout)
      endif                  
      return
 1    format(/,5x,'calculating eigenvalues and eigenvectors for matrix s
     1ize = ',i8)
 2    format(/,5x,'calculating eigenvalues for matrix size = ',i8)
      end       
