c \documentclass{article}
c \ublkackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{VTRBLK: Compute Davidson Trial Vectors}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck vtrblk.f
c***begin prologue     vtrblk
c***date written       010828   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            guess vectors based on blocked hamiltonian
c***                   for davidson routine.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vtrblk
      subroutine vtrblk(eig,eigtr,vec,ind,n,ntrial,prnt)
      implicit integer (a-z)
      logical prnt
      real*8 eig, eigtr, vec
      character*24 str
      character*2 itoc
      character*8 code
      character*80 title
      dimension eig(n), eigtr(ntrial), vec(n,ntrial)
      dimension ind(n)
      common/io/inp, iout
      data code / 'trial:' /
      call iosys('rewind vectors on ham read-and-write',0,0,0,' ')
      do 10 i=1,n
         ind(i)=i
 10   continue   
      do 20 ii=2,n
         i=ii-1
         k=i
         tmp=eig(i)
         i1=ind(i)
         do 30 j=ii,n
            if(eig(j).lt.tmp) then
               k=j
               tmp=eig(j)
            endif   
 30      continue
         if(k.ne.i) then
            ind(i)=ind(k)
            ind(k)=i1
            eig(k) = eig(i)
            eig(i) = tmp
         endif
         call copy(eig,eigtr,ntrial)
         do 40 trl=1,ntrial
            itr=ind(trl)
            offset=(itr-1)*n
            call iosys('read real vector from ham',n,vec(1,trl),
     1                  offset,' ')
 40      continue
 20   continue      
      if(prnt) then
         title='guess eigenvalues'
         call prntfm(title,eig,nroots,1,nroots,1,iout)
      endif   
      do 50 i=1,nroots
         call iosys('write real "'//code(1:6)//itoc(i)//'" to rwf',
     1               n,vec(1,i),0,' ')
 50   continue   
      if(prnt) then
         title='guess eigenvectors'
         call prntfm(title,vec,n,nroots,n,nroots,iout)
      endif            
      return
      end       

















