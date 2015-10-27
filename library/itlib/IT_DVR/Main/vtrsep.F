c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{VTRSEP: Compute Davidson Trial Vectors}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck vtrsep.f
c***begin prologue     vtrsep
c***date written       010828   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            guess vectors based on separable hamiltonian
c***                   for davidson routine.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vtrsep
      subroutine vtrsep(u01,u02,u03,eig01,eig02,eig03,
     1                  vec,eig,ind,n1,n2,n3,dim,n,nroots,prnt)
      implicit integer (a-z)
      logical prnt
      real*8 u01, u02, u03, eig01, eig02, eig03, vec
      real*8 eig, tmp, temp
      character*80 title, cpass, chrkey
      character*240 card
      character*5 itoc
      character*8 code
      logical dollar
      dimension u01(n1,n1), u02(n2,n2), u03(n3,n3)
      dimension eig01(n1), eig02(n2), eig03(n3)
      dimension vec(n,nroots), eig(nroots), ind(n,*)
      common/io/inp, iout
      data code / 'trial:' /
c
      if( dollar('$trial',card,cpass,inp) ) then
          title=chrkey(card,'type-trials','unperturbed',' ')
      endif
c
      if(title.eq.'unit') then
         call rzero(vec,n*nroots)
	 do 10 i=1,nroots
	    vec(i,i)=1.d0
 10      continue	          	     
      elseif(title.eq.'unperturbed') then
c
c       the ind array is destroyed by this routine and needs to be copied or
c       regenerated if it is used again.
c        
         if(dim.eq.1) then
            call copy(eig01,eig,nroots)
            call copy(u01,vec,n*nroots)  
         elseif(dim.eq.2) then           
            count=0
            do 100 i=1,n1
               do 200 j=1,n2
                  count=count+1
                  eig(count) = eig01(i) + eig02(j)
 200           continue   
 100        continue
            do 1000 ii=2,n
               i=ii-1
               k=i
               tmp=eig(i)
               i1=ind(i,1)
               j1=ind(i,2)
               do 2000 j=ii,n
                  if(eig(j).lt.tmp) then
                     k=j
                     tmp=eig(j)
                  endif   
 2000          continue
               if(k.ne.i) then
                  ind(i,1)=ind(k,1)
                  ind(i,2)=ind(k,2)
                  ind(k,1)=i1
                  ind(k,2)=j1
                  eig(k) = eig(i)
                  eig(i) = tmp
                  endif
 1000       continue                 
            do 3000 i=1,nroots
               count=0
               do 4000 j=1,n1
                  do 5000 k=1,n2
                     count=count+1
                     vec(count,i)=u01(j,ind(i,1))*u02(k,ind(i,2))
 5000             continue
 4000          continue
 3000       continue
         else
            count=0
            do 15000 i=1,n1
                  do 16000 j=1,n2
                     do 17000 k=1,n3
                        count=count+1
                        eig(count) = eig01(i) + eig02(j) + eig03(k)
17000                continue
16000             continue   
15000          continue
            do 18000 ii=2,n
               i=ii-1
               k=i
               tmp=eig(i)
               i1=ind(i,1)
               j1=ind(i,2)
               k1=ind(i,3)
               do 19000 j=ii,n
                  if(eig(j).lt.tmp) then
                     k=j
                     tmp=eig(j)
                  endif   
19000          continue
               if(k.ne.i) then
                  ind(i,1)=ind(k,1)
                  ind(i,2)=ind(k,2)
                  ind(i,3)=ind(k,3)
                  ind(k,1)=i1
                  ind(k,2)=j1
                  ind(k,3)=k1
                  eig(k) = eig(i)
                  eig(i) = tmp
               endif
18000       continue   
            do 20000 i=1,nroots
               count=0
               do 21000 j=1,n1
                  do 22000 k=1,n2
                     do 23000 l=1,n3
                        count=count+1
                        vec(count,i) = u01(j,ind(i,1)) * 
     1                                 u02(k,ind(i,2)) * 
     2                                 u03(l,ind(i,3))
23000                continue     
22000             continue
21000          continue
20000       continue
         endif
	 if(prnt) then
            title='guess eigenvalues'
            call prntfm(title,eig,nroots,1,nroots,1,iout)
	 endif   
      else
         call lnkerr('error in trial type')
      endif
      do 20 i=1,nroots
         call iosys('write real "'//code(1:6)//itoc(i)//'" to rwf',
     1               n,vec(1,i),0,' ')
 20   continue   
      if(prnt) then
         title='guess eigenvectors'
         call prntfm(title,vec,n,nroots,n,nroots,iout)
      endif            
      return
      end       

















