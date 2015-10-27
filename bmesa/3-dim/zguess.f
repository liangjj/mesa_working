*deck zguess.f
c***begin prologue     zguess
c***date written       960615   (yymmdd)
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
c***end prologue       zguess
      subroutine zguess(eig1,eig2,eig3,t1,t2,t3,ttmp,vec,hold,root,
     1                  ind,scr,nply,dim,matdim,nroots,pottyp,prnt)
      implicit integer (a-z)
      real*8 t1, t2, t3, eig1, eig2, eig3
      real*8 ttmp
      real*8 hold, vec, root, scr, tmp
      character*(*) pottyp
      character*80 title
      logical prnt
      dimension eig1(nply), eig2(nply), eig3(nply)
      dimension t1(nply,nply), t2(nply,nply), t3(nply,nply)
      dimension vec(matdim,nroots), root(nroots), scr(nply,4)
      dimension ind(matdim,*), hold(*), ttmp(nply,nply,3)
      common/io/inp, iout
      call copy(t1,ttmp(1,1,1),nply*nply)
      call copy(t2,ttmp(1,1,2),nply*nply)
      call tred2(nply,nply,ttmp(1,1,1),scr(1,1),scr(1,4),ttmp(1,1,1)) 
      call tql2(nply,nply,scr(1,1),scr(1,4),ttmp(1,1,1),ierr)
      call tred2(nply,nply,ttmp(1,1,2),scr(1,2),scr(1,4),ttmp(1,1,2)) 
      call tql2(nply,nply,scr(1,2),scr(1,4),ttmp(1,1,2),ierr)
      if(dim.eq.3) then
         call copy(t3,ttmp(1,1,3),nply*nply)
         call tred2(nply,nply,ttmp(1,1,3),scr(1,3),scr(1,4),ttmp(1,1,3))
         call tql2(nply,nply,scr(1,3),scr(1,4),ttmp(1,1,3),ierr)
      endif
      if(dim.eq.2) then           
         count=0
         do 10 i=1,nply
            do 20 j=1,nply
               count=count+1
               hold(count) = scr(i,1) + scr(j,2)
 20         continue   
 10      continue
         do 30 ii=2,matdim
            i=ii-1
            k=i
            tmp=hold(i)
            i1=ind(i,1)
            j1=ind(i,2)
            do 40 j=ii,matdim
               if(hold(j).lt.tmp) then
                  k=j
                  tmp=hold(j)
               endif   
 40         continue
            if(k.ne.i) then
               ind(i,1)=ind(k,1)
               ind(i,2)=ind(k,2)
               ind(k,1)=i1
               ind(k,2)=j1
               hold(k) = hold(i)
               hold(i) = tmp
            endif
 30      continue                 
         do 50 i=1,nroots
            root(i)=hold(i)
            count=0
            do 60 j=1,nply
               do 70 k=1,nply
                  count=count+1
                  vec(count,i)=ttmp(j,ind(i,1),1)*ttmp(k,ind(i,2),2)
 70            continue
 60         continue
 50      continue
      else
         count=0
         do 80 i=1,nply
            do 90 j=1,nply
               do 100 k=1,nply
                  count=count+1
                  hold(count) = scr(i,1) + scr(j,2) + scr(k,3)
 100           continue
 90         continue   
 80      continue
         do 200 ii=2,matdim
            i=ii-1
            k=i
            tmp=hold(i)
            i1=ind(i,1)
            j1=ind(i,2)
            k1=ind(i,3)
            do 210 j=ii,matdim
               if(hold(j).lt.tmp) then
                  k=j
                  tmp=hold(j)
               endif   
 210        continue
            if(k.ne.i) then
               ind(i,1)=ind(k,1)
               ind(i,2)=ind(k,2)
               ind(i,3)=ind(k,3)
               ind(k,1)=i1
               ind(k,2)=j1
               ind(k,3)=k1
               hold(k) = hold(i)
               hold(i) = tmp
            endif
 200     continue   
         do 300 i=1,nroots
            root(i)=hold(i)
            count=0
            do 310 j=1,nply
               do 320 k=1,nply
                  do 330 l=1,nply
                     count=count+1
                     vec(count,i) = ttmp(j,ind(i,1),1) * 
     1                                               ttmp(k,ind(i,2),2)
     2                                                 *
     3                                               ttmp(l,ind(i,3),3)
 330              continue     
 320           continue
 310        continue
 300     continue   
      endif
      if(prnt) then
         title='guess eigenvalues'
         call prntrm(title,hold,nroots,1,nroots,1,iout)
         title='guess eigenvectors'
         call prntrm(title,vec,matdim,nroots,matdim,nroots,iout)
      endif                     
      return
      end       






