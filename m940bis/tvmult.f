*deck tvmult.f
c***begin prologue     tvmult
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            test.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       tvmult
      subroutine tvmult(hbuf,ibuf,diag,vecs,temp,rep,fzcore,
     1                  lenbuf,n,nvec,title)
      implicit integer (a-z)
      real*8 hbuf, diag, vecs, temp, energy, sdot, rep, fzcore
      character*(*) title
      character*3 itoc
      logical nodsk
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n), vecs(n,nvec)
      dimension temp(n,nvec), title(3)
      common/io/inp, iout
      data nodsk/.false./
      do 10 i=1,nvec
         call iosys('read real "ci root '//itoc(i)
     1                          //'" from rwf',n,vecs(1,i),0,' ') 
 10   continue
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
      if(n.gt.0) then
         call iosys('read real '//title(3)//' from hamiltonian',
     1               n,diag,0,' ')
      endif
      call iosys('read integer '//title(2)//' from hamiltonian',
     1            1,nonzro,0,' ')
      if(nonzro.eq.0) then
         return
      endif
      call hmult(ibuf,hbuf,lenbuf,title(1),nonzro,vecs,temp,n,
     1           nvec,nodsk)
      do 20 j=1,n
         do 30 i=1,nvec
            temp(j,i) = temp(j,i) + diag(j)*vecs(j,i)
 30      continue
 20   continue
      do 40 i=1,nvec
         energy = sdot(n,vecs(1,i),1,temp(1,i),1) + rep + fzcore
	 write(iout,1) i, energy
 40   continue	    
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
      return
 1    format(/,1x,'root = ',i3,2x,'energy = ',e15.8)      
      end       


















