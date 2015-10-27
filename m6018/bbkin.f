*deck @(#)bbkin.f	1.1 9/8/91
c***begin prologue     bbkin
c***date written       930419   (yymmdd)
c***revision date               (yymmdd)
c***keywords           kohn integrals
c***author             schneider barry (lanl)
c***source             m6018
c***purpose            read bound-bound kinetic energy integrals and
c***                   form channel bound-bound kinetic energy matrix
c***
c***                   nbtot  = total no. bound orbitals associated with
c***                            a channel. 
c***                   orblst = list of bound orbitals with the variational
c***                            orbitals first and then the target orbitals. 
c***                   kein   = input kinetic energy matrix.
c***                   keout  = output matrix containing transformation matrix
c***                            matrix.
c***                   eig    = eigenvalues of kinetic energy operator. 
c***                   
c***references       
c
c***routines called
c***end prologue       bbkin  
      subroutine bbkin(kein,keout,eig,dum,echan,nchan,nmo,nbtot,
     1                 orblst,dimmo,dimc,prnt)
      implicit integer (a-z)
      real*8 kein, keout, echan, eig, dum
      character*80 title
      character*3 itoc
      logical prnt
      common /io/ inp, iout
      dimension kein(nmo,nmo), keout(*), nbtot(dimc)
      dimension echan(dimc), orblst(dimmo,dimc), eig(*), dum(*)
      call iosys ('read real "mo kinetic energy integrals" from kohndt',
     1             nmo*nmo,kein,0,' ')
      loch=1
      do 10 ch=1,nchan
         call filmat(kein,keout(loch),nmo,nmo,nbtot(ch),nbtot(ch),
     1               orblst(1,ch))
c-----------------------------------------------------------------------c
c             diagonalize kinetic energy operator for this channel      c
c-----------------------------------------------------------------------c
         call tred2(nbtot(ch),nbtot(ch),keout(loch),eig,dum,
     1              keout(loch))
         call tql2(nbtot(ch),nbtot(ch),eig,dum,keout(loch),ierr)
         title='"bound kinetic energy eigenvalues for channel-'//       
     1           itoc(ch)//'"'
         call iosys ('write real '//title//' to kohnint',nbtot(ch),eig,
     1                0,' ')
         write (iout,1) title
         write (iout,2) (eig(ii), ii=1,nbtot(ch))
         title='"bound-bound kinetic energy transformation matrix '//
     1         'for channel-'//itoc(ch)//'"'
         call iosys ('write real '//title//' to kohnint',
     1                 nbtot(ch)*nbtot(ch),keout(loch),0,' ')               
         if (prnt) then
             call prntrm (title,keout(loch),nbtot(ch),nbtot(ch),
     1                    nbtot(ch),nbtot(ch),iout)  
         endif
         loch=loch+nbtot(ch)*nbtot(ch)
   10 continue    
      return
    1 format (a80)
    2 format ((/,5x,5e15.8))
      end


