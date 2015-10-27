      subroutine rho(nstate,labels,allrho,nbf,maxp,maxnbf,d)
c
c reads density matrices in ao basis from mesa m950
c they come in full square and are indexed precisely the same
c way that the contracted basis functions are indexed on geobas
c which was read by subroutine inputs in this code.  The only check
c that can be done on the indexing at this stage is to test that
c the number of ao's (contracted gaussians) is the same as the dimension
c of the density matrices.
c
c
      implicit real*8 (a-h,o-z)
c      dimension irow(maxnbf)
      dimension allrho((maxnbf*(maxnbf+1))/2,maxp),labels(maxp)
      dimension d(maxnbf,maxnbf)
      mx = 100
c
c  read number of channels and check for # of basis functions
c
      read(12) nchan, nbfden
      write(66,*)" nchan, nbfden",nchan, nbfden
      if(nbf.ne.nbfden) then
       write(66,202)
  202  format(' Error size of density matrices .ne. #basis functions')
       stop
      endif
c
c  open loop on channel pairs
c      we read a density matrix and load it after triangularization in allrho
c
      do 50 i=1,nchan
      do 50 j=1,i
      read(12) ichan, jchan
      write(66,*)" ichan, jchan",ichan, jchan
      if(i.ne.ichan.or.j.ne.jchan) then
       write(66,201)
  201  format(' Error channel indices do not match in rho')
       stop
      endif
      write(77,"(a10,i5,a10,i5)")"iroot= ", i, "jroot= ",j
      do 9 l=1,nbf
      read(12) (d(m,l),m=1,nbf)
 9    write(77,777)l,(d(m,l),m=1,nbf)
 777  format(i3,8f12.6/(3x,8f12.6))
c
c  put this density matrix into allrho
c
      istate = i*(i-1)/2 + j
      labels(istate) = 1000*ichan + jchan
      do 15 ia=1,nbf
      iam1 = ia-1
      if(iam1.ne.0) then
       do 14 ib=1,iam1
       indx = ia*(ia-1)/2 + ib
  14   allrho(indx,istate) = d(ia,ib)+d(ib,ia)
      endif
      indx = ia*(ia+1)/2
      allrho(indx,istate) = d(ia,ia)
  15  continue
c
c  close loop on pairs of states
c
c*****
c      write(66,9991) i,j,istate
c 9991 format(' i,j,istate',3i5)
c      do 9992 junk=1,nbf
c         write(66,9993) junk,(d(kunk,junk),kunk=1,nbf)
c 9993 format(1x,i3,5f12.5,/,(4x,5f12.5))
c 9992  continue
  50  continue
      nstate=nchan*(nchan+1)/2
      return
      end
