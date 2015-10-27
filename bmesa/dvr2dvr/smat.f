*deck smat.f
c***begin prologue     smat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            overlap matrix of two DVR.
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       smat
      subroutine smat(pi,pj,over,nfini,nfouti,nfinj,nfoutj,npts,
     1                ntot,sets)
      implicit integer (a-z)
      real*8 pi, pj, over, thresh
      character*80 title
      character*(*) sets
      character*2 itoc
      dimension pi(npts,0:nfini-1), pj(npts,0:nfinj-1)
      dimension over(nfini,nfinj)
      common/io/inp, iout
      data thresh/1.d-08/
      call  schmq2p(pi,pj,thresh,npts,nfini,nfouti,nfinj,
     1              nfoutj,.true.,sets)
      write(iout,1) nfini, nfouti, nfinj, nfoutj
      if(sets.eq.'same') then
         ntot=ntot+nfoutj
      endif
      do 40 i=1,nfouti
         do 50 j=1,nfoutj
            over(i,j)=0.d0
            do 60 k=1,npts
               over(i,j)=over(i,j) +pi(k,i-1)*pj(k,j-1)
 60         continue
 50      continue
 40   continue                     
      title='overlap of DVR basis sets ni ='//itoc(nfouti)//
     1      ' nj ='//itoc(nfoutj)
      call prntrm(title,over,nfouti,nfoutj,nfini,nfinj,iout)
      return
 1    format(/,1x,'number of initial functions in set 1 = ',i4,/,1x,
     1            'number of final functions in set 1   = ',i4,/,1x,
     2            'number of initial functions in set 2 = ',i4,/,1x,
     3            'number of final functions in set 2   = ',i4)
      end       
