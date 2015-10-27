*deck filodg
c***begin prologue     filodg
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           hamiltonian, matrix
c***author             schneider, barry (nsf)
c***source             
c***purpose            fill diagonal and off diagonal channel blocks
c***                   of the hamiltonian matrix.
c***description        four types of channel blocks can appear
c***                   corresponding to the closed-closed, closed-open,
c***                   open-closed and open-open portions of the hamiltonian
c***                   matrix. the hamiltonian matrix is passed for all four
c***                   parts of this block in the calling routine.  the v
c***                   matrix for the block is passed in its entirety.  it is
c***                   assumed that the v matrix is blocked in the same fashion
c***                   with all closed orbitals first in the list.  this
c***                   matrix is not hermitian when the channel indices are
c***                   unequal.
c***                   
c
c***routines called
c***end prologue       filodg
      subroutine filodg(hcc12,hcc21,hco12,hco21,hoc12,hoc21,hoo12,
     1                  hoo21,v,ci,cj,nic,nio,njc,njo,ma,mb,n)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 hcc12, hcc21, hco12, hco21, hoc12, hoc21, hoo12, hoo21
      real*8 v
      character*80 title
      dimension hcc12(n,n), hcc21(n,n), hco12(n,n), hco21(n,n)
      dimension hoc12(n,n), hoc21(n,n), hoo12(n,n), hoo21(n,n)
      dimension v(ma,mb)
      nbegi=nic+1
      nfnali=nic+nio
      nbegj=njc+1
      nfnalj=njc+njo
      if(nfnali.ne.ma) then
         call lnkerr('error in indexing in filodg')
      endif
      if(nfnalj.ne.mb) then
         call lnkerr('error in indexing in filodg')
      endif
c      write(iout,*) 'channel = ',ci,' channel = ',cj
c      title='v'
c      call prntrm(title,v,ma,mb,ma,mb,iout)
c      title='hcc12'
c      call prntrm(title,hcc12,nic,njc,n,n,iout)
c      title='hcc21'
c      call prntrm(title,hcc21,njc,nic,n,n,iout)
c      title='hco12'
c      call prntrm(title,hco12,nic,njo,n,n,iout)
c      title='hco21'
c      call prntrm(title,hco21,njc,nio,n,n,iout)
c      title='hoc12'
c      call prntrm(title,hoc12,nio,njc,n,n,iout)
c      title='hoc21'
c      call prntrm(title,hoc21,njo,nic,n,n,iout)
c      title='hoo12'
c      call prntrm(title,hoo12,nio,njo,n,n,iout)
c      title='hoo21'
c      call prntrm(title,hoo21,njo,nio,n,n,iout)
      if (ci.eq.cj) then
          if (nic.ne.0) then
              do 10 i=1,nic
                 do 20 j=1,i
                    hcc12(i,j)=hcc12(i,j)+v(i,j)
                    hcc12(j,i)=hcc12(i,j)
   20            continue
   10         continue
              if (nio.ne.0) then
                  do 30 i=1,nic
                     cntj=0
                     do 40 j=nbegi,nfnali
                        cntj=cntj+1
                        hco12(i,cntj)=hco12(i,cntj)+v(i,j)
                        hoc12(cntj,i)=hco12(i,cntj)
   40                continue
   30             continue
              endif
          endif
          if( nio.ne.0) then
              cnti=0
              do 50 i=nbegi,nfnali
                 cnti=cnti+1
                 cntj=0
                 do 60 j=nbegi,i
                    cntj=cntj+1
                    hoo12(cnti,cntj)=hoo12(cnti,cntj)+v(i,j)
                    hoo12(cntj,cnti)=hoo12(cnti,cntj)
   60            continue                                                     
   50         continue
          endif
      else          
         if (nic.ne.0.and.njc.ne.0) then      
             do 70 i=1,nic
                do 80 j=1,njc
                   hcc12(i,j)=hcc12(i,j)+v(i,j)
                   hcc21(j,i)=hcc12(i,j)
   80           continue
   70        continue
         endif
         if (nic.ne.0.and.njo.ne.0) then         
             do 90 i=1,nic
                cntj=0
                do 100 j=nbegj,nfnalj
                   cntj=cntj+1
                   hco12(i,cntj)=hco12(i,cntj)+v(i,j)
                   hoc21(cntj,i)=hco12(i,cntj)
  100           continue
   90        continue
         endif
         if (nio.ne.0.and.njc.ne.0) then         
             cnti=0
             do 110 i=nbegi,nfnali
                cnti=cnti+1 
                do 120 j=1,njc
                   hoc12(cnti,j)=hoc12(cnti,j)+v(i,j)
                   hco21(j,cnti)=hoc12(cnti,j)
  120           continue
  110        continue
         endif   
         if (nio.ne.0.and.njo.ne.0) then
              cnti=0
              do 200 i=nbegi,nfnali
                 cnti=cnti+1
                 cntj=0
                 do 210 j=nbegj,nfnalj
                    cntj=cntj+1
                    hoo12(cnti,cntj)=hoo12(cnti,cntj)+v(i,j)
                    hoo21(cntj,cnti)=hoo12(cnti,cntj)
  210            continue
  200         continue
         endif
      endif                                               
      return
      end
