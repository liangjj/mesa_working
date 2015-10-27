*deck fildag
c***begin prologue     fildag
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           hamiltonian, matrix
c***author             schneider, barry (nsf)
c***source             
c***purpose            fill diagonal blocks of the hamiltonian matrix.
c***description        four types of channel blocks can appear
c***                   corresponding to the closed-closed, closed-open,
c***                   open-closed and open-open portions of the hamiltonian
c***                   matrix. the hamiltonian matrix is passed for all four
c***                   parts of this block in the calling routine.  the t
c***                   matrix for the block is passed in its entirety.  it is
c***                   assumed that the t matrix is blocked in the same fashion
c***                   with all closed orbitals first in the list.
c***                   
c
c***routines called
c***end prologue       fildag
      subroutine fildag(hcc,hco,hoc,hoo,t,ec,nic,nio,ma,n)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 hcc, hco, hoc, hoo, t, ec
      dimension hcc(n,n), hco(n,n), hoc(n,n), hoo(n,n)
      dimension t(ma,*)
      nbeg=nic+1
      nfnal=nic+nio
      if (nic.ne.0) then
          do 10 i=1,nic
             hcc(i,i)=hcc(i,i)+ec
             do 20 j=1,i
                hcc(i,j)=hcc(i,j)+t(i,j)
                hcc(j,i)=hcc(i,j)
   20        continue
   10     continue
          if (nio.ne.0) then
              do 30 i=1,nic
                 cntj=0
                 do 40 j=nbeg,nfnal
                    cntj=cntj+1
                    hco(i,cntj)=hco(i,cntj)+t(i,j)
                    hoc(cntj,i)=hco(i,cntj)
   40            continue
   30         continue
          endif
      endif
      if( nio.ne.0) then
          cnti=0
          do 50 i=nbeg,nfnal
             cnti=cnti+1
             hoo(cnti,cnti)=hoo(cnti,cnti)+ec
             cntj=0
             do 60 j=nbeg,i
                cntj=cntj+1
                hoo(cnti,cntj)=hoo(cnti,cntj)+t(i,j)
                hoo(cntj,cnti)=hoo(cnti,cntj)
   60        continue                                                     
   50     continue
      endif
      return
      end
