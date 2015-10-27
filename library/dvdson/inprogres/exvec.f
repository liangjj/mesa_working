*deck exvec.f
c***begin prologue     exvec
c***date written       970531   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           davidson, trial, vector
c***author             schneider, barry (nsf)
c***source             
c***purpose            set up new trial vector based an exact or approximate
c***                   solution to the residual equation.
c***references         
c
c***routines called    
c***end prologue       exvec
      subroutine exvec(ham,eig,hbuf,ibuf,diag,resid,vec,
     1                              temp,ipvt,n,m,nel,iter,prnt)
      implicit integer (a-z)
      real*8 ham, eig, hbuf, diag, resid, vec, temp
      logical prnt
      character*4 itoc
      character*80 title  
      dimension ham(n,n), eig(n), hbuf(*), ibuf(2,*), diag(n)
      dimension resid(n,m), vec(n,m), temp(n,n), ipvt(n)
      common/io/inp, iout
      call rzero(ham,n*n)
      do 10 i=1,n
         ham(i,i) = ham(i,i) + diag(i)
 10   continue
      do 20 i=1,nel
         iel = ibuf(1,i)
         jel = ibuf(2,i)
         ham(iel,jel) = ham(iel,jel) + hbuf(i)
         ham(jel,iel) = ham(jel,iel) + hbuf(i)                  
 20   continue
      do 1000 i=1,m
         call rzero(temp,n*n)
          do 30 j=1,n
             temp(j,j) = ham(j,j) - eig(i)
 30       continue
          temp(1,2)=ham(1,2)
          temp(n,n-1)=ham(n,n-1)
          do 40 j=2,n-1
             temp(j,j+1)=ham(j,j+1)
             temp(j,j-1)=ham(j,j-1)
 40       continue   
          call vneg(temp,temp,n*n)
          call sgefa(temp,n,n,ipvt,info)
          if(info.ne.0) then
             call lnkerr('error in sgefa')
          else
             call sgesl(temp,n,n,ipvt,resid(1,i),0)
             call copy(resid(1,i),vec(1,i),n)
          endif
 1000 continue                                  
      if(prnt) then
         title='new trial vectors iteration = '//itoc(iter)
         call prntrm(title,vec,n,m,n,m,iout)
      endif
      return
      end       





