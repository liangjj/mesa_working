*deck ihamxt.f
c***begin prologue     ihamxt
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            calculate and store non-zero matrix elements
c***                   of hamiltonian and their indices in dvr
c***                   representation.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       ihamxt
      subroutine ihamxt(hr,ht,eigr,eigt,psi0,hbuf,ibuf,ind,diag,rhs,
     1                  hbar,omegat,omega,lenbuf,nply,dim,matdim,
     2                  pottyp,prnt,ntot,incore)
      implicit integer (a-z)
      real*8 hr, eigr, eigt, psi0, hbar, omegat, omega
      complex*16 ht, hbuf, diag, ham, rhs, tmp
      character*(*) pottyp
      character*1 itoc
      logical prnt, incore
      dimension nply(2)
      dimension eigt(nply(1)), eigr(nply(2))
      dimension hr(nply(2),nply(2)), ht(nply(1),nply(1)), psi0(nply(2))
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(matdim)
      dimension ind(matdim,*), rhs(matdim)
      common/io/inp, iout
      call iosys('create integer "'//itoc(dim)//'d buffers" '//
     1           'on ham',-1,0,0,' ')
      call czero(diag,matdim)
      call vscale(hr,hr,hbar*omegat,nply(2)*nply(2))
      ntot=0
      count=0
      do 10 i=1,matdim
            i1=ind(i,1)
            j1=ind(i,2)
            indi=ind(i,3)
            do 20 j=1,i-1
               ham=(0.d0,0.d0)
               i2=ind(j,1)
               j2=ind(j,2)
               indj=ind(j,3)
               if(j1.eq.j2) then
                  ham = ham + ht(i1,i2)
                  if(i1.eq.i2) then
                     ham = ham - hr(j1,j1)
                  endif
               else
                  if(i1.eq.i2) then
                     ham = ham - hr(j1,j2)
                  endif
               endif
               if(ham.ne.(0.d0,0.d0)) then
                  count=count+1
                  if(count.gt.lenbuf) then
                     ntot=ntot+lenbuf
                     call iosys('write integer "'//itoc(dim)//
     1                          'd buffers" to ham without rewinding',
     2                           2*lenbuf,ibuf,0,' ')
                     call iosys('write integer "'//itoc(dim)//
     1                          'd buffers" to ham without rewinding',
     2                           wptoin(2*lenbuf),hbuf,0,' ')
                     count=1
                  endif
                  ibuf(1,count)=indi
                  ibuf(2,count)=indj
                  hbuf(count)=ham
               endif                        
 20         continue
            diag(indi) = diag(indi) + ht(i1,i1) - hr(j1,j1)
 10      continue
c     now add in the potential energy contribution.
      do 30 i=1,dim
         i1=ind(i,1)
         j1=ind(i,2)
         indi = ind(i,3)
         tmp=eigr(j1)*cos(omega*eigt(i1))   
         diag(indi)=diag(indi) - tmp
         rhs(indi) = tmp*psi0(j1)       
 30   continue
      if(count.gt.0) then
         ntot=ntot+count
         call iosys('write integer "'//itoc(dim)//
     1              'd buffers" to ham without rewinding',2*count,
     2               ibuf,0,' ')
         call iosys('write integer "'//itoc(dim)//
     1              'd buffers" to ham without rewinding',
     2               wptoin(2*count),hbuf,0,' ')
      endif
      call iosys('endfile "'//itoc(dim)//'d buffers" on ham',0,0,0,' ')
      call iosys('write integer "'//itoc(dim)//'d number of elements '
     1           //'" to ham',1,ntot,0,' ')
      call iosys('write real "'//itoc(dim)//'d diagonals '
     1          //'" to ham',2*matdim,diag,0,' ')
      call iosys('rewind all on ham read-and-write',0,0,0,' ')
      write(iout,2) ntot, matdim
      incore=.false.
      if(ntot.le.lenbuf) then
         incore=.true.
         write(iout,3)
      endif         
      if(prnt) then
         call rdham(hbuf,ibuf,diag,dim,lenbuf,matdim)
      endif
      return
 1    format(a80)     
 2    format(/,1x,'number non-zero,off-diagonal matrix elements = ',i6,
     1       /,1x,'number of diagonal matrix elements = ',i5)
 3    format(/,1x,'all non-zero matrix elements and indices kept in'
     1            ' core')     
      end       






