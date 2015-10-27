*deck @(#)symorb.f	5.1  11/6/94
      subroutine symorb(cvec,salc,temp,mosave,aosave,
     1                  hold,tt,proj, s,xm, mov,aov, eval,tsym,
     2                  mosym,numsym, nsym,num, smllmo,tstmo,
     3                  smleig, symtre,maxerr,nsmall,kw,notest,pvec)
      implicit integer(a-z)
      real*8 cvec(num,num),smllmo,tstmo,xx,temp(num,*),
     1 salc(num,num),proj(num*num,nsym),tt(num*num),eval(num),smleig,
     2 tsym(nsym,nsym),hold(num,*),s(num*num),xm(num*num),maxerr
      integer aosave(num,nsym),mosave(num,nsym),mov(num),aov(num)
      integer numsym(nsym),mosym(nsym)
      character*16 bflabl(300)
      character *1 itoc
      logical pvec
c
      common /io/ inp, iout
c
      call iosys('read character "basis function labels" from rwf',
     $              -1,0,0,bflabl)
c
c   orthonormalize the salcs
c
      nnp=num*(num+1)/2
      call iosys('read real "overlap integrals" from rwf',
     $ nnp,xm,0,' ')
      call schmdt(salc,xm,s,temp,tt,num,num,nnp,maxerr)
      if(pvec) then
         write(iout,*)' input vectors '
         call wvec(cvec,eval,num,num,bflabl,' ')
      end if
c
      eval(num+1)=-100000.d+00
      eval(num+2)=-100000.d+00
c
c
c   construct the symmetry projection operators,  ps
c
c    ps = |cs><cs|s  where cs are symmetry vectors and
c                          s  is  the overlap matrix
c
      ip=1
      do 301 i=1,nsym
         numsym(i)=0
         mm=mosym(i)
         if(mm.ne.0) then
            call ebct(temp,salc(1,ip),salc(1,ip),num,mm,num)
            call ebc(proj(1,i),temp,s,num,num,num)
         end if
         ip=ip+mm
 301  continue
c
      do 203 i=1,num
         mov(i)=0
 203  continue
c
       ij=1
1000   continue
       i=ij
       i1=i+1
       i2=i+2
       if(abs(eval(i)-eval(i1)).lt.smleig) then
          ndeg=2
          ij=ij+2
          if(abs(eval(i)-eval(i2)).lt.smleig) then
             ndeg=3
             ij=ij+1
          end if
       else
          ndeg=1
          ij=ij+1
       end if
c
       kk=i
       call rzero(tt,num)
       do 1 k=1,ndeg
          do 2 j=1,num
             tt(j)=tt(j)+cvec(j,kk)
  2       continue
          kk=kk+1
  1    continue
 
       do 3 j=1,nsym
          if(mosym(j).ne.0) then
             call ebc(temp(1,j),proj(1,j),tt,num,num,1)
          end if
  3    continue
c
       call ebc(hold,s,temp,num,num,nsym)
       call ebtc(tsym,temp,hold,nsym,num,nsym)
c
       ns=0
       do 4 j=1,nsym
          if(mosym(j).ne.0) then
             if(tsym(j,j).gt.smllmo) then
                mov(ns+i)=j
                ns=ns+1
                numsym(j)=numsym(j)+1
                xx=1.d0/sqrt(tsym(j,j))
                do 5 k=1,num
                   hold(k,ns)=temp(k,j)*xx
  5             continue
             end if
          end if
  4    continue
c
c
        if(ns.ne.ndeg) then
           write(iout,*)' ndeg ne ns .. stop ',ns,ndeg
           write(iout,*)' i i1 i2 ',i,i1,i2
           call lnkerr(' m604: symmetry projection ')
        end if
c
        if(notest.ne.0) then
           call tstsym(hold,cvec(1,i),s,temp,tt,num,ndeg,tstmo)
        end if
c
        call scopy(ndeg*num,hold,1,cvec(1,i),1)
c
       if(ij.le.num) go to 1000
c
       do 15 i=1,nsym
          if(mosym(i).ne.numsym(i)) then
             write(iout,*)' symmetry error ',mosym(i),numsym(i)
             call lnkerr(' m604: sym error ')
          end if
  15   continue
c
       ks=0
       do 100 i=1,nsym
c
          numsym(i)=0
          if(mosym(i).ne.0) then
             do 7 j=1,num
                aov(j)=0
  7          continue
c
             do 9 j=ks+1,ks+mosym(i)
                do 6 k=1,num
                   if(abs(salc(k,j)).gt.smllmo) then
                      aov(k)=1
                   end if
  6             continue
  9          continue
             ks=ks+mosym(i)
             jj=0
             do 8 j=1,num
               if(aov(j).gt.0) then
                  jj=jj+1
                  aosave(jj,i)=j
               endif
 8           continue
             numsym(i)=jj
          end if
c
100    continue
c
c
c
      do 201 i=1,nsym
         numsym(i)=0
 201  continue
c
      if(symtre.eq.0) then
         do 200 i=1,num
            is=mov(i)
            numsym(is)=numsym(is)+1
            mosave(numsym(is),is)=i
 200  continue
      else
         do 210 i=1,nsmall
            is=mov(i)
            aov(i)=mov(i)
            tt(i)=eval(i)
            numsym(is)=numsym(is)+1
            mosave(numsym(is),is)=i
 210     continue
         call scopy(nsmall*num,cvec,1,hold,1)
         ii=nsmall
         do 400 i=nsmall+1,num
            is=mov(i)
            if(is.eq.symtre) then
               numsym(is)=numsym(is)+1
               ii=ii+1
               mosave(numsym(is),is)=ii
               aov(ii)=is
               tt(ii)=eval(i)
               call scopy(num,cvec(1,i),1,hold(1,ii),1)
            end if
 400     continue
         do 410 j=1,nsym
            if(j.ne.symtre) then
               do 420 i=nsmall+1,num
                  is=mov(i)
                  if(is.eq.j) then
                     numsym(is)=numsym(is)+1
                     ii=ii+1
                     mosave(numsym(is),is)=ii
                     aov(ii)=is
                     tt(ii)=eval(i)
                     call scopy(num,cvec(1,i),1,hold(1,ii),1)
                  end if
  420          continue
            end if
 410     continue
         do 430 i=1,num
            mov(i)=aov(i)
            eval(i)=tt(i)
 430     continue
         call scopy(num*num,hold,1,cvec,1)
      end if
c
      write(iout,106)
 106  format(/,' kohnopt input: written to rwf and chk ')
      call iosys ('write integer "kohn mo numsym" to rwf',nsym,mosym,
     1             0,' ')
      call iosys ('write integer "kohn mo numsym" to chk',nsym,mosym,
     1             0,' ')
      write(iout,*) ' '
      write(iout,*) 'mo orbital symmetry information'
      do 107 i=1,nsym
         if(numsym(i).ne.mosym(i)) then
            write(iout,*)' symmetry error in mo section '
            call lnkerr(' m604: sym error ')
         end if
         if(mosym(i).gt.0) then
            write (iout,*) ' '
            write (iout,*) 'symmetry', i
            write(iout,108) mosym(i)
            write(iout,109) (mosave(jj,i),jj=1,mosym(i))
            call iosys ('write integer "symmos-'//itoc(i)//'" to rwf',
     1                   mosym(i),mosave(1,i),0,' ')
            call iosys ('write integer "symmos-'//itoc(i)//'" to chk',
     1                   mosym(i),mosave(1,i),0,' ')
         endif
 107  continue
 103  format(10(2x,i3))
 108  format(/,10x,'number of symmetry orbitals in this irrep',1x,i4)
 109  format(/,5x,'list of mos',(/,5x,10(i4,1x)))
c
      call schmdt(cvec,xm,s,temp,tt,num,num,nnp,maxerr)
c
      write(iout,110)
 110  format(/,' final vectors ')
      call wvec(cvec,eval,num,num,bflabl,' ')
c
      call iosys ('rewind all on chk read-and-write',0,0,0,' ')
      call iosys ('close chk',0,0,0,' ')
      return
      end
