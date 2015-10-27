*deck @(#)symorb.f	1.3  9/3/91
      subroutine symorb(movec,salc,symsalc,s,eval,hold,
     1                  t1,t2,t3,tol,maxerr,mosym,nummo,
     2                  molst,nsym,num,symtre,nsmall,prn)
      implicit integer(a-z)
      real*8 movec, salc, symsalc, s, eval, hold, t1, t2, t3
      real*8 tol, maxerr
      logical prn, symtre
      character*8 symtyp
      character*4 subgrp
      dimension movec(num,num), salc(num,num), symsalc(num,num)
      dimension s(num,num), eval(num,2)
      dimension t1(num,num), t2(num,num)
      dimension t3(num,num), mosym(nsym,2), nummo(nsym)
      dimension molst(num,nsym,2), hold(num,num), prn(*)
      dimension symtyp(14)
      character*16 bflabl(300)
      character *1 itoc
      character*80 title
c
      common /io/ inp, iout
c
      call iosys('read character "basis function labels" from rwf',
     $              -1,0,0,bflabl)
      call iosys('read character "group symbol" from rwf',0,0,0,subgrp)
      call grpsym(subgrp,symtyp,nirrep)
      call iosys('read integer "number of irreducible representations"'
     $             //' from rwf',1,nirrep,0,' ')
c       call iosys('read character "labels of irreducible '//
c     $            'representations" from rwf',nirrep*len(symtyp(1)),
c     $             0,0,lirrep)
      nnp=num*(num+1)/2
      call iosys('read real "overlap integrals" from rwf',nnp,t1,0,' ')
      call trtosq(s,t1,num,nnp)
      call rzero(symsalc,num*num)
      if(prn(1)) then
         write(iout,*) ' input molecular orbitals '
         call wvec(movec,eval,num,num,bflabl,' ')
      endif
c
      loc=1
      cnt=0
      do 10 i=1,nsym      
         mosym(i,2)=0
         nummo(i)=0
         if(mosym(i,1).ne.0) then
c
c        orthonormalize the salcs for this symmetry
c
            write(iout,1) symtyp(i)
            call lowdin(s,salc(1,loc),symsalc(1,loc),t3,t1,t2,
     1                  tol,mosym(i,1),mosym(i,2),num,prn(2))    
            if(mosym(i,1).ne.mosym(i,2)) then
               write(iout,2) symtyp(i), mosym(i,1),
     #                       symtyp(i), mosym(i,2)
               call lnkerr('linear dependence error')
            endif 
c
c        form the projection operator
c
            call ebctxx(t1,symsalc(1,loc),symsalc(1,loc),num,
     1                  mosym(i,1),num,num,num,num)
c
c        calculate the coefficients of the projected mo's
c
            call ebc(t2,t1,s,num,num,num)
            call ebc(t1,t2,movec,num,num,num)
c
c        calculate the overlap matrix of the projected set
c
            call ebtc(t2,t1,s,num,num,num)
            call ebc(t3,t2,t1,num,num,num)
c
c        a diagonal element of "zero" indicates that this mo has no
c        projection in this symmetry and should not be counted.
c        
            if(prn(3)) then
               title='projected vector overlap matrix'
               call prntrm(title,t3,num,num,num,num,iout)
            endif   
            do 20 j=1,num
               if(abs(t3(j,j)).gt.maxerr) then
                  nummo(i)=nummo(i)+1
                  molst(nummo(i),i,1)=j
               endif
 20         continue
            if(nummo(i).ne.0) then
               do 30 j=1,nummo(i)
                  cnt=cnt+1
                  eval(cnt,2)=eval(molst(j,i,1),1)
                  call copy(t1(1,molst(j,i,1)),hold(1,cnt),num)
 30            continue   
            endif
            loc=loc+mosym(i,1)
         endif
 10   continue 
      call copy(hold,movec,num*num)
      call copy(eval(1,2),eval(1,1),num)
      call icopy(molst(1,1,1),molst(1,1,2),num*nsym)
c
c     the orbitals are now symmetry ordered.  molst contains the original location
c     of each orbital in the transformation matrix.
c
      if(symtre) then
         write(iout,6)  
         cnt=0 
         cntvrt=nsmall
         do 50 i=1,nsym
            if(nummo(i).ne.0) then
               cntsym=0
               do 60 j=1,nummo(i)
                  flg=0
                  do 70 k=1,nsmall
                     if(molst(j,i,2).eq.k) then
c
c                       this mo goes up front
c
                        flg=1
                        cnt=cnt+1
                        call copy(hold(1,cnt),movec(1,k),num)
                        eval(k,1)=eval(cnt,2)
                        cntsym=cntsym+1
                        molst(cntsym,i,1)=k
                     endif                     
 70               continue
                  if(flg.eq.0) then
                     cnt=cnt+1
                     cntvrt=cntvrt+1
                     call copy(hold(1,cnt),movec(1,cntvrt),num)
                     eval(cntvrt,1)=eval(cnt,2)
                     cntsym=cntsym+1
                     molst(cntsym,i,1)=cntvrt
                  endif
 60            continue   
            endif
 50      continue
      endif
      do 100 i=1,nsym
         if(nummo(i).ne.0) then       
            write(iout,3) symtyp(i), nummo(i)
            write(iout,4)
            write(iout,5) (molst(j,i,1),j=1,nummo(i))            
            call iosys ('write integer "symmos-'//itoc(i)//'" to rwf',
     1                   nummo(i),molst(1,i,1),0,' ')
            call iosys ('write integer "symmos-'//itoc(i)//'" to chk',
     1                   nummo(i),molst(1,i,1),0,' ')
         endif
 100  continue   
c
c      perform a final schmidt orthogonalization to clean up any 
c      non-orthonormality issues.
c
      call wvec(movec,eval,num,num,bflabl,' ')
      call sqtotr(t1,s,num,nnp)
      call schmdt(movec,t1,s,t2,t3,num,num,nnp,maxerr)
c      if(prn(1)) then
         call wvec(movec,eval,num,num,bflabl,' ')
c      endif
      call iosys ('write integer "kohn mo numsym" to rwf',nsym,nummo,
     1             0,' ')
      call iosys ('write integer "kohn mo numsym" to chk',nsym,nummo,
     1       0,' ')
      call iosys ('write integer "scattering mo numsym" to rwf',
     #             nsym,nummo,0,' ')
      call iosys ('write integer "scattering mo numsym" to chk',
     #             nsym,nummo,0,' ')
      return
 1    format(/,1x,'Lowdin orthogonalization for symmetry = ',a3)
 2    format(/,1x,'number of ',a3,' orbitals input  = ',i3,/,1x,
     1            'number of ',a3,' orbitals output = ',i3)
 3    format(/,1x,'number of molecular orbitals of ',a3,
     #            ' symmetry = ',i3)
 4    format(/,1x,'list of orbitals')
 5    format((/,1x,( 10(i3,1x) ) ))
 6    format(/,1x,'Reordering Molecular Orbitals')
      end
