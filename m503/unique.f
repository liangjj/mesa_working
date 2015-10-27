*deck @(#)unique.f	1.3  7/30/91
      subroutine unique(c,eigval,f,jmat,kmat,values,d,u,t1,t2,h,
     $                  nbf,nnp,nshell,ntriang,iout,nmat,
     $                  shlmin,shlmax,ops,s,smhalf,fcoef,
     $                  triang,t5,temp,usesym,nirrep,numso,lambda,
     $                  orbsym,salc)
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf),eigval(nbf),f(nnp),jmat(nnp,nmat),kmat(nnp,nmat)
      real*8 values(nnp,ntriang),d(nnp,nmat),u(nbf,nbf),t1(nbf,nbf)
      real*8 t2(nbf,nbf,nmat),h(nnp),s(nnp),smhalf(nnp),fcoef(nshell)
      real*8 triang(nnp),t5(nbf,nbf),salc(nbf,nbf)
      integer shlmin(nshell),shlmax(nshell)
      integer temp(nbf),orbsym(nbf)
      integer numso(nirrep),lambda(nirrep)
      character*(*) ops
      logical logkey,usesym
c
c     ----- form the density matrix or matrices -----
c
      do 124 i=1,nshell-1
         call gdmat(d(1,i),c,nbf,nnp,shlmin(i),shlmax(i))
  124 continue
c
c     ----- and the coulomb and exchange matrices -----
c
      call jandks(values,d,nnp,nbf,jmat,kmat,nmat,nmat,ntriang,
     $            t1,nmat,t2)
c
c     ----- transform the j and k matrices to the current mo basis ----
c
      do 1 i=1,nmat
         call trtosq(t1,jmat(1,i),nbf,nnp)
         call ebc(t2,t1,c,nbf,nbf,nbf)
         call ebtc(t1,c,t2,nbf,nbf,nbf)
         call sqtotr(jmat(1,i),t1,nbf,nnp)
c
         call trtosq(t1,kmat(1,i),nbf,nnp)
         call ebc(t2,t1,c,nbf,nbf,nbf)
         call ebtc(t1,c,t2,nbf,nbf,nbf)
         call sqtotr(kmat(1,i),t1,nbf,nnp)
    1 continue
c
c     ----- and construct the averaged fock matrix
c                f=h+f(i)[2j(i)-k(i)]
c
      call trtosq(t1,h,nbf,nnp)
      call ebc(t2,t1,c,nbf,nbf,nbf)
      call ebtc(t1,c,t2,nbf,nbf,nbf)
      call sqtotr(f,t1,nbf,nnp)
c
      do 3 shell=1,nshell-1
         do 2 i=1,nnp
            f(i)=f(i)+fcoef(shell)*(2*jmat(i,shell)-kmat(i,shell))
    2    continue
    3 continue
c
c     ----- zero out the sub-blocks mixing orbital shells -----
c
      do 7 ishell=2,nshell
         do 6 jshell=1,ishell-1
            do 5 i=shlmin(ishell),shlmax(ishell)
               ia=i*(i-1)/2
               do 4 j=shlmin(jshell),shlmax(jshell)
                  f(ia+j)=0.0d+00
    4          continue
    5       continue
    6    continue
    7 continue
c
c     ----- block and diagonalize the matrix by symmetry -----
c     ----- then put eigenvectors back in proper place   -----
c
      call rzero(u,nbf*nbf)       
      if (usesym) then
          begin=0
          do 9 irrep=1,nirrep
             tstnum=0
             do 10 i=1,nbf
                if ( orbsym(i).eq.irrep) then
                     tstnum=tstnum+1
                     temp(tstnum)=i
                endif
   10        continue
             if (tstnum.ne.numso(irrep)) then
                 call lnkerr ('error in symorb count in unique')
             endif
             if (tstnum.ne.0) then
                 count=0
                 do 20 i=1,tstnum
                    do 30 j=1,i
                       ij=temp(i)*(temp(i)-1)/2 +temp(j)
                       count=count+1
                       h(count)=f(ij)
   30               continue
   20            continue
                 call degrsp(tstnum,count,h,eigval(begin+1),1,t5,t1,t2)
                 begin=begin+tstnum
                 call tmptou(u,t5,temp,nbf,tstnum)
             endif
    9     continue
      else
             call degrsp(nbf,nnp,f,eigval,1,u,t1,t2)
      endif
c
c     ----- transform back to the ao basis -----
      call ebc(t1,c,u,nbf,nbf,nbf)
      call vmove(c,t1,nbf**2)
c
      if (logkey(ops,'print=scf=core-lagrangian',.false.,' ')) then
         write (iout,1000)
 1000    format (5x,'mo core lagrangian')
         call matout(eigval,nbf,1,nbf,1,iout)
      end if
c
      call iosys('write real "scf mo core lagrangian" to rwf',
     $           nbf,eigval,0,' ')
c
c
c
      return
      end
