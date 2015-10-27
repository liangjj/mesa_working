*deck @(#)srtv.f	3.1  11/20/92
       subroutine srtv(c,ct,nbf,norb,symtyp,eig,teig,id2h,mompt,
     $ localz,prtd2h,thrshv)
       implicit integer (a-z)
       real*8 c(nbf,norb),ct(nbf,norb),eig(nbf),teig(nbf),thrshv
       real*8 xtest,phase,root2
       dimension symtyp(norb),mompt(nbf)
       dimension nsym(8),symstr(4),symorb(8)
c
       common/io/inp,iout
c
c
       call izero(nsym,8)
c
       do 1 i=1,norb
        nsym(symtyp(i))=nsym(symtyp(i))+1
   1   continue
c
       symstr(1)=0
       do 2 i=2,4
        symstr(i)=symstr(i-1)+nsym(i-1)
   2   continue
c
       do 3 i=1,norb
       js=symstr(symtyp(i))+1
       teig(js)=eig(i)
         do 4 j=1,nbf
          ct(j,js)=c(j,i)
   4     continue
       symstr(symtyp(i))=js
   3   continue
c
       write(iout,5)(nsym(i),i=1,4)
   5   format(/,'  symmetry orbitals detected in this class ',/,
     # 4i6)
c
       numsym=4
c
       if(id2h.ne.0) then
c
c.       write(iout,*)' c2v orbitals '
c.       call matout(ct,nbf,nbf,nbf,nbf,iout)
c
       numsym=8
       write(iout,6)
   6   format(/,'  d2h symmetry orbitals will be generated  ',/)
c
       call iosys('read integer mompt from rwf',nbf,mompt,0,' ')
c
       norb2=norb/2
       nbf2=nbf/2
       ioff=0
       iu=norb2
       ig=0
       do 10 isym=1,4
        nobs=nsym(isym)
c.        write(iout,*)' symmetry nobs ',isym,nobs
        ius=0
        igs=0
        if(nobs.eq.0) go to 10
         do 9 i=1,nobs
          is=ioff+i
c
c find the element with the largest absolute value..scilib call
c
          jtest=isamax(nbf2,ct(1,is),1)
c
          phase=(-1.0)**mompt(jtest)
c23456
c.       write(iout,*)' orbital jtest phase ',i,jtest,phase,
c.     $  mompt(jtest)
c.       write(iout,*)'terms ',ct(jtest,is),ct(jtest+nbf2,is)
c
         xtest=abs(ct(jtest,is)-phase*ct(jtest+nbf2,is))
     $         / abs(ct(jtest,is))
c
c.       write(iout,*)'  xtest  thrshv ',xtest,thrshv
c
        if(xtest.gt.thrshv) then
c.         write(iout,*)' u-symmetry '
         iu=iu+1
         ius=ius+1
         eig(iu)=teig(is)
         call scopy(norb,ct(1,is),1,c(1,iu),1)
        else
c.         write(iout,*)' g-symmetry '
         ig=ig+1
         igs=igs+1
         eig(ig)=teig(is)
         call scopy(norb,ct(1,is),1,c(1,ig),1)
        end if
c
 9     continue
       ioff=ioff+nobs
       if(igs.ne.ius) then
        write(iout,11) isym,igs,ius
        call lnkerr(' srtv: d2h error ')
       else
        write(iout,11) isym,igs,ius
       end if
c
       symorb(isym)=igs
       symorb(isym+4)=ius
c
10     continue
c
       call scopy(norb,eig,1,teig,1)
       call scopy(nbf*norb,c,1,ct,1)
       call scopy(8,symorb,1,nsym,1)
c
11     format(' in symmetry ',i4,' g-orbs=  ',i6,'  u-orbs=  ',i6)
c
       end if
c
       write(iout,*)' writing bhlnsym ',numsym
       call iosys('write integer bhlnsym to rwf',1,numsym,0,' ')
       call iosys('write integer bhlnobs to rwf',numsym,nsym,0,' ')
       call iosys('write real "symmetry orbitals" to rwf',
     $  nbf*nbf,ct,0,' ')
c
       if(localz.ne.0) then
       write(iout,*)'  localizing the orbitals '
       if(prtd2h.ne.0) then
        write(iout,*)' delocalized orbitals in d2h '
        call matout(ct,nbf,nbf,nbf,nbf,iout)
       end if
       ix=0
       do 9011 i=1,nbf2
        jt1=isamax(nbf2,ct(1,i),1)
        if(ct(jt1,i).lt.0.) then
         do 7012 j=1,nbf
         ct(j,i)=-ct(j,i)
 7012    continue
        end if
        if(ct(jt1,i+nbf2).lt.0.) then
         do 8012 j=1,nbf
         ct(j,i+nbf2)=-ct(j,i+nbf2)
 8012    continue
        end if
       ix=ix+1
       teig(ix)=eig(i)
        do 9012 j=1,nbf
        c(j,ix)=ct(j,i)+ct(j,i+nbf2)
 9012   continue
       ix=ix+1
       teig(ix)=eig(i+nbf2)
        do 9013 j=1,nbf
        c(j,ix)=ct(j,i)-ct(j,i+nbf2)
 9013   continue
 9011  continue
       root2=1.d0/sqrt(2.d0)
       do 9014 i=1,nbf
        do 9015 j=1,nbf
        ct(j,i)=c(j,i)*root2
 9015   continue
 9014  continue
c
c.       write(iout,*)' debug print '
c.       call matout(ct,nbf,4,nbf,4,iout)
c
       end if
c
       return
       end
