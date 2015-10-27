*deck @(#)rdints.f	5.1  11/6/94
      subroutine rdints(valone,valtwo,valx,valy,valz)
c
c  read the mo integrals (from the transformation program) and store
c  them according to their indices
c
      implicit real*8 (a-h,o-z)
      character*32 lab1
c mesa
c
      character*128 tints
c
*mdc*if cray
*      parameter (n2=32, n4=16)
*mdc*else
      parameter (n2=16, n4=8)
*mdc*endif
      common/icntrl/ iciwrt, icipun, intwrt, ihamrd, ksym
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      common /core/ ecore,thresh,ev,mask,nmask,limit1,limit2,intmxo
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /detio/ idfrst, idspc, irt, iau, inextu, intmax, ihso
      common /io/inp,iout
c
c mesa
c
      dimension valone(*),valtwo(*),valx(*),valy(*),valz(*)
c      dimension valone(*),valtwo(*),v(intmxo),lbl(intmxo),lb(4)
c
c  zero out the arrays to be used for the one-electron integrals
c
      do 10 i=1,limit1
        valx(i)=0.d0
        valy(i)=0.d0
        valz(i)=0.d0
   10 valone(i) = 0.0d0
c
c  read in the one-electron integrals from tape iunt1b
c
c mesa read
c
       call iosys('read character "transformed integral filename"'
     $           //' from rwf',0,0,0,tints)
       call iosys('open tints as old',0,0,0,tints)
       call iosys('read real "mo one-electron integrals" from tints',
     $-1,valone,0,' ')
c
c read spin orbit integrals
c
c  lx
c
       call iosys('read real "msox integrals" from tints',
     $-1,valx,0,' ')
c       if(logkey(ops,'print=int=so',.false.,' '))then
c         write(iout,*)'rdints x',(valx(i),i=1,limit1)
c       end if
c
c  ly
c
       call iosys('read real "msoy integrals" from tints',
     $-1,valy,0,' ')
c       if(logkey(ops,'print=int=so',.false.,' '))then
c         write(iout,*)'rdints y',(valy(i),i=1,limit1)
c       end if
c
c  lz
c
       call iosys('read real "msoz integrals" from tints',
     $-1,valz,0,' ')
c       if(logkey(ops,'print=int=so',.false.,' '))then
c         write(iout,*)'rdints z',(valz(i),i=1,limit1)
c       end if
c
c add contributions from valx, valy, and valz to valone array
c
c       do 100 i=1,limit1
c         valone(i)=valone(i)+valx(i)+valy(i)+valz(i)
c100    continue
c
c      if(ihso.eq.0) then
c        intsu=1
c        if(intwrt.ge.1) lab1='h+(core)'
c      else
c        intsu=4
c        if(intwrt.ge.1) lab1='h+(core), hso(x), hso(y), hso(z)'
c      endif
c
c      do 50 ints=1,intsu
c
c   30 read (iunt1b) last, nints, v, lbl
c      do 40 m=1,nints
c      call unpck(lbl(m),lb,n2,2)
c   40 valone((lb(1)*(lb(1)-1))/2+lb(2)) = v(m)
c      if(last.eq.0) go to 30
c
c   50 continue
c
      if(intwrt.ge.1)  lab1='h+(core)'
      if(intwrt.ge.1) call oneout(lab1,valone,nbf,5,3)
c
c  zero out the array to be used for the two-electron integrals
c
      do 70 i=1,limit2
   70 valtwo(i) = 0.0d0
c
c  read the two-electron integrals from tape iunt1b
c
c mesa read
c
      do 61 i=1,nbf
         do 61 j=1,i
            do 61 k=1,nbf
               do 61 l=1,k
       call iosys('read real "mo two-electron integrals"'// 
     $' from tints without rewinding',
     $1,temp,0,' ')
                  valtwo(indexf(i,j,k,l)) = temp
   61 continue
      call iosys('close tints',0,0,0,' ')
c
c   80 read (iunt1b) last, nints, v, lbl
c      do 90 m=1,nints
c      call unpck(lbl(m),lb,n4,4)
c   90 valtwo(indexf(lb(1),lb(2),lb(3),lb(4))) = v(m)
c      if(last.eq.0) go to 80
c
       write(iw,*)'two - elec ints: rdints 1 -4',(valtwo(i),i=1,4)
       write(iw,*)'two - elec ints: rdints 5 -8',(valtwo(i),i=5,8)
       write(iw,*)'two - elec ints: rdints 9 -12',(valtwo(i),i=9,12)
       write(iw,*)'two - elec ints: rdints 13 -16',(valtwo(i),i=13,16)
       if(intwrt.ge.2) call twoout('el. rep.',valtwo,nbf,4)
c
      return
c
c
      end
