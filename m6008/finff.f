*deck @(#)finff.f	1.1 9/8/91
c***begin prologue     finff
c***date written       890605   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           finff, link 6008, kohn variational
c***author             schneider, barry (lanl)  
c***source             m6008
c***purpose            final free-free integrals
c***                   
c***description        the free-free integrals are transformed to a
c***                   free basis orthogonal to all bound molecular
c***                   orbitals.
c***references         schneider and rescigno, physical review
c
c***routines called    iosys, util and mdutil
c***end prologue       finff
      subroutine finff (hpvhp,hpvhm,hpvb,ovbf,ovpb,ovmb,hpb,hmb,hpp,
     1                  hpm,hmm,hambb,scrb,scrc,nbtot,orblst,finlst,
     2                  zeroc,nlm,nchan,ntchn,matbb,maxlm,nstri,nmo,
     3                  dimmo,dimc,prntff,prntbf,prntov,prntfn)
      implicit integer(a-z)
      character *80 title
      logical prntff, prntbf, prntov, prntfn, zeroc
      real*8 ovmb, hmb, hambb, hmm, rowv, scrb
      complex*16 ovpb, hpb, hpp, hpm, scrc
      complex *16 hpvhp, hpvhm, ovbf, hpvb
      character *8 rowt, colt
      dimension hpvhp(1:maxlm,1:maxlm,nstri)
      dimension hpvhm(1:maxlm,1:maxlm,nchan,nchan)
      dimension ovbf(1:maxlm,nmo,nchan), hpvb(1:maxlm,nmo,nchan,nchan)
      dimension ovpb(ntchn,matbb), ovmb(ntchn,matbb)
      dimension hpb(ntchn,matbb), hmb(ntchn,matbb)
      dimension hpp(ntchn,ntchn), hpm(ntchn,ntchn), hmm(ntchn,ntchn)
      dimension hambb(matbb,matbb), nbtot(dimc), orblst(dimmo,dimc)
      dimension finlst(dimmo,dimc), nlm(dimc), zeroc(dimc)
      dimension scrc(ntchn,matbb), scrb(ntchn,matbb)
      common /io/ inp,iout
c----------------------------------------------------------------------c
c          set up matrix elements in channel form                      c
c----------------------------------------------------------------------c    
      rowv=-99.d0
      colv=-99
      call czero(ovpb,ntchn*matbb)
      call rzero(ovmb,ntchn*matbb)
      call czero(hpb,ntchn*matbb)
      call rzero(hmb,ntchn*matbb) 
      cntch1=0
      do 10 ch1=1,nchan
         if (.not.zeroc(ch1)) then
             do 20 bfn=1,nbtot(ch1)
                cntb=finlst(bfn,ch1)
                orb=orblst(bfn,ch1)
                cntf=cntch1
                do 30 nolm1=1,nlm(ch1)
                   cntf=cntf+1
                   ovpb(cntf,cntb)=ovbf(nolm1,orb,ch1)
                   ovmb(cntf,cntb)=imag(ovpb(cntf,cntb))
   30           continue
   20        continue
         endif
         cntch1=cntch1+nlm(ch1)
   10 continue
      cntch1=0
      do 40 ch1=1,nchan
         if (.not.zeroc(ch1)) then 
             do 50 nolm1=1,nlm(ch1)
                cntch1=cntch1+1
                do 60 ch2=1,nchan
                   if (.not.zeroc(ch2)) then
                       do 70 bfn=1,nbtot(ch2)
                          cntb=finlst(bfn,ch2)
                          orb=orblst(bfn,ch2)
                          hpb(cntch1,cntb)=hpvb(nolm1,orb,ch1,ch2)
                          hmb(cntch1,cntb)=imag(hpb(cntch1,cntb))
   70                  continue
                   endif
   60           continue                  
   50        continue
         else
             cntch1=cntch1+nlm(ch1)
         endif  
   40 continue
c----------------------------------------------------------------------c
c               same for free-free integrals                           c
c----------------------------------------------------------------------c
      call czero(hpp,ntchn*ntchn)
      call czero(hpm,ntchn*ntchn) 
      call rzero(hmm,ntchn*ntchn) 
      cntch1=0
      do 100 ch1=1,nchan
         if (.not.zeroc(ch1)) then
             cntch2=0
             do 200 ch2=1,ch1
                if (.not.zeroc(ch2)) then
                    ist=ch1*(ch1-1)/2+ch2
                    cntf1=cntch1
                    do 300 nolm1=1,nlm(ch1)
                       cntf1=cntf1+1
                       cntf2=cntch2
                       do 400 nolm2=1,nlm(ch2)
                          cntf2=cntf2+1            
                          hpp(cntf1,cntf2)=hpvhp(nolm1,nolm2,ist)
                          hpp(cntf2,cntf1)=hpp(cntf1,cntf2)
                          hpm(cntf1,cntf2)=hpvhm(nolm1,nolm2,ch1,ch2)
                          hpm(cntf2,cntf1)=hpvhm(nolm2,nolm1,ch2,ch1)
                          hmm(cntf1,cntf2)=imag(hpm(cntf1,cntf2))
                          hmm(cntf2,cntf1)=imag(hpm(cntf2,cntf1))
  400                  continue
  300               continue
                endif
                cntch2=cntch2+nlm(ch2)
  200        continue
         endif  
         cntch1=cntch1+nlm(ch1)
  100 continue
c----------------------------------------------------------------------c
c                now we have the unique matrices                       c
c                    print if desired                                  c
c----------------------------------------------------------------------c 
      if (prntov) then
          title='ovpb'
          call cmprir(ovpb,rowv,colv,ntchn,matbb,ntchn,matbb,title,
     1                rowt,colt,iout)
          title='ovmb'
          call mprir(ovmb,rowv,colv,ntchn,matbb,ntchn,matbb,title,
     1               rowt,colt,iout)
      endif
      if (prntbf) then
          title='hpb-no'
          call cmprir(hpb,rowv,colv,ntchn,matbb,ntchn,matbb,title,
     1                rowt,colt,iout)
          title='hmb-no'
          call mprir(hmb,rowv,colv,ntchn,matbb,ntchn,matbb,title,
     1                rowt,colt,iout)
      endif
      if (prntff) then
          title='hpp-no'
          call cmprir(hpp,rowv,colv,ntchn,ntchn,ntchn,ntchn,title,
     1                rowt,colt,iout)
          title='hpm-no'
          call cmprir(hpm,rowv,colv,ntchn,ntchn,ntchn,ntchn,title,
     1                rowt,colt,iout)
          title='hmm-no'
          call mprir(hmm,rowv,colv,ntchn,ntchn,ntchn,ntchn,title,
     1               rowt,colt,iout)
      endif
c----------------------------------------------------------------------c
c              calculate the new matrix elements                       c
c                    the form is                                       c
c----------------------------------------------------------------------c
c                                                                      c 
c    m(f,f') = m(f0,f0') -sum <f0,b> m(b,f0') -sum m(f0,b) <b,f0')     c
c                           +                                          c
c                          sum sum <f0,b> m(b,b') <b',f0'>             c
c                                                                      c
c                            m = h - e                                 c
c                                                                      c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c                    hpp first                                         c
c----------------------------------------------------------------------c
c     ----- first sum -----
      call cambct(hpp,ovpb,hpb,ntchn,matbb,ntchn)
c     ----- second sum -----
      call cambct(hpp,hpb,ovpb,ntchn,matbb,ntchn)
c     ----- first matrix multiply of double sum -----
      call ecbc(scrc,ovpb,hambb,ntchn,matbb,matbb)
c     ----- final result back in original matrix -----
      call capbct(hpp,scrc,ovpb,ntchn,matbb,ntchn)
c----------------------------------------------------------------------c
c                    hmm now                                           c
c----------------------------------------------------------------------c
c     ----- first sum -----
      call ambct(hmm,ovmb,hmb,ntchn,matbb,ntchn)
c     ----- second sum -----
      call ambct(hmm,hmb,ovmb,ntchn,matbb,ntchn)
c     ----- first matrix multiply of double sum -----
      call ebc(scrb,ovmb,hambb,ntchn,matbb,matbb)
c     ----- final result back in original matrix -----
      call apbct(hmm,scrb,ovmb,ntchn,matbb,ntchn)
c----------------------------------------------------------------------c
c                    hpm now                                           c
c----------------------------------------------------------------------c
c     ----- first sum -----
      call amcbct(hpm,ovpb,hmb,ntchn,matbb,ntchn)
c     ----- second sum -----
      call amcbct(hpm,hpb,ovmb,ntchn,matbb,ntchn)
c     ----- first matrix multiply of double sum -----
c     ----- already done above  -------------------
c     ----- final result back in original matrix -----
      call apcbct(hpm,scrc,ovmb,ntchn,matbb,ntchn)
      if (prntfn) then
          title='hpp-o'
          call cmprir(hpp,rowv,colv,ntchn,ntchn,ntchn,ntchn,title,
     1                rowt,colt,iout)
          title='hpm-o'
          call cmprir(hpm,rowv,colv,ntchn,ntchn,ntchn,ntchn,title,
     1                rowt,colt,iout)
          title='hmm-o'
          call mprir(hmm,rowv,colv,ntchn,ntchn,ntchn,ntchn,title,
     1               rowt,colt,iout)
      endif
      return
      end
