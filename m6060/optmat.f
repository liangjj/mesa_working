*deck @(#)optmat.f
c***begin prologue     optmat
c***date written       920415   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           optmat, link 6060, kohn variational
c***author             schneider, barry (nsf)  
c***source             m6060
c***purpose            final optical potential matrix elements
c***                   
c***description        matrix elements in channel form
c
c***routines called    iosys, util and mdutil
c***end prologue       optmat
      subroutine optmat (vppin,vpmin,vmmin,vbpin,vbmin,vbbin,vppout,
     1                   vpmout,vmmout,vpbout,vmbout,vbbout,nbscat,
     2                   orblst,finlst,nlm,nchan,ntchn,maxlm,nmo,matbv,
     3                   itri,dimmo,dimc,prnt,zroset)
      implicit integer(a-z)
      character *80 title
      logical prnt, zroset
      real*8 rowv
      complex*16 vppin,vpmin,vmmin,vppout,vpmout,vmmout
      complex *16 vbpin,vbmin,vpbout,vmbout,vbbin,vbbout
      character *8 rowt, colt
      dimension vppin(1:maxlm,1:maxlm,itri), nbscat(dimc)
      dimension vpmin(1:maxlm,1:maxlm,nchan,nchan), orblst(dimmo,dimc)
      dimension vmmin(1:maxlm,1:maxlm,itri), nlm(dimc)
      dimension vbbin(nmo,nmo,itri), vbbout(matbv,matbv)
      dimension vppout(ntchn,ntchn), vpmout(ntchn,ntchn)
      dimension vmmout(ntchn,ntchn), finlst(dimmo,dimc)
      dimension vbpin(nmo,1:maxlm,nchan,nchan), vpbout(ntchn,matbv)
      dimension vbmin(nmo,1:maxlm,nchan,nchan), vmbout(ntchn,matbv)
      common /io/ inp,iout
c----------------------------------------------------------------------c
c          set up matrix elements in channel form                      c
c----------------------------------------------------------------------c    
      rowv=-99.d0
      colv=-99
      call czero(vppout,ntchn*ntchn)
      call czero(vpmout,ntchn*ntchn) 
      call czero(vmmout,ntchn*ntchn) 
      call czero(vpbout,ntchn*matbv) 
      call czero(vmbout,ntchn*matbv) 
      call czero(vbbout,matbv*matbv)
      if (.not.zroset) then
c----------------------------------------------------------------------c
c                   those involving bound-free functions               c
c----------------------------------------------------------------------c
          cntch1=0
          do 40 ch1=1,nchan
             do 50 nolm1=1,nlm(ch1)
                cntch1=cntch1+1
                do 60 ch2=1,nchan
                   do 70 bfn=1,nbscat(ch2)
                      cntb=finlst(bfn,ch2)
                      orb=orblst(bfn,ch2)
                      vpbout(cntch1,cntb)=vbpin(orb,nolm1,ch2,ch1)
                      vmbout(cntch1,cntb)=vbmin(orb,nolm1,ch2,ch1)
   70              continue
   60           continue                  
   50        continue
   40     continue
c----------------------------------------------------------------------c
c                    free-free functions                               c
c----------------------------------------------------------------------c
          cntch1=0
          do 100 ch1=1,nchan
             cntch2=0
             do 200 ch2=1,ch1
                ist=ch1*(ch1-1)/2+ch2
                cntf1=cntch1
                do 300 nolm1=1,nlm(ch1)
                   cntf1=cntf1+1
                   cntf2=cntch2
                   do 400 nolm2=1,nlm(ch2)
                      cntf2=cntf2+1            
                      vppout(cntf1,cntf2)=vppin(nolm1,nolm2,ist)
                      vppout(cntf2,cntf1)=vppout(cntf1,cntf2)
                      vmmout(cntf1,cntf2)=vmmin(nolm1,nolm2,ist)
                      vmmout(cntf2,cntf1)=vmmout(cntf1,cntf2)
                      vpmout(cntf1,cntf2)=vpmin(nolm1,nolm2,ch1,ch2)
                      vpmout(cntf2,cntf1)=vpmin(nolm2,nolm1,ch2,ch1)
  400              continue
  300           continue
                cntch2=cntch2+nlm(ch2)
  200        continue
  100     continue
c----------------------------------------------------------------------c
c                       bound-bound                                    c
c----------------------------------------------------------------------c
          do 500 ch1=1,nchan
             do 600 ch2=1,ch1
                ist=ch1*(ch1-1)/2+ch2
                do 700 bfn1=1,nbscat(ch1)
                   cnt1=finlst(bfn1,ch1)
                   orb1=orblst(bfn1,ch1)
                   do 800 bfn2=1,nbscat(ch2)
                      cnt2=finlst(bfn2,ch2)
                      orb2=orblst(bfn2,ch2)
                      vbbout(cnt1,cnt2)=vbbin(orb1,orb2,ist)
                      vbbout(cnt2,cnt1)=vbbout(cnt1,cnt2)
  800              continue
  700           continue
  600        continue
  500     continue
c----------------------------------------------------------------------c
c                now we have the unique matrices                       c
c                    print if desired                                  c
c----------------------------------------------------------------------c 
          if (prnt) then
              title='vppout'
              call cmprir(vppout,rowv,colv,ntchn,ntchn,ntchn,ntchn,
     1                    title,rowt,colt,iout)
              title='vpmout'
              call cmprir(vpmout,rowv,colv,ntchn,ntchn,ntchn,ntchn,
     1                    title,rowt,colt,iout)
              title='vmmout'
              call cmprir(vmmout,rowv,colv,ntchn,ntchn,ntchn,ntchn,
     1                    title,rowt,colt,iout)
              title='vpbout'
              call cmprir(vpbout,rowv,colv,ntchn,matbv,ntchn,matbv,
     1                    title,rowt,colt,iout)
              title='vmbout'
              call cmprir(vmbout,rowv,colv,ntchn,matbv,ntchn,matbv,
     1                    title,rowt,colt,iout)
              title='vbbout'
              call cmprir(vbbout,rowv,colv,matbv,matbv,matbv,matbv,
     1                    title,rowt,colt,iout)
          endif
      endif
      return
      end



