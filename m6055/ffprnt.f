*deck @(#)ffprnt.f	1.1 9/8/91
c***begin prologue     ffprnt
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            print free-free matrix elements
c***references
c
c***routines called
      subroutine ffprnt(hpvhp,hpvhp1,hpvhm,hpvhm1,ovpp,ovpm,nlm,
     1                  nchan,nstri,dimc,maxlm,vdrctv)
      implicit integer (a-z)
      character *80 title
      character *8 colt, rowt
      character *3 itoc
      complex *16 hpvhp, hpvhm, hpvhp1, hpvhm1, ovpp, ovpm
      real *8 rowv
      logical vdrctv
      common /io/ inp, iout
      dimension hpvhp(maxlm,maxlm,nstri), hpvhp1(maxlm,nchan)
      dimension hpvhm(maxlm,maxlm,nchan,nchan), hpvhm1(maxlm,nchan)
      dimension nlm(dimc), ovpp(maxlm,nchan), ovpm(maxlm,nchan)
      if (.not.vdrctv) then
          ist=0
          do 10 ch1=1,nchan
             do 20 ch2=1,ch1
                ist=ist+1
                title='hpvhp matrix ch1-'//itoc(ch1)//' ch2-'//itoc(ch2)
                rowt='lm index'
                colt=rowt
                rowv=-99.d0
                colv=-99
                call cmprir(hpvhp(1,1,ist),rowv,colv,nlm(ch1),nlm(ch2),
     1                      maxlm,maxlm,title,rowt,colt,iout)
   20        continue
   10     continue
          do 30 ch1=1,nchan
             do 40 ch2=1,nchan
                title='hpvhm matrix ch1-'//itoc(ch1)//' ch2-'//itoc(ch2)
                rowt='lm index'
                colt=rowt
                rowv=-99.d0
                colv=-99
                call cmprir(hpvhm(1,1,ch1,ch2),rowv,colv,nlm(ch1),
     1                      nlm(ch2),maxlm,maxlm,title,rowt,colt,iout)
   40        continue
   30     continue
      else
          do 50 ch1=1,nchan
             title='hpvhp1 matrix ch1-'//itoc(ch1)
             rowt='lm index'
             colt=rowt
             rowv=-99.d0
             colv=-99
             call cmprir(hpvhp1(1,ch1),rowv,colv,nlm(ch1),1,
     1                   maxlm,1,title,rowt,colt,iout)
   50     continue
             do 60 ch1=1,nchan
             title='ovpp matrix ch1-'//itoc(ch1)
             rowt='lm index'
             colt=rowt
             rowv=-99.d0
             colv=-99
             call cmprir(ovpp(1,ch1),rowv,colv,nlm(ch1),1,
     1                   maxlm,1,title,rowt,colt,iout)
   60     continue
          do 70 ch1=1,nchan
             title='ovpm matrix ch1-'//itoc(ch1)
             rowt='lm index'
             colt=rowt
             rowv=-99.d0
             colv=-99
             call cmprir(ovpm(1,ch1),rowv,colv,nlm(ch1),1,
     1                   maxlm,1,title,rowt,colt,iout)
   70     continue
      endif
      return
      end
