*deck @(#)srt.f	1.1 9/7/91
c***begin prologue     srt
c***date written       890907   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m6011, link 6011, orbital sort
c***author             schneider, barry (lanl)
c***source             m6011
c***purpose            sort initial numerical orbital file
c***                   into a subset kept for m6005
c***description        
c***                   
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       srt
      subroutine srt (fin,fout,grid,ncon,numin,numout,nkept,list,npts,
     1                nwrin,nwdin,nwrout,nwdout,norm,ireg)
      implicit integer (a-z)
      common /io/ inp, iout
      real *8 fin, fout, grid, norm
      dimension fin (npts,numin), fout(npts,numout), list(200)
      dimension grid(4,npts), norm(numout)
      call iosys ('read real "trns grid" from grid without rewinding',
     1            4*npts,grid,0,' ')
      npass=ncon/numin
      nlast=ncon-npass*numin
      if (nlast.ne.0) then
          npass=npass+1
      else
          nlast=numin
      endif
      conf=0
      conl=0
      count=0
      wrdout=0
      cntout=0
      do 10 pass=1,npass     
         noin=numin
         if (pass.eq.npass) then
             noin=nlast
         endif
         wrdin=noin*npts
         call iosys ('read real "con array" from orbs '//
     1               'without rewinding',wrdin,fin,0,' ') 
         nwrin=nwrin+1
         nwdin=nwdin+wrdin
         conf=conl+1
         conl=conl+noin         
         concnt=0
         do 20 con=conf,conl
            concnt=concnt+1 
            do 30 tstc=1,nkept
               if (con.eq.list(tstc)) go to 40
   30       continue
            go to 50
   40       count=count+1
            cntout=cntout+1
            call copy(fin(1,concnt),fout(1,count),npts)
            do 200 ipt=1,npts            
               norm(cntout)=norm(cntout)+fout(ipt,count)*
     1                      fout(ipt,count)*grid(4,ipt)
  200       continue
            nwdout=nwdout+npts
            wrdout=wrdout+npts 
            if (count.eq.numout) then
                call iosys ('write real "sorted orbs" on orbs '//
     1                      'without rewinding',wrdout,fout,0,' ')
                nwrout=nwrout+1
                count=0
                wrdout=0 
            endif
   50       continue
   20    continue
   10 continue
      if (wrdout.ne.0) then
          call iosys ('write real "sorted orbs" on orbs without '//
     1                'rewinding',wrdout,fout,0,' ')
          nwrout=nwrout+1
      endif
      return
      end
