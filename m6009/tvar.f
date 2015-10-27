*deck @(#)tvar.f	1.1 9/8/91
c***begin prologue     tvar
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           tvar, link 1106, kohn variational
c***author             schneider, barry (lanl), rescigno, tom(llnl)  
c***source             m1106
c***purpose            find variational t-matrix
c***description        
c***                   
c***                   
c***                  
c***references         schneider and rescigno, physical review
c
c***routines called    iosys, util and mdutil
c***end prologue       tvar
      subroutine tvar(vrr,vrrc,vir,tvir,vfrb,vbfr,vbfrc,tmat,seig,svec,
     1                dum,dir,type,ntchn,matbb,fil,typem)
      implicit integer (a-z)
      real *8 vrr, vfrb, vbfr, phase, rtest, dum
      real *8 atan2, epsum, impart, repart
      complex *16 vir, tvir, tmat, cfac, svec, seig, catan
      complex *16 vrrc, vbfrc
      logical dir
      character *80 title
      character *24 fil
      character *(*) type, typem
      dimension vrr(ntchn,ntchn), vfrb(ntchn,matbb), vir(ntchn,ntchn)
      dimension tvir(ntchn,ntchn), vbfr(matbb,ntchn)
      dimension seig(ntchn), svec(ntchn,ntchn), dum(3*ntchn)
      dimension tmat(ntchn,ntchn), vbfrc(matbb,ntchn), vrrc(ntchn,ntchn)
      common /io/ inp, iout
      data cfac/ (0.d+00,2.d+00)/
c----------------------------------------------------------------------------c
c        now form variationally stable t-matrix                              c
c----------------------------------------------------------------------------c
      if (type.eq.'partitioned') then
          if (typem.eq.'real') then
              do 10 i=1,ntchn
                 do 20 j=1,ntchn
                    tmat(j,i)=vrr(i,j)
   20            continue
   10         continue     
              call ebc(vrr,vfrb,vbfr,ntchn,matbb,ntchn)
              do 30 i=1,ntchn
                 do 40 j=1,ntchn
                    tmat(i,j)=tmat(i,j)-vrr(i,j)
   40            continue
   30         continue              
          elseif(typem.eq.'complex') then
              do 50 i=1,ntchn
                 do 60 j=1,ntchn
                    tmat(j,i)=vrr(i,j)
                    vrrc(i,j)=vrr(i,j)
   60            continue              
   50         continue     
              call ebcc(vrrc,vfrb,vbfrc,ntchn,matbb,ntchn)
              do 70 i=1,ntchn
                 do 80 j=1,ntchn
                    tmat(i,j)=tmat(i,j)-vrrc(i,j)
   80            continue
   70         continue              
          endif
          call cambtc(tmat,tvir,vir,ntchn,ntchn,ntchn)
      endif
      do 90 i=1,ntchn
         do 100 j=1,ntchn
            tmat(i,j)=-2.d+00*tmat(i,j)
  100    continue
   90 continue     
      if (.not.dir) then
         title='t-matrix'
      else
         title='k-matrix'
      endif
      write (iout,250)
      call iosys ('write real '//fil//' to tmat',2*ntchn*ntchn,
     1            tmat,0,' ')
      call prntcm(title,tmat,ntchn,ntchn,ntchn,ntchn,iout)     
      if (.not.dir) then
         do 150 i=1,ntchn
            do 140 j=1,ntchn
               tmat(i,j)=cfac*tmat(i,j)
  140       continue
  150    continue
         do 160 j=1,ntchn
            tmat(j,j)=tmat(j,j)+1.d+00
  160    continue   
      endif
      job=1
      call cgeev (tmat,ntchn,ntchn,seig,svec,ntchn,dum,job,info)
      write (iout,200)
      if (.not.dir) then
         epsum=0.d+00
         do 170 i=1,ntchn
            impart=imag(seig(i))
            repart=real(seig(i))
            phase=atan2(impart,repart)
            phase=phase*.5d+00
            rtest=seig(i)*conjg(seig(i))
            write (iout,210) i,phase,rtest
            epsum=epsum+phase
  170    continue
      else
         epsum=0.d+00
         do 180 i=1,ntchn
            phase=catan(seig(i))
            rtest=1.d+00
            write (iout,210) i,phase,rtest
            epsum=epsum+phase
  180    continue
      endif
      write (iout,300) epsum
      return
  200 format (/,10x,' eigenphases of s matrix')
  210 format (/,2x,'phase no.',1x,i3,2x,'phase =',2x,e15.8,2x,'modulus =
     1 ',2x,f10.5)
  250 format (/,10x,'variationally corrected results')
  300 format (/,5x,'eigenphase sum:',1x,e15.8)
      end
