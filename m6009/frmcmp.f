*deck @(#)frmcmp.f	1.1 9/8/91
c***begin prologue     frmcmp
c***date written       880630   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           frmcmp, link 1107, kohn variational
c***author             schneider, barry (lanl), rescigno, tom(llnl)  
c***source             m6009
c***purpose            calculate full complex matrix from scattering
c***                   integrals and then the t matrices.
c***description        assembles complete matrix for kohn solution
c***                   without partitioning off the bound-bound
c***                   component first. t matrices also extracted.
c***                   useful for checks and for a more direct method.
c***references         schneider and rescigno, physical review
c
c***routines called    iosys, util and mdutil
c***end prologue       frmcmp
c
      subroutine frmcmp(vrr,vrrc,vii,vir,vri,vbfr,vbfrc,vbfi,vfrb,
     1                  vfib,hambb,hambbc,matcmp,crhs,clhs,ipvt,tmat,
     2                  seig,svec,dum,dir,matbb,ntchn,n,prh,prs,
     3                  prn,fil,typem)
      implicit integer (a-z)
      character *80 title
      character *24 fil
      logical dir, prh, prs, prn
      real *8 vrr, vbfr, vfrb, hambb, dum
      complex *16 vii, vri, vir, vbfi, vfib, matcmp, crhs, clhs
      complex *16 vrrc, vbfrc, hambbc 
      complex *16 tmat, seig, svec
      character *24 type
      character *(*) typem
      dimension vrr(ntchn,ntchn), vii(ntchn,ntchn), vri(ntchn,ntchn)
      dimension vir(ntchn,ntchn), vbfr(matbb,ntchn), vbfi(matbb,ntchn)
      dimension vfrb(ntchn,matbb), vfib(ntchn,matbb), hambb(matbb,matbb)
      dimension matcmp(n,n), crhs(n,ntchn), clhs(n,ntchn), ipvt(n)
      dimension tmat(ntchn,ntchn), svec(ntchn,ntchn)
      dimension seig(ntchn), dum(3*ntchn), hambbc(matbb,matbb)
      dimension vrrc(ntchn,ntchn), vbfrc(matbb,ntchn)
      common /io/ inp, iout
      call czero(matcmp,n*n)
      call czero(crhs,n*ntchn)
c----------------------------------------------------------------------c
c         fill in the full matrix using the basic integrals            c
c----------------------------------------------------------------------c
c         free-free 
      icnt=0
      do 10 i=matbb+1,n
         icnt=icnt+1
         jcnt=0
         do 20 j=matbb+1,n
            jcnt=jcnt+1
            matcmp(i,j)=vii(icnt,jcnt)
            crhs(i,jcnt)=vir(icnt,jcnt)
   20    continue
   10 continue
c      free-bound
      icnt=0
      do 30 i=matbb+1,n
         icnt=icnt+1
         do 40 j=1,matbb
            matcmp(i,j)=vfib(icnt,j)
   40    continue
   30 continue
c     bound-free
      if (typem.eq.'complex') then
          do 50 i=1,matbb
             jcnt=0
             do 60 j=matbb+1,n
                jcnt=jcnt+1
                matcmp(i,j)=vbfi(i,jcnt)
                crhs(i,jcnt)=vbfrc(i,jcnt)
   60        continue
   50     continue
      else
          do 55 i=1,matbb
             jcnt=0
             do 65 j=matbb+1,n
                jcnt=jcnt+1
                matcmp(i,j)=vbfi(i,jcnt)
                crhs(i,jcnt)=vbfr(i,jcnt)
   65        continue
   55     continue
      endif
c     bound-bound
      if (typem.eq.'complex') then
          do 70 i=1,matbb
             do 80 j=1,matbb
                matcmp(i,j)=hambbc(i,j)
   80        continue
   70     continue
      else
          do 75 i=1,matbb
             do 85 j=1,matbb
                matcmp(i,j)=hambb(i,j)
   85        continue
   75     continue
      endif
      if (prh) then
          title='full matrix'
          call prntcm(title,matcmp,n,n,n,n,iout)
          title='full complex rhs'
          call prntcm(title,crhs,n,ntchn,n,ntchn,iout)
      endif
      call cc2opy(crhs,clhs,n*ntchn)
      call cgefa (matcmp,n,n,ipvt,info)
      do 90 i=1,ntchn
         call cgesl(matcmp,n,n,ipvt,crhs(1,i),0)
   90 continue
      if (prs) then
          title='full complex solution'
          call prntcm(title,crhs,n,ntchn,n,ntchn,iout)
      endif
      call matm(crhs,n*ntchn)
      icnt=0      
      do 100 i=matbb+1,n
         icnt=icnt+1
         jcnt=0
         do 110 j=matbb+1,n
            jcnt=jcnt+1
            vir(icnt,jcnt)=crhs(i,jcnt)
  110    continue
  100 continue     
c----------------------------------------------------------------------c
c                extract non variational t matrices                    c
c----------------------------------------------------------------------c
      call tnonvr (vir,tmat,seig,svec,dum,dir,ntchn)
c----------------------------------------------------------------------c
      call cebtc(vri,clhs,crhs,ntchn,n,ntchn)
      if (prn) then
          title='full complex numerator'
          call prntcm(title,vri,ntchn,ntchn,ntchn,ntchn,iout)
      endif
      if(typem.eq.'complex') then
         do 200 i=1,ntchn
            do 210 j=1,ntchn
               tmat(i,j)=vrr(i,j)+vri(i,j)
  210       continue
  200    continue    
      else
         do 300 i=1,ntchn
            do 310 j=1,ntchn
               tmat(i,j)=vrr(i,j)+vri(i,j)
  310       continue
  300    continue    
      endif
c----------------------------------------------------------------------c
c              extract variational t matrices                          c
c----------------------------------------------------------------------c
      type='not partitioned'
      call tvar(vrr,vrrc,vir,vri,vfrb,vbfr,vbfrc,tmat,seig,svec,dum,
     1          dir,type,ntchn,matbb,fil,typem)
c----------------------------------------------------------------------c
      return
      end

