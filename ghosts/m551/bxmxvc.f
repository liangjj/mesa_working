*deck %W%  %G%
      subroutine bxmxvc(b,t,cv,c0,len,lok,mix,iblock,ipnt,
     $     buf,cr,icr)
C
C***Begin prologue     bxmxvc
C***Date written       871022   (yymmdd)
C***Revision date      871116   (yymmdd)
C
C  16 November 1987    BHL at LANL
C   core allocation bug fixed.
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
C
C***References
C
C***Routines called    (none)
C
C***End prologue       bxmxvc
C
      implicit real*8 (a-h,o-z)
c
      character*8 citype
      character*8 mcscf
c
      dimension nbf(2),nob(2),ncob(2),naob(2),cv(2),c0(2)
      dimension locsym(2),len(2),lok(2),mix(2),nfob(2)
      dimension iblock(14,2),ipnt(2)
      dimension buf(2),cr(2),icr(2)
      dimension b(2),t(2),ihd(35)
      real*8 bufix(lbufso),sg(*),tg(*)
      real*8 rabcx(*),rabix(*),raibx(*),hss(*)
      integer temp
c
      integer wpadti
c
      common /number/ zero,pt5,one,two,four,eight
      common /io/ inp,iout
c
c
      return
c
      entry enmxvc(nsym,nbf,nob,nfob,ncob,naob,cv,c0,
     $     locsym,len,lok,mix,
     $     lij,lijkl,
     $     nda1,lda1,nf41,lenb,nf35,nf36,nf37,nf16,
     $     nfsm,ndshd1,ndshd2,nhd,
     $     buf,cr,icr,ncore)
c
      return
c
c
      entry mcmxvc(b,nmix,mdim,energy,noci,iter,t,bufix,lbufso,
     $     rabcx,rabix,raibx,hss,incor,sg,tg)
c
      call iosys('read character mcscf_type from mcscr',-1,0,0,mcscf)
      call iosys('read character "mcscf: ci used" from rwf',
     $     -1,0,0,citype)
c
      ncsf = mdim - nmix
c-----------------------------------------------------------------------
c
c     b * delta
c
c-----------------------------------------------------------------------
c
      do 210 i=1,mdim
         t(i)=0.d0
 210  continue
      if(incor.lt.2) then
         call hesvec(t,b,hss,lenb,nmix)
      else
         call ebc(t,hss,b,nmix,nmix,1)
      end if
c
      if(noci.eq.1) return
c
c-----------------------------------------------------------------------
c
c     a * gamma
c
c-----------------------------------------------------------------------
      zero=0.d0
      intape=2
      iflag=1
      junk=0
      if (citype.eq.'m902') then
         call mn820(cr,cr,ncore,'mcscr','mcscf','ci')
         call mn902(cr,cr,ncore,'h times c',b(nmix+1),t(nmix+1),
     $              'mcscr','ci',junk,'mc')
      else if (citype.eq.'m903') then
         call mn903(cr,cr,ncore,'h times c',b(nmix+1),t(nmix+1),
     $        'mcscr','ci',junk,'mc')
      end if
c
      do 310 i = 1, ncsf
 310  t(nmix+i) = t(nmix+i) - energy * b(nmix+i)
c
c-----------------------------------------------------------------------
c
c     c * delta
c
c-----------------------------------------------------------------------
 400  lstor1 = 1
      ncor1 = ncore - lstor1
c-------------------------------
c     generate updated integrals
c-------------------------------
      mmda1=lda1
      intape=50
      call mcgvca(nsym,nbf,nob,nfob,ncob,naob,cv,
     $     locsym,len,lok,mix,b,lijkl,
     $     nda1,lda1,nf35,nf36,nf37,ltrb,
     $     cr(lstor1),icr(lstor1),ncor1,intape,bufix,lbufso,
     $     rabcx,rabix,raibx,incor,sg)
      lda1=mmda1
c
c------------------------------------
c     call h*c with updated integrals
c------------------------------------
      ncor2 = ncore-ncsf
      nc1=ncsf+1
      if (citype.eq.'m902') then
         call mn820(cr(nc1),cr(nc1),ncor2,'mcscr','mcscf','mcscf')
         call mn902(cr(nc1),cr(nc1),ncor2,'h times c',c0,
     $              cr,'mcscr','mcscf',junk,'mc')
      else if (citype.eq.'m903') then
         call mn903(cr(nc1),cr(nc1),ncor2,'h times c',c0,cr,
     $        'mcscr','mcscf',junk,'mc')
      end if
c
      call vadd(t(nmix+1),t(nmix+1),cr,ncsf)
c
c-----------------------------------------------------------------------
c
c     c(t) * gamma
c
c-----------------------------------------------------------------------
c
      lg1 = lstor1
      lg2 = lg1 + lij
      lstor2 = lg2 + lijkl
      ncor2 = ncore - lstor2
c-----------------------------------------
c     generate transition density matrices
c-----------------------------------------
c
      if (citype.eq.'m902') then
         call mn902(cr(lstor2),cr(lstor2),ncor2,
     $              'transition density',c0,b(nmix+1),'mcscr','mcscf',
     $              junk,'mc')
         call mn830(cr(lstor2),cr(lstor2),ncor2,'mcscr','mcscf',
     $        'guga transition density matrix','mo 1ptdm',
     $        'mo 2ptdm')
      else if (citype.eq.'m903') then
         call mn903(cr(lstor2),cr(lstor2),ncor2,'transition density',
     $        c0,b(nmix+1),'mcscr',' ',junk,'mc')
      end if
c
c---------------------------------------------
c     call cg with transition density matrices
c---------------------------------------------
c
c
c  ----- first .. read transition density
c
      temp=lg2+lijkl
      need=wpadti(temp+lij**2)
      call getscm(need,icr,junk,'DCAN',0)
c
      call iosys('read real "mo 1ptdm" from mcscr',lij,cr(lg1),0,' ')
      call iosys('read real "mo 2ptdm" from mcscr',
     $     lij**2,cr(temp),0,' ')
c
      call dcan(junk,cr(temp),cr(lg2),naob(1),lij,(lij+1)*lij/2)
c
c    1pdm from guga needs to be scaled
c
      call scdm1(cr(lg1),naob(1),lij)
c
      lg10 = lg1 - 1
      lg20 = lg2 - 1
c
      call mccga(nsym,nbf,nob,nfob,ncob,naob,cv,
     $     locsym,len,lok,mix,cr(lg1),cr(lg2),
     $     nda1,lda1,nf35,nf36,
     $     t,cr(lstor2),icr(lstor2),ncor2,bufix,lbufso,
     $     rabcx,rabix,raibx,incor,sg,tg)
c
      lda1=mmda1
c
c
      return
      end
