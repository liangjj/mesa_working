*deck @(#)ciout.f	5.1  11/6/94
      subroutine ciout(dtsefi,dtsefj,number,detdet,detsef,sefsef,
     1                 row,cimat,msymtp)
c
c  transform blocks of determinant matrix elements to
c  double-group-function matrix elements.  write out
c  rows of the hamiltonian matrix.
c
      implicit real*8 (a-h,o-z)
      character*80 rtitle, ctitle, blabel
c
      common /headng/ rtitle, ctitle, blabel
      common/icntrl/ iciwrt, icipun, intwrt, ihamrd, ksym
      common /core/ ecore,thresh,ev,mask,nmask,limit1,limit2,intmxo
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      common /c4/iex(4),jex(4),jsymi,jsymj,nopeni,nopenj,ndeti,ndetj
      common /c5/ i,j,nsefi,nsefj,isef,jsef,nzero,kmax,kmaxb
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /detio/ idfrst, idspc, irt, iau, inextu, intmax, ihso
      dimension dtsefi(maxdet,*),dtsefj(maxdet,*),number(*),
     1  detdet(ndeti,*),detsef(ndeti,*),sefsef(nsefi,*),row(*),
     2  cimat(nsefi,*)
      integer and,or
c
*mdc*if cray
*      equivalence (next,rnext)
*mdc*else
      dimension nexta(2)
      equivalence (nexta(1),rnext), (nexta(2),next)
      save rnext, nexta
*mdc*endif
c
      dimension isoch(4,4), mult(4,4)
      data isoch/2*0,1,0,1,0,3*1,3*0,3*1,0/
      data mult/1,2,3,4, 2,1,4,3, 3,4,1,2, 4,3,2,1/
c
c  variables determined in this routine:
c
c     detdet     the hamiltonian matrix between all determinants
c                arising from configurations i and j
c     detsef     the half-transformed hamiltonian matrix over
c                double-group functions
c     sefsef     the hamiltonian matrix between all double-group
c                functions arising from configurations i and j
c
c  variables redefined in this routine:
c
c     dtsefi     coefficients for the expansion of the
c                double-group functions in terms of determinants for
c                configuration i
c     dtsefj     coefficients for the expansion of the
c                double-group functions in terms of determinants for
c                configuration j
c
c  transform the hamiltonian matrix block to the
c  double-group-adapted functions of configuration j
c
      if(mel.eq.1) go to 200
c
      if(nopenj.eq.0) go to 120
c
      do 100 jjj=1,nsefj
      do 100 ii=1,ndeti
  100 detsef(ii,jjj)=0.0d0
      do 110 jjj=1,nsefj
      do 110 jj=1,ndetj
      do 110 ii=1,ndeti
  110 detsef(ii,jjj)=detsef(ii,jjj)+dtsefj(jj,jjj)*detdet(ii,jj)
c
c  transform the hamiltonian matrix block to the
c  double-group-adapted functions of configuration i
c
  120 if(nopeni.eq.0) go to 150
c
      do 130 jjj=1,nsefj
      do 130 iii=1,nsefi
  130 sefsef(iii,jjj)=0.0d0
      do 140 jjj=1,nsefj
      do 140 ii=1,ndeti
      do 140 iii=1,nsefi
  140 sefsef(iii,jjj)=sefsef(iii,jjj)+dtsefi(ii,iii)*detsef(ii,jjj)
c
c  put in phase factors for spin-orbit matrix elements.
c
  150 if(msymtp.eq.1) then
        do 160 jjj=1,nsefj
        do 160 iii=1,nsefi
  160   if(isoch(mod(iii-1,4)+1,mod(jjj-1,4)+1).eq.1)
     1    sefsef(iii,jjj)=-sefsef(iii,jjj)
      elseif(msymtp.eq.2) then
        if(jsymi.eq.ksym) then
          iadd=1
        else
          iadd=3
        endif
        if(jsymj.eq.ksym) then
          jadd=1
        else
          jadd=3
        endif
        do 170 jjj=1,nsefj
        do 170 iii=1,nsefi
  170   if(isoch(mod(iii-1,2)+iadd,mod(jjj-1,2)+jadd).eq.1)
     1    sefsef(iii,jjj)=-sefsef(iii,jjj)
      else
        mji=mod(jsymi-1,4)+1
        mjj=mod(jsymj-1,4)+1
        mk =mod(ksym-1,4)+1
        if(isoch(mult(mji,mk),mult(mjj,mk)).eq.1) then 
          do 180 jjj=1,nsefj
          do 180 iii=1,nsefi
  180     sefsef(iii,jjj)=-sefsef(iii,jjj)
        endif
      endif
c
c  assemble the hamiltonian matrix row-by-row from the sefsef
c  matrix element blocks.
c
  200 do 220 ii=1,nsefi
c      write(iw,*) 'nsefi',nsefi
c      write(iw,*) 'row',(row(irlm),irlm=1,3)
c
c         write(iw,*)'ciout: isef,jsef',isef,jsef
c
         if(isef.eq.jsef) then
            jmax=ii
         else
            jmax=nsefj
         endif
         ncol=jsef
c
         do 210 jj=1,jmax
            rnext=sefsef(ii,jj)
            if(abs(rnext).ge.thresh.or.(isef.eq.jsef.and.ii.eq.jj)) then
               nzero = nzero + 1
               number(ii)=number(ii)+1
               if( number(ii).gt.maxne0 ) then
                  iauu = iau+max(kmax,number(ii))
     $                      +nsefi*number(ii)+1-intmax
                  write (iw,900) iauu, i, ii, j, jj
                  stop
               endif
*mdc*if cray
*        cimat(ii,number(ii))=or(and(rnext,nmask),ncol)
*mdc*elseif sun
               next= or( and(next,nmask),ncol)
               cimat(ii,number(ii))=rnext
c
c               write(iw,*)'ciout: ncol,next,ii,number(ii),cimat',
c     >                     ncol,next,ii,number(ii),cimat(ii,number(ii))
c
*mdc*else ibm vax stellar
*        next=ior(iand(next,nmask),ncol)
*        cimat(ii,number(ii))=rnext
*mdc*endif
            endif
  210    ncol=ncol+1
c
  220 continue
c      write(iw,*) 'row',(row(irlm),irlm=1,3)
c
      jsef=jsef+nsefj
c      write(iw,*) 'i,j',i,j
      if(j.lt.i) return
      
c
*mdc*if ibm vax sun stellar
      nexta(1)=0
*mdc*endif
      if(isef.ne.1) then
         next=number(1)+1
         nlast=kmax+1
         row(nlast)=rnext
c
c         write(iw,*)'ciout: number(1),next,kmax,nlast,row(nlast)',
c     >            number(1),next,kmax,nlast,row(nlast)
c
         call hout(row,nlast)
c
c         write(iw,*)'ciout: hout 1'
      endif
c
  230 if(nsefi.ne.1) then
         nn=nsefi-1
         do 250 ii=1,nn
            kmax=number(ii)
            kmaxb = max(kmaxb,kmax)
            next=number(ii+1)+1
            do 240 kk=1,kmax
  240          row(kk)=cimat(ii,kk)
            nlast=kmax+1
            row(nlast)=rnext
  250    call hout(row,nlast)
c
c         write(iw,*)'ciout: hout 2'
      endif
c
  260 kmax=number(nsefi)
      kmaxb = max(kmaxb,kmax)
      do 270 kk=1,kmax
c         write(iw,*) 'nsefi,isef,jsef,kk',nsefi,isef,jsef,kk
c         write(iw,*) 'cimat(nsfi,kk',cimat(nsefi,kk)
  270    row(kk)=cimat(nsefi,kk)
c
      isef=isef+nsefi
      if (i.ge.ntotfg) then
         nlast=kmax+1
         row(nlast)=rnext
         call hout(row,nlast)
c
c         write(iw,*)'ciout: hout 3'
      endif
c
c
c      write(iw,*) 'leaving ciout: row',(row(irlm),irlm=1,3)
      return
c
  900 format(/' more space needed in ciout.  iaup at least',i12//
     1  10x,'matrix element being processed:'/
     2  2(/10x,'spatial configuration',i5,5x,'spin function',i5))
c
      end
