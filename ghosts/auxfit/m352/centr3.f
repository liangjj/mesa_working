deck  %W% %G%
      subroutine centr3(c,ex,z,iz,cont,s,ptprim,noprim,nocont,ptcont,
     #                  nat,nprim,maxcor,ntypes,nbtype,nnp,ncont,
     #                  start,nbasis,zan,nocart,nobf,maxmom,mintyp,
     #                  nx,ny,nz,minmom,ops,exaux,npaux,
     $                  xptprm,xnoprm,xnocon,xptcon,xstart,xcont,
     $                  nbasx,ncontx)
c***begin prologue     %M%
c***date written       840723  
c***revision date      %G%      
c
c***keywords           
c***author             saxe, paul and martin,richard(lanl) 
c***source             %W%   %G%
c***purpose            
c   module to form the three-center overlap integrals over generally
c   gaussian basis sets.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       %M% 
c
      implicit integer (a-z)
c
c     ----- input arrays(unmodified) -----
      character*(*) ops
      real*8 c(3,nat),ex(nprim),exaux(npaux)
      real*8 cont(ncont),xcont(ncontx)
      real*8 zan(nat)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes)
      integer nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(*),ny(*),nz(*)
      integer xptprm(nat,ntypes),xnoprm(nat,ntypes),xnocon(nat,ntypes)
      integer xptcon(nat,ntypes),xstart(nat,ntypes)
c
c     ----- output arrays -----
      real*8 s(nnp,nbasx)
c
c     ----- scratch arrays -----
      real*8 z(maxcor)
      integer iz(*)
c
c     ----- local variables -----
      real*8 pi,one,four
      logical debug
c
      parameter (debug=.false.)
      parameter (one=1.0d+00,four=4.0d+00)
c
      common/io/inp,iout
c
c
 1000 format(5x,'3-center overlap integrals:')
c
c     ----- form the overlap integrals -----
c
      pi=four*atan(one)
      do 100 iatom=1,nat
         do 90 jatom=1,nat
            do 80 katom=1,nat
               do 70 itype=1,nbtype
                  if (noprim(iatom,itype).le.0) go to 70
                  if (iatom.ne.jatom) then
                     jtypmx=nbtype
                  else
                     jtypmx=itype
                  end if
                  do 60 jtype=1,jtypmx
                     if (noprim(jatom,jtype).le.0) go to 60
c
                     imax=maxmom(itype)
                     jmax=maxmom(jtype)
                     nprimi=noprim(iatom,itype)
                     nprimj=noprim(jatom,jtype)
                     nconti=nocont(iatom,itype)
                     ncontj=nocont(jatom,jtype)
                     npint=nprimi*nprimj
                     lenblk=nocart(itype)*nocart(jtype)
                     do 50 ktype=1,nbtype
                        if(xnoprm(katom,ktype).le.0) go to 50
c
                        kmax=maxmom(ktype)
                        nprimk=xnoprm(katom,ktype)
                        ncontk=xnocon(katom,ktype)
                        lenaux=nocart(ktype)
c     ----- allocate core for temporary vectors, etc. -----
c
                        alpha=1
                        ainv=alpha+3*npint*nprimk
                        xyza=ainv+npint*nprimk
                        expon=xyza+3*npint*nprimk
                        kappa=expon+npint*nprimk
                        prmint=1
                        xyz=max(kappa+npint*nprimk,
     $                          prmint+lenblk*npint*lenaux*nprimk)
                        top1=xyz+3*npint*nprimk
     $                            *(imax+1)*(jmax+1)*(kmax+1)
c
                        conint=prmint+npint*nprimk*lenblk*lenaux
c
                        tmp1=conint+nconti*ncontj*ncontk*lenblk*lenaux
                        tmp2=tmp1+nconti*nprimj
                        top2=tmp2+nconti*ncontj*nprimk 
c
                        if (top1.gt.maxcor.or.top2.gt.maxcor) then
                           call lnkerr('m352: '
     $                            //'centr3..not enough core for s')
                        end if
                        if(debug) then
                           write(iout,*) 'iat,jat,kat,ityp,jtyp,ktyp'
                           write(iout,*) iatom,jatom,katom,
     $                                   itype,jtype,ktype
                        endif
c
c     ----- form the two-dimensional integrals -----
c
                        call sint3(iatom,jatom,katom,imax,jmax,kmax,
     $                             nprimi,nprimj,nprimk,
     $                             ptprim(iatom,itype),
     $                             ptprim(jatom,jtype),
     $                             xptprm(katom,ktype),ex,exaux,c,
     $                             c,z(alpha),z(ainv),z(xyza),z(expon),
     $                             z(kappa),z(xyz),npint*nprimk,nprim,
     $                             npaux,nat,pi)
c
c     ----- form the primitive integrals -----
c
                        call form3(z(prmint),z(xyz),npint*nprimk,
     $                         lenblk*lenaux,imax,jmax,kmax,
     #                         mintyp(itype),
     $                         mintyp(itype)+nocart(itype)-1,
     #                         mintyp(jtype),
     $                         mintyp(jtype)+nocart(jtype)-1,
     #                         mintyp(ktype),
     $                         mintyp(ktype)+nocart(ktype)-1,
     #                         nx,ny,nz)
c
c     ----- transform to contracted functions -----
c
                        call trans3(z(prmint),z(conint),nprimi,nprimj,
     #                              nprimk,nconti,ncontj,ncontk,
     $                              cont(ptcont(iatom,itype)),
     #                              cont(ptcont(jatom,jtype)),
     #                              xcont(xptcon(katom,ktype)),
     #                              z(tmp1),z(tmp2),lenblk,
     $                              minmom(itype),maxmom(itype),
     #                              minmom(jtype),maxmom(jtype),
     $                              minmom(ktype),maxmom(ktype),
     $                              nocart)
c
c     ----- transfer integrals to total array -----
c
                        call put3(s,z(conint),start,xstart,
     $                            iatom,jatom,katom,
     $                            itype,jtype,ktype,
     $                            nconti,ncontj,ncontk,
     $                            nnp,lenblk,nat,nbtype,nobf)
c
   50                continue
   60             continue
   70          continue
   80       continue
   90    continue
  100 continue
c
c
      return
      end
