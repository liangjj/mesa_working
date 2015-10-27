*deck @(#)detin.f	5.1  11/6/94
      subroutine detin(joutfg,iopn,ms,ifiga,ifigb,isuma,isumb,
     1  nphase,idstor,ij,jsym,nsef,nopen,ndet,ncp)
c
c  read the array idstor from iunts1 and unpack the configurational
c  and determinantal data.
c
      common/icntrl/ iciwrt, icipun, intwrt, ihamrd, ksym
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      common /detio/ idfrst, idspc, irt, iau, inextu, intmax, ihso
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
c
      dimension joutfg(nbf,2),iopn(mxopn,2),ms(maxdet,2),
     1  ifiga(mxocc,maxdet,2),ifigb(mxocc,maxdet,2),isuma(maxdet,2),
     2  isumb(maxdet,2),nphase(maxdet,2),idstor(idspc)
c
*mdc*if cray
*      parameter (n1=64,n1m=63,n2=32,n2m=31,n4=16,n4m=15,
*mdc*else
      parameter (n1=32,n1m=31,n2=16,n2m=15,n4=8,n4m=7,
*mdc*endif
     1  n8=n1/8,n8m=n8-1,n16=n1/16,n16m=n16-1)
c
      read (iunts1) idstor
c
      jsym = idstor(1)
      nsef = idstor(2)
c
      call unpck(idstor(3),joutfg(1,ij),2,nbf)
      idim = 3 + (nbf+n2m)/n2
c
      ndet = idstor(idim)
      ncp = idstor(idim+1)
      nopen = idstor(idim+2)
      if(ihamrd.ge.1) go to 20
      idim = idim + 3
c
      if(nopen.ge.1) then
        call unpck(idstor(idim),iopn(1,ij),8,nopen)
        idim = idim + (nopen+n8m)/n8
      endif
c
      call unpck(idstor(idim),ms(1,ij),4,ndet)
      idim = idim + (ndet+n4m)/n4
c
      call unpck(idstor(idim),isuma(1,ij),16,ndet)
      idim = idim + ((ndet+n16m)/n16)
      call unpck(idstor(idim),isumb(1,ij),16,ndet)
      idim = idim + ((ndet+n16m)/n16)
c
      call unpck(idstor(idim),nphase(1,ij),1,ndet)
      idim = idim + (ndet+n1m)/n1
c
      iadd = (ndet*mxocc+n8m)/n8
      call unpck(idstor(idim),ifiga(1,1,ij),8,ndet*mxocc)
      idim = idim + iadd
      call unpck(idstor(idim),ifigb(1,1,ij),8,ndet*mxocc)
      idim = idim + iadd
c
   20 idspc = idstor(idspc)
c
      return
      end
