*deck @(#)detout.f	5.1  11/6/94
      subroutine detout(joutfg,iopn,ms,ifiga,ifigb,isuma,isumb,
     1  nphase,idstor,ispc)
c
c  pack the configurational data into the array idstor and write the
c  array onto iunts1.  the last word of idstor contains the number
c  of elements in the next array.
c
      common/icntrl/ iciwrt, icipun, intwrt, ihamrd, ksym
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      common /detio/ idfrst, idspc, irt, iau, inextu, intmax, ihso
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /c2/ msbas,jsym,nopen,kdbl,ndeti,nsefi,ncpi
c
      dimension joutfg(*),iopn(*),ms(*),ifiga(mxocc,*),
     1  ifigb(mxocc,*),isuma(*),isumb(*),nphase(*),idstor(idspc)
c
*mdc*if cray
*      parameter (n1=64,n1m=63,n2=32,n2m=31,n4=16,n4m=15,
*mdc*else
      parameter (n1=32,n1m=31,n2=16,n2m=15,n4=8,n4m=7,
*mdc*endif
     1  n8=n1/8,n8m=n8-1,n16=n1/16,n16m=n16-1)
c
      if(ihamrd.ge.1) then
        inext = 6 + (nbf+n2m)/n2
      else
        inext=6+(nbf+n2m)/n2+(nopen+n8m)/n8+(ndeti+n4m)/n4+2*
     1    ((ndeti+n16m)/n16)+(ndeti+n1m)/n1+2*((ndeti*mxocc+n8m)/n8)
      endif
c
      if(inext.gt.inextu) then
        iau = iau + (inext - inextu + irt -1)/irt
        write (iw,20) iau
   20     format(/'more space needed in detout.  iaup at least',i8)
        stop
      endif
c
      if(ispc.ne.1) then
        idstor(idspc) = inext
        write (iunts1) idstor
      endif
c
      idstor(1) = jsym
      idstor(2) = nsefi
c
      call pck(idstor(3),joutfg,2,nbf)
      idim = 3 + (nbf+n2m)/n2
c
      idstor(idim) = ndeti
      idstor(idim+1) = ncpi
      idstor(idim+2) = nopen
      idim = idim + 3
c
      if(ihamrd.ge.1) go to 50
c
      if(nopen.ne.0) then
        call pck(idstor(idim),iopn,8,nopen)
        idim = idim + (nopen+n8m)/n8
      endif
c
      call pck(idstor(idim),ms,4,ndeti)
      idim = idim + (ndeti+n4m)/n4
c
      call pck(idstor(idim),isuma,16,ndeti)
      idim = idim + (ndeti+n16m)/n16
      call pck(idstor(idim),isumb,16,ndeti)
      idim = idim + (ndeti+n16m)/n16
c
      call pck(idstor(idim),nphase,1,ndeti)
      idim = idim + (ndeti+n1m)/n1
c
      iadd = (ndeti*mxocc+n8m)/n8
      call pck(idstor(idim),ifiga,8,ndeti*mxocc)
      idim = idim + iadd
      call pck(idstor(idim),ifigb,8,ndeti*mxocc)
      idim = idim + iadd
c
   50 idspc = idim
c
      return
c
      end
