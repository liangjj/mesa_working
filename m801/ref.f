*deck @(#)ref.f	5.1  11/6/94
      subroutine ref(arc,wght,levnr,levpt,orbtbf,bfcode,bfkey,
     #csav,rowsv,wtsav,refwt,mxref)
c
c
c
      implicit integer (a-z)
      character*1 multrf,valenc,bfkey
      character*3 codes,words*18
c
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     #,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     #,               valvir,opensh,multi,speshl
      common /drtchr/ codes(9),words(9),multrf,valenc
      common /dimens/ nbf,nsym,norbs,nrowsp,nrows4p,nrows,nrows4
     #,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     #,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc,nvref
     #,               nijvir
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     #,               levocc,spec,sspesh,val
      common /tapes/  out,errout,input,drttap
      common /cases/  casev(9)
c
      dimension orbtbf(norbs),bfcode(nrefs,nbf),bfkey(nbf)
      dimension arc(4,nrows)
      dimension levnr(nlevs),levpt(nlevs)
      dimension wght(4,nrows)
      dimension rowsv(nlevs),csav(nlevs),wtsav(nlevs),refwt(mxref)
c
      if (levval.le.10000) go to 110
      do 100 refer=1,nrefs
         wt=1
         row=1
         do 50 orb=norbs,1,-1
            lev=orb+1
            case=casev(bfcode(refer,orbtbf(orb)))
            wt=wt+wght(case,levpt(lev)+row)
            row=arc(case,levpt(lev)+row)
   50    continue
c
         refwt(refer)=wt
         write (out,51) refer,wt
   51    format ('  reference #',i4,' is walk',i8)
  100 continue
      nvref=nrefs
c
      return
c
  110 continue
c
      row=1
      do 150 orb=norbs,levval-1,-1
         lev=orb+1
         row=arc(casev(bfcode(1,orbtbf(orb))),levpt(orb+1)+row)
  150 continue
c
      write (out,151) orb+1,row
  151 format (' levval pt',3i5)
c
      bottom=lev-1
      rowm1=row
      rowsv(bottom)=rowm1
      nvref=0
c
      wt=1
      do 155 i=levval-2,1,-1
         case=casev(bfcode(1,orbtbf(i)))
         wt=wt+wght(case,levpt(i+1)+row)
         row=arc(case,levpt(i+1)+row)
  155 continue
c
      minc=1
  159 continue
      min=levpt(lev)+1
      max=levpt(lev)+levnr(lev)
  160 continue
      do 165 row=min,max
      do 164 cs=minc,4
         if (arc(cs,row).eq.rowm1) go to 170
  164 continue
      minc=1
  165 continue
  166 continue
      lev=lev-1
      wt=wtsav(lev)
      if (lev.le.bottom) return
      rowm1=rowsv(lev-1)
      min=levpt(lev)+rowsv(lev)
      max=levpt(lev)+levnr(lev)
      minc=csav(lev)+1
      go to 160
c
  170 continue
      csav(lev)=cs
      rowsv(lev)=row-levpt(lev)
      wtsav(lev)=wt
      wt=wt+wght(cs,row)
      if (lev.ge.nlevs) go to 200
      rowm1=row-levpt(lev)
      lev=lev+1
      go to 159
c
  200 continue
      nvref=nvref+1
      refwt(nvref)=wt
      write (out,201) nvref,wt
  201 format (' valence reference #',i5,' is walk ',i8)
      go to 166
c
      end
