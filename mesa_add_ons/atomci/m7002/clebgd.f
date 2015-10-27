*deck clebgd
      subroutine clebgd (ia,ib,ic,id,ie,if,kselc,rac)
      implicit real *8 (a-h,o-z)
c  calculates clebsch-gordan coeffs.  parameters ia thru if are 2j1, 2j2
c  2j3, 2m1, 2m2, 2m3, resp.  rac is the value returned.
c  kselc=1 checks the selection rules.  kselc=0 does not.
c  t tamura, comp phys comm 1 (1970) 337.
c  to increase size of angular momenta that can be handled, just increas
c  the dimension on faclog and the initial do loop which calculates it.
      common /racah/ faclog(500), check
      data check /0.0d0/
c  construct logs of factorials the first time thru.
      if (check.ne.0.0d0) go to 20
      faclog(1)=0.0d0
      faclog(2)=0.0d0
      fn=1.0d0
      do 10 n=3,500
      fn=fn+1.0d0
   10 faclog(n)=faclog(n-1)+log(fn)
      check=1.0d0
   20 continue
      rac=0.0
      ig4=ia+ib+ic
      if (kselc.eq.0) go to 30
      if (id+ie.ne.if) go to 90
      if (mod(ig4,2).ne.0) go to 90
      if (ia+ib-ic.lt.0.or.ic-iabs(ia-ib).lt.0) go to 90
      if (min0(ia-iabs(id),ib-iabs(ie),ic-iabs(if)).lt.0) go to 90
      if (mod(ib+ie,2).ne.0.or.mod(ic+if,2).ne.0) go to 90
   30 if (ia.eq.0.or.ib.eq.0) go to 40
      if (ic) 90,50,60
   40 rac=1.0
      go to 90
   50 fb=ib+1
      rac=((-1.0d0)**((ia-id)/2))/sqrt(fb)
      go to 90
   60 if (id.ne.0.or.ie.ne.0) go to 70
      ig2=ig4/2
      if (mod(ig2,2).ne.0) go to 90
      ig=ig2/2
      i1=(ia+ib-ic)/2+1
      i2=(ic+ia-ib)/2+1
      i3=(ic+ib-ia)/2+1
      i4=ig2+2
      i5=ig+1
      i6=(ig2-ia)/2+1
      i7=(ig2-ib)/2+1
      i8=(ig2-ic)/2+1
      f1=exp(.5d0*(faclog(i1)+faclog(i2)+faclog(i3)-faclog(i4))+
     1            (faclog(i5)-faclog(i6)-faclog(i7)-faclog(i8)))
      f2=ic+1
      f2=sqrt(f2)
      s1=1-2*mod((ig2+ic)/2,2)
      rac=s1*f1*f2
      go to 90
   70 fc2=ic+1
      iabcp=ig4/2+1
      iabc=iabcp-ic
      icab=iabcp-ib
      ibca=iabcp-ia
      iapd=(ia+id)/2+1
      iamd=iapd-id
      ibpe=(ib+ie)/2+1
      ibme=ibpe-ie
      icpf=(ic+if)/2+1
      icmf=icpf-if
      sqfclg=0.5d0*(log(fc2)-faclog(iabcp+1)+faclog(iabc)+faclog(icab)
     1 +faclog(ibca)+faclog(iapd)+faclog(iamd)+faclog(ibpe)+faclog(ibme)
     2 +faclog(icpf)+faclog(icmf))
      nzmic2=(ib-ic-id)/2
      nzmic3=(ia-ic+ie)/2
      nzmi=max0(0,nzmic2,nzmic3)+1
      nzmx=min0(iabc,iamd,ibpe)
      s1=(-1.0d0)**(nzmi-1)
      do 80 nz=nzmi,nzmx
      nzm1=nz-1
      nzt1=iabc-nzm1
      nzt2=iamd-nzm1
      nzt3=ibpe-nzm1
      nzt4=nz-nzmic2
      nzt5=nz-nzmic3
      termlg=sqfclg-faclog(nz)-faclog(nzt1)-faclog(nzt2)-faclog(nzt3)
     1 -faclog(nzt4)-faclog(nzt5)
      ssterm=s1*exp(termlg)
      rac=rac+ssterm
   80 s1=-s1
      if (abs(rac).lt.0.000001d0) rac=0.0d0
   90 return
      end
