*deck @(#)trssss.f	5.1  11/6/94
      subroutine trssss(half,temp1,temp,test,a,b,nijkl,ia,fa,ib,fb,
     #                  numkl,minkl,nbatch)
c
      implicit integer (a-z)
c
      real*8 half(*),temp1(nijkl),temp(nijkl),a(ia,fa)
      real*8 b(ib,fb)
      integer test(nijkl)
c
c     ----- start timing -----
c
c
      call rzero(temp1,nijkl)
c
c      call xpandm(nijkl,test,1,temp1,1,temp,1,junk)
      do 101 ipt=1,nbatch
         temp1(test(ipt))=temp(ipt)
  101 continue
c
      call rzero(temp,fa*ib*numkl)
c
      do 30 fapt=1,fa
         do 29 iapt=1,ia
            if (a(iapt,fapt).ne.0.0d+00) then
               call saxpy(ib*numkl,a(iapt,fapt),temp1(iapt),ia,
     #                                temp(fapt),fa)
            end if
   29    continue
   30 continue
c
      do 40 fbpt=1,fb
         do 39 ibpt=1,ib
            if (b(ibpt,fbpt).ne.0.0d+00) then
               do 38 fapt=1,fa
                  abkln=fapt+fa*(fbpt-1+fb*(minkl-1))
                  abkl=fapt+fa*(ibpt-1)
                  do 37 kl=1,numkl
                     half(abkln)=half(abkln)+temp(abkl)*b(ibpt,fbpt)
                     abkln=abkln+fa*fb
                     abkl=abkl+fa*ib
   37             continue
   38          continue
            end if
   39    continue
   40 continue
c
c     ----- end timing -----
c
c
c
      return
      end
