*deck incmem
c***prologue           incmem
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           memory
c***author             schneider, barry (nsf)
c***source             m942
c***purpose            
c***                   
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       incmem
      subroutine incmem(rbuf,ibuf,diag,hqq,hqp,xqp,opt,tmp,ipvt,
     #                  header,lenbuf,nwks,npvec,ntri,nonzqq,nonzqp,
     #                  nodsk,need,more,maxcor)
      implicit integer (a-z)
      character*(*) header
      logical nodsk
      dimension header(10)
      common/io/inp, iout 
      call iosys('read integer '//header(10)//' from hamiltonian',
     #            1,nonzqq,0,' ')
      call iosys('read integer '//header(6)//' from hamiltonian',
     #            1,nonzqp,0,' ')
      nonzmx=max(nonzqq,nonzqp)
      nodsk=.false.
      if(nonzmx.le.lenbuf) then
         nodsk=.true.
         lenbuf=nonzmx
         write(iout,1) lenbuf
      end if
      ntri=nwks*(nwks+1)/2
      rbuf=1
      ibuf=wpadti(rbuf+lenbuf)
      diag=iadtwp(ibuf+lenbuf+lenbuf) 
      need=wpadti(diag+nwks)
      maxcor=maxcor-need
      hqq=1
      hqp=hqq+ntri
      xqp=hqp+npvec*nwks
      opt=xqp+npvec*nwks
      ipvt=wpadti(opt+max(npvec*npvec,nwks))
      tmp=iadtwp(ipvt+npvec)
      more=wpadti(tmp+ntri)
      return
 1    format(/,10x,'all of the non-zero matrix elements can be '
     #            'stored in main memory',
     #       /,10x,'set buffer = ',i8)
      end       
