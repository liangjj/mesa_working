*deck outmem
c***prologue           outmem
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
c***end prologue       outmem
      subroutine outmem(rbuf,ibuf,diag,hqp,opt,vec,hvec,resid,t,mat,
     #                  mattmp,b,btmp,list,header,lenbuf,nwks,nrhs,
     #                  nonzqq,nonzqp,nodsk,mxdvd,need,more,maxcor)
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
      if(nonzmx.le.lenbuf) then
         lenbuf=nonzmx
         nodsk=.true.
      endif
      rbuf=1
      ibuf=wpadti(rbuf+lenbuf)
      diag=iadtwp(ibuf+lenbuf+lenbuf)
      need=wpadti(diag+nwks)
      maxcor=maxcor-need
      hqp=1
      opt=hqp+nrhs*nwks
      tonow=opt+nrhs*nrhs
      maxcor = maxcor - wptoin(tonow+nwks*nrhs)
c
c     try to figure out the maximum number of vectors
c     you can get in memory.
c
      rcore=iadtwp(maxcor)
      maxvec=rcore/nwks
      maxvec = maxvec - nrhs
      mxdvd=min(maxvec,nwks)
      mxdvd=max(mxdvd,nrhs)
      vec=tonow
      hvec=vec+nwks*mxdvd
      resid=hvec+nwks*mxdvd
      t=resid+nwks*mxdvd
      mat=t+nwks*nrhs
      mattmp=mat+mxdvd*mxdvd
      b=mattmp+mxdvd*mxdvd
      btmp=b+mxdvd*nrhs
      list=wpadti(btmp+mxdvd*nrhs)
      more=list+nrhs
      write(iout,1) mxdvd
      return
 1    format(/10x,'available memory will allow a maximum of ',/,10x,
     #            i5,' davidson vectors to be held in core')
      end       
