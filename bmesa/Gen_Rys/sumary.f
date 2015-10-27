*deck sumary.f
c***begin prologue     sumary
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           summary of polynomial information
c***author             schneider, barry (nsf)
c***source
c***purpose            
c***                   
c***                   
c***description      
c***                 
c***                                                                       
c***                                                          
c***references         
c
c***routines called    
c***end prologue       sumary

      subroutine sumary (nrun,maxdim,typwt)
      implicit integer (a-z)
      parameter ( runmax=200 )
      real*8 alf, rys
      character*(*) typwt
      character*2 itoc
      character*5 label
      character*80 title
      dimension alf(runmax), nord(runmax)
      common/io/inp, iout 
      pointer (ps,rys(1))
      wds=4*maxdim
      call memory(wptoin(nrun*wds),ps,ngot,'summary',0)
      a0=1
      b0=a0+nrun*maxdim
      pt0=b0+nrun*maxdim
      wt0=pt0+nrun*maxdim
      do 10 i=1,nrun
         call iosys('read integer "quad size-'//itoc(i)
     1              //'" from kohn',1,nord(i),0,' ')
         if(typwt.eq.'rys') then
            call iosys('read real "X-'//itoc(i)
     1                 //'" from kohn',1,alf(i),0,' ')
         endif 
         call iosys('read real "alpha-'//itoc(i)
     1              //'" from kohn',nord(i),rys(a0),0,' ')
         call iosys('read real "beta-'//itoc(i)
     1              //'" from kohn',nord(i),rys(b0),0,' ')
         call iosys('read real "points-'//itoc(i)
     1             //'" from kohn',nord(i),rys(pt0),0,' ')
         call iosys('read real "weights-'//itoc(i)
     1             //'" from kohn',nord(i),rys(wt0),0,' ')         
         a0=a0+maxdim
         b0=b0+maxdim
         pt0=pt0+maxdim
         wt0=wt0+maxdim
 10   continue 
      a0=1
      b0=a0+nrun*maxdim
      pt0=b0+nrun*maxdim
      wt0=pt0+nrun*maxdim  
      ncmx=5
      label='Rys X'
      title='          alpha recursion coefficients'
      call matprnt(title,label,alf,rys(a0),maxdim,nrun,maxdim,nrun,ncmx)
      title='          beta recursion coefficients'
      call matprnt(title,label,alf,rys(b0),maxdim,nrun,maxdim,nrun,ncmx)
      title='          points'
      call matprnt(title,label,alf,rys(pt0),maxdim,nrun,maxdim,
     1             nrun,ncmx)
      title='          weights'
      call matprnt(title,label,alf,rys(wt0),maxdim,nrun,maxdim,
     1             nrun,ncmx)
      return
      end       
