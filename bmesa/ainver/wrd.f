*deck wrd.f
c***begin prologue     wrd
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            memory outlay
c***                   
c***references         
c
c***routines called    
c***end prologue       wrd
      subroutine wrd(ham,eig,v,u,tpr,dim,words)
      implicit integer (a-z)
      common/io/inp, iout
      ibeg=words
      ham=ibeg
      eig=ham+dim*dim
      v=eig+dim
      u=v+dim*dim
      tpr=u+dim*dim
      words=tpr+dim*dim
      return
      end       
