*deck etimes
      subroutine etimes (tims)
      dimension tims(3)
      data cpu0 /0.0/, sys0 /0.0/, aio0 /0.0/
*
*       return elapsed cpu,sys,aio times since the last call
*
      call timing (cpu,sys,aio)
      tims(1)=cpu-cpu0
      tims(2)=sys-sys0
      tims(3)=aio-aio0
      cpu0=cpu
      sys0=sys
      aio0=aio
*
      return
      end
