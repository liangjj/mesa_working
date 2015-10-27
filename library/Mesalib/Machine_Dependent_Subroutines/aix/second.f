*deck @(#)second.f	5.1  11/6/94
       function second()
       real*8 second
       external readrtc
       dimension itime(2)
       call readrtc(itime)
       second = itime(1) + 1.0d-9 * itime(2)
       return
       end
