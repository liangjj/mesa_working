*deck @(#)indexf.f	5.1  11/6/94
      function indexf(ii,jj,kk,ll)
c
c  compute the canonical position specified by (ii,jj,kk,ll)
c
      ij     = (((max(ii,jj)) - 1)*(max(ii,jj)))/2 + min(ii,jj)
      kl     = (((max(kk,ll)) - 1)*(max(kk,ll)))/2 + min(kk,ll)
      indexf = (((max(ij,kl)) - 1)*(max(ij,kl)))/2 + min(ij,kl)
      return
      end
