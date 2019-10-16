c$Id:$
      subroutine datrim(al,au,ad,jp,neqr,neqs,neqt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'dlog10' to 'log10'                       17/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Triangular decomposition of a matrix stored in profile
c               form

c      Notes:

c         Equations '   1  ' to 'neqr' are reduced
c         Remainder is a "reduced" tangent for export.

c         Equations '   1  ' to 'neqs' are symmetric.
c         Equations 'neqs+1' to 'neqt' are unsymmetric.

c         Use:
c          a.) All equations are unsymmetric       : neqs = 1 (or 0)
c              N.B.  The top 1 x 1 submatrix is always symmetric.
c                    Both 'al' and 'au' must be provided.

c          b.) All equations are symmetric         : neqs = neqt
c              N.B.  In this case array 'al' is not used.

c          c.) First 'neqs' equations are symmetric: 1 < neqs < neqt
c              N.B.  Storage of 'al' for unsymmetric equations only.

c      Inputs:
c         al(*)  - Unfactored lower triangular terms
c         au(*)  - Unfactored upper triangular terms
c         ad(*)  - Unfactored diagonal terms
c         jp(*)  - Pointer to row/column ends in 'al' and 'au'.
c         neqr   - Number of equations to reduce by elimination.
c         neqs   - Number of symmetric equations
c         neqt   - Number of equations in A.

c      Outputs:
c         al(*)  - Factored lower triangular terms
c         au(*)  - Factored upper triangular terms
c         ad(*)  - Factored diagonal terms
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'complx.h'
      include  'counts.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'setups.h'

      integer   neqs,neqt,neq, n0,n1,n2,n3, neq1,neq2,neqr, jp(*)
      integer   i, id, ie, ih, is, idh, ifig, j, jd, je, jh, jr, jrh
      real*4    etime, tt, tary(2)
      real*8    tol, dd, daval, dfig, dimx, dimn, dot, dsred, dured
      real*8    al(*),au(*),ad(*)

      save

c     N.B.  tol should be set to approximate half-word precision.

      data      tol/0.5d-07/

c     Factor for complex equations

      if(cplxfl) then
        is = jp(neqt) + 1
        id = neqt     + 1
        call cdatri(al,al(is),au,au(is),ad,ad(id),jp,neqs,neqt)
        return
      endif

c     Use two column elimination

      if(soltyp.eq.2 .and. neqr.ge.neqt) then
        call d4triu(al,au,ad,jp, neqs, neqt)
        return
      endif

c     Set error counters and initial values for conditioning check

      n0   = 0
      n1   = 0
      n2   = 0
      n3   = 0
      dfig = 1.0d0
      dimx = 0.0d0
      dimn = 0.0d0
      do j = 1,neqt
        dimn = max(dimn,abs(ad(j)))
      end do ! j

c     Loop through columns to perform triangular decomposition

      neq = max(1,neqs)

c     Reduce symmetric equations

      neq1 = neqr + 1
      neq2 = neqr + 2
      jd   = 1
      do j = 1,neq
        jr = jd + 1
        jd = jp(j)
        jh = jd - jr
        if(jh.gt.0) then
          is = j - jh
          ie = j - 1

c         if diagonal is zero compute a norm for singularity test

          if(ad(j).eq.0.0d0) call datest(au(jr),jh,daval)

          do i = is,min(ie,neq)
            jr = jr + 1
            id = jp(i)
            ih = min(id-jp(i-1),i-is+1)
            if(ih.gt.0) then
              jrh = jr - ih
              idh = id - ih + 1
              ih  = min(ih,neq2-is)
              au(jr) = au(jr) - dot( au(jrh), au(idh), ih )
            endif
          end do ! i
        endif

c       Reduce diagonal

        dd = ad(j)
        if(jh.ge.0) then
          jr    = jd - jh
          jrh   = j - jh - 1
          ih    = min(jh,neq1-is)
          ad(j) = ad(j) - dsred(au(jr),ad(jrh),ih+1)

c         Count errors

          if(abs(ad(j)).lt.tol*abs(dd))      n0 = n0 + 1
          if(dd.lt.0.0d0.and.ad(j).gt.0.0d0) n1 = n1 + 1
          if(dd.gt.0.0d0.and.ad(j).lt.0.0d0) n1 = n1 + 1
          if(ad(j) .eq.  0.0d0)              n2 = n2 + 1

c         Complete rank test for a zero diagonal case

          if(dd.eq.0.0d0.and.jh.gt.0) then
            if(abs(ad(j)).lt.tol*daval)      n3 = n3 + 1
          endif
        endif

c       Store reciprocal of diagonal

        if(ad(j).ne.0.0d0) then
          dimx  = max(dimx,abs(ad(j)))
          dimn  = min(dimn,abs(ad(j)))
          dfig  = max(dfig,abs(dd/ad(j)))
          if(j.le.neqr) ad(j) = 1.0d0/ad(j)
        endif

c       Output interactive equation processing

        if(ior.lt.0 .and. mod(j,1000).eq.0 ) then
          tt = etime(tary)
          write(*,3000) j,tary
        endif
      end do ! j

c     Reduce unsymmetric equations

      je  = jp(neq)
      jd  = je
      do j = neq+1,neqt
        jr = jd + 1
        jd = jp(j)
        jh = jd - jr
        if(jh.gt.0) then
          is = j - jh
          ie = j - 1

c         If diagonal is zero compute a norm for singularity test

          if(ad(j).eq.0.0d0) call datest(au(jr),jh,daval)
          do i = is,min(ie,neq)
            jr = jr + 1
            id = jp(i)
            ih = min(id-jp(i-1),i-is+1)
            if(ih.gt.0) then
              jrh = jr - ih
              idh = id - ih + 1
              ih  = min(ih,neq2-is)
              au(jr)    = au(jr   ) - dot( au(jrh   ), au(idh), ih )
              al(jr-je) = al(jr-je) - dot( al(jrh-je), au(idh), ih )
            end if
          end do ! i

          do i = max(neq+1,is),ie
            jr = jr + 1
            id = jp(i)
            ih = min(id-jp(i-1),i-is+1)
            if(ih.gt.0) then
              jrh = jr - ih
              idh = id - ih + 1
              ih  = min(ih,neq1-is+1)
              if(i.gt.neqr) then
                ih = min(ih,neq1+id-jp(i-1)-i)
              end if
              au(jr)    = au(jr)    - dot( au(jrh), al(idh-je), ih )
              al(jr-je) = al(jr-je) - dot( al(jrh-je), au(idh), ih )
            endif
          end do ! i
        endif

c       Reduce diagonal

        dd = ad(j)
        if(jh.ge.0) then
          jr    = jd - jh
          jrh   = j - jh - 1
          ad(j) = ad(j) - dured(al(jr-je),au(jr),ad(jrh),
     &                          min(jh,neq1-is)+1)

c         Count errors

          if(abs(ad(j)).lt.tol*abs(dd))      n0 = n0 + 1
          if(dd.lt.0.0d0.and.ad(j).gt.0.0d0) n1 = n1 + 1
          if(dd.gt.0.0d0.and.ad(j).lt.0.0d0) n1 = n1 + 1
          if(ad(j) .eq.  0.0d0)              n2 = n2 + 1

c         Complete rank test for a zero diagonal case

          if(dd.eq.0.0d0.and.jh.gt.0) then
            if(abs(ad(j)).lt.tol*daval)      n3 = n3 + 1
          endif
        endif

c       Store reciprocal of diagonal

        if(ad(j).ne.0.0d0) then
          dimx  = max(dimx,abs(ad(j)))
          dimn  = min(dimn,abs(ad(j)))
          dfig  = max(dfig,abs(dd/ad(j)))
          if(j.le.neqr) ad(j) = 1.0d0/ad(j)
        endif

c       Output interactive equation processing

        if(ior.lt.0 .and. mod(j,1000).eq.0 ) then
          tt = etime(tary)
          write(*,3000) j,tary
        endif
      end do ! j

c     Print error summary and conditioning information

      dd = 0.0d0
      if(dimn.ne.0.0d0) dd = dimx/dimn
      ifig = nint(log10(dfig) + 0.6d0)
      if(n0.gt.0) write(iow,2000) n0,nstep,niter
      if(n1.gt.0) write(iow,2001) n1,nstep,niter
      if(n2.gt.0) write(iow,2002) n2,nstep,niter
      if(n3.gt.0) write(iow,2003) n3,nstep,niter
      if(pfr) write(iow,2004) dimx,dimn,dd,ifig
      if(ior.lt.0) then
        if(n0.gt.0) write(*,2000) n0,nstep,niter
        if(n1.gt.0) write(*,2001) n1,nstep,niter
        if(n2.gt.0) write(*,2002) n2,nstep,niter
        if(n3.gt.0) write(*,2003) n3,nstep,niter
        if(pfr) write(*,2004) dimx,dimn,dd,ifig
      endif

c     Formats

2000  format(' *WARNING* Lost at least 7 digits in',
     &       ' reducing diagonals of',i5,' equations.'/
     &       ' Step:',i7,' Iteration:',i6)

2001  format(' *WARNING* sign of diagonal changed when reducing',
     &    i5,' equations.'/' Step:',i7,' Iteration:',i6)

2002  format(' *WARNING* reduced diagonal is zero in',
     &    i5,' equations.'/' Step:',i7,' Iteration:',i6)

2003  format(' *WARNING* rank failure for zero unreduced diagonal in',
     &    i5,' equations.'/' Step:',i7,' Iteration:',i6)

2004  format(' Condition check: D-max',e11.4,'; D-min',e11.4,
     & '; Ratio',e11.4/' Maximum no. diagonal digits lost:',i3)

3000  format(' ---> Equation = ',i5,' , Time: CPU = ',f12.2,
     &       ' , System = ',f12.2)

      end


