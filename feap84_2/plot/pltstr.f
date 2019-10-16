c$Id:$
      subroutine pltstr(dt,sp,st,numnp,ndm,spflag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Replace 'third' by 'one3' from 'pconstant.h'       14/11/2006
c     2. Add modification for history projections           09/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute nodal stress values from integrals of
c               element values.

c      Inputs:
c         dt(*)      - Integral of mesh volume (diagonal weights)
c         st(numnp,*)- Integral of element values x volume
c         numnp      - Number of nodes in mesh
c         ndm        - Spatial dimension of mesh
c         spflag     - Stress projection flag: true if lumped projection

c      Outputs:
c         sp(numnp,*)- Principal values of stresses.
c                      N.B. Assumes st(*,1) = sig_11
c                      N.B. Assumes st(*,2) = sig_22
c                      N.B. Assumes st(*,3) = sig_33
c                      N.B. Assumes st(*,4) = sig_12
c                      N.B. Assumes st(*,5) = sig_23
c                      N.B. Assumes st(*,6) = sig_31
c         st(numnp,*)- Nodal values of plot quantities described
c                      by elements
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldatp.h'
      include  'pconstant.h'
      include  'pdata3.h'
      include  'prstrs.h'
      include  'strnum.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   spflag
      integer   numnp,ndm, i,ii
      real*8    dh,press
      real*8    dt(numnp),sp(numnp,*),st(numnp,*),sig(9)

      save

      if(spflag) then

        do ii = 1,numnp

          if(dt(ii).ne.0.0d0) then

c           Error estimator

            dh = 1.d0/dt(ii)
            hr(ner+ii-1) = hr(ner+ii-1)*dh

c           Stress projections

            do i = 1,npstr-1
              st(ii,i) = st(ii,i)*dh
            end do ! i

          endif ! dh > 0

        end do ! ii

c       History projection scaling

        if(histpltfl) then
          call plthis(dt,hr(np(305)),numnp)
        endif

      endif ! spflag

c     Position quantities

      do ii = 1,numnp

        if(istv.gt.0) then

c         Three-dimensional

          if(ndm.eq.3 .or. istp.eq.8) then
            sig(1)   = st(ii,1)
            sig(2)   = st(ii,2)
            sig(3)   = st(ii,3)
            sig(4)   = st(ii,4)
            sig(5)   = st(ii,5)
            sig(6)   = st(ii,6)
            call pstr3d(sig,sig(7))
            sp(ii,1) = sig(7)
            sp(ii,2) = sig(8)
            sp(ii,3) = sig(9)

c         Two-dimensional

          elseif(ndm.eq.2) then
            sig(1)   = st(ii,1)
            sig(2)   = st(ii,2)
            sig(3)   = st(ii,3)
            sig(4)   = st(ii,4)
            call pstr2d(sig,sig(7))
            sp(ii,1) = sig(7)
            sp(ii,2) = sig(8)
            sp(ii,3) = sig(3)
            sp(ii,4) = sig(9)
          endif

c         Compute mean stress and mises stress

          press    = (sp(ii,1) + sp(ii,2) + sp(ii,3))*one3
          sp(ii,5) = press
          sp(ii,6) = sqrt(1.5d0*((sp(ii,1) - press)**2
     &                         + (sp(ii,2) - press)**2
     &                         + (sp(ii,3) - press)**2))
          sp(ii,7) =        one3*(sp(ii,1) - press)**3
     &                          *(sp(ii,2) - press)**3
     &                          *(sp(ii,3) - press)**3
          if(sp(ii,7).lt.0.0d0) then
            sp(ii,7) = -(abs(sp(ii,7))**one3)
          else
            sp(ii,7) =  (abs(sp(ii,7))**one3)
          endif

        endif

      end do ! ii

      end

      subroutine plthis(dt,st,numnp)

c     Divide history projections by projection area

      implicit   none

      include   'eldatp.h'

      integer    numnp, n,i
      real*8     dh
      real*8     dt(*),st(numnp,*)

      do n = 1,numnp
        if(dt(n).ne.0.0d0) then
          dh = 1.0d0/dt(n)
          do i = 1, plhmax
            st(n,i) = st(n,i)*dh
          end do ! i
        endif
      end do ! n

      end
