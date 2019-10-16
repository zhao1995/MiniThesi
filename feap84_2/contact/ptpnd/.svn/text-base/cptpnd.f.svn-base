c$Id:$
      subroutine cptpnd (ndm,ndf,x,u,csw,npair,cp0,ix1,ix2,
     &                   cm1,ch1,ch2,ch3,w1,w3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add ip(1) for calls to pltnod                    01/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor        September 1, 1996            1.0

c      Acronym: Contact DRIVER for ND PTP

c      Purpose: Management of specific contact formulation

c      Inputs :
c         ndm     - Space dimension of mesh
c         ndf     - Number dof/node
c         x(*)    - Nodal coordinates
c         u(*)    - Current nodal solution vectors
c         csw     - Contact switch
c         npair   - # of current pair
c         cp0(*)  - Contact pair control data
c         ix1(*)  - Element nodal connection list for surface 1
c         ix2(*)  - Element nodal connection list for surface 2
c         cm1(*)  - Material data for surface 1
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)
c         ch3(*)  - Contact history variables (static)
c         w1(*)   - Dictionary of variables for CH1 & CH2
c         w3(*)   - Dictionary of variables for CH3

c      Outputs:
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)
c         ch3(*)  - Contact history variables (static)
c                 - Perform requested activity
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'c_pair.h'
      include  'augdat.h'
      include  'cdata.h'
      include  'compas.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'print.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   ifplt,ifprt,once, setvar,palloc
      character w1(*)*(*),w3(*)*(*)
      integer   ndm,ndf, csw,npair, ix1(*), ix2(*)
      integer   i,j,k,k1,k2,kset,fel,lel, ifsurf,idir,istgn
      integer   ida(3), idx(2), ilm(3), ip(1)
      real*8    dist,penn,pent,gapn, fn,ftn,ftrn, mu, augfn,augft
      real*8    x(ndm,*),u(ndf,*), cp0(nr0,n0c3:*), cm1(*)
      real*8    ch1(lh1,*), ch2(lh1,*), ch3(lh3,*)
      real*8    res(8),tan(8,8),r(3),d(3,3),ftr(2),nn(2),kts(2)

      save

      data      ida  /1,2,3/
      data      ilm  /0,0,0/

      call cdebug0 ('    cptpnd',csw)

c-----[--.----+----.----+----.-----------------------------------------]
      if (csw.eq.0) then
        once = .true.

c     Called from PCONTR for activation of requested history variables

      elseif (csw.eq.1) then

c       Load dictionary of history variables for first call only.

c       WARNING - IF accessed only once, even if more than one contact
c                 pair uses driver.

c                 Place here all activities to bE be done only once

        if (once) then

          once = .false.

c         Define history variable names used by driver

          call defhvptp (w1,w3)

c         Print warning if unsymmetric solver is needed

          if (iffric.eq.1) then
            write (*,3000)
          endif

        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     History variables initialization             csw    :  14

      elseif (csw.eq.14) then

c       Set matched values for master contact node

        idir = max(1,min(nint(abs(cp0(4,0))),ndm))
        call setsurfs(ix1,ix2, ch3, idir)

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from FORMFE to compute stiffness and residual:   3
c     Called from FORMFE to compute only                  :   6
c     Called from PTIMPL to compute residual              : 206

      elseif (csw.eq.3 .or. csw.eq.6 .or. csw.eq.206 ) then

c       Set solution parameters

        if(iffric.eq.1) then
          mu = cm1(1)
        else
          mu = 0.0d0
        endif

c       Set material parameters

        penn = cp0(3,2)
        pent = cp0(4,2)

        idir = nint(cp0(4,0))
        k    = abs(idir)

        do i = 1,nint(cp0(5,0))
          if(ch2(p1(4),i).ne.0.0d0) then

            do k2 = 1,ndm
              r(k2) = 0.0d0
              do k1 = 1,ndm
                d(k1,k2) = 0.0d0
              end do
            end do

            augft = 0.0d0
            if(ifaugm.gt.1) then
              augfn = ch2(p1(151),i)
              if(iffric.eq.1) then
                augft = ch2(p1(152),i)
              endif
            else
              augfn = 0.0d0
            endif

            idx(1)       = nint(ch3(p3(1),i))
            idx(2)       = nint(ch3(p3(2),i))
            gapn         = ch3(p3(4),i)
            fn           = penn*gapn + augfn
            r(k)         = fn
            d(k,k)       = penn

c           Add Lagrange multipler force

            k2           = 2*ndm + 1
            do j = 1,k2
              tan(j,k2) = 0.0
              tan(k2,j) = 0.0
            end do ! i

            if( ifsolm.eq.2 ) then
              fn         =  fn + ch2(p1(21),i)
              k1         =  ndm + k
              r(k)       =  fn
              res(k2)    = -gapn
              if(idir.gt.0) then
                dist     =  1.d0
              else
                dist     = -1.d0
              endif
              tan(k2,k ) = -dist
              tan(k2,k1) =  dist
              tan(k ,k2) = -dist
              tan(k1,k2) =  dist
            endif

c           Friction case

            if(mu.gt.0.0d0) then

              k1      = mod(k ,ndm) + 1
              k2      = mod(k1,ndm) + 1
              kts(1)  = ch1(p1(54)  ,i)
              kts(2)  = ch1(p1(54)+1,i)

c             Compute trial values for friction forces

              ftr(1)  = pent*(x(k1,idx(2)) - x(k1,idx(1))
     &                +      (u(k1,idx(2)) - u(k1,idx(1)))) - kts(1)

              if(ndm.eq.3) then
                ftr(2)  = pent*(x(k2,idx(2)) - x(k2,idx(1))
     &                  +      (u(k2,idx(2)) - u(k2,idx(1)))) - kts(2)
              else
                ftr(2)  = 0.0d0
              endif

              if(idir.lt.0) then
                ftr(1) = -ftr(1)
                ftr(2) = -ftr(2)
              endif

c             Compute norms for trial and final friction values

              ftrn    = sqrt(ftr(1)**2 + ftr(2)**2)
              ftn     = mu*abs(fn)

c             Stick case

              if(ftrn .lt. ftn) then

                r(k1)    = ftr(1)
                d(k1,k1) = pent
                if(ndm.eq.3) then
                  r(k2)    = ftr(2)
                  d(k2,k2) = pent
                endif

c             Slip case

              else
                nn(1)    = ftr(1)/ftrn
                nn(2)    = ftr(2)/ftrn
                r(k1)    = ftn*nn(1)
                r(k2)    = ftn*nn(2)
                d(k1,k)  = -mu*penn*nn(1)
                d(k2,k)  = -mu*penn*nn(2)
                ftn      =  penn*ftn/ftrn
                d(k1,k1) = ftn - ftn*nn(1)*nn(1)
                d(k1,k2) =     - ftn*nn(1)*nn(2)
                d(k2,k2) = ftn - ftn*nn(2)*nn(2)
                d(k2,k1) = d(k1,k2)

c               Update slip

                kts(1)   = kts(1) + (ftrn - ftn)*nn(1)
                kts(2)   = kts(2) + (ftrn - ftn)*nn(2)

              endif
              ch2(p1(54)  ,i) = kts(1)
              ch2(p1(54)+1,i) = kts(2)
            endif

c           Fill in stiffness and residual

            do k2 = 1,ndm
              do k1 = 1,ndm
                tan(k1    ,k2    ) =  d(k1,k2)
                tan(k1+ndm,k2    ) = -d(k1,k2)
                tan(k1    ,k2+ndm) = -d(k1,k2)
                tan(k1+ndm,k2+ndm) =  d(k1,k2)
              end do
              if(idir.gt.0) then
                res(k2    ) =  r(k2)
                res(k2+ndm) = -r(k2)
              else
                res(k2    ) = -r(k2)
                res(k2+ndm) =  r(k2)
              endif
            end do

c           Assemble stiffness and residual

            if(ifsolm.eq.1) then
              ilm(1) = 0
              call constass(idx,ida,2,ndm,ilm,0,0,8,tan,res)
            else
              ilm(1) = idx(1)
              call constass(idx,ida,2,ndm,ilm,1,1,8,tan,res)
            endif

          endif
        end do

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PMACR1 for Augmented variables update   :  10

      elseif (csw.eq.10) then

        if(ifaugm.gt.1) then

c         Normal augmented update

          do i = 1,nint(cp0(5,0))
            if(ch2(p1(4),i).ne.0.0d0) then
              ch2(p1(151),i) = cp0(3,2)*ch3(p3(4),i) + ch2(p1(151),i)
              augg           = max(augg,abs(ch3(p3(4),i)))
            endif
          end do ! i

        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PMACR1 to determine contact geometry    : 103
c     Called from PMACR3 to determine contact geometry    : 304

      elseif ((csw.eq.103 .and.      ifistgn) .or.
     &        (csw.eq.304 .and. .not.ifistgn) ) then

        penn = cp0(3,2)
        idir = nint(cp0(4,0))
        k1   = abs(idir)
        do i = 1,nint(cp0(5,0))
          j  = nint(ch3(p3(1),i))
          k  = nint(ch3(p3(2),i))
          if(j.gt.0 .and. k.gt.0) then
            gapn = (x(k1,k) - x(k1,j)) + (u(k1,k) - u(k1,j))
            if(idir.lt.0) then
              gapn = - gapn
            endif
            ch3(p3(4),i) = gapn

c           Set lagrange multiplier contact state parameter

            if( ifsolm.eq.2 ) then
              if(ch2(p1(21),i).lt.0.0d0) then
                istgn = 1
              elseif(ch2(p1(21),i).gt.0.0d0) then
                istgn = 0
              elseif(gapn.le.0.0d0) then
                istgn = 1
              endif

c           Set penalty multiplier contact state parameter

            else
              if(ifaugm.gt.1) then
                fn = penn*gapn + ch2(p1(151),i)
              else
                fn = penn*gapn
              endif
              if(fn.le.0.0d0) then
                istgn = 1
              else
                istgn = 0
              endif
            endif

c           Set contact state

            if(istgn.eq.1) then
              ch2(p1(4),i) = 1.d0
              if(ifsolm.eq.2) then
                ch3(p3(21)  ,i) = max(mr(np(89)+j-1),mr(np(89)+k-1))
                ch3(p3(21)+1,i) = 0.0d0
              endif
            else
              ch2(p1(4),i) = 0.d0
              if(ifsolm.eq.2) then
                ch3(p3(21)  ,i) = 0.0d0
                ch3(p3(21)+1,i) = 0.0d0
              endif
            endif
          endif
        end do

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PMACR5 to show element information

      elseif (csw.eq.200) then

        write (*,2000)

c-----[--.----+----.----+----.-----------------------------------------]
c     Printout of contact status

      elseif(csw.eq.204) then

c      Get output flag and range

       call setcprt (ifprt,fel,lel)

c       Print title

        if (ifprt) then

          write (iow,2001) npair

c         Loop over requested surface 1 elements

        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PPLOTF for plot of contact geometry

      elseif (csw.eq.305) then

        call setcplt(ifplt,ifsurf)
        if(ifplt) then
          fp(1) = np(53) - 3
          do i = 1,nint(cp0(5,0))
            if(ifsurf.ne.2) then
              call pppcol(2,1)
              j = nint(ch3(p3(1),i))
              ip(1) = 1
              call pltnod(hr(fp(1)+3*j),ip,ndm,1,0,1,1)
            endif
            if(ifsurf.ne.1) then
              call pppcol(8,1)
              k = nint(ch3(p3(2),i))
              ip(1) = 1
              call pltnod(hr(fp(1)+3*k),ip,ndm,1,0,1,1)
            endif
          end do
        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PPLOTF to set profile and range to plot variable
c     Called from CONTACT for plot contours of a contact variable

      elseif ((csw.eq.308) .or.
     &        (csw.eq.408)    ) then

c       call c3varplt (ix1,ch1,ch2,ch3,npair,csw)

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PCONTR to activate all history variables

      elseif (csw.eq.313) then

c       Activate history variables: Compute number slave pts

        call csurfmax(ix1,dnope1,nope1,neps1, kset)
        setvar   = palloc(136,'CTEM1',kset,1)
        call cptsurf(ix1,dnope1,nope1,neps1, mr(np(136)), kset)
        setvar   = palloc(136,'CTEM1',0,1)
        cp0(5,0) = kset
        call acthvptp(kset)

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from UPDATE to update lagrange multiplers

      elseif (csw.eq.314) then

        if( ifsolm.eq.2 ) then
          do i = 1,nint(cp0(5,0))
            if (nint(ch2(p1(4),i)).ge.0) then
              ilm(1) = nint(ch3(p3(1),i))
              call getlagm(ilm(1),1,1,ch2(p1(21),i))
            else
              ch2(p1(21),i) = 0.0d0
            endif
          end do ! i
        endif


c-----[--.----+----.----+----.-----------------------------------------]
c     Reset profile for contacts

      elseif (csw.eq.403) then

c       Set degree-of-freedoms to assemble

        do i = 1,nint(cp0(5,0))
          if(ch2(p1(4),i).ne.0.0d0) then
            idx(1) = nint(ch3(p3(1),i))
            idx(2) = nint(ch3(p3(2),i))
            if(ifsolm.eq.1) then
              call modprof (idx,ida,2,ndm)
            elseif(ifsolm.eq.2) then
              ilm(1) = idx(1)
              call modprofl(idx,ida,2,ndm,ilm,1,1)
            endif
          endif
        end do

      endif

c-----[--.----+----.----+----.-----------------------------------------]

2000  format(/10x,'CPTPND: Point to Point Frictional Contact'/1x)

2001  format (/5x,'C o n t a c t   O u t p u t   f o r   P a i r ',i5)

3000  format (/'  *WARNING* Unsymmetric solver required for friction'/)

      end
