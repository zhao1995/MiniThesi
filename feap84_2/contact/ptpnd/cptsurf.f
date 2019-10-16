c$Id:$
      subroutine cptsurf(ixc,dnope,nope,neps, icx, maxn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Determine contact nodal pairs

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    dnope,nope,neps,maxn, ixc(dnope,*), icx(*), ie,in

c     Mark location of unique node numbers

      maxn = 0
      do ie = 1,neps
        do in = 1,nope
          if(ixc(in,ie).gt.0) then
            maxn            = max(maxn,ixc(in,ie))
            icx(ixc(in,ie)) = 1
          endif
        end do ! in
      end do ! ie

c     Compress list

      in = 0
      do ie = 1,maxn
        if(icx(ie).ne.0) then
          in = in + 1
          icx(in) = ie
          if(in.ne.ie) icx(ie) = 0
        endif
      end do ! ie
      maxn = in

      end

      subroutine csurfmax(ixc,dnope,nope,neps,  maxn)

      implicit   none

      integer    dnope,nope,neps,maxn, ixc(dnope,*), ie,in

      maxn = 0
      do ie = 1,neps
        do in = 1,nope
          maxn  = max(maxn,ixc(in,ie))
        end do ! in
      end do ! ie

      end

      subroutine setsurfs(ix1,ix2, cp3, idir)

      implicit   none

      include   'c_geom.h'
      include   'sdata.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    setvar, palloc
      integer    ix1(dnope1,*),ix2(dnope2,*), idir, nmax1,nmax2
      real*8     cp3(*)

      call csurfmax(ix1,dnope1,nope1,neps1, nmax1)
      call csurfmax(ix2,dnope2,nope2,neps2, nmax2)

      setvar = palloc(136,'CTEM1',nmax1,1)
      setvar = palloc(137,'CTEM2',nmax2,1)

      call cptsurf(ix1,dnope1,nope1,neps1, mr(np(136)), nmax1)
      call cptsurf(ix2,dnope2,nope2,neps2, mr(np(137)), nmax2)

      call cptsetp(mr(np(136)),mr(np(137)), cp3, hr(np(43)),
     &             ndm, idir, nmax1,nmax2)

      setvar = palloc(137,'CTEM2',0,1)
      setvar = palloc(136,'CTEM1',0,1)

      end

      subroutine cptsetp(ix1,ix2,ch3, x, ndm,idir,neps1,neps2)

      implicit   none

      include   'c_0.h'
      include   'c_keyh.h'

      integer    ndm,idir,neps1,neps2, i,j, jset, ix1(*),ix2(*)
      real*8     ch3(lh3,*),x(ndm,*), dist,penn,pent

c     Set initial values for master contact node

      pent = 0.0d0
      do i = 1,neps1
        do j = i+1,neps1
          penn = (x(1,ix1(i)) - x(1,ix1(j)))**2
     &         + (x(2,ix1(i)) - x(2,ix1(j)))**2
          if(ndm.eq.3) then
            penn = penn + (x(3,ix1(i)) - x(3,ix1(j)))**2
          endif
          pent = max(pent,penn)
        end do ! j
      end do ! i
      if(pent.gt.0.0d0) then
        pent = 1.d-02*pent
      else
        pent = 0.01d0
      endif

c     Set matching pairs

      do i = 1,neps1
        jset = 1
        dist = (x(   1,ix1(i)) - x(   1,ix2(1)))**2
     &       + (x(   2,ix1(i)) - x(   2,ix2(1)))**2
     &       - (x(idir,ix1(i)) - x(idir,ix2(1)))**2
        if(ndm.eq.3) then
          dist = dist + (x(3,ix1(i)) - x(3,ix2(1)))**2
        endif
        do j = 2,neps2
          penn = (x(1,ix1(i)) - x(1,ix2(j)))**2
     &         + (x(2,ix1(i)) - x(2,ix2(j)))**2
     &         - (x(idir,ix1(i)) - x(idir,ix2(j)))**2
          if(ndm.eq.3) then
            penn = penn + (x(3,ix1(i)) - x(3,ix2(j)))**2
          endif
          if(penn.lt.dist) then
            jset = j
            dist = penn
          endif
        end do ! j

c       Set points which match

        if(dist.gt.pent) jset = 0
        ch3(p3(1),i) = ix1(i)
        ch3(p3(2),i) = ix2(jset)

      end do ! i

      end
