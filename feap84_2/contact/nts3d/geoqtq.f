c$Id:$
      subroutine geoqtq (kset,x,u,ix1,ix2,ns,xs,ch1,ch2,ch3,ids,
     &                   surpoin,insegt,conseg,cpoint,xipatch,mue,csw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'conseg' and 'cpoint' to integer arrayw   01/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Anna Haraldsson             February 1998            1.0

c      Acronym:

c      Purpose: Search for neighbour to edge defined by node knoten
c               sucht nachbar zu kante entsprechend Knoten

c      Inputs:
c         kset    - Segment number
c         x(*)    - Nodal coordinate array
c         u(*)    - Nodal displacement array
c         ix1(*)  - 1-Surface facet node connections
c         ix2(*)  - 2-Surface facet node connections
c         ns      - Slave node number
c         xs(*)   - Slave node position at t_n+1
c         ch1(*)  - History parameters at t_n
c         ch2(*)  - History parameters at t_n+1
c         ch3(*)  - History parameters
c         ids(*)  - Equation/BC indicator for slave node
c         surpoin - Surface points
c         insegt  - In segment flag
c         conseg  - Contact segment flag ( 0 = false, 1 = true)
c         cpoint  - Point flag           ( 0 = false, 1 = true)
c         xipatch - Surface coordinates
c         mue     - Friction flag
c         csw     - Solution option switch (csw =  14: History init.;
c                                               = 304: Geometric state.

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
c              4           3
c               o---------o
c               |    xi2  |
c               |    |    |
c               |    +-xi1|
c               |         |
c               |         |
c               o---------o
c              1           2

c                    s                istgn = 0
c               o------------o
c                ////////////

c                                    istgn = 1
c               o------------o
c                //// s /////

c                   o----------o
c                  / s  ///////       istgn = 2
c                 /
c                /                                   a1
c               o                                o------->o
c                                               nr     nrneben

c               o                    istgn = 3
c                \
c                 \
c                  \
c                   o----------o
c                     /////////
c                 s

c               o                    istgn = 4: unter Knoten
c                \
c                 \
c                  \
c                   @----------o
c                  /                                  o
c                 /                                  nr
c                o
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_contac.h'
      include   'c_geom.h'
      include   'c_keyh.h'
      include   'c_pair.h'
      include   'c_tole.h'
      include   'sdata.h'
      include   'pointer.h'
      include   'comblk.h'

      include   'p_int.h'

      logical    flag2,flag,ueberall
      integer    csw,kset,nel2,nod2,masts,istgn,istgt,mastsold,inseg
      integer    nr,zahl,node,nodes,anzseg,seg,anzkn, ns
      integer    nrneben,nachb,seg2,knoten,nel1,zweiter,istgnold
      integer    iter,imax,nrold,i,j,nocont
      integer    ix1(dnope1,*),ix2(dnope2,*),im(4),surpoin(*)
      integer    insegt(*),ids(3),conseg(*),cpoint(*)
      real*8     aslave,a11,a12,a22,det,dist,gapp,gapn,gn,nbar,maxxi
      real*8     temp1,temp2,xi1,xi2,mue,aij(2,2),a1(3),a2(3), g(3)
      real*8     ni(4),niab(4,2),niabab(4,2,2),xnew(3),nvec(3)
      real*8     t1(3),t2(3),n(3),xp(3),xm(3,4),xi(2),xold(3)
      real*8     x(ndm,*),u(ndf,*),xs(3),ch1(*),ch2(*),ch3(*)
      real*8     xipatch(2,*)

      save

      data       imax / 25 /

      call cdebug0 ('    geoqtq',-1)

c-----[--.----+----.----+----.-----------------------------------------]

      if(csw.eq.14) return  ! HACK HACK

      call as(kset,x,u,ix1,surpoin,aslave, nodes)

      ch3(p3(16)) = aslave

      gapp  = 1.0d8
      istgn = 0
      inseg = 1

c     Search closest master-point

      fp(2) = np(191)+ surpoin(nsurf2)
      anzkn = mr(fp(2))

c     Loop over all surface 2 nodes: Coordinate master = xnew(.)

      nocont = 0
      do nod2 = 1,anzkn
         do j = 1,3
           xnew(j) = x(j,mr(fp(2)+nod2)) + u(j,mr(fp(2)+nod2))
           g(j)    = xnew(j) - xs(j)
         end do ! j

         dist = sqrt(g(1)*g(1) + g(2)*g(2) + g(3)*g(3))

         if (dist.le.gapp) then
           gapp = dist
           nr   = mr(fp(2)+nod2)
           node = nod2
         endif

      end do ! nod2

      mastsold    = nint(ch2(p1( 1)))
      nrold       = nint(ch2(p1(25)))
      ch2(p1(25)) = nr

c-----[--.----+----.----+----.-----------------------------------------]

      gapp   = 1.0d8
      anzseg = mr(np(191)+mr(fp(2)+node+anzkn))

      zahl = 0

c-----[--.----+----.----+----.-----------------------------------------]

c     Loop over elements of patch of closest node

      maxxi  =  10.0d0
      do seg = 1,anzseg

        conseg(seg)    =  0
        xipatch(1,seg) =  0.0d0
        xipatch(2,seg) =  0.0d0

        nel2 = mr(np(191)+mr(fp(2)+node+anzkn)+seg)

        do nod2 = 1,4
          im(nod2) = ix2(nod2,nel2)
          do j = 1,3
            xm(j,nod2) = x(j,im(nod2)) + u(j,im(nod2))
          end do ! j
        end do ! nod2

c       Make sure that slave node is not on facet

        do nod2 = 1,4
          if(im(nod2).eq.ns) then
            nocont = nocont + 1
            go to 110
          endif
        end do ! nod2

c       Perform constrained closest point projection

        call cnproj( ids, xs,xm, xi,xp,t1,t2,g, iter)

        if(iter.ge.imax) then
          if(ifdb) then
            write(*,*) ' No convergence in GEOQTQ, slave node =', nodes
          end if
          nocont = nocont + 1
          go to 110
        end if

c       Compute normal vector to facet (times twice area)

        n(1) = t1(2)*t2(3) - t1(3)*t2(2)
        n(2) = t1(3)*t2(1) - t1(1)*t2(3)
        n(3) = t1(1)*t2(2) - t1(2)*t2(1)
        gapn = g(1)*n(1) + g(2)*n(2) + g(3)*n(3)
        nbar = 1.0d0/(n(1)*n(1) + n(2)*n(2) + n(3)*n(3))

c       If gap negative check postion

        if( gapn*sqrt(nbar) .le. tlopen ) then

c         CONTACT!!!

          conseg(seg)= 1

          aij(1,1) = t1(1)*t1(1) + t1(2)*t1(2) + t1(3)*t1(3)
          aij(1,2) = t1(1)*t2(1) + t1(2)*t2(2) + t1(3)*t2(3)
          aij(2,2) = t2(1)*t2(1) + t2(2)*t2(2) + t2(3)*t2(3)

c         Contact point in Segment

c         if (max(abs(xi(1)),abs(xi(2))).le.1.d0+tlouts) then
          if (max(abs(xi(1)),abs(xi(2))).le.maxxi) then
            maxxi = max(abs(xi(1)),abs(xi(2)))
            nbar =  sqrt(nbar)
            gapn = -gapn*nbar

c           Found contact with a segment before => istgn=2

            if (istgn.eq.1) then
              inseg = inseg + 1
c             istgn = 2
            endif
            insegt(inseg) =  seg
            gn            = -abs(gapn)
            masts         =  nel2
            a11           =  aij(1,1)
            a12           =  aij(1,2)
            a22           =  aij(2,2)
            xi1           =  xi(1)
            xi2           =  xi(2)
            gapp          =  abs(gapn)
            do j = 1,3
              nvec(j) = nbar * n(j)
              a1(j)   = t1(j)
              a2(j)   = t2(j)
              xnew(j) = xp(j)
            end do ! j
            xipatch(1,seg) = xi(1)
            xipatch(2,seg) = xi(2)
            istgn          = 1

c         Contact point not in this segment

          else

c           if no contact to segment found before, save xi & segment
c           to find out if contact to edge or contact to point
c           (and to which edge)

            if((istgn.ne.1).and.(istgn.ne.2)) then
              nbar           =  sqrt(nbar)
              gapn           = -gapn*nbar
              xipatch(1,seg) =  xi(1)
              xipatch(2,seg) =  xi(2)
              istgn          =  3
            endif
          endif
        end if
110     continue
      end do ! seg   Loop over elements of patch of closest node

c-----[--.----+----.----+----.-----------------------------------------]
       if(nocont.eq.anzkn) return
c-----[--.----+----.----+----.-----------------------------------------]

c     Check that contact is within a facet

      if (istgn.eq.1) then

        if(maxxi.gt.1.d0+tlouts) then
          istgn = 0
        endif

c     Contact point in more than 1 Segment

      elseif (istgn.eq.2) then

c       Contact point in more than 2 Segments => contact to point

        if (inseg.ge.3) then

          istgn = 4

c       Contact point in 2 Segments => contact to edge

        else

          if (inseg .ne. 2) then
            write(*,*) 'Somethings wrong in geoqtq: inseg=',inseg
          endif
          seg  = insegt(2)
          nel2 = mr(np(191)+mr(np(191)+surpoin(nsurf2)+node+anzkn)+seg)
          do nod2 = 1,4
            im(nod2) = ix2(nod2,nel2)
          end do ! nod2
          seg  = insegt(1)
          nel1 = mr(np(191)+mr(np(191)+surpoin(nsurf2)+node+anzkn)+seg)
          if ((im(1).eq.nr).or.(im(3).eq.nr)) then
            do nod2 = 1,4
              if ((im(2).eq.ix2(nod2,nel1)).or.
     &            (im(4).eq.ix2(nod2,nel1)))
     &        nrneben = ix2(nod2,nel1)
            end do ! nod2
          elseif ((im(2).eq.nr).or.(im(4).eq.nr)) then
            do nod2 = 1,4
              if ((im(1).eq.ix2(nod2,nel1)).or.
     &            (im(3).eq.ix2(nod2,nel1)))
     &          nrneben = ix2(nod2,nel1)
            end do ! nod2
          endif
        endif

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Out-of-segment-case

      elseif (istgn.eq.3) then

        seg  = 1
        zahl = 0

c       Initialize point array

        do i = 1, anzseg
          cpoint(i) = 0
        end do ! i

        fp(2) = np(191)+surpoin(nsurf2)

 111    if(conseg(seg).eq.1) then
          call suchkant(knoten,flag,ix2,xipatch(1,seg),xipatch(2,seg),
     &                  mr(np(191)+mr(fp(2)+node+anzkn)+seg),nr,zweiter)
          if(.not.flag) then
            call suchnach(seg,knoten,nachb,node,ix2,surpoin)
            if(conseg(nachb).eq.1) then
              zahl = zahl + 1
              call suchkant(knoten,flag2,ix2,
     &                      xipatch(1,nachb),xipatch(2,nachb),
     &                      mr(np(191)+mr(fp(2)+node+anzkn)+nachb),nr,
     &                      zweiter)
              if(.not.flag2) then
                call suchnach(nachb,knoten,seg2,node,ix2,surpoin)
                if (seg2.eq.seg) then
                  goto 400
                endif
              endif
            endif
          else
c           cpoint(i) =  1
            cpoint(seg) = 1
          endif
        endif

        if (zahl.ge.anzseg) then
          ueberall = .true.
          do i = 1,anzseg
            if(cpoint(i).eq.0) ueberall = .false.
          end do ! i

c         Found right node

          if (ueberall) then

            istgn = 4

c         No contact!

          else

            istgn = 0

          endif

          goto 500

        endif
        seg  = seg  + 1
        zahl = zahl + 1
        goto 111
 400    nrneben = knoten
 500    continue
      endif

c-----[--.----+----.----+----.-----------------------------------------]

c     Contact to edge: Compute xi, nvec, a, gap

      if ((istgn.eq.2).or.(istgn.eq.3)) then

        do j = 1,3
          xm(j,1) = x(j,nr)      + u(j,nr)
          xm(j,2) = x(j,nrneben) + u(j,nrneben)
          a1(j)    = xm(j,2) - xm(j,1)
          a2(j)    = xm(j,1) - xs(j)
        end do ! j

c       Projection of slave node to edge

        a11 =   a1(1)*a1(1) + a1(2)*a1(2) + a1(3)*a1(3)
        a11 = -(a1(1)*a2(1) + a1(2)*a2(2) + a1(3)*a2(3))/a11
        do i=1,3
          nvec(i) = a2(i)   + xi1*a1(i)
          xnew(i) = xm(1,j) + xi1*a1(j)
        end do ! i
        gn = sqrt(nvec(1)**2 + nvec(2)**2 + nvec(3)**2)
        do i = 1,3
          nvec(i) = nvec(i)/gn
        end do ! i
        gn = -gn

c     Contact to point: Compute nvec, gap

      elseif(istgn.eq.4) then

        do j = 1,3
          xnew(j) = x(j,nr) + u(j,nr)
          nvec(j) = xs(j)   - xnew(j)
        end do ! j

        gn = sqrt(nvec(1)**2 + nvec(2)**2 + nvec(3)**2)

        do i = 1,3
          nvec(j) = -nvec(j)/gn
        end do ! i

        gn = -gn

      endif

c-----[--.----+----.----+----.-----------------------------------------]

      istgt      = 0
      ch2(p1(4)) = istgn

c     NO CONTACT

      if(istgn.eq.0) then

        if(mue.gt.0.d0) then
          ch2(p1(3))    = 0.d0
          ch2(p1(10))   = 0.d0
          ch2(p1(11))   = 0.d0
          ch2(p1(26)  ) = 0.d0
          ch2(p1(26)+1) = 0.d0
          ch2(p1(26)+2) = 0.d0
        endif

        ch2(p1( 1)  ) = masts
        ch3(p3( 9)  ) = gapn
        ch3(p3(14))   = 0.d0
        ch2(p1(25)  ) = nr

        if(ifdb .and. masts.ne.mastsold) then
          write(99,3000) istgn,kset,masts,mastsold
        endif
        return

c     NORMAL GAP => NORMAL CONTACT

      elseif (istgn.eq.1) then

        ch2(p1( 1)  ) = masts
        ch3(p3( 9)  ) = gn
        ch3(p3(18)  ) = nvec(1)
        ch3(p3(18)+1) = nvec(2)
        ch3(p3(18)+2) = nvec(3)
        ch2(p1(19)  ) = a1(1)
        ch2(p1(19)+1) = a1(2)
        ch2(p1(19)+2) = a1(3)
        ch3(p3(20)  ) = a2(1)
        ch3(p3(20)+1) = a2(2)
        ch3(p3(20)+2) = a2(3)
        ch2(p1(24)  ) = xi1
        ch2(p1(24)+1) = xi2
        ch2(p1(25)  ) = nr

        if(ifdb .and. masts.ne.mastsold) then
          write(99,3000) istgn,kset,masts,mastsold
        endif

        if(mue.gt.0.d0) then
          ch2(p1(3))    = 0.d0
          ch2(p1(26)  ) = 0.d0
          ch2(p1(26)+1) = 0.d0
          ch2(p1(26)+2) = 0.d0
        endif

      elseif((istgn.eq.2).or.(istgn.eq.3))then

        ch2(p1( 1)  ) = nr
        ch3(p3( 9)  ) = gn
        ch3(p3(18)  ) = nvec(1)
        ch3(p3(18)+1) = nvec(2)
        ch3(p3(18)+2) = nvec(3)
        ch2(p1(19)  ) = a1(1)
        ch2(p1(19)+1) = a1(2)
        ch2(p1(19)+2) = a1(3)
        ch2(p1(24)  ) = xi1
        ch2(p1(25)  ) = nrneben

        if(ifdb .and. nr.ne.nrold) then
          write(99,3000) istgn,kset,nr,nrold
        endif

        if(mue.gt.0.d0) then
          ch2(p1(3))    = 0.d0
          ch2(p1(26)  ) = 0.d0
          ch2(p1(26)+1) = 0.d0
          ch2(p1(26)+2) = 0.d0
        endif

      elseif(istgn.eq.4)then

        ch2(p1( 1)  ) = nr
        ch3(p3( 9)  ) = gn
        ch3(p3(18)  ) = nvec(1)
        ch3(p3(18)+1) = nvec(2)
        ch3(p3(18)+2) = nvec(3)

        if(ifdb .and. nr.ne.nrold) then
          write(99,3000) istgn,kset,nr,nrold
        endif

        if(mue.gt.0.d0) then
          ch2(p1(3))    = 0.d0
          ch2(p1(26)  ) = 0.d0
          ch2(p1(26)+1) = 0.d0
          ch2(p1(26)+2) = 0.d0
        endif

      else

      endif

c     Friction checks

      if (mue.gt.0.d0) then

        istgnold = nint(ch1(p1(4)))
        if(istgnold.eq.1) then
          mastsold = nint(ch1(p1(1)))
          do nod2 = 1,4
            im(nod2) = ix2(nod2,mastsold)
            do j=1,3
              xm(j,nod2) = x(j,im(nod2)) + u(j,im(nod2))
            end do ! j
          end do ! nod2

          xi(1) = ch1(p1(24)  )
          xi(2) = ch1(p1(24)+1)

          call formfkt(xi,ni,niab,niabab)

          do j = 1,3
            xold(j) = ni(1)*xm(j,1) + ni(2)*xm(j,2)
     &              + ni(3)*xm(j,3) + ni(4)*xm(j,4)
          end do ! j
        elseif((istgnold.eq.2).or.(istgnold.eq.3)) then
          nrold      = nint(ch1(p1(1)))
          do j = 1,3
            xm(1,j) = x(j,nrold) + u(j,nrold)
          end do ! j
          xi1 = ch1(p1(24)  )
          do j = 1,3
            xold(j) = xm(1,j) + xi1*ch1(p1(19)+j-1)
          end do ! j
        elseif(istgnold.eq.4) then
          nrold = nint(ch1(p1(1)))
          do j = 1,3
            xold(j) = x(j,nrold) + u(j,nrold)
          end do ! j
        endif

        do i = 1,3
          ch2(p1(26)+i-1) = ch1(p1(26)+i-1) + xnew(i) - xold(i)
        end do ! i

        if (istgn.eq.1) then

          if (istgnold.eq.1) then

c           TANGENTIAL GAP => TANGENTIAL CONTACT

            if((xi1.ne.ch1(p1(24))).and.(xi2.ne.ch1(p1(24)+1)))then

              istgt = 1
              if (masts.ne.ch1(p1(1))) then

                det         =  a11*a22 - a12*a12
                temp1       =  xnew(1)*a1(1)
     &                      +  xnew(2)*a1(2)
     &                      +  xnew(3)*a1(3)
     &                      -  xold(1)*a1(1)
     &                      -  xold(2)*a1(2)
     &                      -  xold(3)*a1(3)
                temp2       =  xnew(1)*a2(1)
     &                      +  xnew(2)*a2(2)
     &                      +  xnew(3)*a2(3)
     &                      -  xold(1)*a2(1)
     &                      -  xold(2)*a2(2)
     &                      -  xold(3)*a2(3)
                ch2(p1(10)) = (a22*temp1 - a12*temp2)/det
                ch2(p1(11)) = (a11*temp2 + a12*temp1)/det

c             Same master segment as in last time step

              else
                ch2(p1(10)) = xi1 - ch1(p1(24))
                ch2(p1(11)) = xi2 - ch1(p1(24)+1)
              endif
            endif

            ch2(p1(3)) = istgt

          elseif((istgnold.eq.2).or.(istgnold.eq.3)) then

            istgt       = 1
            det         =  a11*a22 - a12*a12
            temp1       =  xnew(1)*a1(1)
     &                  +  xnew(2)*a1(2)
     &                  +  xnew(3)*a1(3)
     &                  -  xold(1)*a1(1)
     &                  -  xold(2)*a1(2)
     &                  -  xold(3)*a1(3)
            temp2       =  xnew(1)*a2(1)
     &                  +  xnew(2)*a2(2)
     &                  +  xnew(3)*a2(3)
     &                  -  xold(1)*a2(1)
     &                  -  xold(2)*a2(2)
     &                  -  xold(3)*a2(3)
            ch2(p1(10)) = (a22*temp1 - a12*temp2)/det
            ch2(p1(11)) = (a11*temp2 + a12*temp1)/det
            ch2(p1(3))  = istgt

          elseif(istgnold.eq.4) then

            istgt       = 1
            det         =  a11*a22 - a12*a12
            temp1       =  xnew(1)*a1(1)
     &                  +  xnew(2)*a1(2)
     &                  +  xnew(3)*a1(3)
     &                  -  xold(1)*a1(1)
     &                  -  xold(2)*a1(2)
     &                  -  xold(3)*a1(3)
            temp2       =  xnew(1)*a2(1)
     &                  +  xnew(2)*a2(2)
     &                  +  xnew(3)*a2(3)
     &                  -  xold(1)*a2(1)
     &                  -  xold(2)*a2(2)
     &                  -  xold(3)*a2(3)
            ch2(p1(10)) = (a22*temp1 - a12*temp2)/det
            ch2(p1(11)) = (a11*temp2 + a12*temp1)/det
            ch2(p1(3))  = istgt

          elseif(istgnold.eq.0) then

            istgt = 0

          endif

        elseif((istgn.eq.2).or.(istgn.eq.3))then

          if(istgnold.eq.1) then

            do i = 1,3
              ch2(p1(26)+i-1) = ch1(p1(26)+i-1) + xnew(i) - xold(i)
            end do ! i

            istgt = 1

          elseif((istgnold.eq.2).or.(istgnold.eq.3)) then

            if((nr.eq.ch1(p1(1))).and.(nrneben.eq.ch1(p1(25)))) then

              ch2(p1(10)  ) = xi1 - ch1(p1(24))

              istgt = 1

            else

              do i = 1,3
                ch2(p1(26)+i-1) = ch1(p1(26)+i-1) + xnew(i) - xold(i)
              end do ! i

              istgt = 1

            endif

          elseif(istgnold.eq.4) then

            do i = 1,3
              ch2(p1(26)+i-1) = ch1(p1(26)+i-1) + xnew(i) - xold(i)
            end do ! i

            istgt = 1

          elseif(istgnold.eq.0) then

            istgt = 0

          endif

        elseif(istgn.eq.4) then

          if(istgnold.eq.1) then

            do i = 1,3
              ch2(p1(26)+i-1) = ch1(p1(26)+i-1) + xnew(i) - xold(i)
            end do ! i

            istgt = 1

          elseif((istgnold.eq.2).or.(istgnold.eq.3)) then

            do i = 1,3
              ch2(p1(26)+i-1) = ch1(p1(26)+i-1) + xnew(i) - xold(i)
            end do ! i

            istgt = 1

          elseif(istgnold.eq.4) then

            if (nr.eq.ch1(p1(1))) then

              istgt = 0

            else
              do i = 1,3
                ch2(p1(26)+i-1) = ch1(p1(26)+i-1) + xnew(i) - xold(i)
              end do ! i

              istgt = 1

            endif

          elseif(istgnold.eq.0) then

            istgt = 0

          endif

c       No normal contact => no tangential contact

        elseif(istgn.eq.0) then

          istgt = 0

        endif

        ch2(p1(3)) = istgt

      endif

3000  format(' GEOQTQ: ISTGN:',i2,' Slave node',i5,
     &       ' with master segments',2i5)

      end
