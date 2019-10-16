c$Id:$
      subroutine zzpro2(ix,ib,ip,iq, x, st, sp, sh,
     &                  ndm,nen,nen1,numnp,numel,numst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c     2. Increase dimension of 'ip' to numnp+1 to simplify  05/07/2011
c        coding of loops.
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Zienkiewicz-Zhu Projections

c      Inputs:
c        ix(nen1,*)      - Nodal connection array
c        ip(*)           - Pointer array for element patches
c        x(ndm,*)        - Nodal coordinates
c        ndm             - Mesh spatial dimension
c        nen             - Maximum number of nodes on an element
c        nen1            - Dimension of 'ix' array
c        numnp           - Number of total nodal points
c        numel           - Number of total elements
c        numst           - Dimension of 'st' array

c      Outputs:
c        ib(*)           - Used to average projected node values
c        iq(*)           - List of elements connected to each node
c        st(numnp,numst) - Stress prjection array
c        sp(numnp,3)     - Principal stress projection array
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'eldata.h'
      include   'eldatp.h'
      include   'strnum.h'
      include   'pconstant.h'
      include   'pointer.h'
      include   'comblk.h'
      include   'zzcom1.h'

      logical    fa, patchfl
      integer    ndm, nen,nen1,numnp,numel,numst
      integer    i,ii, j, k,kk, l, nn, nps

      real*8     cd,press

      integer    ix(nen1,numel), ib(numnp), ip(numnp),iq(*)
      real*8     x(ndm,*),st(numnp,numst), sp(numnp,*), sh(numnp,*)
      real*8     hh(10,10),aa(30,10),bb(30,10),poly(10),t(8), sig(9)
      real*8     tmp(10)

      real*8     cc(30,10),gg(30,10)

      integer    ixp(343), npn, ips(1000)

      save

      data       fa /.false./

c     4.) Set patch size parameters

      do i = 1,numst
        st(1:numnp,i) = 0.0d0
      end do ! i

      if(hpltfl) then
        do i = 1,plhmax
          sh(1:numnp,i) = 0.0d0
        end do ! i
      endif

c     5.) Build patches

      do nn = 1,ip(numnp+1)
        iq(nn) = 0
      end do ! nn

      do nn = 1,numel
        do i = 1,nen
          ii = ix(i,nn)
          if(ii.gt.0) then
            do j = ip(ii),ip(ii+1)-1
              if(iq(j).eq.0 .or. iq(j).eq.nn) then
                iq(j) = nn
                go to 300
              end if
            end do ! j
          end if ! ii > 0
 300      continue
        end do ! i
      end do ! nn

c     6.) Loop over patches to compute projected quantities

      ips(1:numnp) = 0

      do nn = 1,numnp
        if(ib(nn).gt.0) then

          call pzero(hh,100)
          call pzero(aa,300)
          call pzero(bb,300)
          if(hpltfl) then
            call pzero(cc,300)
            call pzero(gg,300)
          endif
          call pzeroi(ixp,343)
          do i = 1,ndm
            xnodz(i) = x(i,nn)
          end do ! i

c         Form stresses for patch

          npn = 0
          do j = ip(nn),ip(nn+1)-1

c           Get contribution from each element

            if(iq(j).gt.0) then

c             Set default projection quantities and initialize

              nums = numst
              do i = 1,10
                ek( 1:10  ,i) = 0.0d0
                est(1:nums,i) = 0.0d0
                if(hpltfl) then
                  ehp(1:plhmax,i) = 0.0d0
                endif
              end do ! i

c             Set patch

              do k = 1,nen
                if(ix(k,iq(j)).ge.0) then
                  patchfl = .true.
                  do i = 1,npn
                    if(ixp(i).eq.ix(k,iq(j))) then
                      patch fl = .false.
                      exit
                    endif
                  end do ! i
                  if(patchfl) then
                    npn      = npn + 1
                    ixp(npn) = ix(k,iq(j))
                  endif
                endif
              end do ! k

c             Compute projections

              call formfe(np(40),np(26),np(26),np(26),fa,fa,fa,fa,25,
     &                    iq(j),iq(j),1)

c             Default is linear variation

              nps = ndm + 1

c             Look at element type

              if(ndm.eq.2) then
                if(nel.eq.6 .or. nel.eq.8 .or. nel.eq.9) then
                  nps = 6
                elseif(nel.eq.16) then
                  nps = 10
                else
                  nps = 10
                endif
              else
                nps = 10
              endif

c             Assemble least square patch arrays

              do i = 1,nps
                hh(1:nps ,i) = hh(1:nps ,i) + ek(1:nps,i)
                aa(1:nums,i) = aa(1:nums,i) + est(1:nums,i)
                if(hpltfl) then
                  cc(1:plhmax,i) = cc(1:plhmax,i) + ehp(1:plhmax,i)
                endif
              end do ! i

            endif

          end do ! j

c         Solve

          call zinvert(hh,nps,10,tmp)

c         Stress/strain variables

          do i = 1,nums
            t(1:nps) = 0.0d0
            do k = 1,nps
              t(1:nps) = t(1:nps) + hh(1:nps,k)*aa(i,k)
            end do ! k
            bb(i,1:nps) = bb(i,1:nps) + t(1:nps)
          end do ! i

c         History variables

          if(hpltfl) then
            do i = 1,plhmax
              t(1:nps) = 0.0d0
              do k = 1,nps
                t(1:nps) = t(1:nps) + hh(1:nps,k)*cc(i,k)
              end do ! k
              gg(i,1:nps) = gg(i,1:nps) + t(1:nps)
            end do ! i
          endif

c         Assemble nodal stresses

          do i = 1,npn
            kk = ixp(i)
            if(kk.gt.0) then
c             if(ips(kk).le.0) then
              if(ib(kk).le.0) then

c               Increment ib(kk) to count number of nodal projections

c               ips(kk) = ips(kk) - 1
                ib(kk) = ib(kk) - 1

c               Compute polynomial

                poly(1) = 1.d0
                do k = 1,ndm
                  poly(k+1) = x(k,kk) - x(k,nn)
                end do ! k

c               Compute quadratic and cubic polynomials

                poly( 4) = poly(2)*poly(2)
                poly( 5) = poly(2)*poly(3)
                poly( 6) = poly(3)*poly(3)
                poly( 7) = poly(4)*poly(2)
                poly( 8) = poly(5)*poly(2)
                poly( 9) = poly(5)*poly(3)
                poly(10) = poly(6)*poly(3)

c               Project stress/strain to nodes

                do k = 1,nums
                  do l = 1,nps
                    st(kk,k) = st(kk,k) + bb(k,l)*poly(l)
                  end do ! l
                end do ! k

c               Project history to nodes

                if(hpltfl) then
                  do k = 1,plhmax
                    do l = 1,nps
                      sh(kk,k) = sh(kk,k) + gg(k,l)*poly(l)
                    end do ! l
                  end do ! k
                endif
              end if
            end if
          end do !

c        Assemble interior node patch value

          st(nn,1:nums) = bb(1:nums,1)
          if(hpltfl) then
            sh(nn,1:plhmax) = sh(nn,1:plhmax) + gg(1:plhmax,1)
          endif

        end if

      end do ! nn

c     7.) Boundary averages

      do nn = 1,numnp
c       if(ips(nn).lt.0) then

c         cd = 1.d0/dble(-ips(nn))
        if(ib(nn).lt.0) then

          cd = 1.d0/dble(-ib(nn))
          st(nn,1:nums) = st(nn,1:nums)*cd
          if(hpltfl) then
            sh(nn,1:plhmax) = sh(nn,1:plhmax)*cd
          endif
        end if

c       Three-dimensional Principal Stresses

        if(ndm.eq.3 .or. istp.eq.8) then
          sig(1)   = st(nn,1)
          sig(2)   = st(nn,2)
          sig(3)   = st(nn,3)
          sig(4)   = st(nn,4)
          sig(5)   = st(nn,5)
          sig(6)   = st(nn,6)
          call pstr3d(sig,sig(7))
          sp(nn,1) = sig(7)
          sp(nn,2) = sig(8)
          sp(nn,3) = sig(9)

c       Two-dimensional Principal Stresses

        elseif(ndm.eq.2) then
          sig(1)   = st(nn,1)
          sig(2)   = st(nn,2)
          sig(3)   = st(nn,3)
          sig(4)   = st(nn,4)
          call pstr2d(sig,sig(7))
          sp(nn,1) = sig(7)
          sp(nn,2) = sig(8)
          sp(nn,3) = (sig(7)-sig(8))*0.5d0
          sp(nn,4) = sig(9)
        endif

c       Compute mean stress and mises stress

        press    = (sp(nn,1) + sp(nn,2) + sp(nn,3))*one3
        sp(nn,5) = press
        sp(nn,6) = sqrt(1.5d0*((sp(nn,1) - press)**2
     &                       + (sp(nn,2) - press)**2
     &                       + (sp(nn,3) - press)**2))
        sp(nn,7) =        one3*(sp(nn,1) - press)**3
     &                        *(sp(nn,2) - press)**3
     &                        *(sp(nn,3) - press)**3
        if(sp(nn,7).lt.0.0d0) then
          sp(nn,7) = -(abs(sp(nn,7))**one3)
        else
          sp(nn,7) =  (abs(sp(nn,7))**one3)
        endif

      end do ! nn

      end
