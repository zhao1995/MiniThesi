c$Id:$
      subroutine vblke(nr,ns,nt,x,ix,ni,ne,nf,ndm,nen1,mat,ntyp,
     &                 dlayer,ilr, prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set ninc = 2 for ntyp.le.15                      02/07/2007
c       2. Set ninc = 2 for ntyp.le.18                      04/09/2007
c       3. Replace itl by itq                               05/09/2007
c       4. Add 14/15 node generation                        06/09/2007
c       5. Correct allocation of active node for 14/15 node 20/10/2007
c       6. Add sets for 64-node brick type elements         06/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate a block of 3-d 8-node brick elements

c      Inputs:
c         nr        - Number elements in 1-local coordinate dir.
c         ns        - Number elements in 2-local coordinate dir.
c         nt        - Number elements in 3-local coordinate dir.
c         ni        - Initial node number for block
c         ne        - Initial element number for block
c         ndm       - Spatial dimension of coordinate array
c         nen1      - Dimension of ix array
c         mat       - Material set number for block
c         ntyp      - Element type for generations
c                     10:  8-node hexahedron  elements
c                     11:  4-node tetrahedron elements
c                     12: 27-node hexahedron  elements
c                     13: 10-node tetrahedron elements
c                     14: 20-node hexahedron  elements
c                     15: 11-node tetrahedron elements
c                     17: 14-node tetrahedron elements
c                     18: 15-node tetrahedron elements
c         prt       - Output generatied coordinates if true
c         prth      - Output title/header data if true

c      Outputs:
c         x(ndm,*)  - Nodal coordinates
c         ix(*)     - Element nodal connection list for block
c         nf        - Final   element number for block
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat2.h'
      include  'compac.h'
      include  'iofile.h'
      include  'pconstant.h'
      include  'trdata.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   prt,prth, phd
      integer   ni,nf,ne,nn,ndm,nen1,mat,ma,ntyp, mct
      integer   nr,ns,nt,nrs,i,j,k,l,m,n,dlayer, ninc, nodi
      integer   ix(nen1,*),iq(27),it(10),ilr(*),nd(27),itq(10,6)
      real*8    x(ndm,*)

      save

      data      itq / 1,2,4,5,  9, 21, 12, 17, 25, 23,
     &                2,3,4,8, 10, 11, 21, 27, 26, 20,
     &                2,4,5,8, 21, 23, 25, 27, 20, 16,
     &                2,6,3,8, 18, 24, 10, 27, 22, 26,
     &                3,6,7,8, 24, 14, 19, 26, 22, 15,
     &                5,6,2,8, 13, 18, 25, 16, 22, 27 /

c     Check generation order

      do i = 1,10
        it(i) = i
      end do ! i
      do i = 1,27
        iq(i) = i
      end do ! i
      if    (ntyp.eq.12) then
        nn = 27
      elseif(ntyp.eq.14) then
        nn = 20
      elseif(ntyp.eq.19) then
        nn = 64
      else
        nn = 8
      endif
      if(trdet.lt.0.0d0) then
        do i = 1,4
          iq(i+4) = i
          iq(i  ) = i+4
        end do ! i
        if(ntyp.eq.12 .or. ntyp.eq.14) then
          do i = 9,12
            iq(i+4) = i
            iq(i  ) = i+4
          end do ! i
          iq(21) = 22
          iq(22) = 21
        endif
        i     = it(2)
        it(2) = it(3)
        it(3) = i
        if(ntyp.eq.13 .or. ntyp.eq.15 .or.
     &     ntyp.eq.17 .or. ntyp.eq.18) then
          i      = it( 5)
          it( 5) = it( 7)
          it( 7) = i
          i      = it( 9)
          it( 9) = it(10)
          it(10) = i
        endif
      endif

      nrs = nr*ns
      if(ntyp.le.11) then
        ninc  =  1
        nd(1) = -1
        nd(2) =  0
        nd(3) =      nr
        nd(4) = -1 + nr
        nd(5) = -1      + nrs
        nd(6) =           nrs
        nd(7) =      nr + nrs
        nd(8) = -1 + nr + nrs
      elseif(ntyp.le.18) then
        ninc  =  2
        nd( 1) = -1
        nd( 2) =  1
        nd( 3) =  1 + nr*2
        nd( 4) = -1 + nr*2
        nd( 5) = -1        + nrs*2
        nd( 6) =  1        + nrs*2
        nd( 7) =  1 + nr*2 + nrs*2
        nd( 8) = -1 + nr*2 + nrs*2
        nd( 9) =  0
        nd(10) =  1 + nr
        nd(11) =      nr*2
        nd(12) = -1 + nr
        nd(13) =             nrs*2
        nd(14) =  1 + nr   + nrs*2
        nd(15) =      nr*2 + nrs*2
        nd(16) = -1 + nr   + nrs*2
        nd(17) = -1        + nrs
        nd(18) =  1        + nrs
        nd(19) =  1 + nr*2 + nrs
        nd(20) = -1 + nr*2 + nrs
        nd(21) =      nr
        nd(22) =      nr   + nrs*2
        nd(23) = -1 + nr   + nrs
        nd(24) =  1 + nr   + nrs
        nd(25) =             nrs
        nd(26) =      nr*2 + nrs
        nd(27) =      nr   + nrs
      endif

c     Compute element connections

      nodi = ni + nr*ns*nt - 1     ! Start number for internal nodes
      if(dlayer.ge.0) then
        ma = mat
      endif
      nf = ne - 1
      do k = 1,nt-1,ninc
        if(dlayer.eq.3) then
          ma = ilr(k)
        endif
        do j = 1,ns-1,ninc
          if(dlayer.eq.2) then
            ma = ilr(j)
          endif
          n = nr*(j-1 + ns*(k-1)) + ni
          do i = 1,nr-1,ninc
            if(dlayer.eq.1) then
              ma = ilr(i)
            endif
            n = n + 1

c           8-node hexahedral elements

            if(ntyp.eq.10) then
              nf = nf + 1
              ix(nen1,nf) = ma
              do m = 1,8
                ix(iq(m),nf) = n + nd(m)
              end do ! m

c           4-node tetrahedral elements

            elseif(ntyp.eq.11) then
              do l = 1,6
                nf = nf + 1
                ix(nen1,nf) = ma
                do m = 1,4
                  ix(it(m),nf) = n + nd(itq(m,l))
                end do ! m
              end do ! l

c           20 and 27-node hexahedral elements

            elseif(ntyp.eq.12 .or. ntyp.eq.14) then
              nf = nf + 1
              ix(nen1,nf) = ma
              do m = 1,nn
                ix(iq(m),nf) = n + nd(m)
              end do ! m

c           10 and 14-node tetrahedral elements

            elseif(ntyp.eq.13 .or. ntyp.eq.17) then
              do l = 1,6
                nf = nf + 1
                ix(nen1,nf) = ma
                do m = 1,10
                  ix(it(m),nf) = n + nd(itq(m,l))
                end do ! m
              end do ! l

c           11 and 15-node tetrahedral elements

            elseif(ntyp.eq.15 .or. ntyp.eq.18) then
              do l = 1,6
                nf = nf + 1
                ix(nen1,nf) = ma
                do m = 1,10
                  ix(it(m),nf) = n + nd(itq(m,l))
                end do ! m

c               Internal point

                mr(np(190)+nodi) = 0  ! Activate node
                nodi      = nodi + 1
                if(ntyp.eq.15) then
                  ix(11,nf) = nodi
                else
                  ix(15,nf) = nodi
                endif

                do m = 1,ndm
                  x(m,nodi) = 0.250d0*(x(m,ix(it( 5),nf))
     &                               + x(m,ix(it( 6),nf))
     &                               + x(m,ix(it( 7),nf))
     &                               + x(m,ix(it( 8),nf))
     &                               + x(m,ix(it( 9),nf))
     &                               + x(m,ix(it(10),nf)))
     &                      - 0.125d0*(x(m,ix(it( 1),nf))
     &                               + x(m,ix(it( 2),nf))
     &                               + x(m,ix(it( 3),nf))
     &                               + x(m,ix(it( 4),nf)))
                end do ! m

              end do ! l

            endif

            n = n + ninc - 1

          end do ! i
        end do ! j
      end do ! k

c     64-node hexahedron

      if(ntyp.eq.19) then

        nf = ne
        call pqrblk(3,3,3, nr,ns,nt, nf,ni, ma, ix,nen1)

c     14- and 15-node terahedral element face node generation

      elseif(ntyp.eq.17 .or. ntyp.eq.18) then

c       Do numbers in 'nr' direction

        nn = ne -1
        do k = 1,nt-1,ninc
          do j = 1,ns-1,ninc
            do i = 1,nr-1,ninc
              ix(13,nn+1) = nodi + 1
              ix(12,nn+3) = nodi + 2

              ix(12,nn+1) = nodi + 3
              ix(14,nn+3) = nodi + 3

              ix(11,nn+2) = nodi + 4
              ix(13,nn+4) = nodi + 4

              ix(13,nn+2) = nodi + 5
              ix(11,nn+3) = nodi + 5

              ix(13,nn+3) = nodi + 6
              ix(13,nn+6) = nodi + 6

              ix(11,nn+4) = nodi + 7
              ix(12,nn+6) = nodi + 7

              ix(12,nn+4) = nodi + 8
              ix(11,nn+5) = nodi + 8

              ix(14,nn+4) = nodi + 9
              ix(14,nn+5) = nodi + 10

c             Coordinates for new points

              do m = 0,9
                mr(np(190)+nodi+m) = 0  ! Activate node
              end do ! m
              do m = 1,ndm
                x(m,nodi+ 1) = one3*(x(m,ix(3,nn+1))
     &                             + x(m,ix(1,nn+1))
     &                             + x(m,ix(4,nn+1)))
                x(m,nodi+ 2) = one3*(x(m,ix(2,nn+3))
     &                             + x(m,ix(3,nn+3))
     &                             + x(m,ix(4,nn+3)))
                x(m,nodi+ 3) = one3*(x(m,ix(2,nn+1))
     &                             + x(m,ix(3,nn+1))
     &                             + x(m,ix(4,nn+1)))
                x(m,nodi+ 4) = one3*(x(m,ix(1,nn+2))
     &                             + x(m,ix(2,nn+2))
     &                             + x(m,ix(4,nn+2)))
                x(m,nodi+ 5) = one3*(x(m,ix(3,nn+2))
     &                             + x(m,ix(1,nn+2))
     &                             + x(m,ix(4,nn+2)))
                x(m,nodi+ 6) = one3*(x(m,ix(3,nn+3))
     &                             + x(m,ix(1,nn+3))
     &                             + x(m,ix(4,nn+3)))
                x(m,nodi+ 7) = one3*(x(m,ix(1,nn+4))
     &                             + x(m,ix(2,nn+4))
     &                             + x(m,ix(4,nn+4)))
                x(m,nodi+ 8) = one3*(x(m,ix(2,nn+4))
     &                             + x(m,ix(3,nn+4))
     &                             + x(m,ix(4,nn+4)))
                x(m,nodi+ 9) = one3*(x(m,ix(1,nn+4))
     &                             + x(m,ix(2,nn+4))
     &                             + x(m,ix(3,nn+4)))
                x(m,nodi+10) = one3*(x(m,ix(1,nn+5))
     &                             + x(m,ix(2,nn+5))
     &                             + x(m,ix(3,nn+5)))
              end do ! m

              nn          = nn   + 6
              nodi        = nodi + 8
            end do ! i
            nodi = nodi + 2
          end do ! j
        end do ! k

c       Do nodes in ns direction

        do k = 1,nt-1,ninc
          nn = ne + 3*(nr-1)*(ns-1)*(k-1)/4 - 1
          do i = 1,nr-1,ninc
            do j = 1,ns-1,ninc
              ix(11,nn+1) = nodi + 1
              ix(14,nn+6) = nodi + 2

              ix(12,nn+2) = nodi + 3
              ix(13,nn+5) = nodi + 4

c             Coordinates for new points

              do m = 0,3
                mr(np(190)+nodi+m) = 0  ! Activate node
              end do ! m
              do m = 1,ndm
                x(m,nodi+ 1) = one3*(x(m,ix(1,nn+1))
     &                             + x(m,ix(2,nn+1))
     &                             + x(m,ix(4,nn+1)))
                x(m,nodi+ 2) = one3*(x(m,ix(1,nn+6))
     &                             + x(m,ix(2,nn+6))
     &                             + x(m,ix(3,nn+6)))
                x(m,nodi+ 3) = one3*(x(m,ix(2,nn+2))
     &                             + x(m,ix(3,nn+2))
     &                             + x(m,ix(4,nn+2)))
                x(m,nodi+ 4) = one3*(x(m,ix(3,nn+5))
     &                             + x(m,ix(1,nn+5))
     &                             + x(m,ix(4,nn+5)))
              end do ! m

              nn          = nn + 3*nr - 3
              nodi        = nodi + 2
            end do ! j
            nn   = nn - 3*(nr-1)*(ns-1)/2 + 6
            nodi = nodi + 2
          end do ! i
        end do ! k

c       Do nodes in nt direction

        nn = ne - 1
        do j = 1,ns-1,ninc
          nn = ne + 3*(nr-1)*(j-1)/2 - 1
          do i = 1,nr-1,ninc
            do k = 1,nt-1,ninc
              ix(14,nn+1) = nodi + 1
              ix(14,nn+2) = nodi + 2

              ix(11,nn+6) = nodi + 3
              ix(12,nn+5) = nodi + 4

c             Coordinates for new points

              do m = 0,3
                mr(np(190)+nodi+m) = 0  ! Activate node
              end do ! m
              do m = 1,ndm
                x(m,nodi+ 1) = one3*(x(m,ix(1,nn+1))
     &                             + x(m,ix(2,nn+1))
     &                             + x(m,ix(3,nn+1)))
                x(m,nodi+ 2) = one3*(x(m,ix(1,nn+2))
     &                             + x(m,ix(2,nn+2))
     &                             + x(m,ix(3,nn+2)))
                x(m,nodi+ 3) = one3*(x(m,ix(1,nn+6))
     &                             + x(m,ix(2,nn+6))
     &                             + x(m,ix(4,nn+6)))
                x(m,nodi+ 4) = one3*(x(m,ix(2,nn+5))
     &                             + x(m,ix(3,nn+5))
     &                             + x(m,ix(4,nn+5)))
              end do ! m

              nn          = nn + 3*(nr-1)*(ns-1)/2
              nodi        = nodi + 2
            end do ! k
            nn   = nn - 3*(nr-1)*(ns-1)*(nt-1)/4 + 6
            nodi = nodi + 2
          end do ! i
        end do ! j

      endif ! 14- and 15-node tets

c     Output point

      nn = ni + nr*ns*nt
      if(prt .and. nodi.gt.nn) then
        mct = 0
        do n = nn,nodi
          mct = mct + 1
          phd = mod(mct,50).eq.1
          call prtitl(prth.and.phd)
          if(phd) write(iow,2000) (l,'-Coordinate',l=1,ndm)
          write(iow,2001) n,(x(l,n),l=1,ndm)
          if(ior.lt.0) then
            if(phd) write(*,2000) (l,'-Coordinate',l=1,ndm)
            write(*,2001) n,(x(l,n),l=1,ndm)
          endif
        end do ! n
      endif

c     Formats

2000  format(/'  N o d a l   C o o r d i n a t e s'//
     &        '     Node',5(i4,a11))

2001  format(i8,1x,1p,5e15.4)

      end
