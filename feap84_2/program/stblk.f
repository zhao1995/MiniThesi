c$Id:$
      subroutine stblk(nsblk,ndm, nns,ixn,ixs,xs, ns, xmx)

c      * * F E A P * * A Finite Elemen4 Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate unique normal vectors for master block

c      Inputs:
c         nsblk     - Number of blocks
c         ndm       - Spatial dimension of block
c         nns(5,*)  - ??
c         ixs(9,*)  - List of surface nodes for blocks
c         xs(ndm,*) - Nodal coordinates of block
c         xmx       - Maximum length of ??

c      Outputs:
c         ixn(9,*)  - List of
c         ns(3,*,*) - Normals for block
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   nsblk,ndm, i,j,mi,mj,ni,nj, nd,node
      integer   nns(5,nsblk),ixn(9,nsblk),ixs(9,nsblk)
      real*8    xj1,xj2,xj3, d, xmx, tol
      real*8    xs(ndm,9,nsblk), ns(3,9,nsblk)

      save

      data      tol/1.d-06/

c     Compute normals for each bottom surface master node

      call pzeroi(ixn,9*nsblk)
      node = 0
      ni   = nns(1,1)
      do i = 1,ni
        if(ixs(i,1).gt.0) then
          ixn(i,1) = i
          node = node + 1
        end if
      end do ! i
      node = ni

      do mi = 1,nsblk
        ni = nns(1,mi)
        do mj = mi+1,nsblk
          nj = nns(1,mj)

          do j = 1,nj
            if(ixs(j,mj).gt.0) then
              xj1 = xs(1,j,mj)
              xj2 = xs(2,j,mj)
              xj3 = xs(3,j,mj)

              do i = 1,ni
                if(ixs(i,mi).gt.0) then
                  d = (xs(1,i,mi)-xj1)**2 + (xs(2,i,mi)-xj2)**2
     &              + (xs(3,i,mi)-xj3)**2
                  if(d.lt.tol*xmx*xmx) then
                    ixn(j,mj) = ixn(i,mi)
                  end if
                end if
              end do ! i
            end if
          end do ! j
          do j = 1,nj
            if(ixn(j,mj).eq.0 .and. ixs(j,mj).gt.0) then
              node      = node + 1
              ixn(j,mj) = node
            end if
          end do ! j
        end do ! mj
      end do ! mi

c     Average normals for all blocks to form continuous mesh

      do nd = 1,node
        do mi = 1,nsblk
          ni = nns(1,mi)
          do i = 1,ni
            if( ixn(i,mi).eq.nd ) then

              do mj = mi+1,nsblk
                nj = nns(1,mj)
                do j = 1,nj

                  if( ixn(j,mj).eq.nd )  then

c                   Accumulate normal for node

                    ns(1,i,mi) = ns(1,i,mi) + ns(1,j,mj)
                    ns(2,i,mi) = ns(2,i,mi) + ns(2,j,mj)
                    ns(3,i,mi) = ns(3,i,mi) + ns(3,j,mj)

                    d = 1.d0/sqrt(ns(1,i,mi)**2 + ns(2,i,mi)**2
     &                          + ns(3,i,mi)**2)

                    ns(1,i,mi) = ns(1,i,mi)*d
                    ns(2,i,mi) = ns(2,i,mi)*d
                    ns(3,i,mi) = ns(3,i,mi)*d

                  end if
                end do ! j
              end do ! mj
            end if
          end do ! i
        end do ! mi
      end do ! nd

c     Transfer final averaged normals to all other nodes

      do nd = 1,node
        do mi = 1,nsblk
          ni = nns(1,mi)
          do i = 1,ni
            if( ixn(i,mi).eq.nd ) then

              do mj = mi+1,nsblk
                nj = nns(1,mj)
                do j = 1,nj

                  if( ixn(j,mj).eq.nd )  then

                    ns(1,j,mj) = ns(1,i,mi)
                    ns(2,j,mj) = ns(2,i,mi)
                    ns(3,j,mj) = ns(3,i,mi)

                  end if
                end do ! j
              end do ! mj
            end if
          end do ! i
        end do ! mi
      end do ! nd

      end
