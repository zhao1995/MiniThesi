c$Id:$
      subroutine sxblk(nsblk,ndm, nns, ixs, xs, ts, nml, x, prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generates coordinates for surface block

c      Inputs:
c         nsblk     - Number of surface blocks
c         ndm       - Spatial dimension of mesh
c         nns(5,*)  - Block information for surface block
c         ixs(*)    - Master node numbers for blocks
c         xs(ndm,*) - Surface nodal coordinates for blocks
c         ts(*)     - Director thickness for block
c         nml(3,9,*)- Normal vector to master nodes
c         prt       - Flag, output results if true
c         prth      - Flag, output title/header if true

c      Outputs:
c         x(ndm,*)  - Nodal coordinates of blocks
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      logical   phd, prt, prth
      integer   nsblk,ndm,i,j,k,l,ii,jj,j1,j2,mi,mct,nn,nr,ns,nt,n
      integer   nns(5,nsblk), ixs(9,nsblk), ixl(27), ixlmap(9,2)
      real*8    dr,ds,dt, x(ndm,*), ss(3), xl(3,27)
      real*8    xs(ndm,9,nsblk), nml(3,9,nsblk), ts(9,nsblk)

      save

c     Map from surface to 27-node brick numbers

      data ixlmap/1,2,3,4,13,14,15,16,17,5,6,7,8,18,19,20,21,22/

      mct = 0
      do mi = 1,nsblk
        nn = nns(1,mi)
        nr = nns(2,mi)
        ns = nns(3,mi)
        nt = nns(4,mi)
        n  = nns(5,mi)

        do j = 1,27
          ixl(j) = 0
        end do ! j

c       Set generation increments of natural coordinates

        dr = 2.d0/nr
        ds = 2.d0/ns
        dt = 2.d0/nt

c       3-d generations

        ss(2) = -1.0d0
        do j = 1,ns+1
          ss(1) = -1.0d0
          do i = 1,nr+1
            ss(3) = -1.0d0
            do k = 1,nt+1

c             Compute coordinates of node

              do jj = 1,nn
                if(ixs(jj,mi).ne.0) then
                  j1 = ixlmap(jj,1)
                  j2 = ixlmap(jj,2)
                  ixl(j1) = j1
                  ixl(j2) = j2
                  do ii = 1,ndm
                    xl(ii,j1) = xs(ii,jj,mi)
                    xl(ii,j2) = xs(ii,jj,mi) + nml(ii,jj,mi)*ts(jj,mi)
                  end do ! ii
                end if
              end do ! jj
              call xcor3d(ss,xl,ixl, x(1,n))

c             Output point

              if(prt) then
                mct = mct + 1
                phd = mod(mct,50).eq.1
                if(phd) then
                  call prtitl(prth)
                  write(iow,2000) (l,'-coord',l=1,ndm)
                  write(iow,2001) n,(x(l,n),l=1,ndm)
                  if(ior.lt.0) then
                    write(*,2000) (l,'-coord',l=1,ndm)
                    write(*,2001) n,(x(l,n),l=1,ndm)
                  endif
                endif
              endif
              n = n + 1
              ss(3) = ss(3) + dt
            end do ! k
            ss(1) = ss(1) + dr
          end do ! i
          ss(2) = ss(2) + ds
        end do ! j
      end do ! mi

c     Formats

2000  format(/'  N o d a l   C o o r d i n a t e s'//6x,'Node',5(i7,a6))

2001  format(i10,5f13.4)

      end
