c$Id:$
      subroutine nodew(ndw,msum,ix,ixc,nnac,nnid,nen,nen1,ncen,ncen1,
     &                 numnp,numel,numcels,nstart)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute node and element weights for profile
c               minimizations

c      Inputs:
c         ix(nen1,*)     - Element nodal connection list
c         ixc(ncen1,*)   - Contact element nodal connection list
c         nnac(numel)    - Overlay markers
c         nnid(numnp)    - Number active dof at nodes
c         nen            - Maximum number nodes/element
c         nen1           - Dimension of ix  array
c         ncen           - Maximum number nodes/contact element
c         ncen1          - Dimension of ixc array
c         numnp          - Number of nodes in mesh
c         numel          - Number of elements in mesh
c         numcels        - Number of contact elements in mesh

c      Outputs:
c         ndw(*)         - Node weights
c         msum(*)        - Element weights
c         nstart
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   k, l, m,minw, n,nel,nen,nen1,ncen,ncen1
      integer   numnp,numel,numcels,nstart, nsum
      integer   ndw(numnp),msum(numel),ix(nen1,numel),ixc(ncen1,*)
      integer   nnac(numel),nnid(numnp)

      save

c     Evaluation of node weights

      do n = 1,numnp
        ndw(n) = 1
      end do ! n
      minw = 1

c     Loop five times to evaluate node weights

      do k = 1,5

c       Normalize node weights

        do n = 1,numnp
          ndw(n) = ndw(n) / minw
        end do ! n

c       Compute element weights

        do m = 1,numel
          if(msum(m).gt.0.and.ix(nen1-1,m).ge.0.and.nnac(m).eq.0) then
            nsum = 0
            nel  = 0
            do l = 1,nen
              n = abs(ix(l,m))
              if(n.gt.0) then
                if(nnid(n).gt.0) then
                  nel  = nel + 1
                  nsum = nsum + ndw(n)
                endif
              endif
            end do ! l
            nsum     = min(nsum,numnp)
            msum(m)  = nsum / max(1,nel)
          endif
        end do ! m
        do m = 1,numcels
          if(msum(m+numel).gt.0) then
            nsum = 0
            nel  = 0
            do l = 1,ncen
              n = abs(ixc(l,m))
              if(n.gt.0) then
                if(nnid(n).gt.0) then
                  nel  = nel + 1
                  nsum = nsum + ndw(n)
                endif
              endif
            end do ! l
            nsum          = min(nsum,numnp)
            msum(m+numel) = nsum / max(1,nel)
          endif
        end do ! m

c       Compute node weights

        minw = 0
        do n = 1,numnp
          ndw(n) = 0
        end do ! n
        do m = 1,numel
          if (msum(m).gt.0.and.ix(nen1-1,m).ge.0.and.nnac(m).eq.0) then
            nsum = 0
            do l = 1,nen
              n = abs(ix(l,m))
              if(n.gt.0) then
                if(nnid(n).gt.0) then
                  ndw(n) = ndw(n) + msum(m)
                  minw   = max(minw,ndw(n))
                endif
              endif
            end do ! l
          endif
        end do ! m
        do m = 1,numcels
          if (msum(m+numel).gt.0) then
            nsum = 0
            do l = 1,ncen
              n = abs(ixc(l,m))
              if(n.gt.0) then
                if(nnid(n).gt.0) then
                  ndw(n) = ndw(n) + msum(m+numel)
                  minw   = max(minw,ndw(n))
                endif
              endif
            end do ! l
          endif
        end do ! m

c       Find minimum value

        minw = 32000000
        do n = 1,numnp
          if (ndw(n).gt.0 .and. ndw(n).le.minw) then
            minw   = ndw(n)
            nstart = n
          endif
        end do ! n
      end do ! k

      end
