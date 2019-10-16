c$Id:$
      subroutine sh3finte (xl,ul,xln,ndof,ndm,ndf,nel, dir, xi, xj)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c      1. Add increment of director in dir(:,:,5:7)         04/06/2013
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:   SH3FINTE is subroutine which initializes
c                       appropriately local element director field
c                       according on number of degrees of freedom
c                       of nodes in element; as defined by the
c                       integer array ndof(nel). If a node posseses 6
c                       dof's it indicates presence of a shell
c                       intersection; otherwise, node is associated
c                       with a smooth shell mid-surface.

c        Author:        M.S. Rifai J.C. Simo D.D. Fox

c        Date:          January 1992.
c        Revised:       February 1997: reduce to 1-call
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        xl(ndm,*) ..... Local nodal coordinates of element
c        xln(3,3,9,4) .. local nodal rotation matrices
c                        xln(3,3,9,1): time t_n
c                        xln(3,3,9,2): time t_n+a
c                        xln(3,3,9,3): time t_n+1
c                        xln(3,3,9,4): time t_0
c        ndof(nel) ..... Array containing number of DOF/node
c        ndm ........... Spacial dimension (ndm=3)
c        nel ........... Number of nodes in element

c        Routine Output:
c        ---------------
c        dir(3,9,4) .... Local element director field
c                        dir(3,9,1): time t_n
c                        dir(3,9,2): time t_n+a
c                        dir(3,9,3): time t_n+1
c                        dir(3,9,4): time t_0
c                        dir(3,9,5): time t_n - t_0
c                        dir(3,9,6): time t_n+a - t_0
c                        dir(3,9,7): time t_n+1 - t_0
c        xi(3,3,9) ..... Transformation matrix for nodal DOF at t_n+a
c        xj(3,3,9) ..... Transformation matrix for nodal DOF at t_n+1
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'fdata.h'

      integer   i, j, k, l, m, node, ndm, ndf, nel, ntime
      integer   ii(4), jj(4), kk(4), ndof(nel)
      real*8    xnorm, facn, fac1,tt1(3), tt2(3), tt3(3), dir(3,9,7)
      real*8    xl(ndm,*), xln(3,3,9,4), xi(3,3,9), xj(3,3,9)
      real*8    ul(ndf,*), v3(3), len

      save

c     Scatter Data (nel = 4)

      data      ii / 1 , 2 , 3 , 4 /
      data      jj / 2 , 3 , 4 , 1 /
      data      kk / 3 , 4 , 1 , 2 /

c     Define local element directors

      fac1 = theta(3)
      facn = 1.d0 - fac1

      do j = 1,4
        node = jj(j)

        if    (ndof(node) .eq. 5) then
          do ntime = 1,4
            do l = 1,3
              dir(l,node,ntime) = xln(l,3,node,ntime)
            end do ! l
          end do ! ntime

c         Interpolate director to t_n+alpha
          if(fl(9))then
            do l=1,3
              dir(l,node,2) = facn*dir(l,node,1) + fac1*dir(l,node,3)
            end do ! l
          endif

        elseif(ndof(node) .eq. 6) then

          i = ii(j)
          k = kk(j)

c         Compute reference directors

          do l = 1,3
            tt1(l) = xl(l,k) - xl(l,node)
            tt2(l) = xl(l,i) - xl(l,node)
          end do ! l
          call vecp( tt1, tt2, tt3 )
          xnorm = 1.d0/sqrt( tt3(1)**2 + tt3(2)**2 + tt3(3)**2 )
          do l = 1,3
            dir(l,node,4) = tt3(l) * xnorm
          end do ! l

c         Compute current directors

          do l = 1,3

            tt3(l) = xln(1,l,node,4)*dir(1,node,4)
     &             + xln(2,l,node,4)*dir(2,node,4)
     &             + xln(3,l,node,4)*dir(3,node,4)

            dir(l,node,1) = xln(l,1,node,1)*tt3(1)
     &                    + xln(l,2,node,1)*tt3(2)
     &                    + xln(l,3,node,1)*tt3(3)

            dir(l,node,3) = xln(l,1,node,3)*tt3(1)
     &                    + xln(l,2,node,3)*tt3(2)
     &                    + xln(l,3,node,3)*tt3(3)

          end do ! l
          if(fl(9))then
            do l = 1,3
              dir(l,node,2) = facn*dir(l,node,1) + fac1*dir(l,node,3)
            end do ! l
          else
            do l = 1,3
              dir(l,node,2) = xln(l,1,node,2)*tt3(1)
     &                      + xln(l,2,node,2)*tt3(2)
     &                      + xln(l,3,node,2)*tt3(3)
            end do ! l
          endif
        endif

c       Increment of director

        do ntime = 1,3
          dir(:,node,ntime+4) = dir(:,node,ntime) - dir(:,node,4)
        end do ! ntime

c       Linearized rotation for small displacements

        v3(:) = xln(:,1,node,4)*ul(4,node) + xln(:,2,node,4)*ul(5,node)
        len   = v3(1)**2 +v3(2)**2 +v3(3)**2
        if(len.lt.1.d-04) then
          dir(:,node,7) = v3(:)
        endif

      end do ! j

c     Compute transformation matrix

      do node = 1,4

        if (ndof(node) .eq. 5) then
          do m = 1,3
            do l = 1,3
              xi(l,m,node) = xln(l,m,node,2)
              xj(l,m,node) = xln(l,m,node,3)
            end do ! l
          end do ! m

        elseif (ndof(node) .eq. 6) then
          xi(1,1,node) =  0.d0
          xi(2,2,node) =  0.d0
          xi(3,3,node) =  0.d0
          xi(1,2,node) =  dir(3,node,2)
          xi(2,3,node) =  dir(1,node,2)
          xi(3,1,node) =  dir(2,node,2)
          xi(2,1,node) = -dir(3,node,2)
          xi(3,2,node) = -dir(1,node,2)
          xi(1,3,node) = -dir(2,node,2)

          xj(1,1,node) =  0.d0
          xj(2,2,node) =  0.d0
          xj(3,3,node) =  0.d0
          xj(1,2,node) =  dir(3,node,3)
          xj(2,3,node) =  dir(1,node,3)
          xj(3,1,node) =  dir(2,node,3)
          xj(2,1,node) = -dir(3,node,3)
          xj(3,2,node) = -dir(1,node,3)
          xj(1,3,node) = -dir(2,node,3)
        endif

      end do ! node

      end
