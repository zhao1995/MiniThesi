c$Id:$
      subroutine rnodld(ra,r,x,rcg,ixt,nrk,
     &                  ndm,ndf,rlam,umode,n,imf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Computation of residual for rigid body loads

c      Inputs:
c         x(3)        - Point of nodal load in reference state
c         rcg(3,11,*) - Center of mass locations
c         ixt         - Indicator for rigid node
c         nrk         - Pointer to type of rate terms to use
c         ndm         - Spatial dimension of mesh
c         ndf         - Number dof/node

c      Outputs:
c         ra(3)       - Load at t_n+a
c         r(3)        - Load at t_n+1
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'evdata.h'
      include  'modreg.h'

      integer   i, j, n, ixt, nrk, ndm, ndf
      integer   imf(ndf,*)
      real*8    ra(3),r(3),x(3), rcg(3,11,*), rlam(9,6,*), ybig(3)
      real*8    lamu(3), lamuh(3), lam1(3,3), lamn(3,3), lamh(3,3)
      real*8    umode(neqmf,*),umodeh(3)

      save

c     Compute position of center of mass at t_n+1

      call pzero(r,3)
      do i = 1,ndm
        ybig(i) = x(i) - rcg(i,1,ixt)
      end do ! i
      call quavec(rlam(1,3,ixt),ybig,r)

      if(nmbody.gt.0) then
        do i = 1,nmbody
          if(modbod(i).eq.ixt) then

            do j = 1,ndm
              if(imf(j,n).gt.0) then
                umodeh(j) = umode(imf(j,n),6)
              else
                umodeh(j) = 0.0d0
              endif
            end do ! j

            do j = 1,3
              lamu(j) = 0.0d0
            end do ! j
            call quavec(rlam(1,3,ixt),umodeh,lamu)
            do j = 1,ndm
              r(j) = r(j) + lamu(j)
            end do ! j

          end if
        end do ! i
      end if

c     Compute position of center of mass at t_n+a

      if(nrk.le.0) then

        do i = 1,ndm
          ra(i) = r(i)
        end do ! i

      else

        call quamat(rlam(1,3,ixt),lam1)
        call quamat(rlam(1,1,ixt),lamn)

        do i = 1,ndm
          do j = 1,ndm
            lamh(i,j) = 0.5d0*(lam1(i,j) + lamn(i,j))
          end do ! j
        end do ! i

        call pzero(ra,3)
        call pzero(lamuh,3)
        do i = 1,ndm
          do j = 1,ndm
            ra(i) = ra(i) + lamh(i,j)*ybig(j)
            lamuh(i) = lamuh(i) + lamh(i,j)*umodeh(j)
          end do ! j
        end do ! i

        if(nmbody.gt.0) then
          do i = 1,ndm
            ra(i) = ra(i) + lamuh(i)
          end do ! i
        end if

      endif

      end
