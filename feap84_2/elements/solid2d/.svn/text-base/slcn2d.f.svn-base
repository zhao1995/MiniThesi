c$Id:$
      subroutine slcn2d(ix,sig,eps,p,s,se,nel,nes)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise plot for 6-node triangles                 08/03/2007
c       2. Compute projections for 8-node element           24/06/2009
c       3. Add projection of strains                        01/01/2013
c       4. Correct entries in triangle strain plots         05/02/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Project element variables to nodes

c      Inputs:
c        ix(*)    - Nodal connections for element
c        sig(nes,*) - Stresses at quadrature points
c        eps(  6,*) - Strains  at quadrature points
c        nel          - Number nodes on element
c        nes          - Dimension of stress array

c      Outputs:
c        p(nen)   - Weights for 'lumped' projection
c        s(nen,*) - Integral of variables
c        se(nen)  - Error projectors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'eldatp.h'
      include  'sdata.h'
      include  'prstrs.h'
      include  'strnum.h'
      include  'pconstant.h'
      include  'pointer.h'
      include  'qudshp.h'
      include  'comblk.h'

      integer   nel,nes, i,l, ix(*),ixl(9)
      real*8    p(*),s(nen,*),se(*),sig(nes,*),eps(6,*), shpl(3,9), xg

      integer   j
      real*8    mat(3,3),sigg(3,6),gsig(3,6),epsg(3,6),geps(3,6)

      save

      data      ixl /1,2,3,4,5,6,7,8,9/

      if(nel.eq.8) then
        do l = 1,lint
          call meanx(ix,hr(np(44)),ndm)
          call shp2d(sg2(1,l),hr(np(44)),shpl,jac(l),ndm,9,ixl,.false.)
          do i = 1,nel
            xg   = shpl(3,i)*jac(l)
            p(i) = p(i) + xg

c           Stress projections

            s(i,1) = s(i,1) + sig(1,l)*xg
            s(i,2) = s(i,2) + sig(2,l)*xg
            s(i,3) = s(i,3) + sig(3,l)*xg
            s(i,4) = s(i,4) + sig(4,l)*xg
            s(i,5) = s(i,5) + sig(5,l)*xg
            s(i,6) = s(i,6) + sig(6,l)*xg

c           Strain projections

            s(i, 7) = s(i, 7) + eps(1,l)*xg
            s(i, 8) = s(i, 8) + eps(2,l)*xg
            s(i, 9) = s(i, 9) + eps(3,l)*xg
            s(i,10) = s(i,10) + eps(4,l)*xg
            s(i,11) = s(i,11) + eps(5,l)*xg
            s(i,12) = s(i,12) + eps(6,l)*xg

c           Error estimation projection

            se(i)  = se(i)  + erav*xg
          end do ! i

        end do ! l

c     Linear projection on element

      elseif(nel.eq.6) then

        do j = 1,3
          do i = 1,3
            mat(i,j)    = 0.0d0
          end do ! i
          do i = 1,6
            sigg(j,i) = 0.0d0
            epsg(j,i) = 0.0d0
          end do ! i
        end do ! j
        do l = 1,lint
          do j = 1,3
            xg = el2(j,l)*el2(4,l)
            do i = 1,3
              mat(i,j) = mat(i,j) + el2(i,l)*xg
            end do ! i
            do i = 1,6
              sigg(j,i) = sigg(j,i) + sig(i,l)*xg
              epsg(j,i) = epsg(j,i) + eps(i,l)*xg
            end do ! i
          end do ! j
        end do ! l
        call invert(mat,3,3)
        do i = 1,6
          do j = 1,3
            gsig(j,i) = mat(j,1)*sigg(1,i)
     &                + mat(j,2)*sigg(2,i)
     &                + mat(j,3)*sigg(3,i)
            geps(j,i) = mat(j,1)*epsg(1,i)
     &                + mat(j,2)*epsg(2,i)
     &                + mat(j,3)*epsg(3,i)
          end do ! j
          p( i)   = 1.0d0

c         Project stresses

          s(1,i) =  gsig(1,i)
          s(2,i) =  gsig(2,i)
          s(3,i) =  gsig(3,i)
          s(4,i) = (gsig(1,i) + gsig(2,i))*0.5d0
          s(5,i) = (gsig(2,i) + gsig(3,i))*0.5d0
          s(6,i) = (gsig(3,i) + gsig(1,i))*0.5d0

c         Project strains

          s(1,i+6) =  geps(1,i)
          s(2,i+6) =  geps(2,i)
          s(3,i+6) =  geps(3,i)
          s(4,i+6) = (geps(1,i) + geps(2,i))*0.5d0
          s(5,i+6) = (geps(2,i) + geps(3,i))*0.5d0
          s(6,i+6) = (geps(3,i) + geps(1,i))*0.5d0
        end do ! i

c     Row sum lumped projection routine

      else

        do l = 1,lint

          do i = 1,nel

            xg   = shp2(3,i,l)*jac(l)
            p(i) = p(i) + xg

c           Stress projections

            s(i,1) = s(i,1) + sig(1,l)*xg
            s(i,2) = s(i,2) + sig(2,l)*xg
            s(i,3) = s(i,3) + sig(3,l)*xg
            s(i,4) = s(i,4) + sig(4,l)*xg
            s(i,5) = s(i,5) + sig(5,l)*xg
            s(i,6) = s(i,6) + sig(6,l)*xg

c           Strain projections

            s(i, 7) = s(i, 7) + eps(1,l)*xg
            s(i, 8) = s(i, 8) + eps(2,l)*xg
            s(i, 9) = s(i, 9) + eps(3,l)*xg
            s(i,10) = s(i,10) + eps(4,l)*xg
            s(i,11) = s(i,11) + eps(5,l)*xg
            s(i,12) = s(i,12) + eps(6,l)*xg

c           Error estimation projection

            se(i)  = se(i)  + erav*xg

          end do ! i
        end do ! l
      endif

      iste = 12

c     Do history plots if required

      if(hpltfl) then
        call hlcn2d(hr(np(304)),nel)
      end if

      end

      subroutine meanx(ix,xl,ndm)

      implicit   none
      integer    ndm, ix(*)
      real*8     xl(ndm,*)

c     Check midside of edges

      if(ix(5).eq.0) then
        xl(1,5) = 0.5*(xl(1,1)+xl(1,2))
        xl(2,5) = 0.5*(xl(2,1)+xl(2,2))
      endif
      if(ix(6).eq.0) then
        xl(1,6) = 0.5*(xl(1,2)+xl(1,3))
        xl(2,6) = 0.5*(xl(2,2)+xl(2,3))
      endif
      if(ix(7).eq.0) then
        xl(1,7) = 0.5*(xl(1,3)+xl(1,4))
        xl(2,7) = 0.5*(xl(2,3)+xl(2,4))
      endif
      if(ix(8).eq.0) then
        xl(1,8) = 0.5*(xl(1,4)+xl(1,1))
        xl(2,8) = 0.5*(xl(2,4)+xl(2,1))
      endif

c     Compute center node location

      xl(1,9) = 0.50d0*(xl(1,5)+xl(1,6)+xl(1,7)+xl(1,8))
     &        - 0.25d0*(xl(1,1)+xl(1,2)+xl(1,3)+xl(1,4))
      xl(2,9) = 0.50d0*(xl(2,5)+xl(2,6)+xl(2,7)+xl(2,8))
     &        - 0.25d0*(xl(2,1)+xl(2,2)+xl(2,3)+xl(2,4))
      end
