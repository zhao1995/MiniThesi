c$Id:$
      subroutine ptrans(ia,angl,ul,p,s,nel,ndf,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set transformation data for element computations
c               with sloping boundary conditions

c      Inputs:
c         ia(*)     - Degrees of freedom to rotate
c         angl(*)   - Array of element nodal angles
c         nel       - Number of nodes on element
c         ndf       - Number dof/node
c         nst       - Dimension of element arrays
c         isw       - Switch: rotate ul if isw=1; otherwise element
c                     arrays

c      Outputs:
c         ul(*)     - Element solution variables
c         p(*)      - Element vector
c         s(*,*)    - Element matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'

      integer   nel,ndf,nst,isw, i1,ij1,ij2, i, j, ia(2)
      real*8    cs,sn,tm, angl(*),ul(ndf,nen,4),p(ndf,*),s(nst,nst)

      save

c     Transform displacement quantities to element coordinates

      if(ndf.le.1) then
        return
      elseif(isw.eq.1) then
        do i = 1,nel
          if(angl(i).ne.0.0d0) then
            call pdegree(angl(i), sn,cs)
            do j = 1,6
              tm            = cs*ul(ia(1),i,j) - sn*ul(ia(2),i,j)
              ul(ia(2),i,j) = sn*ul(ia(1),i,j) + cs*ul(ia(2),i,j)
              ul(ia(1),i,j) = tm
            end do ! j
          endif
        end do ! i

c     Transform element arrays to global coordinates

      else
        i1  = 0
        ij1 = ia(1)
        ij2 = ia(2)
        do i = 1,nel
          if(angl(i).ne.0.0d0) then

c           Transform load vector

            call pdegree(angl(i), sn,cs)
            tm       = cs*p(ij1,i) + sn*p(ij2,i)
            p(ij2,i) =-sn*p(ij1,i) + cs*p(ij2,i)
            p(ij1,i) = tm
            if(isw.eq.2) then

c             Postmultiply s by transformation

              do j = 1,nst
                tm         = s(j,i1+ij1)*cs + s(j,i1+ij2)*sn
                s(j,i1+ij2)=-s(j,i1+ij1)*sn + s(j,i1+ij2)*cs
                s(j,i1+ij1)= tm
              end do ! j

c             Premultiply s by transformation

              do j = 1,nst
                tm         = cs*s(i1+ij1,j) + sn*s(i1+ij2,j)
                s(i1+ij2,j)=-sn*s(i1+ij1,j) + cs*s(i1+ij2,j)
                s(i1+ij1,j)= tm
              end do ! j
            endif
          endif
          i1 = i1 + ndf
        end do ! i
      endif

      end
