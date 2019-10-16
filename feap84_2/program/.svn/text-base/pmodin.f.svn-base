c$Id:$
      subroutine pmodin(phi, y, md,mu,id,u,v,ti,tu,nv,neq,tneq,mtyp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Initial Conditions for modal solution of
c               linear transient problems

c      Inputs:
c         phi(neq,*) - Mass orthonormal eigenvectors
c         md(*)      - Diagonal part of mass matrix
c         mu(*)      - Upper part of consistent mass
c         id(*)      - Pointer array to active equations
c         u(*)       - Initial displacements
c         v(*)       - Initial velocities
c         nv         - Number of eigenpairs (converged)
c         neq        - Number of active components in eigenvectors
c         tneq       - Number of d.o.f. in model
c         mtyp       - Mass type (true = consistent,false = lumped)

c      Scratch:
c         ti(*,2)    - Temporary vector
c         tu(*,2)    - Temporary vector

c      Outputs:
c         y(nv,3)    - Eigensolution at time t_0
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none


      include  'pfeapb.h'
      include  'sdata.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   mtyp
      integer   i, nv,neq,tneq, id(*)
      real*8    phi(vneq,*), y(nv,3), md(*),mu(*),u(*),v(*)
      real*8    ti(tneq,*),tu(tneq,*), dot

      save

c     Initial displacement solution

      do i = 1,tneq
        tu(i,1) = 0.0d0
        tu(i,2) = 0.0d0
      end do ! i

      do i = 1,tneq
        if(id(i).gt.0) then
          tu(id(i),1) = u(i)
          tu(id(i),2) = v(i)
        endif
      end do ! i

c     Parallel solution

      if(pfeap_on) then

        call ppmodin(phi, y, u,v,tu,nv,tneq,mtyp)

c     Serial solution

      else
c       Consistent mass

        if(mtyp) then
          do i = 1,neq
            ti(i,1) = 0.0d0
            ti(i,2) = 0.0d0
          end do ! i
          call caprod(md,mu,tu(1,1),ti(1,1),mr(np(90)),mr(np(91)),neq)
          call caprod(md,mu,tu(1,2),ti(1,2),mr(np(90)),mr(np(91)),neq)

c       Lumped mass

        else
          do i = 1,neq
            ti(i,1) = md(i)*tu(i,1)
            ti(i,2) = md(i)*tu(i,2)
          end do ! i
        end if

c       Initial motion participation

        do i = 1,nv
          y(i,1) = dot(ti(1,1),phi(1,i),neq)
          y(i,2) = dot(ti(1,2),phi(1,i),neq)
        end do ! i

      endif

c     Exit

      end
