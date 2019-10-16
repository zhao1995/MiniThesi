c$Id:$
      subroutine pmodis(id,phi,y,base,phib,wb,ub,mp,neq,ndf,numnp,u,ud)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set modal contributions to nodal displacements

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'evdata.h'
      include  'part0.h'
      include  'prlod.h'
      include  'prld1.h'
      include  'pfeapb.h'

      integer   i,j,jj,m,n, mp,neq,ndf,numnp

      real*8    aa,bb,un
      integer   id(ndf,*),base(ndf,*)
      real*8    phi(vneq,*),y(mf,*),phib(neq,2,mp),wb(mp,3)
      real*8    ub(ndf,*), u(ndf,*),ud(ndf,numnp,*)

      save

c     Set the modal contributions to solution state

      call pzero(ud,2*ndf*numnp)

      aa = theta(3)
      bb = 1.d0 - aa
      do n = 1,numnp
        do j = 1,ndf
          if(npart.eq.ndfp(j)) then
            un = u(j,n)
            jj = id(j,n)
            if(jj.gt.0) then
              u(j,n) = 0.0d0
              do i = 1,mf
                u(j,n)    = u(j,n)    + phi(jj,i)*y(i,1)
                ud(j,n,1) = ud(j,n,1) + phi(jj,i)*y(i,2)
                ud(j,n,2) = ud(j,n,2) + phi(jj,i)*y(i,3)
              end do ! i

c             Add base mode contributions

              do m = 1,mp
                u(j,n)    = u(j,n)    + phib(jj,1,m)*wb(m,1)
                ud(j,n,1) = ud(j,n,1) + phib(jj,1,m)*wb(m,2)
                ud(j,n,2) = ud(j,n,2) + phib(jj,1,m)*wb(m,3)
              end do ! m

c           Boundary displacements

            else

              if(mp.gt.0 .and. base(j,n).gt.0) then
                u(j,n)    = ub(j,n)*wb(base(j,n),1)
                ud(j,n,1) = ub(j,n)*wb(base(j,n),2)
                ud(j,n,2) = ub(j,n)*wb(base(j,n),3)
              endif

            endif

            if(nrk.gt.0) then
              ud(j,n,nrk) = bb*un + aa*u(j,n)
            endif

          endif
        end do ! j
      end do ! n

      end
