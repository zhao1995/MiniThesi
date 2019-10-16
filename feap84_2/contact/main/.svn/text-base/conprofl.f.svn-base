c$Id:$
      subroutine conprofl(jp,idl,id,ixl,ida,nnod,ndof,iad,ilm,lnod,nlag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor           May   30, 1996            1.0

c      Acronym: CONtact PRoFile with Lagrange multipliers

c      Purpose: Modify profile for active contacts

c      Inputs :
c        id(ndf,*)    - Equation numbers for nodal dof
c        ixl(*)       - Contact element node connections
c        ida(*)       - Active contact dof list
c        nnod         - Number of nodes on contact element ixl(*)
c        ndof         - Number of degree of freedoms in ida(*)
c        iad(numnp,*) - Lagrange multiplier determination
c        ilm(*)       - Lagrange multiplier equation numbers
c        lnod         - Number multiplier nodes
c        nlag         - Number of lagrange multiplier equations/node

c      Scratch:
c        idl(*)     - Store element active equations, etc.

c      Outputs:
c        jp(*)      - Row/column lengths for each equation
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'part0.h'
      include  'part1.h'
      include  'sdata.h'

      integer   i,ii,j,jj,mm,nad, nnod,ndof,lnod,nlag
      integer   idl(*),ixl(nnod),ida(ndof),ilm(lnod),id(ndf,*),jp(*)
      integer   iad(numnp,*)

      save

c     Check active contact nodes/degree of freedoms

      mm  = 0
      nad = 0

      do i = 1,nnod

        ii = ixl(i)

c       Set element profile

        if(ii.gt.0) then

          do jj = 1,ndof
            j = ida(jj)
            if( ndfp(j).eq.npart .and. id(j,ii).gt.0 ) then
              if(mm.eq.0) mm = id(j,ii)
              mm       = min(mm,id(j,ii))
              nad      = nad + 1
              idl(nad) = id(j,ii)
            end if
          end do ! jj

        endif

      end do ! i

c     Add Lagrange multiplier equations

      do i = 1,lnod
        ii = iad(ilm(i),2)
        jj = iad(ii,1) + iad(ilm(i),3) - nlag
        do j = 1,nlag
          nad      = nad + 1
          idl(nad) = jj + j
        end do ! j
      end do ! i

c     Compute column heights

      do i = 1,nad
        ii = idl(i)
        jp(ii) = max(jp(ii),ii-mm)
      end do

      end
