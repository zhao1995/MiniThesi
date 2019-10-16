c$Id:$
      subroutine defdof (ixl,ida,idl,nnod,ndof,id)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: DEFine contact DOF

c      Purpose: Determine dof involved for a contact element

c      Inputs :
c         ixl(*)  - IX Local vector with involved nodes
c         ida(*)  - ID Active vector with active contact dof
c         nnod    - # of nodes in ixl
c         ndof    - # of dof in ida
c         id(*)   - Equation numbers for degree of freedoms   (iop = 1)

c      Outputs:
c         idl(*)  - Temporary array with DOF of the contact element
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_tanfl.h'
      include  'part0.h'
      include  'part1.h'
      include  'sdata.h'

      integer   nnod,ndof, kl,knod,kdof,node,dof
      integer   ixl(*),ida(*),idl(ndof,*),id(ndf,*)

      save

      call cdebug0 ('        defdof',-1)

c     Set assembly for full vector

      if(ddfl) then
        do knod = 1,nnod
          node = ixl(knod)
          kl   = ndf*node - ndf
          do kdof = 1,ndof
            dof            = ida(kdof)
            idl(kdof,knod) = kl + dof
          end do ! kdof
        end do ! knod

c     Set assembly for active part of compressed vector
c     and check for partitions

      else
        do knod = 1,nnod
          node = ixl(knod)
          do kdof = 1,ndof
            if(npart.eq.ndfp(dof)) then
              dof            = ida(kdof)
              idl(kdof,knod) = id(dof,node)
            else
              idl(kdof,knod) = 0
            endif
          end do ! kdof
        end do ! knod
      endif

      end
