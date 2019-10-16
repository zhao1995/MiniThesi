c$Id:$
      subroutine uasblem(s,p,ld,ns,afl,bfl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Added petscksp.h                                 01/05/2007
c       2. Change 'include/finclude' to 'finclude'          23/01/2009
c       3. Update petsc includes to v3.1                    20/07/2010
c       4. Remove 'isw', assembly moved to parform.F        24/03/2012
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: User function to assemble current mass matrix/vector

c     Inputs:
c       s(ns,ns)  - element consistent mass (matrix)
c       p(ns)     - element lumped mass     (vector)
c       ld(ns,2)  - local/global active equation numbers
c                   1 = rows; 2 = columns.
c       ns        - size of arrays
c       afl       - Assemble s(ns,ns) into global storage: Mmat
c       bfl       - Assemble p(ns)    into global storage: Mdiag

c     Outputs:
c                 - Assembled Mass matrix/vector (in Mmat & Mdiag)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit    none

#     include    "finclude/petsc.h"
#     include    "fdata.h"
#     include    "pfeapb.h"
#     include    "pfeapc.h"

      PetscInt    ierr
      PetscViewer M_view

      logical     afl,bfl
      integer     i, ns, n, m, j
      integer     ld(ns,*)
      real*8      s(ns,ns),p(ns),val

c     Assemble consistent mass matrix

      if(afl) then
        do m = 1,nar
          if(ld(m,5).eq.0) then       ! Assemble active rows only
            i = ld(m,1) - 1
            do n = ld(m,3),ld(m,4)
              if(ld(n,5).eq.0) then   ! Assemble active cols only
                j = ld(n,2) - 1
                val = s(m,n)
                if(val.ne.0.0d0) then
                  call MatSetValue( Mmat, i, j, val, ADD_VALUES, ierr)
                endif
              endif
            end do ! m
          endif
        end do ! n

      endif

c     Assemble lumped mass vector

      if(bfl) then
        val = 0.0d0
        do n = 1, nar
          if(ld(n,5).eq.0) then
            i = ld(n,1) - 1
            call VecSetValue( Mdiag, i, p(n), ADD_VALUES, ierr)
          else
            i = ld(n,1) - 1
            call VecSetValue( Mdiag, i, val , ADD_VALUES, ierr)
          endif
        end do ! n
      endif

      end
