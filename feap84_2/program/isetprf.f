c$Id:$
      subroutine isetprf(jp,idl,id,ix,ie,intel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute interface addition to profile of global arrays

c      Inputs:
c        jp(*)    - Pointer array to row/column ends of profile
c        idl(*)   - Local temporary storage vector
c        id(*)    - Equation numbers for degree of freedoms
c        ix(*)    - Global node numbers on elements
c        ie(*)    - Element properties
c        intel(*) - Interface element pairs

c      Outputs:
c        jp(*)  - Modified pointer array to row/column ends of profile
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'sdata.h'
      include   'part0.h'
      include   'pointer.h'
      include   'comblk.h'

      integer    jp(*),idl(*),id(ndf,*),ix(nen1,*),ie(nie,*),intel(8,*)
      integer    i,ii,j, ml,mm, nn, n1,n2,nad

      save

      nn = 1
      do while(intel(1,nn).ne.0)
        n1 = intel(1,nn)
        n2 = intel(2,nn)
        nn = nn + 1

c       Test for active elements

        if(n1.gt.0 .and. n2.gt.0) then

          mm  = 0
          nad = 0
          do i = 1,nen

c           Set element 1 profile

            ii = ix(i,n1)
            if(ii.gt.0) then
              do j = 1,ndf
                if( ndfp(j).eq.npart .and. id(j,ii).gt.0 ) then
                  if(mm.eq.0) mm = id(j,ii)
                  mm       = min(mm,id(j,ii))
                  nad      = nad + 1
                  idl(nad) = id(j,ii)
                end if
              end do ! j
            endif

c           Set element 2 profile

            ii = ix(i,n2)
            if(ii.gt.0) then
              do j = 1,ndf
                if( ndfp(j).eq.npart .and. id(j,ii).gt.0 ) then
                  if(mm.eq.0) mm = id(j,ii)
                  mm       = min(mm,id(j,ii))
                  nad      = nad + 1
                  idl(nad) = id(j,ii)
                end if
              end do ! j
            endif

          end do ! i

c         Element 1 Lagrange multiplier equations

          if(ie(nie-8,ix(nen1,n1)).gt.0) then
            ml = ix(nen+4,n1)
            if(ml.gt.0) then
              ml = mr(np(211)-1+ml) + ix(nen+5,n1) - 1
              do j = 1,ie(nie-8,ix(nen1,n1))
                if(mm.eq.0) mm = ml + j
                mm       = min(mm,ml + j)
                nad      = nad + 1
                idl(nad) = ml + j
              end do ! j
            endif
          endif

c         Element 2 Lagrange multiplier equations

          if(ie(nie-8,ix(nen1,n2)).gt.0) then
            ml = ix(nen+4,n2)
            if(ml.gt.0) then
              ml = mr(np(211)-1+ml) + ix(nen+5,n2) - 1
              do j = 1,ie(nie-8,ix(nen1,n2))
                if(mm.eq.0) mm = ml + j
                mm       = min(mm,ml + j)
                nad      = nad + 1
                idl(nad) = ml + j
              end do ! j
            endif
          endif

c         Compute column heights

          do i = 1,nad
            ii = idl(i)
            jp(ii) = max(jp(ii),ii-mm)
          end do ! i

        end if

      end do ! while

      end
