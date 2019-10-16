c$Id:$
      subroutine chkrot(ie,ix,nie,nen1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Check for rotation type and auto sets.

c      Inputs:
c         ie(nie,*) - Element control information
c         ix(nen1,*) - Element connection list
c         nie        - Dimension of ie array
c         nen1       - Dimension of ix array

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'cdata.h'
      include 'crotas.h'
      include 'pointer.h'
      include 'comblk.h'

      logical  setvar,palloc
      integer  i,m,n,etyp, nie,nen1
      integer  ie(nie,*),ix(nen1,*)

      save

c     Check for nonzero rotational updates

      do n = 1,numel
        etyp = ix(nen1,n)
        do  m = 1,nummat
          if(ie(nie-2,m).eq.etyp .and. ie(nie-6,m).ne.0) then
            go to 100
          end if
        end do ! m
      end do ! n
      return

c     Allocate storage for rotational updates

100   if(np(81).eq.0) then
        frotas = .true.
        setvar = palloc(81,'MO   ',numnp*2 ,1)
        setvar = palloc(83,'MT   ',numnp   ,2)
        setvar = palloc(82,'MR   ',numnp*54,2)
      endif

c     Set rotational update for unassigned nodes

      do n = 1,numel
        etyp = ix(nen1,n)
        do  m = 1,nummat
          if(ie(nie-2,m).eq.etyp) then
            do i = 1,nen
              if(ix(i,n).gt.0) then
                if(mr(np(81)+ix(i,n)-1).eq.0) then
                   mr(np(81)+ix(i,n)-1) = ie(nie-6,m)
                endif
              endif
            end do ! i
          end if
        end do ! m
      end do ! n

      end
