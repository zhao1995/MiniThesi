c$Id:$
      subroutine psetid(eq,id,ix,ip,ie,iedof,nie,ndf,nen,nen1,
     &                  numel,numnp,nummat)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Separate eq and id from old id numbering         29/04/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Remove unused dof's using ie(nie,*) array

c      Inputs:
c        eq(ndf,numnp) - Active equation numbers
c        id(ndf,numnp) - Boundary condition codes

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,ii,j,ma,n,nen,nen1,nie,ndf,numnp,numel,nummat
      integer   eq(ndf,*), id(ndf,*),ie(nie,*),iedof(ndf,nen,*)
      integer   ip(ndf,*),ix(nen1,*)

      save

      do n = 1,numnp
        do j = 1,ndf
          ip(j,n)   = 0
          eq(j,n) = id(j,n)
        end do ! j
      end do ! n

c     Check nodes on each element for active dof's

      do n = 1,numel

        if(ix(nen1-1,n).ge.0) then

c         Loop over the material sets

          do ma = 1,nummat
            if(ie(nie-2,ma).eq.ix(nen1,n)) then
              do i = 1,nen
                ii = ix(i,n)
                if(ii.gt.0) then
                  do j = 1,ndf
                    if(iedof(j,i,ma).gt.0) then
                      ip(iedof(j,i,ma),ii) = 1
                    endif
                  end do ! j
                endif
              end do ! i
            endif
          end do ! ma
        endif
      end do ! n

c     Set b.c. restraints for unused dof's

      do n = 1,numnp
        do j = 1,ndf
          if(ip(j,n).eq.0) then
            eq(j,n) = 1
          end if
        end do ! j
      end do ! n

      end
