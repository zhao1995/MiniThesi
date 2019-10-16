c$Id:$
      subroutine plfacz(id,ix,ia,numnp,numel,ie,nie,ifc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Determine pointer array for number of faces attached
c               to each node.

c      Inputs:
c         ix(nen1,*)- Element nodal connections
c         ia(*)     - Active plot elements
c         numnp     - Number of nodes
c         numel     - Number of elements
c         ie(nie,*) - Assembly data for material sets
c         nie       - Dimension of ie array

c      Outputs:
c         id(*)     - Pointer to faces for each node
c         ifc       - Number of entries in id array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pbody.h'
      include  'pdata5.h'
      include  'pdata6.h'
      include  'sdata.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   lclip
      integer   numnp,numel,nie,ifc, i,j,k,l,m,n
      integer   ii,jj,kk,ll,iel,iiel
      integer   ix(nen1,numel), ia(*), id(*), il(4,6), ie(nie,*)

      save

      data il/3,2,1,4, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 5,6,7,8/

c     Compute location of faces which face forward

      call pzeroi(id,numnp)

c     Compute location of boundary faces by eliminating repeated faces

      do n = 1,numel

c       Set for continuum element

        if(ix(nen1-1,n).ge.nreg1 .and. ix(nen1-1,n).le.nreg2 .and.
     &     ix(nen1-1,n).ge.0 .and. ia(n).ge.0) then
         iel = ie(nie-1,ix(nen1,n))
         if(iel.gt.0) then
           iiel = inord(iel)
         else
           iiel = exord(-iel)
         endif
         if (iiel .gt. 10 ) then
          if( lclip(ix(1,n),8,hr(npxx),ndm) ) then
            do m = 1,6
              i = 1
              do j = 2,4
                if(ix(il(j,m),n).lt.ix(il(i,m),n)) i = j
              end do ! j
              j = mod(i,4) + 1
              k = mod(j,4) + 1
              l = mod(k,4) + 1
              ii = ix(il(i,m),n)
              jj = min(ix(il(j,m),n),ix(il(l,m),n))
              kk = ix(il(k,m),n)
              ll = max(ix(il(j,m),n),ix(il(l,m),n))
              if(ii.ne.kk .and. jj.ne.ll) then
                if(kk.ne.ll) then
                  kk = 32768*min(kk,ll) + max(kk,ll)
                elseif(ii.ne.jj) then
                  kk = 32768*min(kk,jj) + max(kk,jj)
                else
                  kk = 0
                end if

c               Add to list if new

                if(kk.gt.0) then
                  id(ii) = id(ii) + 1
                end if
              end if

            end do ! m
          end if
         end if
        end if

      end do ! n

c     Accumlate pointer

      j     = 1
      i     = id(1)
      id(1) = 1
      do n = 2,numnp+1
        j     = j + i
        i     = id(n)
        id(n) = j
      end do ! n
      ifc = j

      end
