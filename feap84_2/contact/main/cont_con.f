c$Id:$
      subroutine cont_con(ix1, inod, dnope,nope, neps, numnp, nod, nte)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Determine facets connected to each control point

c      Inputs :
c        ix1(:,:)  - List of master facets
c        dnope     - Dimension for ix1
c        nope      - Number nodes/facet
c        neps      - Number of facets
c        numnp     - Number of nodes/control points

c      Working:
c        inod(:)   - Node length variable

c      Output:
c        nod       - Number of unique nodes on master surface
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'pointer.h'
      include   'comblk.h'

      include   'iofile.h'

      integer    dnope, nope, neps, numnp, nod
      integer    ix1(dnope,neps), inod(numnp)

      logical    setvar, palloc
      integer    n,i, nte

c     Zero node list

      inod = 0

c     Determine master surface nodes and number times on facet

      do n = 1,neps
        do i = 1,nope
          inod(ix1(i,n)) = inod(ix1(i,n)) + 1
        end do ! i
      end do ! n

c     Count entries for data

      nod = 0
      nte = 0
      do n = 1,numnp
        if(inod(n).gt.0) then
          nod = nod + 1
          nte = nte + inod(n)
        endif
      end do ! n

c     Allocate number of nodes and total connections

      setvar = palloc(137,'CTEM2',nod  , 1)
      setvar = palloc(138,'CTEM3',nod+1, 1)
      setvar = palloc(139,'CTEM4',nte  , 1)

      call cont_list(ix1, inod, dnope, nope, neps, numnp,
     &                    mr(np(137)),mr(np(138)),mr(np(139)),nod,nte)

      end

      subroutine cont_list(ix1,inod, dnope, nope, neps, numnp,
     &                          cp_nd,cp_pt,cp_el,nod,nte)

      implicit   none

      integer    dnope, nope, neps, numnp, nod, nte
      integer    ix1(dnope,neps), inod(numnp)
      integer    cp_nd(nod), cp_pt(0:nod),cp_el(nte)

      integer    nd,n,i,p

c     Set list of master facet nodes in array; establish pointer

      cp_nd = 0
      cp_pt = 0
      nd    = 0
      do n = 1,numnp
        if(inod(n).gt.0) then
          nd = nd + 1
          cp_nd(nd) = n
          cp_pt(nd) = inod(n)
          inod(n)   = nd
        endif
      end do ! n

c     Check that all nodes are found

      if(nd.ne.nod) then
        write(*,*) ' --> ERROR: ND =',nd,' NOD =',nod
      endif

c     Convert nd_pt to ptr

      do n = 1,nod
        cp_pt(n) = cp_pt(n) + cp_pt(n-1)
      end do ! n

c     Load facets into array

      cp_el = 0
      do n = 1,neps
        do i = 1,nope
          if(ix1(i,n).gt.0) then
            nd = inod(ix1(i,n))
            do p = cp_pt(nd-1)+1,cp_pt(nd)
              if(cp_el(p).eq.0) then
                cp_el(p) = n
                exit
              endif
            end do ! p
          endif
        end do ! i
      end do ! n

      end
