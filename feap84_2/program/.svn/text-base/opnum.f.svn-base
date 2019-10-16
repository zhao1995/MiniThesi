c$Id:$
      subroutine opnum(ix,ixc,nd,ln,ne,ndw,msum,nfrnt,nnac,nnid,
     &                 numnp,numel,numcels,nen,nen1,ncen,ncen1,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Calculation of best order to number equations.
c               Numbers equations for minimum front width/profile.

c               Ref: M. Hoit and E.L. Wilson, 'An Equation
c                    Numbering Algorithm Based on a Minimum
c                    Front Criteria,' Computers & Structures,
c                    v 16, No. 1-4, pp225-239, 1983.
c               Modified: R.L. Taylor; 1 December 2000
c                         Ignores overlayed elements and fixed nodes in
c                         computing weights.

c      Inputs:
c         ix(nen1,*)     - Element nodal connection list
c         ixc(ncen1,*)   - Contact element nodal connection list
c         nnac(numel)    - Overlay markers
c         nnid(numnp)    - Number dof at each node
c         numnp          - Number of nodes in mesh
c         numel          - Number of elements in mesh
c         numcels        - Number of contact elements in mesh
c         nen            - Maximum number of nodes/element
c         nen1           - Dimension of ix  array
c         ncen           - Maximum number of nodes/contact element
c         ncen1          - Dimension of ixc array
c         prt            - Print results if true

c      Scratch:
c         nd             - Nodes currently in front
c         ln             - Pointer array for element number
c         ne             - Element numbers connected to nodes
c         ndw            - Node weights
c         msum           - Element weights

c      Outputs:
c         nfrnt(i)       - Original node for new number i
c         msum(i)        - New element number for element i
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      logical   prt
      integer   ie, l, m,minw,ml,mm,mh, n,nf,nn,nsum
      integer   numnp,numel,numcels,nen,nen1,ncen,ncen1
      integer   node,nume,nstart,numb
      integer   ix(nen1,*),ixc(ncen1,*),nd(numnp),ln(numnp),ne(*)
      integer   ndw(numnp),msum(numel),nfrnt(numnp),nnac(numel)
      integer   nnid(numnp)

      save

c     Initialization

      if(prt) write (iow,2000)
      node = 0
      nume = 0

c     Set start lists

      call nodel(ix,ixc,nnac,nnid,nd,ln,ne,numnp,numel,numcels,
     &           nsum,nen,nen1,ncen,ncen1)

c     Initialize the element weight

      do m = 1,numel + numcels
        msum(m) = 1
      end do ! m

c     Locate starting node

  100 nstart  = 0
      call nodew(ndw,msum,ix,ixc,nnac,nnid,nen,nen1,ncen,ncen1,
     &           numnp,numel,numcels,nstart)
      if(nstart.ne.0) then
        if(prt) write (iow,2004) nstart
        nfrnt(1) = nstart
        numb = 1

c       Find next element to be added to front

  110   ie = 0
        minw = 32000000

c       Loop over existing nodes on front

        do nn = 1,numb
          nf = nfrnt(nn)
          mh = ln(nf)
          ml = 1
          if(nf.ne.1) ml = ln(nf-1) + 1

c         For each node on front check elements ---

          do mm = ml,mh

            m = ne(mm)
            if(m.gt.0 .and.m.le.numel) then
              if(msum(m).gt.0 .and. ix(nen1-1,m).ge.0
     &                        .and. nnac(m).eq.0) then

c               Evaluate increase or decrease in front --

                nsum = 0
                do l = 1,nen
                  n = abs(ix(l,m))
                  if(n.gt.0) then
                    if(nnid(n).gt.0) then
                      if(ndw(n).eq.msum(m)) nsum = nsum - 1
                      if(nd(n).ge.0)        nsum = nsum + 1
                    endif
                  endif
                end do ! l

c               Compare with previous minimum

                if(nsum.lt.minw) then
                  minw = nsum
                  ie   = m
                endif
              endif
            elseif(m.gt.numel) then
              if(msum(m).gt.0) then

c               Evaluate increase or decrease in front --

                nsum = 0
                do l = 1,ncen
                  n = abs(ixc(l,m-numel))
                  if(n.gt.0) then
                    if(nnid(n).gt.0) then
                      if(ndw(n).eq.msum(m)) nsum = nsum - 1
                      if(nd(n).ge.0)               nsum = nsum + 1
                    endif
                  endif
                end do ! l

c               Compare with previous minimum

                if(nsum.lt.minw) then
                  minw = nsum
                  ie   = m
                endif
              endif
            endif
          end do ! mm

        end do ! nn

c       Subtract element sums from nodal values

        m = ie
        if(prt) write (iow,2001) m

        if(m.gt.0 .and. m.le.numel .and. nnac(m).eq.0) then
          do l = 1,nen

c           Reduce node sums by element weights

            n = abs(ix(l,m))
            if(n.gt.0) then
              if(nnid(n).gt.0) then
                ndw(n) = ndw(n) - msum(m)
                call front(nfrnt,nd,n,numb,1)

c               Check if equation is to be numbered

                if(ndw(n).eq.0) then
                  node = node + 1
                  nd(n) = node
                  call front(nfrnt,nd,n,numb,2)
                  if(prt) write (iow,2002) n
                endif

              endif
            endif
          end do ! l

        elseif(m.gt.numel) then
          do l = 1,ncen

c           Reduce node sums by element weights

            n = abs(ixc(l,m-numel))
            if(n.gt.0) then
              if(nnid(n).gt.0) then
                ndw(n) = ndw(n) - msum(m)
                call front(nfrnt,nd,n,numb,1)

c               Check if equation is to be numbered

                if(ndw(n).eq.0) then
                  node = node + 1
                  nd(n) = node
                  call front(nfrnt,nd,n,numb,2)
                  if(prt) write (iow,2002) n
                endif

              endif
            endif
          end do ! l

        endif

c       Remove element from system

        nume = nume + 1
        if(m.gt.0) then
          msum(m) = - nume
        endif

c       Check if front has been reduced to zero

        if(numb.eq.0) go to 100
        go to 110
      endif

c     Put in final order - add inactive nodes

      do n = 1,numnp

        if(nd(n).eq.0) then
          node = node + 1
          nfrnt(node) = n
          if(prt) write (iow,2003) n
        else
          nn = nd(n)
          nfrnt(nn) = n
        endif

      end do ! n

c     Formats

 2000 format(/1x,'Calculation of Minimum Front'/)

 2001 format( 5x,'Next element on front=',i6)

 2002 format( 5x,'Next node numbered =',i6)

 2003 format( 5x,'Node with no element attached =',i7)

 2004 format( 5x,'Starting node number =',i7/)

      end
