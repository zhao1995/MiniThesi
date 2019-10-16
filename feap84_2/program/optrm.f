c$Id:$
      subroutine optrm(ix,ixc,ic,ip, nen,nen1,ncen1,numel,numcel,
     &                 numnp,nummat,nnel,nnid,nnac,nnrm,imat,imrm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Check for elements which are 'overlays' to other
c               elements.  Used for profile optimizer to avoid high
c               weight values.

c     Inputs:
c       ix(nen1,*)   - Element connection list
c       icx(ncen1,*) - Contact element connection list
c       ic(numnp)    - Pointer to element connections
c       ip(*)        - Elements connected to each node
c       nen          - Maximum number of nodes of elements
c       nen1         - First dimension for 'ix'
c       ncen1        - First dimension for 'ixc'
c       numel        - Number of elements in problems
c       numcel       - Number of contact elements in problems
c       numnp        - Number of nodes in problems
c       nummat       - Number of material sets
c       nnel(*)      - Number of nodes attached to each element
c       nnid(*)      - Number of active dof at node

c     Output:
c       nnac(*)      - Marker for overlay elements: 1 = overlay
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'

      integer    nen,nen1,ncen1,numel,numcel,numnp,nummat
      integer    i,ii, j,jj, n,ni,nj, mi,mj
      integer    ix(nen1,*),ixc(ncen1,*),ic(0:*),ip(*),nnel(*),nnid(*)
      integer    nnac(*),nnrm(*),imat(nummat),imrm(nummat)

      save

c     Loop over nodes

      do n = 1,numel+numcel
        nnac(n) = 0
      end do ! n

      do n = 1,numnp

c       Check elements attached to node 'n'

        do ii = ic(n-1)+1,ic(n)
          ni = ip(ii)

c         Check elements which have less than maximum number nodes

          if(ni.le.numel) then
            if(nnel(ni).lt.nen .and. nnac(ni).eq.0 ) then
              do jj = ic(n-1)+1,ic(n)
                nj = ip(jj)

c               Only check elements with greater/equal number nodes

                if(nj.le.numel) then

                  if(nnel(nj).ge.nnel(ni) .and. nnac(nj).eq.0
     &                                    .and. ni.ne.nj ) then
c                 if(nnel(nj).ge.nnel(ni) .and. ni.ne.nj ) then
c                 if(nnel(nj).ge.nnel(ni) .and. nnac(nj).eq.0) then

c                   Check for match of nodes

                    do i = 1,nnel(ni)
                      if(nnid(ix(i,ni)).eq.0) go to 100
                      do j = 1,nnel(nj)
                        if(ix(i,ni).eq.ix(j,nj)) go to 100
                      end do ! j

c                     No match of all nodes - exit

                      go to 300
100                   continue
                    end do ! i

c                   All nodes on element 'ni' match those on 'nj'

                    nnac(ni) = 1
                    nnrm(ni) = nj

                  endif

c               Check contact elements with greater/equal number nodes

                else

                  mj = nj - numel
                  if(nnel(nj).ge.nnel(ni) .and. nnac(nj).eq.0
     &                                    .and. ni.ne.nj ) then

c                   Check for match of nodes

                    do i = 1,nnel(ni)
                      if(nnid(ix(i,ni)).eq.0) go to 200
                      do j = 1,nnel(nj)
                        if(ix(i,ni).eq.ixc(j,mj)) go to 200
                      end do ! j

c                     No match of all nodes - exit

                      go to 300
200                   continue
                    end do ! i

c                   All nodes on element 'ni' match those on 'nj'

                    nnac(ni) = 1
                    nnrm(ni) = nj

                  endif
                endif
300             continue
              end do ! jj
            endif

c         Check contact elements with less than maximum number nodes

          else

            mi = ni - numel
            if(nnel(ni).lt.nen .and. nnac(ni).eq.0 ) then
              do jj = ic(n-1)+1,ic(n)
                nj = ip(jj)

c               Only check elements with greater/equal number nodes

                if(nj.le.numel) then

                  if(nnel(nj).ge.nnel(ni) .and. nnac(nj).eq.0
     &                                    .and. ni.ne.nj ) then

c                   Check for match of nodes

                    do i = 1,nnel(ni)
                      if(nnid(ixc(i,mi)).eq.0) go to 400
                      do j = 1,nnel(nj)
                        if(ixc(i,mi).eq.ix(j,nj)) go to 400
                      end do ! j

c                     No match of all nodes - exit

                      go to 600
400                   continue
                    end do ! i

c                   All nodes on element 'ni' match those on 'nj'

                    nnac(ni) = 1
                    nnrm(ni) = nj

                  endif

c               Check contact elements with greater/equal number nodes

                else


                  mj = nj - numel
                  if(nnel(nj).ge.nnel(ni) .and. nnac(nj).eq.0
     &                                    .and. ni.ne.nj ) then

c                   Check for match of nodes

                    do i = 1,nnel(ni)
                      if(nnid(ixc(i,mi)).eq.0) go to 500
                      do j = 1,nnel(nj)
                        if(ixc(i,mi).eq.ixc(j,mj)) go to 500
                      end do ! j

c                     No match of all nodes - exit

                      go to 600
500                   continue
                    end do ! i

c                   All nodes on element 'ni' match those on 'nj'

                    nnac(ni) = 1
                    nnrm(ni) = nj

                  endif
                endif
600             continue
              end do ! jj
            endif

          endif
        end do ! ii
      end do ! n

c     Report number of faces removed for each material set

      do n = 1,nummat
       imat(n) = 0
       imrm(n) = 0
      end do ! n
      do n = 1,numel
        i       = ix(nen1,n)
        imat(i) = imat(i) + 1
        if(nnac(n).ne.0) then
          imrm(i) = imrm(i) + 1
        endif
      end do ! n
      ii = 0
      jj = 0
      write(iow,2000)
      if(ior.lt.0) then
        write(*,2000)
      endif
      do n = 1,nummat
        if(imat(n).gt.0) then
          if(ior.lt.0) then
            write(*,2001) ' Material',n,' No. Elements',imat(n),
     &                    ' No. Overlayed',imrm(n)
          endif
          write(iow,2001) ' Material',n,' No. Elements',imat(n),
     &                    ' No. Overlayed',imrm(n)
          ii = ii + imat(n)
          jj = jj + imrm(n)
        endif
      end do ! n
      ni = 0
      do n = 1,numnp
        if(nnid(n).eq.0) then
          ni = ni + 1
        endif
      end do ! n

      if(numcel.gt.0) then
        nj = 0
        do n = numel+1,numel+numcel
          if(nnac(n).gt.0) then
            nj = nj + 1
          endif
        end do ! n
        jj = jj + nj
        if(ior.lt.0) then
          write(*,2002) ' Number contact elements   ',numcel,
     &                  ' No. Overlayed',nj
        endif
        write(iow,2002) ' Number contact elements   ',numcel,
     &                  ' No. Overlayed',nj
      endif

      if(ior.lt.0) then
        write(*,2003) ' Total Number of Elements',ii,
     &                ' No. Overlayed',jj,' No. Fixed Nodes',ni
      endif
      write(iow,2003) ' Total Number of Elements',ii,
     &                ' No. Overlayed',jj,' No. Fixed Nodes',ni

c     Output formats

2000  format( 5x,'O v e r l a y   S u m m a r y   i n   O p t i m i z e'
     &      /)
2001  format( 9x,a,i5,a,i10,a,i10)
2002  format( 9x,a,i10,a,i10)
2003  format(/9x,a,i12,a,i10/46x,a,i8)

      end
