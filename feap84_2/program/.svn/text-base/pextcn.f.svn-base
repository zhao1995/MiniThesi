c$Id:$
      subroutine pextcn(id,ic,ie,ix,ip,numnp,numel,nen1,nie)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comblk.h'

      logical   ifl,iend
      integer   numnp,numel,nen1,nie
      integer   i, j, k, ii, jj, ij, iju, n, n1, n2, ni
      integer   id(*),ic(*),ie(nie,*),ix(nen1,*),ip(*),jplt(30)

      save

c     Initialize connection array

      do i = 1,numnp
        ip(i) = 0
      end do ! i

      do i = 1,id(numnp+1)
        ic(i) = 0
      end do ! i

c     Loop through elements to set up list

      do n = 1,numel
        ii = ix(nen1,n)
        call pltord(ix(1,n),ie(nie-1,ii), iju,jplt)

c       Look up element nodes

        ii = abs(ix(jplt(1),n))
        do ij = 2,iju
          j = jplt(ij)
          jj = abs(ix(j,n))
          if(jj.ne.ii) then
            n1 = min(ii,jj)
            n2 = max(ii,jj)
            do k = id(n1),id(n1+1)-1
              if(ic(k).eq.0) then
                ic(k) =  n2
                go to 100
              elseif(abs(ic(k)).eq.n2) then
                ic(k) = -abs(n2)
                go to 100
              endif
            end do ! k
100         ii = jj
          endif
        end do ! ij
      end do ! n

c     Compute outline nodes

      do ni = 1,numnp
        iend = .true.
        do n = 1,numnp
          ifl = .true.
          n1  =  n
101       do i = id(n1),id(n1+1)-1
            if(ic(i).gt.0) then
              go to 102
            elseif(ic(i).eq.0) then
              go to 103
            endif
          end do ! i
          go to 103
102       iend = .false.
          if(ifl) then
            ip(n1) = 1
            ifl = .false.
          endif
          n2     =  ic(i)
          ic(i)  = -n2
          ip(n2) = 1
          n1     = n2
          go to 101
103       continue
        end do ! n
        if(iend) go to 104
      end do ! ni

104   continue

      end
