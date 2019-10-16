c$Id:$
      subroutine pnodcn(ix,iext,ic,ielc,indc,ndm,nen,nen1,numnp,in2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Find nodes for smoothing

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   flag
      integer   ndm,nen,nen1,numnp, i,j,k,n,nn,nel, ic1,ic2,in1,in2
      integer   ix(nen1,*),iext(*),ic(*),ielc(*),indc(*), ie(0:30)

      save

      ic1 = 0
      in1 = 0
      in2 = 0
      do n = 1,numnp
        ic2 = ic(n)
        if(iext(n).eq.0) then
          do i = ic1+1,ic2

c           Do only adjacent nodes for 2-d

            if(ndm.eq.2) then

              nel = 0
              do j = 1,nen
                if(ix(j,ielc(i)).gt.0) then
                  nel = nel + 1
                  ie(nel) = ix(nel,ielc(i))
                endif
              end do ! j
              ie(0)     = ie(nel)
              ie(nel+1) = ie(1)

              do j = 1,nel
                if(ie(j).eq.n) then

c                 Previous node

                  nn = ie(j-1)
                  flag = .true.
                  do k = in1+1,in2
                    if(nn.eq.indc(k)) flag = .false.
                  end do ! k
                  if(flag) then
                    in2 = in2 + 1
                    indc(in2) = nn
                  endif

c                 Next node

                  nn   = ie(j+1)
                  flag = .true.
                  do k = in1+1,in2
                    if(nn.eq.indc(k)) flag = .false.
                  end do ! k
                  if(flag) then
                    in2 = in2 + 1
                    indc(in2) = nn
                  endif

                endif
              end do ! j

c           Do all nodes

            else

              do j = 1,nen
                nn = ix(j,ielc(i))
                if(nn.gt.0 .and. nn.ne.n) then
                  flag = .true.
                  do k = in1+1,in2
                    if(nn.eq.indc(k)) flag = .false.
                  end do ! k
                  if(flag) then
                    in2 = in2 + 1
                    indc(in2) = nn
                  endif
                endif
              end do ! j

            endif

          end do ! i
        endif
        ic(n) = in2
        ic1   = ic2
        in1   = in2
      end do ! n

      end
