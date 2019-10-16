c$Id:$
      subroutine setrlk(id,rlink,rixt,ndf,numnp,neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Reset equation numbers for master/slave effects

c      Inputs:
c        id(ndf,*)    - Equation numbers before master/slave
c        rixt(*)      - List of master nodes
c        rlink(ndf,*) - Linking indicators
c        ndf          - Number dof/node
c        numnp        - Number of nodes in mesh

c      Outputs:
c        id(ndf,*)    - Equation numbers after master/slave
c        neq          - Number of active equations in problem
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndf,numnp,neq, i,j, m,n, master,slave
      integer   id(ndf,numnp),rlink(ndf,numnp),rixt(*)

      save

      do n = 1,numnp
        master = rixt(n)
        if(master.gt.0) then
          do i = 1,ndf
            if(rlink(i,n).eq.0 .and. id(i,n).gt.0) then
              slave   = id(i,n)
              id(i,n) = id(i,master)
              do m = 1,numnp
                do j = 1,ndf
                  if(id(j,m).gt.slave) then
                    id(j,m) = id(j,m) - 1
                  endif
                end do ! j
              end do ! m
            end if
          end do ! i
        end if
      end do ! n

c     Check for maximum equation number

      neq = 0
      do n = 1,numnp
        do i = 1,ndf
          neq = max(neq,id(i,n))
        end do ! i
      end do ! n

      end
