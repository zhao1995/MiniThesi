c$Id:$
      subroutine pclink(id,ip,ndf,numnp,neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:

c      Inputs:
c         id(*)  - Number of equations
c         ip(*)  - Pointer array
c         ndf    - Number dof/node
c         numnp  - Number of nodes
c         neq    - Number of active equations

c      Outputs:
c         id(*)  - Number of equations after links
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'mxsiz.h'
      include   'part0.h'

      logical    errck
      integer    ndf,numnp,neq, nmax,n1,n2,i,ii,j
      integer    id(ndf,numnp),ip(numnp)

      save

      do n1 = 1,numnp-1
        if(ip(n1).ne.n1) then
          do n2 = 1,numnp
            if(ip(n1).eq.n2) then
              do j = 1,ndf
                if(ndfp(j).eq.npart .and. clnk(j).eq.0) then
                  if(id(j,n1).gt.0 .and. id(j,n2).gt.0) then

c                   Select node to renumber dof

                    if(id(j,n1).lt.id(j,n2)) then
                      nmax     = id(j,n2)
                      id(j,n2) = id(j,n1)
                    else
                      nmax     = id(j,n1)
                      id(j,n1) = id(j,n2)
                    endif
                    do ii = 1,numnp
                      if(id(j,ii).eq.nmax) then
                        id(j,ii) = id(j,n1)
                      end if
                    end do ! ii

c                   Loop through all nodes to reduce equation numbers

                    errck = .false.
                    do i = 1,ndf
                      if(ndfp(i).eq.npart) then
                        do ii = 1,numnp
                          if(id(i,ii).gt.nmax) then
                            id(i,ii) = id(i,ii) - 1
                            errck    = .true.
                          endif
                        end do ! ii
                      endif
                    end do ! i
                    if(errck) neq = neq - 1

                  endif

                end if
              end do ! j
            endif
          end do ! n2
        endif
      end do ! n1

      end
