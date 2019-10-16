c$Id:$
      subroutine elimrbm(eval,evec,neq,mf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Eliminate rigid body modes for appending to finite
c               motion rigid body.

c      Inputs:
c          eval(mf)     - Eigenvalues
c          evec(neq,mf) - Eigenvectors
c          neq          - Number active dof for eigenproblem
c          mf           - Number of modes (including RB's)

c      Outputs:
c          eval(mf)     - Eigenvalues (zeros removed)
c          evec(neq,mf) - Eigenvectors(zero vectors removed
c          mf           - Number of modes (without RB's)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   neq,mf, i,j,n
      real*8    eval(*),evec(neq,*),maxeval

      save

c     Purpose: Eliminates rigid body modes

      maxeval = 0.d0

      do j = 1,mf
        maxeval = max(maxeval,abs(eval(j)))
      end do ! j
      maxeval = maxeval*1.d-6

      j = 1
100   continue

      if(abs(eval(j)).lt.maxeval) then
        do i = j+1,mf
          eval(i-1) = eval(i)
          do n = 1,neq
            evec(n,i-1) = evec(n,i)
          end do ! n
        end do ! i
        mf = mf - 1
      else
        j = j + 1
      endif

      if(j.lt.mf) then
        go to 100
      elseif(j.eq.mf) then
        if(abs(eval(j)).lt.maxeval) then
          mf = mf - 1
        end if
      endif

      end
