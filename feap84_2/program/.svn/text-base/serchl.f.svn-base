c$Id:$
      subroutine serchl(gtol,id,prsd,pu,d,stol,t,neq, step)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Line search driver program

c      Inputs:
c         gtol      - Line search reference value
c         id(*)     - Equation numbers for each dof
c         prsd      - Pointer to current residual
c         pu        - Pointer to current solution
c         stol      - Line search convergence tolerance
c         t(*)      - Working solution storage
c         neq       - Number of equations

c      Outputs:
c         d(*)      - Line search vector * step
c         step      - Step size for line search
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'counts.h'
      include  'endata.h'
      include  'iofile.h'
      include  'hdatam.h'
      include  'print.h'

      include  'p_serchl.h'

      integer   j, neq,linmax, id(*)
      real*8    gamma1, g, gtol,stol,step, sa,sb, ga,gb, d(*),t(*)

      save

      data      linmax /10/

c     Set history update flag false for line searches only

      hflgu  = .false.
      h3flgu = .false.

c     Line search in direction 'd' and return step size in 'step'

      sa = 1.0d0
      ga = gamma1(id,pu,prsd,d,t,sa)
      sb = 0.0d0
      gb = aengy

c     Find bracket on zero

      if(ga*gb.gt.0.0d0) then
         if(prnt) then
           write(iow,3000)
           if(ior.lt.0) write(*,3000)
         endif

c     Perform 'linmax' steps of line-search

      else
        j = 0
10      j = j + 1
        if (j.le.linmax) then

          step = sa - ga*(sa-sb)/(ga-gb)
          g    = gamma1(id,pu,prsd,d,t,step)

c         Output line-search parameters

          if(prnt) then
            write(iow,3001) j,step,g,sa,sb
            if(ior.lt.0) write(*,3001) j,step,g,sa,sb
          endif

c         Update postions for next iteration

          gb = 0.5d0*gb
          if (g*ga.lt.0.0d0) then
            sb = sa
            gb = ga
          endif
          sa = step
          ga = g

c         Check convergence

          if(abs(g).gt.stol*abs(gtol)) go to 10
        endif

c       Update counter for RHS forms

        nform = nform + iform
        iform = 0

c       Multiply solution increment by step size computed

        do j = 1,neq
          d(j) = step*d(j)
        end do ! j

c       Update histories

        hflgu  = .true.
        h3flgu = .false.
        g      = gamma1(id,pu,prsd,d,t,1.d0)

      endif

c     Formats

3000  format(' -> No line search - Energy positive at full step.')

3001  format(' -> Iteration',i3,' Step Size =',1p,e12.5,
     &        ' Energy =',1p,e12.5,' sa=',0p,f6.3,' sb=',0p,f6.3)

      end
