c$Id:$
      subroutine pmacr7 (lct,ct,j)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove argument 'prt7'                           13/12/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Command language instruction subprogram: Part 7

c      Inputs:
c         lct      - Command option for current command
c         ct(3)    - Command parameters for current command
c         j        - Number of command to execute

c      Outputs:
c         Depends on value of command j: None for serial version
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    j
      character  text*15,lct*15
      real*8     ct(3)

      save

c     Solution command instruction subprogram - part 7.

      text = lct

      if(nint(ct(1)).ge.0) then

c       [graph,<node>,d] d = number of domains
c       [graph,file]     Input graph partion from existing file

        if(j.eq.1) then

          write(*,2000) ' GRAPh '

c       [outd,<bloc>k] - Output of domain meshes

        elseif(j.eq.2) then

          write(*,2000) ' OUTDomain '

c       [pets]c,<on,off,view,novi> - Parallel solution using PETSc

        elseif(j.eq.3) then

          write(*,2000) ' PETSc '

c       Subspace eigencomputations (for: mass,iden,geom)

        elseif(j.eq.4) then

          write(*,2000) ' PSUBspace '

c       [u_in]<inpu,outp> - Initial solution inputs/outputs

        elseif(j.eq.5) then

          write(*,2000) ' U_IN '

c       [glis]t,,# -- Input list of global node numbers to output

        elseif(j.eq.6) then

          write(*,2000) ' GLISt '

c       [gplo]t,disp,#  -- Output displacement comp. # global plot ndat
c       [gplo]t,stre,#  -- Output nodal stress comp. # global plot ndat
c       [gplo]t,pstr,#  -- Output nodal stress comp. # global plot ndat

        elseif(j.eq.7) then

          write(*,2000) ' GPLOt '

        endif
      endif

c     Formats

2000  format('  *ERROR*',a,'command only available in parallel version')

      end
