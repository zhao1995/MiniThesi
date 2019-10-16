c$Id:$
      subroutine crel01 (nope,dnope,ics,emax,nelmn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Assign length to cc of 4 characters              17/04/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact Read ELEMents

c      Purpose: Input of contact surface facet nodes
c                with generation of missing data

c      Inputs:
c         nope    - # of NOde Per Element
c         dnope   - Dimension of NOde Per Element

c      Outputs:
c         ics(*)  - Contact element nodal connection array
c         emax    - Element max number found
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_dict.h'
      include  'iofile.h'
      include  'print.h'

      logical   errck,gettd,gettxtd,pcomp,whfl,whf1
      character tx(2)*15,cc*4
      integer   nope,dnope,ics(dnope,*),emax
      integer   nelo,incno,nel,nelmn,ng,kn,ke,ks,scp
      real*8    td(16)

      save

c     Check sub-command option and subcommand data not needed

      call cdebug0 ('      crel01',-1)

c     Set defaults

      emax  = 0
      nelo  = 0
      incno = 0

      errck = gettd(td,nope+2,'skip')
      nel   = nint(td(1))
      if(nelmn.eq.0) then
        nelmn = nel
      else
        nelmn = min(nelmn,nel)
      endif

c    Storage

      whfl = .true.
      do while (whfl)
        do kn = 1,nope
          ics(kn,nel) = nint(td(kn+2))
        end do ! kn

c       Automatic generation

        if((nel.gt.(nelo+1)) .and. (nelo.ne.0) .and.
     &     (incno.ne.0)                             ) then
          ng = nel-nelo-1
          do ke = 1,ng
            do kn = 1,nope
              ics(kn,nelo+ke) = ics(kn,nelo)+incno*ke
            end do ! kn
            nelmn = min(nelmn,nelo+ke)
          end do ! ke
        endif

c       Update values

        nelo  = nel
        incno = nint(td(2))
        emax  = max(emax,nel)

c       Check new data line

        errck = gettxtd(tx(1),1,td,0,'skip')
        cc    = tx(1)(1:4)

c       Stop if new string or blank line found

        if (pcomp(cc,'    ',4)) then
          whfl = .false.

c       Stop if another sub-command found

        else
          ks = 0
          whf1 = .true.
          do while (whf1)
            ks = ks+1
            if (ks.gt.nsc(1)) then

c             Sub-command not found, read line as facet data

              whf1  = .false.
              backspace (ior)
              errck = gettd(td,nope+2,'skip')
              nel   = nint(td(1))
              nelmn = min(nelmn,nel)

c           Sub-command found

            elseif (pcomp(cc,cis(scp(1,ks)),4))then
              whfl = .false.
              whf1 = .false.
              backspace (ior)
            endif
          end do ! while
        endif
      end do ! while

c     Output

      if (prt) write (iow,2001)

c     Formats

2001  format ('         ELEM:      Manual Input'/)

      end
