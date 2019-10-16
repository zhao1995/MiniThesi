c$Id:$
      subroutine pactivate(lct,ct,iact)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    11/07/2007
c       1. Separate equation and b.c. on psetid call        29/04/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Controls activation and deactivation of regions

c      Inputs:
c         lct       = Type of activation
c         ct(3)     = Command parameters
c         iact      = 1: Activate; -1: Deactivate region

c      Outputs:
c         none      = All through common blocks
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'comblk.h'
      include   'elacts.h'
      include   'fdata.h'
      include   'hdatam.h'
      include   'idptr.h'
      include   'iofile.h'
      include   'pointer.h'
      include   'region.h'
      include   'sdata.h'

      logical    pcomp, palloc, setvar
      character  lct*15
      integer    iact, k1,k2,k3, i,nnn
      real*8     ct(3)

      pltact = .true.

c     Set all regions to be active/deactive

      if(pcomp(lct,'all',3)) then
        k1 = 1
        k2 = mxreg
        k3 = 1

c     Set specified region only

      else
        k1 = nint(ct(1))

c       Output status of active/deactive regions

        if(k1.le.0) then
          do i = 1,mxreg
            do nnn = nen1-2,numel*nen1,nen1
              if(mr(np(33)+nnn).eq.i) then
                write(iow,2003) i
                if(ior.lt.0) write(*,2003) i
                go to 100
              elseif(-mr(np(33)+nnn).eq.i) then
                write(iow,2004) i
                if(ior.lt.0) write(*,2004) i
                go to 100
              end if
            end do ! nnn
            write(iow,2005) i
            if(ior.lt.0) write(*,2005) i
100         continue
          end do ! i
          return

c       Error: Specified region too large

        elseif(k1.gt.mxreg) then
          write(ilg,3000) k1
          write(iow,3000) k1
          if(ior.lt.0) write(*,3000) k1
          return

c       Set region range to active/deactivate

        else
          k2 = max(k1,min(nint(ct(2)),mxreg))
          k3 = max(1,nint(ct(3)))
        endif
      end if

c     Set 'first' activation indicators in 'ix' array to zero

      do nnn = nen1-3,numel*nen1,nen1
        mr(np(33)+nnn) = 0
      end do ! nnn

      do i = k1,k2,k3
        if(iact.eq.1) then
          if(ior.lt.0) write(*,2000) i
          write(iow,2000) i
          numact = i
          itract = 0
        else
          if(ior.lt.0) write(*,2001) i
          write(iow,2001) i
          numdea = i
          itrdea = 0
        endif

c       Set region indicators on IX array

        do nnn = nen1-2,numel*nen1,nen1
          if(abs(mr(np(33)+nnn)).eq.i) then
            if(iact.eq.1 .and. mr(np(33)+nnn).lt.0) then
              mr(np(33)+nnn-1) =  i                   ! new activation
            endif
            mr(np(33)+nnn) = iact*abs(mr(np(33)+nnn)) ! active/deactive
          end if
        end do ! nnn

c       Initialize strains in activated elements

        if(pcomp(lct,'init',4)) then
          if(ior.lt.0) write(*,2002) i
          write(iow,2002) i
          nreg   = i
          hflgu  = .false.
          h3flgu = .true.
          if(iact.eq.1) then
            call formfe(np(40),np(26),np(26),np(26),
     &                 .false.,.false.,.false.,.false.,17,1,numel,1)
          else
            call formfe(np(40),np(26),np(26),np(26),
     &                 .false.,.false.,.false.,.false.,18,1,numel,1)
          endif
        endif

      end do ! i

      nreg = -1

c     Reset ID equation numbers

      setvar = palloc(111,'TEMP1',nneq,  1)
      call psetid(mr(id31),mr(np(31)+nneq),mr(np(33)),mr(np(111)),
     &            mr(np(32)),mr(np(240)),nie,ndf,nen,nen1,
     &            numel,numnp,nummat)
      setvar = palloc(111,'TEMP1',   0,  1)

c     Set projection flag false (forces recompute of projections)

      fl(11) = .false.

c     Reset matrix profile

      call pnewpr()

c     Formats

2000  format('   Activate   region',i4)
2001  format('   Deactivate region',i4)

2002  format('   Setting initial strain for region =',i5)

2003  format(10x,'Region',i4,' activated')
2004  format(10x,'Region',i4,' deactivated')
2005  format(10x,'Region',i4,' does not exist')

3000  format(/' *ERROR* Cannot ACTIVATE/DEACTIVTE region number',i4/)

      end
