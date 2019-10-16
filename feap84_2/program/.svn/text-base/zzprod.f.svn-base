c$Id:$
      subroutine zzprod(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct allocation of NDNS by adding *           17/05/2007
c       2. Add set of storage for HSELM and HDNP hist plots 05/01/2012
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Driver routine for Zienkiewicz-Zhu SPR projection

c      Inputs:
c         lct     - Character array for options
c                   'off' - return to lumped nodal projection
c         ctl     - Real parameter for material number to project

c      Outputs:
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'bdata.h'
      include   'cdata.h'
      include   'eldata.h'
      include   'eldatp.h'
      include   'fdata.h'
      include   'iofile.h'
      include   'pbody.h'
      include   'pdata3.h'
      include   'prstrs.h'
      include   'sdata.h'
      include   'strnum.h'

      include   'pointer.h'
      include   'comblk.h'

      logical    pcomp,setvar,palloc
      character  lct*15
      integer    numst,ipmax
      real*8     ctl

      save

c     Reset to lumped projection

      if(pcomp(lct,'off ',4)) then

        fl(11) = .false.
        if(ior.lt.0) then
          write(*,2000)
        endif
        write(iow,2000)

c     Assemble patches for ZZ-projection

      else

c       Set parameters

        numst = npstr - 1
        ma    = nint(ctl)
        maplt = ma
        if(ma.eq.0) then
          if(ior.lt.0) then
            write(*,2001)
          endif
          write(iow,2001)
        else
          if(ior.lt.0) then
            write(*,2002) ma
          endif
          write(iow,2002) ma
        endif

c       Set arrays to find patches

        if (plfl) then
          setvar = palloc( 58,'NDNP ',numnp*npstr, 2)
          setvar = palloc( 57,'NDER ',numnp*8    , 2)
          setvar = palloc( 60,'NDNS ',max(nen*npstr,nst*nst),2)
          setvar = palloc(207,'NSCR ',numel      , 2)
          nper   = np(57)
          npnp   = np(58)
          plfl = .false.
          if(histpltfl) then
            setvar = palloc(304,'HSELM',nen*hplmax  ,2)
            setvar = palloc(305,'HDNP ',numnp*hplmax,2)
          endif
        endif
        nph  = npnp
        ner  = nper

        setvar = palloc(218,'ZZIB ',numnp    , 1)
        setvar = palloc(219,'ZZIP ',numnp+1  , 1)

c       Call routines to compute Zienkiewicz-Zhu projection

        call zzpro1(mr(np(33)),mr(np(218)),mr(np(219)),
     &              ma,ndm,nen,nen1,numnp,numel,ipmax)

        setvar = palloc(119,'TEMP9 ',ipmax, 1)

        call zzpro2(mr(np(33)),mr(np(218)),mr(np(219)),mr(np(119)),
     &              hr(np(43)),hr(nph+numnp),hr(ner+numnp),
     &              hr(np(305)),ndm,nen,nen1,numnp,numel,numst)

        setvar = palloc(119,'TEMP9 ',0, 1)

        fl(11) = .true.

      endif

c     Formats

2000  format(/'      Disable ZZ-projections - Use lumped method'/)
2001  format(/'      Perform ZZ-projections for all materials'/)
2002  format(/'      Perform ZZ-projections for material',i4/)

      end
