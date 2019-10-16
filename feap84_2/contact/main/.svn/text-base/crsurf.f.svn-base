c$Id:$
      subroutine crsurf (nsurf,cs0,lcs1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Add 'nurb' and 'tspline' surface type outputs      02/08/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact Read SURfaces data

c      Purpose: Input of contact surface data

c      Inputs :
c         nsurf   - Number of SURface

c      Outputs:
c         cs0(*)  - Contact surfaces control data
c         lcs1    - Length of facet array ICS
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include   'c_0.h'
      include   'c_contac.h'
      include   'c_comnd.h'
      include   'c_dict.h'
      include   'cdata.h'
      include   'chdata.h'
      include   'iofile.h'
      include   'print.h'
      include   'p_ptname.h'
      include   'pointer.h'
      include   'sdata.h'
      include   'comblk.h'

      include   'p_int.h'

      logical    whfl,polfl, outfac, setvar,palloc, elmnf, pcomp
      character  getlab*70,label*70,type*4,cc*4,c2*4
      integer    nsurf,lcs1,nelmn, ncom,ntype,nfeat,nopti,nsubc,nsopt
      integer    emax,kr,nope,dnope,neps,kn,ke,k,labl, regn, nblk
      real*8     cs0(nr0,n0c1:*), tydat(15),td(14)

      save

      call cdebug0 ('    crsurf',-1)

c     Print command title

      if (prt) then
        write(iow,2100) nsurf
        label = getlab(2,labl)
        if (labl.ne.0) write (iow,2101) label
      endif

c     Store surface # and offset for surface facets

      ncom      = 1
      cs0(1,-1) = nsurf
      cs0(2,-1) = ofssurf

c     TYPE declaration

      call crtype (ncom,nsurf,type,ntype,tydat)
      cs0(1,0) = ntype

c     Set default number of nodes/facet

      elmnf = .false.
      nelmn = 0
      dnope = 0
      if (ntype.eq.1) then                     ! Type: line
        tydat(1) = max(tydat(1),2.d0)
        dnope    = 2
      elseif (ntype.eq.2) then                 ! Type: tria
        tydat(1) = max(tydat(1),3.d0)
      elseif (ntype.eq.3) then                 ! Type: quad
        tydat(1) = max(tydat(1),4.d0)
      elseif (ntype.eq.4) then                 ! Type: beam
        tydat(1) = max(tydat(1),2.d0)
        dnope    = 2
      elseif (ntype.eq.5) then                 ! Type: point
        tydat(1) = max(tydat(1),1.d0)
        elmnf    = .true.
      elseif (ntype.eq.6) then                 ! Type: rigid
        tydat(1) = max(tydat(1),1.d0)
      elseif (ntype.eq.7) then                 ! Type: NURBs
        tydat(1) = max(tydat(1),1.d0)
        nblk     = max(nint(tydat(2)),1)
      elseif (ntype.eq.8) then                 ! Type: T-spline
        tydat(1) = max(tydat(1),1.d0)
        nblk     = max(nint(tydat(2)),1)
      endif

c     Set value for dimensioning array

      nope      = nint(tydat(1))
      dnope     = nope + dnope
      cs0(4,-1) = dnope

c     Save data in command dictionary

      do kr = 1,15
        cs0(kr+1,0) = tydat(kr)
      end do

c     Output type of facet information

      outfac = .true.

      if (ntype.eq.1) then          ! Line facet

        if(prt) then
          write (iow,2110) type,nope
          write (iow,2111)
        endif

      elseif (ntype.eq.2) then      ! Triangular facet

        if(prt) then
          write (iow,2110) type,nope
          write (iow,2112)
        endif

      elseif (ntype.eq.3) then      ! Quadrilateral facet

        if(prt) then
          write (iow,2110) type,nope
          write (iow,2113)
        endif

      elseif (ntype.eq.4) then      ! Beam facet

        if(prt) then
          write (iow,2110) type,nope
          write (iow,2114)
        endif

      elseif (ntype.eq.5) then      ! Point facet

        if(prt) then
          write (iow,2110) type,nope
          write (iow,2115)
        endif

      elseif(ntype.eq.6) then       ! Rigid facet

        if(prt) then
          write (iow,2110) type,nope
          write (iow,2116)
        endif

      elseif(ntype.eq.7) then       ! NURBS facet

        if(prt) then
          write (iow,2110) type,nope
          write (iow,2117) 'NURBS'
        endif

      elseif(ntype.eq.8) then       ! T-spline facet

        if(prt) then
          write (iow,2110) type,nope
          write (iow,2117) 'T-spline'
        endif

      elseif (ntype.eq.9) then        ! Facets from part number

        call acheck(xxx,yyy,15,80,80)
        do kr = 1,ipart
          if(pcomp(yyy(16:30),partname(kr),15)) then
            nope = ptelm(3,kr)
            if (prt) then
              write (iow,2110) type,nope
              write (iow,2119) partname(kr),ptelm(1,kr),ptelm(2,kr)
            endif

c           Set memory for storing facets

            dnope  = nope + 1
            emax   = ptelm(2,kr) - ptelm(1,kr) + 1
            setvar = palloc(140,'CTEM5',dnope*emax,1)
            call crel09(nope,dnope,mr(np(33)),ptelm(1,kr),ptelm(2,kr),
     &                  mr(np(140)) )
c           Update and store total # of facets & dimensions for facets
            neps      = max(int(cs0(3,-1)),emax)
            cs0(3,-1) = neps
            cs0(4,-1) = dnope
            cs0(2, 0) = nope
            cs0(3, 0) = ptelm(1,kr)
            cs0(4, 0) = ptelm(2,kr)

            exit
          endif
        end do ! kr

      endif

c     Read data

      polfl = .false.
      whfl  = .true.
      do while (whfl)
        call crdata (ncom,cc,c2,nfeat,nopti,nsubc,nsopt,td)

c       FEATURE - deposit values in command table

        if (nfeat.gt.0) then
          cs0(1,nfeat) = nfeat
          cs0(2,nfeat) = nopti
          do k = 1,14
            cs0(k+2,nfeat) = td(k)
          end do
          if (prt) then
            write (iow,2120) cc
            if ((nfeat.eq.1)) then
c              write (iow,2121) c2
            endif
          endif

c       SUB-COMMAND - perform subcommand

        elseif (nsubc.gt.0) then
          if (prt) then
            write (iow,2130) cc,c2
          endif

c         Set last facet number (required for automatic generation)

          emax = nint(cs0(3,-1))

c         Facet direct reading/generation

          if (nsubc.eq.1) then
            setvar = palloc(140,'CTEM5',dnope*max(numnp,numel),1)
            call crel01 (nope,dnope,mr(np(140)),emax,nelmn)
c           Update and store total # of facets
            neps      = max(int(cs0(3,-1)),emax)
            cs0(3,-1) = neps
            if(elmnf .and. nelmn.gt.1) then
              write(iow,3003) nelmn
              write(ilg,3003) nelmn
              call plstop()
            endif

c         Block generation

          elseif (nsubc.eq.2) then
            setvar = palloc(140,'CTEM5',(dnope+3)*numel+numnp+1,1)
            call crel02 (nsopt,td,nope,dnope,mr(np(140)),emax,polfl)
c           Update and store total # of facets
            if(nsopt.eq.2) then
              neps      = max(int(cs0(3,-1)),emax)
              cs0(3,-1) = neps
            endif

c         Blend generation

          elseif (nsubc.eq.3) then
            setvar = palloc(140,'CTEM5',dnope*numel+numnp+1,1)
            call crel03 (nsopt,td,nope,dnope,mr(np(140)),emax)
c           Update and store total # of facets
            if(nsopt.eq.2) then
              neps      = max(int(cs0(3,-1)),emax)
              cs0(3,-1) = neps
            endif

c         Region generation

          elseif (nsubc.eq.4) then
            setvar = palloc(140,'CTEM5',dnope*numel,1)
            regn   = nint(td(1))
            call crel04(mr(np(140)),mr(np(33)),mr(np(150)),
     &                  regn,dnope,emax)
c           Update and store total # of facets
            neps      = max(int(cs0(3,-1)),emax)
            cs0(3,-1) = neps

c         Function generation

          elseif (nsubc.eq.5) then
            setvar = palloc(140,'CTEM5',dnope*numnp,1)
            call crfunc(nsopt,td,emax)
c           Update and store total # of facets
            neps      = max(int(cs0(3,-1)),emax)
            cs0(3,-1) = neps
            cs0(2, 0) = nsopt
            do k = 3,c_nr0
              cs0(k,0) = td(k-2)
            end do ! k
            outfac = .false.

          endif

c       TYPE DATA --> read type data
c         N.B.  Currently no type data exists

        elseif (nfeat.eq.-3) then
          if (prt) then
            write (iow,2140)
          endif
c         call cr.... ()
          write (ilg,3001) nsurf,xxx
          write (iow,3001) nsurf,xxx
          write (  *,3001) nsurf,xxx
          call coutcis(1)
          call plstop()

c       Blank line found - search again

        elseif (nfeat.eq.-1) then
          continue

c       New command found - go back to PCONT

        elseif (nfeat.eq.-2) then
          backspace (ior)
          whfl = .false.
        endif
      end do

c     Add search variables to surface facets

      if (ntype.eq.1) then                    ! Node to surface, 2d
        call csurflr(mr(np(140)),dnope,neps)
c     elseif (ntype.eq.2) then                ! Smooth node to surf, 2d
c       call csurflr(mr(np(140)),dnope,neps)
      elseif (ntype.eq.4) then                ! Beam to beam, 3d
        call csurflr(mr(np(140)),dnope,neps)
      endif

c     Output of surface facets

      if (prt .and. outfac) then
        if(dnope.gt.nope) then
          write (iow,2151) (kn,' Node',kn=1,nope),
     &                     (kn-nope,' Side',kn=nope+1,nope+2)
          do ke = 1,neps
            fp(1) = np(140) - 1 + dnope*(ke-1)
            write (iow,2152) ke,(mr(fp(1)+kn),kn=1,nope),
     &                          (mr(fp(1)+kn),kn=dnope-1,dnope)
          end do ! ke
        else
          write (iow,2151) (kn,' Node',kn=1,nope)
          do ke = 1,neps
            fp(1) = np(140) - 1 + dnope*(ke-1)
            write (iow,2152) ke,(mr(fp(1)+kn),kn=1,nope)
          end do
        endif
      endif

c     Check for multiple surfaces

      if(dnope.gt.nope) then
        ke = -1
        kr = -1
        do k = 1,neps
          fp(1) = np(140) -1 + dnope*(k-1)
          if(mr(fp(1)+dnope  ).eq.0) kr = kr + 1
          if(mr(fp(1)+dnope-1).eq.0) ke = ke + 1
        end do
        if(kr+ke.gt.0) then
          write(ilg,3002) ke,kr
          write(iow,3002) ke,kr
          write(  *,3002) ke,kr
          call plstop()
        endif
      endif

c     Compute new offset for next surface

      if(ntype.eq.6) then
        lcs1  = 0
      else
        lcs1  = dnope*neps
      endif
      ofssurf = ofssurf + lcs1

c     Formats

2100  format (/5x,'C o n t a c t   S u r f a c e   D a t a'//
     &         5x,'Data Set for Surface Number       ',i4/)

2101  format (/5x,'Surface Label: ',a)

2110  format ( 5x,'Contact Facet Type                ',a/
     &         5x,'Max # of Nodes per Facet          ',i4)
2111  format (/5x,'Two D line facet with ',
     &            '2 or more nodes')
2112  format (/5x,'Three D triangular facet with ',
     &            '3 or more nodes')
2113  format (/5x,'Three D quadrilateral facet with ',
     &            '4 or more nodes')
2114  format (/5x,'Three D beam facet with ',
     &            '2 or more nodes')
2115  format (/5x,'Three D point facet with 1 node.')

2116  format (/5x,'Rigid contact surface.')

2117  format (/5x,a,' contact surface.')

2119  format (/5x,'PART contact surface: ',a/
     &        10x,'Elements =',i8,' to',i8)

2120  format (/5x,'  Surface Feature                 ',a)

2130  format ( 5x,'  Surface Sub-command             ',a,' ',a)

2140  format ( 5x,'  Type Declaration Data:          ')

2151  format (/5x,'S u r f a c e    F a c e t    C o n n e c t i o n s'
     &       //5x,'Facet',8(i3,a5):/(10x,8(i3,a5)))
2152  format (i10,8i8:/(10x,8i8:))

3001  format (/' *ERROR* CRSURF: Reading data for SURFACE ',i5/
     &         '  Unrecognized data for command: '/
     &         1x,a)

3002  format (/' *ERROR* CRSURF: Multiple segments on surfaces.'/
     &         '         Correct by using TIE for surface area;'/
     &         '         or redefine surface data to avoid gaps.'/
     &         '         Number Side 1 Surfaces ',i5/
     &         '         Number Side 2 Surfaces ',i5/)

3003  format (/' *ERROR* CRSURF: Elements for point input do not',
     &         ' start from 1 (one).'/
     &         '         MINIMUM = ',i5)
      end
