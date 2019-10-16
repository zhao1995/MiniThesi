c$Id:$
      subroutine crpair (npair,ncdim,cp0)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1.  Add set of inititial PENEtration option         03/05/2007
c           Format 2328
c       2.  Add output of axisymmetry for pairs             09/07/2010
c       3.  Add option to use time function                 14/03/2011
c       4.  Activate equations for perturbed solutions      09/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact Read SURfaces data

c      Purpose: Input of contact surface data

c      Inputs :
c         nsurf   - Number of the SURface
c         ncdim   - Number Contact DIMensions

c      Outputs:
c         cp0(*)  - Contact pair control data
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_dict.h'
      include  'c_tanfl.h'
      include  'chdata.h'
      include  'iofile.h'
      include  'print.h'
      include  'sdata.h'

      logical   whfl, pcomp, tlfl
      character getlab*70,label*70,type*4,cc*4,c2*4, cl*4
      integer   npair,ncdim
      integer   ncom,ntype,nfeat,nopti,nsubc,nsopt,kr,k,labl
      real*8    cp0(nr0,n0c3:*), tydat(15),td(14), tl(3)

      save

      call cdebug0 ('    crpair',0)

c     Print command title

      if (prt) then
        if(ncdim.eq.0) then
          write(iow,2300) npair
        else
          write(iow,2300) npair,ncdim
        endif
        label=getlab(2,labl)
        if (labl.ne.0) write (iow,2310) label
      endif

c     Store surface #

      ncom       = 3
      cp0( 1,-1) = npair
      cp0(13,-1) = ncdim

c     TYPE declaration

      call crtype (ncom,npair,type,ntype,tydat)

      if(prt) write (iow,2320) type

c     User Contact Elements 1 to 20

      if (    ntype.eq.1) then
        call cpair01(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.2) then
        call cpair02(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.3) then
        call cpair03(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.4) then
        call cpair04(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.5) then
        call cpair05(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.6) then
        call cpair06(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.7) then
        call cpair07(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.8) then
        call cpair08(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.9) then
        call cpair09(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.10) then
        call cpair10(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.11) then
        call cpair11(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.12) then
        call cpair12(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.13) then
        call cpair13(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.14) then
        call cpair14(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.15) then
        call cpair15(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.16) then
        call cpair16(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.17) then
        call cpair17(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.18) then
        call cpair18(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.19) then
        call cpair19(npair,ncdim,cp0,tydat)
      elseif (ntype.eq.20) then
        call cpair20(npair,ncdim,cp0,tydat)

c     Node to Segment Contact

      elseif (ntype.eq.21) then

        if(ndm.eq.3) then
          cp0(12,-1) = 1.0d0  ! force surface generation for sets
        endif
        if(prt) write (iow,2411) (int(tydat(kr)),kr=1,2)

c     Node to Node Contact

      elseif (ntype.eq.22 .or. ntype.eq.23) then

        if(prt) write (iow,2412) (int(tydat(kr)),kr=1,3)

c     Node to Rigid Contact

      elseif (ntype.eq.24 .or. ntype.eq.25) then

        if(prt) write (iow,2414) (int(tydat(kr)),kr=1,2)

c     Tied surface to surface

      elseif (ntype.eq.26) then

        if(prt) write (iow,2415) (int(tydat(kr)),kr=1,2)

      endif

c     Store command parameters in control table

      cp0(1,0) = ntype
      do kr = 2,c_nr0
        cp0(kr,0) = tydat(kr-1)
      end do ! kr

c     Read data

      whfl = .true.
      tlfl = .false.
      do while (whfl)
        call crdata (ncom,cc,c2,nfeat,nopti,nsubc,nsopt,td)

c       FEATURE - deposit values in command table

        if (nfeat.gt.0) then
          cp0(1,nfeat) = nfeat
          cp0(2,nfeat) = nopti

          if ((nfeat.eq.1)) then                ! SWITch data

            if(prt) write (iow,2321) cc,c2
            if (pcomp(c2,'on',2)) then
              do k = 1,3
                td(k) = max(0,min(1,nint(td(k))))
              end do ! k
              if(prt) write (iow,23211) (nint(td(k)),k=1,3)
            elseif(pcomp(c2,'timf',4) .or. pcomp(c2,'time',4)) then
              td(4) = nint(td(1))
              td(1) = 1.0d0
              if(prt) write (iow,23212) nint(td(4))
            endif

          elseif (nfeat.eq.2) then              ! SOLM data

            if(prt) then
              if(nopti.le.4 .or. nopti.eq.7) then
                write (iow,2322) cc,c2,(td(k),k=1,3)
              else
                write (iow,2322) cc,c2
              endif
            endif
            if(nopti.eq.2) lagrm   = .true.
            if(nopti.eq.7) lagrm   = .true.
            if(nopti.eq.5) shakefl = .true.
            if(nopti.eq.6) rattlfl = .true.

          elseif (nfeat.eq.3) then              ! DETA data

            if(prt) write (iow,2323) cc,c2

          elseif (nfeat.eq.4) then              ! MATE data

            if(nint(td(1)).eq.0) then
              td(1) = 1.d0
            endif
            if(nint(td(2)).eq.0) then
              td(2) = td(1)
            endif
            if(prt) write (iow,2324) nint(td(1)),nint(td(2))

          elseif (nfeat.eq.5) then              ! AUGM data

            if(pcomp(c2,'    ',4)) c2 = 'Basi'
            if(prt) write (iow,2325) cc,c2

          elseif (nfeat.eq.6) then              ! TOLE data

            cl   = cc
            tlfl = .true.
            if(nopti.eq.2) then
              if(cp0(4,nfeat).eq.0.0d0) then
                td(2) = 1.d-08
              else
                td(2) = cp0(4,nfeat)
              endif
              if(cp0(5,nfeat).eq.0.0d0) then
                td(3) = 1.d-05
              else
                td(3) = cp0(5,nfeat)
              endif
            elseif(nopti.eq.3) then
              td(2) = td(1)
              if(cp0(3,nfeat).eq.0.0d0) then
                td(1) = 1.d-08
              else
                td(1) = cp0(4,nfeat)
              endif
              if(cp0(5,nfeat).eq.0.0d0) then
                td(3) = 1.d-05
              else
                td(3) = cp0(5,nfeat)
              endif
            elseif(nopti.eq.4) then
              td(3) = td(1)
              if(cp0(3,nfeat).eq.0.0d0) then
                td(1) = 1.d-08
              else
                td(1) = cp0(3,nfeat)
              endif
              if(cp0(4,nfeat).eq.0.0d0) then
                td(2) = 1.d-08
              else
                td(2) = cp0(4,nfeat)
              endif
            endif
            do k = 1,3
              tl(k) = td(k)
            end do ! k

          elseif (nfeat.eq.7) then              ! ADHE data

            if(nopti.eq.1) then
              if(prt) write (iow,2327) cc,c2
            else
              if(prt) write (iow,2327) cc,c2,td(1)
            endif

          elseif (nfeat.eq.8) then              ! PENE data
            if(nopti.eq.1) then
              if(prt) write (iow,2328) cc,'ON'
              td(1) = 1.d0
            else
              if(prt) write (iow,2328) cc,'OFF'
              td(1) = 0.d0
            endif

          elseif (nfeat.eq.9) then              ! AXISymetric data
            if(ndm.le.2) then
              if(prt) write (iow,2329) 'Axisymmetric'
              td(1) = 1.d0
            else
              write (iow,3002)
              write (ilg,3002)
            endif

          elseif (nfeat.eq.10) then              ! INTErpolation method
            if(prt) write (iow,23210) cc,c2

          elseif (nfeat.eq.11) then              ! QUADrature order
            if(prt) write (iow,23213) cc,(k,td(k),k=1,ndm-1)

          endif

c         Store parameters in feature table

          do k = 1,14
            cp0(k+2,nfeat) = td(k)
          end do ! k

c         SUB-COMMAND - perform subcommand
c         No subcommand defined till now

        elseif (nsubc.gt.0) then
          if (prt) then
            write (iow,2330) cc
          endif

c         TYPE DATA --> read type data
c         Actually no pair data exist

        elseif (nfeat.eq.-3) then
          if (prt) then
            write (iow,2340)
          endif
c         call cr.... ()
          write (ilg,3001) npair,xxx
          write (iow,3001) npair,xxx
          write (  *,3001) npair,xxx
          call coutcis(3)
          call plstop()

c         Blank line found - search again

        elseif (nfeat.eq.-1) then
          continue

c         New command found - go back to PCONT

        elseif (nfeat.eq.-2) then
          backspace (ior)
          whfl = .false.
        endif
      end do ! while

c     Output tolerance values if input

      if(tlfl .and. prt) write (iow,2326) cl,(tl(k),k=1,3)

2300  format (/5x,'C o n t a c t   P a i r   D a t a'//
     &         5x,'Data Set for Pair Number          ',i4:/
     &         5x,'Contact Surface Dimension         ',i4/)

2310  format (/5x,'Pair Label: ',a)

2320  format ( 5x,'Contact Pair Type:                ',a)

2321  format (/5x,'  P a i r   F e a t u r e:        ',a/
     &         5x,'  Switch on/off                   ',a/
     &         5x,'    OFF  - Switch Off All Pairs   '/
     &         5x,'    ON   - Switch On Pair List    '/
     &         5x,'    TIMF - Use Life Time Function ')

23211 format (/5x,'    Normal (on = 1, off = 0)      ',i5/
     &         5x,'    Tangential (Friction)         ',i5/
     &         5x,'    Thermal                       ',i5)

23212 format ( 5x,'      Proportional load number    ',i5)

2322  format (/5x,'  P a i r   F e a t u r e:        ',a/
     &         5x,'  Solution Method:                ',a/
     &         5x,'    PENA - Penalty                '/
     &         5x,'    LAGM - Lagrangian Multiplier  '/
     &         5x,'    PERT - Perterbed Tangent      ':/
     &         5x,'  Normal     Penalty             ',1p,e12.5:/
     &         5x,'  Tangential Penalty             ',1p,e12.5:/
     &         5x,'  Normal Force Limit             ',1p,e12.5)

2323  format (/5x,'  P a i r   F e a t u r e:        ',a/
     &         5x,'  Contact Detection Method:       ',a/
     &         5x,'    BASI - Search Each Iteration  '/
     &         5x,'    SEMI - Search to Iteration 1  '/
     &         5x,'    RIGI - Not Available          ')

2324  format (/5x,'  P a i r   F e a t u r e:        ',a/
     &         5x,'  Contact Material Model:         '/
     &         5x,'    Material for Surface 1        ',i4/
     &         5x,'    Material for Surface 2        ',i4)

2325  format (/5x,'  P a i r   F e a t u r e:        ',a/
     &         5x,'  Augmentation Method:            ',a)

2326  format (/5x,'  P a i r   F e a t u r e:        ',a/
     &         5x,'  Contact Tolerances:'/
     &         5x,'    Initial Penetration Tolerance ',1p,e12.5/
     &         5x,'    Open Gap Tolerance            ',1p,e12.5/
     &         5x,'    Out of Facet Tolerance        ',1p,e12.5)

2327  format (/5x,'  P a i r   F e a t u r e:        ',a/
     &         5x,'  Contact Adhesion:               ',a:/
     &         5x,'    Adhesion Stress               ',1p,e12.5)

2328  format (/5x,'  P a i r   F e a t u r e:        ',a/
     &         5x,'  Contact Initial Penetrate Check ',a)

2329  format (/5x,'  P a i r   F e a t u r e:        '/
     &         5x,'  ',a12,' contact            ','ON')

23210 format (/5x,'  P a i r   F e a t u r e:        ',a/
     &         5x,'  Interpolation Method:           ',a)

23213 format (/5x,'  P a i r   F e a t u r e:        ',a/
     &         5x,'  Quadrature Order:               '/
     &       (10x,'Direction',i2,' = '1p,1e12.5:))

2330  format (/5x,'  Pair Sub-command                ',a)

2340  format ( 5x,'  Type Declaration Data:          ')

2411  format (/5x,'Node-to-Segment Contact Model'/
     &         5x,'  Slave  Surface                  ',i4/
     &         5x,'  Master Surface                  ',i4)

2412  format (/5x,'Point-to-Point Contact Model'/
     &         5x,'  Slave  Surface                  ',i4/
     &         5x,'  Master Surface                  ',i4/
     &         5x,'  Direction for contact           ',i4)

2414  format (/5x,'Point-to-Rigid Contact Model'/
     &         5x,'  Slave  Surface                  ',i4/
     &         5x,'  Master Surface                  ',i4)

2415  format (/5x,'Tied Surface-to-Surface Model'/
     &         5x,'  Slave  Surface                  ',i4/
     &         5x,'  Master Surface                  ',i4)

3001  format (/' *ERROR* CRPAIR: Reading data for PAIR ',i5/
     &         '  Unrecognized data for command: '/
     &         3x,a)

3002  format (/' *ERRORE CRPAIR: Axisymmetric not allowed in 3-d')

      end
