c$Id:$
      subroutine contlib (csw,npair,cs0,cm0,cp0,ics,cm,ch1,ch2,ch3,hic)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change fp(1) to point                            15/01/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact DRIVers LIBrary

c      Purpose: Switch to selected driver to perform action

c      Inputs :
c         csw     - Contact switch
c         npair   - # of current pair
c         cs0(*)  - Contact surfaces control data
c         cm0(*)  - Contact material control data
c         cp0(*)  - Contactpair control data
c         ics(*)  - Contact element nodal connection array
c         cm(*)   - Contact materials data storage
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)
c         ch3(*)  - Contact history variables (static)
c         hic(*)  - HIstory Correspondence table

c      Outputs:
c                 - Perform requested activity
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_dict.h'
      include  'c_geom.h'
      include  'c_mate.h'
      include  'c_pair.h'
      include  'cdata.h'
      include  'ddata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'comblk.h'

      include  'p_point.h'  ! Get temporary pointer

      logical   ifscan
      integer   csw,npair,ics(*), ndv
      integer   hic((c_lp1+c_lp3),*)
      real*8    cs0(nr0,n0c1:nc01,*),cm0(nr0,n0c2:nc02,*)
      real*8    cp0(nr0,n0c3:nc03,*),cm(*),ch1(*),ch2(*),ch3(*)

      save

      call cdebug0 ('  contlib',csw)

c     Manual sel. performed to show contact characteristics

      if ((csw.eq.200) .or. (csw.eq.400)) then
        ndv    =  npair
        ifscan = .true.

c     Automatic sel. performed in normal operations

      else
        call setcomp (npair,cs0,cm0,cp0,hic,csw)
        ndv    =  ndrv
        ifscan = .false.
      endif

c     Contactpair not active -> scan driver YES, normal operations NO

      if ((ifon.eq.1) .or. (ifscan)) then

c       Set pointer for displacement array

        if(nrk.eq.0) then
          point = np(40)
        elseif(nrk.gt.0) then
          point = np(42) + ndf*numnp*(nrk-1)
        else
          write(  *,3000) nrk
          write(ilg,3000) nrk
          call plstop()
        endif

c       User Contact Element Type 1

        if (ndv.eq.1) then
          call celmt01 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 2

        elseif (ndv.eq.2) then
          call celmt02 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 3

        elseif (ndv.eq.3) then
          call celmt03 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 4

        elseif (ndv.eq.4) then
          call celmt04 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 5

        elseif (ndv.eq.5) then
          call celmt05 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 6

        elseif (ndv.eq.6) then
          call celmt06 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 7

        elseif (ndv.eq.7) then
          call celmt07 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 8

        elseif (ndv.eq.8) then
          call celmt08 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 9

        elseif (ndv.eq.9) then
          call celmt09 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 10

        elseif (ndv.eq.10) then
          call celmt10 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 11

        elseif (ndv.eq.11) then
          call celmt11 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 12

        elseif (ndv.eq.12) then
          call celmt12 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 13

        elseif (ndv.eq.13) then
          call celmt13 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 14

        elseif (ndv.eq.14) then
          call celmt14 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 15

        elseif (ndv.eq.15) then
          call celmt15 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 16

        elseif (ndv.eq.16) then
          call celmt16 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 17

        elseif (ndv.eq.17) then
          call celmt17 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 18

        elseif (ndv.eq.18) then
          call celmt18 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 19

        elseif (ndv.eq.19) then
          call celmt19 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       User Contact Element Type 20

        elseif (ndv.eq.20) then
          call celmt20 (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf1),cs0(1,n0c1,nsurf2),
     &                  cm0(1,n0c2,nmat1), cm0(1,n0c2,nmat2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),cm(ofm1),cm(ofm2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       Program Contact Element Type NToS

        elseif (ndv.eq.21) then
          if(ndm.eq.2 .or. cndm.eq.2) then
            call cnts2d (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                   cs0(1,n0c1,nsurf2),
     &                   cp0(1,n0c3,npair),
     &                   ics(ofs1),ics(ofs2),cm(ofm1),
     &                   ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                   w1(1,ndv),w3(1,ndv))
          elseif(ndm.eq.3) then
            call cnts3d (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                   cp0(1,n0c3,npair),
     &                   ics(ofs1),ics(ofs2),cm(ofm1),
     &                   ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                   w1(1,ndv),w3(1,ndv))
          endif

c       Program Contact Element Type PToP or NToN Problems

        elseif (ndv.eq.22 .or. ndv.eq.23) then
          call cptpnd (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                 cp0(1,n0c3,npair),ics(ofs1),ics(ofs2),
     &                 cm(ofm1),ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                 w1(1,ndv),w3(1,ndv))

c       Program Contact Element Type PToR or NToR Problems

        elseif (ndv.eq.24 .or. ndv.eq.25) then
          call cntrnd  (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))

c       Program Contact Element Type TIED

        elseif (ndv.eq.26) then
          call ctied2d (ndm,ndf,hr(np(43)),hr(point),csw,npair,
     &                  cs0(1,n0c1,nsurf2),
     &                  cp0(1,n0c3,npair),
     &                  ics(ofs1),ics(ofs2),
     &                  ch1(ofh1),ch2(ofh1),ch3(ofh3),
     &                  w1(1,ndv),w3(1,ndv))
        endif

      endif

c     Format

3000  format(' *ERROR* Incorrect solution vector "u" location:',
     &       ' nrk =',i10)

      end
