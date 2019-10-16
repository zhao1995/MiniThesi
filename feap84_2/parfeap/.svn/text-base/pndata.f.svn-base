c$Id:$
      subroutine pndata(tx,ct,nxd,nxn,nne,labl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change np33 to npix in call to pltcon            08/04/2011
c       2. Change npix to plix for plots                    16/11/2011
c       3. Update call to rprint                            09/01/2013
c          Add nne to argument list
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Collected nodal data for plots
c               Parallel version

c      Inputs:
c         tx(2)     - Text identifier data
c         ct(3)     - Plot command parameters

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comblk.h'
      include  'cdata.h'
      include  'cdat1.h'
      include  'chdata.h'
      include  'comfil.h'
      include  'fdata.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'plcapt.h'
      include  'pointer.h'
      include  'prange.h'
      include  'tdata.h'

      character tx(2)*15, gplf*128,fext*5, c*1
      logical   labl
      logical   pcomp,palloc, setvar
      logical   ndatafl
      integer   i,icp,k4
      integer   ii,jj,kk
      integer   n,nxd,nxn,nne
      integer   inode, indx
      real*8    ct(3)
      real*8    value

      save

      if(pcomp('disp',tx(2),4) .or. pcomp('stre',tx(2),4) .or.
     &   pcomp('pstr',tx(2),4) ) then

        if(pcomp('disp',tx(2),4)) c = 'd'
        if(pcomp('stre',tx(2),4)) c = 's'
        if(pcomp('pstr',tx(2),4)) c = 'p'

        inquire(unit=ios,opened = ndatafl)
        if(ndatafl) then
          close(ios)
        endif

        setvar = palloc(113,'TEMP3',numnp,2)

        kk = 1
        jj = 1
        do while (kk .le. jj)
          call adomnam(fplt,gplf,kk)
          gplf(1:1) = 'G'
          fext(1:5) = ' '
          ii = int(ct(1))

          if(ii.lt.10) then
            write(fext,'(a1,a3,i1)') c,'000',ii
          elseif(ii.lt.100) then
            write(fext,'(a1,a2,i2)') c,'00',ii
          elseif(ii.lt.1000) then
            write(fext,'(a1,a1,i3)') c,'0',ii
          elseif(ii.lt.10000) then
            write(fext,'(a1,i4)') c,ii
          endif
          call addext(gplf,fext,128,5)
          inquire(file=gplf, exist=ndatafl )
          if(ndatafl) then
            open(unit=ios,file=gplf,form='formatted',status='old')
            rewind ios
            caption(1:15) = ' '
            read(ios,*) jj
            read(ios,*) ttim,caption, ii
            write(caption,'(a13,i2)'),caption(1:13),ii
            ncapt = 1
c           Read data
            do  indx = 1, numnp
              read(ios,*,end=55) inode,value
              hr(np(113) + inode - 1) = value
            end do
55          close(ios)
          else
            write(*,*) 'File ',gplf,' not found'
          endif
          kk = kk + 1
        end do

c       Plot nodal data
        icp = 1
        i   = 1
        call rprint(mr(plix),nxn,nxd,nne,hr(np(113)),icp,k4)

        call pltcon(hr(np(53)),mr(np(32)),mr(plix),mr(np(62)),
     &              hr(np(113)),nie,3,icp,nxd,nxn,i,n,i,2,labl)

c       Destroy nodal data array
        setvar = palloc(113,'TEMP3',0,2)
      else
        write(*,*) 'Allowed NDATa options disp, stre, pstr'
      endif ! acceptable options

      end
