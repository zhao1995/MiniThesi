c$Id:$
      subroutine ploadi(prt,prth,error,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Input data for surface loading from elements

c      Inputs:
c         prt       - Output generated data if true
c         prth      - Output title/header if true
c         isw       - (1) determine data amount

c      Outputs:
c         error     - Flag, true if error detected during inputs
c                     Data stored in pointers for later use
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'allotd.h'
      include  'cdata.h'
      include  'cdat2.h'
      include  'chdata.h'
      include  'comfil.h'
      include  'conval.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'sldata.h'
      include  'trdata.h'
      include  'comblk.h'

      logical   prt,prth,error, pinput,setvar,palloc
      integer   i,k,n, ns,nv,nl, nsv,nrec
      integer   isw, iorsv
      integer   ixl(16)
      real*8    valu(64),td(52)

      save

c     Input of surface loading data to compute required data

      if(isw.eq.1) then
        if(numsl.eq.0) then
          sldfil = finp
          sldfil(1:1) = 'J'
          open (unit = ios, file = sldfil)
        endif

        write(ios,1002) prt,prth
        write(ios,1001) vvv
        write(ios,1001) tr,xr,trdet,x0

100     if(ior.lt.0) write(*,3000)
        setvar = pinput(td,4)
        if(setvar) go to 100
        ns = nint(td(2))
        nv = nint(td(3))
        nl = nint(td(4))
        if(nl.eq.0) nl = 1

c       Set parameters

        numsl         = numsl + 1
        iels(1,numsl) = nint(td(1))
        iels(2,numsl) = nim
        iels(3,numsl) = nre
        iels(4,numsl) = nl
        inods(numsl)  = ns
        ivals(numsl)  = nv

c       Input records for current loading state

        nrec  = 0
        error = .false.
1       nsv = ns+nv
110     do i = 1,nsv,8
          k = min(nsv-i+1,8)
          setvar = pinput(td(i),k)
          if(setvar) go to 110
          write(ios,'(a)') xxx
        end do ! i
        if(nint(td(1)).ne.0) then
          nim  = nim + ns
          nre  = nre + nv
          nrec = nrec + 1
          go to 1
        endif
        iels(5,numsl) = nrec

c     Allocate memory

      elseif(isw.eq.2 .and. numsl.gt.0) then

        close(ios)
        setvar = palloc( 37, 'SLODI', nim, 1 )
        setvar = palloc(241, 'SLODR', nre, 2 )

c       Input of surface loading data

        open (unit = ios, file = sldfil)

        error = .false.
        iorsv = ior
        ior   = ios
        rewind ior

        do n = 1,numsl

          read(ios,1002) prt,prth
          read(ios,1001) vvv
          read(ios,1001) tr,xr,trdet,x0

          nim = iels(2,n)
          nre = iels(3,n)
          nl  = iels(4,n)
          ns  = inods(n)
          nv  = ivals(n)

c         Input records for current loading state

          if(prt) then
            call prtitl(prth)
            write(iow,2000) (i,i=1,ns)
            write(iow,2001) (i,i=1,nv)
            if(ior.lt.0) then
              write(*,2000) (i,i=1,ns)
              write(*,2001) (i,i=1,nv)
            endif
          endif
2         nsv   = ns + nv
210       do i = 1,nsv,8
            k = min(nsv-i+1,8)
            setvar = pinput(td(i),k)
            if(setvar) go to 210
          end do ! i
          if(nint(td(1)).ne.0) then
            do i = 1,ns
              ixl(i) = nint(td(i))
            end do ! i
            do i = 1,nv
              valu(i) = td(i+ns)
            end do ! i
            do i = 1,ns
              if(ixl(i).gt.numnp) then
                write(ilg,3001) (ixl(k),k=1,ns)
                write(iow,3001) (ixl(k),k=1,ns)
                if(ior.lt.0) write(*,3001) (ixl(k),k=1,ns)
                error = .true.
                go to 2
              endif
            end do ! i
            if(prt) then
              if(ior.gt.0) then
                write(iow,2002) (ixl(i),i=1,ns)
                write(iow,2003) (valu(i),i=1,nv)
              else
                write(*,2002) (ixl(i),i=1,ns)
                write(*,2003) (valu(i),i=1,nv)
              endif
            endif
            call pstore(ns,nv,ixl,valu,mr(np(37)+nim),hr(np(241)+nre))
            nim = nim + ns
            nre = nre + nv
            go to 2
          endif
        end do ! n

        close(ios,status = 'delete')
        ior           = iorsv
      endif

c     Formats

1001  format (1p,4e20.12)
1002  format (2l5)

2000  format('     S u r f a c e   L o a d i n g'//4x,10(i3,'-node'))

2001  format(5x,6(i5,'-value '))

2002  format(/8i8)

2003  format(4x,1p,6e12.4)

3000  format(' Input: el-type, #-nodes, #-vals, ld-type'/'   >',$)

3001  format(' *ERROR* PLOADI: Node too large'/(8i8))

      end
