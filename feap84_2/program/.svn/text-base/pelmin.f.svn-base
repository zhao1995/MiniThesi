c$Id:$
      subroutine pelmin(tx,idl,ix,rben,nen1,prt,prth,error,elabel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'elabel' to describe current generation type 20/07/2007
c       2. Add set of last element number to last_elm       29/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Data input routine for element connections

c      Inputs:
c         tx        - Option identifier
c         nen1      - Dimension for ix array
c         prt       - Flag, print input data if true
c         prth      - Flag, print title/header if true
c         rben      - Rigid generation flag

c      Scratch:
c         idl(*)    - Local degree of freedom integer data

c      Outputs:
c         ix(*)     - Element nodal connection lists
c         rben(*)   - Rigid body numbers for elements
c         error     - True if error occurs during input
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cblktr.h'
      include  'cdata.h'
      include  'chdata.h'
      include  'dstars.h'
      include  'iofile.h'
      include  'p_ptname.h'
      include  'region.h'
      include  'rigid2.h'

      logical   prt,prth,error,errck,fprt,genfl
      logical   pcomp,pinput,vinput,cksep
      character elabel*(*),tx*15,mtype*200
      integer   i,j,k,l,n,nen1, ii,il,is,ilx,lg,lk,llx,ma,ng
      integer   ixg(16),idl(*),ix(nen1,*),rben(*)
      real*8    td(16)

      save

c     Check for input format descriptor

      if(nen.le.16) then
        do j = 3,80
          if(cksep(xxx(j:j))) then
            do i = j+1,80
              if(cksep(xxx(i:i))) go to 200
            end do ! i
          endif
        end do ! j
        i = 199
200     mtype = xxx(i+1:200)
        fprt = vinput(mtype,200,td,nen)
        i = 0
        do n = 1,nen
          ixg(n) = nint(td(n))
          i      = max(i,abs(ixg(n)))
        end do ! n
        if(i.gt.0) then
          if(prt) call iprint(ixg,1,nen,1,'Element Generation Array')
          genfl = .true.
        else
          genfl = .false.
        endif
      else
        genfl = .false.
      endif

c     Perform input of data for element connections

      ilx = 0
      l   = 0
      ma  = 0
      if(ior.lt.0) write(*,2003)
      do i = 1, numel, 50
        j = min(numel,i+49)
        do n = i, j
          fprt = .false.
          if (l .lt. n) then
            llx = ilx

c           Input element records - N.B. limit is 16 nos. / record

            il = min(nen+3,16)
201         errck = pinput(td,il)
            if(errck) go to 201

c           Discontinue input on blank record or negative number.

            l  = nint(td(1))
            if(l .le. 0) then
              last_elm = n
              return
            else
              l   = l + starel
              neo = max(neo,l)
            endif

c           Set element group number

            is = nint(td(2))

c           Transfer to element integer data

            do k = 3,il
              idl(k-2) = nint(td(k))
            end do ! k

c           Input additional records

            do ii = 14,nen,16
              il = min(nen+1-ii,16)
202           errck = pinput(td,il)
              if(errck) go to 202

c             Transfer to element integer data

              do k = 1,il
                idl(ii+k) = nint(td(k))
              end do ! k
            end do ! ii

c           Check if old sequence

            if(pcomp(tx,'old',3)) then

              lg = idl(nen+1)
              lk = is

c           Else change order

            else

              lg = is
              lk = idl(1)
              do k = 1,nen
                idl(k) = idl(k+1)
              end do ! k

            endif

c           Add star node number to inputs

            do k = 1,nen
              if(idl(k).gt.0) idl(k) = idl(k) + starnd
            end do ! k

            if (lg .eq. 0) then
              lg = 1
            endif
            ilx = lg
          endif

c         Error in input data

          if (l .lt. n) then
            write(ilg,3001) l,n
            write(iow,3001) l,n
            if(ior.lt.0) then
              write(*,3001) l,n
            endif
            error = .true.

c         Generate missing elements

          else if ((l .gt. n) .and. (llx .ne. 0)) then
            do k = 1, nen
              if(genfl) then
                ix(k,n) = ix(k,n-1) + ixg(k)
              else
                ix(k,n) = ix(k,n-1) + ng
              endif
              if (ix(k,n-1) .eq. 0) then
                ix(k,n) = 0
              endif
              if ((ix(k,n) .gt. numnp) .or. (ix(k,n) .lt. 0)) then
                write(ilg,3002) n
                write(iow,3002) n
                if(ior.lt.0) then
                  write(*,3002) n
                endif
                error = .true.
              endif
            end do ! k
            ix(nen1,n)   = ix(nen1,n-1)
            ix(nen1-1,n) = nreg
            rben(n)      = nrigid
            fprt         = .true.

c         Transfer input to current element

          else if (l.eq.n) then
            ng = lg
            do k = 1, nen
              if ((idl(k) .gt. numnp) .or. (idl(k) .lt. 0)) then
                write(ilg,3002) n
                write(iow,3002) n
                if(ior.lt.0) then
                  write(*,3002) n
                endif
                error = .true.
              endif
              ix(k,l) = idl(k)
            end do ! k
            ix(nen1,l)   = lk
            ix(nen1-1,l) = nreg
            rben(l)      = nrigid
            fprt         = .true.
          endif

c         Output element list

          if ((prt) .and. (.not. error) .and. fprt) then
            if(mod(ma,50).eq.0) then
              call prtitl(prth)
              write(iow,2001) elabel,(k,k=1,nen)
              if(ior.lt.0) then
                write(*,2001) elabel,(k,k=1,nen)
              endif
            endif
            ma = ma + 1
            write(iow,2002) n,ix(nen1,n),ix(nen1-1,n),(ix(k,n),k=1,nen)
            if(ior.lt.0) then
              write(*,2002) n,ix(nen1,n),ix(nen1-1,n),(ix(k,n),k=1,nen)
            endif
          endif
        end do ! n
      end do ! i

c     Formats

2001  format(5x,a//3x,'Elmt Mat Reg',
     &           8(i3,' Node'):/(15x,8(i3,' Node')))

2002  format(i7,2i4,8i8:/(15x,8i8))

2003  format(' Input: Elmt#, Matl#, (ix(i),i=1,nen), inc'/3x,'>')

3001  format(' *ERROR* PELMIN: Element',i5,' appears after element',i5)

3002  format(' *ERROR* PELMIN: Element',i5,' has illegal nodes')

      end
