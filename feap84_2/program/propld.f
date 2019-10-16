c$Id:$
      function propld(t,lctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set values for Type 4. Limit ixp(i) to 1 or 2    08/06/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Proportional load table (maximum 50)

c      Load types:
c         1)
c            tmin = a(1)
c            tmax = a(2)
c            prop = a(3) + a(4)*t + a(5)*sin(a(6)*(t-tmin)+a(7))**iex
c                   for: tmin <= t <= tmax

c         2.) Piecewise linear interpolation
c            prop = vm(i)+(vm(i+1)-vm(i))*(t-tm(i))/(tm(i+1)-tm(i)
c                   tm(i) < t < tm(i+1)

c         3.) Roll and slide
c            tmin = a(1)
c            tmax = a(2)

c         4.) Sawtooth - cyclic loading

c         5.) Polynomial
c            prop - a(3) + a(4)*t + a(5)*t^2 + a(6)*t^3 + a(7)*t^4
c            tmin = a(1)
c            tmax = a(2)

c         6.) Exponential -
c            prop = a(3) + exp(a(5)*(t-tmin))*(a(4)
c                 + a(6)*(sin(a(7)*(t-tmin)))**iex)
c            tmin = a(1)
c            tmax = a(2)

c      Individual proportional factors also stored in /prld1/
c      Total is returned as value of function

c      Inputs:
c         t         - Time for lookup
c         lctl(2)   - Number of proportional load to define.
c                     N.B. Only used for inputs
c      Outputs:
c         propld    - Value of total proportional load
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'prlod.h'
      include  'prld1.h'
      include  'pointer.h'
      include  'tdata.h'
      include  'umac1.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   pcomp, errck, pinput, setvar, palloc
      integer   i,j,k,l,m,lunit,nprtot, lctl(2), ilast(50),itime(50)
      real*8    propld,propl,prop2i,prop2p,prop2,prop3,upropld
      real*8    t, tt, td(16)

      save

c     Initialize propld

      propld = 0.0d0

c     Delete proportional load

      if(pcomp(uct,'dele',4)) then
        do i=lctl(1),lctl(2)

          if(ik(i).eq.2 ) then
            j = mr(np(121) + 2*i - 1)
            l = mr(np(121) + 2*i) - j
            k = 0
            do m = 1,npld
              if(ik(m).eq.2 .and. mr(np(121)+2*m-1).gt.j) then
                mr(np(121)+2*m-1) = mr(np(121)+2*m-1) - l
                mr(np(121)+2*m  ) = mr(np(121)+2*m  ) - l
                k                 = max(k,mr(np(121)+2*m))
              endif
            end do ! m
            mr(np(121))       = k               ! Reset last entry value
            mr(np(121)+2*i-1) = 0               ! Zero deleted entry
            mr(np(121)+2*i)   = 0
            do m = j,k-1
              hr(np(122)+m)   = hr(np(122)+m+l) ! Pack values
              hr(np(122)+m+1) = 0.0d0           ! Zero moved value
            end do ! m
          else
            do m = 1,7
              ap(m,i) = 0.0d0
            end do
          end if
          prldv(i) = 0.0d0
          ik(i)    = 0
          j        = npld
          do m = 1,j
            if(ik(m).gt.0) npld = m
          end do
          write(iow,3002) i,npld
          if(ior.lt.0) then
            write(*,3002) i,npld
          endif
        end do
        return
      end if

c     Input table of proportional loads

      if(lctl(1).gt.0) then
        do i=lctl(1),lctl(2)
          write(iow,2000) i
          if(ior.lt.0) then
            write(*,2000) i
            write(*,2002)
            call pprint('   >')
          endif
          if(pcomp(uct,'load',4)) then
            ap(1,i) = 0.d0
            ap(2,i) = 1.d+20
            ap(3,i) = prldv(i) - t
            ap(4,i) = 1.d0
            ap(5,i) = 0.d0
            ap(6,i) = 0.d0
            ap(7,i) = 0.d0
            ik(i)   = 1
          elseif(pcomp(uct,'relo',4)) then
            ap(1,i) = 0.d0
            ap(2,i) = 1.d+20
            ap(3,i) =-t
            ap(4,i) = 1.d0
            ap(5,i) = 0.d0
            ap(6,i) = 0.d0
            ap(7,i) = 0.d0
            ik(i)   = 1
          elseif(pcomp(uct,'unlo',4)) then
            ap(1,i) = 0.d0
            ap(2,i) = 1.d+20
            ap(3,i) = prldv(i) + t
            ap(4,i) =-1.d0
            ap(5,i) = 0.d0
            ap(6,i) = 0.d0
            ap(7,i) = 0.d0
            ik(i)   = 1
          elseif(pcomp(uct,'ramp',4)) then
            ap(1,i) = 0.d0
            ap(2,i) = 1.d+20
            ap(3,i) = 1.d0 + t
            ap(4,i) =-1.d0
            ap(5,i) = 0.d0
            ap(6,i) = 0.d0
            ap(7,i) = 0.d0
            ik(i)   =-1
          elseif(pcomp(uct,'cons',4)) then
            ap(1,i) = 0.d0
            ap(2,i) = 1.d+20
            ap(3,i) = prldv(i)
            ap(4,i) = 0.d0
            ap(5,i) = 0.d0
            ap(6,i) = 0.d0
            ap(7,i) = 0.d0
            ik(i)   = 1
          elseif(pcomp(uct,'user',4)) then
201         errck = pinput(td,9)
            if(errck) go to 201
            j        =  abs(nint(td(1)))
            ik(i)    = -j - 100
            iexp(i)  = nint(td(2))
            ap(1,i)  = td(3)
            ap(2,i)  = td(4)
            ap(3,i)  = td(5)
            ap(4,i)  = td(6)
            ap(5,i)  = td(7)
            ap(6,i)  = td(8)
            ap(7,i)  = td(9)
            call uprop(j,iexp(i),ap(1,i),ap(2,i),ap(3,i), t,upropld, 1)
          else
101         errck = pinput(td,9)
            if(errck) go to 101
            ik(i) = nint(td(1))
            if(ik(i).gt.0) then

c             Type 1.

              if(ik(i).eq.1) then
                iexp(i)  = nint(td(2))
                ap(1,i)  = td(3)
                ap(2,i)  = td(4)
                ap(3,i)  = td(5)
                ap(4,i)  = td(6)
                ap(5,i)  = td(7)
                ap(6,i)  = td(8)
                ap(7,i)  = td(9)

c             Type 2.

              elseif(ik(i) .eq. 2) then
                setvar = palloc(121,'PROP0',2*npld+1,1) ! allocate table
c               Keep table clean
                do m = 1,npld
                  if(ik(m).ne.2) then
                    mr(np(121)+2*m-1) = 0
                    mr(np(121)+2*m  ) = 0
                  endif
                end do ! m
                if(mr(np(121)+2*i).gt.0) then
                  if(ior.lt.0) then
                    write(*,3001) i
                    return
                  endif
                  write(ilg,3001) i
                  write(iow,3001) i
                  call plstop()
                else
                  l      = max(1,int(td(2)))
                  propl  = prop2i(lunit,l,ilast(i))
                  setvar = palloc(123,'PROP2',2*ilast(i),2)
                  propl  = prop2p(lunit,l,hr(np(123)),itime(i))

c                 Set offset pointer into array PROP0

                  nprtot            = mr(np(121))
                  mr(np(121)+2*i-1) = nprtot
                  mr(np(121)+2*i  ) = nprtot + 2*ilast(i)
                  mr(np(121))       = nprtot + 2*ilast(i)

c                 Move  PROP2 to PROP1 and destroy PROP2

                  setvar   = palloc(122,'PROP1',nprtot+2*ilast(i),2)
                  call pmove(hr(np(123)),hr(np(122)+nprtot),2*ilast(i))
                  setvar   = palloc(123,'PROP2',0,2)

                endif

c             Type 3.  Slide and roll

              elseif(ik(i) .eq. 3) then

                ap(1,i)  = td(2)
                ap(2,i)  = td(3)
                ap(3,i)  = td(4)
                ap(4,i)  = td(5)
                ap(5,i)  = td(6)
                ap(6,i)  = td(7)
                ap(7,i)  = td(8)

c             Type 4.

              elseif(ik(i).eq.4) then

                iexp(i)  = max(1,min(2,nint(td(2))))
                ap(1,i)  = td(3)
                ap(2,i)  = td(4)
                ap(3,i)  = td(5)
                ap(4,i)  = td(6)
                ap(5,i)  = td(7)
                ap(6,i)  = td(8)
                ap(7,i)  = td(9)

                if(ap(4,i).eq.0.0d0) then
                  write(iow,2004)
                  if(ior.lt.0) then
                    write(*,2004)
                    go to 101
                  endif
                  call plstop()
                endif

c             Type 5.  Polynomial

              elseif(ik(i) .eq. 5) then

                ap(1,i)  = td(2)
                ap(2,i)  = td(3)
                ap(3,i)  = td(4)
                ap(4,i)  = td(5)
                ap(5,i)  = td(6)
                ap(6,i)  = td(7)
                ap(7,i)  = td(8)

c             Type 6.  Exponential

              elseif(ik(i) .eq. 6) then

                iexp(i)  = max(1,nint(td(2)))
                ap(1,i)  = td(3)
                ap(2,i)  = td(4)
                ap(3,i)  = td(5)
                ap(4,i)  = td(6)
                ap(5,i)  = td(7)
                ap(6,i)  = td(8)
                ap(7,i)  = td(9)

c             Error on input type

              else
                write(ilg,3000) ik(i)
                write(iow,3000) ik(i)
                if(ior.lt.0) then
                  write(*,3000) ik(i)
                  go to 101
                endif
                call plstop()
              endif

c           Default: Ramp loading with unit slope.

            elseif(ik(i).eq.0) then
              ap(1,i) = 0.d0
              ap(2,i) = 1.d+20
              ap(3,i) = 0.d0
              ap(4,i) = 1.d0
              ap(5,i) = 0.d0
              ap(6,i) = 0.d0
              ik(i)   = 1
            endif

          endif

c         Output proportional parameters

          if(ik(i).eq.1 .or.ik(i).ge.4.and.ik(i).le.6) then
            write(iow,2001) i,ik(i),ap(1,i),ap(2,i),
     &                      (m,m=1,5),(ap(m,i),m=3,7),iexp(i)
            if(ior.lt.0) then
              write(*,2001) i,ik(i),ap(1,i),ap(2,i),
     &                      (m,m=1,5),(ap(m,i),m=3,7),iexp(i)
            endif
          elseif(ik(i).eq.3) then
            write(iow,2003) i,ik(i),(ap(m,i),m=1,7)
            if(ior.lt.0) then
              write(*,2003) i,ik(i),(ap(m,i),m=1,7)
            endif
          endif
        end do
      endif

c     Compute individual proportional factors at time t

      propl = 0.0d0
      do i = 1,npld

        if(t.ge.ap(1,i).and.t.le.ap(2,i) .and. ik(i).ne.2 ) then

c         Type 1. Functional proportional loading

          if(abs(ik(i)).eq.1) then
            tt       = t - ap(1,i)
            l        = max(iexp(i),1)
            prldv(i) = ap(3,i) + ap(4,i)*tt
     &               + ap(5,i)*((sin(ap(6,i)*tt+ap(7,i)))**l)

c           For 'unload' option prevent negative values

            if(ik(i).lt.0) prldv(i) = max(0.0d0,prldv(i))
            propl = propl + prldv(i)

c         Type 3. Slide and roll

          elseif(ik(i).eq.3) then

            tt       = t  - ap(1,i)
            prldv(i) = prop3(i,tt,ap(3,i),hr(np(43)),hr(np(27)),
     &                       mr(np(29)))

c         Type 4. Sawtooth proportional loading

          elseif(ik(i).eq.4) then
            l        = int((t+ap(4,i))/(ap(4,i)+ap(4,i)))
            m        = iexp(i)
            m        = int((t + (m-1)*ap(4,i))/(m*ap(4,i))) + 1
            prldv(i) = (1 - 2*mod(m,2))*ap(3,i)*(l+l-t/ap(4,i))
            propl    = propl + prldv(i)

c         Type 5.  Polynomial

          elseif(ik(i).eq.5) then

            tt       = t  - ap(1,i)
            prldv(i) = ap(3,i)     + tt*(ap(4,i)
     &                             + tt*(ap(5,i)
     &                             + tt*(ap(6,i)
     &                             + tt* ap(7,i))))
            propl    = propl + prldv(i)

c         Type 6.  Exponential

          elseif(ik(i).eq.6) then

            tt       = t  - ap(1,i)
            l        = max(iexp(i),1)
            prldv(i) = ap(3,i)     + exp(ap(5,i)*tt)*(ap(4,i)
     &                             + ap(6,i)*(sin(ap(7,i)*tt))**l)
            propl    = propl + prldv(i)

c         User proportional load

          elseif(ik(i) .le. -100) then

            j = abs(ik(i) + 100)
            call uprop(j,iexp(i),ap(1,i),ap(2,i),ap(3,i), t,upropld, 2)
            prldv(i) = upropld
            propl    = propl + upropld

          endif

c       Type 2. Linear interpolation loading

        elseif(ik(i).eq.2) then

          l        = 0
          fp(1)    = np(122) + mr(np(121) + 2*i-1)
          prldv(i) = prop2(t,l,hr(fp(1)),ilast(i),itime(i))
          propl    = propl + prldv(i)

c       Outside range of loadings

        else

          prldv(i) = 0.0d0

        endif

      end do

      propld = propl

c     Formats

2000  format(30x,'Proportional Load Table Number',i3)

2001  format(/,' Number    Type      Tmin',10x,'Tmax',/i3,i10,6x,g11.4,
     & 3x,g11.4,/6x,5('a(',i1,')',10x),'exp',/4x,5(g11.4,3x),i5/)

2002  format(' Input: type, exponent, tmin, tmax, a(i),i=1,4')

2003  format(/,' Number    Type      Tmin',10x,'Tmax',/i3,i10,6x,g11.4,
     & 3x,g11.4,/5x,'1-vel',9x,'2-vel',9x,'r-vel',9x,'1-ref',9x,'2-ref'
     &/4x,5(g11.4,3x)/)

2004  format( ' Zero a(2) in sawtooth function - reinput data')

3000  format(' *ERROR* PROPLD: Option,',i3,' does not exist:',
     &       ' Reinput data')

3001  format(' *ERROR* PROPLD: Table',i3,' already defined')

3002  format(' *WARNING* Deleted proportional load number',i3/
     &       '           Maximum proportional load number',i3)

      end
