c$Id:$
      subroutine sld3d3(d,ul,xl,tl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add augmenting for element using d(185);         14/03/2007
c          Format as an if-then-else structure
c       2. Add cylindrical material property option         15/06/2007
c       3. Add 'psig' & 'peps'  to output principal values  04/03/2008
c       4. Add 'zn' to outputs, change format to be same    12/06/2008
c          for all sld3d*.f files.
c       5. Revise augment: augfp -> augf                    14/04/2009
c       6. Add 'elcoor' to use xref                         25/07/2009
c       7. Move computation of psil to material models      26/07/2009
c       8. Add 'l' to modlsd call                           05/01/2012
c       9. Add average of density for multiscale            10/05/2012
c      10. Add eps to slcn3d call                           01/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 3-D linear elastic enhanced strain element

c     Output records:

c     Prints in element: sig-11, sig-22, sig-33, sig-12, sig-23, sig-31
c                        eps-11, eps-22, eps-33, eps-12, eps-23, eps-31
c                        sig-1 , sig-2 , sig-3 , eps-1 , eps-2 , eps-3

c     Prints at nodes:   1=sig-11, 2=sig-22, 3=sig-33
c                        1=sig-12, 2=sig-23, 3=sig-31
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'
      include  'bdata.h'
      include  'cdata.h'
      include  'elcoor.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'incshp3.h'
      include  'iofile.h'
      include  'oelmt.h'
      include  'part0.h'
      include  'prstrs.h'
      include  'qudshp.h'
      include  'rdata.h'
      include  'comblk.h'


      logical   noconv
      integer   ndf,ndm,nst,isw
      integer   i,i1,j,j1,l,istrt, nhv,nn, ninc,nhi, icount,ncount
      integer   is(24)
      real*8    alam, ha, ta ,  dvs, dvm, xn ,  yn,  zn, rr
      real*8    a11, a12, a13, a21, a22, a23, a31, a32, a33, am
      real*8    a41, a42, a43, a51, a52, a53, a61, a62, a63
      real*8    cfac,lfac,mass,unorm, uinorm, duinrm, utol
      real*8    d(*), xl(ndm,*), ul(ndf,nen,*), tl(*), s(nst,*),r(ndf,*)
      real*8    shpi(3,4,9),sig(10,8),eps(9,3),epsl(6,8)
      real*8    gg(9,24),hh(9,9),bb(9),hg(9,24),dd(6,6,5),epsm(6)
      real*8    dui(9),epsv(9), psig(3), peps(3), dot

      save

      data      utol   / 1.d-10 /
      data      nhi    /14/, ninc  /9/,  ncount / 25 /
      data      alam,ha / 2*0.0d0 /

c     Augmented update

      if(isw.eq.10) then

        hr(nh2) = hr(nh2) + augf*d(185)*hr(nh2+1)

c     Initialize constitutive history data

      elseif(isw.eq.14) then

c       Get quadrature information

        l = 2
        call int3d(l,lint,sg3)

c       Initialize history variables only

        nhv   = nint(d(15))
        istrt = nint(d(84))
        nn   = nhi
        alam = hr(nh2)
        call shpm3d(shpi,xl,ndm,.true.)
        do l = 1,lint

c         Compute strains and coordinates

          call strnen(d,xl,ul,tl,shp3(1,1,l),shpi(1,1,l),
     &                ndm,ndf,nel, eps,ta)

c         Compute constitution

          call modlsd(l,d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                dd,sig(1,l),alam,ha,isw)
          nn = nn + nhv
        end do ! l

c     ARRAY COMPUTATIONS

      elseif(isw.eq. 3 .or. isw.eq. 4 .or. isw.eq. 6 .or.
     &       isw.eq. 8 .or. isw.eq.16) then

c       Set flag to compute enhanced strains and material model #

c       Get quadrature information

        l = 2
        call int3d(l,lint,sg3)

c     Extract enhanced mode displacements

        nhv   = nint(d(15))
        istrt = nint(d(84))

        do i = 1,ninc
          ui(i,1,1) = hr(nh2+i+1)
          ui(i,1,2) = hr(nh2+i+1) - hr(nh1+i+1)
        end do ! i
        unorm = 0.0d0
        do i = 1,8
          unorm = max(unorm,abs(ul(1,i,1)),
     &                      abs(ul(2,i,1)),
     &                      abs(ul(3,i,1)))
        end do ! i

c       Compute shape functions

        call shpm3d(shpi,xl,ndm,.true.)

        noconv = .true.
        icount = 0
        do while(noconv)

c         Initialize for enhanced modes

          do i = 1,9
            bb(i) = 0.0d0
            do j = 1,9
              hh(j,i) = 0.0d0
            end do ! j
          end do ! i

          nn   = nhi
          alam = hr(nh2)
          do l = 1,lint

            dvs = jac(l)*ctan(1)

c           Compute stress and strain at point

            call strnen(d,xl,ul,tl,shp3(1,1,l),shpi(1,1,l),
     &                  ndm,ndf,nel, eps,ta)

c           Compute constitution

            call modlsd(l,d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                  dd,sig(1,l),alam,ha,isw)

c           Residual computations

            do j = 1,3
              xn = shpi(1,j,l)*jac(l)
              yn = shpi(2,j,l)*jac(l)
              zn = shpi(3,j,l)*jac(l)

c             Compute enhanced stress divergence loads

              bb(3*j-2) = bb(3*j-2) - sig(1,l)*xn
     &                              - sig(4,l)*yn
     &                              - sig(6,l)*zn
              bb(3*j-1) = bb(3*j-1) - sig(2,l)*yn
     &                              - sig(4,l)*xn
     &                              - sig(5,l)*zn
              bb(3*j  ) = bb(3*j  ) - sig(3,l)*zn
     &                              - sig(5,l)*yn
     &                              - sig(6,l)*xn
            end do ! j

            do j = 1,3
              xn = shpi(1,j,l)*dvs
              yn = shpi(2,j,l)*dvs
              zn = shpi(3,j,l)*dvs

c             Compute d * b matrix = a

              a11 = dd(1,1,1)*xn + dd(1,4,1)*yn + dd(1,6,1)*zn
              a21 = dd(2,1,1)*xn + dd(2,4,1)*yn + dd(2,6,1)*zn
              a31 = dd(3,1,1)*xn + dd(3,4,1)*yn + dd(3,6,1)*zn
              a41 = dd(4,1,1)*xn + dd(4,4,1)*yn + dd(4,6,1)*zn
              a51 = dd(5,1,1)*xn + dd(5,4,1)*yn + dd(5,6,1)*zn
              a61 = dd(6,1,1)*xn + dd(6,4,1)*yn + dd(6,6,1)*zn

              a12 = dd(1,2,1)*yn + dd(1,4,1)*xn + dd(1,5,1)*zn
              a22 = dd(2,2,1)*yn + dd(2,4,1)*xn + dd(2,5,1)*zn
              a32 = dd(3,2,1)*yn + dd(3,4,1)*xn + dd(3,5,1)*zn
              a42 = dd(4,2,1)*yn + dd(4,4,1)*xn + dd(4,5,1)*zn
              a52 = dd(5,2,1)*yn + dd(5,4,1)*xn + dd(5,5,1)*zn
              a62 = dd(6,2,1)*yn + dd(6,4,1)*xn + dd(6,5,1)*zn

              a13 = dd(1,3,1)*zn + dd(1,5,1)*yn + dd(1,6,1)*xn
              a23 = dd(2,3,1)*zn + dd(2,5,1)*yn + dd(2,6,1)*xn
              a33 = dd(3,3,1)*zn + dd(3,5,1)*yn + dd(3,6,1)*xn
              a43 = dd(4,3,1)*zn + dd(4,5,1)*yn + dd(4,6,1)*xn
              a53 = dd(5,3,1)*zn + dd(5,5,1)*yn + dd(5,6,1)*xn
              a63 = dd(6,3,1)*zn + dd(6,5,1)*yn + dd(6,6,1)*xn

              do i = 1,3
              xn = shpi(1,i,l)
                yn = shpi(2,i,l)
                zn = shpi(3,i,l)
                hh(3*i-2,3*j-2) = hh(3*i-2,3*j-2) + xn*a11 + yn*a41
     &                                            + zn*a61
                hh(3*i-2,3*j-1) = hh(3*i-2,3*j-1) + xn*a12 + yn*a42
     &                                            + zn*a62
                hh(3*i-2,3*j  ) = hh(3*i-2,3*j  ) + xn*a13 + yn*a43
     &                                            + zn*a63
                hh(3*i-1,3*j-2) = hh(3*i-1,3*j-2) + yn*a21 + xn*a41
     &                                            + zn*a51
                hh(3*i-1,3*j-1) = hh(3*i-1,3*j-1) + yn*a22 + xn*a42
     &                                            + zn*a52
                hh(3*i-1,3*j  ) = hh(3*i-1,3*j  ) + yn*a23 + xn*a43
     &                                            + zn*a53
                hh(3*i  ,3*j-2) = hh(3*i  ,3*j-2) + zn*a31 + yn*a51
     &                                            + xn*a61
                hh(3*i  ,3*j-1) = hh(3*i  ,3*j-1) + zn*a32 + yn*a52
     &                                            + xn*a62
                hh(3*i  ,3*j  ) = hh(3*i  ,3*j  ) + zn*a33 + yn*a53
     &                                            + xn*a63
              end do ! i
            end do ! j
            nn = nn + nhv
          end do ! l
          ha = hr(nh2+1)  ! Save augmented function

c         Eliminate enhanced modes

          call invert(hh,9,9)

c         Compute enhanced parameters

          do i = 1,9
            dui(i) = 0.0d0
            do j = 1,9
              dui(i) = dui(i) + hh(i,j)*bb(j)
            end do ! j
          end do ! i

c         Check convergence

          uinorm = unorm
          duinrm = 0.0d0
          do i = 1,9
            uinorm = max(uinorm,abs(ui(i,1,1)))
            duinrm = max(duinrm,abs(dui(i)))
          end do !

          icount = icount + 1
          noconv = duinrm.gt.utol*uinorm .and. icount.lt.ncount

          if(noconv) then
            do i = 1,9
              ui(i,1,1) = ui(i,1,1) + dui(i)
              ui(i,1,2) = ui(i,1,2) + dui(i)
            end do ! i
          else
            do i = 1,9
              hr(nh2+i+1) = ui(i,1,1)
            end do ! i
          end if ! noconv

        end do ! while

c       Stiffness and/or Residual computation

        if(isw.eq.3 .or. isw.eq.6) then

c         Set mass factors

          rr    = d(4)
          if((d(7).ge.0.0 .or. d(183).ne.0.0d0) .and.
     &               (ndfo(1).gt.0 .or. shflg)) then
            cfac = d(7)
            lfac = 1.d0 - cfac
          else
            cfac = 0.0d0
            lfac = 0.0d0
          endif

c         Initialize coupling matrix

          do i = 1,24
            do j = 1,9
              gg(j,i) = 0.0d0
            end do ! j
          end do ! i

          nn   = nhi
          alam = hr(nh2)
          do l = 1,lint

c           Compute stress and strain at point

            call strnen(d,xl,ul,tl,shp3(1,1,l),shpi(1,1,l),
     &                  ndm,ndf,nel, eps,ta)

c           Compute constitution

            call modlsd(l,d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                  dd,sig(1,l),alam,ha,isw)

c           Compute parameters for multiscale plane strain

            v_avg = v_avg + jac(l)
            v_rho = v_rho + jac(l)*d(4)

c           Store time history plot data for element

            i = 6*(l-1)
            do j = 1,6
              tt(j+i) = sig(j,l)
            end do ! j

c           Add stiffness part of Rayleigh damping to stress

            if(d(78).ne.0.0d0) then
              call rays3d(d,shp3(1,1,l),shp3(1,1,l),sig(1,l),dd,
     &                    ul(1,1,4),ndf,nel,.false.)
            endif

c           Residual computations

            call resid3d(jac(l),shp3(1,1,l),sig(1,l),d, xl,ul(1,1,4),
     &                   ul(1,1,5),r,ndm,ndf, .true.)

c           Stiffness computations

            if(isw.eq.3) then

              dvs = (ctan(1) + d(78)*ctan(2))*jac(l)
              dvm = (ctan(3) + d(77)*ctan(2))*jac(l)*rr

              j1  = 1
              do j = 1,nel

c               Compute d * b matrix = a

                xn  = shp3(1,j,l)*dvs
                yn  = shp3(2,j,l)*dvs
                zn  = shp3(3,j,l)*dvs

                a11 = dd(1,1,1)*xn + dd(1,4,1)*yn + dd(1,6,1)*zn
                a21 = dd(2,1,1)*xn + dd(2,4,1)*yn + dd(2,6,1)*zn
                a31 = dd(3,1,1)*xn + dd(3,4,1)*yn + dd(3,6,1)*zn
                a41 = dd(4,1,1)*xn + dd(4,4,1)*yn + dd(4,6,1)*zn
                a51 = dd(5,1,1)*xn + dd(5,4,1)*yn + dd(5,6,1)*zn
                a61 = dd(6,1,1)*xn + dd(6,4,1)*yn + dd(6,6,1)*zn

                a12 = dd(1,2,1)*yn + dd(1,4,1)*xn + dd(1,5,1)*zn
                a22 = dd(2,2,1)*yn + dd(2,4,1)*xn + dd(2,5,1)*zn
                a32 = dd(3,2,1)*yn + dd(3,4,1)*xn + dd(3,5,1)*zn
                a42 = dd(4,2,1)*yn + dd(4,4,1)*xn + dd(4,5,1)*zn
                a52 = dd(5,2,1)*yn + dd(5,4,1)*xn + dd(5,5,1)*zn
                a62 = dd(6,2,1)*yn + dd(6,4,1)*xn + dd(6,5,1)*zn

                a13 = dd(1,3,1)*zn + dd(1,5,1)*yn + dd(1,6,1)*xn
                a23 = dd(2,3,1)*zn + dd(2,5,1)*yn + dd(2,6,1)*xn
                a33 = dd(3,3,1)*zn + dd(3,5,1)*yn + dd(3,6,1)*xn
                a43 = dd(4,3,1)*zn + dd(4,5,1)*yn + dd(4,6,1)*xn
                a53 = dd(5,3,1)*zn + dd(5,5,1)*yn + dd(5,6,1)*xn
                a63 = dd(6,3,1)*zn + dd(6,5,1)*yn + dd(6,6,1)*xn

c               Add diagonal mass effects

                am           = shp3(4,j,l)*dvm
                mass         = am*lfac
                s(j1  ,j1  ) = s(j1  ,j1  ) + mass
                s(j1+1,j1+1) = s(j1+1,j1+1) + mass
                s(j1+2,j1+2) = s(j1+2,j1+2) + mass

                i1 = 1
                do i = 1,nel
                  xn   = shp3(1,i,l)
                  yn   = shp3(2,i,l)
                  zn   = shp3(3,i,l)
                  mass = shp3(4,i,l)*am*cfac

                  s(i1  ,j1  ) = s(i1  ,j1  ) + xn*a11 + yn*a41 + zn*a61
     &                                        + mass
                  s(i1  ,j1+1) = s(i1  ,j1+1) + xn*a12 + yn*a42 + zn*a62
                  s(i1  ,j1+2) = s(i1  ,j1+2) + xn*a13 + yn*a43 + zn*a63

                  s(i1+1,j1  ) = s(i1+1,j1  ) + yn*a21 + xn*a41 + zn*a51
                  s(i1+1,j1+1) = s(i1+1,j1+1) + yn*a22 + xn*a42 + zn*a52
     &                                        + mass
                  s(i1+1,j1+2) = s(i1+1,j1+2) + yn*a23 + xn*a43 + zn*a53

                  s(i1+2,j1  ) = s(i1+2,j1  ) + zn*a31 + yn*a51 + xn*a61
                  s(i1+2,j1+1) = s(i1+2,j1+1) + zn*a32 + yn*a52 + xn*a62
                  s(i1+2,j1+2) = s(i1+2,j1+2) + zn*a33 + yn*a53 + xn*a63
     &                                        + mass
                  i1 = i1 + ndf
                end do ! i

                do i = 1,3
                  xn = shpi(1,i,l)
                  yn = shpi(2,i,l)
                  zn = shpi(3,i,l)
                  gg(3*i-2,j1  ) = gg(3*i-2,j1  ) + xn*a11+yn*a41+zn*a61
                  gg(3*i-2,j1+1) = gg(3*i-2,j1+1) + xn*a12+yn*a42+zn*a62
                  gg(3*i-2,j1+2) = gg(3*i-2,j1+2) + xn*a13+yn*a43+zn*a63
                  gg(3*i-1,j1  ) = gg(3*i-1,j1  ) + yn*a21+xn*a41+zn*a51
                  gg(3*i-1,j1+1) = gg(3*i-1,j1+1) + yn*a22+xn*a42+zn*a52
                  gg(3*i-1,j1+2) = gg(3*i-1,j1+2) + yn*a23+xn*a43+zn*a53
                  gg(3*i  ,j1  ) = gg(3*i  ,j1  ) + zn*a31+yn*a51+xn*a61
                  gg(3*i  ,j1+1) = gg(3*i  ,j1+1) + zn*a32+yn*a52+xn*a62
                  gg(3*i  ,j1+2) = gg(3*i  ,j1+2) + zn*a33+yn*a53+xn*a63
                end do ! i
                j1 = j1 + ndf
              end do ! j
            end if ! isw.eq.3

            nn = nn + nhv
          end do ! l

c         Construct static condensation

          if(isw.eq.3) then
            do i = 1,3
              is(i   ) = i
              is(i+ 3) = i + ndf
              is(i+ 6) = i + ndf*2
              is(i+ 9) = i + ndf*3
              is(i+12) = i + ndf*4
              is(i+15) = i + ndf*5
              is(i+18) = i + ndf*6
              is(i+21) = i + ndf*7
            end do ! i
            do i = 1,9
              do j = 1,24
                hg(i,j) = dot(hh(1,i),gg(1,j),9)
              end do ! j
            end do ! i
            do j = 1,24
              do i = 1,24
                s(is(i),is(j)) = s(is(i),is(j)) - dot(gg(1,i),hg(1,j),9)
              end do ! i
            end do ! j
          endif

c       Compute and output element variables

        elseif(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16) then

c         Set initial counter for history terms in stress/strain

          call shpm3d(shpi,xl,ndm,.true.)

          nn = 0
          do l = 1,lint

c           Compute stress and strain at point

            call strnen(d,xl,ul,tl,shp3(1,1,l),shpi(1,1,l),
     &                  ndm,ndf,nel, eps,ta)
            epsv(l) = eps(1,1) + eps(2,1) + eps(3,1)

c           Compute constitution

            call modlsd(l,d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                  dd,sig(1,l),alam,ha,isw)

c           Compute coordinates

            if(isw.eq.4) then

              xn = 0.0
              yn = 0.0
              zn = 0.0
              do j = 1,nel
                xn = xn + shp3(4,j,1)*xl(1,j)
                yn = yn + shp3(4,j,1)*xl(2,j)
                zn = zn + shp3(4,j,1)*xl(3,j)
              end do ! j

c             Output stress values

              do j = 1,3
                epsm(j  ) = eps(j  ,1)
                epsm(j+3) = eps(j+3,1)*0.5d0
              end do ! j
              call pstr3d(sig(1,l),psig)
              call pstr3d(epsm    ,peps)

              mct = mct - 5
              if(mct.le.0) then
                write(iow,2010) o,head
                if(ior.lt.0) write(*,2010) o,head
                mct = 50
              endif
              write(iow,2011) n,ma,xn,yn,zn,
     &                        (sig(i,l),i=1,6),(epsm(i),i=1,6),
     &                        (psig(i),i=1,3),(peps(i),i=1,3)
              if(ior.lt.0) then
                write(*,2011) n,ma,xn,yn,zn,
     &                        (sig(i,l),i=1,6),(epsm(i),i=1,6),
     &                        (psig(i),i=1,3),(peps(i),i=1,3)
              endif
            elseif(isw.eq.8) then
              do j = 1,6
                epsl(j,l) = eps(j,1)
              end do ! j
            end if

            nn = nn + nhv
          end do ! l

c         Plot stress values

          if(isw.eq.8) then
            call slcn3d(sig,epsl, r,s, nel)

c         Compute J-integrals and material forces

          elseif(isw.eq.16) then
            call pjint3d(ul,tl,epsv,sig,r,ndf,ndm)
          endif
        endif
      endif

c     Formats

2010  format(a1,20a4//5x,'Element Stresses'//'     Elmt Mat',
     &   '    1-coord    2-coord    3-coord'/12x,
     &   '  11-stress  22-stress  33-stress  12-stress',
     &   '  23-stress  31-stress'/12x,
     &   '  11-strain  22-strain  33-strain  12-strain',
     &   '  23-strain  31-strain'/12x,
     &   '   1-stress   2-stress   3-stress',
     &   '   1-strain   2-strain   3-strain')

2011  format(/i9,i4,1p,3e11.3/(12x,1p,6e11.3))

      end
