c$Id:$
      subroutine sld3d1tm(d,ul,xl,tl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    06/02/2010
c       1. Add isw to call of iner3d                        27/01/2011
c       2. Add 'l' to modlsd call                           05/01/2012
c       3. Introduce epsp(6,125) to plot strain contours    15/08/2013
c          also pass to slcn3d.
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 3-D linear coupled displacment element for feap

c     Output records:

c     Prints in element: sig-11, sig-22, sig-33, sig-12, sig-23, sig-31
c                        eps-11, eps-22, eps-33, eps-12, eps-23, eps-31
c                         1-sig,  2-sig,  3-sig

c     Prints at nodes:   1=sig-11, 2=sig-22, 3=sig-33,
c                        4=sig-12  5=sig-23, 6=sig-31
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'complx.h'
      include  'elcoor.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'prstrs.h'
      include  'ptdat6.h'
      include  'qudshp.h'
      include  'rdata.h'
      include  'comblk.h'

      integer   i,j,l,nn,i1,j1,ndf,ndm,nst,isw,nhv,istrt
      real*8    sig(10,125),eps(9,3),dd(6,6,5),vv(3)
      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),tl(*),s(nst,*),r(*)
      real*8    th(125),epsv(125),epsp(6,125)
      real*8    mom(3,64), psig(3),epp(6),peps(3)
      real*8    dvm, xn, yn, zn, ta, lfac, cfac, sfac
      real*8    a11, a12, a13, a14, a15, a16
      real*8    a21, a22, a23, a24, a25, a26
      real*8    a31, a32, a33, a34, a35, a36
      real*8    bt1, bt2, bt3, sh3
      real*8    alam, ha

      save

      data      alam,ha / 2*0.0d0 /

c     Input modifications

      if(isw.eq.1) then
        return
      endif

c     Compute element tangent array

      nhv   = nint(d(15))
      istrt = nint(d(84))

c     Get quadrature information

      call quadr3d(d,.true.)

c     Compute shape functions

      do l = 1,lint
        call interp3d(l, xl, ndm,nel)
      end do ! l

c     Compute the residual and tangent arrays or energy

      if(isw.eq.3 .or. isw.eq.6 .or. isw.eq.13) then

c       Compute intertial effects

        if(ctan(3).ne.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
          cfac = d(7)
          lfac = 1.0d0 - cfac
          call iner3d(d,xl,ul(1,1,4),ul(1,1,5),s,r,
     &                nel,ndf,ndm,nst, isw)
        else
          cfac = 0.0d0
          lfac = 0.0d0
        endif

c       Initialize momenta

        if(isw.eq.13) then
          do j = 1,nel
            do i = 1,3
              mom(i,j) = 0.0d0
            end do ! i
          end do ! j
        endif

c       Compute residual and tangent arrays

        nn = 0
        do l = 1,lint

c         Compute strain at point

          call strn3d(d,xl,ul,th,shp3(1,1,l),3,ndf,nel, eps,ta)

c         Compute stress at point

          call modlsd(l,d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                dd,sig(1,l),alam,ha,isw)

c         Residual and tangent computation

          if(isw.eq.3 .or. isw.eq.6) then

c           Residual computations

c           Store time history plot data for element

            i = 6*(l-1)
            do j = 1,6
              tt(j+i) = sig(j,l)
            end do ! j

c           Add stiffness part of Rayleigh damping to stress

            if(d(78).ne.0.0d0) then
              call rays3d(d,shp3(1,1,l),shp3(1,1,l),sig(1,l),
     &                    dd(1,1,5),ul(1,1,4),ndf,nel,.false.)
              sfac = d(78)*ctan(2)
            else
              sfac = 0.0d0
            endif

c           Form residual

            call resid3d(jac(l),shp3(1,1,l),sig(1,l),d,xl,
     &                   ul(1,1,4),ul(1,1,5),r,ndm,ndf, .false.)

c           Thermal Velocity

            do i = 1,ndm
              vv(i) = 0.0d0
              do j = 1,nel
                vv(i) = vv(i) + shp3(4,j,l)*ul(i,j,4)
              end do ! j
            end do ! i

            j1 = 4
            do j = 1,nel

c             Compute B_trans * D * alpha * j * w

              xn  = shp3(1,j,l)*jac(l)*d(129)
              yn  = shp3(2,j,l)*jac(l)*d(129)
              zn  = shp3(3,j,l)*jac(l)*d(129)
              bt1 = xn*dd(1,1,2) + yn*dd(4,1,2) + zn*dd(6,1,2)
              bt2 = xn*dd(4,1,2) + yn*dd(2,1,2) + zn*dd(5,1,2)
              bt3 = xn*dd(6,1,2) + yn*dd(5,1,2) + zn*dd(3,1,2)

              r(j1) = r(j1) - bt1*vv(1) - bt2*vv(2) - bt3*vv(3)

              j1  = j1 + ndf
            end do ! j

c           Stiffness computations

            if(isw.eq.3) then

c             Modify tangent for stiffness Rayleigh damping

              do j = 1,6
                do i = 1,6
                  dd(i,j,1) = dd(i,j,1)*ctan(1) + dd(i,j,5)*sfac
                end do ! i
              end do ! j

              i1 = 1
              do i = 1,nel

c               Compute d * b matrix = a

                xn  = shp3(1,i,l)*jac(l)
                yn  = shp3(2,i,l)*jac(l)
                zn  = shp3(3,i,l)*jac(l)

                a11 = xn*dd(1,1,1) + yn*dd(4,1,1) + zn*dd(6,1,1)
                a12 = xn*dd(1,2,1) + yn*dd(4,2,1) + zn*dd(6,2,1)
                a13 = xn*dd(1,3,1) + yn*dd(4,3,1) + zn*dd(6,3,1)
                a14 = xn*dd(1,4,1) + yn*dd(4,4,1) + zn*dd(6,4,1)
                a15 = xn*dd(1,5,1) + yn*dd(4,5,1) + zn*dd(6,5,1)
                a16 = xn*dd(1,6,1) + yn*dd(4,6,1) + zn*dd(6,6,1)

                a21 = xn*dd(4,1,1) + yn*dd(2,1,1) + zn*dd(5,1,1)
                a22 = xn*dd(4,2,1) + yn*dd(2,2,1) + zn*dd(5,2,1)
                a23 = xn*dd(4,3,1) + yn*dd(2,3,1) + zn*dd(5,3,1)
                a24 = xn*dd(4,4,1) + yn*dd(2,4,1) + zn*dd(5,4,1)
                a25 = xn*dd(4,5,1) + yn*dd(2,5,1) + zn*dd(5,5,1)
                a26 = xn*dd(4,6,1) + yn*dd(2,6,1) + zn*dd(5,6,1)

                a31 = xn*dd(6,1,1) + yn*dd(5,1,1) + zn*dd(3,1,1)
                a32 = xn*dd(6,2,1) + yn*dd(5,2,1) + zn*dd(3,2,1)
                a33 = xn*dd(6,3,1) + yn*dd(5,3,1) + zn*dd(3,3,1)
                a34 = xn*dd(6,4,1) + yn*dd(5,4,1) + zn*dd(3,4,1)
                a35 = xn*dd(6,5,1) + yn*dd(5,5,1) + zn*dd(3,5,1)
                a36 = xn*dd(6,6,1) + yn*dd(5,6,1) + zn*dd(3,6,1)

c               Compute B_trans * D * alpha * j * w

                bt1  = xn*dd(1,1,2) + yn*dd(4,1,2) + zn*dd(6,1,2)
                bt2  = xn*dd(4,1,2) + yn*dd(2,1,2) + zn*dd(5,1,2)
                bt3  = xn*dd(6,1,2) + yn*dd(5,1,2) + zn*dd(3,1,2)

                j1 = 1
                do j = 1,nel

c                 Compute stiffness matrix

                  xn = shp3(1,j,l)
                  yn = shp3(2,j,l)
                  zn = shp3(3,j,l)

                  s(i1  ,j1  ) = s(i1  ,j1  ) + a11*xn + a14*yn + a16*zn
                  s(i1  ,j1+1) = s(i1  ,j1+1) + a14*xn + a12*yn + a15*zn
                  s(i1  ,j1+2) = s(i1  ,j1+2) + a16*xn + a15*yn + a13*zn

                  s(i1+1,j1  ) = s(i1+1,j1  ) + a21*xn + a24*yn + a26*zn
                  s(i1+1,j1+1) = s(i1+1,j1+1) + a24*xn + a22*yn + a25*zn
                  s(i1+1,j1+2) = s(i1+1,j1+2) + a26*xn + a25*yn + a23*zn

                  s(i1+2,j1  ) = s(i1+2,j1  ) + a31*xn + a34*yn + a36*zn
                  s(i1+2,j1+1) = s(i1+2,j1+1) + a34*xn + a32*yn + a35*zn
                  s(i1+2,j1+2) = s(i1+2,j1+2) + a36*xn + a35*yn + a33*zn

c                 Add thermo-mechanical coupling term

                  sh3          = shp3(4,j,l)*ctan(1)
                  s(i1  ,j1+3) = s(i1  ,j1+3) + bt1*sh3
                  s(i1+1,j1+3) = s(i1+1,j1+3) + bt2*sh3
                  s(i1+2,j1+3) = s(i1+2,j1+3) + bt3*sh3

c                 Rate terms

                  if(d(129).ne.0.0d0) then
                    sh3            = shp3(4,j,l)*d(129)*ctan(2)
                    s(i1+3,j1  ) = s(i1+3,j1  ) + bt1*sh3
                    s(i1+3,j1+1) = s(i1+3,j1+1) + bt2*sh3
                    s(i1+3,j1+2) = s(i1+3,j1+2) + bt3*sh3
                  endif

                  j1 = j1 + ndf
                end do ! j
                i1 = i1 + ndf
              end do ! i
            endif
            if(cplxfl) then
              call cst3d1(shp3(1,1,l),jac(l),eps,ul(1,1,8),
     &                    dd,dd(1,1,3),s(1,nst+1),r,r(nst+1))
            endif

c         Compute energy

          elseif(isw.eq.13) then

            dvm  = d(4)*jac(l)
            call sengy(ul,shp3(1,1,l), jac(l),dvm, lfac,cfac, 3,4, mom)

          endif ! isw
          nn = nn + nhv
        end do ! l

c       Accumulate momentum

        if(isw.eq.13) then
          do j = 1,nel
            do i = 1,3
              epl(i) = epl(i) + mom(i,j)
            end do ! i
            epl(4) = epl(4) + xl(2,j)*mom(3,j) - xl(3,j)*mom(2,j)
            epl(5) = epl(5) + xl(3,j)*mom(1,j) - xl(1,j)*mom(3,j)
            epl(6) = epl(6) + xl(1,j)*mom(2,j) - xl(2,j)*mom(1,j)
          end do ! j
        endif

c     Compute and output element variables

      elseif(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16) then

c       Set initial counter for history terms in stress/strain

        nn = 0
        do l = 1,lint

c         Compute strain at point

          call strn3d(d,xl,ul,th,shp3(1,1,l),ndm,ndf,nel, eps,ta)
          epsv(l) = eps(1,1) + eps(2,1) + eps(3,1)

c         Compute stress at point

          call modlsd(l,d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                dd,sig(1,l),alam,ha,isw)

C         Output values

          if(isw.eq.4) then

c           Compute coordinates

            xn = 0.0
            yn = 0.0
            zn = 0.0
            do j = 1,nel
              xn = xn + shp3(4,j,l)*xl(1,j)
              yn = yn + shp3(4,j,l)*xl(2,j)
              zn = zn + shp3(4,j,l)*xl(3,j)
            end do ! j

c           Compute principal stress values

            do i = 1,3
              epp(i  ) = eps(i  ,1)
              epp(i+3) = eps(i+3,1)*0.5d0
            end do ! i
            call pstr3d(sig(1,l),psig)
            call pstr3d(epp     ,peps)

            mct = mct - 5
            if(mct.le.0) then
              write(iow,2010) o,head
              if(ior.lt.0) write(*,2010) o,head
              mct = 50
            endif
            write(iow,2011) n,ma,xn,yn,zn,
     &                      (sig(i,l),i=1,6),(epp(i),i=1,6),
     &                      (psig(i),i=1,3),(peps(i),i=1,3)
            if(ior.lt.0) then
              write(*,2011) n,ma,xn,yn,zn,
     &                      (sig(i,l),i=1,6),(epp(i),i=1,6),
     &                      (psig(i),i=1,3),(peps(i),i=1,3)
            endif
          elseif(isw.eq.8) then
            do i = 1,6
              epsp(i,l) = eps(i,1)
            end do ! i
          endif
          nn = nn + nhv
        end do ! l

c       Plot stress values

        if(isw.eq.8) then
          call slcn3d(sig,epsp, r,s, nel)

c       Compute J-integrals and material forces

        elseif(isw.eq.16) then
          call pjint3d(ul,tl,epsv,sig,r,ndf,ndm)
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
     &   '   1-strain   2-strain   3-strain',/39(' -'))

2011  format(/i9,i4,1p,3e11.3/(12x,1p,6e11.3))

      end
