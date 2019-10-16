c$Id:$
      subroutine sld3d1(d,ul,xl,tl,s,r,ndf,ndm,nst,isw, ther)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct computation of 4-node tet inertia        29/03/2007
c       2. Add cylindrical material property option         15/06/2007
c       3. Add 14/15 node tet computation                   06/07/2007
c       4. Change 'ord' to 'nel' on 'tetshp' calls          05/11/2007
c       5. Modify quadrature for interia, add 'iner3d'      28/12/2007
c       6. Add 'psig' & 'peps'  to output principal values  04/03/2008
c       7. Add output of 'zn', change formats to be same    12/06/2008
c          for all sld3d*.f files.
c       8. Add direct call to quadr3d and interp3d          11/11/2008
c       9. Remove 'nel' from call to 'quadr3d'              23/01/2009
c      10. Increase arrays to store 64 node brick           03/02/2009
c      11. Add thermal coupling terms                       30/04/2009
c      12. Add 'elcoor' to use xref                         25/07/2009
c      13. Move computation of psil to material models      26/07/2009
c      14. Increase arrays to store 125 node brick          20/12/2010
c      15. Move inertial computation after stiffness        01/01/2011
c          Add isw to iner3d to control tangent computation
c      16. Add 'l' to modlsd call                           05/01/2012
c      17. Add isw.eq.14 to execution option                04/05/2012
c      18. Add average of density for multiscale            10/05/2012
c      19. Add eps to slcn3d call                           01/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 3-D linear elastic displacment element for feap

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
      include  'debugs.h'
      include  'elcoor.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'oelmt.h'
      include  'part0.h'
      include  'prstrs.h'
      include  'ptdat6.h'
      include  'qudshp.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   ther
      integer   i,j,l,nn,i1,j1,ndf,ndm,nst,isw,nhv,istrt,tdof,tdf
      real*8    sig(10,125),eps(9,3),dd(6,6,5),vv(3),epsp(6,125)
      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),tl(*),s(nst,*),r(*)
      real*8    th(125),epsv(125)
      real*8    mom(3,125), psig(3),epp(6),peps(3)
      real*8    dvm, xn, yn, zn, ta, lfac, cfac, sfac
      real*8    a11, a12, a13, a14, a15, a16
      real*8    a21, a22, a23, a24, a25, a26
      real*8    a31, a32, a33, a34, a35, a36
      real*8    bt1, bt2, bt3, sh3
      real*8    alam, ha

      save

      data      alam,ha / 2*0.0d0 /

c     Set nodal temperatures: Can be specified or computed

      if(isw.gt.1) then
        tdof = nint(d(19))
        if(tdof.le.0) then
          do i = 1,nel ! {
            th(i) = tl(i)
          end do ! i     }
        else
          do i = 1,nel ! {
            th(i) = ul(tdof,i,1)
          end do ! i     }
        endif
      endif

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

      if(isw.eq.3 .or. isw.eq.6 .or. isw.eq.13 .or. isw.eq.14) then

c       Compute intertial effects

        if(ctan(3).ne.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
          cfac = d(7)
          lfac = 1.0d0 - cfac
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

c           Compute parameters for multiscale plane strain

            v_avg = v_avg + jac(l)
            v_rho = v_rho + jac(l)*d(4)

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

c           Thermal coupling property for rates

            if(ther .and. d(129).ne.0.0d0) then
              tdof = nint(d(19))

              if(tdof.gt.0) then

c               Velocity

                do i = 1,ndm
                  vv(i) = 0.0d0
                  do j = 1,nel
                    vv(i) = vv(i) + shp3(4,j,l)*ul(i,j,4)
                  end do ! j
                end do ! i

                j1 = 0
                do j = 1,nel

c                 Compute B_trans * D * alpha * j * w

                  xn  = shp3(1,j,l)*jac(l)*d(129)
                  yn  = shp3(2,j,l)*jac(l)*d(129)
                  zn  = shp3(3,j,l)*jac(l)*d(129)
                  bt1 = xn*dd(1,1,2) + yn*dd(4,1,2) + zn*dd(6,1,2)
                  bt2 = xn*dd(4,1,2) + yn*dd(2,1,2) + zn*dd(5,1,2)
                  bt3 = xn*dd(6,1,2) + yn*dd(5,1,2) + zn*dd(3,1,2)

                  r(j1+tdof) = r(j1+tdof)-bt1*vv(1)-bt2*vv(2)-bt3*vv(3)

                  j1  = j1 + ndf
                end do ! j
              endif
            endif

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

                  if(ther .and. tdof.gt.0) then
                    tdf          = tdof - 1
                    sh3          = shp3(4,j,l)*ctan(1)
                    s(i1  ,j1+tdf) = s(i1  ,j1+tdf) + bt1*sh3
                    s(i1+1,j1+tdf) = s(i1+1,j1+tdf) + bt2*sh3
                    s(i1+2,j1+tdf) = s(i1+2,j1+tdf) + bt3*sh3

c                   Rate terms

                    if(d(129).ne.0.0d0) then
                      sh3            = shp3(4,j,l)*d(129)*ctan(2)
                      s(i1+tdf,j1  ) = s(i1+tdf,j1  ) + bt1*sh3
                      s(i1+tdf,j1+1) = s(i1+tdf,j1+1) + bt2*sh3
                      s(i1+tdf,j1+2) = s(i1+tdf,j1+2) + bt3*sh3
                    endif

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

c       Compute intertial effects

        if(ctan(3).ne.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
          call iner3d(d,xl,ul(1,1,4),ul(1,1,5),s,r,
     &                nel,ndf,ndm,nst, isw)
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

      subroutine cst3d1(shp,xsj,epsr,ul,dr,di, si,pr,pi)

      implicit   none
      include   'cdata.h'
      include   'eldata.h'
      include   'eltran.h'
      include   'pmod2d.h'
      include   'sdata.h'

      integer    j,jj,j1,k,kk,k1
      real*8     shp(4,*),epsr(*),epsi(6),ul(ndf,*)
      real*8     dr(6,6),di(6,6), sigr(6),sigi(6),si(nst,*),pr(*),pi(*)
      real*8     xsj, aj1,aj2,aj3, bd(3,6)

c     Compute imaginary strains

      do j = 1,6
        epsi(j) = 0.0d0
      end do ! j
      do j = 1,nel
        epsi(1) = epsi(1) + shp(1,j)*ul(1,j)
        epsi(2) = epsi(2) + shp(2,j)*ul(2,j)
        epsi(3) = epsi(3) + shp(3,j)*ul(3,j)
        epsi(4) = epsi(4) + shp(1,j)*ul(2,j) + shp(2,j)*ul(1,j)
        epsi(5) = epsi(5) + shp(2,j)*ul(3,j) + shp(3,j)*ul(2,j)
        epsi(6) = epsi(6) + shp(3,j)*ul(1,j) + shp(1,j)*ul(3,j)
      end do ! j

c     Compute imaginary stress

      do j = 1,6
        sigr(j) = 0.0d0
        sigi(j) = 0.0d0
        do k = 1,6
          sigr(j) = sigr(j) + di(j,k)*epsi(k)
          sigi(j) = sigi(j) + dr(j,k)*epsi(k) + di(j,k)*epsr(k)
        end do ! k
      end do ! j

      j1 = 0
      do jj = 1,nel

        aj1 = shp(1,jj)*xsj
        aj2 = shp(2,jj)*xsj
        aj3 = shp(3,jj)*xsj

c       Compute B_trans * D * j * w

        do k = 1,6
          bd(1,k) = (aj1*di(1,k) + aj2*di(4,k) + aj3*di(6,k))*ctan(1)
          bd(2,k) = (aj2*di(2,k) + aj1*di(4,k) + aj3*di(5,k))*ctan(1)
          bd(3,k) = (aj3*di(3,k) + aj2*di(5,k) + aj1*di(6,k))*ctan(1)
        end do ! k

c       Loop over columns (symmetry noted)

        k1 = 0
        do kk = 1,nel
          do j = 1,3
            si(j1+j,k1+1) = si(j1+j,k1+1) + bd(j,1)*shp(1,kk)
     &                                    + bd(j,4)*shp(2,kk)
     &                                    + bd(j,6)*shp(3,kk)

            si(j1+j,k1+2) = si(j1+j,k1+2) + bd(j,4)*shp(1,kk)
     &                                    + bd(j,2)*shp(2,kk)
     &                                    + bd(j,5)*shp(3,kk)

            si(j1+j,k1+3) = si(j1+j,k1+3) + bd(j,6)*shp(1,kk)
     &                                    + bd(j,5)*shp(2,kk)
     &                                    + bd(j,3)*shp(3,kk)
          end do ! j
          k1 = k1 + ndf
        end do ! kk

c       Residual for added real part

        pr(j1+1) = pr(j1+1) + aj1*sigr(1) + aj2*sigr(4) + aj3*sigr(6)
        pr(j1+2) = pr(j1+2) + aj2*sigr(2) + aj1*sigr(4) + aj3*sigr(5)
        pr(j1+3) = pr(j1+3) + aj3*sigr(3) + aj2*sigr(5) + aj1*sigr(6)

c       Residual for imaginary part

        pi(j1+1) = pi(j1+1) - aj1*sigi(1) - aj2*sigi(4) - aj3*sigi(6)
        pi(j1+2) = pi(j1+2) - aj2*sigi(2) - aj1*sigi(4) - aj3*sigi(5)
        pi(j1+3) = pi(j1+3) - aj3*sigi(3) - aj2*sigi(5) - aj1*sigi(6)

        j1 = j1 + ndf
      end do ! jj

      end
