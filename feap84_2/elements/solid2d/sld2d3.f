c$Id:$
      subroutine sld2d3(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use d(185) for augmenting factor                 14/03/2007
c       2. Revise augment: augfa -> augf                    14/04/2009
c       3. Move computation of psil to material models      26/07/2009
c       4. Set quadrature for case when d(5).eq.0           23/11/2010
c          Add return for isw.eq.10,
c          Remove call shp2d for isw.eq.4,
c          Correct index from 5 to 7 for pstr2d call
c       5. Compute v_avg and sig_33 for multiscale use      03/12/2010
c       6. Add 'l' to modlsd call                           05/01/2012
c       7. Add average of density for multiscale            10/05/2012
c       8. Add eps on call to slcn2d                        01/01/2013
c       9. Pass strains to stcn2z for z-zhu projections     01/01/2014
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Plane enhanced strain element for FEAP

c      Output records:
c      Prints in element: sig-11, sig-22, sig-33, sig-12, sig-1 sig-2
c                         eps-11, eps-22, eps-33, eps-12

c      Prints at nodes:   1=sig-11, 2=sig-22, 3=sig-33, 4=sig-12
c                         psig-1  , psig-2    (computed by FEAP)

c      History Variable Storage (relative to nh1 or nh2)

c      Start           Description             Variable  Length
c      hr(0)           Enhanced displacement      ui(*,1)   5
c      hr(5)           Stress history at point-1    -      nhv
c      hr(5+  nhv)     Stress history at point-2    -      nhv
c      hr(5+2*nhv)     Stress history at point-3    -      nhv
c      hr(5+3*nhv)     Stress history at point-4    -      nhv

c      Total number of words / element is 5 + 4*nhv
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'
      include  'bdata.h'
      include  'cdata.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'incshp.h'
      include  'iofile.h'
      include  'oelmt.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'prld1.h'
      include  'qudshp.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   noconv
      integer   ndm,ndf,nst,isw, ix(*)
      real*8    alam,ha, ta,tol1,tolu
      real*8    cfac,lfac,sfac,dmas,lms,cms,jac0 ,rr,zz
      integer   i,j,l,nhv,ni,nn,nenit,i1,j1,i2,j2,ii,jj
      integer   istrt, ncp,ncs,nci
      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),tl(*),s(nst,*),p(ndf,*)
      real*8    sig(10,9),eps(9,3),epsv(9), epsp(6,64)
      real*8    gg(5,12),hh(5,5),bb(5),hg(5,12),dui(5),dd(6,6,5,9)
      real*8    ss(12,12),shpr(9,9),sigl(6),aa(6,6,9)
      real*8    bd(6,3),dvol(9),r0(9),epsm(9)

      save

      data      ni /7/, eps / 27*0.0d0 /
      data      alam,ha     /  2*0.0d0 /

c     Data inputs

      if( isw.eq.1 ) then

c       Adjust history storage

        nh1 = nh1 + 5 + 2
        return

      elseif(isw.eq.10) then

        hr(nh2) =  hr(nh2) + augf*d(185)*hr(nh2+1)
        return

      endif

c     Set for torsion case

      if(stype.eq.8) then
        ncp = 3
        ncs = 6
        nci = 12
      else
        ncp = 2
        ncs = 4
        nci = 8
      endif

c     Recover enhanced modes (saved in last iteration)

      nhv   = nint(d(15))
      istrt = nint(d(84))
      do i = 1,5
        ui(i,1)   = hr(nh2+1+i)
        ui(i,2)   = 0.0d0
      end do ! i

c     Compute quadrature and weights

      if(nint(d(182)).gt.0) then
        call int2dn(nel,lint,sg2)
      else
        l = min(3,nint(d(5)))
        if(l.eq.0) then
          l = 2
        endif
        call int2d(l,lint,sg2)
      endif

c     Initialize history variables only

      if(isw.eq.14) then
        if(stype.eq.1) then
          nn = ni + lint
        else
          nn = ni
        endif
        alam = hr(nh2)
        do l = 1,lint
          call modlsd(l,d,ta,eps,hr(nn+nh1),hr(nn+nh2),nhv,istrt,
     &                dd(1,1,1,l),sig(1,l),alam,ha,isw)
          nn = nn + nhv
        end do ! l
        return
      endif ! isw.eq.14

c     Compute shape functions

      do l = 1,lint
        call shp2d(sg2(1,l),xl,shp2(1,1,l),jac(l),ndm,nel,ix,.false.)

c       Axisymmetry

        if(stype.eq.3 .or. stype.eq.8) then
          r0(l)   = shp2(3,1,l)*xl(1,1) + shp2(3,2,l)*xl(1,2)
     &            + shp2(3,3,l)*xl(1,3) + shp2(3,4,l)*xl(1,4)
          dvol(l) = jac(l)*sg2(3,l)*r0(l)
          do j = 1,nel
            shpr(j,l) = shp2(3,j,l)/r0(l)
          end do ! j

c       Plane

        else
          dvol(l) = jac(l)*sg2(3,l)
          do j = 1,nel
            shpr(j,l) = 0.0d0
          end do ! j
        endif
      end do ! l

c     Compute enhanced modes by local iteration

      tolu   =  1.d-03*tol*rnmax/dble(numel)
      nenit  =  0
      noconv = .true.
      do while(noconv)

        do j = 1,5
          bb(j) = 0.0d0
          do i = 1,5
            hh(i,j) = 0.0d0
          end do ! i
        end do ! j

c       Set initial counter for history terms in stress/strain

        if(stype.eq.1) then
          nn = ni + lint
        else
          nn = ni
        endif
        alam = hr(nh2)
        do l = 1,lint

c         Compute enhanced strain shape functions

          call shpi2d(sg2(1,l),jac(l),xl,ndm)

c         Compute stress and strain at point

          call strn2d(d,xl,ul,tl,shp2(1,1,l),ndf,ndm,nel,rr,zz,ta,eps)

c         Plane stress modification

          if(stype.eq.1) then
            eps(3,1) = hr(nh2-1+ni+l)
            eps(3,2) = hr(nh1-1+ni+l)
          endif

          call modlsd(l,d,ta,eps,hr(nn+nh1),hr(nn+nh2),nhv,istrt,
     &                dd(1,1,1,l),sig(1,l),alam,ha,isw)

c         Store plane stress strain

          if(stype.eq.1) then
            hr(nh2-1+ni+l) = eps(3,1)
          endif

c         Scale moduli and stresses

          do i = 1,ncs
            do j = 1,ncs
              aa(j,i,l) = dd(j,i,1,l)*dvol(l)*ctan(1)
            end do ! j
            sigl(i) = sig(i,l)*dvol(l)
          end do ! i

c         Compute parameters for multiscale plane strain

          v_avg  = v_avg  + dvol(l)
          v_rho  = v_rho  + dvol(l)*d(4)
          sig_33 = sig_33 + sigl(3)

c         Store time history plot data for element

          i = 12*(l-1)
          do j = 1,6
            tt(j+i)   = sig(j,l)
            tt(j+i+6) = eps(j,1)
          end do ! j

c         Enhanced residual computations

          do j = 1,2
            bb(2*j-1) = bb(2*j-1) - sigl(1)*shpi(1,j)
     &                            - sigl(4)*shpi(2,j)
            bb(2*j  ) = bb(2*j  ) - sigl(2)*shpi(2,j)
     &                            - sigl(4)*shpi(1,j)
          end do ! j

          if(stype.eq.3 .or. stype.eq.8) then
            shpi(3,3) = shpi(3,3)/r0(l)
            bb(5)     = bb(5) - sigl(3)*shpi(3,3)
            do j = 1,2
              bd(3,1)     = aa(3,1,l)*shpi(1,j) + aa(3,4,l)*shpi(2,j)
              bd(3,2)     = aa(3,2,l)*shpi(2,j) + aa(3,4,l)*shpi(1,j)
              hh(5,2*j-1) = hh(5,2*j-1) + shpi(3,3)*bd(3,1)
              hh(5,2*j  ) = hh(5,2*j  ) + shpi(3,3)*bd(3,2)
              bd(1,1)     = aa(1,3,l)*shpi(3,3)
              bd(2,1)     = aa(2,3,l)*shpi(3,3)
              bd(4,1)     = aa(4,3,l)*shpi(3,3)
              hh(2*j-1,5) = hh(2*j-1,5) + shpi(1,j)*bd(1,1)
     &                                  + shpi(2,j)*bd(4,1)
              hh(2*j  ,5) = hh(2*j  ,5) + shpi(2,j)*bd(2,1)
     &                                  + shpi(1,j)*bd(4,1)
            end do ! j
            hh(5,5) = hh(5,5) + shpi(3,3)*aa(3,3,l)*shpi(3,3)
          else
            shpi(3,3) = 0.0d0
            hh(5,5)   = 1.0d0
          endif

c         Stiffness computations

          do j = 1,2

c           Compute d * b matrix = a

            bd(1,1) = aa(1,1,l)*shpi(1,j) + aa(1,4,l)*shpi(2,j)
            bd(2,1) = aa(2,1,l)*shpi(1,j) + aa(2,4,l)*shpi(2,j)
            bd(4,1) = aa(4,1,l)*shpi(1,j) + aa(4,4,l)*shpi(2,j)
            bd(1,2) = aa(1,2,l)*shpi(2,j) + aa(1,4,l)*shpi(1,j)
            bd(2,2) = aa(2,2,l)*shpi(2,j) + aa(2,4,l)*shpi(1,j)
            bd(4,2) = aa(4,2,l)*shpi(2,j) + aa(4,4,l)*shpi(1,j)
            do i = 1,2
              hh(2*i-1,2*j-1) = hh(2*i-1,2*j-1) + shpi(1,i)*bd(1,1)
     &                                          + shpi(2,i)*bd(4,1)
              hh(2*i-1,2*j  ) = hh(2*i-1,2*j  ) + shpi(1,i)*bd(1,2)
     &                                          + shpi(2,i)*bd(4,2)
              hh(2*i  ,2*j-1) = hh(2*i  ,2*j-1) + shpi(2,i)*bd(2,1)
     &                                          + shpi(1,i)*bd(4,1)
              hh(2*i  ,2*j  ) = hh(2*i  ,2*j  ) + shpi(2,i)*bd(2,2)
     &                                          + shpi(1,i)*bd(4,2)
            end do ! i
          end do ! j
          nn = nn + nhv
        end do ! l
        ha = hr(nh2+1)

        call invert(hh,5,5)

c       Compute incremental enhanced displacements enhanced modes

        do i = 1,5
          dui(i)  = hh(i,1)*bb(1) + hh(i,2)*bb(2) + hh(i,3)*bb(3)
     &            + hh(i,4)*bb(4) + hh(i,5)*bb(5)
          ui(i,1) = ui(i,1) + dui(i)
          ui(i,2) = ui(i,2) + dui(i)
        end do ! i

c       Check convergence

        tol1 = abs(bb(1)*dui(1) + bb(2)*dui(2) + bb(3)*dui(3)
     &           + bb(4)*dui(4) + bb(5)*dui(5))

        if(tol1.le.tolu .and. nenit.ge.1) then
          noconv = .false.
        endif
        nenit = nenit + 1
        if(nenit.ge.3 .or. tolu.eq.0.0d0) then
          noconv = .false.
        endif
      end do ! while

c     Save enhanced modes

      do i = 1,5
        hr(nh2+1+i) = ui(i,1)
      end do ! i

c     Set initial counter for history terms in stress/strain

      if(isw.eq.3. or. isw.eq.6) then

c       Time integration order set to static or dynamic

        if((d(7).ge.0.0 .or. d(183).ne.0.0d0) .and.
     &             (ndfo(1).gt.0 .or. shflg)) then
          cfac = d(7)
          lfac = 1.0d0 - cfac
        else
          cfac = 0.0d0
          lfac = 0.0d0
        endif

        do i = 1,nci
          do j = 1,5
            gg(j,i) = 0.0d0
          end do ! j
        end do ! i

        do l = 1,lint

          call shpi2d(sg2(1,l),jac(l),xl,ndm)

          if(stype.eq.3 .or. stype.eq.8) then
            shpi(3,3) = shpi(3,3)/r0(l)
            jac0      = jac(l)*sg2(3,l)
          else
            shpi(3,3) = 0.0d0
            jac0      = 0.0d0
          endif

c         Rayleigh damping effects

          if(d(78).ne.0.0d0) then
            call rays2d(d,shp2(1,1,l),sig(1,l),dd(1,1,5,l),ul(1,1,4),
     &                  xl,ndf,ndm,nel)
            sfac = d(78)*ctan(2)
          else
            sfac = 0.0d0
          endif

c         Residual computations

          call resid2d(cfac,lfac,dvol(l),jac0,rr,zz,shp2(1,1,l),
     &                 sig(1,l),d,ul(1,1,4),ul(1,1,5), p,ndf)

c         Stiffness computations

          if(isw.eq.3) then

            dmas = d(4)*(ctan(3) + d(77)*ctan(2))*dvol(l)

            do j = 1,ncs
              do i = 1,ncs
                aa(i,j,l) = aa(i,j,l) + dd(i,j,5,l)*sfac
              end do ! i
            end do ! j

            j1 = 0
            j2 = 0
            do j = 1,nel

c             Compute d * b matrix = a

              do i = 1,ncs
                bd(i,1) = aa(i,1,l)*shp2(1,j,l) + aa(i,4,l)*shp2(2,j,l)
     &                  + aa(i,3,l)*shpr(j,l)
                bd(i,2) = aa(i,2,l)*shp2(2,j,l) + aa(i,4,l)*shp2(1,j,l)
              end do ! i

c             Torsion

              if(stype.eq.8) then
                do i = 1,ncs
                  bd(i,3) = aa(i,5,l)*shp2(2,j,l)
     &                    + aa(i,6,l)*(shp2(1,j,l) - shpr(j,l))
                end do ! i
              endif

c             Lumped mass effects

              lms = shp2(3,j,l)*dmas
              do jj = 1,ncp
                s(j1+jj,j1+jj) = s(j1+jj,j1+jj) + lms*lfac
              end do ! jj

              i1 = 0
              do i = 1,nel
                cms = lms*shp2(3,i,l)*cfac
                do jj = 1,ncp
                  s(i1+jj,j1+jj) = s(i1+jj,j1+jj) + cms
                  s(i1+1 ,j1+jj) = s(i1+1 ,j1+jj) + shp2(1,i,l)*bd(1,jj)
     &                                            + shp2(2,i,l)*bd(4,jj)
     &                                            + shpr(i,l) *bd(3,jj)
                  s(i1+2 ,j1+jj) = s(i1+2 ,j1+jj) + shp2(2,i,l)*bd(2,jj)
     &                                            + shp2(1,i,l)*bd(4,jj)
                end do ! jj
                if(stype.eq.8) then
                  do jj = 1,ncp
                    s(i1+3,j1+jj) = s(i1+3,j1+jj) + shp2(2,i,l)*bd(5,jj)
     &                              + (shp2(1,i,l) - shpr(i,l))*bd(6,jj)
                  end do ! jj
                endif
                i1 = i1 + ndf
              end do ! i

c             Enhanced coupling array

              do jj = 1,ncp
                do i = 1,2
                  gg(2*i-1,j2+jj) = gg(2*i-1,j2+jj) + shpi(1,i)*bd(1,jj)
     &                                              + shpi(2,i)*bd(4,jj)

                  gg(2*i  ,j2+jj) = gg(2*i  ,j2+jj) + shpi(2,i)*bd(2,jj)
     &                                              + shpi(1,i)*bd(4,jj)
                end do ! i
                gg(5,j2+jj) = gg(5,j2+jj) + shpi(3,3)*bd(3,jj)
              end do ! jj
              j1 = j1 + ndf
              j2 = j2 + ncp
            end do ! j
          endif
        end do ! l

c       Eliminate enhanced modes

        do i = 1,5
          do j = 1,nci
            hg(i,j) = hh(1,i)*gg(1,j) + hh(2,i)*gg(2,j)
     &              + hh(3,i)*gg(3,j) + hh(4,i)*gg(4,j)
     &              + hh(5,i)*gg(5,j)
          end do ! j
        end do ! i

        if(isw.eq.3) then
          do j = 1,nci
            do i = 1,nci
              ss(i,j) = gg(1,i)*hg(1,j) + gg(2,i)*hg(2,j)
     &                + gg(3,i)*hg(3,j) + gg(4,i)*hg(4,j)
     &                + gg(5,i)*hg(5,j)
            end do ! i
          end do ! j

c         Construct static condensation

          j1 = 0
          j2 = 0
          do j = 1,4
            i1 = 0
            i2 = 0
            do i = 1,4
              do jj = 1,ncp
                do ii = 1,ncp
                  s(i1+ii,j1+jj) = s(i1+ii,j1+jj) - ss(i2+ii,j2+jj)
                end do ! ii
              end do ! jj
              i1 = i1 + ndf
              i2 = i2 + ncp
            end do ! i
            j1 = j1 + ndf
            j2 = j2 + ncp
          end do ! j

c         Compute reduced residual

          j2 = 0
          do j = 1,4
            do jj = 1,ncp
              p(jj,j) = p(jj,j) - hg(1,j2+jj)*bb(1) - hg(2,j2+jj)*bb(2)
     &                          - hg(3,j2+jj)*bb(3) - hg(4,j2+jj)*bb(4)
     &                          - hg(5,j2+jj)*bb(5)
            end do ! jj
            j2     = j2 + ncp
          end do ! j
        endif

c       Multiply by thickness if not unity

        if((isw.eq.3 .or. isw.eq.6) .and. d(14).ne.1.d0) then

          do j = 1,nst
            do i = 1,nst
              s(i,j) = s(i,j)*d(14)
            end do ! i
          end do ! j
          do j = 1,nel
            do i = 1,ndf
              p(i,j) = p(i,j)*d(14)
            end do ! i
          end do ! j

        endif

c     Compute and output element variables

      elseif(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16 .or. isw.eq.25) then

c       Set initial counter for history terms in stress/strain

        if(stype.eq.1) then
          nn = ni + lint
        else
          nn = ni
        endif
        alam = hr(nh2)
        do l = 1,lint
          call shpi2d(sg2(1,l),jac(l),xl,ndm)

c         Compute stress and strain at point

          call strn2d(d,xl,ul,tl,shp2(1,1,l),ndf,ndm,nel,rr,zz,ta,eps)

c         Plane stress modification

          if(stype.eq.1) then
            eps(3,1) = hr(nh2-1+ni+l)
            eps(3,2) = hr(nh1-1+ni+l)
          endif

          epsv(l) = eps(1,1) + eps(2,1) + eps(3,1)

          call modlsd(l,d,ta,eps,hr(nn+nh1),hr(nn+nh2),nhv,istrt,
     &                dd,sig(1,l),alam,ha,isw)

c         Compute principal stress values

          if(isw.eq.4) then
            do i = 1,3
              epsm(i  ) = eps(i  ,1)
              epsm(i+3) = eps(i+3,1)*0.5d0
            end do ! i
            mct = mct - 4
            if(stype.eq.8) then
              call pstr3d(sig(1,l),sig(7,l))
              call pstr3d(epsm    ,epsm(7))
              if(mct.le.0) then
                write(iow,2001) o,head
                if(ior.lt.0) write(*,2001) o,head
                mct = 50
              endif
              write(iow,2002) n,ma,rr,zz,(sig(i,l),i=7,9),
     &                       (sig(i,l),i=1,6),(eps(i,1),i=1,6),
     &                       (epsm(i),i=7,9)
              if(ior.lt.0) then
                write(*,2002) n,ma,rr,zz,(sig(i,l),i=7,9),
     &                       (sig(i,l),i=1,6),(eps(i,1),i=1,6),
     &                       (epsm(i),i=7,9)
              endif
            else
              call pstr2d(sig(1,l),sig(7,l))
              call pstr3d(epsm    ,epsm(7))
              if(mct.le.0) then
                write(iow,2003) o,head
                if(ior.lt.0) write(*,2003) o,head
                mct = 50
              endif
              write(iow,2004) n,ma,rr,zz,(sig(i,l),i=7,9),
     &                       (sig(i,l),i=1,4),(eps(i,1),i=1,4),
     &                       (epsm(i),i=7,9)
              if(ior.lt.0) then
                write(*,2004) n,ma,rr,zz,(sig(i,l),i=7,9),
     &                       (sig(i,l),i=1,4),(eps(i,1),i=1,4),
     &                       (epsm(i),i=7,9)
              endif
            endif

c         Store strains for plots

          elseif(isw.eq.8) then

            do j = 1,6
              epsp(j,l) = eps(j,1)
            end do ! j

          endif
          nn = nn + nhv
        end do ! l

c       Plot stress values

        if(isw.eq.8) then

          call slcn2d(ix,sig,epsp,p,s,p(nen+1,1),nel,10)

c       Compute J-integrals and material forces

        elseif(isw.eq.16) then

          call pjint2d(d,ul,tl,shp2,jac,epsv,sig,p,ndf,ndm,lint,10)

c       Compute Z-Z projections

        elseif(isw.eq.25) then

          call stcn2z(xl,sig,epsp,shp2,jac,lint,ndm,nel,10)

        endif

      endif

c     Formats

2001  format(a1,20a4//5x,'Element Stresses'//'     Elmt Mat',
     &    4x,'1-coord    2-coord   1-stress   2-stress   3-stress'/
     &    2x,'11-stress  22-stress  33-stress  12-stress',
     &    2x,'23-stress  31-stress'/'11-strain  22-strain  33-strain',
     &    2x,'12-strain  23-strain  31-strain'/
     &   15x,' 1-strain   2-strain   3-strain'/39(' -'))
2002  format(i9,i4,1p,2e11.3,1p,3e11.3/13x,1p,6e11.3/13x,1p,6e11.3/
     &       13x,1p,3e11.3)

2003  format(a1,20a4//5x,'Element Stresses'//'     Elmt Mat',
     &    4x,'1-coord    2-coord   1-stress   2-stress      Angle'/
     &   15x,'11-stress  22-stress  33-stress  12-stress'/
     &   15x,'11-strain  22-strain  33-strain  12-strain'/
     &   15x,' 1-strain   2-strain   3-strain'/39(' -'))
2004  format(i9,i4,1p,2e11.3,1p,2e11.3,0p,1f11.3/
     &       13x,1p,4e11.3/13x,1p,4e11.3/
     &       13x,1p,3e11.3)

      end
