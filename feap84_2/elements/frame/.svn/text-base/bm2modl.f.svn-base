c$Id:$
      subroutine bm2modl (d,ul,forca,a,energy,shp,hn,h1,ndf,nen,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set constitutive model for 2-d frame element

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm2com.h'
      include  'ddata.h'
      include  'eldata.h'
      include  'pmod2d.h'
      include  'tdata.h'

      integer   ndf,nen, i,ii,isw,imat,nh
      real*8    dpsi,dpsih, energy, denr ,psi, cn,sn
      real*8    theta1,thetan,rtheta, phia1,phia2,phi11,phi12
      real*8    d(*),ul(ndf,nen,*),fa(3),df(3),cc(3,3,2),hn(*),h1(*)
      real*8    forca(4),forc1(2),a(4,4),b(2,2)
      real*8    defn(3),def1(3),shp(2,3),cay(2,2),dcay(2,2)
      real*8    lamn(2,2),lam1(2,2),lama(2,2),dlam(2,2)

      save

c     Compute constitutive model type

      imat = int(d(20))
      nh   = int(d(15))

c     Compute gradient

      do i = 1,3
        fa(i) = 0.0d0
        df(i) = 0.0d0
        do ii = 1,nel
          fa(i) = fa(i) + ul(i,ii,1)*shp(1,ii)
          df(i) = df(i) + ul(i,ii,2)*shp(1,ii)
        end do ! ii
      end do ! i
      fa(1) = fa(1) + 1.d0

c     Conserving model form

      if(imat.eq.6) then

c       Compute rotation increment: Psi_n+1 - Psi_n

        if(isw.eq.13) then
          theta1 = 1.d0
          thetan = 0.d0
          rtheta = 0.d0
        else
          theta1 = theta(3)
          thetan = 1.d0 - theta1
          rtheta = 1.d0/theta1 - 1.d0
        endif
        thetan = 1.d0 - theta1
        dpsi   = 0.d0
        do ii = 1,nel
          dpsi = dpsi + ul(3,ii,2)*shp(2,ii)
        end do ! ii
        dpsi   = dpsi/theta(3)
        dpsih  = 0.5d0*dpsi

c       Compute entries into rotation array

        lamn(1,1) =  hn(1)
        lamn(2,1) =  hn(2)
        lamn(1,2) = -hn(2)
        lamn(2,2) =  hn(1)

        denr      =  1.d0/(1.d0 + dpsih*dpsih)
        cay(1,1)  =  (1.d0 - dpsih*dpsih)*denr
        cay(2,1)  =  dpsi*denr
        cay(1,2)  = -cay(2,1)
        cay(2,2)  =  cay(1,1)

        dcay(1,1) =  denr*cay(1,2)
        dcay(2,1) =  denr*cay(1,1)
        dcay(1,2) = -dcay(2,1)
        dcay(2,2) =  dcay(1,1)

c       Spatial form of update

        lam1(1,1) =  cay(1,1)*lamn(1,1) + cay(1,2)*lamn(2,1)
        lam1(2,1) =  cay(2,1)*lamn(1,1) + cay(2,2)*lamn(2,1)
        lam1(1,2) =  cay(1,1)*lamn(1,2) + cay(1,2)*lamn(2,2)
        lam1(2,2) =  cay(2,1)*lamn(1,2) + cay(2,2)*lamn(2,2)

        dlam(1,1) =  dcay(1,1)*lamn(1,1) + dcay(1,2)*lamn(2,1)
        dlam(2,1) =  dcay(2,1)*lamn(1,1) + dcay(2,2)*lamn(2,1)
        dlam(1,2) =  dcay(1,1)*lamn(1,2) + dcay(1,2)*lamn(2,2)
        dlam(2,2) =  dcay(2,1)*lamn(1,2) + dcay(2,2)*lamn(2,2)

        lama(1,1) =  thetan*lamn(1,1) + theta1*lam1(1,1)
        lama(2,1) =  thetan*lamn(2,1) + theta1*lam1(2,1)
        lama(1,2) =  thetan*lamn(1,2) + theta1*lam1(1,2)
        lama(2,2) =  thetan*lamn(2,2) + theta1*lam1(2,2)

c       Save terms in LAMBDA_n+1

        h1(1) =  lam1(1,1)
        h1(2) =  lam1(2,1)

c       Compute strains and forces in Gaussian frame

c       t_n   deformation gradient evaluation

        defn(1) = lamn(1,1)*(fa(1)-df(1))+lamn(2,1)*(fa(2)-df(2))-1.d0
        defn(2) = lamn(1,2)*(fa(1)-df(1)) + lamn(2,2)*(fa(2)-df(2))
        defn(3) = (fa(3) - df(3))

c       t_n+1 deformation gradient evaluation

        do i = 1,3
          df(i) = df(i)*rtheta
        end do ! i

        def1(1) = lam1(1,1)*(fa(1)+df(1))+lam1(2,1)*(fa(2)+df(2))-1.d0
        def1(2) = lam1(1,2)*(fa(1)+df(1))+lam1(2,2)*(fa(2)+df(2))
        def1(3) = (fa(3) + df(3))

c       Mid-point evaluation

        defa(1,1) = thetan*defn(1) + theta1*def1(1)
        defa(2,1) = thetan*defn(2) + theta1*def1(2)
        defa(3,1) = thetan*defn(3) + theta1*def1(3)

c       Constitution for stresses

        call bm2con (d,hn(3),h1(3),nh,cc,strs,defa,def1,isw)

c       Compute stored energy density

        if(isw.eq.13) then

          energy = strs(1,1)*def1(1)
     &           + strs(2,1)*def1(2)
     &           + strs(3,1)*def1(3)

c       Compute tangent tensor

        elseif(isw.ne.12) then

c         Compute first Piola-material frame

          forca(1) = lama(1,1)*strs(1,1) + lama(1,2)*strs(2,1)
          forca(2) = lama(2,1)*strs(1,1) + lama(2,2)*strs(2,1)
          forca(3) = strs(3,1)
          forca(4) = fa(2)*forca(1) - fa(1)*forca(2)

          forc1(1) = dlam(1,1)*strs(1,1) + dlam(1,2)*strs(2,1)
          forc1(2) = dlam(2,1)*strs(1,1) + dlam(2,2)*strs(2,1)

c         Compute tangent elastic tensor

          b(1,1) = lama(1,1)*cc(1,1,1) + lama(1,2)*cc(2,1,1)
          b(1,2) = lama(1,1)*cc(1,2,1) + lama(1,2)*cc(2,2,1)
          b(2,1) = lama(2,1)*cc(1,1,1) + lama(2,2)*cc(2,1,1)
          b(2,2) = lama(2,1)*cc(1,2,1) + lama(2,2)*cc(2,2,1)

          a(1,1) = b(1,1)*lam1(1,1)    + b(1,2)*lam1(1,2)
          a(1,2) = b(1,1)*lam1(2,1)    + b(1,2)*lam1(2,2)
          a(1,3) = lama(1,1)*cc(1,3,1) + lama(1,2)*cc(2,3,1)

          a(2,1) = b(2,1)*lam1(1,1)    + b(2,2)*lam1(1,2)
          a(2,2) = b(2,1)*lam1(2,1)    + b(2,2)*lam1(2,2)
          a(2,3) = lama(2,1)*cc(1,3,1) + lama(2,2)*cc(2,3,1)

          a(3,1) = cc(3,1,1)*lam1(1,1) + cc(3,2,1)*lam1(1,2)
          a(3,2) = cc(3,1,1)*lam1(2,1) + cc(3,2,1)*lam1(2,2)
          a(3,3) = cc(3,3,1)

c         Compute remainder

          phia1  = fa(1)
          phia2  = fa(2)
          phi11  = dlam(1,1)*(fa(1)+df(1)) + dlam(2,1)*(fa(2)+df(2))
          phi12  = dlam(1,2)*(fa(1)+df(1)) + dlam(2,2)*(fa(2)+df(2))

          a(1,4) = b(1,1)*phi11    + b(1,2)*phi12
          a(2,4) = b(2,1)*phi11    + b(2,2)*phi12
          a(3,4) = cc(3,1,1)*phi11 + cc(3,2,1)*phi12

          a(4,1) = phia2*a(1,1) - phia1*a(2,1)
          a(4,2) = phia2*a(1,2) - phia1*a(2,2)
          a(4,3) = phia2*a(1,3) - phia1*a(2,3)

          a(4,4) = phia2*a(1,4) - phia1*a(2,4)

c         Add geometric tangent part

          if(gflag) then
            a(1,4) = a(1,4) + forc1(1)
            a(2,4) = a(2,4) + forc1(2)

            a(4,1) = a(4,1) - forca(2)
            a(4,2) = a(4,2) + forca(1)

            a(4,4) = a(4,4) + phia2*forc1(1) - phia1*forc1(2)
          endif

        endif

c     Generalized mid-point formulation

      else

c       Compute matrix LAMBDA

        psi = 0.d0
        do ii = 1,nel
          psi = psi + ul(3,ii,1)*shp(2,ii)
        end do ! ii
        cn = cos(psi)
        sn = sin(psi)

c       Compute strains and forces in gaussian frame

        defa(1,1) = fa(1)*cn + fa(2)*sn - 1.d0
        defa(2,1) = fa(2)*cn - fa(1)*sn
        defa(3,1) = fa(3)

        call bm2con (d,hn(3),h1(3),nh,cc,strs,defa,def1,isw)

c       Compute stored energy density

        if(isw.eq.13) then

          energy = strs(1,1)*defa(1,1)
     &           + strs(2,1)*defa(2,1)
     &           + strs(3,1)*defa(3,1)

        elseif(isw.ne.12) then

c         Compute first Piola-material frame

          forca(1) = cn*strs(1,1) - sn*strs(2,1)
          forca(2) = sn*strs(1,1) + cn*strs(2,1)
          forca(3) = strs(3,1)
          forca(4) = fa(2)*forca(1) - fa(1)*forca(2)

c         Compute tangent elastic tensor

          b(1,1) = cn*cc(1,1,1) - sn*cc(2,1,1)
          b(1,2) = cn*cc(1,2,1) - sn*cc(2,2,1)
          b(2,1) = sn*cc(1,1,1) + cn*cc(2,1,1)
          b(2,2) = sn*cc(1,2,1) + cn*cc(2,2,1)

          a(1,1) = b(1,1)*cn    - b(1,2)*sn
          a(1,2) = b(1,1)*sn    + b(1,2)*cn
          a(1,3) = cn*cc(1,3,1) - sn*cc(2,3,1)

          a(2,1) = a(1,2)
          a(2,2) = b(2,1)*sn    + b(2,2)*cn
          a(2,3) = sn*cc(1,3,1) + cn*cc(2,3,1)

          a(3,1) = a(1,3)
          a(3,2) = a(2,3)
          a(3,3) = cc(3,3,1)

          phia1  = defa(1,1) + 1.d0
          phia2  = defa(2,1)

          a(1,4) = b(1,1)*phia2    - b(1,2)*phia1
          a(2,4) = b(2,1)*phia2    - b(2,2)*phia1
          a(3,4) = cc(1,3,1)*phia2 - cc(2,3,1)*phia1

          a(4,4) = phia2*(cc(1,1,1)*phia2 - cc(1,2,1)*phia1)
     &           + phia1*(cc(2,2,1)*phia1 - cc(1,2,1)*phia2)

          if(gflag) then
            a(1,4) = a(1,4) - forca(2)
            a(2,4) = a(2,4) + forca(1)
            a(4,4) = a(4,4) - fa(1)*forca(1) - fa(2)*forca(2)
          endif

          a(4,1) = a(1,4)
          a(4,2) = a(2,4)
          a(4,3) = a(3,4)

        endif

      endif

      end
