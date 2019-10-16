      subroutine mate3d01(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c
c     Purpose: calculate S and C for
c               a linear elastic isotropic material law
c
c     Inputs:
c         h1(nh)    - history array h1
c         d(md)     - local d-array
c         Eps       - strains
c         isw       - solution option from element
c
c     Input material parameter:
c         E         - Young's modulus
c         v         - Poisson's ratio
c
c     Outputs:
c         md = 2    - number of used data for control of d-array
c         nh = 0    - number of history parameter at Gauss-Point
c         sig       - stresses
c         cmat      - tangent modulus
c         plout(10) - plot data
c
c     Allocation of d-array:
c         d(1) = E  - Young's modulus
c         d(2) = v  - Poisson's ratio
c
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(3),eps(6),sig(6),cmat(6,6),plout(10)

      if(isw.eq.1) then
c....   input data
        call mati3d01(d,md,nh)
c
      else
c....   elastic tangent modul
        call matm3d01(cmat,d(1),d(2))
c
c....   stresses sig(i)=cmat(i,j)*eps(j)
        sig = matmul(cmat,eps)
c
c        call mprint(cmat,6,6,6,'Cmat')
c        call mprint( sig,6,1,6,' Sig')

      end if

      return
      end
c
      subroutine mati3d01(d,md,nh)
c-----------------------------------------------------------------------
c
c     mate3d01 Linear elastic isotropic
c     input material parameter
c
c-----------------------------------------------------------------------
c
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(2)

      md=2
      if(ior.lt.0) write(*,1001)
1001  format(
     + ' Input: E,nu')
      nh   = 0
      call dinput(d,md)
      if(d(2).le.-1.d0.or.d(2).ge.0.5d0)stop ' -1 < nue < 0.5  ! '
      g=0.5d0*d(1)/(1.d0+d(2))
                  write(iow,1002) nh,(d(i),i=1,md),g
      if(ior.lt.0)write(*  ,1002) nh,(d(i),i=1,md),g
1002  format(5x,'Linear elastic isotropic material data',/,
     +  5x,'length nh of h1,h2 ........',i12,/,
     +  5x,'elastic modulus E .........',g12.5,/,
     +  5x,'poissons ratio  v .........',f12.5,/,
     +  5x,'  shear modulus G .........',g12.5,/)

      return
      end
c
      subroutine matm3d01(cmat,e1,xnu)
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension cmat(6,6)
c.... linear isotropic elasticity matrix
c
      cc = e1/((1.0d0+xnu)*(1.0d0-2.0d0*xnu))
      tt = (1.0d0-2.0d0*xnu)/2.0d0
      cc1xnu = cc*(1.0d0-xnu)
      ccxnu = cc*xnu
      cctt  = cc*tt
c
      cmat(1,1) = cc1xnu
      cmat(1,2) = ccxnu
      cmat(1,3) = ccxnu
      cmat(2,2) = cc1xnu
      cmat(2,3) = ccxnu
      cmat(3,3) = cc1xnu
      cmat(4,4) = cctt
      cmat(5,5) = cctt
      cmat(6,6) = cctt
c
c.... symmetry
      cmat(2,1) = cmat(1,2)
      cmat(3,1) = cmat(1,3)
      cmat(3,2) = cmat(2,3)
c
      return
      end