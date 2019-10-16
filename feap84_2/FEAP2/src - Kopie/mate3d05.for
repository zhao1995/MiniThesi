      subroutine mate3d05(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c     Purpose: calculate S and IC for 
c              a finite strain elastoplastic material law 
c              F = Fe * Fp 
c
c     Inputs:
c         h1(nh)     - history array h1
c         nh = 7     - length of history arrays, to define in element
c         d(md)      - local d-array
c         Eps        - strains  
c         isw        - solution option from element 
c
c     Input material parameter:                
c         E          - Young's modulus
c         v          - Poisson's ratio
c         y_o        - initial yield stress
c         y_i        - yield stress at t=infty
c         xk         - linear hardening modul
c         xd         - exponential hardening
c
c     Outputs:
c         md = 6     - number of used data for control of d-array 
c         nh = 7     - number of history parameter at Gauss-Point
c         h2(nh)     - history array h2
c         Sig        - 2nd Piola-Kirchhoff Stress
c         Cmat       - algorithimic consistent tangent modulus
c         plout(10)  - plot data    
c
c     Allocation of d-array:
c         d(1) = E   - Young's modulus
c         d(2) = v   - Poisson's ratio
c         d(3) = y_o - initial yield stress
c         d(4) = y_i - yield stress at t=infty
c         d(5) = xk  - linear hardening modul
c         d(6) = xd  - exponential hardening
c         yield condition: f = g - y_o+ xk*a + (y_i-y_o)*(1-exp(-xd*a))
c                          g = sqrt(3/2 Sd:Sd) - sd deviatoric stress
c                          a = equivalent plastic strain
c
c     References:
c     [1] Simo,J.C.; Algorithm for static and dynamic multiplicative
c         plasticity that preserve the classical return mapping schemes
c         of infinitesimal theory;Comp.Meth.Appl.Mech.Eng.99(1992)61-112
c     [2] Klinkel,S.; Theorie und Numerik eines Volumen-Schalen-
c         Elementes bei finiten elastischen und plastischen Verzerrungen
c         Diss. Baustatik Uni KA Bericht-Nr.7(2000)
c
c     (c) s.klinkel                                             May,2001
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(*),Eps(*),Sig(*),Cmat(6,6),plout(10),
     +          C_p(6)
c
      if(isw.eq.1) then
c....   input data
        call mati3d05(d,md,nh)
c
      else
c....   Store history variables at time n (h1-array) in plastic variables
        do i = 1,6
           C_p(i) = h1(i)
        end do
        a = h1(7)
c
c       compute elastoplastic tangent modul and stresses
        call matm3d05(d,Eps,C_p,a,sig,Cmat,ngp)
c
c....   store history variables at time n+1 in h2-array
        do i = 1,6
           h2(i) = C_p(i)
        end do
        h2(7) = a
        plout(1) = a    
      end if
c
      return
      end
c
      subroutine mati3d05(d,md,nh)
c-----------------------------------------------------------------------
c
c     mate3d05 finite strain elastoplastic 
c     input material parameter 
c     
c-----------------------------------------------------------------------
c
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(6)

      md=6
      if(ior.lt.0) write(*,1001)
1001  format(
     + ' Input: E,nu,y_0,y_i,xk,xd')
      nh = 7
      call dinput(d,md)
                  write(iow,1002) nh,(d(i),i=1,md)
      if(ior.lt.0)write(*  ,1002) nh,(d(i),i=1,md)
1002  format(5x,'Finite strain elastoplastic material data',/,
     +  5x,'length nh of h1,h2 ...........',i12,/,
     +  5x,'Youngs modulus E .............',g12.5,/,
     +  5x,'Poissons ratio v .............',g12.4,/,
     +  5x,'Initial yield stress    y_0 ..',g12.5,/,
     +  5x,'Yield stress at t=infty y_i ..',g12.5,/,
     +  5x,'Linear hardening modul xk ....',g12.5,/,
     +  5x,'Exponential hardening  xd ....',g12.5,/)
c
      return
      end 
c
      subroutine matm3d05(d,Eps,C_p,a,S,Cmat,ngp)
c-----------------------------------------------------------------------
c.... stresses and tangent matrix for plane stress J-2 plasticity with
c     linear isotropic hardening  (conv. coordinates!)
c
c.... Input parameters     in h1
c      epn(5)    - plastic strains           at t_n
c      alpha     - effective plastic strain  at t_n
c
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      double precision L1(3,3),L2(3,3)
      dimension d(*),S(6),cmat(6,6)
      dimension G(3,3),
     +   C(3,3),C_n(3,3),C_p(6),Cp(3,3),Cp_n(3,3),
     +   xlam(3),xN(3,3),fv1(3),fv2(3),
     +   ee(3),ep(3),e_tr(3),
     +   Ce(3,3),T_tr(3),Td_tr(3),sig(3),Eps(*),
     +   T1(6,3),T2(6,3),T1L1(6,3),T2L2(6,3),T1L1T1(6,6),T2L2T2(6,6)

      data G / 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0 /
c
      xkappa  = d(1)/(3.d0*(1.d0-2.d0*d(2)))
      xmu  = d(1)/(2.d0*(1.d0+d(2)))
      y_o  = d(3)
      y_i  = d(4)
      xk   = d(5)
      xd   = d(6)
c
c.....Right Cauchy Green tensor
      do i=1,3
         C(i,i) = 2.d0*Eps(i) + 1.d0
      end do
      C(1,2) = Eps(4) 
      C(1,3) = Eps(5) 
      C(2,3) = Eps(6) 
      C(2,1) = Eps(4) 
      C(3,1) = Eps(5) 
      C(3,2) = Eps(6) 

c     Plastic Metric tensor Cp
      do i = 1,3
         Cp(i,i) = C_p(i)
      end do
      Cp(1,2) = C_p(4)
      Cp(1,3) = C_p(5)
      Cp(2,3) = C_p(6)
      Cp(2,1) = C_p(4)
      Cp(3,1) = C_p(5)
      Cp(3,2) = C_p(6)
      detCp =  Cp(1,1)*(Cp(2,2)*Cp(3,3)-Cp(2,3)*Cp(3,2))
     +        -Cp(1,2)*(Cp(2,1)*Cp(3,3)-Cp(2,3)*Cp(3,1))
     +        +Cp(1,3)*(Cp(2,1)*Cp(3,2)-Cp(2,2)*Cp(3,1))
      detC  =  C(1,1)* (C(2,2)*C(3,3)-C(2,3)*C(3,2))
     +        -C(1,2)* (C(2,1)*C(3,3)-C(2,3)*C(3,1))
     +        +C(1,3)* (C(2,1)*C(3,2)-C(2,2)*C(3,1))
      if (detC.lt.1.e-14) then
         write(*,*) 'In Elas det C <= 0 ', detC
         write(iow,*) 'In Elas det C <= 0 '
c         write(*,*) 'S(3)=', S(3)
c         write(*,*) 'e(i)=',(e(i),i=1,6)
c         write(*,*) 'C(i,j)=',(C(1,i),i=1,3)
c         write(*,*) 'C(i,j)=',(C(2,i),i=1,3)
c         write(*,*) 'C(i,j)=',(C(3,i),i=1,3)
         stop
      end if
      do i=1,3
         do j=1,3
            if (abs(detCp).lt.1.e-14 ) then
               Cp(i,j) = G(i,j)
            end if
            if(dabs(C(i,j)).lt.1.e-13) C(i,j) = 0.d0
            if(dabs(Cp(i,j)).lt.1.e-13) Cp(i,j) = 0.d0
            Cp_n(i,j) = Cp(i,j)
            C_n(i,j) = C(i,j)
         end do
      end do
c
c.....Eigenvalue problem [C - lam_tr Cp] N = 0
      call pzero(xlam,3)
      call pzero(xN,3*3)
      call pzero(fv1,3)
      call pzero(fv2,3)
      ierr = 0
      call rsg(3,3,C,Cp,xlam,3,xN,fv1,fv2,ierr)
      do i=1,3
         if (xlam(i).lt.1.e-14) then
            write(*,*) 'xlam(',i,')=',xlam(i)
            stop
         end if
      end do
      if (abs(xlam(1)-xlam(2)).lt.1.e-15) xlam(2)=xlam(1)+1.e-15
      if (abs(xlam(1)-xlam(3)).lt.1.e-15) xlam(3)=xlam(1)+1.e-15
      if (abs(xlam(2)-xlam(3)).lt.1.e-15) xlam(3)=xlam(2)+1.e-15
c
c.....test of the eigenvectors
c     xN * (Cp xN) = 1
      do k=1,3
         xx = 0.d0
         do i=1,3
            xx = xx + xN(i,k)*
     +          (Cp_n(i,1)*xN(1,k)+Cp_n(i,2)*xN(2,k)+Cp_n(i,3)*xN(3,k))
         end do
         if (dabs(xx)-1.d0.gt.1.e-14) then
            write(*,*) 'N(',k,') Cp N(',k,') >1'
            stop
         end if
      end do
c
c.....logarithmic strains
      do i=1,3
         e_tr(i) = 0.5d0 * dlog(xlam(i))
      end do
c
c.....trial values
      do i=1,3
         do j=1,3
            fac=2.d0
            if(i.ne.j) fac=-1.d0
            Ce(i,j) = xkappa + 2.d0 * xmu * fac/3.d0
         end do
      end do
      do i=1,3
         T_tr(i) = 0.d0
         do j=1,3
            T_tr(i) = T_tr(i) + Ce(i,j) * e_tr(j)
         end do
      end do
      Td_tr(1)= 2.d0/3.d0*T_tr(1)-1.d0/3.d0*T_tr(2)-1.d0/3.d0*T_tr(3)
      Td_tr(2)=-1.d0/3.d0*T_tr(1)+2.d0/3.d0*T_tr(2)-1.d0/3.d0*T_tr(3)
      Td_tr(3)=-1.d0/3.d0*T_tr(1)-1.d0/3.d0*T_tr(2)+2.d0/3.d0*T_tr(3)
      g_tr = dsqrt(3.d0/2.d0*
     +       (Td_tr(1)*Td_tr(1)+Td_tr(2)*Td_tr(2)+Td_tr(3)*Td_tr(3)))
      y = y_o + xk*a + (y_i-y_o)*(1-exp(-xd*a))
      f = g_tr - y
c
c.....plastic corrector
      if (y_o.gt.1e-10.and.f.gt.1.d-10*y) then
         call plas3d05(d,T_tr,Td_tr,g_tr,f,Ce,ep,a,ngp)
c.....   compute Cp_n+1
         do i=1,3
            ee(i) = e_tr(i) - ep(i)
         end do
         call pzero(Cp,3*3)
         do k=1,3
            fac = dexp(2.d0*ee(k)) / dexp(2.d0*e_tr(k))
            do i=1,3
               do j=1,3
               Cp(i,j) = Cp(i,j) + fac * xN(i,k)*xN(j,k)
               end do
            end do
         end do
         call invert(Cp,3,3)
      end if
      C_p(1) = Cp(1,1)
      C_p(2) = Cp(2,2)
      C_p(3) = Cp(3,3)
      C_p(4) = Cp(1,2)
      C_p(5) = Cp(1,3)
      C_p(6) = Cp(2,3)
c
c.....transformation of the stresses
      do i=1,3
         sig(i) = T_tr(i)/xlam(i)
         T1(1,i) = xN(1,i)*xN(1,i)
         T1(2,i) = xN(2,i)*xN(2,i)
         T1(3,i) = xN(3,i)*xN(3,i)
         T1(4,i) = xN(1,i)*xN(2,i)
         T1(5,i) = xN(1,i)*xN(3,i)
         T1(6,i) = xN(2,i)*xN(3,i)
      end do
      do i = 1,6
         S(i) = 0.d0
         do j=1,3
            S(i) = S(i) + T1(i,j) * sig(j)
         end do
      end do
c
c.....transformation of the elastoplastic tangent
      T2(1,1) = xN(1,1)*xN(1,2) + xN(1,2)*xN(1,1)
      T2(2,1) = xN(2,1)*xN(2,2) + xN(2,2)*xN(2,1)
      T2(3,1) = xN(3,1)*xN(3,2) + xN(3,2)*xN(3,1)
      T2(4,1) = xN(1,1)*xN(2,2) + xN(1,2)*xN(2,1)
      T2(5,1) = xN(1,1)*xN(3,2) + xN(1,2)*xN(3,1)
      T2(6,1) = xN(2,1)*xN(3,2) + xN(2,2)*xN(3,1)
      T2(1,2) = xN(1,1)*xN(1,3) + xN(1,3)*xN(1,1)
      T2(2,2) = xN(2,1)*xN(2,3) + xN(2,3)*xN(2,1)
      T2(3,2) = xN(3,1)*xN(3,3) + xN(3,3)*xN(3,1)
      T2(4,2) = xN(1,1)*xN(2,3) + xN(1,3)*xN(2,1)
      T2(5,2) = xN(1,1)*xN(3,3) + xN(1,3)*xN(3,1)
      T2(6,2) = xN(2,1)*xN(3,3) + xN(2,3)*xN(3,1)
      T2(1,3) = xN(1,2)*xN(1,3) + xN(1,3)*xN(1,2)
      T2(2,3) = xN(2,2)*xN(2,3) + xN(2,3)*xN(2,2)
      T2(3,3) = xN(3,2)*xN(3,3) + xN(3,3)*xN(3,2)
      T2(4,3) = xN(1,2)*xN(2,3) + xN(1,3)*xN(2,2)
      T2(5,3) = xN(1,2)*xN(3,3) + xN(1,3)*xN(3,2)
      T2(6,3) = xN(2,2)*xN(3,3) + xN(2,3)*xN(3,2)
c
      do i=1,3
         L1(i,i) = ( Ce(i,i) - 2.d0*T_tr(i) ) / (xlam(i)*xlam(i))
      end do
      L1(1,2) = Ce(1,2) / (xlam(1)*xlam(2))
      L1(1,3) = Ce(1,3) / (xlam(1)*xlam(3))
      L1(2,3) = Ce(2,3) / (xlam(2)*xlam(3))
      L1(2,1) = L1(1,2)
      L1(3,1) = L1(1,3)
      L1(3,2) = L1(2,3)
c
      call pzero(L2,3*3)
      L2(1,1) = (Sig(1)-Sig(2))/(xlam(1)-xlam(2))
      L2(2,2) = (Sig(1)-Sig(3))/(xlam(1)-xlam(3))
      L2(3,3) = (Sig(2)-Sig(3))/(xlam(2)-xlam(3))
c
      do i=1,6
         do j=1,3
            T1L1(i,j) = 0.d0
            T2L2(i,j) = 0.d0
            do k=1,3
              T1L1(i,j) = T1L1(i,j) + T1(i,k)*L1(k,j)
              T2L2(i,j) = T2L2(i,j) + T2(i,k)*L2(k,j)
            end do
         end do
      end do
      do i=1,6
         do j=1,6
            T1L1T1(i,j) = 0.d0
            T2L2T2(i,j) = 0.d0
            do k=1,3
               T1L1T1(i,j) = T1L1T1(i,j) + T1L1(i,k)*T1(j,k)
               T2L2T2(i,j) = T2L2T2(i,j) + T2L2(i,k)*T2(j,k)
            end do
            cmat(i,j) = T1L1T1(i,j) + T2L2T2(i,j)
         end do
      end do
c
      return
      end
c
      subroutine plas3d05(d,T_tr,Td_tr,g_tr,f,Ce,ep,a,ngp)
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(*),T_tr(3),Td_tr(3),Ce(3,3),ep(3)
c
      xkappa  = d(1)/(3.d0*(1.d0-2.d0*d(2)))
      xmu  = d(1)/(2.d0*(1.d0+d(2)))
      y_o  = d(3)
      y_i  = d(4)
      xh   = d(5)
      xd   = d(6)
c.....Startvalues
      gam = 0.d0
      newton = 0
c
c.....begin local Newton iteration.......................
100   continue
      newton = newton + 1
c      if (n.eq.1) write(iow,*) 'tau=',T_tr(1),T_tr(2),T_tr(3)
      if (newton.gt.99) then
         write( * ,*) 'no converg. in local newton'
         write( * ,*) 'f,n',f,ngp
         write(iow,*) 'no converg. in local newton'
         stop
      endif
      dgam = f/(3.d0*xmu + xh + (y_i-y_o)*xd*exp(-xd*a))
      gam = gam + dgam
      a = a + dgam
      y = y_o + xh*a + (y_i-y_o)*(1.d0 - exp(-xd*a))
      f = g_tr - 3.d0*xmu*gam - y
      if (dabs(f).gt.1.d-10*y) go to 100
      gam = gam - 1.e-9
c.....end local Newton iteration.........................
c
c.....update stresses and plastic strains
      g = g_tr-3.d0*xmu*gam
      trT = (T_tr(1)+T_tr(2)+T_tr(3))/3.d0
      do i=1,3
         T_tr(i) = g/g_tr * Td_tr(i) + trT
         ep(i) = 3.d0/2.d0*gam/g_tr*Td_tr(i)
      end do
c
c.....elasto-plastic tangent operator
      dyda = xh + (y_i-y_o)*xd*exp(-xd*a)
      xnt = Td_tr(1)*Td_tr(1)+Td_tr(2)*Td_tr(2)+Td_tr(3)*Td_tr(3)
csk      beta1 = 2.d0*xmu*g/g_tr
      beta1 = 2.d0*xmu /(1.d0 + 3.d0 * gam * xmu / y)
      beta2 =((1.d0/(dyda/(3.d0*xmu)+1.d0) + g/g_tr - 1.d0)*2.d0*xmu)/
     +       xnt
      Ce(1,1) = xkappa + beta1*2.d0/3.d0
      Ce(1,2) = xkappa - beta1*1.d0/3.d0
      Ce(1,3) = xkappa - beta1*1.d0/3.d0
      Ce(2,1) = xkappa - beta1*1.d0/3.d0
      Ce(2,2) = xkappa + beta1*2.d0/3.d0
      Ce(2,3) = xkappa - beta1*1.d0/3.d0
      Ce(3,1) = xkappa - beta1*1.d0/3.d0
      Ce(3,2) = xkappa - beta1*1.d0/3.d0
      Ce(3,3) = xkappa + beta1*2.d0/3.d0
      do i=1,3
         do j=1,3
            Ce(i,j) = Ce(i,j) - beta2*Td_tr(i)*Td_tr(j)
         end do
      end do
c
      return
      end
