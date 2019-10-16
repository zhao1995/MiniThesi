      subroutine mate3d06(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c     Purpose: calculate Sig and Cmat for 
c              a finite strain elastic material law of Ogden-type
c
c     Inputs:
c         h1(nh)  - history array h1
c         d(md)   - local d-array
c         Eps     - strains  
c         isw     - solution option from element 
c
c     Input material parameter:                
c         mue_1   - shear modulus
c         mue_2   - shear modulus
c         mue_3   - shear modulus
c         alpha_1 - sum(mue_i*alpha_i) = 2mue
c         alpha_2
c         alpha_3
c         lame    - Lamé const.
c
c     Outputs:
c         md = 7  - number of used data for control of d-array 
c         nh = 0  - number of history parameter at Gauss-Point
c         Sig     - 2nd Piola-Kirchhoff Stress
c         Cmat    - algorithimic consistent tangent modulus
c         plout(10)- plot data    
c
c     Allocation of d-array:
c         d(1) = mue_1   - shear modulus
c         d(2) = mue_2   - shear modulus
c         d(3) = mue_3   - shear modulus
c         d(4) = alpha_1 - sum(mue_i*alpha_i) = 2mue
c         d(5) = alpha_2
c         d(6) = alpha_3
c         d(7) = lame    - Lamé const.
c
c     References:
c     [1] Ogden,R.W.; Large deformation isotropic Elasticity - On the 
c         Correlation of Theory and Experiment for compressible 
c         rubberlike solids; Proc. Royal Soc. London A328(1972)565-584
c     [2] Klinkel,S.; Theorie und Numerik eines Volumen-Schalen-
c         Elementes bei finiten elastischen und plastischen Verzerrungen
c         Diss. Baustatik Uni KA Bericht-Nr.7(2000)
c
c     (c) s.klinkel                                             May,2001
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(*),Eps(6),Sig(6),Cmat(6,6),plout(10)
c
      if(isw.eq.1) then
c....   input data
        call mati3d06(d,md,nh)
c
      else
c....   compute elastic tangent moduli and stresses
        call matm3d06(d,Eps,sig,Cmat)
      end if
c
      return
      end
c
      subroutine mati3d06(d,md,nh)
c-----------------------------------------------------------------------
c
c     mate05 finite elastic strain Ogden material 
c     input material parameter 
c     
c-----------------------------------------------------------------------
c
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(7)

      md=7
      if(ior.lt.0) write(*,1001)
1001  format(
     + ' Input:mue1,mue2,mue3,alpha1,alpha2,alpha3,lame')
      nh = 0
      call dinput(d,md)
                  write(iow,1002) nh,(d(i),i=1,md)
      if(ior.lt.0)write(*  ,1002) nh,(d(i),i=1,md)
1002  format(5x,'Finite elastic strain Ogden material data',/,
     +  5x,'length nh of h1,h2 ...........',i12,/,
     +  5x,'Shear moduli mue_1 ...........',g12.5,/,
     +  5x,'............ mue_2 ...........',g12.5,/,
     +  5x,'............ mue_3 ...........',g12.5,/,
     +  5x,'Exponents    alpha_1 .........',g12.5,/,
     +  5x,'............ alpha_2 .........',g12.5,/,
     +  5x,'............ alpha_3 .........',g12.5,/,
     +  5x,'Lamé constant ................',g12.5,/)
c
      return
      end 
c
      subroutine matm3d06(d,Eps,S,Cmat)
c-----------------------------------------------------------------------
c
c     mate3d06 finite elastic strains, Ogden model 
c
c-----------------------------------------------------------------------
      USE eldata
      USE iofile
      implicit double precision (a-h,o-z)
      double precision L1(3,3),L2(3,3)
      dimension d(*),S(6),cmat(6,6)
      dimension G(3,3),C(3,3),
     +   xmue(3),alpha(3),
     +   xlam(3),xN(3,3),fv1(3),fv2(3),
     +   Eps(*),sig(3),
     +   T1(6,3),T2(6,3),T1L1(6,3),T2L2(6,3),T1L1T1(6,6),T2L2T2(6,6)

      data G / 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0 /
c    
c
      xmue(1)  = d(1)
      xmue(2)  = d(2)
      xmue(3)  = d(3) 
      alpha(1) = d(4)  
      alpha(2) = d(5)  
      alpha(3) = d(6)
      xlame    = d(7)    

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
c
      detC  =  C(1,1)* (C(2,2)*C(3,3)-C(2,3)*C(3,2))
     +        -C(1,2)* (C(2,1)*C(3,3)-C(2,3)*C(3,1))  
     +        +C(1,3)* (C(2,1)*C(3,2)-C(2,2)*C(3,1))  
      if (detC.lt.1.e-14) then 
         write(*,*) 'In matm3d06 det C <= 0 '
         write(iow,*) 'In matm3d06 det C <= 0 '
      end if
      do i=1,3
         do j=1,3
            if(dabs(C(i,j)).lt.1.e-13) C(i,j) = 0.d0
         end do
      end do
c
c.....Eigenvalue problem [C - lam_tr G] N = 0
      call pzero(xlam,3)
      call pzero(xN,3*3)
      call pzero(fv1,3)
      call pzero(fv2,3)
      ierr = 0
      call rsg(3,3,C,G,xlam,3,xN,fv1,fv2,ierr)
cfg      do i=1,3
cfg         if (xlam(i).lt.1.e-14) then
cfg            write(*,*) 'xlam(',i,')=',xlam(i)
cfg         end if
cfg      end do
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
     +         (G(i,1)*xN(1,k)+G(i,2)*xN(2,k)+G(i,3)*xN(3,k))
         end do
         if (dabs(xx)-1.d0.gt.1.e-14) then 
            write(*,*) 'N(',k,') G N(',k,') <>1'
            stop
         end if
      end do
c
c.....Spannungsberechnung
      xJ       = xlam(1)*xlam(2)*xlam(3)  
      do i=1,3
         sig(i)=0.d0
         do k=1,3
          sig(i) = sig(i)+
     +     1.d0/xlam(i)*(xmue(k)*(xlam(i)**(alpha(k)/(2.d0))-1.d0))
         end do
       sig(i) = sig(i) + xlame/(2.d0*xlam(i))*(xJ-1.d0)
      end do
c
c.....transformation of the stresses      
      do i=1,3
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
c.....transformation of the elastic tangent
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
         L1(i,i) = 0.d0
         do k=1,3
           L1(i,i) = L1(i,i)+
     +               xmue(k)/(xlam(i)**2.d0)*
     +               ((alpha(k)-2.d0)*xlam(i)**(alpha(k)/2.d0)+2.d0) 
         end do
         L1(i,i) = L1(i,i) + xlame/(xlam(i)**2.d0)  
      end do
      L1(1,2) = xlame*xlam(3)
      L1(1,3) = xlame*xlam(2)
      L1(2,3) = xlame*xlam(1)
      L1(2,1) = xlame*xlam(3)
      L1(3,1) = xlame*xlam(2) 
      L1(3,2) = xlame*xlam(1)
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
