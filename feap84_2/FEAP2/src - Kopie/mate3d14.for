      subroutine mate3d14(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +            xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c     Purpose: calculate sig and cmat for 
c              a finite strain elastic material law of Blatz-Ko
c
c     Inputs:
c         h1(nh)  - history array h1
c         d(md)   - local d-array
c         eps     - strains  
c         isw     - solution option from element 
c
c     Input material parameter:                
c         xmu     - Shear modulus
c         xnu     - Poisson's ratio
c         xf      - Interpolation parameter
c        
c
c     Outputs:
c         md = 3  - number of used data for control of d-array 
c         nh = 0  - number of history parameter at Gauss-Point
c         sig     - Stress
c         cmat    - Tangent modulus
c         plout(10)- plot data    
c
c     Allocation of d-array:
c         d(1) = xmu
c         d(2) = xnu
c         d(3) = xf
c
c     References:
c     [1] Holzapfel,G.A.; Nonlinear solid mechanics - A continuum approach 
c         for engineering; John Wiley & Sons Ltd(2004)
c     [2] Klinkel,S.; Theorie und Numerik eines Volumen-Schalen-
c         Elementes bei finiten elastischen und plastischen Verzerrungen
c         Diss. Baustatik Uni KA Bericht-Nr.7(2000)
c
c     (c) m.romero, s.klinkel                                   Aug,2011
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)      
      dimension d(3),Eps(6),Sig(6),Cmat(6,6),plout(10)
      dimension C(3,3),Cinv(3,3),xgp(3)
c
      if(isw.eq.1) then
c....   input data
        call mati3d14(d,md,nh)
c
      else
c....   compute elastic tangent modul and stresses
        call elas3d14(d,Eps, C,Cinv,Sig,Cmat)
      end if
c
      return
      end
c
      subroutine mati3d14(d,md,nh)
c-----------------------------------------------------------------------
c
c     mate14 finite elastic strain blatz-Ko material 
c     input material parameter 
c     
c-----------------------------------------------------------------------
c
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(3)
      md=3
      if(ior.lt.0) write(*,1001)
1001  format(
     + ' Input:xmu_o,xnu_o,f')
      nh = 0
      call dinput(d,md)
                  write(iow,1002) nh,(d(i),i=1,md)
      if(ior.lt.0)write(*  ,1002) nh,(d(i),i=1,md)
1002  format(5x,'Finite elastic strain Blatz-Ko material data',/,
     +  5x,'length nh of h1,h2 ............',i12,/,
     +  5x,'Shear modulus xmu .............',g12.5,/,
     +  5x,'Poisson ratio xnu .............',g12.5,/,
     +  5x,'Interpolation parameter xf ....',g12.5,/)
c
      return
      end 
c
      subroutine elas3d14(d,eps, C,Cinv,Sig,Cmat)
c-----------------------------------------------------------------------
c
c     mate3d14 finite elastic strains, Blatz-Ko 
c     calculate C, sig and cmat
c
c-----------------------------------------------------------------------    
      USE eldata
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(3),Sig(6),Cmat(6,6),eps(6)
      dimension x1o1(6,6),CinvoCinv(6,6),CinvoC(6,6),CoCinv(6,6)
      dimension Cinvob2(6,6),dCinvdC(6,6)
      dimension b2oCinv(6,6),S(3,3),xII(6,6)
      dimension C(3,3),Cinv(3,3),b2(3,3)
c    
c
c.....Material Parameter
      xmu = d(1)
      xnu = d(2)
      xf  = d(3)
c.....Cauchy-Green Strain tensor
      C=0.0d0
      C(1,1)= (2.0d0*eps(1))+1.0d0
      C(2,2)= (2.0d0*eps(2))+1.0d0
      C(3,3)= (2.0d0*eps(3))+1.0d0
      C(1,2) = eps(4)
      C(1,3) = eps(5)
      C(2,3) = eps(6)
      C(2,1) = eps(4)
      C(3,1) = eps(5)
      C(3,2) = eps(6)
c
c.....compute Invariants      
      xI1 = C(1,1) + C(2,2) + C(3,3)
      xI2 = C(1,1)*C(2,2) + C(1,1)*C(3,3) + C(2,2)*C(3,3)
     +     -C(1,2)*C(1,2) - C(1,3)*C(1,3) - C(2,3)*C(2,3)
      xI3 = C(1,1)* (C(2,2)*C(3,3)-C(2,3)*C(3,2))
     +     -C(1,2)* (C(2,1)*C(3,3)-C(2,3)*C(3,1))  
     +     +C(1,3)* (C(2,1)*C(3,2)-C(2,2)*C(3,1))  
c
c.....compute coefficients of stresses
c     S = a0 I + a1 C + a2 C^-1
      chi =  2.0d0*xmu*(1.0d0-xf)/xI3
      a0 =  (xmu*xf)+(chi*xI1/2.0d0)
      a1 =  -chi/2
      beta = xnu/(1.0d0-(2.0d0*xnu))
      a2 = (-xmu*xf*(xI3**(-beta)))+(xmu*(1-xf)*xI3**beta)
     +     -((chi/2)*xI2)
c
c.....identity matrix 3x3
      b2=0.0d0
      do i=1,3
         b2(i,i)=1.0d0
      enddo
c.....compute Cinverse
      Cinv=0.0d0
      do i=1,3
         do j=1,3
            if (dabs(C(i,j)).lt.1.d-14) C(i,j)=0.0d0
            Cinv(i,j) = C(i,j)
         enddo
      enddo 
      call pivot(C,3,3,Cinv)
c.....compute stresses in matrix form
      S=0.0d0
      S = a0*b2 + a1*C + a2*Cinv
c.....compute stresses in vector notation
      Sig=0.0d0
      do i=1,3
         Sig(i) =  S(i,i)
      enddo
      Sig(4) = S(1,2)
      Sig(5) = S(1,3)
      Sig(6) = S(2,3)
c
c.....compute coefficients of elastic tangent
c     Cmat = g0 (IoI)+ g1 (dC^-1/dC)  + g2 (C^-1oC^-1) + g3 (C^-1oC+CoC^-1)
c          + g4 (C^-1oI+IoC^-1) + g5 (xII)
      g0 =  chi
      g1 = -2.0d0*xmu*xf*(xI3**(-beta)) 
     +     +2.0d0*xmu*(1.0d0-xf)*(xI3**(beta))
     +     -2.0d0*xmu*(1.0d0-xf)*(xI2/xI3)
      g2 = +2.0d0*xmu*xf*beta*(xI3**(-beta)) 
     +     +2.0d0*xmu*beta*(1.0d0-xf)*(xI3**beta)
     +     +chi*xI2
      g3 =  chi
      g4 = -chi*xI1
      g5 = -chi
c
c     compute (1o1)
      x1o1=0.0d0
      call AotimesB(b2,b2,x1o1)
c
c     compute dC^-1/dC  
c...  derivative dCi_dC(ijkl) = -Ci(ik)*Ci(jl)
      dCinvdC=0.0d0
      do i=1,3
        do j=1,3
           dCinvdC(i,j) = -Cinv(i,j)*Cinv(i,j)
        enddo
      enddo
      dCinvdC(1,4) = -(Cinv(1,1)*Cinv(1,2) + Cinv(1,2)*Cinv(1,1))/2.d0
      dCinvdC(1,5) = -(Cinv(1,1)*Cinv(1,3) + Cinv(1,3)*Cinv(1,1))/2.d0
      dCinvdC(1,6) = -(Cinv(1,2)*Cinv(1,3) + Cinv(1,3)*Cinv(1,2))/2.d0
      dCinvdC(2,4) = -(Cinv(2,1)*Cinv(2,2) + Cinv(2,2)*Cinv(2,1))/2.d0
      dCinvdC(2,5) = -(Cinv(2,1)*Cinv(2,3) + Cinv(2,3)*Cinv(2,1))/2.d0
      dCinvdC(2,6) = -(Cinv(2,2)*Cinv(2,3) + Cinv(2,3)*Cinv(2,2))/2.d0
      dCinvdC(3,4) = -(Cinv(3,1)*Cinv(3,2) + Cinv(3,2)*Cinv(3,1))/2.d0
      dCinvdC(3,5) = -(Cinv(3,1)*Cinv(3,3) + Cinv(3,3)*Cinv(3,1))/2.d0
      dCinvdC(3,6) = -(Cinv(3,2)*Cinv(3,3) + Cinv(3,3)*Cinv(3,2))/2.d0
c
      dCinvdC(4,1) = -(Cinv(1,1)*Cinv(2,1) + Cinv(2,1)*Cinv(1,1))/2.d0
      dCinvdC(4,2) = -(Cinv(1,2)*Cinv(2,2) + Cinv(2,2)*Cinv(1,2))/2.d0
      dCinvdC(4,3) = -(Cinv(1,3)*Cinv(2,3) + Cinv(2,3)*Cinv(1,3))/2.d0
      dCinvdC(5,1) = -(Cinv(1,1)*Cinv(3,1) + Cinv(3,1)*Cinv(1,1))/2.d0
      dCinvdC(5,2) = -(Cinv(1,2)*Cinv(3,2) + Cinv(3,2)*Cinv(1,2))/2.d0
      dCinvdC(5,3) = -(Cinv(1,3)*Cinv(3,3) + Cinv(3,3)*Cinv(1,3))/2.d0
      dCinvdC(6,1) = -(Cinv(2,1)*Cinv(3,1) + Cinv(3,1)*Cinv(2,1))/2.d0
      dCinvdC(6,2) = -(Cinv(2,2)*Cinv(3,2) + Cinv(3,2)*Cinv(2,2))/2.d0
      dCinvdC(6,3) = -(Cinv(2,3)*Cinv(3,3) + Cinv(3,3)*Cinv(2,3))/2.d0
c
      dCinvdC(4,4) = -(Cinv(1,1)*Cinv(2,2) + Cinv(1,2)*Cinv(2,1)
     +               + Cinv(2,1)*Cinv(1,2) + Cinv(2,2)*Cinv(1,1))/4.d0
      dCinvdC(4,5) = -(Cinv(1,1)*Cinv(2,3) + Cinv(1,3)*Cinv(2,1)
     +               + Cinv(2,1)*Cinv(1,3) + Cinv(2,3)*Cinv(1,1))/4.d0
      dCinvdC(4,6) = -(Cinv(1,2)*Cinv(2,3) + Cinv(1,3)*Cinv(2,2)
     +               + Cinv(2,2)*Cinv(1,3) + Cinv(2,3)*Cinv(1,2))/4.d0
      dCinvdC(5,4) = -(Cinv(1,1)*Cinv(3,2) + Cinv(1,2)*Cinv(3,1)
     +               + Cinv(3,1)*Cinv(1,2) + Cinv(3,2)*Cinv(1,1))/4.d0
      dCinvdC(5,5) = -(Cinv(1,1)*Cinv(3,3) + Cinv(1,3)*Cinv(3,1)
     +               + Cinv(3,1)*Cinv(1,3) + Cinv(3,3)*Cinv(1,1))/4.d0
      dCinvdC(5,6) = -(Cinv(1,2)*Cinv(3,3) + Cinv(1,3)*Cinv(3,2)
     +               + Cinv(3,2)*Cinv(1,3) + Cinv(3,3)*Cinv(1,2))/4.d0
      dCinvdC(6,4) = -(Cinv(2,1)*Cinv(3,2) + Cinv(2,2)*Cinv(3,1)
     +               + Cinv(3,1)*Cinv(2,2) + Cinv(3,2)*Cinv(2,1))/4.d0
      dCinvdC(6,5) = -(Cinv(2,1)*Cinv(3,3) + Cinv(2,3)*Cinv(3,1)
     +               + Cinv(3,1)*Cinv(2,3) + Cinv(3,3)*Cinv(2,1))/4.d0
      dCinvdC(6,6) = -(Cinv(2,2)*Cinv(3,3) + Cinv(2,3)*Cinv(3,2)
     +               + Cinv(3,2)*Cinv(2,3) + Cinv(3,3)*Cinv(2,2))/4.d0
c
c
c     compute (C^-1 o C^-1)
      CinvoCinv=0.0d0
      call AotimesB(Cinv,Cinv,CinvoCinv)
c
c     compute (C^-1 o C + C o C^-1)
      CinvoC=0.0d0
      CoCinv=0.0d0
      call AotimesB(Cinv,C,CinvoC)
      call AotimesB(C,Cinv,CoCinv)
c      
c     compute (I o C^-1+ C^-1 o I)
c     btwo(i) = b2(i,j)
      Cinvob2=0.0d0
      b2oCinv=0.0d0
      call AotimesB(Cinv,b2,Cinvob2)
      call AotimesB(b2,Cinv,b2oCinv)
c
c.....compute xII
      xII=0.0d0
      xII(1,1)=1.0d0
      xII(2,2)=1.0d0
      xII(3,3)=1.0d0
      xII(4,4)=0.5d0
      xII(5,5)=0.5d0
      xII(6,6)=0.5d0
c
c     compute Cmat
c     isotropic elastic tangent 
      Cmat=0.0d0 
      Cmat = g0*(x1o1) 
     +     + g1*(dCinvdC)
     +     + g2*(CinvoCinv) 
     +     + g3*(CoCinv + CinvoC)
     +     + g4*(b2oCinv + Cinvob2) 
     +     + g5*(xII)  
c
      return
      end
c
c
      subroutine AotimesB(A,B,AoB)
c-----------------------------------------------------------------------
c
c     AoB(ijkl) = A(ij)*B(kl)
c
c-----------------------------------------------------------------------
       implicit double precision (a-h,o-z)
       dimension A(3,3),B(3,3),AoB(6,6)
c      
       do i=1,3
         do j=1,3
         AoB(i,j) = A(i,i)*B(j,j)
        enddo
       enddo
       AoB(1,4) = (A(1,1)*B(1,2) + A(1,1)*B(2,1))/2.d0
       AoB(1,5) = (A(1,1)*B(1,3) + A(1,1)*B(3,1))/2.d0
       AoB(1,6) = (A(1,1)*B(2,3) + A(1,1)*B(3,2))/2.d0
       AoB(2,4) = (A(2,2)*B(1,2) + A(2,2)*B(2,1))/2.d0
       AoB(2,5) = (A(2,2)*B(1,3) + A(2,2)*B(3,1))/2.d0
       AoB(2,6) = (A(2,2)*B(2,3) + A(2,2)*B(3,2))/2.d0
       AoB(3,4) = (A(3,3)*B(1,2) + A(3,3)*B(2,1))/2.d0
       AoB(3,5) = (A(3,3)*B(1,3) + A(3,3)*B(3,1))/2.d0
       AoB(3,6) = (A(3,3)*B(2,3) + A(3,3)*B(3,2))/2.d0
c
       AoB(4,1) = (A(1,2)*B(1,1) + A(2,1)*B(1,1))/2.d0
       AoB(4,2) = (A(1,2)*B(2,2) + A(2,1)*B(2,2))/2.d0
       AoB(4,3) = (A(1,2)*B(3,3) + A(2,1)*B(3,3))/2.d0
       AoB(5,1) = (A(1,3)*B(1,1) + A(3,1)*B(1,1))/2.d0
       AoB(5,2) = (A(1,3)*B(2,2) + A(3,1)*B(2,2))/2.d0
       AoB(5,3) = (A(1,3)*B(3,3) + A(3,1)*B(3,3))/2.d0
       AoB(6,1) = (A(2,3)*B(1,1) + A(3,2)*B(1,1))/2.d0
       AoB(6,2) = (A(2,3)*B(2,2) + A(3,2)*B(2,2))/2.d0
       AoB(6,3) = (A(2,3)*B(3,3) + A(3,2)*B(3,3))/2.d0
c
       AoB(4,4) = (A(1,2)*B(1,2) + A(1,2)*B(2,1)
     +           + A(2,1)*B(1,2) + A(2,1)*B(2,1))/4.d0
       AoB(4,5) = (A(1,2)*B(1,3) + A(1,2)*B(3,1)
     +           + A(2,1)*B(1,3) + A(2,1)*B(3,1))/4.d0
       AoB(4,6) = (A(1,2)*B(2,3) + A(1,2)*B(3,2)
     +           + A(2,1)*B(2,3) + A(2,1)*B(3,2))/4.d0
       AoB(5,4) = (A(1,3)*B(1,2) + A(1,3)*B(2,1)
     +           + A(3,1)*B(1,2) + A(3,1)*B(2,1))/4.d0
       AoB(5,5) = (A(1,3)*B(1,3) + A(1,3)*B(3,1)
     +           + A(3,1)*B(1,3) + A(3,1)*B(3,1))/4.d0
       AoB(5,6) = (A(1,3)*B(2,3) + A(1,3)*B(3,2)
     +           + A(3,1)*B(2,3) + A(3,1)*B(3,2))/4.d0
       AoB(6,4) = (A(2,3)*B(1,2) + A(2,3)*B(2,1)
     +           + A(3,2)*B(1,2) + A(3,2)*B(2,1))/4.d0
       AoB(6,5) = (A(2,3)*B(1,3) + A(2,3)*B(3,1)
     +           + A(3,2)*B(1,3) + A(3,2)*B(3,1))/4.d0
       AoB(6,6) = (A(2,3)*B(2,3) + A(2,3)*B(3,2)
     +           + A(3,2)*B(2,3) + A(3,2)*B(3,2))/4.d0
c
       return
       end



