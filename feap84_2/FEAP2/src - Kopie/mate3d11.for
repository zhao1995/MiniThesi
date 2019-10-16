      subroutine mate3d11(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c     Purpose: calculate S and C for 
c              a small strain visco-elastic material law 
c
c     Inputs:
c         h1(nh)         - history array h1
c         d(md)          - local d-array
c         Eps            - strains  
c         isw            - solution option from element 
c
c     Input material parameters:   
c      record 1                
c         E              - Young's modulus
c         nu             - Poisson ratio
c         alpha          - Thermal expansion coefficient
c         nv             - Number of viscoelastic terms 
c         mn             - 0 = deviatoric strains, 1 = total strains 
c      record i=1,nv  for each element  
c         mu_i           - Viscoelastic shear parameter
c         tau_i i=1,nv   - Viscoelastic relaxation time 
c                      
c     Outputs:         
c         md = 5+2*nv    - number of used data for control of d-array 
c         nh = 6*(nv+1)  - number of history parameters at Gauss point
c         h2(nh)         - history array h2
c         Sig            - 2nd Piola-Kirchhoff Stress
c         Cmat           - consistent tangent matrix
c         plout(10)      - plot data    
c
c     Allocation of d-array:
c         d( 1)  E       - Young's modulus
c         d( 2)  nu      - Poisson ratio
c         d( 3)  alpha   - Thermal expansion coefficient
c         d( 4)  nv      - Number of viscoelastic terms 
c         d( 5)  mn      - Material number: 0 = deviatoric strains, 1 = total strains 
c         d( 6)  mu_1    - Viscoelastic shear parameter
c         d( 7)  tau_1   - Viscoelastic relaxation time 
c         d( 8)  mu_2    - Viscoelastic shear parameter 
c         d( 9)  tau_2   - Viscoelastic relaxation time 
c         ... 
c
c     References: 
c      J.C. Simo and T.J.R. Hughes: Computational Inelasticity, Springer 
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),h1(*),h2(*),eps(*),sig(*),cmat(6,6),plout(10)
c
      if(isw.eq.1) then
c....   input data
        call mati3d11(d,md,nh)
c
      else
c....   Get internal variables from h1 at time t_n 
        do i = 1,nh
          h2(i) = h1(i)
        end do
c       compute stresses and tangent matrix
c        
        mn=d(5)
        if (mn.eq.0) then            
          call viscoe (d,tgp,eps,h2(1),h2(7),6,sig,cmat) ! only deviatoric parts
        else  
          call viscoe1(d,tgp,eps,h2(1),h2(7),6,sig,cmat) ! deviatoric + volumetric parts
        end if
c        
        do i = 1,10 
         plout(i) = h2(6+i)
        end do
c
      end if
c
      return
      end
c
      subroutine mati3d11(d,md,nh)
c-----------------------------------------------------------------------
c
c     isotropic small strain viscoelastic model
c     input material parameters 
c     
c-----------------------------------------------------------------------
c
      USE cdat1
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(*)
c
      if(ior.lt.0) write(*,1001)
1001  format('Input: E,nu,alpha_t,nv,mn >')
      call dinput(d,5) 
      nv = d(4)
      mn = d(5) 
      md = 5 + 2*nv
      nh = 6 + 6*nv   ! 6 en(i) + 6*nv qi(i,n) 

                  write(iow,1002) nh,(d(i),i=1,3),nv,mn
      if(ior.lt.0)write(*  ,1002) nh,(d(i),i=1,3),nv,mn
1002  format(5x,'isotropic small strain viscoelastic model',/,
     +  5x,'length nh of h1,h2 ..................',i12,/,
     +  5x,'Youngs modulus E.....................',g12.5,/,
     +  5x,'Poisson ratio v......................',g12.5,/,
     +  5x,'Thermal expansion coefficient........',g12.5,/,
     +  5x,'Number of viscoelastic models........',i10,/,
     +  5x,'0=deviatoric strains 1=total strains.',i10,/,
     +  5x,'El.No...shear parameter mu...relaxation time tau')
      
c.... length d-array
      if(md.gt.ndd) then
                      write(iow,1003) ndd,md
         if(ior.lt.0) write(iow,1003) ndd,md
1003     format(1x,'darray is set to  ',i4,' values',/,
     1          1x,'darray needs      ',i4,' values')
         stop
      end if

      if(ior.lt.0) write(*,1004)
1004  format('Input: mui,taui,i=1,nv >')
      do i=1,nv   
        call dinput(d(2*i+4),2)
                    write(iow,1005) i,(d(2*i+k),k=4,5)
        if(ior.lt.0)write(*  ,1005) i,(d(2*i+k),k=4,5)
      end do
1005  format(5x,i4,g12.5,1x,g12.5)

      return
      end 
c
      subroutine viscoe(d,ta,eps,en,qi,ntm,sig,cmat)
c-----------------------------------------------------------------------
c      Purpose: Linear viscoelastic model (shear only)
c              
c     (c)  R.L. Taylor  UCB 
c-----------------------------------------------------------------------

c      Inputs:
c         d(*)    - Material parameters
c         ta      - Temperature
c         eps(*)  - Strains at t_n+1
c         en(*)   - Strains at t_n
c         qi(*)   - Viscoelastic strain
c         ntm     - Number of terms

c      Outputs:
c         sig(6)  - Stresses
c         cmat(6,6) - Viscoelastic moduli
c-----------------------------------------------------------------------
      USE tdata
      implicit  none


      integer   i,j, n,nv,ntm
      real*8    G,Gg,K,Kg,Kth, gfac,exp_n,mu_0,mu_n,dq_n,dtau, theta
      real*8    alpha,ta, d(*),eps(*),en(*),qi(ntm,*),ee(6)
      real*8    sig(6),cmat(6,6), hvisc

c     Set elastic parameters for G (mu) and lambda

      G     = d(1)/(2.d0*(1.d0 + d(2)))
      K     = d(1)/(3.d0*(1.d0 - 2.d0*d(2)))
      alpha = d(3)

c     Compute volumetric strain and deviatoric components

      theta = (eps(1) + eps(2) + eps(3))/3.d0
      do i = 1,3
        ee(i  ) = eps(i) - theta
      end do ! i
      do i = 4,ntm
        ee(i) = eps(i)*0.5d0
      end do ! i

c     Set properties for integrating the q_i terms

      mu_0 = 0.0d0
      gfac = 0.0d0
      sig  = 0.0d0

      nv   = d(4)
      do n = 1,nv 
        mu_n  =    d(2*n+4)
        dtau  = dt/d(2*n+5)
        exp_n = dexp(-dtau)

        if(dtau.lt.1.d-04) then
          hvisc = 1.d0 - 0.5d0*dtau*(1.d0 - dtau/3.d0*(1.d0
     +                 - 0.25d0*dtau*(1.d0 - 0.2d0*dtau)))
        else
          hvisc = (1.d0 - exp_n)/dtau
        end if

        dq_n = mu_n * hvisc
        gfac = gfac + dq_n
        mu_0 = mu_0 + mu_n

c       Update history and compute viscoelastic deviatoric stress

        do i = 1,ntm
          qi(i,n) = exp_n*qi(i,n) + dq_n*(ee(i) - en(i))
          sig(i)  = sig(i) + qi(i,n)
        end do ! i
      end do ! n

c     Finish updates and save the strains

      mu_0 = 1.d0 - mu_0
      gfac = gfac + mu_0
      do i = 1,ntm
        sig(i) = 2.d0*G*(mu_0*ee(i) + sig(i))
        en(i)  = ee(i)
      end do ! i

c     Add elastic bulk term

      Kth = K*(theta*3.0d0 - alpha*ta)
      do i = 1,3
        sig(i) = sig(i) + Kth
      end do ! i

c     Set tangent parameters

      Gg = G*gfac
      Kg = K - 2.d0*Gg/3.d0
      do j =1,3
        do i = 1,3
          cmat(i,j) = Kg
        end do ! i
        cmat(j,j) = cmat(j,j) + 2.d0*Gg
      end do ! i

      do i = 4,ntm
        cmat(i,i) = Gg
      end do ! i

      return
      end
c
      subroutine viscoe1(d,ta,eps,en,qi,ntm,sig,cmat)
c-----------------------------------------------------------------------
c      Purpose: Linear viscoelastic model
c     
c-----------------------------------------------------------------------

c      Inputs:
c         d(*)    - Material parameters
c         ta      - Temperature
c         eps(*)  - Strains at t_n+1
c         en(*)   - Strains at t_n
c         qi(*)   - Viscoelastic strain
c         ntm     - Number of terms

c      Outputs:
c         sig(6)  - Stresses
c         cmat(6,6) - Viscoelastic moduli
c
c      Comments:
c      alpha not implemented!!
c-----------------------------------------------------------------------
      USE tdata
      implicit double precision (a-h,o-z)
      dimension d(*),sig(6),eps(6),cmat(6,6),en(6),qi(ntm,*)
      

      call pzero(cmat,6*6)
      call pzero(sig,6)

c     Set elastic parameters 

      E     = d(1)
      xnue  = d(2)
c      alpha = d(3)
      call matm3d01(cmat,E,xnue)

c     Integrating the q_i terms

      rmu_0 = 0.0d0

      nv  = d(4)
      do n = 1,nv
        rmu_n = d(2*n+4)
        dtau  = dt/d(2*n+5)
        exp_n = dexp(-dtau)

        if(dtau.lt.1.d-04) then
          hvisc = 1.d0 - 0.5d0*dtau*(1.d0 - dtau/3.d0*(1.d0
     +                 - 0.25d0*dtau*(1.d0 - 0.2d0*dtau)))
        else
          hvisc = (1.d0 - exp_n)/dtau
        end if

        dq_n  = rmu_n * hvisc
        rmu_0 = rmu_0 + rmu_n

        do i = 1,ntm
          qi(i,n) = exp_n*qi(i,n) + dq_n*(eps(i) - en(i))
          sig(i)  = sig(i) + qi(i,n)
        end do ! i
      end do ! n

      rmu_0 = 1.d0 - rmu_0

c     compute stresses and moduli

      sig  = matmul(cmat,(rmu_0*eps+sig))
      cmat = (rmu_0 + dq_n)*cmat
      en   = eps
c
      return
      end
