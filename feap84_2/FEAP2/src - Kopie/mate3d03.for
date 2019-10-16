      subroutine mate3d03(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c
c     Purpose: calculate S and C for 
c               a linear elastic transversal isotropic material law 
c
c     Inputs:
c         h1(nh) - history array h1
c         d(md)  - local d-array
c         Eps    - strains  
c         isw    - solution option from element 
c
c     Input material parameter:                
c         E_11   - Young's modulus
c         E_22   - Young's modulus
c         v_12   - Poisson's ratio 12=13
c         G_12   - Shear   modulus  G_12 = G_13
c         G_23   - Shear   modulus
c         phi    - angle phi
c
c     Outputs:
c         md = 6 - number of used data for control of d-array 
c         nh = 0 - number of history parameter at Gauss-Point
c         sig    - stresses
c         cmat   - tangent modulus
c         plout(10)- plot data    
c
c     Allocation of d-array:
c         d(1) = E_11 - Young's modulus
c         d(2) = E_22 - Young's modulus
c         d(3) = v_12 - Poisson's ratio 12=13
c         d(4) = G_12 - Shear   modulus  G_12 = G_13
c         d(5) = G_23 - Shear   modulus
c         d(6) = phi  - angle phi
c
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(8),eps(6),sig(6),cmat(6,6),plout(10)

      if(isw.eq.1) then
c....   input data
        call mati3d03(d,md,nh)
c
      else 
c....   elastic tangent modul
        call matm3d03(d,cmat)
c
c....   stresses sig(i)=cmat(i,j)*eps(j)
        call mvmul(cmat,eps,6,6,sig)
c      
      end if
      
      return
      end
c
      subroutine mati3d03(d,md,nh)
c-----------------------------------------------------------------------
c
c     mate3d03 Linear elastic transversal isotropic 
c     input material parameter 
c     
c-----------------------------------------------------------------------
c
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(6)

      md=6
      if(ior.lt.0) write(*,1001)
1001  format(' Input: E_1,E_2,nu_12,G_12,G_23,phi,')
      nh   = 0
      call dinput(d,md)

c.... v_23= 0.5*E_2/G_23 - 1
      v23   = 0.5d0*d(2)/d(5)-1.d0     
     
                  write(iow,1002) nh,(d(i),i=1,md-1),v23,d(md)
      if(ior.lt.0)write(*  ,1002) nh,(d(i),i=1,md-1),v23,d(md)
1002  format(5x,'Linear elastic transversal isotropic material data',/,
     +  5x,'length nh of h1,h2 ........',i12,/,
     +  5x,'elastic modulus E_1 .......',g12.5,/,
     +  5x,'elastic modulus E_2 .......',g12.5,/,
     +  5x,'poissons ratio v_12 .......',f12.5,/,
     +  5x,'shear modulus  G_12 .......',g12.5,/,
     +  5x,'shear modulus  G_23 .......',g12.5,/,
     +  5x,'poissons ratio v_23 .......',g12.5,/,
     +  5x,'angle phi to x_1 ..........',f12.5,/)

c.... test

      e1    = d(1)
      e2    = d(2)
      v12   = d(3)
      v21   = v12*e2/e1 

c      
c.... test 1     
      t1 = 1.d0 - v12*v21 
c.... test 2     
      t2 = 1.d0-v23*v23 
c.... test 3     
      t3 = 1.d0-v23 - 2*v12*v21 
      
                   write(iow,1003) v21,t1,t2,t3
      if(ior.lt.0) write(*  ,1003) v21,t1,t2,t3
1003  format(
     +  5x,'poissons ratio v_21 .......',g12.5,/,
     +  5x,'cond 1: 1-v12*v21>0 .......',g12.5,/,
     +  5x,'cond 2: 1-v23*v23>0 .......',g12.5,/,
     +  5x,'cond 3: 1-v23-2v12*v21>0 ..',g12.5)
      
      if(t1.le.0.d0) stop 
     +'Warning: problem material data input: choose 1-v_12^2*E_2/E_1 >0'
      if(t2.le.0.d0) stop 
     +'Warning: mat.data input: choose 1-v_23^2 >0 or  G_23 > 0.25 E_2'
      if(t3.le.0.d0) stop 
     +'Warning: problem material data input: 1-v23-2v12*v21>0'

      return
      end 
c
      subroutine matm3d03(d,cmat)
c-----------------------------------------------------------------------
c
c     mate3d03 Linear elastic transversal isotropic 
c     calculate C 
c     
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(*),cmat(6,6)
    
      e1    = d(1)
      e2    = d(2)
      v12   = d(3)
      v21   = v12*e2/e1 
      c44   = d(4)
      c55   = c44
      c66   = d(5)
      v23   = 0.5d0*d(2)/d(5)-1.d0     
      phi   = d(6)

c.... local elasticity matrix, coefficient
      cc  = 1.0d0/(1.0d0+v23)/(1.0d0-v23-2.0d0*v12*v21)
      c11 = cc*e1*(1.0-v23*v23)
      c22 = cc*e2*(1.d0-v12*v21)
      c33 = c22
      c12 = cc*e2*v12*(1.0d0+v23)
      c13 = c12
      c23 = cc*e2*(v23+v12*v21)

c.... rotation angle
      phi = phi*datan(1.0d0)/45.0
      s   = dsin(phi)
      s2  = s*s
      s3  = s*s*s
      s4  = s*s*s*s
      c   = dcos(phi)
      c2  = c*c
      c3  = c*c*c
      c4  = c*c*c*c
c.... global material matrix (x_1--x_2)
      cmat(1,1) = c4*c11 + 2.0d0*c2*s2*(c12+2.0d0*c44) + s4*c22
      cmat(1,2) = c2*s2*(c11+c22-4.0d0*c44) + (c4+s4)*c12
      cmat(1,4) = c3*s*(c11-c12-2.0d0*c44)+s3*c*(c12-c22+2.0d0*c44)
      cmat(2,2) = s4*c11 + 2.0d0*c2*s2*(c12+2.0d0*c44) + c4*c22
      cmat(2,4) = s3*c*(c11-c12-2.0d0*c44)+c3*s*(c12-c22+2.0d0*c44)
      cmat(3,4) = c*s*(c13-c23)
      cmat(4,4) = c2*s2*(c11+c22-2.0d0*c12-2.0d0*c44) + (c4+s4)*c44
c.... terms in 3-direction
      cmat(1,3) = c2*c13 + s2*c23
      cmat(2,3) = s2*c13 + c2*c23
      cmat(3,3) = c33
c.... shear terms
      cmat(5,5) = c2*c55 + s2*c66
      cmat(5,6) = c*s*(c55-c66)
      cmat(6,6) = s2*c55 + c2*c66
c.... symmetry
      cmat(2,1) = cmat(1,2)
      cmat(3,1) = cmat(1,3)
      cmat(3,2) = cmat(2,3)
      cmat(4,1) = cmat(1,4)
      cmat(4,2) = cmat(2,4)
      cmat(4,3) = cmat(3,4)
      cmat(6,5) = cmat(5,6)
c
c.... test
      do i =1,6
        if(cmat(i,i) .lt.0.d0 ) write(*,2000) i,i
      end do

2000  format(1x,'neg. val. of cmat at ',i2,i2)
      return
      end
c
