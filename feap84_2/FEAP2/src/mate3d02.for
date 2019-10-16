      subroutine mate3d02(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c
c     Purpose: calculate S and C for 
c               a linear elastic orthotropic material law 
c               with a rotation phi in 1-2 plane
c
c     Inputs:
c         h1(nh)   - history array h1
c         d(md)    - local d-array
c         Eps      - strains  
c         isw      - solution option from element 
c
c     Input material parameter:                
c         E_11     - Young's modulus
c         E_22     - Young's modulus
c         E_33     - Young's modulus
c         v_12     - Poisson's ratio 12
c         v_13     - Poisson's ratio 13
c         v_23     - Poisson's ratio 23
c         G_12     - Shear   modulus    
c         G_13     - Shear   modulus    
c         G_23     - Shear   modulus
c         phi      - angle phi
c
c     Outputs:
c         md = 10  - number of used data for control of d-array 
c         nh = 0   - number of history parameter at Gauss-Point
c         sig      - stresses
c         cmat     - tangent modulus
c         plout(10)- plot data    
c
c     Allocation of d-array:
c        d(1)  = E_11  - Young's modulus
c        d(2)  = E_22  - Young's modulus
c        d(3)  = E_33  - Young's modulus
c        d(4)  = v_12  - Poisson's ratio 12
c        d(5)  = v_13  - Poisson's ratio 13
c        d(6)  = v_23  - Poisson's ratio 23
c        d(7)  = G_12  - Shear   modulus    
c        d(8)  = G_13  - Shear   modulus    
c        d(9)  = G_23  - Shear   modulus
c        d(10) = phi   - angle phi
c
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(11),eps(6),sig(6),cmat(6,6),plout(10)

      if(isw.eq.1) then
c....   input data
        call mati3d02(d,md,nh)
c
      else 
c....   elastic tangent modul C
        call matm3d02(d,cmat)
c
c....   stresses sig(i)=cmat(i,j)*eps(j)
        call mvmul(cmat,eps,6,6,sig)
c      
      end if

      return
      end
c
      subroutine mati3d02(d,md,nh)
c-----------------------------------------------------------------------
c
c     mate3d02 Linear elastic orthotropic, rotation with alpha in 1-2 plane
c     input material parameter 
c     
c-----------------------------------------------------------------------
c
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(10)

      md=10
      if(ior.lt.0) write(*,1001)
1001  format(
     + ' Input: E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,phi')
      nh   = 0
      call dinput(d,md)
                  write(iow,1002) nh,(d(i),i=1,md)
      if(ior.lt.0)write(*  ,1002) nh,(d(i),i=1,md)
1002  format(5x,'Linear elastic orthotropic material data',/,
     +  5x,'length nh of h1,h2 ........',i12,/,
     +  5x,'elastic modulus E_1 .......',g12.5,/,
     +  5x,'elastic modulus E_2 .......',g12.5,/,
     +  5x,'elastic modulus E_3 .......',g12.5,/,
     +  5x,'poissons ratio v_12 .......',f12.5,/,
     +  5x,'poissons ratio v_13 .......',f12.5,/,
     +  5x,'poissons ratio v_23 .......',f12.5,/,
     +  5x,'shear modulus  G_12 .......',g12.5,/,
     +  5x,'shear modulus  G_13 .......',g12.5,/,
     +  5x,'shear modulus  G_23 .......',g12.5,/,
     +  5x,'angle phi to x_1 ..........',f12.5,/)

c.... test
      e1  = d(1)
      e2  = d(2)
      e3  = d(3)
      v12 = d(4)
c     v21 = v12*e2/e1 
      v13 = d(5)
      v23 = d(6)

c.... test 1     
      t1 = e2 - v23*v23*e3 
c.... test 2     
      t2 = e1 - v13*v13*e3 
c.... test 3     
      t3 = e1 - v12*v12*e2 
c.... test 4
      t4  = 1.0d0/(e1*e2-v13**2*e2*e3-v23**2*e1*e3-v12**2*e2**2
     +            -2.d0*e2*v12*v13*v23*e3)
      
                   write(iow,1003) t1,t2,t3,t4
      if(ior.lt.0) write(*  ,1003) t1,t2,t3,t4
1003  format(
     +  5x,'cond 1: e2-v23^2*e3>0 .....',g12.5,/,
     +  5x,'cond 2: e1-v13^2*e3>0 .....',g12.5,/,
     +  5x,'cond 3: e1-v12^2*e2>0......',g12.5,/,
     +  5x,'cond 4: Determinant>0......',g12.5)
      
      if(t1.le.0.d0) stop    
     +'Warning: problem material data input: choose e2-v23^2*e3 >0'
      if(t2.le.0.d0) stop 
     +'Warning: problem material data input: choose e1-v13^2*e3 >0'
      if(t3.le.0.d0) stop 
     +'Warning: problem material data input: e1-v12^2*e2>0'
      if(t4.le.0.d0) stop 
     +'Warning: problem material data input: Determinant>0'

      return
      end 
c
      subroutine matm3d02(d,cmat)
c-----------------------------------------------------------------------
c
c     mate3d02 Linear elastic orthotropic, rotation with alpha in 1-2 plane
c     calculate C 
c     
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),cmat(6,6)

c.... material data
      e1  = d(1)
      e2  = d(2)
      e3  = d(3)
      v12 = d(4)
c     v21 = v12*e2/e1 
      v13 = d(5)
      v23 = d(6)
      c44 = d(7)
      c55 = d(8)
      c66 = d(9)
      phi = d(10)

c.... local material matrix
      cc  = 1.0d0/(e1*e2-v13**2*e2*e3-v23**2*e1*e3-v12**2*e2**2
     +            -2.d0*e2*v12*v13*v23*e3)
      
      c11 =  (e2-v23**2*e3)*e1**2*cc
      c12 =  (v12*e2+v13*v23*e3)*e1*e2*cc
      c13 =  (v12*v23+v13)*e1*e2*e3*cc
      c22 =  (e1-v13**2*e3)*e2**2*cc
      c23 =  (v23*e1+v12*v13*e2)*e2*e3*cc
      c33 =  (e1-v12**2*e2)*e2*e3*cc

c.... rotation in 1-2 plane
      phi = phi*datan(1.0d0)/45.0
      s   = dsin(phi)
      s2  = s*s
      s3  = s*s*s
      s4  = s*s*s*s
      c   = dcos(phi)
      c2  = c*c
      c3  = c*c*c
      c4  = c*c*c*c

c.... global material matrix 
      cmat(1,1) = c4*c11 + 2.0d0*c2*s2*(c12+2.0d0*c44) + s4*c22
      cmat(1,2) = c2*s2*(c11+c22-4.0d0*c44) + (c4+s4)*c12
      cmat(1,3) = c2*c13 + s2*c23
      cmat(1,4) = c3*s*(c11-c12-2.0d0*c44)+s3*c*(c12-c22+2.0d0*c44)
      cmat(2,2) = s4*c11 + 2.0d0*c2*s2*(c12+2.0d0*c44) + c4*c22
      cmat(2,3) = s2*c13 + c2*c23
      cmat(2,4) = s3*c*(c11-c12-2.0d0*c44)+c3*s*(c12-c22+2.0d0*c44)
      cmat(3,3) = c33
      cmat(3,4) = c*s*(c13-c23)
      cmat(4,4) = c2*s2*(c11+c22-2.0d0*c12-2.0d0*c44) + (c4+s4)*c44
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

c.... test
      do i =1,6
        if(cmat(i,i) .lt.0.d0 ) then
          write(*,2000) i,i
        end if 
      end do
      
2000  format(1x,'neg. val. of Cmat at ',i2,i2)

c     call mprint(cmat,6,6,6,'cmat')

      return
      end

