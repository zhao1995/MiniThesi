      subroutine mate3d12(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c
c     Purpose: calculate Sig(S,-D) and Cmat for 
c               a linear piezoelectric material law
c
c               | S  |   | IC    -ie | | E  |
c               |    | = |           |*|    |
c               |-vD |   |-ie^T  -eps| | vE |
c
c     Inputs:
c         h1(nh)    - history array h1
c         d(md)     - local d-array
c         Eps       - gradient fields
c         isw       - solution option from element 
c
c     Input material parameter:                
c         E         - Young's modulus
c         v         - Poisson's ratio
c         e13       - piezo moduli
c         e33       - 
c         e15       - 
c         eps       - dielectric constant eps
c
c      Outputs:
c         md = 6   - number of used data for control of d-array 
c         nh = 0   - length of history array at Gauss-Point
c         sig      - generlized stresses
c         cmat     - tangent modulus
c         plout(10)- plot data    
c
c     Allocation of d-array:
c         d(1) = E  - Young's modulus
c         d(2) = v  - Poisson's ratio
c         d(3) = e13 - piezo moduli
c         d(4) = e33 - 
c         d(5) = e15 - 
c         d(6) = eps - dielectric constant eps
c
c     (c) s.klinkel                                           March,2011
c-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(*),eps(nsig),sig(nsig),cmat(nsig,nsig),
     +          plout(10),xgp(3)
c
      if(isw.eq.1) then
c....   input data
        call mati3d12(d,md,nh)
      else
c....   elastic tangent modul
        call matm3d12(d,cmat)
c....   stresses sig(i)=cmat(i,j)*eps(j)
        call mvmul(cmat,eps,9,9,sig)
      end if
c      
      return
      end
c
      subroutine mati3d12(d,md,nh)
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(6)

      md=6
      if(ior.lt.0) write(*,1001)
1001  format(
     + ' Input:E,nu,e_13,e_33,e_15,epsilon')
      nh = 0
      call dinput(d,md)
                  write(iow,1002) nh,(d(i),i=1,md)
      if(ior.lt.0)write(*  ,1002) nh,(d(i),i=1,md)
1002  format(5x,'linear piezolectricity',/,
     +  5x,'length nh of h1,h2 ........',i12,/,
     +  5x,'elastic modulus E .........',g12.5,/,
     +  5x,'Poissons ratio  v .........',g12.5,/,
     +  5x,'Piezo moduli e_13..........',g12.5,/,
     +  5x,'............ e_33..........',g12.5,/,
     +  5x,'............ e_15..........',g12.5,/,
     +  5x,'Dielectric constant eps....',g12.5)
c
      return
      end 
c
      subroutine matm3d12(d,Cmat)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(6),Cmat(9,9)
      Cmat = 0.d0
c...  elasticity matrix
c.... local iostropy elasticity matrix
      xnu= d(2)
      cc = d(1)/((1.0d0+xnu)*(1.0d0-2.0d0*xnu))
      tt = (1.0d0-2.0d0*xnu)/2.0d0
c.... material parameters 
      Cmat(1,1) =cc*(1.0d0-xnu)  
      Cmat(1,2) =cc*xnu          
      Cmat(1,3) =cc*xnu          
      Cmat(2,2) =cc*(1.0d0-xnu)  
      Cmat(2,3) =cc*xnu          
      Cmat(3,3) =cc*(1.0d0-xnu)  
      Cmat(4,4) =cc*tt           
      Cmat(5,5) =cc*tt           
      Cmat(6,6) =cc*tt           
c.... symmetry
      Cmat(2,1) = cmat(1,2)
      Cmat(3,1) = cmat(1,3)
      Cmat(3,2) = cmat(2,3)
c
c...  couple matrix -e^t
      Cmat(1,9) = -d(3)
      Cmat(2,9) = -d(3)
      Cmat(3,9) = -d(4)
      Cmat(5,7) = -d(5)
      Cmat(6,8) = -d(5)
      Cmat(9,1) = Cmat(1,9)
      Cmat(9,2) = Cmat(2,9)
      Cmat(9,3) = Cmat(3,9)
      Cmat(7,5) = Cmat(5,7)
      Cmat(8,6) = Cmat(6,8)
c
c...  dielectric matrix -eps
      Cmat(7,7) = -d(6)
      Cmat(8,8) = -d(6)
      Cmat(9,9) = -d(6)
c
      return
      end
