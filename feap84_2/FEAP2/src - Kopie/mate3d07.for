      subroutine mate3d07(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c
c     Purpose: calculate S and C for fiber reinforced material
c              Method of Cells (Aboudi)
c
c     Inputs:
c         h1(nh) - history array h1
c         d(md)  - local d-array
c         eps    - strains  
c         isw    - solution option from element 
c
c     Input material parameter:                
c       card 1  
c         E_11 fiber
c         E_22 fiber
c         E_33 fiber
c         v_12 fiber
c         v_13 fiber
c         v_23 fiber
c         G_12 fiber
c         G_13 fiber
c         G_23 fiber
c         alpha_11 fiber
c         alpha_22 fiber
c         alpha_33 fiber
c         e matrix
c         poissons ratio matrix
c         shear modulus matrix
c         alpha matrix
c          
c       card 2          
c         spring stiffness in 1 dir
c         spring stiffness in 2 dir
c         spring stiffness in 3 dir
c         temperature change
c         fiber volume fraction
c         tensile strength-fiber
c         compressive strength-fiber
c         ultimate matrix tensile strength
c         ultimate matrix compressive strength
c         ultimate matrix shear strength
c         number of cells  (4)
c         fiber angle phi
c
c     Outputs:
c         md = 125 - number of used data for control of d-array 
c         nh = 30  - number of history parameters at Gauss point
c         h2(nh)   - history array h2
c         sig      - stresses
c         cmat     - tangent modulus
c         plout(10)- plot data    
c
c     Allocation of d-array:
c         d( 1) = E_11 fiber
c         d( 2) = E_22 fiber
c         d( 3) = E_33 fiber
c         d( 4) = v_12 fiber
c         d( 5) = v_13 fiber
c         d( 6) = v_23 fiber
c         d( 7) = G_12 fiber
c         d( 8) = G_13 fiber
c         d( 9) = G_23 fiber
c         d(10) = alpha_11 fiber
c         d(11) = alpha_22 fiber
c         d(12) = alpha_33 fiber
c         d(13) = e matrix
c         d(14) = poissons ratio matrix
c         d(15) = shear modulus matrix
c         d(16) = alpha matrix
c         d(17) = spring stiffness in 1 dir
c         d(18) = spring stiffness in 2 dir
c         d(19) = spring stiffness in 3 dir
c         d(20) = temperature change
c         d(21) = fiber volume fraction
c         d(22) = tensile strength-fiber
c         d(23) = compressive strength-fiber
c         d(24) = ultimate matrix tensile strength
c         d(25) = ultimate matrix compressive strength
c         d(26) = ultimate matrix shear strength
c         d(27) = number of cells  (4)
c         d(28) = fiber angle phi
c
c-----------------------------------------------------------------------
c              
      USE debugs
      USE iofile
      implicit double precision(a-h,o-z)
      dimension h1(*), h2(*), d(*), eps(6), sig(6), cmat(6,6),
     +          crds(3,24), rhs(18), amatrx(18,18), smatrx(6,6),  
     +          stiff(72,72), frhs(72),  svars(12),  
     +          astar(48,6), stiffred(48,48), lm(18), te(6,6),
     +          u(18), ub(48), ubkib(24), ui(24), lhol(4,4), s(6),
     +          epsl(6),sigl(6),cmatl(6,6),
     +          plout(10) 
      
      data lhol  / 1,3,2,4,  3,1,2,4,  1,3,4,2, 3,1,4,2 /
c              
c              
      if(isw.eq.1) then          
        call mati3d07(d,md,nh)    
c              
      else

       ndofel = 18
       mcrd   = 3   
       nnode  = 24
       lstop  = 24
       ndof   = 72
       ndn    = 28
c      
c        transform eps to local coordinate system
c      
       phi   = d(28)
       phi = phi*datan(1.0d0)/45.d0
       sn  = dsin(phi)
       s2  = sn*sn
       c   = dcos(phi)
       c2  = c*c
       sc  = sn*c
       sc2 = 2.d0*sc
       call pzero(te,6*6)
       te(1,1) =  c2  
       te(1,2) =  s2  
       te(1,4) =  sc
       te(2,1) =  s2  
       te(2,2) =  c2
       te(2,4) = -sc  
       te(3,3) =  1.d0
       te(4,1) = -sc2  
       te(4,2) =  sc2  
       te(4,4) =  c2-s2  
       te(5,5) =  c
       te(5,6) =  sn
       te(6,5) = -sn
       te(6,6) =  c
       
       epsl = matmul(te,eps)
       
c      
c       compute the size of the cells based on the input volume fraction     
c      
       numel1 = d(27)
       volf = d(21)
       volm = 1.d0-volf
       if (numel1.eq.4) then 
        df  = 1.d0
        rl1 = dsqrt(volf/df)
        rad = dsqrt((4.d0*rl1*rl1)+(4.d0*volm))
        rl2 = ((-2.d0*rl1)+rad)/2.d0
       endif
c      
       call pzero (astar,6*48) 
       call pzero (frhs,72) 
       call pzero (stiff,72*72) 
c      
c.....  loop through the elements
c      
       kk = ndn
       do m = 1,numel1
       
c       topology 
        do i = 1,18
         lm(i) =  d(kk+i)
        enddo
       
        if (m.eq.1) idofcon = lm(9)
c      
c....    coordinate matrices
c       
        call coor(m, lstop, rl1, rl2, df, lm, crds, astar)
c       
c        get the stiffness matrix for this element
c       
        jelem=0
        call pzero (u,18) 
        call aboudi(rhs, amatrx, svars, ndofel, d, crds, mcrd,
     +              nnode, u, jelem, m)
c       
c        assemble the global stiffness matrix and force vector
c       
        do j=1,18
         do i=1,18
          stiff(lm(i),lm(j)) = stiff(lm(i),lm(j)) + amatrx(i,j)
         enddo 
         frhs(lm(j)) = frhs(lm(j)) + rhs(j)
        enddo
       
        kk = kk + 18
       enddo   ! end loop over elements
c      
c       put in spring stiffness (set up for only 4)
c      
       kk = ndn+18*4
       do m = 1,4
       
        do i = 1,6
         lm(i) = d(kk+i)
        enddo
         
        call pzero(u,6)
        call spring(rhs, smatrx, d, u)
c       
        do j=1,6
         frhs(lm(j)) = frhs(lm(j)) + rhs(j)
         do i=1,6
          stiff(lm(i),lm(j)) = stiff(lm(i),lm(j)) + smatrx(i,j)
         enddo
        enddo
       
        kk = kk + 6
       enddo
c      
c        constrain out dof 3 for node 3,
c        row and column are set to 0 with 1 on the diagonal
c      
       do i = 1,ndof
        stiff(i,idofcon) = 0.d0
        stiff(idofcon,i) = 0.d0
       enddo 
       stiff(idofcon,idofcon) = 1.d0
       frhs(idofcon)=0.d0
c      
c      condense out the internal degrees of freedom by gauss elimination
c      
       do j=1,lstop
        stiffj=stiff(j,j)
        do i=j,ndof
         stiff(j,i)=stiff(j,i)/stiffj
        enddo
        frhs(j)=frhs(j)/stiffj
        do jj=j+1, ndof
         stfjj=stiff(jj,j)
         do i=j,ndof
          stiff(jj,i)=stiff(jj,i)-stiff(j,i)*stfjj
         enddo
         frhs(jj)=frhs(jj)-frhs(j)*stfjj
        enddo
       enddo
c      
       do i = 1,48
        do j = 1,48
         stiffred(i,j) = stiff(lstop+i,lstop+j)
        enddo
       enddo 
c      
c       compute the final stiffness matrix:  a^T  k  a
c      
       cmatl = matmul(transpose(astar),matmul(stiffred,astar))
       cmat  = matmul(transpose(te),matmul(cmatl,te))
c      
       do i=1,48
        frhs(i+lstop) = -frhs(i+lstop)
       enddo
c      
c       compute the boundary displacements from the averaged strains
c
       ub = matmul(astar,epsl) 
c      
c       compute the boundary forces
c      
       do j=1,48
        do i = 1+lstop,ndof
         frhs(i) = frhs(i) + stiff(i,lstop+j)*ub(j)
        enddo
       enddo
c      
c       compute the averaged stresses from the boundary forces
c      
       do i=1,6
        sigl(i)=0.d0
        do j=1,48
         sigl(i) = sigl(i) + astar(j,i)*frhs(lstop+j)
        enddo
       enddo
       sig = matmul(transpose(te),sigl)
c      
       if(debug.eq.1.and.ngp.eq.1.and.lgp.eq.1) then
        if(dabs(phi).gt.1.d-2)write(iow,*)' !!! phi is not zero !!! '
        cmatl = cmat
        call invert(cmatl,6,6)
        epsl = -matmul(cmatl,sig)
        call mprint(cmat,6,6,6,'cmat ')
        call mprint(epsl,6,1,6,'epst ')
        e11 =  1.d0/cmatl(1,1)
        e22 =  1.d0/cmatl(2,2)  
        v12 = -cmatl(1,2)*e11
        v23 = -cmatl(2,3)*e22
        g12 =  1.d0/cmatl(4,4)
        g23 =  1.d0/cmatl(6,6)
        write(iow,*)'e1=',e11
        write(iow,*)'e2=',e22
        write(iow,*)'v1=',v12
        write(iow,*)'v2=',v23
        write(iow,*)'g1=',g12
        write(iow,*)'g2=',g23
       endif
c      
c       compute the internal displacements
c      
       call pzero(ubkib,lstop)
       do j=1,48
        do i=1,lstop
         ubkib(i) = ubkib(i) + stiff(i,lstop+j)*ub(j)
        enddo
       enddo
       do i=1,lstop
        ui(i) = frhs(i)-ubkib(i)
       enddo
       
       ui(lstop) = ui(lstop)/stiff(lstop,lstop)
       do i = lstop-1, 1, -1
        do j = i+1, lstop
         ui(i) = ui(i) - ui(j)*stiff(i,j)
        enddo
       enddo
c      
c       call the element routines again passing in the displacements
c       to retrieve the element stresses and strains 
c      
       do m=1,numel1
        m1 = m-1
        m6 = m1*6
        m12 = m1*12
       
        do i=1,6
         u(i) = ub(m12+i)
        enddo
c        
        do i=1,4
         if (lhol(i,m).lt.3) then
          u((i-1)*3+7) = ui(1+(lhol(i,m)-1)*3+m6)
          u((i-1)*3+8) = ui(2+(lhol(i,m)-1)*3+m6)
          u((i-1)*3+9) = ui(3+(lhol(i,m)-1)*3+m6)
         else            
          u((i-1)*3+7) = ub(7+(lhol(i,m)-3)*3+m12)
          u((i-1)*3+8) = ub(8+(lhol(i,m)-3)*3+m12)
          u((i-1)*3+9) = ub(9+(lhol(i,m)-3)*3+m12)
         endif
        enddo
c       
        jelem=1
        call aboudi(rhs, amatrx, svars, ndofel, d, crds, mcrd,
     +              nnode, u, jelem, m)
       
        do i=1,6
         h2(i+m6) = svars(i)
        enddo
c      
       enddo !  end loop over elements
       
       
       if (h2(1).gt.0.0) then
        h2(25) = h2(1)/d(22)
       else 
        h2(25) = -h2(1)/d(23)
       endif
       h2(26) = 0.0
c      
c      von mises stress - fiber cell
c      
       s12 = h2(1) - h2(2)
       s23 = h2(2) - h2(3)
       s31 = h2(3) - h2(1)
       
       h2(27) = dsqrt( 0.5d0 * (s12*s12  + s23*s23 + s31*s31) 
     +               +  3.d0 * (h2(4)*h2(4)+h2(5)*h2(5)+ h2(6)*h2(6)) )
       
       do i=1,3
        
        do j=1,6
         s(j) = h2(i*6+j)
        enddo
        
        oper = dsqrt(((s(2)-s(3))/2.d0)**2+(s(6)*s(6)))
        spr1 = (s(2)+s(3))/2+oper
        spr2 = (s(2)+s(3))/2-oper
        
        if (dabs(spr2).gt.dabs(spr1)) then
         spr1=spr2
        endif
        
        if (spr1.le.0.0) then
         stren = d(25)
        else
         stren = d(24)
        endif
        
        oper = ((spr1*spr1)/(stren*stren))+(((s(4)*s(4))+
     +          (s(5)*s(5)))/(d(26)*d(26)))
        
        if (oper.gt.h2(26)) h2(26) = oper
c       
c       von mises stress  - matrix cells
c       
        s12 = s(1) - s(2)
        s23 = s(2) - s(3)
        s31 = s(3) - s(1)
       
        h2(i+27) = dsqrt(0.5 * (s12*s12 + s23*s23 + s31*s31)
     +                + 3.d0 * (s(4)*s(4)+s(5)*s(5)+s(6)*s(6)))
       enddo
c      
      endif
c
      return
      end
c
      subroutine aboudi(rhs, amatrx, svars, ndofel, d, crds, mcrd,
     +                  nnode, u, jelem, m)
      implicit real*8(a-h,o-z)
      dimension rhs(ndofel), amatrx(ndofel,ndofel), svars(*), d(*),
     +          crds(mcrd, nnode), u(ndofel), b(6,18), corth(6,6), 
     +          strain(6), tstr(6), thstr(6), stress(6)
c
c      6 node cell element
c      the vectors are arranged in the following manner:
c         dof 1 at node 1,
c         dof 2 at node 1,
c         dof 3 at node 1,
c         dof 1 at node 2,
c         dof 2 at node 2,
c         dof 3 at node 2,  and so on.  
c
c      the program is called in parts:
c      the first part calls the element to obtain the element stiffness 
c      in the second part, the element stresses and strains are returned   
c      given the displacement field, or state of strain, that exists
c      at the boundary    
c
c      svars(1) to svars(6) = stresses at integration points
c      svars(7) to svars(12) = strains at integration points
c
c     setup b-matrix
c
      m1=m-1
      d1 = 1.d0/(dabs(crds(1, 1+m1*6)-crds(1,2+m1*6)))
      d2 = 1.d0/(dabs(crds(2, 3+m1*6)-crds(2,4+m1*6)))
      d3 = 1.d0/(dabs(crds(3, 5+m1*6)-crds(3,6+m1*6)))
      prod = 1.d0/(d1*d2*d3)

      call pzero(b,6*18)

      b(1,1)  =  d1
      b(1,4)  = -d1
      b(2,8)  =  d2
      b(2,11) = -d2
      b(3,15) =  d3
      b(3,18) = -d3
      b(4,2)  =  d1
      b(4,5)  = -d1
      b(4,7)  =  d2
      b(4,10) = -d2
      b(5,3)  =  d1
      b(5,6)  = -d1
      b(5,13) =  d3
      b(5,16) = -d3
      b(6,9)  =  d2
      b(6,12) = -d2
      b(6,14) =  d3
      b(6,17) = -d3
c
c      constitutive matrix corth
c
      call orthotropic(d, corth, m)
          
      if (jelem.eq.1) go to 100

      amatrx = matmul(transpose(b),matmul(corth,b))*prod

c     calculate the total strain from the displacements
100   continue
      tstr = matmul(b,u) 
c
c     calculate the thermal strain in the element
c
      call pzero(thstr,6)
      if (m.eq.1) then
       do j=1,3
        thstr(j) = d(j+9)*d(20)
       enddo
      else
       do j=1,3
        thstr(j) = d(16)*d(20)
       enddo
      endif
c
c     substract the thermal strains from the total strain
c
       strain = tstr - thstr
c
c     calculate the stresses from the strain
c
       stress = matmul(corth,strain)
c
c     compute the force vector, rhs
c
       rhs = - matmul(transpose(b),stress)*prod
c
c     set state variables
c
      do i=1,6
       svars(i) = stress(i)
       svars(i+6) = strain(i)
      enddo
c      
      return
      end
c
      subroutine coor(m, lstop, rl1, rl2, df, lm, crds, astar)
      implicit real*8(a-h,o-z)
      dimension crds(3,24), astar(48,6), lm(18)
c
       m1=m-1
c
       if (m.lt.5) then
        crds(1,1+m1*6)=0.0
        crds(1,2+m1*6)=-df
        do i=3,6
         crds(1,i+m1*6)=-0.5d0*df
        enddo
       else
        crds(1,1+m1*6)=-df
        crds(1,2+m1*6)=-2.d0*df
        do i=3,6
         crds(1,i+m1*6)=-1.5d0*df
        enddo
       endif
       if (mod(m,2).gt.0) then
        iodd=0
       else
        iodd=1
       endif
       if ((m.eq.1).or.(m.eq.2).or.(m.eq.5).or.(m.eq.6)) then
        ihigh=0
       else
        ihigh=1
       endif
       crds(2,1+m1*6) = 0.5d0*rl1+iodd*(0.5d0*rl1+0.5d0*rl2)
       crds(2,2+m1*6) = crds(2,1+m1*6)
       crds(2,5+m1*6) = crds(2,1+m1*6)
       crds(2,6+m1*6) = crds(2,1+m1*6)
       crds(2,3+m1*6) = rl1+iodd*rl2
       crds(2,4+m1*6) = iodd*rl1
       crds(3,1+m1*6) = 0.5d0*rl1+ihigh*(0.5d0*rl1+0.5*rl2)
       do i=2,4
        crds(3,i+m1*6) = crds(3,1+m1*6)
       enddo
       crds(3,5+m1*6) = rl1+ihigh*(rl2)
       crds(3,6+m1*6) = ihigh*rl1
c      
c      setup a matrix, which converts the boundary displacements 
c      into strains
c      
       do i=1,18,3
        if (lm(i).gt.lstop) then
         li = lm(i)-lstop
         ii = 1+i/3
         astar(li,1)   = crds(1,ii+m1*6)
         astar(li,4)   = 0.5d0*crds(2,ii+m1*6)
         astar(li,5)   = 0.5d0*crds(3,ii+m1*6)
         astar(li+1,2) = crds(2,ii+m1*6)
         astar(li+1,4) = 0.5d0*crds(1,ii+m1*6)
         astar(li+1,6) = 0.5d0*crds(3,ii+m1*6)
         astar(li+2,3) = crds(3,ii+m1*6)
         astar(li+2,5) = 0.5d0*crds(1,ii+m1*6)
         astar(li+2,6) = 0.5d0*crds(2,ii+m1*6)
        endif
       enddo
c
      return
      end
c
      subroutine spring(rhs, smatrx, d, u)
      implicit real*8(a-h,o-z)
      dimension rhs(*), smatrx(6,6), d(*), u(6)
c
      call pzero(smatrx,6*6)
c
      smatrx(1,1) =  d(17)
      smatrx(1,4) = -d(17)
      smatrx(2,2) =  d(18)
      smatrx(2,5) = -d(18)
      smatrx(3,3) =  d(19)
      smatrx(3,6) = -d(19)
      smatrx(4,1) = smatrx(1,4)
      smatrx(4,4) = smatrx(1,1)
      smatrx(5,2) = smatrx(2,5)
      smatrx(5,5) = smatrx(2,2)
      smatrx(6,3) = smatrx(3,6)
      smatrx(6,6) = smatrx(3,3)
c
c     compute the forc vector, rhs
c
      do i = 1,6
       rhs(i) = 0.d0
       do j = 1,6 
        rhs(i) = rhs(i) - smatrx(i,j)*u(j)
       enddo
      enddo  
c
      return
      end
c
      subroutine orthotropic(d, corth, m)
      implicit real*8(a-h,o-z)
      dimension d(*), corth(6,6)
c
c     orthotropic material
c
      if(m.eq.1) then      ! orthotropic fiber
        e11 = d(1)
        e22 = d(2)
        e33 = d(3)
        v12 = d(4)
        v13 = d(5)
        v23 = d(6)
        v21 = v12*(d(2)/d(1))
        v31 = v13*(d(3)/d(1))
        v32 = v23*(d(3)/d(2))
        g12 = d(7)
        g13 = d(8)
        g23 = d(9)
      else                ! isotropic matrix 
        e11 = d(13)
        e22 = e11
        e33 = e11
        v12 = d(14)
        v13 = v12
        v23 = v12
        g12 = d(15)
        g13 = g12
        g23 = g12
        v21 = v12
        v31 = v13
        v32 = v23
      endif
      det = 1.d0-(v12*v21)-(v23*v32)-(v13*v31)-(2.0*v12*v23*v31)
      
      call pzero(corth,6*6)
       
      corth(1,1) = (1-(v23*v32))*e11/det
      corth(1,2) = (v12+(v13*v32))*e22/det
      corth(1,3) = (v13+(v23*v12))*e33/det
      corth(2,2) = (1-(v31*v13))*e22/det
      corth(2,3) = (v23+(v12*v13))*e33/det
      corth(3,3) = (1-(v12*v21))*e33/det
      corth(2,1) = corth(1,2)
      corth(3,1) = corth(1,3)
      corth(3,2) = corth(2,3)
      corth(4,4) = g12
      corth(5,5) = g13
      corth(6,6) = g23
c      
      return
      end
c
      subroutine mati3d07(d,md,nh)
c-----------------------------------------------------------------------
c
c     input material parameters and compute topology
c     
c-----------------------------------------------------------------------
c
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(*), lhol(4,4), lhols(2,4), lm(18)
      data lhol    / 1,3,2,4,  3,1,2,4,  1,3,4,2,  3,1,4,2 /
      data lhols  / 0,6,  3,15,  9,21, 12,18 /
c
      ndn = 28
      md = ndn + 18*4 + 6*4 
      if(ior.lt.0) write(*,1001)
1001  format(' Input: material parameters > ')
      call dinput(d(1),16)
      call dinput(d(17),12)
      
      numel1 = d(27)  
      volf = d(21)
      if ((numel1.ne.4).and.(numel1.ne.8)) stop 'numel1 = 4 or 8 '
      if ((volf.lt.0.0).or.(volf.gt.1.0)) stop ' 0 <  volf < 1  '

      nh   = numel1*6 + 6
c....  compute lm-array (topology) 
c     for cell elements      
      kk = ndn
      do m = 1,4
       m1 = m-1
       do i=1,3
        lm(i)   = i+12+(m*12)
        lm(i+3) = i+15+(m*12)
       enddo
       do j=1,4
        do i=1,3
         if (lhol(j,m).eq.1) then
          lm(i+3+j*3) = i+(m1*6)
         elseif (lhol(j,m).eq.2) then
          lm(i+3+j*3) = i+3+(m1*6)
         elseif (lhol(j,m).eq.3) then
          lm(i+3+j*3) = i+18+(m*12)
         else 
          lm(i+3+j*3) = i+21+(m*12)
         endif
        enddo
       enddo
c      store in d-array      
       do i = 1,18
        d(kk+i) = lm(i)
       enddo
       kk = kk + 18
      enddo  

c     for spring elements
      kk = ndn+18*4
      do m = 1,4
       do i=1,3
        lm(i)   = i+lhols(1,m)
        lm(i+3) = i+lhols(2,m)
       enddo
c      store in d-array      
       do i = 1,6
        d(kk+i) = lm(i)
       enddo
       kk = kk + 6
      enddo
c      
c     output material properties
c
                  write(iow,1002) (d(i),i=1,26),numel1,d(28)
      if(ior.lt.0)write(*  ,1002) (d(i),i=1,26),numel1,d(28)
1002  format(5x,'material data for cell model ',/,
     +  5x,'01: e11 fiber...........................',e12.3,/,
     +  5x,'02: e22 fiber...........................',e12.3,/,
     +  5x,'03: e33 fiber...........................',e12.3,/,
     +  5x,'04: v12 fiber...........................',f12.3,/,
     +  5x,'05: v13 fiber...........................',f12.3,/,
     +  5x,'06: v23 fiber...........................',f12.3,/,
     +  5x,'07: g12 fiber...........................',f12.3,/,
     +  5x,'08: g13 fiber...........................',f12.3,/,
     +  5x,'09: g23 fiber...........................',f12.3,/,
     +  5x,'10: alpha11 fiber.......................',f12.3,/,
     +  5x,'11: alpha22 fiber.......................',f12.3,/,
     +  5x,'12: alpha33 fiber.......................',f12.3,/,
     +  5x,'13: e matrix............................',e12.3,/,
     +  5x,'14: poissons ratio matrix...............',f12.3,/,
     +  5x,'15: shear modulus matrix................',f12.3,/,
     +  5x,'16: alpha matrix........................',f12.3,/,
     +  5x,'17: spring stiffness in 1 dir...........',e12.3,/,
     +  5x,'18: spring stiffness in 2 dir...........',e12.3,/,
     +  5x,'19: spring stiffness in 3 dir...........',e12.3,/,
     +  5x,'20: temperature change..................',f12.3,/,
     +  5x,'21: fiber volume fraction...............',f12.3,/,
     +  5x,'22: tensile strength-fiber..............',e12.3,/,
     +  5x,'23: compressive strength-fiber..........',e12.3,/,
     +  5x,'24: ultimate matrix tensile strength....',e12.3,/,
     +  5x,'25: ultimate matrix compressive strength',e12.3,/,
     +  5x,'26: ultimate matrix shear strength......',e12.3,/,
     +  5x,'27: number of cells.....................',i4/,
     +  5x,'28: fiber angle phi.....................',e12.3/)
c
      return
      end 
