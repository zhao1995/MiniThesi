      subroutine elmt31(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-----------------------------------------------------------------------
c
c.... 4 node nonlinear flat shell element (Membrane + Bathe/Dvorkin)
c.... with moderate rotations                                       
c
c
c     LÄUFT nicht unter INTEL????????? WW08/03/2010
c
c
c-----------------------------------------------------------------------
c     open: 
c     ################################################################
c     Material data only with respect to local base system(diagonals)
c     only rectangular elements are allowed!! 
c     ################################################################
c     Base Vector t_i only for postprocessing and error analysis
c.....  transform stresses compare elmt36
c.....  here T^T(!) 


c     ACTUAL: Igeo used for Base vectors under ityp=2


c
c       
c     test rayleigh damping  
c 
c     energy conserving algorithm only for membrane part tested, 
c     geht nicht mehr wg. damage!! 
c
c-----------------------------------------------------------------------
c.... history 
c....     06/04: calculate dme,dmb,db,ds if necessary for isw=4
c.....           removed for damage 03/07
c....     03/06: Shear correction factor                   
c....     03/06: external loads QLOA
c....     03/06: energy conserving time integration BETA,ENER
c....     01/07: print/plot stresses layerwise
c....     03/07: failure tsai-wu model included  
c....     03/07: degradation     model included  
c....     04/07: calculate shear stresses gamma=J^-1*gamma_conv
c-...            instead of ds = J^T-1*ds*J^-1   
c                                                                   
c-----------------------------------------------------------------------
c     Macros included:                                                
c      1=Input                                                        
c      3=TANG                                                         
c      4=STRE(Print)                                                  
c      5=LMAS/CMAS                                                    
c      6=REAC                                                         
c      8=STRE(Plot)                                                   
c     11=EXT                                                          
c     12=DAMP only consistent                                        
c     13=STR1 but here: Plot local base vectors!                     
c     22=QLOA INPUT: q_z, q_n
c-----------------------------------------------------------------------
c
c.... Input data   1. line                                            
c
c     d( 1) = nlay: no. of layers                                     
c     d( 2) = E_1                          E  in case of isotropy     
c     d( 3) = E_2                          E                          
c     d( 4) = v_12                         v                          
c     d( 5) = G_12 = G13                   G                          
c     d( 6) = G_23                         G                          
c     d( 7) =  --    calculated  v_21 (not used)                                 
c     d( 8) = rho=gamma/g   per volume                                
c     d( 9) = q_z do not use!! use QLOA                                                      
c     d(10) = q_n do not use!! use QLOA                                                    
c-----------------------------------------------------------------------
c
c     Input data 2. line                                              
c
c     d(11) = 0=lin,1=w,x**2,2=u,x**2+v,x**2+w,x**2                   
c     d(12) = 0/1 warping correction with offset coordinates          
c     d(13) = 0/1 shear correction factor                             
c     d(14) = relative stiffness on 6.local dof: suggest: (1.e-4)     
c     d(15) = alpha Rayleigh damping  alpha*M + beta*K                 
c     d(16) = beta  Rayleigh damping                                  
c-----------------------------------------------------------------------
c
c     Input data 3. line - basis for stress output, see pdirec3       
c     d(17) type                | d(18)| d(19)| d(20)                 
c     1  basic via diagonals    |      |      |                       
c     2  t_3(1)+given vector t_1|  t1x |  t1y |  t1z                  
c     3 - 17    see macro BASE                       
c                                                                     
c     control correct base system via plot,forc,1                     
c     plates may lead to a wrong orientation (det J<0)                
c     use in that case type 2                                         
c-----------------------------------------------------------------------
c
c     Input data 4. line - data for failure models d(21)-d(36)      
c     read failure values        d(21)-d(28)
c     set          values        d(29)-d(36)
c     d(21)= ityp (damage criterion)
c
c     Input data 5. line - data for degradation models d(37)-d(46)      
c     read degradation values    d(37)-d(44)
c         set          values    d(45)-d(46)
c     d(37)= ityp (damage criterion)
c
c-----------------------------------------------------------------------
c     set
c     d(47) = Total thickness of shell 
c     d(48) = ivor = number of elements in d-field before ply values  
c-----------------------------------------------------------------------
c
c     Input data 6.line                                               
c
c     Input data for each layer (on separate input line)              
c     phi       angle between global and local coordinate system      
c     zs....... coordinate to midpoint of layer                        
c     h........ thickness of layer                                    
c
c     angles with respect to local coordinate system via diagonals!!
c     makes only sense for rectangular elements!!
c
c-----------------------------------------------------------------------
c     displacements and rotations in global coordinate directions     
c                                                                     
c-----------------------------------------------------------------------
c.... stress plot           stress print (dir from igeo)              
c    1   n_11              1   n_11                                   
c    2   n_12              2   n_22  !!    reason: Main stresses       
c    3   n_22              3   n_12  !!                               
c    4   m_11              4   m_11                                    
c    5   m_12              5   m_22  !!                               
c    6   m_22              6   m_12  !!                               
c    7   q_13              7   q_13                                   
c    8   q_23              8   q_23                                   
c
c-----------------------------------------------------------------------
c
c.... energy conserving algorithm nop=5 in /ddata/
c.... vl|n     = vl1  =  ul-dul  
c.... vl|n+1/2 = vlm  =  ul -0.5*dul
c.... vl|n+1   = vl2  =  ul  
c
c-----------------------------------------------------------------------
c                                                                  
      USE bdata       
      USE cdata 
      USE cdat1 
      USE ddata 
      USE eldata 
      USE hdata 
      USE iofile 
      USE plslay 
      USE prlod 
      USE qload       
      USE strnam 

      implicit double precision(a-h,o-z)
      dimension d(*),s(nst,nst),p(nst),tl(*),
     +          gp(3),pgs(2,4),wgs(4),
     +          tr(3,3),ts(3,3),dummy(24),     
     +          xl(ndm,*),yl(3,4),shp(3,4),xjaci(2,2),deriv(2,4),
     +          gr1(10),grm(10),gr2(10),
     +          gr1w(4),grmw(4),gr2w(4),
     +          ul(ndf,*),vl1(6,4),vlm(6,4),vl2(6,4),
     +          bmi(3,3),bmj(3,3),bbi(3,3),bbj(3,3),bsi(2,3),bsj(2,3),
     +          bmtd(3,3),bbtd(3,3),bmbtd(3,3),bbmtd(3,3),bstd(3,2),
     +          dme(3,3),dmb(3,3),db(3,3),ds(2,2),
     +          sig1(8),sig2(8),eps(8),
     +          sij1(3,3),sij2(3,3),sij3(3,3),sij4(3,3),sij5(3,3),
     +          ri1(3),ri2(3),ri3(3),
     +          ql(3),qld(2)
      dimension h1(*),h2(*),h3(*)
c
c.... go to correct array processor
      go to(1,2,3,4,5,3,2,4,2,2,3,3,13,2,2,2,2,2,2,2,2,22) isw
      return
c
c.... mat.properties 1:nlay,E_1,E_2,v_12,G_12=G_13,G_23,rho,q_z,q_n
1     call dinput(d,9)
      d(10) = d(9)
      d(9)  = d(8)
      d(8)  = d(7)
      d(7)  = d(4)*d(3)/d(2) ! v_21 = v_12*E_2/E_1
                   write(iow,1001) (d(i),i=1,10)
      if(ior.lt.0) write(*  ,1001) (d(i),i=1,10)
1001  format(/,5x,'Material data for nonl. flat shell:',/
     + 5x,'Number of layers................',f8.0,/,
     + 5x,'Elastic modulus E_1 ........',g12.5,/,
     + 5x,'Elastic modulus E_2 ........',g12.5,/,
     + 5x,'Poisson ratio   v_12........',g12.5,/,
     + 5x,'Shear modulus   G_12 = G_13.',g12.5,/,
     + 5x,'Shear modulus   G_23 .......',g12.5,/,
     + 5x,'Poisson ratio   v_21 .......',g12.5,/,
     + 5x,'rho=gamma/g ................',g12.5,/,
     + 5x,'q_z load in z-direction.....',g12.5,/,
     + 5x,'q_n load in transv. dir.....',g12.5) 

c.... test of correct length of d-array
      nlay = d(1)
      d(48)= 48 !set
      ivor = d(48)
      nmax =  ivor + 3 * nlay
      if(nmax.gt.ndd) then
         write(*,1010) ndd,nmax
1010     format(1x,'darray is set to  ',i4,' values',/,
     1          1x,'darray needs      ',i4,' values')
         stop
      end if
c
c.... mat. properties 2: lin,warp,shear,c6dof,alpha,beta(Rayleigh damping)                
      call dinput(d(11),6)
      if(d(13).eq.0.d0) d(13)=5.d0/6.d0  ! Def. cappa
                   write(iow,1002) (d(i),i=11,16)
      if(ior.lt.0) write(*  ,1002) (d(i),i=11,16)
1002  format(5x,'Type of calculation:',/
     + 5x,'lin=0, nonlin w=1, nonlin uvw=2......',f8.0,/,
     + 5x,'warping correction 0/1...............',f8.0,/,
     + 5x,'shear   correction 0/1...........',g12.5,/,
     + 5x,'relative stiffness on 6.local dof',g12.5,/,
     + 5x,'alpha for Rayleigh damping       ',g12.5,/,
     + 5x,'beta  for Rayleigh damping       ',g12.5)

c.... mat. properties 3: type,val1,val2,val3
      call dinput(d(17),4)
                   write(iow,1003) (d(i),i=17,20)
      if(ior.lt.0) write(*  ,1003) (d(i),i=17,20)
1003  format(5x,'Type of base vector + add. values',/
     1 5x,'1 via diagonals, 2=t_3(1)+given vector t_1 with t1x,t1y,t1z'
     1 ,/,5x,'Others: see Macro DREC: ',f4.0,3f12.5)
c
c.... mat. properties 4: data for failure models iftyp,f1,f2,....
      call damage31(d,1)

c.... mat. properties 5: data for degradation models idtyp,d1,d2,....
      call damage31(d,2)

c.... mat. properties 6: phi,zs,h for each layer
                   write(iow,1005)
      if(ior.lt.0) write(*  ,1005)
1005     format(5x,'material data for layers',/,
     +   5x,'layer  angle phi-x1 coordinate z_s thickness  h')

c.... input of data for each layer (3 values)
      do ilay =1,nlay
        ia = (ilay - 1) * 3 + ivor
        ii = ia + 1
        call dinput(d(ii),3)
                     write(iow,1006) ilay,(d(i),i=ii,ii+2)
        if(ior.lt.0) write(*  ,1006) ilay,(d(i),i=ii,ii+2)
      end do
1006  format(5x,i5,3(2x,g12.5))
c
c...  thickness of shell
      hs = 0.0
      do ii = 1, nlay
        ia  = (ii - 1) * 3 + ivor
        hs  = hs + d(ia + 3)
      end do
      d(47) = hs  ! control
                   write(iow,1007) d(47)
      if(ior.lt.0) write(*  ,1007) d(47)
1007  format(5x,'Total shell thickness ............',g12.5)

c
c.... set values for h-array
      idtyp=d(37)

      ngs = 2
      if(idtyp.eq.0)then  ! dummy
        nh1=ngs*ngs*nlay*1  
      end if 

      if(idtyp.eq.1)then  
        nh1=ngs*ngs*nlay*1  
      end if 

2     return
c
c.... stiffness matrix and residual
3     nlay  = d(1)
      lin   = d(11)
      iwarp = d(12)
      call trans31(d,xl,yl,ul,vl1,vlm,vl2,tr,ndm,ndf,nen)

c...  Gauss-Points 
      ngs = 2
      lgs = ngs*ngs
      call gauss31(ngs,pgs,wgs)
      
c.....dvorkin/bathe values + shear correction 
      call bdval31(yl,ndm)
      call shearfac31(d,xl,cappa)

c.... loop Gauss points
      nn = 0 ! counter h-array 
      do 30 igs = 1,lgs
        iko = (nn+1)*nlay
        xsi = pgs(1,igs)
        eta = pgs(2,igs)

        call  sfr31(shp,xsi,eta)                       ! shape functions
        call jaco31(shp,djacb,xjaci,deriv,yl,n,nel,ndm)! jacobian 
        da = djacb*wgs(igs)

        call grad31(shp,gr2,gr2w,yl,vl2,ndf,iwarp)     ! displacement gradients
        call stra31(vl2,gr2,gr2w,eps,xjaci,xsi,eta,lin)! strains

        if(nop.eq.5) then 
          call grad31(shp,gr1,gr1w,yl,vl1,ndf,iwarp)
          call grad31(shp,grm,grmw,yl,vlm,ndf,iwarp)
        end if

        call mod31(dme,dmb,db,ds,cappa,d,eps,h1(1+iko),h2(1+iko),
     +             idam,lin)                            ! Mat.matrix    

        call stre31(dme,dmb,db,ds,sig2,eps)             ! stresses

        c05 = 1.d0
        if(nop.eq.5) then
          call stra31(vl1,gr1,gr1w,eps,xjaci,xsi,eta,lin)
          call stre31(dme,dmb,db,ds,sig1,eps)         ! das ist falsch, da dme,db.... nicht an richtiger Stelle!!!   
c....     stresses at t_n+1/2, store now at ...2 
          do i=1,8
            sig2(i) = 0.5d0*(sig1(i)+sig2(i))
          end do
          c05 = 0.5d0
        end if

c....  loop over node i
        do 31 ino =1,nel
          ia = (ino-1)*ndf
c.....    external load vector                               
          if(d(9).eq.0.0.and.d(10).eq.0.0) goto 36
          call pzero(ql,3)
          call qload31(qld,aqloa,d,numel,n,mqloa,propq,prop,isw)

c.....    external constant load qz in z-direction
          do i = 1,3
            ql(i) = tr(i,3)*qld(1)
          end do         
c.....    external constant load qn transversal to element
          ql(3) = ql(3) + qld(2)
c.....    add loads to load vector (prop/propq in qload!)
          do i = 1,3
            p(ia+i) = p(ia+i) + shp(3,ino)*ql(i)*da 
          end do
36        continue

c.....    b-matrices
          if(nop.eq.5) then 
            call bmat31(bmi,bbi,bsi,shp,deriv,ino,grm,xjaci,lin)
          else  
            call bmat31(bmi,bbi,bsi,shp,deriv,ino,gr2,xjaci,lin)
          end if

c.....    residual
          call mulp31(bmi,sig2(1),ri1,da,3,3)
          call mulp31(bbi,sig2(4),ri2,da,3,3)
          call mulp31(bsi,sig2(7),ri3,da,3,2) 
          p(ia+1) = p(ia+1) - ri1(1) 
          p(ia+2) = p(ia+2) - ri1(2) 
          p(ia+3) = p(ia+3) - ri1(3) - ri2(1) - ri3(1) 
          p(ia+4) = p(ia+4)          - ri2(2) - ri3(2) 
          p(ia+5) = p(ia+5)          - ri2(3) - ri3(3)
          if(isw.eq.6) goto 31
c.....    btd
          call muls31(bmi,dme,bmtd,da,3,3,3)
          call muls31(bbi, db,bbtd,da,3,3,3)
          call muls31(bsi, ds,bstd,da,3,2,2)  
c....     coupling terms          
          call muls31(bmi,dmb,bmbtd,da,3,3,3)
          call muls31(bbi,dmb,bbmtd,da,3,3,3) 

c.....    loop over node j 
          jnos = ino
          if(nop.eq.5) jnos=1 ! unsymmetry    
          do 32 jno = jnos,nel
            ja = (jno-1)*ndf
c.....      b-matrices
            call bmat31(bmj,bbj,bsj,shp,deriv,jno,gr2,xjaci,lin)
c....       stiffness matrix
c....       membrane + bending + shear terms
            call mulk31(bmtd,bmj,sij1,3,3,3)
            call mulk31(bbtd,bbj,sij2,3,3,3)
            call mulk31(bstd,bsj,sij5,3,2,3)
            do ii = 1,3
              do jj = 1,3
                s(ia+ii,  ja+jj  )=s(ia+ii,  ja+jj  ) + sij1(ii,jj)*c05
                s(ia+ii+2,ja+jj+2)=s(ia+ii+2,ja+jj+2) + sij2(ii,jj)*c05 
     +                                                + sij5(ii,jj)*c05
              end do
            end do     
c....       coupling terms          
            call mulk31(bmbtd,bbj,sij3,3,3,3)
            call mulk31(bbmtd,bmj,sij4,3,3,3)
            do ii = 1,3
              do jj = 1,3
                s(ia+ii  ,ja+jj+2)=s(ia+ii  ,ja+jj+2)+sij3(ii,jj)*c05
                s(ia+ii+2,ja+jj  )=s(ia+ii+2,ja+jj  )+sij4(ii,jj)*c05
              end do
            end do     
            if(lin.gt.0) then
c....         nonlinear stiffness matrix
              call ksig31(aksig,shp,sig2,ino,jno)
              s(ia+3,ja+3)   = s(ia+3,ja+3) + aksig*da
              if(lin.eq.2) then
                s(ia+1,ja+1) = s(ia+1,ja+1) + aksig*da
                s(ia+2,ja+2) = s(ia+2,ja+2) + aksig*da
              end if
            end if
32        continue
31      continue
30    continue
c.... stiffness on 6.local dof on diagonal K_66 = d(14)*K_ii(bending), suggestion: d(14)=10**-4  
      sk66 = 0.d0
      do ii = 1,4
        jj = (ii-1)*6
        do kk = 3,5
          sk66 = sk66+ s(jj+kk,jj+kk)    
        end do
      end do
      sk66=sk66/12.d0 ! average value, 4 nodes * 3 bending terms 
      s( 6, 6) = d(14)*sk66
      s(12,12) = d(14)*sk66
      s(18,18) = d(14)*sk66
      s(24,24) = d(14)*sk66

c.... warping: modify for offset projections (Taylor, Mafelap87)
      if(iwarp.eq.1) then
        i1 = 1
        if(yl(3,1).ne.0.0d0) then
          do i = 1,4
            p(i1+3) = p(i1+3) + yl(3,i)*p(i1+1)  
            p(i1+4) = p(i1+4) - yl(3,i)*p(i1)
            j1 = i1
            if(nop.eq.5) j1=1 
            do j = i,4
              call proj31(s(i1,j1),yl(3,i),yl(3,j),nst)
              j1 = j1 + ndf
            end do    
            i1 = i1 + ndf
          end do   
        end if
      end if 
c
      if(isw.eq.12) goto 12  ! damping matrix for Rayleigh-Damping C=alpha*M+beta*K
c
c.... symmetry+transformation 
      call transl31(s,p,tr,nst,ndf,nop)

      return
c
c.... print and plot stresses
c...  Gauss-Points
4     ngs = 2
      lgs = ngs*ngs
      call gauss31(ngs,pgs,wgs)
      nlay  = d(1)
      lin   = d(11)
      iwarp = d(12) 
c.... description of stresses  
      call plsn31(klay)
c.... local arbitrary basis tr 
      call trans31(d,xl,yl,ul,vl1,vlm,vl2,tr,ndm,ndf,nen)
c.....dvorkin/bathe values + shear factor
      call bdval31(yl,ndm)
      call shearfac31(d,xl,cappa)
c.... local basis ts with respect to igeo
cww      call transs31(d,xl,tr,ts,ndm)
c
      nn = 0 ! for history data
      do 40 igs = 1,lgs
        iko = (nn+1)*nlay
        xsi = pgs(1,igs)
        eta = pgs(2,igs)
        call  sfr31(shp,xsi,eta)
        call jaco31(shp,djacb,xjaci,deriv,yl,n,nel,ndm)
        da = djacb*wgs(igs)
        call grad31(shp,gr2,gr2w,yl,vl2,ndf,iwarp)
        call stra31(vl2,gr2,gr2w,eps,xjaci,xsi,eta,lin)
        if(klay.eq.0) then  ! stress resultants
          call mod31(dme,dmb,db,ds,cappa,d,eps,h1(1+iko),h2(1+iko),
     +               idam)
          call stre31(dme,dmb,db,ds,sig2,eps)
        else ! stresses at layer
          call strlay31(d,sig2,eps,cappa,klay,h1(1+iko),h2(1+iko)
     +                 ,idam)
        end if 

c....   transform stresses with respect to igeo tr->ts
cww        call transsig31(tr,ts,sig2)

        if(isw.eq.4) then
c....     Print stresses
          do idm = 1,3
            gp(idm) = 0.0
            do ino = 1,nel
              gp(idm) = gp(idm) + xl(idm,ino) * shp(3,ino)
            end do
          end do
          mct = mct - 1
          if(mct.gt.0) go to 41
          if(klay.eq.0) then  ! resultants
                         write(iow,4000)
            if(ior.lt.0) write(  *,4000)
4000       format(2x,'STRESS RESULTANTS Nonlinear flat shell element',/,  
     *    1x,' El ',' Mat ',' X-Coord ',' Y-coord ',' Z-coord ',
     *    '      ',   '    N_11     ','    N_22     ','    N_12     ',
     *                '    M_11     ','    M_22     ','    M_12     ',
     *                '    Q_13     ','    Q_23     ','damage',' warp ')

          else  ! stresses
                         write(iow,4001)
            if(ior.lt.0) write(  *,4001)
4001       format(2x,'STRESSES Nonlinear flat shell element',/,  
     *    1x,' El ',' Mat ',' X-Coord ',' Y-coord ',' Z-coord ',
     *    ' Layer',   '    S_11     ','    S_22     ','    S_12     ',
     *                '             ','             ','             ',
     *                '    T_13     ','    T_23     ','damage',' warp ')
          end if
          mct = 50
41        warp = 0.d0
          if(yl(3,1).ne.0.d0) warp = yl(3,1)          
                       write(iow,4002)n,ma,gp,klay,sig2,idam,warp 
          if(ior.lt.0) write(  *,4002)n,ma,gp,klay,sig2,idam,warp 
4002      format(1x,i4,i4,3f9.3,i6,1x,8g13.5,i6,g13.5)
c
        else if(isw.eq.8) then
c....     calculate nodal stresses
          if(iplma(ma).eq.0)  return ! only if MATN
          call plot31(ix,strea,strea(1+numnp),shp,sig2,idam,da,numnp)
        end if
   40 continue
      return
c
c.... mass matrix   
5     call trans31(d,xl,yl,ul,vl1,vlm,vl2,tr,ndm,ndf,nen) ! necessary?->for da!
c...  Gauss-Points 
      ngs = 2
      lgs = ngs*ngs
      call gauss31(ngs,pgs,wgs)
c
c...  gamma, thickness of shell
      rho = d(8)
      hs  = d(47)
c      
      do 50 igs = 1,lgs
        xsi = pgs(1,igs)
        eta = pgs(2,igs)
        call  sfr31(shp,xsi,eta)
        call jaco31(shp,djacb,xjaci,deriv,yl,n,nel,ndm)
        da  = djacb*wgs(igs)
        do 51 ino =1,nel
c
c......   lmas
          iq  = (ino-1)*ndf
          pma = shp(3,ino)*rho*da*hs
          pmi = pma*hs*hs/12.d0
          p(iq+1) = p(iq+1) + pma
          p(iq+2) = p(iq+2) + pma
          p(iq+3) = p(iq+3) + pma
c......   lmas incl. rot terms for all terms!! (6.dof=?)
          p(iq+4) = p(iq+4) + pmi
          p(iq+5) = p(iq+5) + pmi
          p(iq+6) = p(iq+6) + pmi ! ??
c
c......   cmas including rot. terms (6.dof=?)
          do 52 jno = ino,nel
            iq = (ino-1)*ndf
            jq = (jno-1)*ndf
            s(iq+1,jq+1)=s(iq+1,jq+1)+pma*shp(3,jno)
            s(iq+2,jq+2)=s(iq+2,jq+2)+pma*shp(3,jno)
            s(iq+3,jq+3)=s(iq+3,jq+3)+pma*shp(3,jno)
            s(iq+4,jq+4)=s(iq+4,jq+4)+pmi*shp(3,jno)
            s(iq+5,jq+5)=s(iq+5,jq+5)+pmi*shp(3,jno)
            s(iq+6,jq+6)=s(iq+6,jq+6)+pmi*shp(3,jno) ! ??
52        continue
51      continue
50    continue
c.... symmetry, transformation not allowed for p(negative masses)->dummy
      call pzero(dummy,24)      
      call transl31(s,dummy,tr,nst,ndf,nop)
      return
c
c.... damping matrix for Rayleigh-Damping C=alpha*M+beta*K  only consistent   
12    alpha=d(15)
      beta =d(16)
c.... term beta*K
      do i = 1,nst
        do j = i,nst
          s(i,j) = beta*s(i,j)
        end do
      end do
c.... term alpha*M
c...  thickness of shell
      nlay = d(1)
      ivor = d(48)
      hs = 0.0
      do ii = 1, nlay
        ia  = (ii - 1) * 3 + ivor
        hs  = hs + d(ia + 3)
      end do
c.... gamma      
      rho = d(8)
c      
      do 120 igs = 1,lgs
        xsi = pgs(1,igs)
        eta = pgs(2,igs)
        call  sfr31(shp,xsi,eta)
        call jaco31(shp,djacb,xjaci,deriv,yl,n,nel,ndm)
        da  = djacb*wgs(igs)
        do 121 ino =1,nel
          iq  = (ino-1)*ndf
          pma = alpha*shp(3,ino)*rho*da*hs
          pmi = pma*hs*hs/12.d0
cwwc......   lmas incl. rot terms for all terms!! (6.dof=?) on Diag[K]!!
cww       s(iq+1,iq+1) = s(iq+1,iq+1) + pma
cww       s(iq+2,iq+2) = s(iq+2,iq+2) + pma
cww       s(iq+3,iq+3) = s(iq+3,iq+3) + pma
cww       s(iq+4,iq+4) = s(iq+4,iq+4) + pmi
cww       s(iq+5,iq+5) = s(iq+5,iq+5) + pmi
cww       s(iq+6,iq+6) = s(iq+6,iq+6) + pmi
c
c......   cmas including rot. terms (6.dof=?)
          do 122 jno = ino,nel
            iq = (ino-1)*ndf
            jq = (jno-1)*ndf
            s(iq+1,jq+1)=s(iq+1,jq+1)+pma*shp(3,jno)
            s(iq+2,jq+2)=s(iq+2,jq+2)+pma*shp(3,jno)
            s(iq+3,jq+3)=s(iq+3,jq+3)+pma*shp(3,jno)
            s(iq+4,jq+4)=s(iq+4,jq+4)+pmi*shp(3,jno)
            s(iq+5,jq+5)=s(iq+5,jq+5)+pmi*shp(3,jno)
            s(iq+6,jq+6)=s(iq+6,jq+6)+pmi*shp(3,jno)
122       continue
121     continue
120   continue
c
c.... symmetry, transformation not allowed for p(negative masses)->dummy ???
      call transl31(s,dummy,tr,nst,ndf,nop)
      return
c
c.... plot local base system at center of element via FORC,1
13    continue
c     in case of matn
      if(iplma(ma).eq.0) return
      call trans31(d,xl,yl,ul,vl1,vlm,vl2,tr,ndm,ndf,nen)
cww      call trans31(xl,yl,ul,vl2,tr,ndm,ndf)
      call pppcol(2)
      call plloco31(tr,xl)
      return
c
c.... external loads QLOA
22    lin   = d(11)
      iwarp = d(12)
      call trans31(d,xl,yl,ul,vl1,vlm,vl2,tr,ndm,ndf,nen)
c
c...  Gauss-Points
      ngs = 2
      lgs = ngs*ngs
      call gauss31(ngs,pgs,wgs)
c
c.....dvorkin/bathe values
c      call bdval31(yl,ndm)???

c.... external load vector (only conservative loads, no temperature) 
      do 220 igs = 1,lgs
        xsi = pgs(1,igs)
        eta = pgs(2,igs)
        call  sfr31(shp,xsi,eta)
        call jaco31(shp,djacb,xjaci,deriv,yl,n,nel,ndm)
        da = djacb*wgs(igs)

        do 221 ino =1,nel
          ia = (ino-1)*ndf
c.....    external load vector                               
          call pzero(ql,3)
          call qload31(qld,aqloa,d,numel,n,mqloa,propq,prop,isw)

c.....    external constant load qz in z-direction
          do i = 1,3
            ql(i) = tr(i,3)*qld(1)
          end do         
c.....    external constant load qn transversal to element
          ql(3) = ql(3) + qld(2)
c.....    add loads to load vector
          do i = 1,3
            p(ia+i) = p(ia+i) + shp(3,ino)*ql(i)*da
          end do
221     continue
220   continue

c.... warping: modify for offset projections (Taylor, Mafelap87)
      if(iwarp.eq.1) then
        i1 = 1
        if(yl(3,1).ne.0.0d0) then
          do i = 1,4
            p(i1+3) = p(i1+3) + yl(3,i)*p(i1+1)  
            p(i1+4) = p(i1+4) - yl(3,i)*p(i1)
            i1 = i1 + ndf
          end do   
        end if
      end if 
c
c.... symmetry+transformation 
      call transl31(s,p,tr,nst,ndf,nop)


      return

      end
c
      subroutine jaco31(shp,djacb,xjaci,deriv,xl,n,nel,ndm)
c----------------------------------------------------------------------
c.... Jacobi-Matrix and  cartesian derivatives of shape functions
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension cartd(2,4),shp(3,4),xl(ndm,nel),xjaci(2,2),xjacm(2,2),
     *          deriv(2,4)

c.... Jacobi-matrix xjacm
      call pzero(xjacm,4)
      do 4 idm = 1,2
         do 4 jdm = 1,2
            do 4 ino = 1,nel
    4          xjacm(idm,jdm)=xjacm(idm,jdm)+shp(idm,ino)*xl(jdm,ino)

c.... Determinant and Inverse of Jacobi-Matrix
      djacb = xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1)
      if (djacb) 6,6,8
    6 write(*,600) n
      stop
    8 xjaci(1,1) =  xjacm(2,2)/djacb
      xjaci(2,2) =  xjacm(1,1)/djacb
      xjaci(1,2) = -xjacm(1,2)/djacb
      xjaci(2,1) = -xjacm(2,1)/djacb

c.... cartesian derivatives
      do 10 idm = 1,2
         do 10 ino = 1,nel
            cartd(idm,ino) = 0.0
            do 10 jdm = 1,2
   10       cartd(idm,ino)=cartd(idm,ino)+xjaci(idm,jdm)*shp(jdm,ino)

c.... save convective derivatives
      do 12 i = 1,2
        do 12 j = 1,4
12        deriv(i,j) = shp(i,j)

c.... copy in shp
      do 14 i = 1,2
        do 14 j = 1,4
14        shp(i,j) = cartd(i,j)
  600 format(/,' Det J negativ in elmt31s,x',/,
     +      'wrong coordinates or problem with transformation type',/,
     +      'element number',i5)
      return
      end
c
      subroutine bmat31(bm,bb,bs,shp,deriv,i,gr,xjaci,lin)
c----------------------------------------------------------------------
c.... B-Matrices for flat-Shell membran+bending+shear(B/D)
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      common /bdval1/ xxsim(4),yxsim(4),xetal(4),yetal(4)
      dimension bm(3,3),bb(3,3),bs(2,3),bsc(2,3),shp(3,4),gr(10),
     +          xsii(4),etai(4),deriv(2,4),xjaci(2,2)
      data xsii /-1.d0, 1.d0,1.d0,-1.d0/
      data etai /-1.d0,-1.d0,1.d0, 1.d0/
      dndx = shp(1,i)
      dndy = shp(2,i)
      dnds = deriv(1,i)
      dndt = deriv(2,i)
c
      if(lin.ne.0) then
        ukx = gr(1)
        vkx = gr(2)
        wkx = gr(3)
        uky = gr(6)
        vky = gr(7)
        wky = gr(8)
      end if
c       
c.... form bm
      call pzero(bm,9)
      bm(1,1) =  dndx
      bm(2,2) =  dndy
      bm(3,1) =  dndy
      bm(3,2) =  dndx

c.... nonlinear part
      if(lin.eq.2) then
        bm(1,1) =  bm(1,1) + ukx*dndx
        bm(2,1) =            uky*dndy
        bm(3,1) =  bm(3,1) + uky*dndx + ukx*dndy
c                          
        bm(1,2) =            vkx*dndx
        bm(2,2) =  bm(2,2) + vky*dndy
        bm(3,2) =  bm(3,2) + vky*dndx + vkx*dndy
c
      end if
      if(lin.gt.0) then
        bm(1,3) =            wkx*dndx
        bm(2,3) =            wky*dndy
        bm(3,3) =            wky*dndx + wkx*dndy
      end if
c       
c.... form bb
      call pzero(bb,9)
      bb(1,3) =  dndx
      bb(2,2) = -dndy
      bb(3,2) = -dndx
      bb(3,3) =  dndy
c
c.... form bs convective
      bsc(1,1) =  dnds                    
      bsc(1,2) = -dnds*xsii(i)*yxsim(i)
      bsc(1,3) =  dnds*xsii(i)*xxsim(i)
      bsc(2,1) =  dndt                 
      bsc(2,2) = -dndt*etai(i)*yetal(i)
      bsc(2,3) =  dndt*etai(i)*xetal(i)
c.... form bs cartesian
      call pzero(bs,6)
      do ii=1,2
        do k=1,3
          do j=1,2    
            bs(ii,k)=bs(ii,k)+xjaci(ii,j)*bsc(j,k)
          end do
        end do  
      end do 
c
      return
      end
c
      subroutine ksig31(aksig,shp,sig,i,j)
c-----------------------------------------------------------------------
c.... Initial stress  matrix ksigma = g T * sig * g
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension shp(3,4),sig(8)
      dndxi = shp(1,i)
      dndyi = shp(2,i)
      dndxj = shp(1,j)
      dndyj = shp(2,j)
      aksig = dndxi * sig(1) * dndxj + dndxi * sig(3) * dndyj +
     +        dndyi * sig(2) * dndyj + dndyi * sig(3) * dndxj
      return
      end
c
      subroutine grad31(shp,gr,grw,yl,vl,ndf,iwarp)
c----------------------------------------------------------------------
c.... displacement gradients for flat shell
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension shp(3,4),gr(10),grw(4),vl(ndf,*),yl(3,4)
      call pzero(gr,10)
      call pzero(grw,4)
      do ino = 1,4
        dndx = shp(1,ino)
        dndy = shp(2,ino)
        do idf = 1,5
          idfn = 5 + idf
          const = vl(idf,ino)
          gr(idf)  = gr(idf) + dndx * const
          gr(idfn) = gr(idfn)+ dndy * const
        end do

        if(iwarp.eq.1) then
          z = yl(3,ino)
          if(z.ne.0.d0) then ! warping
            grw(1) = grw(1) + z*dndx*vl(4,ino)
            grw(2) = grw(2) + z*dndx*vl(5,ino)
            grw(3) = grw(3) + z*dndy*vl(4,ino)
            grw(4) = grw(4) + z*dndy*vl(5,ino)
          end if
        end if

      end do
      return
      end
c
      subroutine mod31(dme,dmb,db,ds,cappa,d,eps,h1,h2,idam)
c----------------------------------------------------------------------
c     Elasticity matrix
c     with failure criterion and degradation model 
c....     if failure criterion was fullfilled earlier: E1=E2, E2=0.1*E2

c----------------------------------------------------------------------
      USE iofile 
      implicit double precision(a-h,o-z)
      dimension dme(3,3),dmb(3,3),db(3,3),ds(2,2),d(*),
     +  eps(8),em(3),eb(3),es(2),h1(*),h2(*)
      call pzero(dme,9)
      call pzero(dmb,9)
      call pzero(db,9)
      call pzero(ds,4)

c.... global strains   
c.... Membrane
      em(1) = eps(1)
      em(2) = eps(2)
      em(3) = eps(3)

c.... Bending
      eb(1) = eps(4)
      eb(2) = eps(5)
      eb(3) = eps(6)

c.... Shear
      es(1) = eps(7)
      es(2) = eps(8)

c.... layer data
      nlay = d(1)
      ivor = d(48)
      pi = 4.0d0 * datan(1.0d0)
c.... material parameter input
      e1  = d(2)
      e2  = d(3)
      v12 = d(4)
      g12 = d(5)
      g23 = d(6)
c.... damage parameter input
      iftyp  = d(21)
      idtyp  = d(37)
      idam   = 0
c
c...  loop over all layers
      do 10 ii = 1, nlay
c...    data for actual layer
        ia  = (ii - 1) * 3 + ivor
        phi = d(ia + 1) / 180.0d0 * pi
        zs  = d(ia + 2)
        h   = d(ia + 3)
        if ( h.eq.0.0) goto 10
        zs2 = zs  * zs
        h3  = h * h * h
        c  = dcos(phi)
        c2 = c * c
        c3 = c * c * c
        c4 = c**4
        s  = dsin(phi)
        s2 = s * s
        s3 = s * s * s
        s4 = s**4
        
c....   local elasticity matrix of layer including damage
        call elast31(h11l,h22l,h12l,h33l,q11l,q22l,e1,e2,v12,g12,g23,
     +     cappa,em,eb,es,zs,d,h1(ii),h2(ii),iftyp,idtyp,c2,s2,c,s,idam)

c....   calculate global material parameters (x-y)
        h11g = c4*h11l + 2*c2*s2*(h12l+2*h33l) + s4*h22l
        h12g = c2*s2*(h11l+h22l-4*h33l) + (c4+s4)*h12l
        h13g = c3*s*(h11l-h12l-2*h33l) - s3*c*(h22l-h12l-2*h33l)
        h22g = s4*h11l + 2*c2*s2*(h12l+2*h33l) + c4*h22l
        h23g = s3*c*(h11l-h12l-2*h33l) - c3*s*(h22l-h12l-2*h33l)
        h33g = c2*s2*(h11l+h22l-2*h12l-2*h33l) + (c4+s4)*h33l
        q11g = c2*q11l + s2*q22l
        q12g = c*s*(q11l-q22l)
        q22g = s2*q11l + c2*q22l
c
c....   store into Elasticity matrix
c....   Membrane part
        dme(1,1) = dme(1,1) + h11g * h
        dme(1,2) = dme(1,2) + h12g * h
        dme(1,3) = dme(1,3) + h13g * h
        dme(2,2) = dme(2,2) + h22g * h
        dme(2,3) = dme(2,3) + h23g * h
        dme(3,3) = dme(3,3) + h33g * h

c....   Coupling Membrane-Bending part
        hzs = h * zs
        dmb(1,1) = dmb(1,1) + h11g * hzs
        dmb(1,2) = dmb(1,2) + h12g * hzs
        dmb(1,3) = dmb(1,3) + h13g * hzs
        dmb(2,2) = dmb(2,2) + h22g * hzs
        dmb(2,3) = dmb(2,3) + h23g * hzs
        dmb(3,3) = dmb(3,3) + h33g * hzs

c....   Bending part
        hzs2 =  h3 / 12.0d0 + h * zs2
        db (1,1) = db (1,1) + h11g * hzs2
        db (1,2) = db (1,2) + h12g * hzs2
        db (1,3) = db (1,3) + h13g * hzs2
        db (2,2) = db (2,2) + h22g * hzs2
        db (2,3) = db (2,3) + h23g * hzs2
        db (3,3) = db (3,3) + h33g * hzs2

c....   Shear part linear
        ds(1,1) = ds(1,1) + q11g * h
        ds(1,2) = ds(1,2) + q12g * h
        ds(2,2) = ds(2,2) + q22g * h
   10 continue

c.... symmetrize elasticity matrices
      dme(2,1) = dme(1,2)
      dme(3,1) = dme(1,3)
      dme(3,2) = dme(2,3)
c
      dmb(2,1) = dmb(1,2)
      dmb(3,1) = dmb(1,3)
      dmb(3,2) = dmb(2,3)
c
      db(2,1) = db(1,2)
      db(3,1) = db(1,3)
      db(3,2) = db(2,3)
c
      ds(2,1) = ds(1,2)
c
      return
      end
c
      subroutine mulp31(a,b,c,da,ni,nj)
c----------------------------------------------------------------------
c.... Calculate C = A T * B * da  with A(nj,ni),B(nj)  
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension a(nj,ni),b(nj),c(ni)
      call pzero(c,ni)
      do i = 1,ni
        do j = 1,nj
          c(i) = c(i)+a(j,i)*b(j)*da
        end do
      end do
      return
      end
c
      subroutine muls31(a,b,c,da,ni,nj,nk)
c----------------------------------------------------------------------
c.... Calculate C = A T * B * da  with A(nj,ni),B(nj,nk)  
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension a(nj,ni),b(nj,nk),c(ni,nk)
      call pzero(c,ni*nk)
      do i = 1,ni
        do k = 1,nk
          do j = 1,nj
            c(i,k) = c(i,k)+a(j,i)*b(j,k)*da
          end do
        end do
      end do     
      return
      end
c
      subroutine mulk31(a,b,c,ni,nj,nk)
c----------------------------------------------------------------------
c.... Calculate C = A * B   with A(ni,nj),B(nj,nk)  
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension a(ni,nj),b(nj,nk),c(ni,nk)
      call pzero(c,ni*nk)
      do i = 1,ni
        do k = 1,nk
          do j = 1,nj
            c(i,k) = c(i,k)+a(i,j)*b(j,k)
          end do
        end do
      end do     
      return
      end
c
      subroutine stra31(vl,gr,grw,eps,xjaci,xsi,eta,lin)
c----------------------------------------------------------------------
c.... Calculate strains for flat shell
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      common /bdval1/ xxsim(4),yxsim(4),xetal(4),yetal(4)
      dimension vl(6,*),gr(10),grw(4),xjaci(2,2),
     1          em(3),eb(3),es(2),eps(8)

c.... membrane
      em(1) = gr(1)
      em(2) = gr(7)
      em(3) = gr(2) + gr(6)

c.... add terms from warping       
      em(1) = em(1) - grw(2)
      em(2) = em(2) + grw(3)
      em(3) = em(3) - grw(4) + grw(1)

      if(lin.gt.0) then
        em(1) = em(1) + 0.5d0*gr(3)*gr(3)
        em(2) = em(2) + 0.5d0*gr(8)*gr(8)
        em(3) = em(3) +       gr(3)*gr(8)
      end if

      if(lin.eq.2) then
        em(1) = em(1) + 0.5d0*gr(1)*gr(1) + 0.5d0*gr(2)*gr(2)
        em(2) = em(2) + 0.5d0*gr(6)*gr(6) + 0.5d0*gr(7)*gr(7)
        em(3) = em(3) +       gr(1)*gr(6) +       gr(2)*gr(7)
      end if

c.... Bending
      eb(1) =  gr(5)
      eb(2) = -gr(9)
      eb(3) = -gr(4) + gr(10)

c.....compatible shear strains in M,L
      esxsiB = (xxsim(1)*(vl(5,1) + vl(5,2)) 
     +       -  yxsim(1)*(vl(4,1) + vl(4,2)) + vl(3,2)-vl(3,1) )/2.d0  
      esxsiD = (xxsim(3)*(vl(5,3) + vl(5,4))
     +       -  yxsim(3)*(vl(4,3) + vl(4,4)) + vl(3,3)-vl(3,4) )/2.d0  
      esetaA = (xetal(1)*(vl(5,1) + vl(5,4))
     +       -  yetal(1)*(vl(4,1) + vl(4,4)) + vl(3,4)-vl(3,1) )/2.d0  
      esetaC = (xetal(2)*(vl(5,2) + vl(5,3))
     +       -  yetal(2)*(vl(4,2) + vl(4,3)) + vl(3,3)-vl(3,2) )/2.d0  

c.... convective shear strains Es
      esxsi = 0.5d0*( (1.d0-eta)*esxsiB + (1.d0+eta)*esxsiD )
      eseta = 0.5d0*( (1.d0-xsi)*esetaA + (1.d0+xsi)*esetaC )

c.....cartesian shear strains Ex=J-1*Es  
      es(1) = xjaci(1,1)*esxsi + xjaci(1,2)*eseta
      es(2) = xjaci(2,1)*esxsi + xjaci(2,2)*eseta

c.... store strains
      do i = 1,3
         eps(i)   = em(i)
         eps(i+3) = eb(i)
      end do
      eps(7) = es(1)   
      eps(8) = es(2)   
        
      return
      end
c

      subroutine stre31(dme,dmb,db,ds,sig,eps)
c----------------------------------------------------------------------
c.... Calculate stress resultants for flat shell
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension dme(3,3),dmb(3,3),db(3,3),ds(2,2),sig(8),eps(8),
     1          em(3),eb(3),es(2)
      call pzero(sig,8)

c.... membrane
      em(1) = eps(1)
      em(2) = eps(2)
      em(3) = eps(3)

c.... Bending
      eb(1) = eps(4)
      eb(2) = eps(5)
      eb(3) = eps(6) 

c.....shear strains  
      es(1) = eps(7)
      es(2) = eps(8)

c.... Normal forces
      sig(1)  = dme(1,1)*em(1)+dme(1,2)*em(2)+dme(1,3)*em(3)
      sig(2)  = dme(2,1)*em(1)+dme(2,2)*em(2)+dme(2,3)*em(3)
      sig(3)  = dme(3,1)*em(1)+dme(3,2)*em(2)+dme(3,3)*em(3)

c.... Bending moments
      sig(4) = db(1,1)*eb(1)+ db(1,2)*eb(2)+ db(1,3)*eb(3)
      sig(5) = db(2,1)*eb(1)+ db(2,2)*eb(2)+ db(2,3)*eb(3)
      sig(6) = db(3,1)*eb(1)+ db(3,2)*eb(2)+ db(3,3)*eb(3)

c.... coupling terms 
c.... Normal forces
      sig(1) = sig(1) + dmb(1,1)*eb(1)+dmb(1,2)*eb(2)+dmb(1,3)*eb(3)
      sig(2) = sig(2) + dmb(2,1)*eb(1)+dmb(2,2)*eb(2)+dmb(2,3)*eb(3)
      sig(3) = sig(3) + dmb(3,1)*eb(1)+dmb(3,2)*eb(2)+dmb(3,3)*eb(3)

c.... Bending moments
      sig(4) = sig(4) + dmb(1,1)*em(1)+dmb(1,2)*em(2)+dmb(1,3)*em(3)
      sig(5) = sig(5) + dmb(2,1)*em(1)+dmb(2,2)*em(2)+dmb(2,3)*em(3)
      sig(6) = sig(6) + dmb(3,1)*em(1)+dmb(3,2)*em(2)+dmb(3,3)*em(3)

c.....Shear forces   
      sig(7)= ds(1,1)*es(1) + ds(1,2)*es(2) 
      sig(8)= ds(2,1)*es(1) + ds(2,2)*es(2)

      return
      end
c

      subroutine strlay31(d,sig,eps,cappa,klay,h1,h2,idam)
c-----------------------------------------------------------------------
c     Calculation of stresses in global directions at center of layer  
c     klay with failure criterion and degradation model 
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      dimension d(*),sig(8),eps(8),em(3),eb(3),es(2),eg(3),h1(*),h2(*)
      call pzero(sig,8) 

c.... global strains   
c.... Membrane
      em(1) = eps(1)
      em(2) = eps(2)
      em(3) = eps(3)

c.... Bending
      eb(1) = eps(4)
      eb(2) = eps(5)
      eb(3) = eps(6)

c.... Shear
      es(1) = eps(7)
      es(2) = eps(8)

      nlay = d(1)
      ivor = d(48)
      pi = 4.0d0 * datan(1.0d0)
c.... material parameter 
      e1  = d(2)
      e2  = d(3)
      v12 = d(4)
      g12 = d(5)
      g23 = d(6)
c.... damage parameter input
      iftyp  = d(21)
      idtyp  = d(37)
      idam   = 0
c
c...  loop over all layers
      do 10 ii = 1, nlay
        if(ii.ne.klay) goto 10

c...    data for actual layer
        ia  = (ii - 1) * 3 + ivor
        phi = d(ia + 1) / 180.0d0 * pi
        zs  = d(ia + 2)
        h   = d(ia + 3)
        if ( h.eq.0.0) goto 10
      
        c  = dcos(phi)
        c2 = c * c
        c3 = c * c * c
        c4 = c**4
        s  = dsin(phi)
        s2 = s * s
        s3 = s * s * s
        s4 = s**4

c....   local elasticity matrix of layer including damage
        call elast31(h11l,h22l,h12l,h33l,q11l,q22l,e1,e2,v12,g12,g23,
     +  cappa,em,eb,es,zs,d,h1(ii),h2(ii),iftyp,idtyp,c2,s2,c,s,idam)

c....   global material parameters (x-y)
        h11g = c4*h11l + 2*c2*s2*(h12l+2*h33l) + s4*h22l
        h12g = c2*s2*(h11l+h22l-4*h33l) + (c4+s4)*h12l
        h13g = c3*s*(h11l-h12l-2*h33l) - s3*c*(h22l-h12l-2*h33l)
        h22g = s4*h11l + 2*c2*s2*(h12l+2*h33l) + c4*h22l
        h23g = s3*c*(h11l-h12l-2*h33l) - c3*s*(h22l-h12l-2*h33l)
        h33g = c2*s2*(h11l+h22l-2*h12l-2*h33l) + (c4+s4)*h33l

        q11g = c2*q11l + s2*q22l
        q12g = c*s*(q11l-q22l)
        q22g = s2*q11l + c2*q22l

c....   global strains at center of layer
        eg(1) = em(1) + zs*eb(1)
        eg(2) = em(2) + zs*eb(2)
        eg(3) = em(3) + zs*eb(3)

c....   global stresses 
        sig(1) =  h11g*eg(1) + h12g*eg(2) + h13g*eg(3)
        sig(2) =  h12g*eg(1) + h22g*eg(2) + h23g*eg(3)
        sig(3) =  h13g*eg(1) + h23g*eg(2) + h33g*eg(3)
        sig(7) =  q11g*es(1) + q12g*es(2)
        sig(8) =  q12g*es(1) + q22g*es(2)
	
10    continue
      return
      end 

c
      subroutine plsn31 (i)
c-----------------------------------------------------------------------
c.... define name for stresses in plot-output 
c     i = 0   : stress resultants
c         else: stresses and history values at defined layer
c-----------------------------------------------------------------------
      USE strnam 
      implicit  double precision (a-h,o-z)

      istv = -9

      do is =1,25
        strsus(is) = '               '
      end do
      
      if (i.eq.0) then
        strsus( 1) =  '  N-FORCE N_11 '
        strsus( 2) =  '  N-FORCE N_12 '
        strsus( 3) =  '  N-FORCE N_22 '
        strsus( 4) =  '  MOMENT  M_11 '
        strsus( 5) =  '  MOMENT  M_12 '
        strsus( 6) =  '  MOMENT  M_22 '
        strsus( 7) =  '  Q-FORCE Q_13 '
        strsus( 8) =  '  Q-FORCE Q_23 '
        strsus( 9) =  'Damage Par(0-1)'
      else           
        strsus( 1) = '  Sigma_11     '
        strsus( 2) = '  Sigma_22     '
        strsus( 3) = '  Sigma_12     '
        strsus( 7) = '  Sigma_13     '
        strsus( 8) = '  Sigma_23     '
        strsus( 9) = 'Damage Par(0-1)'
      end if
      return
      end
c
      subroutine plot31(ix,dt,st,shp,sig,idam,da,numnp)
c----------------------------------------------------------------------
c.... Plot stress resultants
c.... N11(1),N12(3),N22(2),M11(4),M12(6),M22(5),Q13(7),Q23(8),
c     Dam(9)= damage somewhere in layers (0-1), no details    
c
c     or   stresses
c.... S11(1),S12(3),S22(2),                    ,T13(7),T23(8),
c     Dam(9) = damage in layer (0-1) 
c.... 
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension dt(numnp),st(numnp,*),shp(3,4),sig(8),ix(*)
      do 10 j = 1,4
        xsji = da*shp(3,j)
        ii = ix(j)
        if (ii.gt.0) then
           dt(ii) = dt(ii) + xsji
c....      stress/stress resultants
             st(ii,1) = st(ii,1)+sig(1)*xsji
             st(ii,2) = st(ii,2)+sig(3)*xsji
             st(ii,3) = st(ii,3)+sig(2)*xsji
             st(ii,4) = st(ii,4)+sig(4)*xsji
             st(ii,5) = st(ii,5)+sig(6)*xsji
             st(ii,6) = st(ii,6)+sig(5)*xsji
             st(ii,7) = st(ii,7)+sig(7)*xsji
             st(ii,8) = st(ii,8)+sig(8)*xsji
             st(ii,9) = st(ii,9)+idam  *xsji

        end if
10    continue
      return
      end
c
      subroutine sfr31(shp,s,t)
c----------------------------------------------------------------------
c.... Shape functions and derivatives for bilinear 2-d elements
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension shp(3,4)
      d1 = 1.0d0
      d4 = 1.d0/4.d0
      st = s * t

c.... shape functions for 4-node element
         shp(3,1) = (d1-t-s+st) * d4
         shp(3,2) = (d1-t+s-st) * d4
         shp(3,3) = (d1+t+s+st) * d4
         shp(3,4) = (d1+t-s-st) * d4

c.... shp functions derivatives
         shp(1,1) = (-d1+t) * d4
         shp(1,2) = (+d1-t) * d4
         shp(1,3) = (+d1+t) * d4
         shp(1,4) = (-d1-t) * d4
         shp(2,1) = (-d1+s) * d4
         shp(2,2) = (-d1-s) * d4
         shp(2,3) = (+d1+s) * d4
         shp(2,4) = (+d1-s) * d4
         return
      end
c
      subroutine gauss31(ng,pg,wg)
c----------------------------------------------------------------------
c.... Position and weight of Gauss points  ngs = 1,2
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension pg(2,4),wg(4)
      call pzero(pg,8)
      call pzero(wg,4)
      pos = sqrt(3.0d0)/3.0d0
      if(ng.eq.1) then
        pg(1,1) = 0.0d0
        pg(2,1) = 0.0d0
        wg(1)   = 4.0d0
      else if(ng.eq.2) then
        pg(1,1) = -pos
        pg(2,1) = -pos
        pg(1,2) = -pos
        pg(2,2) =  pos
        pg(1,3) =  pos
        pg(2,3) = -pos
        pg(1,4) =  pos
        pg(2,4) =  pos
        wg(1) =  1.0d0
        wg(2) =  1.0d0
        wg(3) =  1.0d0
        wg(4) =  1.0d0
      end if
      return
      end
c
      subroutine trans31(d,xl,yl,ul,vl1,vlm,vl2,t,ndm,ndf,nen)
c---------------------------------------------------------------------------
c
c.... compute the transformation array and surface coords.
c     (see RLT-shell)
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension   x0(3),xl(ndm,*),yl(3,4),t(3,3),
     *            ul(ndf,*),vl1(6,4),vlm(6,4),vl2(6,4),d(*)
c
c.... compute the center (0,0)
      do i = 1,3
        x0(i) = 0.25*(xl(i,1) + xl(i,2) + xl(i,3) + xl(i,4))
      end do
      call pzero(t,9)

c     basic via diagonals        

c.... compute the inplane direction cosines (bisect diagonals)
      do i = 1,3
        t(1,i) = xl(i,3) - xl(i,1)
        t(2,i) = xl(i,2) - xl(i,4)
      end do    
      
      dl1 = sqrt(t(1,1)**2 + t(1,2)**2 + t(1,3)**2)
      dl2 = sqrt(t(2,1)**2 + t(2,2)**2 + t(2,3)**2)

      do i = 1,3
        v1 = t(1,i)/dl1
        v2 = t(2,i)/dl2
        t(1,i) = v1 + v2
        t(2,i) = v1 - v2
      end do    

      dl1 = sqrt(t(1,1)**2 + t(1,2)**2 + t(1,3)**2)
      dl2 = sqrt(t(2,1)**2 + t(2,2)**2 + t(2,3)**2)

c.... vectors t_1 and t_2
      do i = 1,3
        t(1,i) = t(1,i)/dl1
        t(2,i) = t(2,i)/dl2
      end do       

c.... compute the normal to the surface
      t(3,1) = t(1,2)*t(2,3) - t(2,2)*t(1,3)
      t(3,2) = t(1,3)*t(2,1) - t(2,3)*t(1,1)
      t(3,3) = t(1,1)*t(2,2) - t(2,1)*t(1,2)
c
c#####################
      igeo = d(17)
c.... modify base vector t_1, input must not be a unit vector  
      if(igeo.eq.2) then 
       tt = dot(d(18),d(18),3)
       if(tt.gt.0.d0) then 
         tt=sqrt(tt)
         t(1,1) = d(18)/tt
         t(1,2) = d(19)/tt
         t(1,3) = d(20)/tt

c....     calculate t_2 = t_3 x t_1
         t(2,1) = t(3,2)*t(1,3)-t(3,3)*t(1,2)
         t(2,2) =-t(3,1)*t(1,3)+t(3,3)*t(1,1)
         t(2,3) = t(3,1)*t(1,2)-t(3,2)*t(1,1)
       
c....     calculate t_1 = t_2 x t_3 new
         t(1,1) = t(2,2)*t(3,3)-t(2,3)*t(3,2)
         t(1,2) =-t(2,1)*t(3,3)+t(2,3)*t(3,1)
         t(1,3) = t(2,1)*t(3,2)-t(2,2)*t(3,1)
       end if
      end if
c#####################
c
c.... compute the projected middle surface coordinates
      do 140 i = 1,4
        do 140 j = 1,3
          yl(j,i) = 0.0
          do 130 k = 1,3
            yl(j,i) = yl(j,i) + t(j,k)*(xl(k,i) - x0(k))
130       continue
140   continue

c.... set offset coordinates to zero if small compared to plan size
      htol =  0.0d0
      do 150 i = 1,4
        htol = max(htol,abs(yl(1,i)),abs(yl(2,i)))
150   continue
      htol = htol*1.e-7
      do 160 i = 1,4
        if(abs(yl(3,i)) .le. htol) yl(3,i) = 0.0d0
160   continue

c.... compute the transformation of displacements
      call pzero(vl1,24) 
      call pzero(vlm,24) 
      call pzero(vl2,24) 

      do 170 i = 1,4
       do 170 j = 1,3
        do 180 k = 1,3

c....    at t_n      (only for nop=5)
         vl1(  j,i) = vl1(  j,i)+t(j,k)*(ul(  k,i)-ul(  k,i+nen))
         vl1(3+j,i) = vl1(3+j,i)+t(j,k)*(ul(3+k,i)-ul(3+k,i+nen))

c....    at t_n+1/2  (only for nop=5)       
         vlm(  j,i) = vlm(  j,i)+t(j,k)*(ul(  k,i)-ul(  k,i+nen)*0.5d0)
         vlm(3+j,i) = vlm(3+j,i)+t(j,k)*(ul(3+k,i)-ul(3+k,i+nen)*0.5d0)
          
c....    at t_n+1         
         vl2(  j,i) = vl2(  j,i)+t(j,k)*ul(  k,i)
         vl2(3+j,i) = vl2(3+j,i)+t(j,k)*ul(3+k,i)
180     continue
170   continue
      return
      end
c
      subroutine trans31_original(xl,yl,ul,vl,t,ndm,ndf)
c---------------------------------------------------------------------------
c.... compute the transformation array and surface coords.
c     (see RLT-shell)
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension   x0(3),xl(ndm,*),yl(3,4),t(3,3),ul(ndf,*),vl(6,4)
c
c.... compute the center (0,0)
      do i = 1,3
        x0(i) = 0.25*(xl(i,1) + xl(i,2) + xl(i,3) + xl(i,4))
      end do
      call pzero(t,9)

c     basic via diagonals        

c.... compute the inplane direction cosines (bisect diagonals)
      do i = 1,3
        t(1,i) = xl(i,3) - xl(i,1)
        t(2,i) = xl(i,2) - xl(i,4)
      end do    
      dl1 = sqrt(t(1,1)**2 + t(1,2)**2 + t(1,3)**2)
      dl2 = sqrt(t(2,1)**2 + t(2,2)**2 + t(2,3)**2)
      do i = 1,3
        v1 = t(1,i)/dl1
        v2 = t(2,i)/dl2
        t(1,i) = v1 + v2
        t(2,i) = v1 - v2
      end do    
      dl1 = sqrt(t(1,1)**2 + t(1,2)**2 + t(1,3)**2)
      dl2 = sqrt(t(2,1)**2 + t(2,2)**2 + t(2,3)**2)

c.... vectors t_1 and t_2
      do i = 1,3
        t(1,i) = t(1,i)/dl1
        t(2,i) = t(2,i)/dl2
      end do       

c.... compute the normal to the surface
      t(3,1) = t(1,2)*t(2,3) - t(2,2)*t(1,3)
      t(3,2) = t(1,3)*t(2,1) - t(2,3)*t(1,1)
      t(3,3) = t(1,1)*t(2,2) - t(2,1)*t(1,2)
c

c.... compute the projected middle surface coordinates
      do 140 i = 1,4
        do 140 j = 1,3
          yl(j,i) = 0.0
          do 130 k = 1,3
            yl(j,i) = yl(j,i) + t(j,k)*(xl(k,i) - x0(k))
130       continue
140   continue

c.... set offset coordinates to zero if small compared to plan size
      htol =  0.0d0
      do 150 i = 1,4
        htol = max(htol,abs(yl(1,i)),abs(yl(2,i)))
150   continue
      htol = htol*1.e-7
      do 160 i = 1,4
        if(abs(yl(3,i)) .le. htol) yl(3,i) = 0.0d0
160   continue

c.... compute the transformation of displacements
      do 170 i = 1,4
        do 170 j = 1,3
          vl(j,i)   = 0.0
          vl(j+3,i) = 0.0
          do 180 k = 1,3
            vl(j,i)    =  vl(j,i)   + t(j,k)* ul(k,i)
            vl(j+3,i)  =  vl(j+3,i) + t(j,k)* ul(k+3,i)
180       continue
170   continue
      return
      end
c
      subroutine transs31(d,xl,tr,ts,ndm)
c---------------------------------------------------------------------------
c.... compute the basis ts due to igeo for stress output
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),xl(ndm,*),x0(3),tr(3,3),ts(3,3),tt(3,3),add(3)
c
      igeo = d(17) 

      add(1)=d(18)
      add(2)=d(19)
      add(3)=d(20)

      call pzero(ts,9)

      if(igeo.le.2) then ! arbitrary
        do i=1,3
          do k=1,3 
            ts(i,k)=tr(i,k)
          end do
        end do
        if(igeo.eq.2) then ! t_3(1) + given vector t_1
          ts(1,1) = add(1)
          ts(1,2) = add(2)
          ts(1,3) = add(3)

c....     calculate t_2 = t_3 x t_1
          ts(2,1) = ts(3,2)*ts(1,3)-ts(3,3)*ts(1,2)
          ts(2,2) =-ts(3,1)*ts(1,3)+ts(3,3)*ts(1,1)
          ts(2,3) = ts(3,1)*ts(1,2)-ts(3,2)*ts(1,1)

c....     calculate t_1 = t_2 x t_3  new (if |t_1| .ne.1)
          ts(1,1) = ts(2,2)*ts(3,3)-ts(2,3)*ts(3,2)
          ts(1,2) =-ts(2,1)*ts(3,3)+ts(2,3)*ts(3,1)
          ts(1,3) = ts(2,1)*ts(3,2)-ts(2,2)*ts(3,1)
        end if

      else if(igeo.gt.2) then ! analytical solutions         

c....   compute the center (0,0)
        do i = 1,3
          x0(i) = 0.25*(xl(i,1) + xl(i,2) + xl(i,3) + xl(i,4))
        end do
        call pdirec3(tt,add,x0,igeo) ! ts=tt^T
c....   transpose tt
        do i=1,3
          do k=1,3 
            ts(i,k)=tt(k,i)
          end do
        end do
      end if
      
      return
      end
c
      subroutine transsig31(tr,ts,sig)
c---------------------------------------------------------------------------
c.... transform the stresses to basis from igeo
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension tr(3,3),ts(3,3),sig(8),sigs(8),t(3,3)

c.... Transformation matrix  T_ij=ts(i)*tr(j)
      call pzero(t,9) 
      do i = 1,3
        do j=1,3
          do k=1,3
            t(i,j) = t(i,j) + ts(i,k)*tr(j,k) 
          end do
        end do
      end do  

c.... Transform stresses

c.... normal forces 
      sigs(1) =      t(1,1)*t(1,1)*sig(1) + t(1,2)*t(1,2)*sig(2)
     +        + 2.d0*t(1,1)*t(1,2)*sig(3)
     +        + 2.d0*t(1,1)*t(1,3)*sig(7)
     +        + 2.d0*t(1,2)*t(1,3)*sig(8)
      sigs(2) =      t(2,1)*t(2,1)*sig(1) + t(2,2)*t(2,2)*sig(2)
     +        + 2.d0*t(2,1)*t(2,2)*sig(3)
     +        + 2.d0*t(2,1)*t(2,3)*sig(7)
     +        + 2.d0*t(2,2)*t(2,3)*sig(8)
      sigs(3) =      t(1,1)*t(2,1)*sig(1) + t(1,2)*t(2,2)*sig(2)
     +        +      t(1,1)*t(2,2)*sig(3) + t(1,2)*t(2,1)*sig(3)
     +        +     (t(1,1)*t(2,3) + t(1,3)*t(2,1))*sig(7)
     +        +     (t(1,2)*t(2,3) + t(1,3)*t(2,2))*sig(8)


c.... bending moments 
      sigs(4) =      t(1,1)*t(1,1)*sig(4) + t(1,2)*t(1,2)*sig(5)
     +        + 2.d0*t(1,1)*t(1,2)*sig(6)
      sigs(5) =      t(2,1)*t(2,1)*sig(4) + t(2,2)*t(2,2)*sig(5)
     +        + 2.d0*t(2,1)*t(2,2)*sig(6)
      sigs(6) =      t(1,1)*t(2,1)*sig(4) + t(1,2)*t(2,2)*sig(5)
     +        +      t(1,1)*t(2,2)*sig(6) + t(1,2)*t(2,1)*sig(6)

c.... shear forces
      sigs(7) =       t(1,1)*t(3,1)*sig(1) + t(1,2)*t(3,2)*sig(2)
     +        +      (t(1,1)*t(3,2) + t(1,2)*t(3,1))*sig(3)
     +        +      (t(1,1)*t(3,3) + t(1,3)*t(3,1))*sig(7)
     +        +      (t(1,2)*t(3,3) + t(1,3)*t(3,2))*sig(8)
      sigs(8) =       t(2,1)*t(3,1)*sig(1) + t(2,2)*t(3,2)*sig(2)
     +        +      (t(2,1)*t(3,2) + t(2,2)*t(3,1))*sig(3)
     +        +      (t(2,1)*t(3,3) + t(2,3)*t(3,1))*sig(7)
     +        +      (t(2,2)*t(3,3) + t(2,3)*t(3,2))*sig(8)

c.... store
      do i = 1,8
        sig(i) = sigs(i) 
      end do

      return
      end
c      
      subroutine transl31(s,p,t,nst,ndf,nop)
c---------------------------------------------------------------------------
c.... transform the loads and stiffness to global coords.
c     including symmetry terms
c     nop=5: no symmetry! 
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension s(nst,nst),p(nst),t(3,3),a(3,3),b(6)
c
      i0 = 0
      
      do 170 ir = 1,4
        do ii = 1,3
          b(ii  ) = dot(t(1,ii),p(i0+1),3)
          b(ii+3) = dot(t(1,ii),p(i0+4),3)
        end do
      
        do ii = 1,6
          p(i0+ii) = b(ii)
        end do    
      
        j0  = i0
        jc1 = ir
        if(nop.eq.5) then ! for unsymmetry
          j0  = 0
          jc1 = 1  
        end if
      
        do 160 jc = jc1,4
          i1 = i0
          do 150 i = 1,2
            j1 = j0
            do 140 j = 1,2
              do 120 ii = 1,3
                do 120 jj = 1,3
                  a(jj,ii) = dot(t(1,ii),s(i1+1,jj+j1),3)
120           continue
      
              do 130 ii = 1,3
                do 130 jj = 1,3
                  s(ii+i1,jj+j1) = dot(a(1,ii),t(1,jj),3)
130           continue

140         j1 = j1 + 3

150       i1 = i1 + 3

c.... compute the symmetric block
        if(nop.ne.5) then
          if(ir.ne.jc) then
            do 155 i = 1,6
             do 155 j = 1,6
              s(j0+j,i0+i) = s(i0+i,j0+j)
155         continue
          end if
        end if  

160     j0 = j0 + ndf

170   i0 = i0 + ndf
      return
      end
c
      subroutine plloco31(tr,xl)
c--------------------------------------------------------+
c.... Plot local basis
c--------------------------------------------------------+
      implicit double precision (a-h,o-z)
      dimension tr(3,3),x0(3),xl(3,*),d1(3),d2(3)

c...  plot at midpoint of element (global coordinates!)
      x0(1) = 0.25*(xl(1,1) + xl(1,2) + xl(1,3) + xl(1,4))
      x0(2) = 0.25*(xl(2,1) + xl(2,2) + xl(2,3) + xl(2,4))
      x0(3) = 0.25*(xl(3,1) + xl(3,2) + xl(3,3) + xl(3,4))

c.... length for plot

c.... compute the diagonals
      do i = 1,3
        d1(i) = xl(i,3) - xl(i,1)
        d2(i) = xl(i,2) - xl(i,4)
      end do      
      dl1 = sqrt(d1(1)**2 + d1(2)**2 + d1(3)**2)
      dl2 = sqrt(d2(1)**2 + d2(2)**2 + d2(3)**2)
      xm = 0.25d0 * min(dl1,dl2)

c...  plot axis 
      call pltaxs31(tr,x0,xm)
      return
      end
c
      subroutine pltaxs31(tr,x0,xm)
c----------------------------------------------------------------------
c.... draw vectors for axes, including rot-macro 
c.... t_1=black
c.... t_2=red
c.... t_3=blue
c----------------------------------------------------------------------
      USE pdata2 
      USE pltran 
      implicit double precision (a-h,o-z)
      logical zoom
      dimension xx(3,5),x0(3),x0z(3),tr(3,3)

c.... test if vectors are in plot-region
      x0z(1) = x0(1)
      x0z(2) = x0(2)
      x0z(3) = x0(3)
      call plxtrn(x0z,tra,vr,3,1)
      if(zoom(x0z,3,1)) then

c....   perspective projecton of axes
        do 120 m = 1,3  
          call pppcol(m)
          call pzero(xx,15)
          do 100 n = 1,3  
            fac1    = tr(m,n)*xm
            xx(n,1) = x0(n)
            xx(n,2) = xx(n,1) + fac1
100         xx(n,5) = xx(n,2)
          fac1 = tr(m,1)*xm
          fac2 = tr(m,2)*xm
          fac3 = tr(m,3)*xm
          xx(1,3) = xx(1,2) -.3*fac1 - .1*(fac2+fac3)
          xx(2,3) = xx(2,2) -.3*fac2 + .1*(fac1+fac3)
          xx(3,3) = xx(3,2) -.3*fac3 + .1*(fac1+fac2)
          xx(1,4) = xx(1,2) -.3*fac1 + .1*(fac2+fac3)
          xx(2,4) = xx(2,2) -.3*fac2 - .1*(fac1+fac3)
          xx(3,4) = xx(3,2) -.3*fac3 - .1*(fac1+fac2)

c.......  transform vector if necessary
          call plxtrn(xx,tra,vr,3,5)

c....     plot the vector
          call plotl(xx(1,1),xx(2,1),xx(3,1),3)
          call plotl(xx(1,2),xx(2,2),xx(3,2),2)
          call plotl(xx(1,2),xx(2,2),xx(3,2),ipgl)
          call plotl(xx(1,3),xx(2,3),xx(3,3),2)
          call plotl(xx(1,4),xx(2,4),xx(3,4),2)
          call plotl(xx(1,2),xx(2,2),xx(3,2),2)
          if(ipgl.eq.1) call clpan
          call plotl(xx(1,2),xx(2,2),xx(3,2),3)
          call plabl(m)
120     continue
      end if
      return
      end
c
      subroutine bdval31(xl,ndm)
c----------------------------------------------------------------------
c.....Bathe/Dvorkin-values
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /bdval1/ xxsim(4),yxsim(4),xetal(4),yetal(4)
      dimension xl(ndm,4)

c.... at M: B,D for all nodes BBDD
      xxsim(1)= 0.5d0*(xl(1,2)-xl(1,1))
      xxsim(2)= xxsim(1)
      xxsim(3)= 0.5d0*(xl(1,3)-xl(1,4))
      xxsim(4)= xxsim(3)
c
      yxsim(1)= 0.5d0*(xl(2,2)-xl(2,1))
      yxsim(2)= yxsim(1)
      yxsim(3)= 0.5d0*(xl(2,3)-xl(2,4))
      yxsim(4)= yxsim(3)
c

c.... at L: A,C for all nodes ACCA
      xetal(1)= 0.5d0*(xl(1,4)-xl(1,1))
      xetal(2)= 0.5d0*(xl(1,3)-xl(1,2))
      xetal(3)= xetal(2)
      xetal(4)= xetal(1)
c
      yetal(1)= 0.5d0*(xl(2,4)-xl(2,1))
      yetal(2)= 0.5d0*(xl(2,3)-xl(2,2))
      yetal(3)= yetal(2)
      yetal(4)= yetal(1)
c
      return
      end
c
      subroutine proj31(s,zi,zj,nst)
c----------------------------------------------------------------------
c.... modify stiffness for offset projections
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension   s(nst,*)

c.... postmultiply by transformation
      do 100 i = 1,5
        s(i,4) = +zj*s(i,2) + s(i,4)
        s(i,5) = -zj*s(i,1) + s(i,5)
100   continue

c.... premultiply using modified terms from postmultiplication
      do 200 i = 1,5
        s(4,i) = +zi*s(2,i) + s(4,i)
        s(5,i) = -zi*s(1,i) + s(5,i)
200   continue
      return
      end
c
      subroutine shearfac31(d,xl,cappa)
c-----------------------------------------------------------------------
c.... shear correction factor due to size
c.... shear correction factor
c     d(13) <0: cappa = input value*SCF(Tessler) 
c     d(13) >0: cappa = input value 
c     d(13) =0: cappa = 5/6 (default)   
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  d(*),xl(3,*),yl(3)

      cappa = d(13)
     
      if(cappa.lt.0.d0)then

c....   AE... element area
c       ae = 4.d0*detj0 
       
c....   AE... or square of longest element side
        hs = d(47)
        xnu = 0.d0   
        call pzero(yl,3) 
        ae = 0.d0
        do i = 1,4
          k = i+1
          if(i.eq.4) k=1
          do j = 1,3
            yl(j) = xl(j,k)-xl(j,i)
          end do 
          sl2 = dot(yl,yl,3)       
          ae = max(ae,sl2)
        end do  
       
        cappa  = dabs(d(13))/(1.d0+ae/(2.d0*hs*hs*(1.d0+xnu))) ! Tessler

      end if 
c      
      return
      end  
c
      subroutine qload31(qld,q,d,numel,n,mqloa,propq,prop,isw)
c----------------------------------------------------------
c.... calculate loads from macro qloa
c
c     qld(01) = qz 
c     qld(02) = qn
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),d(*),qld(2)
      call pzero(qld,2)
      if(isw.eq.4.or.isw.eq.8) then   ! STRE
        qld(1) = d( 9)*prop
        qld(2) = d(10)*prop
        if(mqloa.ne.1) then
          qld(1) = qld(1) + q(n,1)*propq 
          qld(2) = qld(2) + q(n,2)*propq 
        end if 
      else if(isw.eq.22) then         ! only from QLOA
        if(mqloa.ne.1) then
          qld(1) = q(n,1)*propq 
          qld(2) = q(n,2)*propq 
        end if 
      else                            ! only from MATE
        qld(1) = d( 9)*prop
        qld(2) = d(10)*prop
      end if

      return
      end

      subroutine damage31(d,isw)
c----------------------------------------------------------------------
c     isw =1:  read failure     values    d(21)-d(28)
c              set              values    d(29)-d(36)
c     isw =2:  read degradation values    d(37)-d(44)
c              set              values    d(45)-d(46)

c----------------------------------------------------------------------
      USE iofile 
      implicit double precision(a-h,o-z)
      dimension d(48)

      if(isw.eq.1) then

c....   read failure data       
        call pzero(d(21),16)
        call dinput(d(21),8)
        iftyp=d(21)

        if(iftyp.eq.1) then
c....     data for Tsai-Wu-criterion        
c  21     iftyp 1
c  22     X_11  Bruchspannung bei Zug in Faserrichtung  ....values without sign 
c  23     Y_11  Bruchspannung bei Druck in Faserrichtung
c  24     X_22  Bruchspannung bei Zug quer zur Faserrichtung
c  25     Y_22  Bruchspannung bei Druck quer zur Faserrichtung
c  26     X_12  Bruchschubspannung 12
c  27     X_13  Bruchschubspannung 13
c  28     X_23  Bruchschubspannung 23
c......   calculate parameters for tsaiwu criterion
          d(29)=1.d0/d(22)-1.d0/d(23) ! q1 = 1/X_11 - 1/Y_11
          d(30)=1.d0/(d(22)*d(23))    ! a1 = 1/(X_11 * Y_11)
          d(31)=1.d0/d(24)-1.d0/d(25) ! q2 = 1/X_22 - 1/Y_22
          d(32)=1.d0/(d(24)*d(25))    ! a2 = 1/(X_22 * Y_22)
          d(33)= sqrt(d(28)*d(30))    ! a12= dsqrt(f11 * f22)
          d(34)=1.d0/(d(26)*d(26))    ! a4 = 1/(X_12 * X_12)
          d(35)=1.d0/(d(27)*d(27))    ! a5 = 1/(X_13 * X_13)
          d(36)=1.d0/(d(28)*d(28))    ! a6 = 1/(X_23 * X_23)
          write(iow,2001) (d(k),k=22,36)


        else
          write(iow,2010)
        end if

      else if(isw.eq.2) then

c....   read degradation data       
        call pzero(d(37),10)
        call dinput(d(37),8)
        idtyp=d(37)

        if(idtyp.eq.1) then
c....     data for degradation model 1        
c  37     idtyp 1
          write(iow,3001) (d(k),k=37,37)
        else
          write(iow,3010)
        end if

      end if
c.... formats
2001  format(
     + 5x,'TSAI-WU failure criterion (values without sign)  ',/,
     + 5x,'Versagenskriterium                               ',f12.1,/,
     + 5x,'X_11  Bruchspannung Zug in Faserrichtung         ',e12.5,/,
     + 5x,'Y_11  Bruchspannung Druck in Faserrichtung       ',e12.5,/,
     + 5x,'X_22  Bruchspannung Zug quer zur Faserrichtung   ',e12.5,/,
     + 5x,'Y_22  Bruchspannung Druck quer zur Faserrichtung ',e12.5,/,
     + 5x,'X_12  Bruchschubspannung 12                      ',e12.5,/,
     + 5x,'X_13  Bruchschubspannung 13                      ',e12.5,/,
     + 5x,'X_23  Bruchschubspannung 23                      ',e12.5,/,
     + 5x,'q1 = 1/X_11 - 1/Y_11                             ',e12.5,/,
     + 5x,'a1 = 1/(X_11 * Y_11)                             ',e12.5,/,
     + 5x,'q2 = 1/X_22 - 1/Y_22                             ',e12.5,/,
     + 5x,'a2 = 1/(X_22 * Y_22)                             ',e12.5,/,
     + 5x,'a12= dsqrt(f11 * f22)                            ',e12.5,/,
     + 5x,'a4 = 1/(X_12 * X_12)                             ',e12.5,/,
     + 5x,'a5 = 1/(X_13 * X_13)                             ',e12.5,/,
     + 5x,'a6 = 1/(X_23 * X_23)                             ',e12.5)
2010  format(5x,'no failure criterion used')    
3001  format(
     + 5x,'Degradation model                                ',/,
     + 5x,'Model                                            ',f12.1)

3010  format(5x,'no degradation model used')    
      return
      end  
c
      subroutine elast31(h11l,h22l,h12l,h33l,q11l,q22l,e1,e2,v12,
     +  g12,g23,cappa,em,eb,es,zs,d,h1,h2,iftyp,idtyp,c2,s2,c,s,idam)
c----------------------------------------------------------------------
c     Calculate local elasticity matrix including damage

c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension d(48),epsg(5),epsl(5),em(3),eb(3),es(3)

c.... undamaged
      e1q  = e1
      e2q  = e2
      f    = 1.d0/(1.d0-v12*v12*e2q/e1q)
      h11l = f*e1q
      h12l = f*v12*e2q
      h22l = f*e2q
      h33l = g12
      q11l = g12*cappa
      q22l = g23*cappa                

      idv1=h1 ! failure criterion: 0=undam. 1=dam
      idv2=0

      if(idv1.eq.1) then ! damaged
        if(idtyp.eq.1) then   ! Degradation model 1
          e1q = e2
          e2q = 0.1d0*e2
          f   = 1.d0/(1.d0-v12*v12*e2q/e1q)
c....     modify local elasticity matrix of layer
          h11l = f*e1q
          h12l = f*v12*e2q
          h22l = f*e2q
        end if
      else ! check if now damaged 
c....   global strains at position zs
        epsg(1) = em(1) + zs*eb(1)
        epsg(2) = em(2) + zs*eb(2)
        epsg(3) = em(3) + zs*eb(3)
        epsg(4) = es(1) 
        epsg(5) = es(2) 

c....   local  strains at position zs:  eps_l = t_e * eps_g
        epsl(1) =      c2*epsg(1) +     s2*epsg(2) +   s *c *epsg(3)
        epsl(2) =      s2*epsg(1) +     c2*epsg(2) -   s *c *epsg(3)
        epsl(3) = -2*s*c *epsg(1) + 2*s*c *epsg(2) + (c2-s2)*epsg(3)
        epsl(4) =      c *epsg(4) +     s *epsg(5) 
        epsl(5) =     -s *epsg(4) +     c *epsg(5) 

c....   local stresses
        s11 =h11l*epsl(1) + h12l*epsl(2) 
        s22 =h12l*epsl(1) + h22l*epsl(2) 
        s12 =h33l*epsl(3)                 
        s13 =q11l*epsl(4)                 
        s23 =q22l*epsl(5)                 
        if(iftyp.eq.1) then
c....     Tsai-Wu-criterion        
c         q1*s11+a1*s11*s11+a12*s11*s22+q2*s22+a2*s22*s22
c         +a4*s12*s12+a5*s13*s13+a6*s23*s23<1
          q1 =d(29)
          a1 =d(30)
          q2 =d(31)
          a2 =d(32)
          a12=d(33)
          a4 =d(34)
          a5 =d(35)
          a6 =d(36)

          fdam = q1*s11+a1*s11*s11+a12*s11*s22+q2*s22+a2*s22*s22
     1         + a4*s12*s12 + a5*s13*s13+ a6*s23*s23
          if(fdam.gt.1.d0) then
            if(idtyp.eq.1) then   ! Degradation model 1
c....         modify material data  
              e1q = e2
              e2q = 0.1d0*e2
              f   = 1.d0/(1.d0-v12*v12*e2q/e1q)
            end if
c....       modify local elasticity matrix of layer
            h11l = f*e1q
            h12l = f*v12*e2q
            h22l = f*e2q
c....       set failure undex to 1
            idv2=1
          end if
        end if 
      end if
      ih2= max(idv1,idv2)  ! no undamaging possible  
      idam=max(idam,ih2)
c       
      return
      end
