      subroutine fe2micro(a,g1,g2,gl,ixl,ht,tau,melnk,lelnk,x,u,dr,s,p,
     +                    xl,id,ix,numnp,numel,nen,neq,flnk)
c ----------------------------------------------------------------------
c.... Purpose: calculate C=Sum A^T*K*A-G^T*K_aa^-1*G and Tau=-Sum A^T*P
c              for microscale problem
c
c     natyp=1 (3d)        eps=[E_11,E_22,E_33,2E_12,2E_13,2E_23]           
c     natyp=2 (shell)     eps=[E_11,E_22,2E_12,K_11,K_22,2K_12,2E_13,2E_23]
c     natyp=3 (shell ext) eps=[E_11,E_22,2E_12,K_11,K_22,2K_12,2E_13,2E_23,E_330,E_331]
c     natyp=4 (beam)      eps=[E_11,2E_12,2E_13,K_11,K_22,K_33]
c     natyp=5 (beam ext)  eps=[E_11,2E_12,2E_13,K_11,K_22,K_33,E_220,E_221,E_222,E_330,E_331,E_332,E_230,E_231,E_233]
c
c.... Usage
c            MACRO>EPSQ     
c            Macro>Loop,,N   
c            Macro>TANG,,1     
c            Macro>NEXT      
c            Macro>SIGQ
c
c....  Inputs:
c      a(nst,nss)          - local  matrix A  
c      g1(neq,nss)         - global matrix G1
c      g2(neq,nss)         - global matrix G2
c      gL(nst,nss)         - local  matrix G  
c      ixl(ndf,nen)        - DOF indicators 
c      melnk(ndf,numnp)    - link conditions
c      x(ndm,numnp)        - global coordinates
c      u(*)                - global displacements
c      s(nst,nst)          - element stiffness matrix   
c      p(nst)              - element load vector  
c      xl(ndm,nen)         - element nodal coordinates  
c      id(ndf,numnp)       - equation numbers for each active dof
c      ix(nen1,numel)      - element nodal connections of mesh
c           numnp          - number of nodes in mesh   
c           numel          - number of elements in mesh
c           ndm            - spatial dimension of mesh
c           ndf            - number dof/node
c           nen            - max. number of nodes/element
c           nen1           - dimension for ix array: nen+4
c           nst            - dimension for element array: ndf*nen
c           neq            - number of active dofs in mesh
c           nss            - length of strain vector  
c           ibin           - 0=binary transfer,1=ascii, 2=shared memory
c           flnk           - Flag T=link, F=nolink 
c          
c....  Working:    
c        ixl(ndf,nen)      - DOF indicators: 
c                             0 = free dof
c                             1 = boundary dof
c                            -1 = no node
c        lelnk(ndf,nen)    - link indicators: 
c                             idf.ne. 0 = link to node idf 
c                             0 = no link
c        isflg             - f = element without b.c./link nodes 
c                          - t = element with b.c./link nodes
c
c      Outputs:
c        ht(nss,nss)       - C=Sum A^T*K*A - G1^T*X2  with K_uu*X2=G2    
c        tau(nss)          - Tau=-Sum A^T*P
c
c      To be defined:
c           ibin           - 0=binary transfer,1=ascii, 2=shared memory
c
c ----------------------------------------------------------------------
      USE comfil
      USE epsdh       !for iswm,nss,natyp,skfy,skfz,dxyz
      USE epsd1       ! for it1,it2,resife2 in FE2
      USE hdatam
      USE iofile
      USE pdata2      ! for idev
      USE sdata
      USE shmname
      implicit double precision(a-h,o-z)

c.... global
      dimension a(nst,nss),g1(neq,nss),g2(neq,nss),gl(nst,nss),
     +  ixl(ndf,nen),melnk(ndf,numnp),lelnk(ndf,nen),s(nst,nst),p(nst),
     +  id(ndf,numnp),ix(nen1,numel),u(*),dr(*),xl(ndm,nen),
     +  x(ndm,numnp),ht(nss,nss),tau(nss)

c.... local
      dimension dummy(1)

c.... test
      dimension htp(nss,nss),taup(nss)

      character*229 fcisb
      logical sflg,flnk
      data bl/-999.d0/
      
#ifdef __INTEL_      
      INTEGER*8 :: csmap !handles to store file mappings
#endif

c.... name of transfer file  
      iex=ipos(  finp,229)  
      fcisb = finp(1:iex) 
      call dochar2(fcisb,ip)     ! look for i
      call dochar1(fcisb,'b',ip) ! set i -> b

c...  transfer version
      if(idev.eq.3) ibin=0 ! INTEL  
      if(idev.eq.4) ibin=0 ! SALFORD
 
c.... L = l_xm, A=l_x*l_y, V=A*l_z 
      rlen  = dxyz(1)
      area  = dxyz(1)*dxyz(2)
      volm  = dxyz(1)*dxyz(2)*dxyz(3)
      if(natyp.eq.1)volmr = 1.d0/volm
      if(natyp.eq.2)volmr = 1.d0/area
      if(natyp.eq.3)volmr = 1.d0/area
      if(natyp.eq.4)volmr = 1.d0/rlen
      if(natyp.eq.5)volmr = 1.d0/rlen

c.... loop over elements

      tau=0.d0
      ht =0.d0
      g1 =0.d0
      g2 =0.d0

      hflgu  = .false.
      h3flgu = .false.

      do nn = 1,numel
         
        gl = 0.d0
 
c       Check for b.c.
c       Set 'ixl' array to mark dofs with fixed boundaries
        nel = 0  ! not known, set locally in pform
        ixl = 0 
        sflg = .false.
         
        do i = 1,nen
          il = ix(i,nn)
          if(il.gt.0) then
            nel    = i
            do j = 1,ndf
              if(id(j,il).gt.0) then 
                ixl(j,i) = 0
              else if(id(j,il).lt.0) then 
                ixl(j,i) = 1
                sflg = .true. ! element has b.c. nodes 
              end if
            end do ! j
          else !  No node
            do j = 1,ndf
              ixl(j,i) = -1
            end do ! j
          end if ! il
        end do ! i

c       Check for periodic b.c.
        if(flnk) call ochkp2(ix(1,nn),melnk,ndf,nel,sflg)

        if(sflg) then

c         Get element tangent and residual for elmt nn, No assembly
c
          call formfe(u,dummy,dummy,dummy,
     +    .false.,.false.,.false.,.false.,3,nn,nn,1)

c....     xl,s,p are set in formfe/pform         

          if(flnk) then 
c....       Periodic boundary case
            call uprojpp2(ixl,lelnk,melnk,id,ix(1,nn),x,xl,s,p,
     +          g1,gl,ht,tau,a,skfy,skfz,ndm,ndf,nel,nst,neq,nss,
     +          natyp,iswm)
          else
c....       Displacement boundary case
            call uprojpd2(ixl,id,ix(1,nn),xl,s,p,g1,gl,ht,tau,a,skfy,
     +                    skfz,ndm,ndf,nel,nst,neq,nss,natyp,iswm)
          end if
        end if ! sflg 

      end do ! nn
      
      call mprinte(ht,nss,nss,nss,1.d-8,'Dmat_vor')

c      if(iswm.eq.3) then

c....   copy G1 into G2
        call pmove(g1,g2,neq*nss)

c....   material moduli by static condensation
        call formht(ht,g1,g2,nss,neq)

c....   stresses by static condensation
c       call mprint(tau,nss,1,nss,'Sig_vor')  
c       call formtau(tau,g1,u,dr,id,nss,neq,ndf,numnp)
c       call mprint(tau,nss,1,nss,'Sig_nach')  

c      else 

c      end if ! iswm
    
c.... average over RVE  
      call msmul(ht,nss,nss,volmr)
      call msmul(tau,nss,1,volmr)

c.... C after condensation
      call mprinte(ht,nss,nss,nss,1.d-8,'Dmat_nach')  

c.... Tau after condensation
      call mprinte(tau,nss,1,nss,1.d-8,'Sig_nach')  

c.... Test RVE
      if(iepsd.eq.1) then
       if(natyp.eq.3) then
c....    Condensation nss->8
         htp=ht
         taup=tau
         call conden(htp,taup,8,nss-8,nss,.false.)
         call mprinte(htp,8,8,nss,1.d-8,'Dmat')  
         call mprinte(taup,8,1,nss,1.d-8,'Sig')  
       end if
       
       if(natyp.eq.5) then
c....    Condensation nss->6
         htp=ht
         taup=tau
         call conden(htp,taup,6,nss-6,nss,.false.)
         call mprinte(htp,6,6,nss,1.d-8,'Dmat')  
         call mprinte(taup,6,1,nss,1.d-8,'Sig')  
       end if
      endif

c.... Send C and S to macro problem
      if (iepsd.eq.0) then
       if(ibin.eq.0) then
         open(34,file=fcisb,form='unformatted',status='unknown')
         rewind(34)
         write(34) tau
         write(34) ht
c...     code only for ibin=0! ibin=1,2 open! 
         write(34) it1,it2,resife2 ! it-behavior micro-problem

         close(34)

         write(iow,'(a,i3,i3,e16.5)') 'it1,it2,resife2',it1,it2,resife2

       else if(ibin.eq.2) then

#ifdef __INTEL_
        
         call feap_shmem_OpenCSMapping(csmap,shmbasename)         
         call feap_shmem_writeCS(csmap,tau,ht,nss)
         call feap_shmem_CloseHandle(csmap)

#endif
       
       else 
         open(34,file=fcisb,form='formatted',status='unknown')
         rewind(34)
         write(34,'(8e18.9)') (tau(i),i=1,nss) 
         do i = 1,nss
           write(34,'(8e18.9)') (ht(i,k),k=1,nss) 
         end do      
         close(34)
       
       end if
      end if 

      return
      end
c
      subroutine uprojpp2(ixl,lelnk,melnk,id,ix,x,xl,s,p,g1,gl,ht,tau,a,
     +                    skfy,skfz,ndm,ndf,nel,nst,neq,nss,natyp,iswm)
c-----------------------------------------------------------------------
c      Purpose: Set Coupling and diagonal arrays and compute averaged
c               TAU stress and tangent modulus arrays 
c               for classical and periodic b.c.s
c
c      Inputs:
c        ixl(ndf,nen)    - DOF indicators: -1 = no equation
c                                           0 = free dof
c                                           1 = boundary dof
c        melnk(ndf,numnp)- link conditions
c        id(ndf,numnp)   - Equation numbers at nodes
c        ix(nen1)        - Element connection list for this element
c        x(ndm,numnp)    - Nodal coordinates
c        xl(ndm,nel)     - Element nodal coordinates
c        s(nst,nst)      - Element tangent matrix
c        skfy            - shear correction factor y
c        skfz            - shear correction factor z
c        ndm             - Spatial dimension of mesh
c        ndf             - dofs 
c        nel             - Number of maximum node on element
c        nst             - Dimension of tangent matrix
c        neq             - Number of active equations in mesh
c        nss             - Number of modes to project = 6 or 8
c        isw             - Solution switch 
c
c      Working:
c        a (nst,nss)     - displ.-strain matrix U=A*Eps at b.c.
c        gl(nst,nss)     - GL = K*A Coupling matrix local
c
c      Outputs:
c        lelnk(ndf,nen)  - link indicators: i>0= link to node i 
c                                           0 = no link
c        ht(nss,nss)     - C=C+A^T*K*A
c        tau(nss)        - Tau=Tau-A^T*P
c        g1(neq,nss)     - G1 Coupling matrix  global
c
c-----------------------------------------------------------------------

      implicit double precision(a-h,o-z)
      dimension ixl(ndf,*),lelnk(ndf,*),melnk(ndf,*),id(ndf,*),
     +          ix(*),x(ndm,*),xl(ndm,nel),s(nst,nst),p(nst),
     +          xp(nel),yp(nel),zp(nel),
     +          g1(neq,nss),gl(nst,nss),ht(nss,nss),tau(nss),a(nst,nss)
      data d5/0.5d0/

c     set array lelnk and modify nodal coordinates for linked nodes
c     DX=Xslave=Xslave-Xmaster 
      call pzeroi(lelnk,nel*ndf)

      ndm1 = 3
      if(natyp.eq.2)ndm1 = 2
      
      do ir = 1,nel
        xp(ir) = xl(1,ir)
        yp(ir) = xl(2,ir)
        zp(ir) = xl(3,ir)
        if(ix(ir).gt.0) then ! no node
          do i = 1,ndm1
            ipl = melnk(i,ix(ir))     ! number of master node
            
            if(ipl.gt.0.and.ixl(i,ir).ne.1) then   ! if link
              lelnk(i,ir) = ipl ! set lelnk 
              xl(i,ir)    = xl(i,ir) - x(i,ipl)    ! Delta x
            end if
          end do ! i
        else
          do i = 1,ndm1
            lelnk(i,ir) = 0
          end do ! i
        end if
      end do ! ir

      a = 0.d0
      
      if(natyp.eq.1) then ! 3D
        do i = 1,nel
          ia = (i-1)*ndf 
          if(ixl(1,i).eq.1 .or. lelnk(1,i).gt.0) then  
            a(ia+1,1) = xl(1,i)
            a(ia+1,4) = xl(2,i)*d5
            a(ia+1,5) = xl(3,i)*d5
          end if
          if(ixl(2,i).eq.1 .or. lelnk(2,i).gt.0) then  
            a(ia+2,2) = xl(2,i)
            a(ia+2,4) = xl(1,i)*d5
            a(ia+2,6) = xl(3,i)*d5
          end if
          if(ixl(3,i).eq.1 .or. lelnk(3,i).gt.0) then  
            a(ia+3,3) = xl(3,i)
            a(ia+3,5) = xl(1,i)*d5
            a(ia+3,6) = xl(2,i)*d5
          end if
        end do ! i
      else if(natyp.eq.2) then ! Shell
        do i = 1,nel
          ia = (i-1)*ndf 
          if(ixl(1,i).eq.1 .or. lelnk(1,i).gt.0) then  
            a(ia+1,1) = xl(1,i)
            a(ia+1,3) = xl(2,i)*d5
            a(ia+1,4) = xl(1,i)*zp(i)
            a(ia+1,6) = xl(2,i)*zp(i)*d5
            a(ia+1,7) = xl(3,i)*skfy
          end if
          if(ixl(2,i).eq.1 .or. lelnk(2,i).gt.0) then  
            a(ia+2,2) = xl(2,i)
            a(ia+2,3) = xl(1,i)*d5
            a(ia+2,5) = xl(2,i)*zp(i)
            a(ia+2,6) = xl(1,i)*zp(i)*d5
            a(ia+2,8) = xl(3,i)*skfz
          end if
        end do ! i
      else if(natyp.eq.3) then ! Shell ext.
        do i = 1,nel
          ia = (i-1)*ndf 
          if(ixl(1,i).eq.1 .or. lelnk(1,i).gt.0) then
            a(ia+1,1) = xl(1,i)
            a(ia+1,3) = xl(2,i)*d5
            a(ia+1,4) = xl(1,i)*zp(i)
            a(ia+1,6) = xl(2,i)*zp(i)*d5
            a(ia+1,7) = xl(3,i)*d5*skfy
          end if
          if(ixl(2,i).eq.1 .or. lelnk(2,i).gt.0) then
            a(ia+2,2) = xl(2,i)
            a(ia+2,3) = xl(1,i)*d5
            a(ia+2,5) = xl(2,i)*zp(i)
            a(ia+2,6) = xl(1,i)*zp(i)*d5
            a(ia+2,8) = xl(3,i)*d5*skfz
          end if
          if(ixl(3,i).eq.1 .or. lelnk(3,i).gt.0) then
            xm = xp(i)-xl(1,i)
            ym = yp(i)-xl(2,i)
            a(ia+3,4) = -(xp(i)*xp(i)-xm*xm)*d5
            a(ia+3,5) = -(yp(i)*yp(i)-ym*ym)*d5
            a(ia+3,6) = -(xp(i)*yp(i)-xm*ym)*d5
            a(ia+3,7) =  xl(1,i)*d5*skfy
            a(ia+3,8) =  xl(2,i)*d5*skfz
            a(ia+3,9) =  xl(3,i)
            a(ia+3,10)=  xl(3,i)*zp(i)
          end if
        end do ! i
       
      else if(natyp.eq.4) then ! Beam
        do i = 1,nel
          ia = (i-1)*ndf 
          if(ixl(1,i).eq.1 .or. lelnk(1,i).gt.0) then  
            a(ia+1,1) = xl(1,i)
            a(ia+1,2) = xl(2,i)*d5 !*skfy
            a(ia+1,3) = xl(3,i)*d5 !*skfz
            a(ia+1,4) = 0.d0 ! (yp(i)*xl(3,i)-zp(i)*xl(2,i))*d5
            a(ia+1,5) = xl(1,i)*zp(i)
            a(ia+1,6) =-xl(1,i)*yp(i)
          end if
          if(ixl(2,i).eq.1 .or. lelnk(2,i).gt.0) then  
            a(ia+2,2) = xl(1,i)*d5*skfy
            a(ia+2,4) =-xl(1,i)*zp(i)
          end if
          if(ixl(3,i).eq.1 .or. lelnk(3,i).gt.0) then  
            a(ia+3,3) = xl(1,i)*d5*skfz
            a(ia+3,4) = xl(1,i)*yp(i)
          end if
        end do ! i
      else if(natyp.eq.5) then ! Beam ext.

      end if

c.... calculate Tau 
      tau = tau-matmul(transpose(a),p)  

      if(iswm.eq.3.or.iswm.eq.8) then 
      
c....   calculate C and GL

        ht = ht+matmul(transpose(a),matmul(s,a))
        gl = matmul(s,a)  

c....   store GL  into G 
        do i = 1,nel
          ii = ix(i) ! node number in Mesh
          do j = 1,ndf
            jj = id(j,ii) ! dof number in G1
            kk = (i-1)*ndf+j ! row in GL 
            if(jj.gt.0) then
              do k = 1,nss ! column in GL 
                g1(jj,k) = g1(jj,k) + gl(kk,k)
              end do ! k
            end if ! jj
          end do ! j
        end do ! i
      end if ! iswm

      return
      end
c
      subroutine uprojpd2(ixl,id,ix,xl,s,p,g1,gl,ht,tau,a,skfy,skfz,
     +                    ndm,ndf,nel,nst,neq,nss,natyp,iswm)
c-----------------------------------------------------------------------
c      Purpose: Set Coupling and diagonal arrays and compute averaged
c               TAU stress and tangent modulus arrays 
c               for classical b.c.s
c
c      Inputs:
c        ixl(ndf,nen)    - DOF indicators: -1 = no equation
c                                           0 = free dof
c                                           1 = boundary dof
c        id(ndf,numnp)   - Equation numbers at nodes
c        ix(nen1)        - Element connection list for this element
c        xl(ndm,nel)     - Element nodal coordinates
c        s(nst,nst)      - Element tangent matrix
c        skfy            - shear correction factor y
c        skfz            - shear correction factor z
c        ndm             - Spatial dimension of mesh
c        ndf             - dofs 
c        nel             - Number of maximum node on element
c        nst             - Dimension of tangent matrix
c        neq             - Number of active equations in mesh
c        nss             - Number of modes to project = 6 or 8
c        isw             - Solution switch 
c
c      Working:
c        a (nst,nss)     - displ.-strain matrix U=A*Eps at b.c.
c        gl(nst,nss)     - GL = K*A Coupling matrix local
c
c      Outputs:
c        ht(nss,nss)     - C=C+A^T*K*A
c        tau(nss)        - Tau=Tau-A^T*P
c        g1(neq,nss)     - G1 Coupling matrix  global
c
c-----------------------------------------------------------------------

      implicit double precision(a-h,o-z)
      dimension ixl(ndf,*),id(ndf,*),ix(*),xl(ndm,nel),
     +          s(nst,nst),p(nst),g1(neq,nss),gl(nst,nss),ht(nss,nss), 
     +          tau(nss),a(nst,nss)
      data d5/0.5d0/  

      a = 0.d0
      if(natyp.eq.1) then ! 3D
        do i = 1,nel
          ia = (i-1)*ndf 
          if(ixl(1,i).eq.1) then  
            a(ia+1,1) = xl(1,i)
            a(ia+1,4) = xl(2,i)*d5
            a(ia+1,5) = xl(3,i)*d5
          end if
          if(ixl(2,i).eq.1) then  
            a(ia+2,2) = xl(2,i)
            a(ia+2,4) = xl(1,i)*d5
            a(ia+2,6) = xl(3,i)*d5
          end if
          if(ixl(3,i).eq.1) then  
            a(ia+3,3) = xl(3,i)
            a(ia+3,5) = xl(1,i)*d5
            a(ia+3,6) = xl(2,i)*d5
          end if
        end do ! i
      else if(natyp.eq.2) then ! Shell
        stop ' stop in uprojpd2, part not updated'
        do i = 1,nel
          ia = (i-1)*ndf 
          if(ixl(1,i).eq.1) then  
            a(ia+1,1) = xl(1,i)
            a(ia+1,3) = xl(2,i)*d5
            a(ia+1,4) = xl(1,i)*xl(3,i)
            a(ia+1,6) = xl(2,i)*xl(3,i)*d5
            a(ia+1,7) = xl(3,i)*skfy
          end if
          if(ixl(2,i).eq.1) then  
            a(ia+2,2) = xl(2,i)
            a(ia+2,3) = xl(1,i)*d5
            a(ia+2,5) = xl(2,i)*xl(3,i)
            a(ia+2,6) = xl(1,i)*xl(3,i)*d5
            a(ia+2,8) = xl(3,i)*skfz
          end if
        end do ! i
      else if(natyp.eq.3) then ! Shell ext.
        stop ' stop in uprojpd2, part not updated'
        do i = 1,nel
          ia = (i-1)*ndf 
          if(ixl(1,i).eq.1) then  
            a(ia+1,1) = xl(1,i)
            a(ia+1,3) = xl(2,i)*d5
            a(ia+1,4) = xl(1,i)*xl(3,i)
            a(ia+1,6) = xl(2,i)*xl(3,i)*d5
            a(ia+1,7) = xl(3,i)*skfy
          end if
          if(ixl(2,i).eq.1) then  
            a(ia+2,2) = xl(2,i)
            a(ia+2,3) = xl(1,i)*d5
            a(ia+2,5) = xl(2,i)*xl(3,i)
            a(ia+2,6) = xl(1,i)*xl(3,i)*d5
            a(ia+2,8) = xl(3,i)*skfz
          end if
          if(ixl(3,i).eq.1) then  
            a(ia+3,4) = -xl(1,i)*xl(1,i)*d5
            a(ia+3,5) = -xl(2,i)*xl(2,i)*d5
            a(ia+3,6) = -xl(1,i)*xl(2,i)*d5
            a(ia+3,9) =  xl(3,i)
            a(ia+3,10)=  xl(3,i)*xl(3,i)
          end if
        end do ! i
      else if(natyp.eq.4) then ! Beam
        stop ' stop in uprojpd2, part not updated'
        do i = 1,nel
          ia = (i-1)*ndf 
          if(ixl(1,i).eq.1 ) then  
            a(ia+1,1) = xl(1,i)
            a(ia+1,2) = xl(2,i)*d5 !*skfy
            a(ia+1,3) = xl(3,i)*d5 !*skfz
            a(ia+1,4) = 0.d0
            a(ia+1,5) = xl(1,i)*xl(3,i)
            a(ia+1,6) =-xl(1,i)*xl(2,i)
          end if
          if(ixl(2,i).eq.1) then  
            a(ia+2,2) = xl(1,i)*d5 !*skfy
            a(ia+2,4) =-xl(1,i)*xl(3,i)
          end if
          if(ixl(3,i).eq.1) then  
            a(ia+3,3) = xl(1,i)*d5 !*skfz
            a(ia+3,4) = xl(1,i)*xl(2,i)
          end if
        end do ! i
      else if(natyp.eq.5) then ! Beam ext.

      end if ! natyp

c.... calculate Tau 
      tau = tau-matmul(transpose(a),p)  

      if(iswm.eq.3.or.iswm.eq.8) then
c....   calculate C 
        ht = ht+matmul(transpose(a),matmul(s,a))
      
c....   calculate GL
        gl = matmul(s,a)  

c....   store GL  into G 
        do i = 1,nel
          ii = ix(i) ! node number in Mesh
          do j = 1,ndf
            jj = id(j,ii) ! dof number in G1
            kk = (i-1)*ndf+j ! row in GL 
            if(jj.gt.0) then
              do k = 1,nss ! column in GL 
                g1(jj,k) = g1(jj,k) + gl(kk,k)
              end do ! k
            end if ! jj
          end do ! j
        end do ! i
      end if ! iswm
      return
      end
c
      subroutine ochkp2(ix,melnk,ndf,nel,sflg)
c-----------------------------------------------------------------------
c     Purpose: Check for active link in each element
c
c     Inputs:
c        ix(nel)         - Element nodes
c        melnk(ndf,numnp)- Link information at global nodes
c        ndf             - Number of dofs
c        nel             - Number nodes on element
c     Outputs:
c        sflg            - True if links exist
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      logical    sflg
      dimension  ix(nel),melnk(ndf,*)

      do i = 1,nel
        il = ix(i)
        if(il.gt.0) then
          do j = 1,ndf
            if(melnk(j,il).gt.0) then 
              sflg = .true. ! element has linked nodes 
              return
            end if
          end do ! j
        endif ! il
      end do ! i
      end
c
      subroutine formht(ht,g1,g2,nss,neq)
c-----------------------------------------------------------------------
c      Purpose: Modify C array by static condensation on system level
c
c      Inputs:
c        ht(nss,nss)   - C1 
c         g1(neq,nss)  - G1 coupling matrix
c         g2(neq,nss)  - G2 coupling matrix
c         nss          - Number of constraint equations
c         neq          - Number of active equations

c      Outputs:
c        ht(nss,nss)   - C2 = C1-G1^T*X2  with K_aa*X2=G2 
c
c      Comments
c        G = G1=G2
c        K has to be calculated before
c
c-----------------------------------------------------------------------
      USE mdata
      USE ndata
      implicit double precision (a-h,o-z)
      dimension ht(nss,nss),g1(neq,nss),g2(neq,nss)

c.... Loop through columns of G2  
      do mm = 1,nss
c....   K_uu*X(mm)=G2(mm)  
c                   AL          AU          AD     G2       JD
        call dasol (gstiff(nal),gstiff(nau),gstiff,g2(1,mm),jdt12,neq
     +            ,aengy)
      end do ! mm

c.... static condensation of C
      ht = ht-matmul(transpose(g1),g2)

      return
      end
c
      subroutine formtau(tau,g1,u,ddu,id,nss,neq,ndf,numnp)
c-----------------------------------------------------------------------
c      Purpose: Modify S array by static condensation on system level
c
c      Inputs:
c        tau(nss)      - S1 
c         g1(neq,nss)  - G coupling matrix
c         u(*)         - displacement vector 
c         ddu(nneq)    - incremental displacment vector ddu
c         id(ndf,numnp)- equation numbers for each active dof
c         nss          - Number of constraint equations
c         neq          - Number of active equations
c         ndf          - number dof/node
c         numnp        - number of nodes in mesh   

c      Outputs:
c        tau(nss)      - S2 = S1-G^T*ddV_a  with K_aa*ddVa=Ra
c                        K_aa and ddV_a are known from last solution step 
c
c      Comments
c        K has to be calculated before
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      dimension tau(nss),g1(neq,nss),u(*),id(*),ddu(numnp*ndf)
 
      nneq  = numnp*ndf
      nneq2 = nneq+nneq 

c.... incremental displacments 
      call pmovec(id,u(nneq2+1),ddu,nneq)
c      drnorm1=sqrt(dot(u,u,nneq))
c      drnorm2=sqrt(dot(u(nneq+1),u(nneq+1),nneq))
c      drnorm3=sqrt(dot(u(nneq2+1),u(nneq2+1),nneq))
c      drnorm4=sqrt(dot(ddu,ddu,neq))
c      write(16,'(a15,4e12.5)') 
c     +'Norm U,DU,DDU,X',drnorm1,drnorm2,drnorm3,drnorm4

c.... static condensation of Tau
      do i = 1,nss
        do k = 1,neq
          tau(i) = tau(i)-g1(k,i)*ddu(k) 
        end do ! k
      end do ! i  

      return
      end
c
      subroutine epsq_me(prt)
c----------------------------------------------------------------------
c
c      Purpose: read strain values for later use with MACRO EPSQ V = A*E 
c
c      Inputs:
c         prt         - print flag
c
c      Outputs: in epsd.h
c     natyp=1 (3d)        eps=[E_11,E_22,E_33,2E_12,2E_13,2E_23]           
c     natyp=2 (shell)     eps=[E_11,E_22,2E_12,K_11,K_22,2K_12,2E_13,2E_23]
c     natyp=3 (shell ext) eps=[E_11,E_22,2E_12,K_11,K_22,2K_12,2E_13,2E_23,E_33,K_33]
c     natyp=4 (beam)      eps=[E_11,2E_12,2E_13,K_11,K_22,K_33]
c     natyp=5 (beam)      eps=[E_11,2E_12,2E_13,K_11,K_22,K_33,E_220,E_221,E_222,E_330,E_331,E_332,E_230,E_231,E_233]
c 
c----------------------------------------------------------------------
      USE epsdh
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      dimension td(16)

c.... read strains
      if(ior.lt.0) write(*,3000)
      call dinput(td,5)   
c.... iswm,nss,natyp,factx,facty
      iswm  = td(1)
      nss   = td(2)
      natyp = td(3)
      factx = td(4)
      facty = td(5)
c      if(natyp.eq.3.and.nss.ne.10)stop'natyp=3-->nss=10!'  
      call dinput(td,nss)
c.... Get E from input
      do i = 1,nss
        epsd(i)=td(i)
      end do
c...  control input
      if(prt) then  
          write(iow,2000) iswm,nss,natyp,factx,facty,(epsd(i),i=1,nss) 
        if(ior.lt.0) then 
          write(*  ,2000) iswm,nss,natyp,factx,facty,(epsd(i),i=1,nss) 
        end if
      end if  

      return
2000  format(/'  Input of Strains for DISP Generation with EPSQ',/,
     1        '  iswm           = ',i12,/,   
     2        '  nss            = ',i12,/,   
     3        '  natyp          = ',i12,/,   
     4        '  factx          = ',e12.5,/,   
     5        '  facty          = ',e12.5,/,   
     6        '  Strains values = ',15e12.5)
3000  format(' Input: iswm,nsig,natyp,factx,facty,eps(1-6/8/10/15)>',$)
      end

      subroutine epsq_ma(x,f,u,id,melnk,ndm,ndf,numnp,flnk,prt)
c----------------------------------------------------------------------
c
c      Purpose: 1) Set prescribed values for displacements on boundaries
c                  without setting b.c. conditions of a RVE   
c                  V_bc = A_bc*E = Em*X_bc 
c               2) Set values for displacements on linked (slave) nodes for pbc
c                  V_slave = v_master + A_(Dx)*E = Em*Dx,    Dx=X_s-X_m 
c
c     natyp=1 (3d)        eps=[E_11,E_22,E_33,2E_12,2E_13,2E_23]           
c     natyp=2 (shell)     eps=[E_11,E_22,2E_12,K_11,K_22,2K_12,2E_13,2E_23]
c     natyp=3 (shell ext) eps=[E_11,E_22,2E_12,K_11,K_22,2K_12,2E_13,2E_23,E_33,K_33]
c     natyp=4 (beam)      eps=[E_11,2E_12,2E_13,K_11,K_22,K_33]
c     natyp=5 (beam)      eps=[E_11,2E_12,2E_13,K_11,K_22,K_33,E_220,E_221,E_222,E_330,E_331,E_332,E_230,E_231,E_233]

c
c      Inputs:
c         x(ndm,numnp)     - Nodal coordinates of mesh
c         id(ndf,numnp)    - b.c. code for each node
c         melnk(ndf,numnp) - b.c. code for each node
c         ndm              - Spatial dimension of mesh
c         ndf              - Number dof/node
c         numnp            - Number of nodes in mesh
c         flnk             - Flag T=link, F=nolink 
c         prt              - print flag
c
c      Outputs:
c         f(ndf,numnp)     - Load-vector
c         u(ndf,numnp,*    - Displacement-vector
c
c      To be defined:
c         ibin             - 0=binary transfer,else ascii

c      Comments:
c         iepsd            - 1 from input macro, else from file  
c      input macro boun/ebou necessary for setting b.cs. 
c      only b.c. on boundaries allowed! 
c      only link on boundaries allowed! 
c 
c----------------------------------------------------------------------
      USE comfil
      USE epsdh
      USE errchk
      USE iofile
      USE pdata2      ! for idev
      USE shmname
      USE tdata
      implicit double precision (a-h,o-z)
      logical prt
      
      dimension x(ndm,numnp),f(ndf,numnp),u(ndf,numnp,*),id(ndf,numnp),
     +          melnk(ndf,numnp),fl(3),xl(3),inl(3),
     +          emat(3,3),emate(3,3),ematy(3,3),ematz(3,3),
     +          xmin(3),xmax(3)

      character*229 fcisf
      logical flnk
      data bl/-999.d0/
      
#ifdef __INTEL_      
      INTEGER*8 :: emap, nssmap
      CHARACTER*3 :: slgp
      CHARACTER*6 :: sngp
      INTEGER :: n, cur, fin
#endif
      
      double precision, allocatable, dimension(:) :: eps

      data d5/0.5d0/ , d25/0.25d0/ 

c.... name of transfer file  
      iex=ipos(  finp,229)  
      fcisf = finp(1:iex) 
      call dochar2(fcisf,ip)     ! look for i
      call dochar1(fcisf,'f',ip) ! set i -> f

      if (iepsd.eq.1) then
c....   Get E from input macro EPSQ
        skfy = 1.d0
        skfz = 1.d0
        allocate(eps(nss))
        eps = 0.d0
        do i = 1,nss
          eps(i) = epsd(i)
        end do
      else 
        epsd=0.d0 
c....   Get E from macro problem
c...    transfer version
        if(idev.eq.3) ibin=0 ! INTEL 
        if(idev.eq.4) ibin=0 ! SALFORD
        if(ibin.eq.0) then  
          open(33,file=fcisf,form='unformatted',status='unknown')
          rewind(33)
          read(33) dt,ttim,dvp,skfy,skfz,iswm,nss,natyp
          allocate(eps(nss))
          eps=0.d0

          read(33) eps
c...     code only for ibin=0! ibin=1,2 open! 
          read(33) ngp,lgp
          write(iow,*) 'EL, GP',ngp,lgp

          close(33)
        else if(ibin.eq.2) then
          
#ifdef __INTEL_

c     first blank char
cDH          ifb = 1
cDH         do while (fcisf(ifb:ifb) .ne. ' ')
cDH         ifb = ifb + 1
cDH       end do

c     parse lgp number from input filename
cDH          read(fcisf(ifb-2:ifb-1), '(i2)') lgp

cDH   parse LGP from restart filename instead when using parallelization
cDH   and shared memory, as the input file can have an arbitrary number
cDH   just correlating with the thread ID
          
          ifb = 1 ! the point (.) in rrve_01.0001 here
          do while (fres(ifb:ifb).ne.'.')
            ifb = ifb + 1
          end do
          read(fres(ifb-2:ifb-1), '(i2)') lgp
          
cdh   parse ngp (element number) from restart filename
          ifb = 1
          do while (fres(ifb:ifb).ne.' ')
            ifb = ifb + 1
          end do
          read(fres(ifb-4:ifb-1), '(i4)') ngp
        
          shmbasename = 'Global\\FEAPmNNNTTTEEEEEE'C
          
          WRITE(slgp,'(i3.3)') lgp
          shmbasename(13:15) = slgp(1:3)
          WRITE(sngp,'(i6.6)') ngp
          shmbasename(19:24) = sngp(1:6)
          
c         read nss first to be able to allocate eps
          call feap_shmem_OpenNSSMapping(nssmap,shmbasename)
          call feap_shmem_readNSS(nssmap,nss)

          call feap_shmem_CloseHandle(nssmap)
          
          
          allocate(eps(nss))
          eps = 0
            
          call feap_shmem_OpenEMapping(emap,shmbasename)
          call feap_shmem_readE(emap,eps,dt,ttim,dvp,skfy,skfz,iswm,nss,
     &                          natyp)
          call feap_shmem_CloseHandle(emap)
#endif
          
        else 
          open(33,file=fcisf,form='formatted',status='unknown')
          rewind(33)
          read(33,'(5e18.9,3i2)') dt,ttim,dvp,skfy,skfz,iswm,nss,natyp 

          allocate(eps(nss))
          eps=0

          read(33,'(15e18.9)') (eps(i),i=1,nss) 
          close(33)
        end if 
      end if

c.... calculate l_x-l_y-area and volume of RVE
      xmin = 0.d0
      xmax = 0.d0
      dxyz = 0.d0
      
c.... find start node without tie
      do n = 1,numnp 
        if(x(1,n).ne. bl)goto 12  
      end do
      
12    do i = 1,ndm
        xmin(i) = x(i,n)
        xmax(i) = x(i,n)
      end do
      
c.... extreme values
      do n = 1,numnp
        if(x(1,n).ne. bl) then ! tie
          do i = 1,ndm
            xmin(i) = min(xmin(i),x(i,n))
            xmax(i) = max(xmax(i),x(i,n))
          end do ! i
        end if 
      end do ! n
      
c.... length l_x,l_y,l_z
      do i = 1,ndm
        dxyz(i) = xmax(i)-xmin(i) 
      end do

c.... shear correction factor        !! shell
      rlx  = dxyz(1)
      rly  = dxyz(2)
      rlz2 = dxyz(3)*dxyz(3)
      alphax  = 1.d0+factx*rlx*rly/rlz2
      alphay  = 1.d0+facty*rlx*rly/rlz2
      wcappax = skfy
      wcappay = skfz
      skfy    = dsqrt(alphax)*wcappax
      skfz    = dsqrt(alphay)*wcappay
c
c...  control input
      if(prt) then  
       write(iow,2002) ttim,dt,skfy,skfz,iswm,nss,natyp,(eps(i),i=1,nss)
       if(ior.lt.0) then 
       write(*  ,2002) ttim,dt,skfy,skfz,iswm,nss,natyp,(eps(i),i=1,nss) 
       end if
      end if  
c
c...  set matrix E
      if(natyp.eq.1) then ! 3D

        emat(1,1)  = eps(1)
        emat(1,2)  = eps(4)*d5
        emat(1,3)  = eps(5)*d5
  
        emat(2,1)  = eps(4)*d5
        emat(2,2)  = eps(2)
        emat(2,3)  = eps(6)*d5
 
        emat(3,1)  = eps(5)*d5
        emat(3,2)  = eps(6)*d5
        emat(3,3)  = eps(3)

        emate      = emat
       
      else if(natyp.eq.2) then ! Shell
      
        emate = 0
        emate(1,1) = eps(1)
        emate(1,2) = eps(3)*d5
        emate(2,1) = emate(1,2)
        emate(2,2) = eps(2)
        emate(1,3) = eps(7)*skfy
        emate(2,3) = eps(8)*skfz
        
        ematz = 0 
        ematz(1,1) = eps(4)
        ematz(1,2) = eps(6)*d5
        ematz(2,1) = ematz(1,2)
        ematz(2,2) = eps(5)

      else if(natyp.eq.3) then ! Shell ext.
      
        emate = 0
        emate(1,1) = eps(1)
        emate(1,2) = eps(3)*d5
        emate(2,1) = emate(1,2)
        emate(2,2) = eps(2)
        emate(1,3) = eps(7)*d5*skfy
        emate(2,3) = eps(8)*d5*skfz
        emate(3,1) = eps(7)*d5*skfy
        emate(3,2) = eps(8)*d5*skfz
        emate(3,3) = eps(9)
        
        ematz = 0 
        ematz(1,1) = eps(4)
        ematz(1,2) = eps(6)*d5
        ematz(2,1) = ematz(1,2)
        ematz(2,2) = eps(5)
        ematz(3,3) = eps(10)

      else if(natyp.eq.4) then ! Beam

       emate = 0.d0
       ematy = 0.d0
       ematz = 0.d0

       emate(1,1) = eps(1) 
       emate(1,2) = eps(2)*d5 !*skfy 
       emate(2,1) = emate(1,2)
       emate(1,3) = eps(3)*d5 !*skfz 
       emate(3,1) = emate(1,3)
       ematy(1,1) =-eps(6) 
       ematy(1,3) = eps(4)*d5 
       ematy(3,1) = ematy(1,3)
       ematz(1,1) = eps(5) 
       ematz(1,2) =-eps(4)*d5 
       ematz(2,1) = ematz(1,2)

      else if(natyp.eq.5) then ! Beam ext.    

      end if ! natyp

c.... set displacments for b.c.s in f
      do n = 1,numnp
        if(x(1,n).ne. -999.d0) then
          xl(1) = x(1,n)
          xl(2) = x(2,n)
          xl(3) = x(3,n)
          
          if(natyp.eq.2) emat = emate + xl(3)*ematz               ! Shell
          if(natyp.eq.3) emat = emate + xl(3)*ematz               ! Shell ext.
          if(natyp.eq.4) emat = emate + xl(2)*ematy + xl(3)*ematz ! Beam
          
          fl = matmul(emat,xl)
          
          if(natyp.eq.3)then 
           f1 = 0
           f2 = 0 
           f3 = -(xl(1)*xl(1)*eps(4) + xl(2)*xl(2)*eps(5)
     +          + xl(1)*xl(2)*eps(6))*d5
           fl(1) = fl(1) + f1
           fl(2) = fl(2) + f2
           fl(3) = fl(3) + f3
          endif
          
          if(id(1,n).lt. 0) f(1,n) = fl(1)
          if(id(2,n).lt. 0) f(2,n) = fl(2)
          if(id(3,n).lt. 0) f(3,n) = fl(3)
        end if
      end do ! n

c...  Print
      if(prt)              write(iow,2000)        (i,i=1,ndf)
      if(prt.and.ior.lt.0) write(*  ,2000)        (i,i=1,ndf)

      do n = 1,numnp
        if(x(1,n).ne. -999.d0) then
          id3 = ndf
          do i = 1,ndf
            if(id(i,n).gt.0) id3=id3-1 
          end do 
          if(id3.eq.0) go to 10
          if(prt)              write(iow,2001) n,(f(i,n),i=1,ndf)
          if(prt.and.ior.lt.0) write(  *,2001) n,(f(i,n),i=1,ndf)
10        continue      
        end if
      end do

      if(iswm.eq.3.or.iswm.eq.8) then

c....   set displacments for linked nodes with p.b.c. in u
        ndm1 = 3
        if(natyp.eq.2)ndm1 = 2
        if(flnk) then
          do n = 1,numnp
            if(x(1,n).ne. -999.d0) then
              xl(1) = x(1,n)
              xl(2) = x(2,n)
              xl(3) = x(3,n)
              xp = xl(1)
              yp = xl(2)
              zp = xl(3)

              if(natyp.eq.2) emat = emate + xl(3)*ematz               ! Shell   
              if(natyp.eq.3) emat = emate + xl(3)*ematz               ! Shell ext.   
              if(natyp.eq.4) emat = emate + xl(3)*ematz + xl(2)*ematy !Beam
c#########    if(natyp.eq.5) emat = ........ !Beam ext.

c             check link and modify nodal coordinates of slave nodes
c             DX=Xslave=Xslave-Xmaster 

              do i = 1,ndm1
                inl(i) = melnk(i,n) 
                
                if(inl(i).gt.0) then
                  xl(i) = xl(i) - x(i,inl(i))      ! Delta x 
                end if
              end do ! i

              fl = matmul(emat,xl)

             if(natyp.eq.3)then 
              xm = xp-xl(1)
              ym = yp-xl(2)

              f1 = 0
              f2 = 0
              f3 = -d5*((xp*xp-xm*xm)*eps(4)+(yp*yp-ym*ym)*eps(5)
     +                 +(xp*yp-xm*ym)*eps(6))

              fl(1) = fl(1)+f1
              fl(2) = fl(2)+f2
              fl(3) = fl(3)+f3
             endif
 
c             V_slave = V_master + A_(Dx)*E = Em*Dx 
              do i = 1,ndm1
               if(inl(i).gt. 0) u(i,n,1) = u(i,inl(i),1) + fl(i)
c               if(inl(i).gt. 0) u(i,n,2) = u(i,inl(2),2) + fl(i) !open!!!
              enddo

            end if
          end do  

c...      print
          if(prt)              write(iow,2003)        (i,i=1,ndf)
          if(prt.and.ior.lt.0) write(*  ,2003)        (i,i=1,ndf)

          do n = 1,numnp
            if(x(1,n).ne. -999.d0) then
              id3 = ndf
              do i = 1,ndf
                if(melnk(i,n).eq.0) id3=id3-1
              end do 
              if(id3.eq.0) go to 11
              if(prt)              write(iow,2001) n,(u(i,n,1),i=1,ndf)
              if(prt.and.ior.lt.0) write(  *,2001) n,(u(i,n,1),i=1,ndf)
11            continue      
            end if
          end do
        end if
      end if ! iswm

      return
2000  format(/'  EPSQ    n o d a l    prescribed b.c. values'/
     1       4x,'node',9(i3,'-dof ')/(8x,9(i3,'-b.c. val.')))
2001  format(i8,9f8.4/(8x,9f8.4))
2002  format(/'  DISP Generation from Strains with EPSQ',/,
     1        '  Time/DT/KY/KZ/ISWM/NSS/natyp = ',4e12.5,3i4,/,
     1        '  Strains values   = ',15e12.5)
2003  format(/'  EPSQ    n o d a l    prescribed link values'/
     1       4x,'node',9(i3,'-dof ')/(8x,9(i3,'-link val.')))

      end