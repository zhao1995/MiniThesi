      subroutine elmt21(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-----------------------------------------------------------------------
c.....geometrically nonlinear 3 d brick element 
c     with 4/8/20/20=>21/27/64 nodes     
c
c-----------------------------------------------------------------------
c
c.... Material parameters
c     1. Card
c     d(1) = matn   = material number
c     d(2) = lin    = 0 = linear  1 = nonlinear
c     d(3) = isym   = 0 = symmetric  1 = nonsymmetric
c     d(4) = istr   = 1<0>=2.PK  2=Cauchy 
c     d(5) =<kgs    = No. of GP for K> 
c     d(6) =<mgs    = No. of GP for Sigma> 
c
c     ------------------------------------------------------------------
c
c     2. Card loading
c     d(10) = qx     = volume load in x-direction
c     d(11) = qy     = volume load in y-direction
c     d(12) = qz     = volume load in z-direction
c     d(13) = theta  = const. temperature (load)
c     d(14) = at     = alpha_t
c     d(15) = rho    = gamma/g
c
c     ------------------------------------------------------------------
c
c     3. Card material data with respect to chosen mat. number
c        d(16)       = nh  - length of history array at Gauss-Point
c        d(17)-d(..) = ... see SR matelib3d 
c
c-----------------------------------------------------------------------
c     QLOA
c     q(1) = qx     = volume load in x-direction
c     q(2) = qy     = volume load in y-direction
c     q(3) = qz     = volume load in z-direction
c     q(4) = theta  = const. temperature (load)
c
c     ------------------------------------------------------------------
c
c     isw= 3  tangent matrix+residuum 
c     isw= 4  print stresses
c     isw= 5  lumped mass matrix 
c     isw= 6  residuum 
c     isw= 8  plot stresses 
c     isw=11  extended system 
c     isw=15  update history arrays 
c     isw=22  element loads 
c
c-----------------------------------------------------------------------
c     Elements         shape function GP default for all isw   
c      4-nodes  Tetrah lin            4    GP
c      8-nodes  Brick  lin            2x2x2GP 
c     20-nodes  Brick  quad serendip  3x3x3GP
c     21-nodes  Brick  quad serendip  3x3x3GP 
c     27-nodes  Brick  quad lagrange  3x3x3GP  
c     64-nodes  Brick  cube lagrange  4x4x4GP  
c
c     reduced Integration for isw=4/8+elastic material should be used 
c
c-----------------------------------------------------------------------
c     OPEN Tetraeder
c     # Plot: hide geht nicht, damit sind alle Plots praktisch unbrauchbar
c             mit hids geht es 'einigermaßen' bei vielen Elementen,
c             Ausgang über TECPLOT geht!
c     # Lasten müssen genau auf Oberfläche verteilt werden!! 
c     # isw= 5 nicht geprüft
c     # isw=22 nicht eingebaut und nicht geprüft
c     
c     OPEN 64 node hexahedron
c     # Plot: geht i.W. nicht!
c     # Mesh ok, aber nur 8 Eckknoten, siehe ipord64 
c
c-----------------------------------------------------------------------
c     (c)  W. Wagner   11/2005
c          W. Wagner   01/2010 Tetraeder
c          W. Wagner   09/2012 64-node hexagon 
c          W. Wagner   09/2014 QLOA mit Temperatur
c-----------------------------------------------------------------------
      USE bdata
      USE cdat1
      USE cdata
      USE eldata
      USE fe2mat    ! for matfe2=ma
      USE fe2tran   ! for irtyp
      USE hdata
      USE hdatam
      USE iofile
      USE pdata6
      USE prisdat
      USE prlod
      USE qload
      USE strnam
      implicit double precision(a-h,o-z)

      dimension ix(*),xl(ndm,*),tl(*),d(ndd),s(nst,nst),p(nst),
     +          ul(ndf,*),shp(4,64),sg(4,64),dgrad(9),
     +          ql(4),xgp(3),eps(6),sig(6),dmat(6,6),sig0(3),
     +          fi(3,3),ts(6,6),shpm(3,64),tau(6),dd(6,6),plout(10),
     +          ipord4(9),ipord8(17),ipord20(33),ipord21(33),ipord64(17)

      dimension h1(*),h2(*),h3(*)

      data ipord4  /1,2,3,1,4,3,2,4,1/
      data ipord8  /1,2,3,4,1,5,6,2,6,7,3,7,8,4,8,5,1/
      data ipord20 /1,13,2,14,3,15,4,16,1,9,5,17,6,10,2,10,6,18,7,11,3,
     1               11,7,19,8,12,4,12,8,20,5,9,1/
      data ipord21 /1,13,2,14,3,15,4,16,1,9,5,18,6,10,2,10,6,19,7,11,3,
     1               11,7,20,8,12,4,12,8,21,5,9,1/
      data ipord64  /1,4,16,13,1,49,52,4,52,64,16,64,61,13,61,49,1/

      ielno = 21
c.... go to correct array processor
      go to(1,2,3,3,5,3,2,3,2,2,3,2,2,2,3,2,2,2,3,2,2,22), isw
c.... input material properties
1     continue
c.... card 1
      if(ior.lt.0) write(*,1001)
1001  format(' Input:Matn, lin(0/1),isym(0/1),istr,GP,GPS')
      call dinput(d(1),6)
      matn=d(1)  
      lin =d(2) 
      isym=d(3)
      istr=d(4)
      kgs =d(5)
      mgs =d(6)
      if(istr.eq.0) d(4)=1.d0
      if(kgs.lt.0.or.kgs.gt.4) stop 'Wrong Input KGS' 
      if(mgs.lt.0.or.mgs.gt.4) stop 'Wrong Input MGS' 

c.... node numbering for mesh plot in SR pltord and set ngs
      if(nen.eq.4) then  ! nel unknown!!  check if nen is always ok!!
        ngs=4
        inord(ielno) = 9
        do ii = 1,9
          ipord(ii,ielno) = ipord4(ii)
        end do
      else if(nen.eq.8) then
        ngs=2
        inord(ielno) = 17        
        do ii = 1,17
          ipord(ii,ielno) = ipord8(ii)
        end do
      else if(nen.eq.20) then
        ngs=3
        inord(ielno) = 33
        do ii = 1,33
          ipord(ii,ielno) = ipord20(ii)
        end do
      else if(nen.eq.21.or.nen.eq.27) then
        ngs=3
        inord(ielno) = 33
        do ii = 1,33
          ipord(ii,ielno) = ipord21(ii)
        end do
      else if(nen.eq.64) then
        ngs=4
        inord(ielno) = 17        
        do ii = 1,17
          ipord(ii,ielno) = ipord64(ii)
        end do
      end if
      if(kgs.eq.0) then ! GP for K
        kgs=ngs
        d(5)=kgs
      end if  
      if(mgs.eq.0) then ! Gp for Sigma
        mgs=ngs
        d(6)=mgs
      end if 
                   write(iow,1002) nen,(d(i),i=1,6)
      if(ior.lt.0) write(*  ,1002) nen,(d(i),i=1,6)
1002  format(5x,'Materialdata for 3D-',i2,'-node brick element',/,
     +  5x,'General data',/,
     +  5x,'Material type ...........',f12.5,/,
     +  5x,'lin=0,nonlin=1 ..........',f12.5,/,
     +  5x,'sym=0,nonsym=1 ..........',f12.5,/,
     +  5x,'istr=1 - 2.PK stresses...',/,
     +  5x,'istr=2 - Cauchy stresses.',f12.5,/,
     +  5x,'No.of GP for K...........',f12.5,/,
     +  5x,'No.of GP for Sigma.......',f12.5,/)

c.... card 2
      if(ior.lt.0) write(*,1003)
1003  format(' Input: qx, qy, qz, delta T, alpah_T, rho')
      call dinput(d(10),6)
                   write(iow,1004) (d(i),i=10,15)
      if(ior.lt.0) write(*  ,1004) (d(i),i=10,15)
1004  format(5x,'Loading data',/,
     +  5x,'qx ......................',g12.5,/,
     +  5x,'qy ......................',g12.5,/,
     +  5x,'qz ......................',g12.5,/,
     +  5x,'delta T .................',g12.5,/,
     +  5x,'alpha T .................',g12.5,/,
     +  5x,'rho=gamma/g..............',g12.5,/)

c.... card 3
      matfe2=ma  
      call matelib3d(h1,h2,nh,d(17),mdd,eps,sig,dmat,6,1,plout,
     +               xgp,tgp,dvp,1.d0,1.d0,1.d0,n,l,1,1,matn,isw)

      d(16) = nh 
c
c.... check integration order in case of plasticity etc   
cwwc     (check actual version of matelib3D!!)
cww      if(matn.eq.4.or.matn.eq.5.or.matn.eq.7.or.matn.eq.8
cww     +     .or.matn.eq.10) then       

      if((nh.ne.0).and.(kgs.ne.mgs)) then 
        call drawmess('ELMT21: ngs_K .ne. ngs_S!',1,0)
        stop
      end if

c.... calculate length of h-array
      if(nen.gt.4) kgs=kgs*kgs*kgs ! then kgs is per direction 
      nh1=kgs*nh 

c...  check length of d-array 
      nmax = 16 + mdd
      if(nmax.gt.ndd) then
         write(*,1005) ndd,nmax
1005     format(1x,'darray is set to  ',i4,' values',/,
     1          1x,'darray needs      ',i4,' values')
         stop
      end if

c.... description of stresses  
      strsus( 1) = '  STRESS S_xx  '
      strsus( 2) = '  STRESS S_yy  '
      strsus( 3) = '  STRESS S_zz  '
      strsus( 4) = '  STRESS T_xy  '
      strsus( 5) = '  STRESS T_xz  '
      strsus( 6) = '  STRESS T_yz  '
      strsus( 7) = '  STRESS S_11  '
      strsus( 8) = '  STRESS S_22  '
      strsus( 9) = '  STRESS S_33  '
      strsus(10) = '               '
      strsus(11) = '               '
      strsus(12) = '               '
      strsus(13) = '               '
      strsus(14) = '               '
      strsus(15) = '               '
      strsus(16) = ' INT.VAR.(1)   '
      strsus(17) = ' INT.VAR.(2)   '
      strsus(18) = ' INT.VAR.(3)   '
      strsus(19) = ' INT.VAR.(4)   '
      strsus(20) = ' INT.VAR.(5)   '
      strsus(21) = ' INT.VAR.(6)   '
      strsus(22) = ' INT.VAR.(7)   '
      strsus(23) = ' INT.VAR.(8)   '
      strsus(24) = ' INT.VAR.(9)   '
      strsus(25) = ' INT.VAR.(10)  '
c...  names for principal moments
      nptyp = 4 
c
2     return

3     matn = d(1)  
      lin  = d(2)
      isym = d(3) 
      istr = d(4)
      nh   = d(16) 
      l    = d(5) 
      if(isw.eq.4.or.isw.eq.8) l=d(6)

      if(nel.eq.4) then
        call int3dt(l,lint,sg)
      else
        call int3d(l,lint,sg)
c...    special: GP=corner nodes, without mid nodes (Tau too large!!)
        if(isw.eq.8.and.(nel.eq.20 .or. nel.eq.21)) then ! 
c          sg=0 
          call gausnod21(sg,lint) 
        end if
      end if
      
c.... set loads from qloa/mate
32    ql=0
      call qload21(d,aqloa,ql,numel,n,mqloa,propq,prop,isw)


c.... for FE^2 (matn=8) extract restart files from FRESG (only ibin=0)
      if(matn.eq.8.and.irtyp.eq.1) call matt3d08(1,n,numel,lint)

c.....loop over gauss points
      nn = 1 ! counter h-array 

      do l = 1,lint
c.....  shape functions and derivatives
        if(nel.eq.4) then
          call shp3dt(sg(1,l),xsj,shp,xl,ndm,n)
        else 
          call shp3d21(sg(1,l),shp,nel)
          call jaco21(shp,xsj,xl,n,nel,ndm,iow)
        end if
c....   integration
        dvp = xsj*sg(4,l)
c....   coordinates, temperature, gradients, transformation matrix ts
        call grad21(shp,xl,ul,tl,ndf,ndm,nel,xgp,tgp,dgrad,fi,ts,lin)
c.....  strains and stresses
        call strain21(d,dgrad,ql,eps,lin)
        matfe2=ma  
        call matelib3d(h1(nn),h2(nn),nh,d(17),mdd,eps,sig,dmat,6,1,
     +                plout,xgp,tgp,dvp,1.d0,1.d0,1.d0,n,l,1,1,matn,isw)
        if(isw.eq.4)  goto 4  ! stre,prin
        if(isw.eq.8)  goto 8  ! stre,plot
        if(isw.eq.15.or.isw.eq.19) goto 30 ! updh
        
        tau = matmul(ts,sig)
        dd  = matmul(ts,matmul(dmat,transpose(ts)))
c....   stiffness matrix and residual
        i1 = 1
        do ii=1,nel
          i2 = i1 + 1
          i3 = i1 + 2
c....     spatial derivatives of shape functions
          do k = 1,3
           shpm(k,ii) = dot(fi(1,k),shp(1,ii),3)
          enddo 
          xn = shpm(1,ii)*dvp
          yn = shpm(2,ii)*dvp
          zn = shpm(3,ii)*dvp
c....     residual G = P - Bt*S
          f = shp(4,ii)*dvp
          p(i1) = p(i1) - (xn*tau(1)+yn*tau(4)+zn*tau(5)) + ql(1)*f 
          p(i2) = p(i2) - (yn*tau(2)+xn*tau(4)+zn*tau(6)) + ql(2)*f 
          p(i3) = p(i3) - (zn*tau(3)+xn*tau(5)+yn*tau(6)) + ql(3)*f 
c....     B^t*C
          a11 = dd(1,1)*xn + dd(4,1)*yn + dd(5,1)*zn
          a12 = dd(1,2)*xn + dd(4,2)*yn + dd(5,2)*zn
          a13 = dd(1,3)*xn + dd(4,3)*yn + dd(5,3)*zn
          a14 = dd(1,4)*xn + dd(4,4)*yn + dd(5,4)*zn
          a15 = dd(1,5)*xn + dd(4,5)*yn + dd(5,5)*zn
          a16 = dd(1,6)*xn + dd(4,6)*yn + dd(5,6)*zn
          
          a21 = dd(4,1)*xn + dd(2,1)*yn + dd(6,1)*zn
          a22 = dd(4,2)*xn + dd(2,2)*yn + dd(6,2)*zn
          a23 = dd(4,3)*xn + dd(2,3)*yn + dd(6,3)*zn
          a24 = dd(4,4)*xn + dd(2,4)*yn + dd(6,4)*zn
          a25 = dd(4,5)*xn + dd(2,5)*yn + dd(6,5)*zn
          a26 = dd(4,6)*xn + dd(2,6)*yn + dd(6,6)*zn
          
          a31 = dd(5,1)*xn + dd(6,1)*yn + dd(3,1)*zn
          a32 = dd(5,2)*xn + dd(6,2)*yn + dd(3,2)*zn
          a33 = dd(5,3)*xn + dd(6,3)*yn + dd(3,3)*zn
          a34 = dd(5,4)*xn + dd(6,4)*yn + dd(3,4)*zn
          a35 = dd(5,5)*xn + dd(6,5)*yn + dd(3,5)*zn
          a36 = dd(5,6)*xn + dd(6,6)*yn + dd(3,6)*zn 

          if(isw.eq.6) goto 33
c....     tangent stiffness matrix
          j1 = 1
          je = nel
          if(isym.eq.0)je = ii
          do jj = 1,je
            j2 = j1 + 1
            j3 = j1 + 2
            if(isym.eq.1)then
             do k = 1,3
              shpm(k,jj) = dot(fi(1,k),shp(1,jj),3)
             enddo 
            endif 
            if(lin.eq.0)then 
             f = 0.d0
            else
             f = (shp(1,ii)*sig(1)*shp(1,jj)+shp(1,ii)*sig(4)*shp(2,jj)
     +         +  shp(2,ii)*sig(2)*shp(2,jj)+shp(2,ii)*sig(4)*shp(1,jj)
     +         +  shp(3,ii)*sig(3)*shp(3,jj)+shp(1,ii)*sig(5)*shp(3,jj)
     +         +  shp(3,ii)*sig(5)*shp(1,jj)+shp(2,ii)*sig(6)*shp(3,jj)
     +         +  shp(3,ii)*sig(6)*shp(2,jj))*dvp
            endif
           
            xn  = shpm(1,jj)
            yn  = shpm(2,jj)
            zn  = shpm(3,jj)
            s(i1,j1) = s(i1,j1) + xn*a11 + yn*a14 + zn*a15 + f 
            s(i2,j1) = s(i2,j1) + xn*a21 + yn*a24 + zn*a25
            s(i3,j1) = s(i3,j1) + xn*a31 + yn*a34 + zn*a35
            s(i1,j2) = s(i1,j2) + yn*a12 + xn*a14 + zn*a16
            s(i2,j2) = s(i2,j2) + yn*a22 + xn*a24 + zn*a26 + f
            s(i3,j2) = s(i3,j2) + yn*a32 + xn*a34 + zn*a36
            s(i1,j3) = s(i1,j3) + zn*a13 + xn*a15 + yn*a16 
            s(i2,j3) = s(i2,j3) + zn*a23 + xn*a25 + yn*a26 
            s(i3,j3) = s(i3,j3) + zn*a33 + xn*a35 + yn*a36 + f


            j1 = j1 + ndf
          end do ! jj
33        i1 = i1 + ndf
        end do ! ii
        goto 30
        
c....   Output Stresses
4       call stress21(dgrad,sig,sig0,isw,istr)
        mct = mct - 1
        if(mct.gt.0) goto 41
        if(istr.eq.1) then
                       write(iow,4001)o,head
          if(ior.lt.0) write(*  ,4001)o,head
        else if(istr.eq.2) then
                       write(iow,4002)o,head
          if(ior.lt.0) write(*  ,4002)o,head
        end if
                     write(iow,4004)
        if(ior.lt.0) write(*  ,4005)
        mct = 50
41                 write(iow,4006) n,ma,(xgp(i),i=1,3),(sig(i),i=1,6)
        if(ior.lt.0) write(*,4006) n,ma,(sig(i),i=1,6)
        goto 30
c....   plot stresses
8       istv = -25
        if(iplma(ma).eq.0) return ! only if MATN
        call stress21(dgrad,sig,sig0,isw,istr)
        call stre21(ix,strea,strea(1+numnp),plout,sig,sig0,shp,nel,
     +                numnp,dvp)
   30   nn = nn + nh   

      end do ! l

c.... for FE^2 (matn=8) save restart files in FRESG (only ibin=0)
      if(matn.eq.8.and.irtyp.eq.1) then
        if(hflgu.and.h3flgu) call matt3d08(2,n,numel,lint)
      end if   

      if(isw.eq.4.or.isw.eq.6.or.isw.eq.8) return
c.....upper part of stiffness matrix in case of symmetry
      if(isym.eq.0)call msym(s,nst,nst,1)
c      if(n.le.1)call mprint(s,6,6,nst,'s21 ')
c      if(n.le.1)call mprint(p,6,1,nst,'p21 ')
c
c.... calculate eigenvalues
cww   call elemev(s,nel,ndf,nst,1)
c
      return
c
c.... mass matrix lumped
5     rho = d(15)
c.... integration points  
      l    = d(5) 

      if(nel.eq.4) then
        call int3dt(l,lint,sg)
      else
        call int3d(l,lint,sg)
      end if

c.....loop over gauss points
      do l = 1,lint
c.....  shape functions and derivatives
        if(nel.eq. 4) then
          call shp3dt(sg(1,l),xsj,shp,xl,ndm,n)
        else 
          call shp3d21(sg(1,l),shp,nel)
          call jaco21(shp,xsj,xl,n,nel,ndm,iow)
        end if
c.....  jacobian
        dvp = xsj*sg(4,l)
c....   loop over nodes
        i1=0
        do ii=1,nel
          add = rho*shp(4,ii)*dvp
c....     lumped mass matrix           
          p(i1+1) = p(i1+1) + add
          p(i1+2) = p(i1+2) + add
          p(i1+3) = p(i1+3) + add
          i1 = i1 + ndf
        end do     
      end do
      return
c.... load vector for QLOA
22    matn = d(1)  
      lin  = d(2)
      isym = d(3) 
      istr = d(4)
      nh   = d(16) 
      l    = d(5) 

      if(nel.eq.4) then
        call int3dt(l,lint,sg)
      else
        call int3d(l,lint,sg)
      end if

c.... set loads from qloa
      ql=0
      call qload21(d,aqloa,ql,numel,n,mqloa,propq,prop,isw)

c.....loop over gauss points
      nn = 1 ! counter h-array 
      do l = 1,lint

c.....  shape functions and derivatives
        if(nel.eq.4) then
          call shp3dt(sg(1,l),xsj,shp,xl,ndm,n)
        else 
          call shp3d21(sg(1,l),shp,nel)
          call jaco21(shp,xsj,xl,n,nel,ndm,iow)
        end if

c....   integration
        dvp = xsj*sg(4,l)

c....   coordinates, temperature, gradients, transformation matrix ts
        call grad21(shp,xl,ul,tl,ndf,ndm,nel,xgp,tgp,dgrad,fi,ts,lin)

c.....  temperature strains
        eps=0.d0
        if(ql(4).ne.0.d0) then
          etheta = d(14)*ql(4)
          eps(1) = etheta
          eps(2) = etheta
          eps(3) = etheta
        end if

c....   only D=DMAT 
        matfe2=ma  
        call matelib3d(h1(nn),h2(nn),nh,d(17),mdd,eps,sig,dmat,6,1,
     +                plout,xgp,tgp,dvp,1.d0,1.d0,1.d0,n,l,1,1,matn,isw)

c.....  temperature stresses
        sig = matmul(dmat,eps)
        tau = matmul(ts,sig)

c....   stiffness matrix and residual
        i1 = 1
        do ii=1,nel
          i2 = i1 + 1
          i3 = i1 + 2
c....     spatial derivatives of shape functions
          do k = 1,3
           shpm(k,ii) = dot(fi(1,k),shp(1,ii),3)
          end do 
          xn = shpm(1,ii)*dvp
          yn = shpm(2,ii)*dvp
          zn = shpm(3,ii)*dvp

c....     load = P + B^T*S_T
          f = shp(4,ii)*dvp
          p(i1) = p(i1) + (xn*tau(1)+yn*tau(4)+zn*tau(5)) + ql(1)*f 
          p(i2) = p(i2) + (yn*tau(2)+xn*tau(4)+zn*tau(6)) + ql(2)*f 
          p(i3) = p(i3) + (zn*tau(3)+xn*tau(5)+yn*tau(6)) + ql(3)*f 

          if(ql(4).ne.0.d0) then
c....      tangent stiffness matrix KG(SigmaT)
           if(lin.ne.0)then 
            j1 = 1
            je = nel
            if(isym.eq.0)je = ii
            do jj = 1,je
             j2 = j1 + 1
             j3 = j1 + 2
             f = (shp(1,ii)*sig(1)*shp(1,jj)+shp(2,ii)*sig(2)*shp(2,jj)
     +         +  shp(3,ii)*sig(3)*shp(3,jj))*dvp
           
             s(i1,j1) = s(i1,j1)  - f 
             s(i2,j2) = s(i2,j2)  - f
             s(i3,j3) = s(i3,j3)  - f

             j1 = j1 + ndf
            end do ! jj
           end if ! lin 
          end if ! ql(4) 
          i1 = i1 + ndf
        end do ! ii
        nn = nn + nh   
      end do ! l

      if(isym.eq.0)call msym(s,nst,nst,1)

c
4001      format(a1,20a4,/,2x,'Element Stresses (2.P-K)',/)
4002      format(a1,20a4,/,2x,'Element Stresses (Cauchy)',/)

4004      format(1x,'  EL',1x,' Ma',1x,'   1-COR    ',
     +            1x,'   2-COR    ',1x,'   3-COR    ',
     +            1x,'   S_11     ',1x,'   S_22     ',
     +            1x,'   S_33     ',1x,'   T_xy     ',
     +            1x,'   T_xz     ',1x,'   T_yz     ')
4005      format(1x,'  EL',1x,' Ma',
     +            1x,'   S_11     ',1x,'   S_22     ',
     +            1x,'   S_33     ',1x,'   T_xy     ',
     +            1x,'   T_xz     ',1x,'   T_yz     ')
4006      format(1x,i4,1x,i3,9(1x,g12.5))
      end
c
      subroutine grad21
     +          (shp,xl,ul,tl,ndf,ndm,nel,xgp,tgp,dgrad,fi,ts,lin)
c-----------------------------------------------------------------------
c
c      Purpose: calculate coordinates, temperature, gradients, matrix ts
c
c      Inputs:
c         shp(4,*)  - Shape functions and derivatives at point
c         xl(ndm,*) - nodal coordnates 
c         ul(ndf,*) - nodal dsplacements 
c         tl(*)     - nodal temperatures
c         ndf       - No. of Dofs at node
c         nel       - No. of nodes of element
c         lin       - linear/nonlinear 
c
c      Outputs:
c         xgp       - coordinates of integration point
c         tgp       - temperature at integration point
c         dgrad(9)  - Displacement gradients
c                     [u,x; v,x; w,x; u,y; v,y; w,y; u,z; v,z; w,z]^T
c         fi(3,3)   -  inverse deformation gradient
c         ts(6,6)   -  transformation matrix
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension shp(4,*),dgrad(9),xl(ndm,*),ul(ndf,*),tl(*),xgp(*),
     +          fi(3,3),ts(6,6)
c     
c.... coordinate/temp of point
      tgp = 0.0d0
      do idm = 1,3
        xgp(idm) = 0.0d0
        do inode = 1,nel
          xgp(idm) = xgp(idm) + xl(idm,inode)* shp(4,inode)
          tgp      = tgp      + tl(    inode)* shp(4,inode)
        end do ! inode
      end do ! idm  
c     displacement gradient 
      call pzero(dgrad,9)
      do i = 1,nel
        dndx = shp(1,i)
        dndy = shp(2,i)
        dndz = shp(3,i)
        do ikx = 1,3
          dgrad(ikx)   = dgrad(ikx  )+dndx*ul(ikx,i)
          dgrad(ikx+3) = dgrad(ikx+3)+dndy*ul(ikx,i)
          dgrad(ikx+6) = dgrad(ikx+6)+dndz*ul(ikx,i)
        end do
      end do
c     deformation gradient      
      if(lin.eq.0)then
       do i = 1,3
        do j= 1,3
         fi(i,j) = 0.d0
        enddo
        fi(i,i) = 1.d0
       enddo 
      else
       fi(1,1) = 1.d0 + dgrad(1)
       fi(2,1) =        dgrad(2)
       fi(3,1) =        dgrad(3)
       fi(1,2) =        dgrad(4)
       fi(2,2) = 1.d0 + dgrad(5)
       fi(3,2) =        dgrad(6)
       fi(1,3) =        dgrad(7)
       fi(2,3) =        dgrad(8)
       fi(3,3) = 1.d0 + dgrad(9)
      endif
c     transformation matrix ts:     tau = ts * sig      
      ts(1,1) =       fi(1,1)*fi(1,1)
      ts(1,2) =       fi(1,2)*fi(1,2)
      ts(1,3) =       fi(1,3)*fi(1,3)
      ts(1,4) = 2.0d0*fi(1,1)*fi(1,2)
      ts(1,5) = 2.0d0*fi(1,1)*fi(1,3)
      ts(1,6) = 2.0d0*fi(1,2)*fi(1,3)
      ts(2,1) =       fi(2,1)*fi(2,1)
      ts(2,2) =       fi(2,2)*fi(2,2)
      ts(2,3) =       fi(2,3)*fi(2,3)
      ts(2,4) = 2.0d0*fi(2,1)*fi(2,2)
      ts(2,5) = 2.0d0*fi(2,1)*fi(2,3)
      ts(2,6) = 2.0d0*fi(2,2)*fi(2,3)
      ts(3,1) =       fi(3,1)*fi(3,1)
      ts(3,2) =       fi(3,2)*fi(3,2)
      ts(3,3) =       fi(3,3)*fi(3,3) 
      ts(3,4) = 2.0d0*fi(3,1)*fi(3,2)
      ts(3,5) = 2.0d0*fi(3,1)*fi(3,3)
      ts(3,6) = 2.0d0*fi(3,2)*fi(3,3)
      ts(4,1) =       fi(1,1)*fi(2,1)
      ts(4,2) =       fi(1,2)*fi(2,2)
      ts(4,3) =       fi(1,3)*fi(2,3)
      ts(4,4) =       fi(1,1)*fi(2,2) + fi(1,2)*fi(2,1)
      ts(4,5) =       fi(1,1)*fi(2,3) + fi(1,3)*fi(2,1)
      ts(4,6) =       fi(1,2)*fi(2,3) + fi(1,3)*fi(2,2)
      ts(5,1) =       fi(1,1)*fi(3,1)
      ts(5,2) =       fi(1,2)*fi(3,2)
      ts(5,3) =       fi(1,3)*fi(3,3)
      ts(5,4) =       fi(1,1)*fi(3,2) + fi(1,2)*fi(3,1)
      ts(5,5) =       fi(1,1)*fi(3,3) + fi(1,3)*fi(3,1)
      ts(5,6) =       fi(1,2)*fi(3,3) + fi(1,3)*fi(3,2)
      ts(6,1) =       fi(2,1)*fi(3,1)
      ts(6,2) =       fi(2,2)*fi(3,2)
      ts(6,3) =       fi(2,3)*fi(3,3)
      ts(6,4) =       fi(2,1)*fi(3,2) + fi(2,2)*fi(3,1)
      ts(6,5) =       fi(2,1)*fi(3,3) + fi(2,3)*fi(3,1)
      ts(6,6) =       fi(2,2)*fi(3,3) + fi(2,3)*fi(3,2)
c     inverse deformation gradient  
      call invert(fi,3,3)   
c         
      return
      end
c
      subroutine strain21(d,dgrad,ql,eplan,lin)
c-----------------------------------------------------------------------
c
c      Purpose: calculate strains for geom. nonlinear 3d element
c
c      Inputs:
c         dgrad(9)  - Displacement gradients
c         d         - material data
c         ql        - load data
c         lin       - lin/nonlin 
c
c      Outputs:
c         eplan(6)  - Strains
c
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension d(*),dgrad(9),eplan(6),ql(4)

c     gradients
      uxkx=dgrad(1)
      uykx=dgrad(2)
      uzkx=dgrad(3)
      uxky=dgrad(4)
      uyky=dgrad(5)
      uzky=dgrad(6)
      uxkz=dgrad(7)
      uykz=dgrad(8)
      uzkz=dgrad(9)

c.....strains (linear)
      eplan(1) = uxkx
      eplan(2) = uyky
      eplan(3) = uzkz
      eplan(4) = uxky + uykx
      eplan(5) = uxkz + uzkx
      eplan(6) = uykz + uzky

c.....strains (thermal loading)  d(14)=alpha_t
      if(d(14).ne.0.d0) then
        etheta   = d(14)*ql(4)
        eplan(1) = eplan(1) - etheta
        eplan(2) = eplan(2) - etheta
        eplan(3) = eplan(3) - etheta
      end if
      
      if(lin.eq.0) return

c.....strains (nonlinear)
      eplan(1) = eplan(1) + 0.5d0 *(uxkx*uxkx + uykx*uykx + uzkx*uzkx)
      eplan(2) = eplan(2) + 0.5d0 *(uxky*uxky + uyky*uyky + uzky*uzky)
      eplan(3) = eplan(3) + 0.5d0 *(uxkz*uxkz + uykz*uykz + uzkz*uzkz)
      eplan(4) = eplan(4) +        (uxkx*uxky + uykx*uyky + uzkx*uzky)
      eplan(5) = eplan(5) +        (uxkx*uxkz + uykx*uykz + uzkx*uzkz)
      eplan(6) = eplan(6) +        (uxky*uxkz + uykz*uyky + uzkz*uzky)

      return 
      end
c
      subroutine stress21(dgrad,sig,sig0,isw,istr)
c-----------------------------------------------------------------------
c
c      Purpose: calculate diff. stresses for geom. nonlinear 3d element
c
c      Inputs:
c         dgrad(9)  - Displacement gradients
c         sig(6)    - 2.PK-stresses
c         istr      - type of stress
c
c      Outputs:
c         sig(6)    - new stresses
c         sig0(6)   - new main stresses
c
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension dgrad(9),sig(6),stre(3,3),
     +          sig0(3),ssig(3,3),psig(3,3),tsig(3,3),f(3,3)


      if(isw.eq.4 .or. isw.eq.8) then
c.....  main stresses (copy)
        stre(1,1) = sig(1)
        stre(2,2) = sig(2)
        stre(3,3) = sig(3)
        stre(1,2) = sig(4)
        stre(1,3) = sig(5)
        stre(2,3) = sig(6)
        stre(2,1) = stre(1,2)
        stre(3,1) = stre(1,3)
        stre(3,2) = stre(2,3)

        if(istr.lt.2) goto 21

c....   2.PK stress tensor    
        ssig(1,1)= sig(1)
        ssig(1,2)= sig(4)
        ssig(1,3)= sig(5)
        ssig(2,1)=ssig(1,2)
        ssig(2,2)= sig(2)
        ssig(2,3)= sig(6)
        ssig(3,1)=ssig(1,3)
        ssig(3,2)=ssig(2,3)
        ssig(3,3)= sig(3)
        
c.....  material deformation gradient F
        f(1,1)=1.0d0+dgrad(1)
        f(1,2)=dgrad(4)
        f(1,3)=dgrad(7)
        f(2,1)=dgrad(2)
        f(2,2)=1.0d0+dgrad(5)
        f(2,3)=dgrad(8)
        f(3,1)=dgrad(3)
        f(3,2)=dgrad(6)
        f(3,3)=1.0d0+dgrad(9)
        
c.....  1.PK stress tensor
        do 13 i=1,3
          do 13 k=1,3
            psig(i,k)=0.0d0
            do 13 j=1,3
  13          psig(i,k)=psig(i,k)+f(i,j)*ssig(j,k)
        
c.....  Cauchy-stresses(2)
        djacb1 = f(1,1)*(f(2,2)*f(3,3)-f(2,3)*f(3,2))
        djacb2 =-f(1,2)*(f(2,1)*f(3,3)-f(2,3)*f(3,1))
        djacb3 = f(1,3)*(f(2,1)*f(3,2)-f(2,2)*f(3,1))
        detf= djacb1 + djacb2 + djacb3
        do 14 i=1,3
          do 14 k=1,3
            tsig(i,k)=0.0d0
            do 14 j=1,3
  14           tsig(i,k)=tsig(i,k)+psig(i,j)*f(k,j)/detf

c....   stresses (copy)
        sig(1) = tsig(1,1)
        sig(2) = tsig(2,2)
        sig(3) = tsig(3,3)
        sig(4) = tsig(1,2)
        sig(5) = tsig(1,3)
        sig(6) = tsig(2,3)

c....   main stresses (copy)
        do i=1,3
          do k=1,3
            stre(i,k)= tsig(i,k)
          end do
        end do   
c.....  main stresses
21      call stress3d21(stre,sig0)
      end if

      return 
      end
c
      subroutine jaco21(shp,xsj,xl,n,nel,ndm,iow)
c-----------------------------------------------------------------------
c.....Jacobi-Matrix and  cartesian derivatives of shape functions
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension xl(ndm,*),shp(4,nel),xjaci(3,3),xjacm(3,3),cartd(3,64)
      dimension a(3,3)
c.....Jacobi-matrix xjacm
      do 4 idm = 1,3
        do 4 jdm = 1,3
          xjacm(idm,jdm) = 0.0d0
          xjaci(idm,jdm) = 0.0d0
          do 4 inode = 1,nel
4     xjacm(idm,jdm)=xjacm(idm,jdm)+shp(idm,inode)*xl(jdm,inode)
c
c.....Determinant
      djacb1 = xjacm(1,1)*(xjacm(2,2)*xjacm(3,3)-xjacm(2,3)*xjacm(3,2))
      djacb2 =-xjacm(1,2)*(xjacm(2,1)*xjacm(3,3)-xjacm(2,3)*xjacm(3,1))
      djacb3 = xjacm(1,3)*(xjacm(2,1)*xjacm(3,2)-xjacm(2,2)*xjacm(3,1))
      xsj= djacb1 + djacb2 + djacb3
      if(xsj) 6,6,8
    6 write(iow,600) n
      stop
    8 continue
c.....Inverse of Jacobian
c.....Unterdeterminanten
      do 11 idm = 1,3
        do 11 jdm = 1,3
 11     a(idm,jdm) = 0.0d0
      a(1,1) =  xjacm(2,2)*xjacm(3,3)-xjacm(2,3)*xjacm(3,2)
      a(1,2) =-(xjacm(2,1)*xjacm(3,3)-xjacm(2,3)*xjacm(3,1))
      a(1,3) =  xjacm(2,1)*xjacm(3,2)-xjacm(2,2)*xjacm(3,1)
      a(2,1) =-(xjacm(1,2)*xjacm(3,3)-xjacm(1,3)*xjacm(3,2))
      a(2,2) =  xjacm(1,1)*xjacm(3,3)-xjacm(1,3)*xjacm(3,1)
      a(2,3) =-(xjacm(1,1)*xjacm(3,2)-xjacm(1,2)*xjacm(3,1))
      a(3,1) =  xjacm(1,2)*xjacm(2,3)-xjacm(1,3)*xjacm(2,2)
      a(3,2) =-(xjacm(1,1)*xjacm(2,3)-xjacm(1,3)*xjacm(2,1))
      a(3,3) =  xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1)
      xjaci(1,1)= a(1,1)/xsj
      xjaci(1,2)= a(2,1)/xsj
      xjaci(1,3)= a(3,1)/xsj
      xjaci(2,1)= a(1,2)/xsj
      xjaci(2,2)= a(2,2)/xsj
      xjaci(2,3)= a(3,2)/xsj
      xjaci(3,1)= a(1,3)/xsj
      xjaci(3,2)= a(2,3)/xsj
      xjaci(3,3)= a(3,3)/xsj
c
c.....cartesian derivatives
      do 10 idm = 1,3
        do 10 i = 1,nel
         cartd(idm,i) = 0.0d0
         do 10 jdm = 1,3
10    cartd(idm,i)=cartd(idm,i)+xjaci(idm,jdm)*shp(jdm,i)
      do 20 idm=1,3
        do 20 i=1,nel
20    shp(idm,i) = cartd(idm,i)
c
  600 format(1x,'program stop in jaco3d of Elmt21',/,1x,
     +'zero or negative area for element number',i5)
      return
      end
c
      subroutine shp3d21(sg,shp,nel)
c-----------------------------------------------------------------------
c
c      Purpose: Compute 3-d isoparametric 8/20/21/27/64-node shape 
c               functions and their derivatives wrt. r,s,t
c
c      Inputs:
c         sg(3)     - Natural coordinates of point
c         nel       - No. of nodes of element
c
c      Outputs:
c         shp(4,*)  - Shape functions and derivatives at point
c                     shp(1,i) = dN_i/dr
c                     shp(2,i) = dN_i/ds
c                     shp(3,i) = dN_i/dt
c                     shp(4,i) =  N_i
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension shp(4,64),ri(27),si(27),ti(27),sg(3)
      data ri /-1, 1, 1,-1, -1, 1, 1,-1,  -1, 1, 1,-1,   0, 1, 0,-1, 0,
     1                                  0, 1, 0,-1, 0,   0, 1, 0,-1, 0/
      data si /-1,-1, 1, 1, -1,-1, 1, 1,  -1,-1, 1, 1,  -1, 0, 1, 0, 0,
     1                                 -1, 0, 1, 0, 0,  -1, 0, 1, 0, 0/
      data ti /-1,-1,-1,-1,  1, 1, 1, 1,   0, 0, 0, 0,  -1,-1,-1,-1,-1,
     1                                  1, 1, 1, 1, 1,   0, 0, 0, 0, 0/
      r=sg(1)
      s=sg(2)
      t=sg(3)

      d1=1.0d0
      d2=1.0d0/2.0d0
      d8=1.0d0/8.0d0
      d4=1.0d0/4.0d0
      dz=2.0d0
      rr=r*r
      ss=s*s
      tt=t*t
      call pzero(shp,256)
      if(nel.eq.8) then
c....   shape functions for 8-node element
        shp(4,1) = d8*(d1-r)*(d1-s)*(d1-t)
        shp(4,2) = d8*(d1+r)*(d1-s)*(d1-t)
        shp(4,3) = d8*(d1+r)*(d1+s)*(d1-t)
        shp(4,4) = d8*(d1-r)*(d1+s)*(d1-t)
        shp(4,5) = d8*(d1-r)*(d1-s)*(d1+t)
        shp(4,6) = d8*(d1+r)*(d1-s)*(d1+t)
        shp(4,7) = d8*(d1+r)*(d1+s)*(d1+t)
        shp(4,8) = d8*(d1-r)*(d1+s)*(d1+t)
c....   derivatives
        shp(1,1) = -d8*(d1-s)*(d1-t)
        shp(1,2) =  d8*(d1-s)*(d1-t)
        shp(1,3) =  d8*(d1+s)*(d1-t)
        shp(1,4) = -d8*(d1+s)*(d1-t)
        shp(1,5) = -d8*(d1-s)*(d1+t)
        shp(1,6) =  d8*(d1-s)*(d1+t)
        shp(1,7) =  d8*(d1+s)*(d1+t)
        shp(1,8) = -d8*(d1+s)*(d1+t)
        shp(2,1) = -d8*(d1-r)*(d1-t)
        shp(2,2) = -d8*(d1+r)*(d1-t)
        shp(2,3) =  d8*(d1+r)*(d1-t)
        shp(2,4) =  d8*(d1-r)*(d1-t)
        shp(2,5) = -d8*(d1-r)*(d1+t)
        shp(2,6) = -d8*(d1+r)*(d1+t)
        shp(2,7) =  d8*(d1+r)*(d1+t)
        shp(2,8) =  d8*(d1-r)*(d1+t)
        shp(3,1) = -d8*(d1-r)*(d1-s)
        shp(3,2) = -d8*(d1+r)*(d1-s)
        shp(3,3) = -d8*(d1+r)*(d1+s)
        shp(3,4) = -d8*(d1-r)*(d1+s)
        shp(3,5) =  d8*(d1-r)*(d1-s)
        shp(3,6) =  d8*(d1+r)*(d1-s)
        shp(3,7) =  d8*(d1+r)*(d1+s)
        shp(3,8) =  d8*(d1-r)*(d1+s)
c
      else if(nel.eq.20.or.nel.eq.21) then
c....   shape functions for 20-node element 
        do l = 1,21
          r0 = r*ri(l)
          s0 = s*si(l)
          t0 = t*ti(l)
c         midside nodes top/bottom  r-dir
          if(l.eq.13.or.l.eq.15.or.l.eq.18.or.l.eq.20) then
            shp(4,l) = d4*(d1-rr)*(d1+s0)*(d1+t0)
            shp(1,l) =      -d2*r*(d1+s0)*(d1+t0)
            shp(2,l) = d4*(d1-rr)* si(l) *(d1+t0)
            shp(3,l) = d4*(d1-rr)*(d1+s0)* ti(l)
c         midside nodes top/bottom  s-dir
          else if(l.eq.14.or.l.eq.16.or.l.eq.19.or.l.eq.21) then
            shp(4,l) = d4*(d1+r0)*(d1-ss)*(d1+t0)
            shp(1,l) = d4* ri(l) *(d1-ss)*(d1+t0)
            shp(2,l) =-d2*(d1+r0)* s     *(d1+t0)
            shp(3,l) = d4*(d1+r0)*(d1-ss)* ti(l)
c         corner nodes mid plane
          else if(l.ge.9.and.l.le.12) then
            shp(4,l) = d4*(d1+r0)*(d1+s0)*(d1-tt)
            shp(1,l) = d4* ri(l) *(d1+s0)*(d1-tt)
            shp(2,l) = d4*(d1+r0)* si(l) *(d1-tt)
            shp(3,l) =-d2*(d1+r0)*(d1+s0)* t
c         corner nodes top/bottom
          else if(l.ge.1.and.l.le.8) then
            shp(4,l) = d8*(d1+r0)*(d1+s0)*(d1+t0)*(r0+s0+t0-dz)
            shp(1,l) = d8* ri(l) *(d1+s0)*(d1+t0)*(dz*r0+   s0+   t0-d1)
            shp(2,l) = d8*(d1+r0)* si(l) *(d1+t0)*(   r0+dz*s0+   t0-d1)
            shp(3,l) = d8*(d1+r0)*(d1+s0)* ti(l) *(   r0+   s0+dz*t0-d1)
          end if
        end do
        if(nel.eq.20) then        
c....     shape functions for 20-node element 
c....     modify node numbers nodes 18-21 to 17-20
          do l = 18,21 
            shp(4,l-1) = shp(4,l)
            shp(1,l-1) = shp(1,l)
            shp(2,l-1) = shp(2,l)
            shp(3,l-1) = shp(3,l)
          end do
          shp(4,21) = 0.d0
          shp(1,21) = 0.d0
          shp(2,21) = 0.d0
          shp(3,21) = 0.d0
        else
c....     shape functions for 20-node element with unused midpoint 17
c         node  17 (midpoint bottom) not used
        end if
c
      else if(nel.eq.27) then
c....   shape functions for 27-node element
        do l = 1,27
          r0 = r*ri(l)
          s0 = s*si(l)
          t0 = t*ti(l)
c         corner nodes top/bottom
          if(l.ge.1.and.l.le.8) then
            shp(4,l) = d8*(rr+r0)     *(ss+s0)          *(tt+t0)
            shp(1,l) = d8*(dz*r+ri(l))*(ss+s0)          *(tt+t0)
            shp(2,l) = d8*(rr+r0)     *(dz*s+si(l))     *(tt+t0)
            shp(3,l) = d8*(rr+r0)     *(ss+s0)          *(dz*t+ti(l))
c         corner nodes midside
          else if(l.ge.9.and.l.le.12) then
            shp(4,l) = d4*(rr+r0)     *(ss+s0)     *(d1-tt)
            shp(1,l) = d4*(dz*r+ri(l))*(ss+s0)     *(d1-tt)
            shp(2,l) = d4*(rr+r0)     *(dz*s+si(l))*(d1-tt)
            shp(3,l) = d4*(rr+r0)     *(ss+s0)     *(-dz*t)
c         midside nodes top/bottom  r-dir
          else if(l.eq.13.or.l.eq.15.or.l.eq.18.or.l.eq.20) then
            shp(4,l) = d4*(d1-rr)*(ss+s0)     *(tt+t0)
            shp(1,l) = d4*(-dz*r)*(ss+s0)     *(tt+t0)
            shp(2,l) = d4*(d1-rr)*(dz*s+si(l))*(tt+t0)
            shp(3,l) = d4*(d1-rr)*(ss+s0)     *(dz*t+ti(l))
c         midside nodes top/bottom  s-dir
          else if(l.eq.14.or.l.eq.16.or.l.eq.19.or.l.eq.21) then
            shp(4,l) = d4*(rr+r0)     *(d1-ss)*(tt+t0)
            shp(1,l) = d4*(dz*r+ri(l))*(d1-ss)*(tt+t0)
            shp(2,l) = d4*(rr+r0)     *(-dz*s)*(tt+t0)
            shp(3,l) = d4*(rr+r0)     *(d1-ss)*(dz*t+ti(l))
c         midside nodes mid plane  r-dir
          else if(l.eq.23.or.l.eq.25) then
            shp(4,l) = d2*(d1-rr)*(ss+s0)     *(d1-tt)
            shp(1,l) = d2*(-dz*r)*(ss+s0)     *(d1-tt)
            shp(2,l) = d2*(d1-rr)*(dz*s+si(l))*(d1-tt)
            shp(3,l) = d2*(d1-rr)*(ss+s0)     *(-dz*t)
c         midside nodes mid plane  s-dir
          else if(l.eq.24.or.l.eq.26) then
            shp(4,l) = d2*(rr+r0)     *(d1-ss)*(d1-tt)
            shp(1,l) = d2*(dz*r+ri(l))*(d1-ss)*(d1-tt)
            shp(2,l) = d2*(rr+r0)     *(-dz*s)*(d1-tt)
            shp(3,l) = d2*(rr+r0)     *(d1-ss)*(-dz*t)
c         central nodes top/bottom
          else if(l.eq.17.or.l.eq.22) then
            shp(4,l) = d2*(d1-rr)*(d1-ss)*(tt+t0)
            shp(1,l) = d2*(-dz*r)*(d1-ss)*(tt+t0)
            shp(2,l) = d2*(d1-rr)*(-dz*s)*(tt+t0)
            shp(3,l) = d2*(d1-rr)*(d1-ss)*(dz*t+ti(l))
          else if(l.eq.27) then
            shp(4,l) = (d1-rr)*(d1-ss)*(d1-tt)
            shp(1,l) = (-dz*r)*(d1-ss)*(d1-tt)
            shp(2,l) = (d1-rr)*(-dz*s)*(d1-tt)
            shp(3,l) = (d1-rr)*(d1-ss)*(-dz*t)
          end if
        end do
c
      else if(nel.eq.64) then
c....   shape functions for 64-node element
        call shp3dc(sg, shp)
c
      end if
      return
      end
c
      subroutine stre21(ix,dt,st,plout,sig,sig0,shp,nel,numnp,dvp)
c-----------------------------------------------------------------------
c.... plot stresses for nonlinear 3D-element
c....  1=S_xx      2=S_yy   3=S_zz 4=T_xy 5=T_xz 6=T_yz
c....  7=S_11      8=S_22   9=S_33
c.... 16=Plout(1) 25=Plout(10)
c       
c.....(c) w.wagner aug 86
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ix(*),dt(numnp),st(numnp,*),shp(4,nel),sig(*),sig0(3),
     +          plout(10) 
      do 10 i = 1,nel
        xsji = dvp*shp(4,i)
        ii = abs(ix(i))
        if(ii.eq.0) goto 10
          dt(ii) = dt(ii) + xsji
          st(ii,1) = st(ii,1) + sig (1)*xsji
          st(ii,2) = st(ii,2) + sig (2)*xsji
          st(ii,3) = st(ii,3) + sig (3)*xsji
          st(ii,4) = st(ii,4) + sig (4)*xsji
          st(ii,5) = st(ii,5) + sig (5)*xsji
          st(ii,6) = st(ii,6) + sig (6)*xsji
          st(ii,7) = st(ii,7) + sig0(1)*xsji
          st(ii,8) = st(ii,8) + sig0(2)*xsji
          st(ii,9) = st(ii,9) + sig0(3)*xsji
          do k = 1,10
            st(ii,15+k) = st(ii,15+k) + plout(k)*xsji
          end do
10    continue
      return
      end
c
      subroutine gausnod21(sg,lint)
c-----------------------------------------------------------------------
c.....Gauss points = corner nodes
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension sg(4,*)
c.... xsi,eta
c.... node 1,5 
      sg(1,1) = -1.d0 
      sg(2,1) = -1.d0 
      sg(1,5) = -1.d0 
      sg(2,5) = -1.d0 

c.... node 2,6 
      sg(1,2) =  1.d0 
      sg(2,2) = -1.d0 
      sg(1,6) =  1.d0 
      sg(2,6) = -1.d0 

c.... node 3,7 
      sg(1,3) =  1.d0 
      sg(2,3) =  1.d0 
      sg(1,7) =  1.d0 
      sg(2,7) =  1.d0 

c.... node 4,8 
      sg(1,4) = -1.d0 
      sg(2,4) =  1.d0 
      sg(1,8) = -1.d0 
      sg(2,8) =  1.d0 

c.... zeta node 1-4,5-8 
      do im=1,4
        ip=im+4
        sg(3,im)= -1.d0   
        sg(3,ip)=  1.d0   
      end do

c.... weight
      do i=1,8
        sg(4,i)=  1.d0   
      end do

c.... number of GP
      lint = 8

      return
      end
c
      subroutine stress3d21(sig,sig0)
c----------------------------------------------------------------------+
c.... solv  [a -   lambda 1] x = 0    3*3                              |
c----------------------------------------------------------------------+
c
      implicit double precision (a-h,o-z)
      dimension sig(3,3),b(3,3),sig0(3),z(3,3),d(3,3),fv1(3),fv2(3)
      data ev /0.d0/
      call pzero(sig0,3)
      call pzero(b,9)
      b(1,1) = 1.d0
      b(2,2) = 1.d0
      b(3,3) = 1.d0
      call pzero(z,9)
      call pzero(d,9)
      call pzero(fv1,3)
      call pzero(fv2,3)
c.....solve EV-Problem from EISPACK
      matz = 1
      ierr = 0
      nm = 3
      call rsg(nm,3,sig,b,sig0,matz,z,fv1,fv2,ierr)
      if(matz.ne.0) then
c....   scale to length 1
        do k = 1,3
          ev = 0.d0
          do i = 1,3
            ev = ev + z(i,k)*z(i,k)
          end do
          ev = dsqrt(ev)
          do i = 1,3
            z(i,k) = z(i,k)/ev
          end do
        end do
      end if
c.... sort stresses    2-3
10    if(sig0(2).lt.sig0(3)) then
         shelp   = sig0(3) 
         sig0(3) = sig0(2)
         sig0(2) = shelp
      end if
c.... sort stresses    1-2
      if(sig0(1).lt.sig0(2)) then
         shelp   = sig0(2) 
         sig0(2) = sig0(1)
         sig0(1) = shelp
      end if
c.... sort stresses    2-3
      if(sig0(2).lt.sig0(3)) goto 10
c 
      return
      end
c
      subroutine qload21(d,q,ql,numel,n,mqloa,propq,prop,isw)
c----------------------------------------------------------
c.... set loads from macro qloa/mate
c     MATE
c     d(10) = qx     = volume load in x-direction
c     d(11) = qy     = volume load in y-direction
c     d(12) = qz     = volume load in z-direction
c     d(13) = theta  = const. temperature (load)
c     d(14) = at     = alpha_t
c     QLOA
c     q(1) = qx     = volume load in x-direction
c     q(2) = qy     = volume load in y-direction
c     q(3) = qz     = volume load in z-direction
c     q(4) = theta  = const. temperature (load)
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),d(*),ql(4)
      if(isw.eq.22) then ! QLOA 
        if(mqloa.ne.1) then
          do i=1,4 
            ql(i) = q(n,i)*propq 
          end do
        end if  
      else if(isw.eq.4.or.isw.eq.8) then ! STRE
        ql(1) = d(10)*prop
        ql(2) = d(11)*prop
        ql(3) = d(12)*prop
        ql(4) = d(13)*prop
        if(mqloa.ne.1) then
          do i=1,4 
            ql(i) = ql(i)+q(n,i)*propq 
          end do
        end if  
      else  ! TANG etc
        ql(1) = d(10)*prop
        ql(2) = d(11)*prop
        ql(3) = d(12)*prop
        ql(4) = d(13)*prop
      end if
      return
      end
c
c      subroutine gaus21(ngaus,pg,wg)
cc-----------------------------------------------------------------------
cc.....Gauss points  ngaus = 1,3
cc-----------------------------------------------------------------------
c      implicit double precision(a-h,o-z)
c      dimension pg(3),wg(3)
c      call pzero(pg,3)
c      call pzero(wg,3)
c      goto (1,2,3) ngaus
cc.....1 pt  1x1x1
c1     pg(1) = 0.0d0
c      wg(1) = 2.0d0
c      return
cc.....2 pt  2x2x2
c2     pg(1) = -dsqrt(3.0d0)/3.0d0
c      pg(2) = -pg(1)
c      wg(1) =  1.0d0
c      wg(2) =  wg(1)
c      return
cc.....3 pt  3x3x3
c3     pg(1) = -dsqrt(0.6d0)
c      pg(2) =  0.0d0
c      pg(3) = -pg(1)
c      wg(1) =  5.0d0/9.0d0
c      wg(2) =  8.0d0/9.0d0
c      wg(3) =  wg(1)
c      return
c      end
c
