      subroutine elmt45(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-----------------------------------------------------------------------
c
c.....geometrically nonlinear solid shell element 8 nodes
c     enhanced assumed strain method
c
c-----------------------------------------------------------------------
c.... Input element parameters
c     1. Card  Element switch parameters
c       d( 1)   imat: material type, see matelib3d
c       d( 2)   nlay: number of layers (def.=1)
c       d( 3)   lin : linear = 0  nonlinear = 1
c       d( 4)   ieas: Enhanced Assumed Strains (EAS) = 0,
c                    =  3 bending patch test   (Rah et al. 2011)
c                    =  5 membrane patch test  (sk,fg,ww 1999)
c                    =  8 & dist. mesh         (sk,fg,ww 1999)
c                    = 11 & vol. locking       (sk,fg,ww 1999)
c                         & bending patch test
c                    =  7 bending patch test   (Vu-Quoc 2003)
c                    = 30 enhanced brick       (sk,ww 1997)
c       d( 5)   ibd : Bathe-Dvorkin 1/0
c       d( 6)   ibs : Betsch-Stein 1/0
c       d( 7)   ieas: 7 EAS (Rah et al 2011) instead of (Vu-Quoc 2003)
c                    = 1 vol. locking ok for reg mesh
c                    = 2 dist. mesh bending
c       d( 8)   shell thickness=thickness of 1 Element
c               used for thickness integration
c       d( 9)   <no. of Gauss Points/layer> def.=2
c       d(10)   shear corr.factor:(<0=|value|*FE-SCF, >0=value, def.=1)
c
c     2. Card  Volumetric loading
c       d(11)  qx   : volumetric loading in x-direction
c       d(12)  qy   : volumetric loading in y-direction
c       d(13)  qz   : volumetric loading in z-direction
c       d(14)  theta: const. temperature (load)
c       d(15)  at   : alpha_t
c       d(16)  rho  : gamma/g
c
c       further storage in d
c       d(17)
c       d(18)
c       d(19)  nh   : length of history array at Gauss-Point
c       d(20)  ivor : length of d-array without layer data
c
c     3. Card   material data of material imat
c       d(21)-d(..)   see SR matelib3d
c
c.... 4. card and following: layer data
c     Input data for each layer i (on separate input line)
c     beginning at z=-h/2
c        phi   = angle between global and local coordinate system
c        hi    = thickness of layer
c
c-----------------------------------------------------------------------
c
c     isw= 3  tangent matrix+residuum
c     isw= 4  print stresses
c     isw= 5  lumped mass matrix
c     isw= 6  residuum
c     isw= 7  non-conservative loads
c     isw= 8  plot stresses
c     isw=11  extended system
c     isw=15  update history arrays, not tested!!
c
c-----------------------------------------------------------------------
c     (c)  J. Schuett, S. Klinkel  10/97
c          layered formulation for mat.law             FG/WW 06/2012
c          isw=15                                      FG/WW 02/2013
c          shear correction factor                     FG/WW 03/2014
c-----------------------------------------------------------------------
      USE bdata
      USE cdat1
      USE cdata
      USE eldata
      USE hdata
      USE iofile
      USE pdata6
      USE plslay
      USE prisdat
      USE prlod
      USE strnam
      implicit double precision(a-h,o-z)
      parameter(neas=30, nges=neas+24)
      dimension ix(*),xl(ndm,*),tl(*),ul(ndf,*),d(*),ipordl(17),
     +          sg(16),tg(16),wg(16),pgz(5),wgz(5),wg2(4,9),
     +          gp(3),shp(4,8),s(nst,nst),cmat(6,6),dd(30),
     +          bmat(6,3),btd(3,6),gmat(3,3),p(nst),sig(7),du(24),
     +          glo(3,3),glu(3,3),gu(3,3),eu(3,3),e(6),
     +          ToM(6,6),Te(6,6),xjac(3,3),g0o(3,3),g0u(3,3),g00(3,3),
     +          sl(24,24),pl(24),h1(*),h2(*),h3(*),t(3,3),Tm(6,6),
     +          shpk(3,8,4),pkl(3,3,4),pk(3,3,4),eup(3,3,4),xu(3,8),
     +          shpa(3,4),pa(3,3,4),pal(3,3,4),eua(3,3,4),shpka(3,8,3,4)
     +         ,plout(10)
      dimension Sc(nges,nges),Pc(nges),xmtc(nges,6),xm(6,neas),
     +          alpha(neas),eas(6),dvl(nges)
c
      common /batdvo/ shpk,pk,pkl,eup
!$OMP THREADPRIVATE (/batdvo/)
      common /betste/ pa,pal,shpa,shpka,eua
!$OMP THREADPRIVATE (/betste/)
c
      data ipordl /1,2,3,4,1,5,6,2,6,7,3,7,8,4,8,5,1/
      ielno = 45
c
c.... go to correct array processor
      go to(1,2,3,3,5,3,7,3,2,2,3,2,2,2,3,2,2,2,2,2,2,2), isw
c
c.... input material properties
c.... card 1
1     call dinput(d(1),10)
      if(d( 2).eq.0.d0) d( 2)=1.d0  ! default nlay=1
      if(d( 9).eq.0.d0) d( 9)=2.d0  ! Def.no.GP/lay=2
      if(d(10).eq.0.d0) d(10)=1.d0  ! Def. kappa
c.... card 2
      call dinput(d(11),8)
      write(iow,1002) (d(i),i=1,18)
1002  format(5x,'Element data for EAS/BaDvo/BeSt solid shell',/,
     1  5x,'material type..............',f12.5,/,
     2  5x,'number of layers...........',f12.5,/,
     3  5x,'lin=0,nonlin=1.............',f12.4,/,
     4  5x,'eas parameter:.............',f12.4,/,
     5  5x,'Bathe-Dvorkin, 0/1.........',f12.4,/,
     6  5x,'Betsch-Stein, 0/1..........',f12.4,/,
     7  5x,'1=vol.lock., 2=bending.....',f12.5,/,
     8  5x,'Shell thickness............',f12.5,/,
     9  5x,'no. Gauss Points/layer.....',f12.5,/,
     +  5x,'SCF(<0:|VAL|*FE-SCF,>0:VAL)',f12.5,/,
     1  5x,'qx.........................',f12.5,/,
     2  5x,'qy.........................',f12.5,/,
     3  5x,'qz.........................',f12.5,/,
     4  5x,'delta T ...................',f12.5,/,
     5  5x,'alpha T ...................',f12.5,/,
     6  5x,'mass density.rho=gamma/g...',f12.5,/,
     7  5x,'...........................',f12.5,/,
     8  5x,'...........................',f12.5)
c
c.... card 3
      imat = d(1)
      call matelib3d(h1,h2,nh,d(21),mdx,E,Sig,Cmat,6,1,plout,
     +               gp,tgp,dvp,1.d0,1.d0,1.d0,n,l,1,1,imat,isw)
      d(19) = nh
c...  cheque length of d and dd-array
      ivor = 20 + mdx
      d(20)= ivor
      nlay = d(2)
      nmax = ivor + 3 * nlay
      if(nmax.gt.ndd) then
         write(*,1005) ndd,nmax
1005     format(1x,'d array is set to  ',i4,' values',/,
     1          1x,'d array needs      ',i4,' values')
         stop
      end if
      if(mdx.gt.30) then
         write(*,1011) mdx
1011     format(1x,'dd array is set to 30 values',/,
     1          1x,'dd array need  ',i4,' values')
         stop
      end if
c.... card 4 and following: phi,h for each layer
                   write(iow,1008)
      if(ior.lt.0) write(*  ,1008)
1008  format(5x,'material data for layers',/,
     +5x,'layer  angle phi-x1  coordinate zi  thickness hi')
      hlay2=0.d0
      zi   = -d(8)*0.5d0
      do ilay =1,nlay
        ia = (ilay - 1) * 3 + ivor
        call dinput(dd,3)
        d(ia+1) = dd(1)          ! phi
        zi      = zi+0.5d0*(hlay2+dd(2))
        hlay2   = dd(2)
        d(ia+2) = zi           ! zi
        d(ia+3) = dd(2)        ! hi
                     write(iow,1009) ilay,(d(i),i=ia+1,ia+3)
        if(ior.lt.0) write(*  ,1009) ilay,(d(i),i=ia+1,ia+3)
      end do
1009  format(5x,i5,3(2x,g12.5))
      write(iow,1010)
1010  format(//)
c
c.... description of stresses
c     plot>stre,i       = stress i, if layered: sum of stress (check sense!?)
c     plot>stre,i,,-k.m = stress i for layer k and gp m
      strsus( 1) = ' S_11          '
      strsus( 2) = ' S_22          '
      strsus( 3) = ' S_33          '
      strsus( 4) = ' S_12          '
      strsus( 5) = ' S_13          '
      strsus( 6) = ' S_23          '
      strsus(7:15) = '               '
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
      nptyp = 4
c
c.... define node numbering for plot mesh routine, see pltord
      inord(ielno) = 17
      do ii = 1,17
         ipord(ii,ielno) = ipordl(ii)
      end do
c
c.....dimension h-array
      ngb   = 2
      nlay  = d(2)
      ieas  = d(4)
      ngz   = d(9)
      nh    = d(19)
      iges  = ieas + 24
      nhela = ieas*(iges+2)
      nhpla = nh*ngb*ngb*nlay*ngz
      nh3   = nhela               ! h3--array for eas
      nh1   = nhpla               ! h1/h2--array for plasti
c
2     return
3     ngb  = 2
      imat = d(1)
      nlay = d(2)
      lin  = d(3)
      ieas = d(4)
      ibd  = d(5)
      ibs  = d(6)
      hs   = d(8)
      ngz  = d(9)
      nh   = d(19)
      ivor = d(20)
      iges  = ieas + 24

      call pgauss (ngb,lint,sg,tg,wg)
      call gaus1D (ngz,pgz,wgz)

      do i = 1,ivor-20
       dd(i) = d(20+i)
      end do

c.....enhanced element values
      call pzero(g0o,3*3)
      if (ieas.gt.0) then
         call pzero(Sc,nges*nges)
         call pzero(Pc,nges)
         call pzero(alpha,neas)
         call pzero(dvl,nges)
         call pzero(g00,3*3)
c.....   compute transformation matrix To
         call shp45(0.d0,0.d0,0.d0,xl,ul,g0u,g0o,gu,eu,shp,shpa,
     +              ToM,Te,xjac,dj0,n,iow,xu,g00,ndf)
c.....   update internal dofs
         nsc = 1
         npc = nsc + ieas * (ieas + 24)
         nalp= npc + ieas
         do i=1,nen
            du((i-1)*3+1) = ul(1,2*nen+i)
            du((i-1)*3+2) = ul(2,2*nen+i)
            du((i-1)*3+3) = ul(3,2*nen+i)
         end do
c
         call decond45(h3(nsc),h3(npc),h3(nalp),du,dvl,alpha,
     +               ieas,24,nges)
      end if
c
c.....shape functions for shear strains
      if (ibd.gt.0) call batdvo45(shpk,xl,ul,pkl,pk,eup,ndf)
c
c.....shape functions for thickness strains
      if (ibs.gt.0) call betste45(xl,ul,pa,pal,eua,shpka,ndf)
      call pzero(sl,24*24)
      call pzero(pl,24)
c
c.... shear correction factor
      call shearfac45(d,xl,cappa,wcappa)

c.....loop over gauss points
      nn = 1
      do 30 igaus = 1,lint
       xsi = sg(igaus)
       eta = tg(igaus)

       do ilay = 1,nlay
        ia  = (ilay - 1) * 3 + ivor
        phi = d(ia+1)
        zi  = d(ia+2)
        hi  = d(ia+3)
        hi5 = 0.5d0*hi
        if(imat.eq. 2) dd(10) = phi
        if(imat.eq. 3) dd( 6) = phi
        if(imat.eq. 7) dd(28) = phi
        if(imat.eq.15) dd( 6) = phi

        do igz = 1,ngz
         wz = wgz(igz)*hi/hs
         zetai = pgz(igz)
         zs    = zi + zetai*hi5
         zeta  = 2.d0*zs/hs
c
c.....   shape functions and derivatives
         call shp45(xsi,eta,zeta,xl,ul,glu,glo,gu,eu,shp,shpa,
     +             ToM,Te,xjac,detj,n,iow,xu,g0o,ndf)
         gp  = 0
         tgp = 0
         do inode = 1,nel
           do i = 1,3
             gp(i) = gp(i) + xl(i,inode)* shp(4,inode)
           end do
           tgp = tgp + tl(inode)* shp(4,inode)
         end do
         dv = detj*wg(igaus)*wz
         da = detj/hi
c
c.....   interpolation matrix xM and enhanced strains eas
         if (ieas.gt.0) then
            call mmat45(xm,detj,xsi,eta,zeta,d,ToM,dj0)
            do i=1,6
               eas(i)=0.d0
               do j=1,ieas
                  eas(i) = eas(i) + xm(i,j) * alpha(j)
               end do
            end do
         end if
c
c.....   kinematic Green-Langrangean strains e
         call kine45(d,gu,glu,eu,xsi,eta,eas,e)
c
c....    stresses and tangent matrix at t_n+1
         call trafo45_local(e,glu,Tm,sig,cmat,1) ! trafo to cart. coor.
         call trafo45_shear(cappa,wcappa,e,cmat,sig,1) ! trafo E
         call matelib3d(h1(nn),h2(nn),nh,dd,mdx,e,sig,Cmat,6,1,plout,
     +               gp,tgp,dv,1.d0,1.d0,1.d0,n,igaus,ilay,igz,imat,isw)
         call trafo45_shear(cappa,wcappa,e,cmat,sig,2) ! trafo C,S
         if(isw.eq. 4) goto  4 ! print stre
         if(isw.eq. 8) goto  8 ! plot  stre
         if(isw.eq.15) goto 29 ! updh
         call trafo45_local(e,glu,Tm,sig,cmat,2) ! trafo to conv. coor.
c
c.....   nonlinear stiffness matrix and residual
        if (ieas.gt.0) then
c.....     pc(25-nges) = -h = - int( M^T * Sig * dv)
           do 104 i=1,ieas
              do 104 j=1,6
                 xmtc(i,j)=0.d0
                 pc(24+i) = pc(24+i) - xm(j,i) * sig(j) * dv
                 do 104 k=1,6
                 xmtc(i,j) = xmtc(i,j) + xm(k,i) * cmat(k,j)
104        continue
c.....     sc(25-nges,25-nges)= H = int( M^T * C * M)
           do 105 i=1,ieas
              do 105 j=i,ieas
                 do 105 k=1,6
                 Sc(24+i,24+j)=Sc(24+i,24+j)+xmtc(i,k)*xm(k,j)*dv
105        continue
        end if
c
c.....  nonlinear stiffness matrix and residual
        i1=0
        do 31 ii=1,nel
c.....    external load vector
          qx = d(11) * prop
          qy = d(12) * prop
          qz = d(13) * prop
          iq=i1+1
          pl(iq)   = pl(iq)   + shp(4,ii)*qx*dv
          pl(iq+1) = pl(iq+1) + shp(4,ii)*qy*dv
          pl(iq+2) = pl(iq+2) + shp(4,ii)*qz*dv
c....     B- matrix
          call bmat45(gu,glu,shp,ii,bmat,xsi,eta,d)
c....     residual G = P - Bt*S and matrix Bt*D
          do 32 i = 1,3
            do  33 k = 1,6
              btd(i,k) =0.0d0
              pl(i1+i) = pl(i1+i) - bmat(k,i)*sig(k)*dv
              do  34 j = 1,6
34            btd(i,k) = btd(i,k)+bmat(j,i)*cmat(j,k)
33          continue
32        continue
          if (ieas.gt.0) then
c.....       sc(25-nges,24) = LT = int( M^T C * B )^T
             do 107 i=1,3
                do 107 j=1,ieas
                   do 107 k=1,6
                  sc(i1+i,24+j)=sc(i1+i,24+j)+xmtc(j,k)*bmat(k,i)*dv
107          continue
          end if
c....     tangent stiffness matrix
          j1 = i1
          do 35 jj = ii,nel
            call bmat45(gu,glu,shp,jj,bmat,xsi,eta,d)
            do 36  i = 1,3
              do 37  j = 1,3
                do 38  k = 1,6
38              sl(i1+i,j1+j) = sl(i1+i,j1+j) + btd(i,k)*bmat(k,j)*dv
37            continue
36          continue
c.....      initial stress matrix
            if(lin.ne.0) then
            call gmat45(shp,ii,jj,sig,gmat,xsi,eta,d)
            sl(i1+1,j1+1) = sl(i1+1,j1+1) + gmat(1,1)*dv
            sl(i1+2,j1+2) = sl(i1+2,j1+2) + gmat(2,2)*dv
            sl(i1+3,j1+3) = sl(i1+3,j1+3) + gmat(3,3)*dv
            end if
            j1 = j1 + 3
35        continue
          i1 = i1 + 3
31      continue
        goto 29
c
c....   Output Stresses
4       continue
        if ((ilay.eq.klay.and.igz.eq.mlay).or. klay.eq.0) then

        do idm = 1,3
           gp(idm) = 0.0d0
           do inode = 1,nel
              gp(idm) = gp(idm) + xl(idm,inode)* shp(4,inode)
           end do
        end do
        mct = mct - 1
        if (mct.le.0) then
                        write(iow,4000) o,head
           if(ior.lt.0) write(  *,4000) o,head
           mct = 50
        end if
                     write(iow,4001) n,ma,(gp(i),i=1,3),(sig(i),i=1,6)
        if(ior.lt.0) write(  *,4001) n,ma,(gp(i),i=1,3),(sig(i),i=1,6)

4000    format(a1,20a4,/,2x,'Element Stresses (2.P-K)',/,
     +  1x,'EL',1x,'M',
     +  1x,'    1-COR   ',
     +  1x,'    2-COR   ',
     +  1x,'    3-COR   ',
     +  1x,'    S_11    ',
     +  1x,'    S_22    ',
     +  1x,'    S_33    ',
     +  1x,'    S_12    ',
     +  1x,'    S_13    ',
     +  1x,'    S_23    ')
4001    format(1x,i3,1x,i2,9(1x,g12.5))

      end if

c        write(iow,4001) n,ma,(gp(i),i=1,3)
c        write(iow,4002) (sig(i),i=1,6)
c        write(*,4001) n,ma,(gp(i),i=1,3)
c        write(*,4002) (sig(i),i=1,6)
c4000    format(a1,20a4,/,2x,'Element Stresses (2.P-K)',/,
c     +  1x,'EL',1x,'M',1x,'1-COR',1x,'2-COR',1x,'3-COR',/,
c     +  1x,'S_11',1x,'S_22',1x,'S_33',1x,'S_12',1x,'S_13',1x,'S_23')
c4001    format(1x,i3,1x,i2,3g12.5)
c4002    format(6g12.5)
        goto 29
c
c....    plot stresses for klay.mlay or sum of stresses for klay=0
8        istv = -25
         if ((ilay.eq.klay.and.igz.eq.mlay).or. klay.eq.0) then
          if(iplma(ma).eq.0)       return ! only if MATN
          call strepl45(ix,strea,strea(1+numnp),plout,sig,shp,nel,
     +    numnp,dv)
         end if
c
29       nn = nn + nh

        end do  ! igz
       end do   ! ilay
30    continue ! igaus
      if(isw.eq.4.or.isw.eq.8.or.isw.eq.15) return
c
c.....lower part of stiffness matrix
      do 39 i = 1,24
        do 39 j = 1,i
   39    sl(i,j) =sl(j,i)
c.....lower part of stiffness matrix k enhanced
      if (ieas.gt.0) then
        do 120 i=24+1,iges
           do 120 j=24+1,i
              Sc(i,j) = Sc(j,i)
120     continue
        do 121 i=24+1,iges
           do 121 j=1,24
              Sc(i,j) = Sc(j,i)
121     continue
c.....  copy stiffness matrix sl into Sc and pl into Pc
        do 122 i=1,24
           pc(i) = pl(i)
           do 122 j=1,24
              Sc(i,j) = sl(i,j)
122     continue
c.....  condensation S = xK - zL^T H^-1 zL
c.....               P = fext - fint + L^T H^-1 h
        call conden(Sc,Pc,24,ieas,nges,.false.)
c.....  store matrices for iteration
        call store45(sc,pc,alpha,h3(nsc),h3(npc),h3(nalp),24,ieas,nges)
c.....  copy stiffness matrix Sc into Sl and Pc into Pl
        do 123 i=1,24
           Pl(i) = Pc(i)
           do 123 j=1,24
              sl(i,j) = Sc(i,j)
123     continue
      end if
c
c.....copy stiffness matrix sl into S and pl into P
      i1=0
      do ii=1,nel
         i = (ii-1)*3
         j1=0
         p(i1+1:i1+3) = pl(i+1:i+3)
         do jj=1,nel
            j = (jj-1)*3
            s(i1+1:i1+3,j1+1:j1+3) = sl(i+1:i+3,j+1:j+3)
            j1=j1+ndf
         end do
         i1=i1+ndf
      end do
c
      return
c.... lumped mass matrix
5     rho = d(16)
      call int3d(2,lint,wg2)
c.....loop over gauss points
      do l = 1,lint
c.....  shape functions, derivatives and jacobian
         xsi = wg2(1,l)
         eta = wg2(2,l)
         zeta= wg2(3,l)
         call shp45(xsi,eta,zeta,xl,ul,glu,glo,gu,eu,shp,shpa,
     +             ToM,Te,xjac,detj,n,iow,xu,g0o,ndf)
         dvp = detj*wg2(4,l)
c....   loop over nodes
        i1=0
        do ii=1,8
          add = rho*shp(4,ii)*dvp
          p(i1+1) = p(i1+1) + add
          p(i1+2) = p(i1+2) + add
          p(i1+3) = p(i1+3) + add
          i1 = i1 + ndf
        end do
      end do
      return
7     continue
c.....e.g. non-conservative loads,
c.....internal pressure 1 ele through thickness
      call int3d(2,lint,wg2)
      call pzero(g0o,3*3)
      do igaus =1,lint
         xsi = wg2(1,igaus)
         eta = wg2(2,igaus)
         zeta= wg2(3,igaus)
c.....   shape functions and derivatives
         call shp45(xsi,eta,zeta,xl,ul,glu,glo,gu,eu,shp,shpa,
     +              ToM,Te,xjac,detj,n,iow,xu,g0o,ndf)
         da = wg2(4,igaus)/2.d0    !/2 weil 2 Knoten über die Dicke
         i1=0
         do ii=1,nel
c.....     external load vector
           iq=i1+1
           pz = d(1)
           p1=(gu(2,1)*gu(3,2)-gu(2,2)*gu(3,1))*pz
           p2=(gu(3,1)*gu(1,2)-gu(3,2)*gu(1,1))*pz
           p3=(gu(1,1)*gu(2,2)-gu(1,2)*gu(2,1))*pz
           p(iq)  =p(iq)   + p1*shp(4,ii) * da
           p(iq+1)=p(iq+1) + p2*shp(4,ii) * da
           p(iq+2)=p(iq+2) + p3*shp(4,ii) * da
           j1 = 0
           do jj = 1,nel
             s1=(-gu(3,1)*shp(2,jj)+gu(3,2)*shp(1,jj))*shp(4,ii)*pz
             s2=(+gu(2,1)*shp(2,jj)-gu(2,2)*shp(1,jj))*shp(4,ii)*pz
             s3=(-gu(1,1)*shp(2,jj)+gu(1,2)*shp(1,jj))*shp(4,ii)*pz
             s(i1+1,j1+2) = s(i1+1,j1+2) - s1*da
             s(i1+1,j1+3) = s(i1+1,j1+3) - s2*da
             s(i1+2,j1+3) = s(i1+2,j1+3) - s3*da
             s(i1+2,j1+1) = s(i1+2,j1+1) + s1*da
             s(i1+3,j1+1) = s(i1+3,j1+1) + s2*da
             s(i1+3,j1+2) = s(i1+3,j1+2) + s3*da
             j1 = j1 + ndf
           end do
           i1 = i1 + ndf
         end do
      end do
c
      return
      end
c
c----------------------------------------------------------------------
c                            subroutines
c----------------------------------------------------------------------
c
      subroutine shp45 (xsi,eta,zeta,xl,ul,glu,glo,gu,eu,shp,shpa,
     +                  ToM,Te,xjac,detj,n,iow,xu,g0o,ndf)
c-----------------------------------------------------------------------
c.....compute shape functions and their derivatives for linear,quadratic
c.....lagrangian and serendipity isoparametric  3-d elements
c.....global coordinate system x,y,z
c.....local coordinate system xsi,eta,zeta
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension shp(4,8),t(3,3),e(3,3),cro(3),ToM(6,6),to(3,3),
     +          xl(3,8),ul(ndf,8),xu(3,8),xjac(3,3),Te(6,6),
     +          glu(3,3),glo(3,3),gu(3,3),gko(3,3),eu(3,3),g0o(3,3),
     +          shpa(3,4),xj(3,3)
      od=1.0d0
      ed=8.0d0
c
c.....shape functions for 8-node element
      shp(4,1) = (od-xsi)*(od-eta)*(od-zeta)/ed
      shp(4,2) = (od+xsi)*(od-eta)*(od-zeta)/ed
      shp(4,3) = (od+xsi)*(od+eta)*(od-zeta)/ed
      shp(4,4) = (od-xsi)*(od+eta)*(od-zeta)/ed
      shp(4,5) = (od-xsi)*(od-eta)*(od+zeta)/ed
      shp(4,6) = (od+xsi)*(od-eta)*(od+zeta)/ed
      shp(4,7) = (od+xsi)*(od+eta)*(od+zeta)/ed
      shp(4,8) = (od-xsi)*(od+eta)*(od+zeta)/ed
c
c.....derivatives
c.....shp(a,b) is equivalent to SHPb,a with b=1 to 8  and a=1 to 3
      shp(1,1) = (-od)*(od-eta)*(od-zeta)/ed
      shp(1,2) = ( od)*(od-eta)*(od-zeta)/ed
      shp(1,3) = ( od)*(od+eta)*(od-zeta)/ed
      shp(1,4) = (-od)*(od+eta)*(od-zeta)/ed
      shp(1,5) = (-od)*(od-eta)*(od+zeta)/ed
      shp(1,6) = ( od)*(od-eta)*(od+zeta)/ed
      shp(1,7) = ( od)*(od+eta)*(od+zeta)/ed
      shp(1,8) = (-od)*(od+eta)*(od+zeta)/ed

      shp(2,1) = (od-xsi)*(-od)*(od-zeta)/ed
      shp(2,2) = (od+xsi)*(-od)*(od-zeta)/ed
      shp(2,3) = (od+xsi)*( od)*(od-zeta)/ed
      shp(2,4) = (od-xsi)*( od)*(od-zeta)/ed
      shp(2,5) = (od-xsi)*(-od)*(od+zeta)/ed
      shp(2,6) = (od+xsi)*(-od)*(od+zeta)/ed
      shp(2,7) = (od+xsi)*( od)*(od+zeta)/ed
      shp(2,8) = (od-xsi)*( od)*(od+zeta)/ed

      shp(3,1) = (od-xsi)*(od-eta)*(-od)/ed
      shp(3,2) = (od+xsi)*(od-eta)*(-od)/ed
      shp(3,3) = (od+xsi)*(od+eta)*(-od)/ed
      shp(3,4) = (od-xsi)*(od+eta)*(-od)/ed
      shp(3,5) = (od-xsi)*(od-eta)*(+od)/ed
      shp(3,6) = (od+xsi)*(od-eta)*(+od)/ed
      shp(3,7) = (od+xsi)*(od+eta)*(+od)/ed
      shp(3,8) = (od-xsi)*(od+eta)*(+od)/ed

c.....shape functions for Betsch-Stein interpolation
      f=1.d0/4.d0

      shpa(3,1) = f * (1-xsi)*(1-eta)
      shpa(3,2) = f * (1+xsi)*(1-eta)
      shpa(3,3) = f * (1+xsi)*(1+eta)
      shpa(3,4) = f * (1-xsi)*(1+eta)

c.....convect. coordinates ref.config gl(j,k) ; j=1,2,3=x,y,z ;k=1,2,3
c     covariant basis
      do 10 k=1,3
        do 10 j=1,3
          glu(j,k)=0.0d0
          do 10 ii=1,8
            glu(j,k)=glu(j,k)+shp(k,ii)*xl(j,ii)
10    continue
c
c.....coefficient matrix
      do 15 i=1,3
         do 15 j=1,3
            gko(i,j) = 0.d0
            do 15 k=1,3
            gko(i,j) = gko(i,j) + glu(k,i) * glu(k,j)
15    continue
      call invert(gko,3,3)
c
c.....contravariant basis
      do i=1,3
         glo(i,1)=gko(1,1)*glu(i,1)+gko(1,2)*glu(i,2)+gko(1,3)*glu(i,3)
         glo(i,2)=gko(2,1)*glu(i,1)+gko(2,2)*glu(i,2)+gko(2,3)*glu(i,3)
         glo(i,3)=gko(3,1)*glu(i,1)+gko(3,2)*glu(i,2)+gko(3,3)*glu(i,3)
      end do
c
c.....transform into cart. coordintes t(i)
      betr=dsqrt(glu(1,1)*glu(1,1)+glu(2,1)*glu(2,1)+glu(3,1)*glu(3,1))
      t(1,1) = glu(1,1)/betr
      t(2,1) = glu(2,1)/betr
      t(3,1) = glu(3,1)/betr

      cro(1) =  glu(2,1)*glu(3,2)-glu(3,1)*glu(2,2)
      cro(2) = -glu(1,1)*glu(3,2)+glu(3,1)*glu(1,2)
      cro(3) =  glu(1,1)*glu(2,2)-glu(2,1)*glu(1,2)
      betr = dsqrt(cro(1)*cro(1)+cro(2)*cro(2)+cro(3)*cro(3))
      t(1,3) = cro(1)/betr
      t(2,3) = cro(2)/betr
      t(3,3) = cro(3)/betr

      t(1,2) =  t(2,3)*t(3,1)-t(3,3)*t(2,1)
      t(2,2) = -t(1,3)*t(3,1)+t(3,3)*t(1,1)
      t(3,2) =  t(1,3)*t(2,1)-t(2,3)*t(1,1)

c.....transformation matrix (convect. into midpoint convect. coos)
      call pzero(xj,3*3)
      do i=1,3
         do k=1,3
            do j=1,3
            xj(i,k) = xj(i,k) + glu(j,i) * g0o(j,k)
            end do
         end do
      end do
c     only matrix To
      Tom(1,1) =       xj(1,1)*xj(1,1)
      Tom(1,2) =       xj(1,2)*xj(1,2)
      Tom(1,3) =       xj(1,3)*xj(1,3)
      Tom(1,4) =       xj(1,1)*xj(1,2)
      Tom(1,5) =       xj(1,1)*xj(1,3)
      Tom(1,6) =       xj(1,2)*xj(1,3)
      Tom(2,1) =       xj(2,1)*xj(2,1)
      Tom(2,2) =       xj(2,2)*xj(2,2)
      Tom(2,3) =       xj(2,3)*xj(2,3)
      Tom(2,4) =       xj(2,1)*xj(2,2)
      Tom(2,5) =       xj(2,1)*xj(2,3)
      Tom(2,6) =       xj(2,2)*xj(2,3)
      Tom(3,1) =       xj(3,1)*xj(3,1)
      Tom(3,2) =       xj(3,2)*xj(3,2)
      Tom(3,3) =       xj(3,3)*xj(3,3)
      Tom(3,4) =       xj(3,1)*xj(3,2)
      Tom(3,5) =       xj(3,1)*xj(3,3)
      Tom(3,6) =       xj(3,2)*xj(3,3)
      Tom(4,1) = 2.0d0*xj(1,1)*xj(2,1)
      Tom(4,2) = 2.0d0*xj(1,2)*xj(2,2)
      Tom(4,3) = 2.0d0*xj(1,3)*xj(2,3)
      Tom(4,4) =       xj(1,1)*xj(2,2) + xj(1,2)*xj(2,1)
      Tom(4,5) =       xj(1,1)*xj(2,3) + xj(1,3)*xj(2,1)
      Tom(4,6) =       xj(1,2)*xj(2,3) + xj(1,3)*xj(2,2)
      Tom(5,1) = 2.0d0*xj(1,1)*xj(3,1)
      Tom(5,2) = 2.0d0*xj(1,2)*xj(3,2)
      Tom(5,3) = 2.0d0*xj(1,3)*xj(3,3)
      Tom(5,4) =       xj(1,1)*xj(3,2) + xj(1,2)*xj(3,1)
      Tom(5,5) =       xj(1,1)*xj(3,3) + xj(1,3)*xj(3,1)
      Tom(5,6) =       xj(1,2)*xj(3,3) + xj(1,3)*xj(3,2)
      Tom(6,1) = 2.0d0*xj(2,1)*xj(3,1)
      Tom(6,2) = 2.0d0*xj(2,2)*xj(3,2)
      Tom(6,3) = 2.0d0*xj(2,3)*xj(3,3)
      Tom(6,4) =       xj(2,1)*xj(3,2) + xj(2,2)*xj(3,1)
      Tom(6,5) =       xj(2,1)*xj(3,3) + xj(2,3)*xj(3,1)
      Tom(6,6) =       xj(2,2)*xj(3,3) + xj(2,3)*xj(3,2)
C
c.....transformation matrix (convect. into cart. coos)
      do 20 l=1,3
         do 20 m=1,3
            to(l,m)=0.0d0
            do 20 i=1,3
               to(l,m)=to(l,m)+glo(i,l)*t(i,m)
20    continue
c     matrix Te=to(i,j)*to(k,l)
      Te(1,1) =       to(1,1)*to(1,1)
      Te(1,2) =       to(2,1)*to(2,1)
      Te(1,3) =       to(3,1)*to(3,1)
      Te(1,4) =       to(1,1)*to(2,1)
      Te(1,5) =       to(1,1)*to(3,1)
      Te(1,6) =       to(2,1)*to(3,1)
      Te(2,1) =       to(1,2)*to(1,2)
      Te(2,2) =       to(2,2)*to(2,2)
      Te(2,3) =       to(3,2)*to(3,2)
      Te(2,4) =       to(1,2)*to(2,2)
      Te(2,5) =       to(1,2)*to(3,2)
      Te(2,6) =       to(2,2)*to(3,2)
      Te(3,1) =       to(1,3)*to(1,3)
      Te(3,2) =       to(2,3)*to(2,3)
      Te(3,3) =       to(3,3)*to(3,3)
      Te(3,4) =       to(1,3)*to(2,3)
      Te(3,5) =       to(1,3)*to(3,3)
      Te(3,6) =       to(2,3)*to(3,3)
      Te(4,1) = 2.0d0*to(1,1)*to(1,2)
      Te(4,2) = 2.0d0*to(2,1)*to(2,2)
      Te(4,3) = 2.0d0*to(3,1)*to(3,2)
      Te(4,4) =       to(1,1)*to(2,2) + to(2,1)*to(1,2)
      Te(4,5) =       to(1,1)*to(3,2) + to(3,1)*to(1,2)
      Te(4,6) =       to(2,1)*to(3,2) + to(3,1)*to(2,2)
      Te(5,1) = 2.0d0*to(1,1)*to(1,3)
      Te(5,2) = 2.0d0*to(2,1)*to(2,3)
      Te(5,3) = 2.0d0*to(3,1)*to(3,3)
      Te(5,4) =       to(1,1)*to(2,3) + to(2,1)*to(1,3)
      Te(5,5) =       to(1,1)*to(3,3) + to(3,1)*to(1,3)
      Te(5,6) =       to(2,1)*to(3,3) + to(3,1)*to(2,3)
      Te(6,1) = 2.0d0*to(1,2)*to(1,3)
      Te(6,2) = 2.0d0*to(2,2)*to(2,3)
      Te(6,3) = 2.0d0*to(3,2)*to(3,3)
      Te(6,4) =       to(1,2)*to(2,3) + to(2,2)*to(1,3)
      Te(6,5) =       to(1,2)*to(3,3) + to(3,2)*to(1,3)
      Te(6,6) =       to(2,2)*to(3,3) + to(3,2)*to(2,3)
c
c.....Jacobi-Matrix (convect. into global cart. coordintes)
      call pzero(e,3*3)
      e(1,1)=1.0d0
      e(2,2)=1.0d0
      e(3,3)=1.0d0
      do 40 l=1,3
         do 40 m=1,3
         xjac(l,m)=0.0d0
            do 40 i=1,3
               xjac(l,m)=xjac(l,m)+glu(i,l)*e(i,m)
40    continue

c.....determinant
      detj1 = xjac(1,1)*(xjac(2,2)*xjac(3,3)-xjac(2,3)*xjac(3,2))
      detj2 =-xjac(1,2)*(xjac(2,1)*xjac(3,3)-xjac(2,3)*xjac(3,1))
      detj3 = xjac(1,3)*(xjac(2,1)*xjac(3,2)-xjac(2,2)*xjac(3,1))
      detj= detj1 + detj2 + detj3
      if(detj) 6,6,8
    6 write(iow,600) n
      stop
    8 continue

  600 format(1x,'program stop in shp45 of Elmt45',/,1x,
     +'zero or negative area for element number',i5)

c.....convect. coordinates moment.config gu(j,k) ; j=1,2,3=x,y,z ;k=1,2,3
c     covariant basis

      call pzero(xu,3*8)
      do 45 k=1,3
        do 45 l=1,8
          xu(k,l) = xl(k,l)+ul(k,l)
45    continue

      do 50 k=1,3
        do 50 j=1,3
          gu(j,k)=0.0d0
          do 50 ii=1,8
            gu(j,k)=gu(j,k)+shp(k,ii)*xu(j,ii)
50    continue

      do 55 k=1,3
        do 55 j=1,3
          eu(j,k)=0.0d0
          do 55 ii=1,8
            eu(j,k)=eu(j,k)+shp(k,ii)*ul(j,ii)
55    continue

      return
      end
c
      subroutine batdvo45(shpk,xl,ul,pkl,pk,eup,ndf)
c-----------------------------------------------------------------------
c.....Bathe-Dvorkin-formulation of the shear strains 13 / 23
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension shpk(3,8,4),xl(3,8),xu(3,8),pkl(3,3,4),pk(3,3,4),
     +          eup(3,3,4),ul(ndf,8)
c
      call pzero(xu,3*8)
      do 45 k=1,3
        do 45 l=1,8
          xu(k,l) = xl(k,l)+ul(k,l)
45    continue
c
c.....derivatives of the shape functions for shear strains in the
c     in collocation points
c
c.....point A
      xsia = -1.0d0
      etaa = 0.0d0
      zetaa = 0.0d0
      m = 1
      call bdshp45(xsia,etaa,zetaa,shpk,xl,xu,pkl,pk,ul,eup,m,ndf)
c
c.....point B
      xsib = 0.0d0
      etab = -1.0d0
      zetab = 0.0d0
      m = 2
      call bdshp45(xsib,etab,zetab,shpk,xl,xu,pkl,pk,ul,eup,m,ndf)
c
c.....point C
      xsic = 1.0d0
      etac = 0.0d0
      zetac = 0.0d0
      m = 3
      call bdshp45(xsic,etac,zetac,shpk,xl,xu,pkl,pk,ul,eup,m,ndf)
c
c.....point D
      xsid = 0.0d0
      etad = 1.0d0
      zetad = 0.0d0
      m = 4
      call bdshp45(xsid,etad,zetad,shpk,xl,xu,pkl,pk,ul,eup,m,ndf)
c
      return
      end
c
      subroutine bdshp45(xsix,etax,zeta,shpk,xl,xu,pkl,pk,ul,eup,m,ndf)
c-----------------------------------------------------------------------
c.....shape functions for shear strains
c.....convect. coordinates of the point of interpolation
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension shpk(3,8,4),xl(3,8),xu(3,8),pkl(3,3,4),pk(3,3,4),
     +        ul(ndf,8),eup(3,3,4)

      od=1.0d0
      ed=8.0d0

c.....shpk(a,b,m) is equivalent to shpkb,a,m with b=1 to 8 ; a=1 to 3
c     and m=1 to 4 for point A to C

      shpk(1,1,m) = (-od)*(od-etax)*(od-zeta)/ed
      shpk(1,2,m) = ( od)*(od-etax)*(od-zeta)/ed
      shpk(1,3,m) = ( od)*(od+etax)*(od-zeta)/ed
      shpk(1,4,m) = (-od)*(od+etax)*(od-zeta)/ed
      shpk(1,5,m) = (-od)*(od-etax)*(od+zeta)/ed
      shpk(1,6,m) = ( od)*(od-etax)*(od+zeta)/ed
      shpk(1,7,m) = ( od)*(od+etax)*(od+zeta)/ed
      shpk(1,8,m) = (-od)*(od+etax)*(od+zeta)/ed

      shpk(2,1,m) = (od-xsix)*(-od)*(od-zeta)/ed
      shpk(2,2,m) = (od+xsix)*(-od)*(od-zeta)/ed
      shpk(2,3,m) = (od+xsix)*( od)*(od-zeta)/ed
      shpk(2,4,m) = (od-xsix)*( od)*(od-zeta)/ed
      shpk(2,5,m) = (od-xsix)*(-od)*(od+zeta)/ed
      shpk(2,6,m) = (od+xsix)*(-od)*(od+zeta)/ed
      shpk(2,7,m) = (od+xsix)*( od)*(od+zeta)/ed
      shpk(2,8,m) = (od-xsix)*( od)*(od+zeta)/ed

      shpk(3,1,m) = (od-xsix)*(od-etax)*(-od)/ed
      shpk(3,2,m) = (od+xsix)*(od-etax)*(-od)/ed
      shpk(3,3,m) = (od+xsix)*(od+etax)*(-od)/ed
      shpk(3,4,m) = (od-xsix)*(od+etax)*(-od)/ed
      shpk(3,5,m) = (od-xsix)*(od-etax)*(+od)/ed
      shpk(3,6,m) = (od+xsix)*(od-etax)*(+od)/ed
      shpk(3,7,m) = (od+xsix)*(od+etax)*(+od)/ed
      shpk(3,8,m) = (od-xsix)*(od+etax)*(+od)/ed
c
c.....convect. coordinates
c
      do 10 k=1,3
        do 10 j=1,3
          pkl(j,k,m)=0.0d0
          do 10 ii=1,8
            pkl(j,k,m)=pkl(j,k,m)+shpk(k,ii,m)*xl(j,ii)
10    continue
c
      do 20 k=1,3
        do 20 j=1,3
          pk(j,k,m)=0.0d0
          do 20 ii=1,8
            pk(j,k,m)=pk(j,k,m)+shpk(k,ii,m)*xu(j,ii)
20    continue
c
      do 30 k=1,3
        do 30 j=1,3
          eup(j,k,m)=0.0d0
          do 30 ii=1,8
            eup(j,k,m)=eup(j,k,m)+shpk(k,ii,m)*ul(j,ii)
30    continue
c
      return
      end
c
      subroutine betste45(xl,ul,pa,pal,eua,shpka,ndf)
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension ai(8),pa(3,3,4),pal(3,3,4),eua(3,3,4),shpka(3,8,3,4),
     +          xl(3,8),ul(ndf,8),xu(3,8)
      data ai /-1.d0,1.d0,1.d0,-1.d0,-1.d0,-1.d0,1.d0,1.d0/

c     collocation points for E33
      do i=1,4
         xsic = ai(i)
         etac = ai(i+4)
         zetac= 0.d0
         m=3
         call bsshp45(xsic,etac,zetac,m,i,shpka)
      end do

      call pzero(xu,3*8)
      do 45 k=1,3
        do 45 l=1,8
          xu(k,l) = xl(k,l)+ul(k,l)
45    continue

c.... covariant base vectors in collocation planes
c     g^c_zeta
c     G^c_zeta
c     (du/dxsi)^c
      call pzero(pa,3*3*4)
      call pzero(pal,3*3*4)
      call pzero(eua,3*3*4)
      do 101 n=1,4
         do 101 j=1,3
            do 101 ii=1,8
               pa(j,3,n) = pa(j,3,n) + shpka(3,ii,3,n)*xu(j,ii)
               pal(j,3,n) = pal(j,3,n) + shpka(3,ii,3,n)*xl(j,ii)
               eua(j,3,n) = eua(j,3,n) + shpka(3,ii,3,n)*ul(j,ii)
101   continue
c
      return
      end
c
      subroutine bsshp45(xsi,eta,zeta,m,i,shpka)
c-----------------------------------------------------------------------
c.....shape functions for thickness strains
c.....convect. coordinates of the point of interpolation
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension shpka(3,8,3,4)

      od=1.0d0
      ed=8.0d0

c.....shpka(a,b,m,i)
c     a = 1-3 derivation to xsi, eta, zeta
c     b = 1-8 nodes
c     m = 1-3 collocation plane
c     i = 1-4 collocation points

      shpka(1,1,m,i) = (-od)*(od-eta)*(od-zeta)/ed
      shpka(1,2,m,i) = ( od)*(od-eta)*(od-zeta)/ed
      shpka(1,3,m,i) = ( od)*(od+eta)*(od-zeta)/ed
      shpka(1,4,m,i) = (-od)*(od+eta)*(od-zeta)/ed
      shpka(1,5,m,i) = (-od)*(od-eta)*(od+zeta)/ed
      shpka(1,6,m,i) = ( od)*(od-eta)*(od+zeta)/ed
      shpka(1,7,m,i) = ( od)*(od+eta)*(od+zeta)/ed
      shpka(1,8,m,i) = (-od)*(od+eta)*(od+zeta)/ed

      shpka(2,1,m,i) = (od-xsi)*(-od)*(od-zeta)/ed
      shpka(2,2,m,i) = (od+xsi)*(-od)*(od-zeta)/ed
      shpka(2,3,m,i) = (od+xsi)*( od)*(od-zeta)/ed
      shpka(2,4,m,i) = (od-xsi)*( od)*(od-zeta)/ed
      shpka(2,5,m,i) = (od-xsi)*(-od)*(od+zeta)/ed
      shpka(2,6,m,i) = (od+xsi)*(-od)*(od+zeta)/ed
      shpka(2,7,m,i) = (od+xsi)*( od)*(od+zeta)/ed
      shpka(2,8,m,i) = (od-xsi)*( od)*(od+zeta)/ed

      shpka(3,1,m,i) = (od-xsi)*(od-eta)*(-od)/ed
      shpka(3,2,m,i) = (od+xsi)*(od-eta)*(-od)/ed
      shpka(3,3,m,i) = (od+xsi)*(od+eta)*(-od)/ed
      shpka(3,4,m,i) = (od-xsi)*(od+eta)*(-od)/ed
      shpka(3,5,m,i) = (od-xsi)*(od-eta)*(+od)/ed
      shpka(3,6,m,i) = (od+xsi)*(od-eta)*(+od)/ed
      shpka(3,7,m,i) = (od+xsi)*(od+eta)*(+od)/ed
      shpka(3,8,m,i) = (od-xsi)*(od+eta)*(+od)/ed

      return
      end
c
      subroutine kine45(d,gu,glu,eu,xsi,eta,eas,e)
c----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(*),glu(3,3),gu(3,3),green(3,3),eas(6),e(6),
     +          shpk(3,8,4),pkl(3,3,4),pk(3,3,4),eup(3,3,4),eu(3,3),
     +          shpa(3,4),pa(3,3,4),pal(3,3,4),eua(3,3,4),shpka(3,8,3,4)
      common /batdvo/ shpk,pk,pkl,eup
!$OMP THREADPRIVATE (/batdvo/)
      common /betste/ pa,pal,shpa,shpka,eua
!$OMP THREADPRIVATE (/betste/)
c
      lin  = d(3)
      ieas = d(4)
      ibd  = d(5)
      ibs  = d(6)
c
c.....Green--Lagrange strains
c     nonlinear strains
      if (lin.eq.1) then
      do 10 k=1,3
         do 10 l=1,3
            green(k,l)=0.0d0
            do 10 j=1,3
               green(k,l) = green(k,l)
     +                    + (gu(j,k)*gu(j,l)-glu(j,k)*glu(j,l))*0.5d0
10    continue
      if (ibs.gt.0) then
      green(3,3)=0.0d0
         do 11 i=1,4
            do 11 j=1,3
               green(3,3) = green(3,3) + shpa(3,i) * 0.5d0 *
     +                      (pa(j,3,i)*pa(j,3,i)-pal(j,3,i)*pal(j,3,i))
11       continue
      end if
      if (ibd.gt.0) then
         green(1,3) = 0.0d0
         do i=1,3
            green(1,3)=green(1,3)+0.5d0*
     +          ((1-eta)*(pk(i,1,2)*pk(i,3,2)-pkl(i,1,2)*pkl(i,3,2))+
     +           (1+eta)*(pk(i,1,4)*pk(i,3,4)-pkl(i,1,4)*pkl(i,3,4)))
         end do
         green(1,3) = 0.5d0 * green(1,3)
         green(3,1) = green(1,3)
         green(2,3) = 0.0d0
         do i=1,3
            green(2,3)=green(2,3)+0.5d0*
     +          ((1-xsi)*(pk(i,2,1)*pk(i,3,1)-pkl(i,2,1)*pkl(i,3,1))+
     +           (1+xsi)*(pk(i,2,3)*pk(i,3,3)-pkl(i,2,3)*pkl(i,3,3)))
         end do
         green(2,3) = 0.5d0 * green(2,3)
         green(3,2) = green(2,3)
      end if
c.....linear strains
      else if (lin.eq.0) then
      do 20 k=1,3
         do 20 l=1,3
            green(k,l)=0.0d0
            do 20 j=1,3
               green(k,l) = green(k,l)
     +                    + (eu(j,k)*glu(j,l)+glu(j,k)*eu(j,l))*0.5d0
20    continue
      if (ibs.gt.0) then
      green(3,3)=0.0d0
         do 21 i=1,4
            do 21 j=1,3
               green(3,3) = green(3,3) + shpa(3,i) * 0.5d0 *
     +                     (eua(j,3,i)*pal(j,3,i)+eua(j,3,i)*pal(j,3,i))
21       continue
      end if
      if (ibd.gt.0) then
         green(1,3) = 0.0d0
         do i=1,3
            green(1,3)=green(1,3)+0.5d0*
     +          ((1-eta)*(eup(i,1,2)*pkl(i,3,2)+pkl(i,1,2)*eup(i,3,2))+
     +           (1+eta)*(eup(i,1,4)*pkl(i,3,4)+pkl(i,1,4)*eup(i,3,4)))
         end do
         green(1,3) = 0.5d0 * green(1,3)
         green(3,1) = green(1,3)
         green(2,3) = 0.0d0
         do i=1,3
            green(2,3)=green(2,3)+0.5d0*
     +          ((1-xsi)*(eup(i,2,1)*pkl(i,3,1)+pkl(i,2,1)*eup(i,3,1))+
     +           (1+xsi)*(eup(i,2,3)*pkl(i,3,3)+pkl(i,2,3)*eup(i,3,3)))
         end do
         green(2,3) = 0.5d0 * green(2,3)
         green(3,2) = green(2,3)
      end if
      end if
c
c.....ordering and enhanced strains
      do i=1,3
         e(i)=green(i,i)
      end do
      e(4)=2.d0*green(1,2)
      e(5)=2.d0*green(1,3)
      e(6)=2.d0*green(2,3)
      if (ieas.gt.0) then
         do i=1,6
            e(i)=e(i)+eas(i)
         end do
      end if
c
      return
      end
c
      subroutine bmat45(gu,glu,shp,ii,bm,xsi,eta,d)
c-----------------------------------------------------------------------
c.....B-Matrix
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension d(*),gu(3,3),glu(3,3),shp(4,8),bm(6,3),
     +          shpk(3,8,4),pk(3,3,4),pkl(3,3,4),eup(3,3,4),
     +          pa(3,3,4),pal(3,3,4),shpa(3,4),shpka(3,8,3,4),
     +          eua(3,3,4)
      common /batdvo/ shpk,pk,pkl,eup
!$OMP THREADPRIVATE (/batdvo/)
      common /betste/ pa,pal,shpa,shpka,eua
!$OMP THREADPRIVATE (/betste/)
      call pzero(bm,6*3)
c
c
      lin = d(3)
      ibd = d(5)
      ibs = d(6)
c
C.....standard B-matrix
      if (lin.eq.1) then
         bm(1,1)=gu(1,1)*shp(1,ii)
         bm(1,2)=gu(2,1)*shp(1,ii)
         bm(1,3)=gu(3,1)*shp(1,ii)
         bm(2,1)=gu(1,2)*shp(2,ii)
         bm(2,2)=gu(2,2)*shp(2,ii)
         bm(2,3)=gu(3,2)*shp(2,ii)
         bm(3,1)=gu(1,3)*shp(3,ii)
         bm(3,2)=gu(2,3)*shp(3,ii)
         bm(3,3)=gu(3,3)*shp(3,ii)
         bm(4,1)=gu(1,2)*shp(1,ii)+gu(1,1)*shp(2,ii)
         bm(4,2)=gu(2,2)*shp(1,ii)+gu(2,1)*shp(2,ii)
         bm(4,3)=gu(3,2)*shp(1,ii)+gu(3,1)*shp(2,ii)
         bm(5,1)=gu(1,3)*shp(1,ii)+gu(1,1)*shp(3,ii)
         bm(5,2)=gu(2,3)*shp(1,ii)+gu(2,1)*shp(3,ii)
         bm(5,3)=gu(3,3)*shp(1,ii)+gu(3,1)*shp(3,ii)
         bm(6,1)=gu(1,3)*shp(2,ii)+gu(1,2)*shp(3,ii)
         bm(6,2)=gu(2,3)*shp(2,ii)+gu(2,2)*shp(3,ii)
         bm(6,3)=gu(3,3)*shp(2,ii)+gu(3,2)*shp(3,ii)
      else if(lin.eq.0) then
         bm(1,1)=glu(1,1)*shp(1,ii)
         bm(1,2)=glu(2,1)*shp(1,ii)
         bm(1,3)=glu(3,1)*shp(1,ii)
         bm(2,1)=glu(1,2)*shp(2,ii)
         bm(2,2)=glu(2,2)*shp(2,ii)
         bm(2,3)=glu(3,2)*shp(2,ii)
         bm(3,1)=glu(1,3)*shp(3,ii)
         bm(3,2)=glu(2,3)*shp(3,ii)
         bm(3,3)=glu(3,3)*shp(3,ii)
         bm(4,1)=glu(1,2)*shp(1,ii)+glu(1,1)*shp(2,ii)
         bm(4,2)=glu(2,2)*shp(1,ii)+glu(2,1)*shp(2,ii)
         bm(4,3)=glu(3,2)*shp(1,ii)+glu(3,1)*shp(2,ii)
         bm(5,1)=glu(1,3)*shp(1,ii)+glu(1,1)*shp(3,ii)
         bm(5,2)=glu(2,3)*shp(1,ii)+glu(2,1)*shp(3,ii)
         bm(5,3)=glu(3,3)*shp(1,ii)+glu(3,1)*shp(3,ii)
         bm(6,1)=glu(1,3)*shp(2,ii)+glu(1,2)*shp(3,ii)
         bm(6,2)=glu(2,3)*shp(2,ii)+glu(2,2)*shp(3,ii)
         bm(6,3)=glu(3,3)*shp(2,ii)+glu(3,2)*shp(3,ii)
      end if
c
c.....Betsch-Stein interpolation E_(zeta,zeta)
      if (ibs.gt.0) then
         if(lin.eq.1)then
         bm(3,1) = 0.d0
         bm(3,2) = 0.d0
         bm(3,3) = 0.d0
         do i=1,4
            bm(3,1) = bm(3,1) + shpa(3,i)*pa(1,3,i)*shpka(3,ii,3,i)
            bm(3,2) = bm(3,2) + shpa(3,i)*pa(2,3,i)*shpka(3,ii,3,i)
            bm(3,3) = bm(3,3) + shpa(3,i)*pa(3,3,i)*shpka(3,ii,3,i)
         end do
         else if (lin.eq.0) then
         bm(3,1) = 0.d0
         bm(3,2) = 0.d0
         bm(3,3) = 0.d0
         do i=1,4
            bm(3,1) = bm(3,1) + shpa(3,i)*pal(1,3,i)*shpka(3,ii,3,i)
            bm(3,2) = bm(3,2) + shpa(3,i)*pal(2,3,i)*shpka(3,ii,3,i)
            bm(3,3) = bm(3,3) + shpa(3,i)*pal(3,3,i)*shpka(3,ii,3,i)
         end do
         end if
      end if
c
c.....Bathe-Dvorkin Interploation E_(xsi,zeta), E_(eta,zeta)
      if (ibd.gt.0) then
         if (lin.eq.1) then
         bm(5,1)=((1-eta)*
     +            (pk(1,3,2)*shpk(1,ii,2)+pk(1,1,2)*shpk(3,ii,2))+
     +            (1+eta)*
     +            (pk(1,3,4)*shpk(1,ii,4)+pk(1,1,4)*shpk(3,ii,4)))
     +            *0.5d0
         bm(5,2)=((1-eta)*
     +            (pk(2,3,2)*shpk(1,ii,2)+pk(2,1,2)*shpk(3,ii,2))+
     +            (1+eta)*
     +            (pk(2,3,4)*shpk(1,ii,4)+pk(2,1,4)*shpk(3,ii,4)))
     +            *0.5d0
         bm(5,3)=((1-eta)*
     +            (pk(3,3,2)*shpk(1,ii,2)+pk(3,1,2)*shpk(3,ii,2))+
     +            (1+eta)*
     +            (pk(3,3,4)*shpk(1,ii,4)+pk(3,1,4)*shpk(3,ii,4)))
     +            *0.5d0
         bm(6,1)=((1-xsi)*
     +            (pk(1,3,1)*shpk(2,ii,1)+pk(1,2,1)*shpk(3,ii,1))+
     +            (1+xsi)*
     +            (pk(1,3,3)*shpk(2,ii,3)+pk(1,2,3)*shpk(3,ii,3)))
     +            *0.5d0
         bm(6,2)=((1-xsi)*
     +            (pk(2,3,1)*shpk(2,ii,1)+pk(2,2,1)*shpk(3,ii,1))+
     +            (1+xsi)*
     +            (pk(2,3,3)*shpk(2,ii,3)+pk(2,2,3)*shpk(3,ii,3)))
     +            *0.5d0
         bm(6,3)=((1-xsi)*
     +            (pk(3,3,1)*shpk(2,ii,1)+pk(3,2,1)*shpk(3,ii,1))+
     +            (1+xsi)*
     +            (pk(3,3,3)*shpk(2,ii,3)+pk(3,2,3)*shpk(3,ii,3)))
     +            *0.5d0
         else if(lin.eq.0)then
         bm(5,1)=((1-eta)*
     +            (pkl(1,3,2)*shpk(1,ii,2)+pkl(1,1,2)*shpk(3,ii,2))+
     +            (1+eta)*
     +            (pkl(1,3,4)*shpk(1,ii,4)+pkl(1,1,4)*shpk(3,ii,4)))
     +            *0.5d0
         bm(5,2)=((1-eta)*
     +            (pkl(2,3,2)*shpk(1,ii,2)+pkl(2,1,2)*shpk(3,ii,2))+
     +            (1+eta)*
     +            (pkl(2,3,4)*shpk(1,ii,4)+pkl(2,1,4)*shpk(3,ii,4)))
     +            *0.5d0
         bm(5,3)=((1-eta)*
     +            (pkl(3,3,2)*shpk(1,ii,2)+pkl(3,1,2)*shpk(3,ii,2))+
     +            (1+eta)*
     +            (pkl(3,3,4)*shpk(1,ii,4)+pkl(3,1,4)*shpk(3,ii,4)))
     +            *0.5d0
         bm(6,1)=((1-xsi)*
     +            (pkl(1,3,1)*shpk(2,ii,1)+pkl(1,2,1)*shpk(3,ii,1))+
     +            (1+xsi)*
     +            (pkl(1,3,3)*shpk(2,ii,3)+pkl(1,2,3)*shpk(3,ii,3)))
     +            *0.5d0
         bm(6,2)=((1-xsi)*
     +            (pkl(2,3,1)*shpk(2,ii,1)+pkl(2,2,1)*shpk(3,ii,1))+
     +            (1+xsi)*
     +            (pkl(2,3,3)*shpk(2,ii,3)+pkl(2,2,3)*shpk(3,ii,3)))
     +            *0.5d0
         bm(6,3)=((1-xsi)*
     +            (pkl(3,3,1)*shpk(2,ii,1)+pkl(3,2,1)*shpk(3,ii,1))+
     +            (1+xsi)*
     +            (pkl(3,3,3)*shpk(2,ii,3)+pkl(3,2,3)*shpk(3,ii,3)))
     +            *0.5d0
         end if
      end if
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gmat45(shp,ii,jj,sig,gmat,xsi,eta,d)
c-----------------------------------------------------------------------
c.....G-Matrix
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension shp(4,8),sig(6),gmat(3,3),d(*),
     +          shpk(3,8,4),pkl(3,3,4),pk(3,3,4),eup(3,3,4),
     +          shpa(3,4),pa(3,3,4),pal(3,3,4),eua(3,3,4),shpka(3,8,3,4)
      common /batdvo/ shpk,pk,pkl,eup
!$OMP THREADPRIVATE (/batdvo/)
      common /betste/ pa,pal,shpa,shpka,eua
!$OMP THREADPRIVATE (/betste/)
c
      ibd = d(5)
      ibs = d(6)
      fm1 = 0.d0
      fm2 = 0.d0
c
      fm1 = shp(1,ii)*sig(1)*shp(1,jj) + shp(2,ii)*sig(2)*shp(2,jj)
      fm2 = shp(3,ii)*sig(3)*shp(3,jj)
c
      if (ibs.gt.0) then
         fm2 = 0
         do i=1,4
            fm2 = fm2 + shpa(3,i)*shpka(3,ii,3,i)*shpka(3,jj,3,i)*sig(3)
         end do
      end if
c
      fs1 = (shp(1,ii)*shp(2,jj) + shp(2,ii)*shp(1,jj)) *sig(4)
      fs2 = (shp(1,ii)*shp(3,jj) + shp(3,ii)*shp(1,jj)) *sig(5)
     +    + (shp(2,ii)*shp(3,jj) + shp(3,ii)*shp(2,jj)) *sig(6)
c
      if (ibd.gt.0) then
         fs2 =((1-eta)*(shpk(1,ii,2)*shpk(3,jj,2)
     +                + shpk(3,ii,2)*shpk(1,jj,2))
     +      +  (1+eta)*(shpk(1,ii,4)*shpk(3,jj,4)
     +                + shpk(3,ii,4)*shpk(1,jj,4))) * 0.50d0 * sig(5)
     +      + ((1-xsi)*(shpk(2,ii,1)*shpk(3,jj,1)
     +                + shpk(3,ii,1)*shpk(2,jj,1))
     +      +  (1+xsi)*(shpk(2,ii,3)*shpk(3,jj,3)
     +                + shpk(3,ii,3)*shpk(2,jj,3))) * 0.50d0 * sig(6)
      end if
      f=fm1+fm2+fs1+fs2
c
      call pzero(gmat,3*3)
      gmat(1,1)=f
      gmat(2,2)=f
      gmat(3,3)=f
c
      return
      end
c
      subroutine strepl45(ix,dt,st,plout,sig,shp,nel,numnp,dv)
c-----------------------------------------------------------------------
c.....plot stresses for nonlinear 3D--element
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ix(*),dt(numnp),st(numnp,*),shp(4,8),sig(7),plout(10)
      do i = 1,nel
        xsji = dv*shp(4,i)
        ii = abs(ix(i))
        if (ii.eq.0) goto 10
        dt(ii) = dt(ii) + xsji
        st(ii,1) = st(ii,1) + sig(1)*xsji
        st(ii,2) = st(ii,2) + sig(2)*xsji
        st(ii,3) = st(ii,3) + sig(3)*xsji
        st(ii,4) = st(ii,4) + sig(4)*xsji
        st(ii,5) = st(ii,5) + sig(5)*xsji
        st(ii,6) = st(ii,6) + sig(6)*xsji
        do k = 1,10
          st(ii,15+k) = st(ii,15+k) + plout(k)*xsji
        end do
10      continue
      end do
      return
      end
c
      subroutine decond45(skcr,prc,alphao,dv,dvl,alpha,neas,nseq,nges)
c----------------------------------------------------------------------
c     decondensation
c     function : Dc=-Kcc^-1(Kcr Dr - Rc)
c     (see R.D. Cook, page 187)
c
c     Kcr(neas,nseq)   Kcc(neas,neas)   --> scr(neas,nseq+neas)
c                      prc(neas)
c     dv(nseq)         dalpha(neas)     --> dvl(nseq+neas)
c----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension skcr(neas,neas+nseq),prc(*),dv(*),dvl(*),alphao(*),
     +          alpha(*)

c.... test if first time (K_cc=0)
      sum = 0.0d0
      tol45 = 1.d-10
      do i = 1,neas
         j = nseq+i
         dsum = skcr(i,j)
          sum = sum + dsum*dsum
      end do
      if(sum.lt.tol45)  return

      do i = 1,nseq
         dvl(i) = dv(i)
      end do
c
c.... test if first (linear) step in time (alpha=0)
      sum = 0.0d0
      tol451 = 1.d-50
      do i = 1,nseq
          sum = sum + dvl(i)*dvl(i)
      end do
      if (sum.lt.tol451) then
         do i=1,neas
            alpha(i) = alphao(i)
         end do
         goto 100
      end if
c
      do 10 ieq = 1,neas
       keq = nseq + ieq
       dum = 0.d0
       do 20 irow = 1,keq-1
        dum = dum + skcr(ieq,irow)*dvl(irow)
20     continue
       if(dabs(skcr(ieq,keq)).le.1.e-10) write(*,*)'D=0 in DECOND',keq
c
c...   store into long vector  (only local!)
        dvl(keq) = (prc(ieq) - dum) / skcr(ieq,keq)
c
c...   update for alpha
        alpha(ieq) = alphao(ieq) + dvl(keq)
10    continue
100   continue
      return
      end
c
      subroutine store45(sc,pc,alpha,skcr,prc,alphao,nst,neas,nges)
c----------------------------------------------------------------------
c     Store matrices for iteration
c
c        (Kcr+Kcc), Rc
c
c  Dim:  Krc(neas,nst) Kcc(neas,neas) together!
c        Rc (neas)
c        alpha(neas), alphao(neas)
c----------------------------------------------------------------------
      implicit  double precision (a-h,o-z)
      dimension sc(nges,*),pc(*),skcr(neas,nst+neas),prc(neas),
     1          alphao(*),alpha(*)
      do i = 1,neas
        ii = nst+i
        prc(i) = pc(ii)
        alphao(i) = alpha(i)
        do k = 1, neas+nst
          skcr(i,k) = sc(ii,k)
        end do
      end do
      return
      end
c
      subroutine mmat45(xM,dj,xsi,eta,zeta,d,To,dj0)
c-----------------------------------------------------------------------
c     enhanced strain interpolation matrix xM
c     (c) S.Klinkel 5/97
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension xMxi(6,30),xM(6,*),d(*),To(6,6),ToM(6,30)
c
      neas = d(4)
      ieas = d(7)
      call pzero(xMxi,neas*6)
      call pzero(xM,neas*6)
c
      if(neas.eq.3) then                    ! (Rah et al. 2011) ok for
      xMxi(3,1) = zeta                      ! bending patch test
      xMxi(3,2) = zeta*xsi
      xMxi(3,3) = zeta*eta
c
      else if(neas.eq.5) then               ! (sk,fg,ww 1999) ok for
      xMxi(1,1) = xsi                       ! membrane patch test
      xMxi(2,2) = eta
      xMxi(3,3) = zeta
      xMxi(4,4) = xsi
      xMxi(4,5) = eta
C
      else if (neas.eq.8) then              ! (sk,fg,ww 1999)
      xMxi(1,1) = xsi                       ! & dist. mesh
      xMxi(2,2) = eta
      xMxi(3,3) = zeta
      xMxi(1,4) = xsi*eta
      xMxi(2,5) = xsi*eta
      xMxi(4,6) = xsi
      xMxi(4,7) = eta
      xMxi(4,8) = xsi*eta
c
      else if (neas.eq.11) then             ! (sk,fg,ww 1999)
      xMxi(1,1) = xsi                       ! & bending patch test
      xMxi(2,2) = eta                       ! & vol. locking reg and irreg mesh
      xMxi(3,3) = zeta
      xMxi(1,4) = xsi*eta
      xMxi(2,5) = xsi*eta
      xMxi(3,6) = zeta*xsi
      xMxi(3,7) = zeta*eta
      xMxi(3,8) = zeta*xsi*eta
      xMxi(4,9) = xsi
      xMxi(4,10) = eta
      xMxi(4,11) = xsi*eta
c
      else if (neas.eq.7.and.ieas.eq.0) then ! (Vu-Quoc 2003) ok for
      xMxi(1,1) = xsi                        ! bending patch test
      xMxi(2,2) = eta                        ! but vol. locking
      xMxi(3,3) = zeta
      xMxi(3,4) = xsi*zeta
      xMxi(3,5) = eta*zeta
      xMxi(4,6) = xsi
      xMxi(4,7) = eta
c
      else if (neas.eq.7.and.ieas.eq.1) then ! (Rah et al. 2011) ok for
      xMxi(1,1) = xsi                        ! vol.locking only reg. mesh
      xMxi(1,2) = xsi*eta                    ! but locking for irreg. mesh
      xMxi(2,3) = eta
      xMxi(2,4) = xsi*eta
      xMxi(3,5) = zeta
      xMxi(3,6) = zeta*xsi
      xMxi(3,7) = zeta*eta
c
      else if (neas.eq.7.and.ieas.eq.2) then ! (Rah et al. 2011)
      xMxi(1,1) = xsi*zeta                   ! dist. mesh bending
      xMxi(2,2) = eta*zeta
      xMxi(3,3) = zeta
      xMxi(3,4) = zeta*xsi
      xMxi(3,5) = zeta*eta
      xMxi(4,6) = xsi*zeta
      xMxi(4,7) = eta*zeta
c
      else if (neas.eq.4) then
      xMxi(1,1) = xsi
      xMxi(2,2) = eta
      xMxi(3,3) = xsi
      xMxi(3,4) = eta
C
      else if (neas.eq.9) then
      xMxi(3,1) = zeta
      xMxi(5,2) = xsi
      xMxi(5,3) = zeta
      xMxi(6,4) = eta
      xMxi(6,5) = zeta
      xMxi(5,6) = xsi
      xMxi(5,7) = zeta
      xMxi(6,8) = eta
      xMxi(6,9) = zeta
c
      else if (neas.eq.12) then
      xMxi(1,1) = xsi
      xMxi(2,2) = eta
      xMxi(3,3) = zeta
C
      xMxi(4,4) = xsi
      xMxi(4,5) = eta
      xMxi(4,6) = xsi*eta
      xMxi(5,7) = xsi
      xMxi(5,8) = zeta
      xMxi(5,9) = xsi*zeta
      xMxi(6,10) = eta
      xMxi(6,11) = zeta
      xMxi(6,12) = eta*zeta
c
      else if(neas.eq.30) then

      xMxi(1,1) = xsi
      xMxi(2,2) = eta
      xMxi(3,3) = zeta

      xMxi(1,4) = xsi*eta
      xMxi(1,5) = xsi*zeta
      xMxi(2,6) = xsi*eta
      xMxi(2,7) = eta*zeta
      xMxi(3,8) = xsi*zeta
      xMxi(3,9) = eta*zeta

      xMxi(4,10) = xsi
      xMxi(4,11) = eta
      xMxi(5,12) = xsi
      xMxi(5,13) = zeta
      xMxi(6,14) = eta
      xMxi(6,15) = zeta

      xMxi(4,16) = xsi*zeta
      xMxi(4,17) = eta*zeta
      xMxi(5,18) = xsi*eta
      xMxi(5,19) = zeta*eta
      xMxi(6,20) = eta*xsi
      xMxi(6,21) = zeta*xsi

      xMxi(4,22) = xsi*eta
      xMxi(5,23) = xsi*zeta
      xMxi(6,24) = eta*zeta

      xMxi(4,25) = xsi*eta*zeta
      xMxi(5,26) = xsi*eta*zeta
      xMxi(6,27) = xsi*eta*zeta
      xMxi(1,28) = xsi*eta*zeta
      xMxi(2,29) = xsi*eta*zeta
      xMxi(3,30) = xsi*eta*zeta
c
      end if
c
      Do i=1,6
         Do j=1,neas
            ToM(i,j) = 0.0d0
            Do k=1,6
               ToM(i,j) = ToM(i,j) + To(i,k) * xMxi(k,j)
            end do
            xM(i,j) = ToM(i,j) * dj0/dj
            if (abs(xM(i,j)).lt.1.e-15) xM(i,j)=0.0d0
         end do
      end do
c
      return
      end
c
      subroutine shearfac45(d,xl,cappa,wcappa)
c-----------------------------------------------------------------------
c.... shear correction factor
c     d(10) <0: cappa = input value*SCF(Tessler)
c     d(10) >0: cappa = input value
c     d(10) =0: cappa = 1 (default)
c
c     Assumptions:
c     A=(L_max)^2, L only on bottom calculated
c     h=total thickness, not valid for more than 1 element in 3-dir
c     nu=0
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  d(*),xl(3,*),yl(3)

      cappa = d(10)

      if(cappa.lt.0.d0)then
c....   AE... element area
c       ae = 4.d0*detj0  ! not available
c....   AE... square of longest element side-only on bottom of element
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
        hs    = d(8)
        xnu   = 0.d0
        cappa  = dabs(d(10))/(1.d0+ae/(2.d0*hs*hs*(1.d0+xnu))) ! Tessler
      end if

      wcappa=sqrt(cappa)

      return
      end
c

      subroutine trafo45_shear(cappa,wcappa,eps,cmat,sig,isw)
c-----------------------------------------------------------------------
c
c      Purpose:
c          multiply E_S, C, S_S with shear correction factor
c
c      Input:
c             cmat(6,6) - material matrix (cartesian coor. system)
c             sig(6)    - stress vector   (cartesian coor. system)
c             cappa     - shear correction factor
c             isw       - 1: E   2:C,S
c
c      Output:
c       1     eps(6)    - modified   strain vector
c       2     cmat(6,6) - modified material matrix
c       2     sig(6)    - modified   stress vector
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension eps(6),cmat(6,6),sig(6)

      if(isw.eq.1) then
c....  E-shear part
       eps(5)  = eps(5)*wcappa
       eps(6)  = eps(6)*wcappa

      else if(isw.eq.2) then
c....  C-shear part
       cmat(5,5) = cmat(5,5)*cappa
       cmat(5,6) = cmat(5,6)*cappa
       cmat(6,6) = cmat(6,6)*cappa

c....  C-mixed part
       cmat(1,5) = cmat(1,5)*wcappa
       cmat(1,6) = cmat(1,6)*wcappa
       cmat(2,5) = cmat(2,5)*wcappa
       cmat(2,6) = cmat(2,6)*wcappa
       cmat(3,5) = cmat(3,5)*wcappa
       cmat(3,6) = cmat(3,6)*wcappa
       cmat(4,5) = cmat(4,5)*wcappa
       cmat(4,6) = cmat(4,6)*wcappa

c....  C-Symmetry
       cmat(5,1) = cmat(1,5)
       cmat(6,1) = cmat(1,6)
       cmat(5,2) = cmat(2,5)
       cmat(6,2) = cmat(2,6)
       cmat(5,3) = cmat(3,5)
       cmat(6,3) = cmat(3,6)
       cmat(5,4) = cmat(4,5)
       cmat(6,4) = cmat(4,6)

       cmat(6,5) = cmat(5,6)

c....  S-shear part
       sig(5) = sig(5)*wcappa
       sig(6) = sig(6)*wcappa
      end if
      return
      end
c
      subroutine trafo45_local(Eps,glu,Te,Sig,Cmat,iswm)
c-----------------------------------------------------------------------
c
c      Purpose:
c          iswm=1: transformation of Eps in cartesian coordinates
c          iswm=2: transformation of Sig, Cmat in convective coordinates
c
c      Input:
c          if iswm=1
c             Eps(6)   - strain components (11,22,33,12,13,23)
c             glu(3,3) - base vectors of Eps in columns
c          if iswm=2
c             Te(6,6)  - transformation matrix
c             Sig(6)   - stress vector in cartesian coor. system
c             Cmat(6,6)- material matrix in cartesian coor. system
c
c      Output:
c          if iswm=1
c             Eps(6)   - strain components in cartesian coor. system
c             Te(6,6)  - transformation matrix
c          if iswm=2
c             Sig(6)   - stress vector in convective coor. system
c             Cmat(6,6)- material matrix in convective coor. system
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension Eps(6),glu(3,3),Sig(6),Cmat(6,6),tik(3,3),t(3,3)
      dimension Ecar(6),gko(3,3),glo(3,3),Te(6,6),Ccon(6,6),Scon(6)
c
c.... transform Eps to local cartesian coordinate system
      if (iswm.eq.1) then
c       local orthonormal basis
        call vecp (glu(1,1),glu(1,2),t(1,3))
        call norm (t(1,3),t(1,3),3)
        call norm (t(1,1),glu(1,1),3)
        call vecp (t(1,3),t(1,1),t(1,2))

c       contravariant metric coefficient matrix
        do 15 i=1,3
           do 15 j=1,3
              gko(i,j) = 0.d0
              do 15 k=1,3
              gko(i,j) = gko(i,j) + glu(k,i) * glu(k,j)
15      continue
        call invert(gko,3,3)

c       contravariant basis in columns
        do i=1,3
          glo(i,1)=gko(1,1)*glu(i,1)+gko(1,2)*glu(i,2)+gko(1,3)*glu(i,3)
          glo(i,2)=gko(2,1)*glu(i,1)+gko(2,2)*glu(i,2)+gko(2,3)*glu(i,3)
          glo(i,3)=gko(3,1)*glu(i,1)+gko(3,2)*glu(i,2)+gko(3,3)*glu(i,3)
        end do

c       T_ij = G^i . t_j
        do i = 1,3
         do j = 1,3
          tik(i,j) = dot(glo(1,i),t(1,j),3)
         end do
        end do

c       matrix Te=T_mn(i,k)*T_mn(j,l)
        Te(1,1) =       tik(1,1)*tik(1,1)
        Te(1,2) =       tik(2,1)*tik(2,1)
        Te(1,3) =       tik(3,1)*tik(3,1)
        Te(1,4) =       tik(1,1)*tik(2,1)
        Te(1,5) =       tik(1,1)*tik(3,1)
        Te(1,6) =       tik(2,1)*tik(3,1)
        Te(2,1) =       tik(1,2)*tik(1,2)
        Te(2,2) =       tik(2,2)*tik(2,2)
        Te(2,3) =       tik(3,2)*tik(3,2)
        Te(2,4) =       tik(1,2)*tik(2,2)
        Te(2,5) =       tik(1,2)*tik(3,2)
        Te(2,6) =       tik(2,2)*tik(3,2)
        Te(3,1) =       tik(1,3)*tik(1,3)
        Te(3,2) =       tik(2,3)*tik(2,3)
        Te(3,3) =       tik(3,3)*tik(3,3)
        Te(3,4) =       tik(1,3)*tik(2,3)
        Te(3,5) =       tik(1,3)*tik(3,3)
        Te(3,6) =       tik(2,3)*tik(3,3)
        Te(4,1) = 2.0d0*tik(1,1)*tik(1,2)
        Te(4,2) = 2.0d0*tik(2,1)*tik(2,2)
        Te(4,3) = 2.0d0*tik(3,1)*tik(3,2)
        Te(4,4) =       tik(1,1)*tik(2,2) + tik(2,1)*tik(1,2)
        Te(4,5) =       tik(1,1)*tik(3,2) + tik(3,1)*tik(1,2)
        Te(4,6) =       tik(2,1)*tik(3,2) + tik(3,1)*tik(2,2)
        Te(5,1) = 2.0d0*tik(1,1)*tik(1,3)
        Te(5,2) = 2.0d0*tik(2,1)*tik(2,3)
        Te(5,3) = 2.0d0*tik(3,1)*tik(3,3)
        Te(5,4) =       tik(1,1)*tik(2,3) + tik(2,1)*tik(1,3)
        Te(5,5) =       tik(1,1)*tik(3,3) + tik(3,1)*tik(1,3)
        Te(5,6) =       tik(2,1)*tik(3,3) + tik(3,1)*tik(2,3)
        Te(6,1) = 2.0d0*tik(1,2)*tik(1,3)
        Te(6,2) = 2.0d0*tik(2,2)*tik(2,3)
        Te(6,3) = 2.0d0*tik(3,2)*tik(3,3)
        Te(6,4) =       tik(1,2)*tik(2,3) + tik(2,2)*tik(1,3)
        Te(6,5) =       tik(1,2)*tik(3,3) + tik(3,2)*tik(1,3)
        Te(6,6) =       tik(2,2)*tik(3,3) + tik(3,2)*tik(2,3)

c       transformation from the isopar. in kartes. space
        Ecar = matmul(Te,Eps)
        Eps = Ecar
c.... transform Sig Cmat in convective coordinates
      else if (iswm.eq.2) then
        Ccon = matmul(matmul(transpose(Te),Cmat),Te)
        Scon = matmul(transpose(Te),Sig)
        Cmat = Ccon
        Sig = Scon
      end if
c
      return
      end
c
      subroutine trafo45_global(Eps,glu,Te,Sig,Cmat,iswm)
c-----------------------------------------------------------------------
c
c      Purpose:
c          iswm=1: transformation of Eps in cartesian coordinates
c          iswm=2: transformation of Sig, Cmat in convective coordinates
c
c      Input:
c          if iswm=1
c             Eps(6)   - strain components (11,22,33,12,13,23)
c             glu(3,3) - base vectors of Eps in columns
c          if iswm=2
c             Te(6,6)  - transformation matrix
c             Sig(6)   - stress vector in cartesian coor. system
c             Cmat(6,6)- material matrix in cartesian coor. system
c
c      Output:
c          if iswm=1
c             Eps(6)   - strain components in cartesian coor. system
c             Te(6,6)  - transformation matrix
c          if iswm=2
c             Sig(6)   - stress vector in convective coor. system
c             Cmat(6,6)- material matrix in convective coor. system
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension Eps(6),glu(3,3),Sig(6),Cmat(6,6)
      dimension Ecar(6),gko(3,3),glo(3,3),Te(6,6),Ccon(6,6),Scon(6)
c
c.... transform Eps to global cartesian coordinate system
      if (iswm.eq.1) then
c       contravariant metric coefficient matrix
        do 15 i=1,3
           do 15 j=1,3
              gko(i,j) = 0.d0
              do 15 k=1,3
              gko(i,j) = gko(i,j) + glu(k,i) * glu(k,j)
15      continue
        call invert(gko,3,3)
c       contravariant basis in rows
        do i=1,3
          glo(1,i)=gko(1,1)*glu(i,1)+gko(1,2)*glu(i,2)+gko(1,3)*glu(i,3)
          glo(2,i)=gko(2,1)*glu(i,1)+gko(2,2)*glu(i,2)+gko(2,3)*glu(i,3)
          glo(3,i)=gko(3,1)*glu(i,1)+gko(3,2)*glu(i,2)+gko(3,3)*glu(i,3)
        end do
c       matrix Te=glo(i,j)*glo(k,l)
        Te(1,1) =       glo(1,1)*glo(1,1)
        Te(1,2) =       glo(2,1)*glo(2,1)
        Te(1,3) =       glo(3,1)*glo(3,1)
        Te(1,4) =       glo(1,1)*glo(2,1)
        Te(1,5) =       glo(1,1)*glo(3,1)
        Te(1,6) =       glo(2,1)*glo(3,1)
        Te(2,1) =       glo(1,2)*glo(1,2)
        Te(2,2) =       glo(2,2)*glo(2,2)
        Te(2,3) =       glo(3,2)*glo(3,2)
        Te(2,4) =       glo(1,2)*glo(2,2)
        Te(2,5) =       glo(1,2)*glo(3,2)
        Te(2,6) =       glo(2,2)*glo(3,2)
        Te(3,1) =       glo(1,3)*glo(1,3)
        Te(3,2) =       glo(2,3)*glo(2,3)
        Te(3,3) =       glo(3,3)*glo(3,3)
        Te(3,4) =       glo(1,3)*glo(2,3)
        Te(3,5) =       glo(1,3)*glo(3,3)
        Te(3,6) =       glo(2,3)*glo(3,3)
        Te(4,1) = 2.0d0*glo(1,1)*glo(1,2)
        Te(4,2) = 2.0d0*glo(2,1)*glo(2,2)
        Te(4,3) = 2.0d0*glo(3,1)*glo(3,2)
        Te(4,4) =       glo(1,1)*glo(2,2) + glo(2,1)*glo(1,2)
        Te(4,5) =       glo(1,1)*glo(3,2) + glo(3,1)*glo(1,2)
        Te(4,6) =       glo(2,1)*glo(3,2) + glo(3,1)*glo(2,2)
        Te(5,1) = 2.0d0*glo(1,1)*glo(1,3)
        Te(5,2) = 2.0d0*glo(2,1)*glo(2,3)
        Te(5,3) = 2.0d0*glo(3,1)*glo(3,3)
        Te(5,4) =       glo(1,1)*glo(2,3) + glo(2,1)*glo(1,3)
        Te(5,5) =       glo(1,1)*glo(3,3) + glo(3,1)*glo(1,3)
        Te(5,6) =       glo(2,1)*glo(3,3) + glo(3,1)*glo(2,3)
        Te(6,1) = 2.0d0*glo(1,2)*glo(1,3)
        Te(6,2) = 2.0d0*glo(2,2)*glo(2,3)
        Te(6,3) = 2.0d0*glo(3,2)*glo(3,3)
        Te(6,4) =       glo(1,2)*glo(2,3) + glo(2,2)*glo(1,3)
        Te(6,5) =       glo(1,2)*glo(3,3) + glo(3,2)*glo(1,3)
        Te(6,6) =       glo(2,2)*glo(3,3) + glo(3,2)*glo(2,3)
c       transformation from the isopar. in kartes. space
        Ecar = matmul(Te,Eps)
        Eps = Ecar
c.... transform Sig Cmat in convective coordinates
      else if (iswm.eq.2) then
        Ccon = matmul(matmul(transpose(Te),Cmat),Te)
        Scon = matmul(transpose(Te),Sig)
        Cmat = Ccon
        Sig = Scon
      end if
c
      return
      end

