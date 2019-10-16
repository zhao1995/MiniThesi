      subroutine elmt09(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c---------------------------------------------------------------------+
c
c.... quadrilateral shell element for feap                            |
c                                                                     |
c---------------------------------------------------------------------+
c.... incorporating membrane with normal drilling dof                 |
c.... modified to remove effects of constant strains and              |
c.... including deep shell curvature corrections to the               |
c.... discrete kirchhoff quadrilateral plate bending element          |
c     R.L.Taylor                                                      |
c---------------------------------------------------------------------+
c...  * TANG, LMAS, elastic foundation,                               |      
c     * PLOT_STRE                                                     |
c       klay = 0        klay = 0+d(21).ne.0      klay = -1            |
c       1:  N_x         1:  N_x                  1:  S_11(T)          |            
c       2:  N_xy        2:  N_xy                 2:  S_12(T)          |          
c       3:  N_y         3:  N_y                  3:  S_22(T)          |            
c       4:  N_1         4:  as_x bot             4:  S_1 (T)          |      
c       5:  N_2         5:  as_x top             5:  S_2 (T)          |       
c       6:  M_x         6:  M_x                  6:  S_11(B)          |       
c       7:  M_xy        7:  M_xy                 7:  S_12(B)          |      
c       8:  M_y         8:  M_y                  8:  S_22(B)          |      
c       9:  M_1         9:  as_y bot             9:  S_1 (B)          |      
c      10:  M_2        10:  as_y top            10:  S_2 (B)          |
C      11:  press      11:  press               11:  press            |
c       plot for klay =1                                              |
c       -> plot,stre,x,,-1 or macr,stre,lay,a,e,-1                    |
c                                                                     |
c---------------------------------------------------------------------+
c.... Input parameters set as follows:                                |
c                                                                     |
c         ndm = 3 (x,y,z cartesian coordinates at nodes)              |
c         ndf = 6 (u-x,u-y,u-z,r-x,r-y,r-z at nodes)                  |
c         nen = 4 nodes (counterclockwise around element)             |
c                                                                     |
c.... Material parameters                                             |
c                                                                     |
c       Record 1.                                                     |
c         E   = Young's modulus                                       |
c         nu  = Poisson's ratio                                       |
c         rho = mass density per unit volume                          |
c         h   = thickness of shell                                    |
c                                                                     |
c       Record 2. - loading parameters                                |
c         b-1 = uniform loading in 1-direction (per unit area)        |
c         b-2 = uniform loading in 2-direction (per unit area)        |
c         p   = pressure loading normal to shell (per unit area)      |
c         b-x = gravity loading in x-direction (per unit area)        |
c         b-y = gravity loading in y-direction (per unit area)        |
c         b-z = gravity loading in z-direction (per unit area)        |
c         tc  = load type (0=lumped; 1=consistent)                    |
c         cc  = elastic foundation constant                           |
c         pz  = 0:p=p  1=p=p*z                                        |
c                                                                     |
c       Record 3. - base system                                       |
c       Input data 3. line - basis for stress output, see pdirec3     |        
c       d(17)  type                   | d(4) | d(5) | d(6)            | 
c       1      basic via diagonals    |      |      |                 |  
c       2      t_3(1)+given vector t_1|  t1x |  t1y |  t1z            |  
c       3 - 17 see macro BASE         |      |      |                 |
c                                                                     |
c       control correct base system via plot,forc,1                   |  
c                                                                     |
c       Record 4. - concrete design                                   |
c        beta_R                                                       |
c        beta_S                                                       |
c        h_1                                                          |
c        h_2                                                          |
c                                                                     |
c---------------------------------------------------------------------+
c                                                                     |
c       Use of D-array                                                |
c                                                                     |
c       d( 1) =       |          Eh /(1-nu**2)                        |
c       d( 2) =       |       nu*Eh /(1-nu**2)                        |
c       d( 3) =       | (1-nu)/2*Eh /(1-nu**2)                        |
c       d( 4) =  t1   |                                               |
c       d( 5) =  t2   |                                               |
c       d( 6) =  t3   |                                               |
c       d( 7) =  b1   |                                               |
c       d( 8) =  b2   |                                               |
c       d( 9) =  p=b3 |                                               |
c       d(10) =  bx   |                                               |
c       d(11) =  by   |                                               |
c       d(12) =  bz   |                                               |
c       d(13) =  tc   | rho*h                                         |
c       d(14) =  cc   | rho*h**3/12                                   |
c       d(15) =  pz   | h**2/12                                       |
c       d(16) =       | tc                                            |
c       d(17) =  igeo |                                                |
c       d(18) =       | cc                                            |
c       d(19) =       | pz                                            |
c                                                                     |
c       d(20) =   h                                                   |
c       d(21) = betaR                                                 |
c       d(22) = betaS                                                 |
c       d(23) = h_1                                                   |
c       d(24) = h_2                                                   |
c---------------------------------------------------------------------+
c     Input for MACRO QLOA                                            | 
c        q01 =  b1   | see d( 7)                                      |
c        q02 =  b2   | see d( 8)                                      |
c        q03 =  p=b3 | see d( 9)                                      |
c        q04 =  bx   | see d(10)                                      |
c        q05 =  by   | see d(11)                                      |
c        q06 =  bz   | see d(12)                                      |
c---------------------------------------------------------------------+
c open: #calculate stress at klay.mlay instead of T/B                 |
c       #Lastfall gamma*z ist provisorisch eingebaut über d(19) und   |  
c        nur für Last über MATE, weiterhin ist die Nullinie auf zz=20 |
c        festgelegt!!                                                             |   
c---------------------------------------------------------------------+
c     History                                                         |
c     element from RLT see Mafelap conference                         |
c     WW elastic foundation added                                     | 
c     WW error analysis implemented (only bending terms) 8.9.97       |
c     WW concrete design added                                        | 
c     WW plot stresse at top and bottom                               |
c     WW QLOA added                                                   | 
c     WW choose base system (def from diagonals)                      | 
c---------------------------------------------------------------------+
      USE bdata
      USE cdata
      USE cdat1
      USE eldata
      USE errin2
      USE errin3
      USE iofile
      USE pdata10
      USE plslay
      USE pltran
      USE prisdat
      USE prlod
      USE qload
      USE strnam
      implicit double precision (a-h,o-z)
      character comp*25
      dimension xl(ndm,*),tl(*),yl(3,4),sigi(6),sigb(6),ix(*),
     1        d(*),ul(ndf,*),s(nst,*),p(nst),sg(4),tg(4),wg(4),
     1        btd(3,6),gshp1(3,4),gshp2(3,4),dvl(4),eps(6),
     2        vl(6,4),tr(3,3),ev(8),bl(3),as(4),ql(10)
      dimension h1(*),h2(*),h3(*)
      common /shld09/ bm(3,6,4),ii1,ii2
!$OMP THREADPRIVATE (/shld09/)  
      common /shpf09/ shp(3,8,4),shp1(3,4,4),shp2(3,4,4)
!$OMP THREADPRIVATE (/shpf09/)  
c.... transfer to correct processor
      go to (1,2,3,3,3,3,2,3,3,2,2,2,13,3,2,2,2,2,2,2,2,3), isw
      return
c.... input the material properties
c.... record 1+2
1     if(ior.lt.0) write(*,2004)
2004  format(' Specify material parameters for shell:'/
     1  '  a.) E, nu, rho, thick '/
     2  '  b.) 1-grav, 2-grav, pressure, x-gravity, y-gravity,',
     3  ' z-gravity, load-type, el. found. const.,z')
      call dinput(ev,5)
      call dinput(d(7),9)
      d(16) = d(13) ! tc
c.... elast. foundation always/only for compression = 0/1
      d(18) = d(14)
      if(d(18).ge.0.d0)then
         comp='(tension/compression)'
      else
         comp='(  only  compression)'
      end if   
      d(19) = d(15) ! pz
c
      d(1) = ev(1)/(1.-ev(2)*ev(2))*ev(4)
      d(2) = ev(2)*d(1)
      d(3) = (d(1)-d(2))/2.0d0
      xx1  = ev(4)*ev(4)/12.
      d(15)= xx1
      d(20)= ev(4)
cww      d(4) = d(1)*xx1  used for t_1
cww      d(5) = d(2)*xx1
cww      d(6) = d(3)*xx1
      d(13) = ev(3)*ev(4)
      d(14) = d(13)*xx1
      write(iow,2000) (ev(i),i=1,4),(d(i),i=7,12)
      write(iow,2005) d(16)
      write(iow,2006) d(18),comp
      write(iow,2007) d(19)
      if(ior.lt.0) then
          write(*,2000) (ev(i),i=1,5),(d(i),i=7,12)
          write(*,2005) d(16)
          write(*,2006) d(18),comp
          write(*,2007) d(19)
      end if
c.... record 3
      call dinput(d(21),4)
      d(17)=d(21)
      if(d(17).eq.0.d0) d(17)=1.d0     
      d( 4)=d(22)   
      d( 5)=d(23)   
      d( 6)=d(24)   
      d(21)=0.d0 
      d(22)=0.d0 
      d(23)=0.d0 
      d(24)=0.d0 
                   write(iow,2008) (d(17),d(i),i=4,6)
      if(ior.lt.0) write(*  ,2008) (d(17),d(i),i=4,6)
c
c.... record 4
      call dinput(d(21),4)
                   write(iow,2012) (d(i),i=21,24)
      if(ior.lt.0) write(*  ,2012) (d(i),i=21,24)
c
c...  names for principal moments
      nptyp = 5 
c...  position for principal moments
      nprip(4)=6
      nprip(5)=7
      nprip(6)=8
      return
c.... check element for errors
2     return
c.... compute the element tangent array
3      if(n.eq.1) then
c.... description of stresses  
      if(klay.eq.0) then
        strsus( 1) = ' N-FORCE  n_11 '
        strsus( 2) = ' N-FORCE  n_12 '
        strsus( 3) = ' N-FORCE  n_22 '
        strsus( 4) = ' N-FORCE  n_1  '
        strsus( 5) = ' N-FORCE  n_2  '
        strsus( 6) = ' MOMENT   m_11 '
        strsus( 7) = ' MOMENT   m_12 '
        strsus( 8) = ' MOMENT   m_22 '
        strsus( 9) = ' MOMENT   m_1  '
        strsus(10) = ' MOMENT   m_2  '
        strsus(11) = ' Pressure  p   '
        if(d(21).ne.0.d0) then  ! concrete design
          strsus( 4) = 'STEEL as_x(Bot)'
          strsus( 5) = 'STEEL as_x(Top)'                
          strsus( 9) = 'STEEL as_y(Bot)'                
          strsus(10) = 'STEEL as_y(Top)'                
        end if                          
      else if(klay.ne.0) then
        strsus( 1) = ' STRESS S_11(T)'
        strsus( 2) = ' STRESS S_12(T)'
        strsus( 3) = ' STRESS S_22(T)'
        strsus( 4) = ' STRESS S_1 (T)'
        strsus( 5) = ' STRESS S_2 (T)'
        strsus( 6) = ' STRESS S_11(B)'
        strsus( 7) = ' STRESS S_12(B)'
        strsus( 8) = ' STRESS S_22(B)'
        strsus( 9) = ' STRESS S_1 (B)'
        strsus(10) = ' STRESS S_2 (B)'
        strsus(11) = ' Pressure  p   '
      end if
      do is =12,25
        strsus(is) = '               '
      end do
      end if
c 
      l = 2
      call pgauss(l ,lint,sg,tg,wg)
c.... compute the transformation and midsurface coords
      call tran09(d,xl,yl,tr,ndm)
      call jacq09(yl)
      if(yl(3,1).ne.0.0) then
        ii1 = 1
        ii2 = 6
      else
        ii1 = 3
        ii2 = 5
      end if
c.... construct the integrals of the drilling shape functions
      call pzero(gshp1,12)
      call pzero(gshp2,12)
      dv = 0.0
      do 302 l = 1,lint
        call rshp09(sg(l),tg(l),yl,shp(1,1,l),shp1(1,1,l),
     1  shp2(1,1,l),xsj,3)
        dvl(l) = xsj*wg(l)
        dv     = dv + dvl(l)
        do 300 i = 1,3
        do 300 j = 1,4
          gshp1(i,j) = gshp1(i,j) + shp1(i,j,l)*dvl(l)
          gshp2(i,j) = gshp2(i,j) + shp2(i,j,l)*dvl(l)
300     continue
302   continue
      if(isw.eq. 4) go to  4
      if(isw.eq.14) go to 14
c.... construct the modified drilling shape functions
      do 306 i = 1,3
      do 306 j = 1,4
        dv1 = gshp1(i,j)/dv
        dv2 = gshp2(i,j)/dv
        do 304 l = 1,lint
          shp1(i,j,l) = shp1(i,j,l) - dv1
          shp2(i,j,l) = shp2(i,j,l) - dv2
304     continue
306   continue
      if(isw.eq.5) go to 5
      if(isw.eq.8) go to 8
      if(isw.eq.9) go to 9

c.... calculate and transform the element load vector
      call qload09(ql,d,aqloa,numel,n,mqloa,propq,prop,isw)
      do 308 i = 1,3
        bl(i) = ql(i)
        do 308 j = 1,3
          bl(i) = bl(i) + tr(i,j)*ql(j+3)
308   continue
      tc= d(16)
c.... compute membrane/load parts
      do 350 l = 1,lint
c.... membrane and bending stiffness part
        call dktq09(shp(1,5,l),shp(1,1,l))
        dv = dvl(l)
        dv1 = d(1)*dv
        dv2 = d(2)*dv
        dv3 = d(3)*dv
        dv4 = d(1)*dv*d(15)
        dv5 = d(2)*d(15)*dv
        dv6 = d(3)*d(15)*dv
c....   calculate pressure if elastic foundation
        cb = d(18)
        press = 0.d0
        do ip = 1,4
	        press =  press + dabs(cb)*shp(3,ip,l)*ul(3,ip)
        end do     
c

c....   z-position 
        zz=0.d0 
        do ip = 1,nel
          zz = zz + shp(3,ip,l)*xl(3,ip)
        end do
cww        zz = 20-zz  ! ####### NUR FÜR BSP CTWM WW ##########  

        i1 = 1
        do 340 i = 1,4
c.... recover the previously computed shape functions
          shp1i = shp(1,i,l)
          shp2i = shp(2,i,l)
          shp3i = shp(3,i,l)
          shp11 = shp1(1,i,l)
          shp12 = shp1(2,i,l)
          shp13 = shp1(3,i,l)
          shp21 = shp2(1,i,l)
          shp22 = shp2(2,i,l)
          shp23 = shp2(3,i,l)
          i2 = i1 - 1

c.... compute the loading term
c         multiply with z
          pz=1.d0
          if(d(19).ne.0.d0) pz=zz  

          p(i1  ) = p(i1  ) + shp3i*bl(1)*dv*pz
          p(i1+1) = p(i1+1) + shp3i*bl(2)*dv*pz
          p(i1+2) = p(i1+2) + shp3i*bl(3)*dv*pz
          p(i1+3) = p(i1+3) + shp13*bl(3)*dv*tc*pz
          p(i1+4) = p(i1+4) + shp23*bl(3)*dv*tc*pz
          p(i1+5) = p(i1+5)-(shp13*bl(1)*pz + shp23*bl(2)*pz)*dv*tc
          if(isw.eq.22) go to 340  
c.... form the stress-displacement matrix (Bi-trans * D)
          a11 =  dv1*shp1i
          a12 =  dv2*shp2i
          a13 = -dv1*shp11 - dv2*shp22
          a21 =  dv2*shp1i
          a22 =  dv1*shp2i
          a23 = -dv2*shp11 - dv1*shp22
          a31 =  dv3*shp2i
          a32 =  dv3*shp1i
          a33 = -dv3*(shp12+shp21)
c.... form the plate stress-displacement matrix
          do 310 ii = ii1,ii2
            btd(1,ii) = dv4*bm(1,ii,i) + dv5*bm(2,ii,i)
            btd(2,ii) = dv5*bm(1,ii,i) + dv4*bm(2,ii,i)
            btd(3,ii) = dv6*bm(3,ii,i)
310       continue
c.... loop on columns
          j1 = i1
          do 330 j = i,4
            j2 = j1 - 1
            xn = shp(1,j,l)
            yn = shp(2,j,l)
            x1n = - shp1(1,j,l)
            y2n = - shp2(2,j,l)
            xyn = - shp1(2,j,l) - shp2(1,j,l)
c.... compute the membrane part
            s(i1  ,j1  ) = s(i1  ,j1  ) + (a11*xn + a31*yn)
            s(i1+1,j1  ) = s(i1+1,j1  ) + (a12*xn + a32*yn)
            s(i1  ,j1+1) = s(i1  ,j1+1) + (a21*yn + a31*xn)
            s(i1+1,j1+1) = s(i1+1,j1+1) + (a22*yn + a32*xn)
            s(i1+5,j1  ) = s(i1+5,j1  ) + a13*xn  + a33*yn
            s(i1+5,j1+1) = s(i1+5,j1+1) + a23*yn  + a33*xn
            s(i1  ,j1+5) = s(i1  ,j1+5) + a11*x1n + a21*y2n + a31*xyn
            s(i1+1,j1+5) = s(i1+1,j1+5) + a12*x1n + a22*y2n + a32*xyn
            s(i1+5,j1+5) = s(i1+5,j1+5) + a13*x1n + a23*y2n + a33*xyn
c.... compute the bending part
            do 320 ii = ii1,ii2
            do 320 jj = ii1,ii2
            do 320 kk = 1,3
              s(i2+ii,j2+jj) = s(i2+ii,j2+jj) + btd(kk,ii)*bm(kk,jj,j)
320         continue
c.... compute elastic foundation term, add term always/only for negative values
            if(cb.gt.0.or.(cb.lt.0.and.press.lt.0.0d0)) then
              s(i2+3,j2+3) = s(i2+3,j2+3)+dabs(cb)*shp3i*shp(3,j,l)*dv
            end if
330       j1 = j1 + ndf
340     i1 = i1 + ndf
350   continue
      i1 = 1
      if(yl(3,1).ne.0.0) then
        do 370 i = 1,4
          p(i1+3) = p(i1+3) + yl(3,i)*p(i1+1)
          p(i1+4) = p(i1+4) - yl(3,i)*p(i1)
          j1 = i1
          do 360 j = i,4
            call proj09(s(i1,j1),p(i1),yl(3,i),yl(3,j),nst)
            j1 = j1 + ndf
360       continue
          i1 = i1 + ndf
370     continue
      end if
      call rots09(s,p,tr,nst,ndf)
c.... form residual
      do 380 i = 1,nst
      do 380 j = 1,nst
         p(i) = p(i) - s(i,j)*ul(j,1)
380   continue
c
c.... calculate eigenvalues
cww   call elemev(s,nel,ndf,nst)
c
      return
c.... compute and output the element variables
4     continue
      do 400 i = 1,4
      do 400 j = 1,3
         vl(j  ,i) = 0.0
         vl(j+3,i) = 0.0
         do 390 k = 1,3
            vl(j  ,i) = vl(j  ,i) + tr(j,k)*ul(k  ,i)
            vl(j+3,i) = vl(j+3,i) + tr(j,k)*ul(k+3,i)
390      continue
400   continue
      l = 1
      call pgauss(l,lint,sg,tg,wg)
      do 440 l = 1,lint
        call rshp09(sg(l),tg(l),yl,shp,shp1,shp2,xsj,3)
c.... modify the rotational shape functions
        do 405 i = 1,3
        do 405 j = 1,4
            shp1(i,j,1) = shp1(i,j,1) - gshp1(i,j)/dv
            shp2(i,j,1) = shp2(i,j,1) - gshp2(i,j)/dv
405     continue
        call stre09(d,xl,vl,ndm,nel,1, xx,yy,zz,eps,sigi,sigb,press)
c....   modify to stresses
        if(klay.ne.0) call stre09s(d,sigi,sigb,eps)
        mct = mct - 3
        if(mct.gt.0) go to 430
        if(klay.eq.0) then
                       write(iow,2001) o,head
          if(ior.lt.0) write(*  ,2001) o,head
        elseif(klay.ne.0) then
                       write(iow,2011) o,head
          if(ior.lt.0) write(*  ,2011) o,head
        end if
        mct = 50
430                write(iow,2002) n,xx,sigi,press,ma,yy,sigb,zz,eps
      if(ior.lt.0) write(*  ,2002) n,xx,sigi,press,ma,yy,sigb,zz,eps
440   continue
      return
c.... compute element mass arrays (lumped mass matrix)
5     l = 2
      call pgauss(l,lint,sg,tg,wg)
      do 510 l = 1,lint
        call dktq09(shp(1,5,l),shp(1,1,l))
        dv = dvl(l)
        i1 = 1
        do 500 i = 1,4
c.... recover the previously computed shape functions
          shp1i = shp(1,i,l)
          shp2i = shp(2,i,l)
          shp3i = shp(3,i,l)
          shp12 = shp1(2,i,l)
          shp13 = shp1(3,i,l)
          shp23 = shp2(3,i,l)
c.... compute lumped mass term
          p(i1  ) = p(i1  ) + shp3i*d(13)*dv
          p(i1+1) = p(i1+1) + shp3i*d(13)*dv
          p(i1+2) = p(i1+2) + shp3i*d(13)*dv
          p(i1+3) = p(i1+3) + shp13*d(14)*dv
          p(i1+4) = p(i1+4) + shp23*d(14)*dv
          p(i1+5) = p(i1+5) + (shp13*d(13) + shp23*d(13))*dv
c.... compute shape functions
c        call shape(sg(l),tg(l),xl,shpp,xsj,2,nel,ix,.true.)
c        w11 = wg(l)*xsj*d(13)
c        w12 = wg(l)*xsj*d(14)
cc.... for each node j compute db = rho*shape*xsj
c        j1 = 1
c        do 500 j = 1,nel
c          w1x = w11 * shpp(3,j)
c          w2x = w12 * shpp(3,j)
c          p(j1  ) = p(j1) + w1x
c          p(j1+1) = p(j1)
c          p(j1+2) = p(j1)
c          p(j1+3) = p(j1+3) + w2x
c          p(j1+4) = p(j1+3)
          i1 = i1 + ndf
500     continue
510   continue
      return
c
c.... project stresses 
8     istv = -10
c
      if(iplma(ma).eq.0)  return ! only if MATN
c
c.....xl necessary  not yl 
      call stcn09(ix,d,xl,ul,tr,strea,strea(1+numnp),dvl,ndf,nel,numnp,
     1     klay)
      return
c                
c.... calculate errors
9     igeo = d(17)
      if(igeo.ne.2) stop 't_q must be defined in case of adaptivity!'
      call stcn09e(ix,d,xl,ul,tr,strea,strea(1+numnp),dvl,ndf,nel,
     +             numnp,numel,klay,ndm,e_ome)
      return
c
c.... plot local base system at center of element via FORC,1
13    continue
c     in case of matn
      if(iplma(ma).eq.0) return
      call tran09(d,xl,yl,tr,ndm)
      call pppcol(2)
      call plloco09(tr,xl)
      return
c
c.... calculate stresses sig(i) at center of element
14    continue
      do i = 1,4
        do j = 1,3
          vl(j  ,i) = 0.0
          vl(j+3,i) = 0.0
          do k = 1,3
             vl(j  ,i) = vl(j  ,i) + tr(j,k)*ul(k  ,i)
             vl(j+3,i) = vl(j+3,i) + tr(j,k)*ul(k+3,i)
          end do
        end do  
      end do  
      l = 1
      call pgauss(l,lint,sg,tg,wg)
      call rshp09(sg(l),tg(l),yl,shp,shp1,shp2,xsj,3)
c.... modify the rotational shape functions
      do i = 1,3
        do j = 1,4
          shp1(i,j,1) = shp1(i,j,1) - gshp1(i,j)/dv
          shp2(i,j,1) = shp2(i,j,1) - gshp2(i,j)/dv
        end do
      end do
      call stre09(d,xl,vl,ndm,nel,1, xx,yy,zz,eps,sigi,sigb,press)
c.... modify to stresses
      if(klay.ne.0) call stre09s(d,sigi,sigb,eps)
c.... concrete design
      if(d(21).gt.0.0d0.and.klay.eq.0) call design09(d,sigi,sigb,as)
c
      if(nfp.gt.11.or. nfp.lt. 1) return
      if(nfp.ge. 1.and.nfp.le. 5) strp=sigi(nfp)
      if(nfp.ge. 6.and.nfp.le.10) strp=sigb(nfp-5)
      if(nfp.eq.11)               strp=press
      if(klay.eq.0.and.d(21).ne.0) then ! concrete instead of mean stresses
        if(nfp.ge. 4.and.nfp.le. 5) strp=as(nfp-3)
        if(nfp.ge. 9.and.nfp.le.10) strp=as(nfp-6)
      end if
c
      if(flfp) then
c....   calculate extreme values
        xmaxf = max(xmaxf,strp)
        xminf = min(xminf,strp)
      else
c....   plot stresses 
c....   color      
        call pppcolf(strp) 
c.....  transform ccordinates for rot        
        call plxtrn(xl,tra,vr,3,4)
c...... plot element
        call plot9(iel,ix,xl,ndm,nel,1)
      end if
      return
c
c.... formats
2000  format(/5x,'3-D  S h e l l   E l e m e n t '//
     1 10x,'modulus',e18.5/10x,'poisson ratio',f8.5/
     2 10x,'density',e18.5/10x,'thickness',e16.5/
     4 10x,'1-gravity',e16.5/10x,'2-gravity',e16.5/
     5 10x,'pressure ',e16.5/10x,'x-gravity',e16.5/
     6 10x,'y-gravity',e16.5/10x,'z-gravity',e16.5/1x)
2001  format(a1,20a4//5x,'S h e l l  S t r e s s  R e s u l t a n t s '
     1//' elmt x-coord  N_xx       N_xy       N_yy     ',
     2        '   N_1        N_2       angle_N     Pressure '/
     3  ' matl y-coord  M_xx       M_xy       M_yy     ',
     4        '   M_1        M_2       angle_M'/
     5  '      z-coord  eps_xx     eps_xy     eps_yy   ',
     6        '  kappa_xx   kappa_xy   kappa_yy '/
     7    38(' -'))
2002  format(i5,0p1f8.3,1p5e11.3,0p1f8.2,3x,e11.3/i5,0p1f8.3,1p5e11.3,
     1       0p1f8.2/5x,0p1f8.3,1p6e11.3/1x)
c2003  format(' error * * reinput last data line')
2005  format(10x,' Load constant (0=lump, 1=cons)', f5.1)
2006  format(10x,' Elastic foundation constant   ', e12.5,a)
2007  format(10x,' pz  = 0:p=p  1=p=p*z          ', e12.5)
2008  format(10x,' base vector (1-17) see BASE   ', f5.1,/,
     1       10x,' Coordinates of base vector t_q', 3f12.5)
2011  format(a1,20a4//5x,'S h e l l  S t r e s s e s'
     1//' elmt x-coord  S_xx(t)    S_xy(t)    S_yy(t)  ',
     2        '   S_1(t)     S_2(t)    angle(t)'/
     3  ' matl y-coord  S_xx(b)    S_xy(b)    S_yy(b)  ',
     4        '   S_1(b)     S_2(b)    angle(b)'/
     5  '      z-coord  eps_xx(t)  eps_xy(t)  eps_yy(t)',
     6        '  eps_xx(b)  eps_xy(b)  eps_yy(b)'/
     6    38(' -'))
2012  format(10x,'Values for Concrete, Design due to DIN 1045 '/
     + 10x,'beta_R (B25: 17.5 e3 kN/m**2)       =',e12.5/
     + 10x,'beta_S (Bst: 500/550 500 e3 kN/m**2)=',e12.5/
     + 10x,'Position h_1 ...................... =',e12.5/
     + 10x,'Position h_2 ...................... =',e12.5)
      end
c
      subroutine stre09(d,xl,vl,ndm,nel,l, xx,yy,zz,eps,sigi,sigb,press)
      implicit double precision (a-h,o-z)
      dimension   xl(ndm,*),sigi(6),sigb(6),d(*),vl(6,*), eps(6)
      common /shld09/ bm(3,6,4),ii1,ii2
!$OMP THREADPRIVATE (/shld09/)  
      common /shpf09/ shp(3,8,4),shp1(3,4,4),shp2(3,4,4)
!$OMP THREADPRIVATE (/shpf09/)  
c.... compute the membrane and bending strains
      call dktq09(shp(1,5,l),shp(1,1,l))
      do 410 i = 1,6
        eps(i) = 0.0
410   continue
      xx = 0.0
      yy = 0.0
      zz = 0.0
      press = 0.0
      do 420 j = 1,nel
        xx = xx + shp(3,j,l)*xl(1,j)
        yy = yy + shp(3,j,l)*xl(2,j)
        zz = zz + shp(3,j,l)*xl(3,j)
        eps(1) = eps(1) + shp(1,j,l)*vl(1,j)
     1                  - shp1(1,j,l)*vl(6,j)
        eps(3) = eps(3) + shp(2,j,l)*vl(2,j)
     1                  - shp2(2,j,l)*vl(6,j)
        eps(2) = eps(2) + shp(1,j,l)*vl(2,j) + shp(2,j,l)*vl(1,j)
     1                  - (shp1(2,j,l) + shp2(1,j,l))*vl(6,j)
        do 415 i = ii1,ii2
          eps(4) = eps(4) + bm(1,i,j)*vl(i,j)
          eps(6) = eps(6) + bm(2,i,j)*vl(i,j)
          eps(5) = eps(5) + bm(3,i,j)*vl(i,j)
415     continue
        cb = d(18)
        press = press + dabs(cb)*shp(3,j,l)*vl(3,j)
        if(cb.lt.0.d0.and.press.ge.0.d0) press = 0.d0 
420   continue
      sigi(1) = d(1)*eps(1) + d(2)*eps(3)
      sigi(3) = d(1)*eps(3) + d(2)*eps(1)
      sigi(2) = d(3)*eps(2)
      call pstres(sigi,sigi(4),sigi(5),sigi(6))
      sigb(1) = d(1)*d(15)*eps(4) + d(2)*d(15)*eps(6)   
      sigb(3) = d(2)*d(15)*eps(4) + d(1)*d(15)*eps(6)
      sigb(2) = d(3)*d(15)*eps(5)
      call pstres(sigb,sigb(4),sigb(5),sigb(6))
      return
      end
c
      subroutine jacq09(xl)
      implicit double precision (a-h,o-z)
      dimension xl(3,*)
      common /shld09/ bm(3,6,4),ii1,ii2
!$OMP THREADPRIVATE (/shld09/)  
      common /elcom091/ b(4),c(4),aa(4),bb(4),cc(4),dd(4),ee(4)
!$OMP THREADPRIVATE (/elcom091/)  
c.... form geometric constants for DKQ element
      cxx = 0.
      cyy = 0.
      cxy = 0.
      do 100 i = 1,4
      k = mod(i,4) + 1
      b(i) = xl(2,k) - xl(2,i)
      c(i) = xl(1,i) - xl(1,k)
      dz   = xl(3,k) - xl(3,i)
      b2 = b(i)*b(i)
      c2 = c(i)*c(i)
      cxx = cxx + dz*b2/(b2 + c2)
      cyy = cyy + dz*c2/(b2 + c2)
      cxy = cxy + (dz+dz)*b(i)*c(i)/(b2 + c2)
      sql = (b2 + c2)/0.75
      aa(i) = (c(i)+c(i))/sql
      bb(i) = b(i)*c(i)/sql
      cc(i) = c2/sql
      dd(i) =-(b(i)+b(i))/sql
      ee(i) = b2/sql
      b(i)  = b(i)/8.
100   c(i)  = c(i)/8.
      if(xl(3,1).ne.0.0) then
        x31 = xl(1,3) - xl(1,1)
        y31 = xl(2,3) - xl(2,1)
        x42 = xl(1,4) - xl(1,2)
        y42 = xl(2,4) - xl(2,2)
        a4  = (x31*y42 - x42*y31)*2.0
        ad  = a4*a4/4.0
        x31 = x31/ad
        y31 = y31/ad
        x42 = x42/ad
        y42 = y42/ad
        bm(1,1,1) = -cxx*x42
        bm(2,1,1) = -cyy*x42
        bm(3,1,1) = -cxy*x42
        bm(1,2,1) = -cxx*y42
        bm(2,2,1) = -cyy*y42
        bm(3,2,1) = -cxy*y42
        bm(1,6,1) = -cxx/a4
        bm(2,6,1) = -cyy/a4
        bm(3,6,1) = -cxy/a4
        bm(1,1,2) =  cxx*x31
        bm(2,1,2) =  cyy*x31
        bm(3,1,2) =  cxy*x31
        bm(1,2,2) =  cxx*y31
        bm(2,2,2) =  cyy*y31
        bm(3,2,2) =  cxy*y31
        bm(1,6,2) = -cxx/a4
        bm(2,6,2) = -cyy/a4
        bm(3,6,2) = -cxy/a4
        do 210 i = 1,3
          bm(i,1,3) = -bm(i,1,1)
          bm(i,2,3) = -bm(i,2,1)
          bm(i,6,3) =  bm(i,6,1)
          bm(i,1,4) = -bm(i,1,2)
          bm(i,2,4) = -bm(i,2,2)
          bm(i,6,4) =  bm(i,6,2)
210     continue
      end if
      return
      end
c
      subroutine dktq09(shm,shn)
      implicit double precision (a-h,o-z)
      dimension shm(3,4),shn(3,4)
      common /shld09/bm(3,6,4),ii1,ii2
!$OMP THREADPRIVATE (/shld09/)  
      common /elcom091/b(4),c(4),aa(4),bb(4),cc(4),dd(4),ee(4)
!$OMP THREADPRIVATE (/elcom091/)  
c.... form strain-displacement array for DKQ element
      do 110 i = 1,4
      j = mod(i+2,4) + 1
      bm(1,3,i) = aa(i)*shm(1,i) - aa(j)*shm(1,j)
      bm(1,4,i) = bb(i)*shm(1,i) + bb(j)*shm(1,j)
      bm(1,5,i) = cc(i)*shm(1,i) + cc(j)*shm(1,j) - shn(1,i)
      bm(2,3,i) = dd(i)*shm(2,i) - dd(j)*shm(2,j)
      bm(2,4,i) =-ee(i)*shm(2,i) - ee(j)*shm(2,j) + shn(2,i)
      bm(2,5,i) =-bb(i)*shm(2,i) - bb(j)*shm(2,j)
      bm(3,3,i) = aa(i)*shm(2,i) - aa(j)*shm(2,j)
     1          + dd(i)*shm(1,i) - dd(j)*shm(1,j)
      bm(3,4,i) =-ee(i)*shm(1,i) - ee(j)*shm(1,j) + shn(1,i)
     1          - bm(2,5,i)
      bm(3,5,i) = cc(i)*shm(2,i) + cc(j)*shm(2,j) - shn(2,i)
     1          - bm(1,4,i)
110   continue
      return
      end
c
      subroutine tran09(d,xl,yl,t,ndm)
      implicit double precision (a-h,o-z)
c
c.... compute the transformation array and surface coords.
c
      dimension   d(*),x0(3),xl(ndm,*),yl(3,4),t(3,3),tt(3,3),add(3)
c

      igeo = d(17) 

      add(1)=d(4)
      add(2)=d(5)
      add(3)=d(6)

      if(igeo.le.2) then ! diagonals
c....   compute the inplane direction cosines (bisect diagonals)
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
        do i = 1,3
          t(1,i) = t(1,i)/dl1
          t(2,i) = t(2,i)/dl2
c....     compute the center (0,0) displacement
          x0(i) = 0.25*(xl(i,1) + xl(i,2) + xl(i,3) + xl(i,4))
        end do
c....   compute the normal to the surface
        t(3,1) = t(1,2)*t(2,3) - t(2,2)*t(1,3)
        t(3,2) = t(1,3)*t(2,1) - t(2,3)*t(1,1)
        t(3,3) = t(1,1)*t(2,2) - t(2,1)*t(1,2)
        if(igeo.eq.2) then ! given vector t_1
c....     modify base vector t_1, input must not be a unit vector  
          td = dot(d(4),d(4),3)
          if(td.gt.0.d0) then 
            td=sqrt(td)
            t(1,1) = d(4)/td
            t(1,2) = d(5)/td
            t(1,3) = d(6)/td
c....       calculate t_2 = t_3 x t_1
            t(2,1) = t(3,2)*t(1,3)-t(3,3)*t(1,2)
            t(2,2) =-t(3,1)*t(1,3)+t(3,3)*t(1,1)
            t(2,3) = t(3,1)*t(1,2)-t(3,2)*t(1,1)
c....       calculate t_1 = t_2 x t_3 new
            t(1,1) = t(2,2)*t(3,3)-t(2,3)*t(3,2)
            t(1,2) =-t(2,1)*t(3,3)+t(2,3)*t(3,1)
            t(1,3) = t(2,1)*t(3,2)-t(2,2)*t(3,1)
          end if
        end if
      else if(igeo.gt.2) then ! analytical solutions         
c....   compute the center (0,0)
        do i = 1,3
          x0(i) = 0.25*(xl(i,1) + xl(i,2) + xl(i,3) + xl(i,4))
        end do
        call pdirec3(tt,add,x0,igeo) ! tt=t^T
c....   transpose tt
        do i=1,3
          do k=1,3 
            t(i,k)=tt(k,i)
          end do
        end do
      end if

c.... compute the projected middle surface coordinates
      do 140 i = 1,4
      do 140 j = 1,3
        yl(j,i) = 0.0
        do 130 k = 1,3
          yl(j,i) = yl(j,i) + t(j,k)*(xl(k,i) - x0(k))
130     continue
140   continue
c.... set offset coordinates to zero if small compared to plan size
      htol =  0.0
      do 150 i = 1,4
        htol = max(htol,abs(yl(1,i)),abs(yl(2,i)))
150   continue
      htol = htol*1.e-7
      do 160 i = 1,4
        if(abs(yl(3,i)) .le. htol) yl(3,i) = 0.0
160   continue
      return
      end
c
      subroutine rshp09(ss,tt,x,shp,shp1,shp2,xsj,ndm)
      implicit double precision (a-h,o-z)
c
c.... shape function routine for two dimensional elements
c
      dimension   x(ndm,*),s(4),t(4),
     1       shp(3,8),sx(2,2),shp1(3,4),shp2(3,4)
      data s/-0.5,0.5,0.5,-0.5/,t/-0.5,-0.5,0.5,0.5/
c.... form 4-node quadrilateral shape functions
      shp(1,2) = 0.25*(1.-tt)
      shp(1,3) = 0.25*(1.+tt)
      shp(1,1) = -shp(1,2)
      shp(1,4) = -shp(1,3)
      shp(2,4) = 0.25*(1.-ss)
      shp(2,3) = 0.25*(1.+ss)
      shp(2,2) = -shp(2,3)
      shp(2,1) = -shp(2,4)
      shp(3,1) = shp(1,2)*(1.-ss)
      shp(3,2) = shp(1,2)*(1.+ss)
      shp(3,3) = shp(1,3)*(1.+ss)
      shp(3,4) = shp(1,3)*(1.-ss)
      s2 = (1.-ss*ss)*0.5
      t2 = (1.-tt*tt)*0.5
      shp(3,5) =  s2*(1-tt)
      shp(3,6) =  t2*(1+ss)
      shp(3,7) =  s2*(1+tt)
      shp(3,8) =  t2*(1-ss)
      shp(1,5) = -ss*(1-tt)
      shp(1,6) =  t2
      shp(1,7) = -ss*(1+tt)
      shp(1,8) = -t2
      shp(2,5) = -s2
      shp(2,6) = -tt*(1+ss)
      shp(2,7) =  s2
      shp(2,8) = -tt*(1-ss)
c.... construct jacobian and its inverse
      sx(2,2) =  x(1,1)*shp(1,1)+x(1,2)*shp(1,2)+x(1,3)*shp(1,3)
     1          +x(1,4)*shp(1,4)
      sx(1,2) =  x(1,1)*shp(2,1)+x(1,2)*shp(2,2)+x(1,3)*shp(2,3)
     1          +x(1,4)*shp(2,4)
      sx(2,1) =  x(2,1)*shp(1,1)+x(2,2)*shp(1,2)+x(2,3)*shp(1,3)
     1          +x(2,4)*shp(1,4)
      sx(1,1) =  x(2,1)*shp(2,1)+x(2,2)*shp(2,2)+x(2,3)*shp(2,3)
     1          +x(2,4)*shp(2,4)
      xsj = sx(1,1)*sx(2,2)-sx(1,2)*sx(2,1)
      xsj1 = xsj
      if(xsj.eq.0.0d0) xsj1 = 1.0
      sx(2,2) = sx(2,2)/xsj1
      sx(1,1) = sx(1,1)/xsj1
      sx(1,2) =-sx(1,2)/xsj1
      sx(2,1) =-sx(2,1)/xsj1
c.... form global derivatives
      do 140 i = 1,8
      tp        = shp(1,i)*sx(1,1)+shp(2,i)*sx(2,1)
      shp(2,i)  = shp(1,i)*sx(1,2)+shp(2,i)*sx(2,2)
140   shp(1,i) = tp
c.... form the rotational and 5-th shape functions
      call hshp09(shp(1,5),shp1,shp2)
      return
      end
c
      subroutine hshp09(shp,shp1,shp2)
      implicit double precision (a-h,o-z)
c..... dimension arrays
      dimension shp(3,4),shp1(3,4),shp2(3,4),shx(4),shy(4)
      common /elcom091/b(4),c(4),aa(4),bb(4),cc(4),dd(4),ee(4)
!$OMP THREADPRIVATE (/elcom091/)  
c....  form the shape functions
      do 130 l = 1,3
        do 110 k = 1,4
          shx(k) = shp(l,k)*c(k)
          shy(k) = shp(l,k)*b(k)
110     continue
        j = 4
        do 120 i = 1,4
          shp1(l,i) = shy(i) - shy(j)
          shp2(l,i) = shx(i) - shx(j)
          j = i
120     continue
130   continue
      return
      end
c
      subroutine rots09(s,p,t,nst,ndf)
      implicit double precision (a-h,o-z)
c
c.... transform the loads and stiffness to global coords.
c
      dimension s(nst,nst),p(nst),t(3,3),a(3,3),b(6)
      i0 = 0
      do 170 ir = 1,4
        do 110 ii = 1,3
          b(ii  ) = dot(t(1,ii),p(i0+1),3)
          b(ii+3) = dot(t(1,ii),p(i0+4),3)
110     continue
        do 111 ii = 1,6
          p(i0+ii) = b(ii)
111     continue
        j0 = i0
        do 160 jc = ir,4
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
        if(ir.ne.jc) then
          do 155 i = 1,6
          do 155 j = 1,6
            s(j0+j,i0+i) = s(i0+i,j0+j)
155       continue
        end if
160     j0 = j0 + ndf
170   i0 = i0 + ndf
      return
      end
c
      subroutine stcn09(ix,d,yl,ul,tr,dt,st,dvl,ndf,nel,numnp,klay)
      implicit double precision (a-h,o-z)
      dimension   dt(numnp),st(numnp,*),yl(3,*),sigi(6),sigb(6),
     1       tr(3,3),d(*),ul(ndf,*),vl(6,4),eps(6),dvl(4),ix(*),as(4)
      common /shpf09/ shp(3,8,4),shp1(3,4,4),shp2(3,4,4)
!$OMP THREADPRIVATE (/shpf09/)  
c.... compute membrane and bending stresses for the projection.
      do 400 i = 1,4
        do 400 j = 1,3
            vl(j  ,i) = 0.0
            vl(j+3,i) = 0.0
            do 399 k = 1,3
              vl(j  ,i) = vl(j  ,i) + tr(j,k)*ul(k  ,i)
              vl(j+3,i) = vl(j+3,i) + tr(j,k)*ul(k+3,i)
399         continue
400   continue
      do 440 l = 1,4   ! Gauss-Points
        call stre09(d,yl,vl,3,nel,l, xx,yy,zz,eps,sigi,sigb,press)
c....   sum stresses for error analysis
        call elener09(d,sigi,sigb,dvl(l))
c....   modify to stresses
        if(klay.ne.0) call stre09s(d,sigi,sigb,eps)
c....   concrete design
        if(d(21).gt.0.0d0.and.klay.eq.0) then
          call design09(d,sigi,sigb,as)
        end if
c....   store stresses for the projection.
        do 435 j = 1,4   ! Nodes
          xsji = dvl(l)*shp(3,j,l)
          ii = iabs(ix(j))
          if(ii.le.0) go to 435
          dt(ii) = dt(ii) + xsji
c.......  N_x(1),N_xy(2),N_y(3),N_1(4),N_2(5)
c.......  M_x(6),M_xy(7),M_y(8),M_1(9),M_2(10),press
c         or 
c.......  N_x(1),N_xy(2),N_y(3),as1(B),as1(T)
c.......  M_x(6),M_xy(7),M_y(8),as2(B),as2(T),press
c         or
c.......  S_x(t),S_xy(t),S_y(t),S_1(t),S_2(t)
c.......  S_x(b),S_xy(b),S_y(b),S_1(b),S_2(b),press
          st(ii,1)  = st(ii,1)  + sigi(1)*xsji
          st(ii,2)  = st(ii,2)  + sigi(2)*xsji
          st(ii,3)  = st(ii,3)  + sigi(3)*xsji
          st(ii,6)  = st(ii,6)  + sigb(1)*xsji            
          st(ii,7)  = st(ii,7)  + sigb(2)*xsji
          st(ii,8)  = st(ii,8)  + sigb(3)*xsji
          st(ii,11) = st(ii,11) + press  *xsji
          if(d(21).ne.0.d0.and.klay.eq.0) then
            st(ii,4)  = st(ii,4)  + as(1)*xsji
            st(ii,5)  = st(ii,5)  + as(2)*xsji
            st(ii,9)  = st(ii,9)  + as(3)*xsji
            st(ii,10) = st(ii,10) + as(4)*xsji
          else
            st(ii,4)  = st(ii,4)  + sigi(4)*xsji
            st(ii,5)  = st(ii,5)  + sigi(5)*xsji
            st(ii,9)  = st(ii,9)  + sigb(4)*xsji
            st(ii,10) = st(ii,10) + sigb(5)*xsji
          end if 
435     continue
440   continue
      return
      end
c
      subroutine proj09(s,p,zi,zj,nst)
      implicit double precision (a-h,o-z)
      dimension   s(nst,*),p(*)
c.... modify stiffness for offset projections
c.... postmultiply by transformation
      do 100 i = 1,6
        s(i,4) =  zj*s(i,2) + s(i,4)
        s(i,5) = -zj*s(i,1) + s(i,5)
100   continue
c.... premultiply using modified terms from postmultiplication
      do 200 i = 1,6
        s(4,i) =  zi*s(2,i) + s(4,i)
        s(5,i) = -zi*s(1,i) + s(5,i)
200   continue
      return
      end
c
      subroutine stre09s(d,sigto,sigbo,eps)
c--------------------------------------------------------------------
c     calculate stresses at top and bottom from resultants
c     ww bs uniKA  1/95
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension   sigi(6),sigb(6),d(*),sigto(6),sigbo(6),eps(6),eps1(6)
c.... stress resultants, eps,kappa
      do i = 1,6
         sigi(i) = sigto(i)
         sigb(i) = sigbo(i)
         eps1(i) =   eps(i)
      end do
c...  thickness values
      h = d(20)
      h2   = 0.5d0*h
      h312 = h*h*h/12.d0
c.... stresses at top (+h/2)
      sigto(1) = sigi(1)/h - sigb(1)/h312*h2
      sigto(2) = sigi(2)/h - sigb(2)/h312*h2
      sigto(3) = sigi(3)/h - sigb(3)/h312*h2
      call pstres(sigto,sigto(4),sigto(5),sigto(6))
c.... stresses at bottom (-h/2)
      sigbo(1) = sigi(1)/h + sigb(1)/h312*h2
      sigbo(2) = sigi(2)/h + sigb(2)/h312*h2
      sigbo(3) = sigi(3)/h + sigb(3)/h312*h2
      call pstres(sigbo,sigbo(4),sigbo(5),sigbo(6))
c.... strains at top (+h/2)
      eps(1) = eps1(1) - h2 * eps1(4)
      eps(2) = eps1(2) - h2 * eps1(5)
      eps(3) = eps1(3) - h2 * eps1(6)
c.... strains at bottom (+h/2)
      eps(4) = eps1(1) + h2 * eps1(4)
      eps(5) = eps1(5) + h2 * eps1(5)
      eps(6) = eps1(3) + h2 * eps1(6)
      return
      end
c
      subroutine elener09(d,sigi,sigb,da)
c----------------------------------------------------------------------
c.....sum energy for taylor shell element for error analysis
c----------------------------------------------------------------------
      USE errin1
      implicit double precision (a-h,o-z)
      dimension d(*),sigi(6),sigb(6)
c.... terms in el.matrix and its inverse, bend = memb/(h*h/12) = m/(d(15)
      d11    = d(1)
      d12    = d(2)
      d33    = d(3)
      det    = d11*d11 - d12*d12 
      di11   = d11/det
      di12   =-d12/det
      di33   = 1./d33
c.... membrane energy
      u_om(1) = u_om(1) + da *
     +        (di11*sigi(1)*sigi(1)+di11*sigi(3)*sigi(3)       
     +        +2.d0*di12*sigi(1)*sigi(3)+di33*sigi(2)*sigi(2))
      u_om(2) = u_om(2) + da * dot(sigi(1),sigi(1),3)
c.... bending energy  add with factor 1/d(15)!
      u_om(1) = u_om(1) + da / d(15) / d(15) *
     +        (di11*sigb(1)*sigb(1)+di11*sigb(3)*sigb(3)       
     +        +2.d0*di12*sigb(1)*sigb(3)+di33*sigb(2)*sigb(2))
      u_om(2) = u_om(2) + da / d(15) * dot(sigb(1),sigb(1),3)
      return
      end
c
      subroutine stcn09e(ix,d,yl,ul,tr,dt,st,dvl,ndf,nel,numnp,numel,
     +                   klay,ndm,e_ome)
c----------------------------------------------------------------------
c.....energy of stress differences for Taylor shell element
c----------------------------------------------------------------------
      USE errin1
      USE errin2
      implicit double precision (a-h,o-z)
      dimension   dt(numnp),st(numnp,*),yl(3,*),sigi(6),sigb(6),
     1       tr(3,3),d(*),ul(ndf,*),vl(6,4),eps(6),dvl(4),ix(*),
     +         sigpi(3),dsigi(3),sigpb(3),dsigb(3)
      dimension e_ome(*),e_ome09(numerr)
      common /shpf09/ shp(3,8,4),shp1(3,4,4),shp2(3,4,4)
!$OMP THREADPRIVATE (/shpf09/)  
c.... terms in el.matrix and its inverse, bend = memb/(h*h/12)=m/(d(15)
      d11    = d(1)
      d12    = d(2)
      d33    = d(3)
      det    = d11*d11 - d12*d12 
      di11   = d11/det
      di12   =-d12/det
      di33   = 1./d33
c
c.... displ. projection.
      do 400 i = 1,4
        do 400 j = 1,3
            vl(j  ,i) = 0.0
            vl(j+3,i) = 0.0
            do 399 k = 1,3
              vl(j  ,i) = vl(j  ,i) + tr(j,k)*ul(k  ,i)
              vl(j+3,i) = vl(j+3,i) + tr(j,k)*ul(k+3,i)
399         continue
400   continue
c.... sum element energy
c.... intial values for element errors
      e_ome09 = 0.0d0
      do 440 l = 1,4   ! Gauss-Points
c.....  Gauss Point stresses from mat. law
          call stre09(d,yl,vl,3,nel,l, xx,yy,zz,eps,sigi,sigb,press)
c.....  Gauss Point stresses  from nodal stress projection
        call pzero(sigpi,3)
        call pzero(sigpb,3)
          do 435 j = 1,4   ! Nodes
            fac  = shp(3,j,l)
            ii   = iabs(ix(j))
            if(ii.le.0) go to 435
            sigpi(1) = sigpi(1) + st(ii,1)*fac 
            sigpi(2) = sigpi(2) + st(ii,2)*fac 
            sigpi(3) = sigpi(3) + st(ii,3)*fac 
            sigpb(1) = sigpb(1) + st(ii,6)*fac 
            sigpb(2) = sigpb(2) + st(ii,7)*fac 
            sigpb(3) = sigpb(3) + st(ii,8)*fac 
435       continue
c.....  Diff. stresses  at gauss point
        do i = 1,3
          dsigi(i) = sigpi(i)-sigi(i)
          dsigb(i) = sigpb(i)-sigb(i)
        end do    
c.... membrane energy of stress differences
      e_ome09(1) = e_ome09(1) + dvl(l) * 
     +        (di11*dsigi(1)*dsigi(1)+di11*dsigi(3)*dsigi(3)       
     +        +2.*di12*dsigi(1)*dsigi(3)+di33*dsigi(2)*dsigi(2))
      e_ome09(2) = e_ome09(2) + dvl(l) * dot(dsigi(1),dsigi(1),3)
c.... bending energy of stress differences, add with factor 1/d(15)!
      e_ome09(1) = e_ome09(1) + dvl(l) / d(15) / d(15) * 
     +        (di11*dsigb(1)*dsigb(1)+di11*dsigb(3)*dsigb(3)       
     +        +2.*di12*dsigb(1)*dsigb(3)+di33*dsigb(2)*dsigb(2))
      e_ome09(2) = e_ome09(2) + dvl(l)/d(15) * dot(dsigb(1),dsigb(1),3)
440   continue
c.... plot/print stress errors
      call elmterr(ix,xl,ndm,numel,e_ome09,e_ome)
      return
      end
c
      subroutine design09(d,sigi,sigb,as)
c--------------------------------------------------------------------+
c     copy from elmt07s but with normal forces                       |
c     concrete design  (from general design diagram)                 |
c     mu     = gamma*m - gamma*n*zs                                  |
c     nu     = gamma*n                                               |
c     ms    = mu / (b*h**2*beta_R)                                   |
c     mbu   = mu <=0.337                                             |
c     meu   = mu - 0.337 >= 0                                        |
c     k_z   = 0.5 + sqrt(0.25 - 0.5*mbu           (Approximation)    |
c     asz   = mbu /(k_z*h*beta_S)+nu/beta_S+[meu /((h*)*beta_S)]     |
c     asd   =                               [meu /((h*)*beta_S)]     |
c     b     = 1                                                      |
c     gamma = 1.75                                                   |
c     h     = d-hx(y)                                                |
c     h*    = h-hx(y)                                                |
c--------------------------------------------------------------------+
      implicit double precision (a-h,o-z)
      dimension d(*),sigi(6),sigb(6),as(4)
      smg    = 0.337d0
      thickd = d(20)
      betar  = d(21)
      betas  = d(22)
      thickx = d(23)
      thicky = d(24)
      gamma  = 1.75d0
      asxz   = 0.d0
      asyz   = 0.d0
      asxd   = 0.d0
      asyd   = 0.d0
      call pzero(as,4)
c.... design
c.....m_x
c.... distance for normal force
      thickzs = thickx-0.5d0*thickd   
c.... res. moment 
      smx = sigb(1) - sigi(1)*thickzs
      if(smx .ge. 0.0d0) smxu =      smx  * gamma
      if(smx .lt. 0.0d0) smxu = dabs(smx) * gamma
      smxb     = smxu /(thickx*thickx*betar)
      if(smxb  .gt. smg) then
c....   add. steel part 
        dsmx   = (smxb  - smg)
        smxeu  =  dsmx*thickx*thickx*betar
        smxu   =  smxu - smxeu
        smxb   =  smg 
c....   approx. distance for steel h-h'  
        dhx   = thickx - (thickd-thickx)
c....   tension
        asxz  = asxz + smxeu/(betas*dhx)
c....   compression
        asxd  = asxd + smxeu/(betas*dhx)
      end if
c...  find k_z
      val   = 0.25d0 - 0.5d0*smxb
      zk    = 0.5d0 + dsqrt(val)
c....   steel for concrete part 
      asxz  = asxz  + smxu / (zk*thickx*betas) + sigi(1)*gamma/betas
      if(smx .ge. 0.0d0) then
        as(1) = as(1) + max(asxz,0.d0)
        as(2) = as(2) + max(asxd,0.d0)
      elseif(smx .lt. 0.0d0) then
        as(2) = as(2) + max(asxz,0.d0)
        as(1) = as(1) + max(asxd,0.d0)
      end if
c.....m_y
c.... distance for normal force
      thickzs = thicky-0.5d0*thickd   
c.... res. moment 
      smy = sigb(3) - sigi(3)*thickzs
      if(smy .ge. 0.0d0) smyu =      smy  * gamma
      if(smy .lt. 0.0d0) smyu = dabs(smy) * gamma
      smyb    = smyu /(thicky*thicky*betar)
      if(smyb .gt. smg) then
c....   add. steel part 
        dsmy  = (smyb - smg)
        smyeu =  dsmy*thicky*thicky*betar
        smyu  =  smyu - smyeu
        smyb  =  smg
c....   approx. distance for steel h-h'  
        dhy   = thicky - (thickd-thicky)
c....   tension
        asyz  = asyz + smyeu/(betas*dhy)
c....   compression
        asyd  = asyd + smyeu/(betas*dhy)
      end if
c...  find k_z
      val   = 0.25d0 - 0.5d0*smyb
      zk    = 0.5d0 + dsqrt(val)
c....   steel for concrete part 
      asyz  = asyz  + smyu / (zk*thicky*betas) + sigi(3)*gamma/betas
      if(smy .ge. 0.0d0) then
        as(3) = as(3) + max(asyz,0.d0)
        as(4) = as(4) + max(asyd,0.d0)
      elseif(smy .lt. 0.0d0) then
        as(4) = as(4) + max(asyz,0.d0)
        as(3) = as(3) + max(asyd,0.d0)
      end if
      return
      end
c
      subroutine qload09(ql,d,q,numel,n,mqloa,propq,prop,isw)
c----------------------------------------------------------
c.... add loads from macro qloa
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),d(*),ql(10)
      call pzero(ql,10) 
      if(isw.eq.22) then
        if(mqloa.ne.1) then
          do i = 1,6 
            ql(i) = q(n,i)*propq 
          end do 
        end if
      else 
        do i = 1,6 
          ql(i) = d(6+i)*prop 
        end do 
      end if
      return
      end
c
      subroutine plloco09(tr,xl)
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
      call pltaxs09(tr,x0,xm)
      return
      end
c
      subroutine pltaxs09(tr,x0,xm)
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
