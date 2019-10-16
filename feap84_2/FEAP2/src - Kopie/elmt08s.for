      subroutine elmt08(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-------------------------------------------------------+
c     2-4 node shell of revolution element              |
c        "mindlin theory"                               |
c     linear elastic element                            |
c     1991-06  basic             K. Knebel              |
c     1997-06  thermal forces    ww/fg                  |
c     1997-10  curvatures chgd.  ww/fg                  |
c     1997-11  snow loads        kub                    |
c     1999-01  optimization      ww                     |
c     2010-12  elast. foundation ww                     |
c     2014-12  kappa             ww                     |
c-------------------------------------------------------+
c
      USE bdata
      USE cdata
      USE eldata
      USE fornam
      USE iofile
      USE pdata6
      USE pdata7
      USE pdata10
      USE strnam
      implicit real*8 (a-h,o-z)
      character wd(3)*8
      integer ix(*),ixl(5)
      dimension xl(ndm,*),tl(*),pel1(12),pel2(12),shp(4),dsh(4),xll(3,4)
      dimension d(*),ul(*),s(nst,*),p(*),eps(5)
      dimension bmat(5,12),cbmat(12),cb(5,12),sig(5),sigr(5),dd(6)
      dimension h1(*),h2(*),h3(*)
      dimension sg(5),wg(5)
      common /ec08a08/   cmat(5,5),hyp,co,si,cob,sib
!$OMP THREADPRIVATE (/ec08a08/)  
c
      data wd/'2 - node','3 - node','4 - node'/
      data ixl/1,2,3,4,1/
      pi     = 4.d0*datan(1.d0)
      igausb = d(7)
      igauss = d(8)
      ityp   = d(4)
      ielno  = 8
c
c.... go to correct array processor
c
      go to(1,2,3,3,2,3,2,2,2,2,2,2,3), isw
      return
c----------------------------------------------------------
c            isw = 1     material data input
c----------------------------------------------------------
1     continue
c.... input material properties
c----------------
c1     1=e
c1     2=xnu
c1     3=h=thick
c      4=ityp studvers=1   
c1     5=cb
c1     6=rho
c1    15=cappa 
c---------------
c      7=igausb studvers=1
c      8=igauss studvers=1
c---------------
c2     9=s-load 
c2    10=g-load
c2    11=pc 
c2    12=po
c2    13=zniv
c2    18=dir
c---------------
c3    14=alp
c3    16=tb
c3    17=t

c---------------------------------------------- 1 material data
      if(ior.lt.0) write(*,5001) !e,nu,h,rho,cb,cappa
      call dinput(dd,6)
      d(1)  = dd(1)
      d(2)  = dd(2)
      d(3)  = dd(3)
      d(6)  = dd(4)
      d(5)  = dd(5)
      d(15) = dd(6)
c------------------------------------------------------ 2 loads
      if(ior.lt.0) write(*,5002) !g,pc,p0,z0,idir,s
      call dinput(dd,6)
      d(10) = dd(1)
      d(11) = dd(2)
      d(12) = dd(3)
      d(13) = dd(4)
      d(18) = dd(5)
      d(9)  = dd(6)
      if(d(18).ne. 1.0d0) d(18)=2.0d0
c---------------------------------------------- 3 thermal loads
      if(ior.lt.0) write(*,5003)! alpa,temp-bot,temp-top
      call dinput(dd,3)
      d(14) = dd(1)
      d(16) = dd(2)
      d(17) = dd(3)
c
c...  fixed for stud.version
      d(4)  = 1.0d0
      d(7)  = 1.0d0
      d(8)  = 1.0d0
c.... modification for linear pressure load: 2 point integration
      if(d(13).ne.0.d0) d(7) = 2.0d0
c.... set material parameter type and flags
      write(iow,2000) wd(1),d(1),d(2),d(6),d(5),d(15),int(d(7)), 
     1     int(d(8)),d(3),d(10),d(11),d(12),d(13),int(d(18)),d(9),
     +      d(14),d(16),d(17)
      if(ior.lt.0) then
        write(*,2000) wd(1),d(1),d(2),d(6),d(5),d(15),int(d(7)), 
     1      int(d(8)),d(3),d(10),d(11),d(12),d(13),int(d(18)),d(9),
     +      d(14),d(16),d(17),d(5)
      end if
      ipla = 2
c.... description of stresses  
      forsus( 1) =  '  N-FORCE n_s  '
      forsus( 2) =  '  N-FORCE n_t  '
      forsus( 3) =  '  MOMENT  m_s  '
      forsus( 4) =  '  MOMENT  m_t  '
      forsus( 5) =  '  Q-FORCE q_s  '
      do i = 6,11         
         forsus(i) =  ' '
      end do
c.... define node numbering for plot mesh routine, see pltord
      inord(ielno)   = 4
      ipord(1,ielno) = 1
      ipord(2,ielno) = 2
      ipord(3,ielno) = 2
      ipord(4,ielno) = 1
      return
c----------------------------------------------------------
c     isw = 2  element check (not in use)
2     return
c----------------------------------------------------------
c     isw = 3  stiffness matrix and element-load-vector
c----------------------------------------------------------
3     continue
      dx  = xl(1,nel) - xl(1,1)
      dy  = xl(2,nel) - xl(2,1)
      hyp = dsqrt ( dx * dx + dy * dy)
c.....direction cosine for transformation
      si  = dy / hyp
      co  = dx / hyp
c.....direction cosine for b-matrix
      sib = co
      cob = si
c.....set material-matrix cmat(5x5)
      call shearfac08(d,hyp,cappa)
      call cmat08(cmat,d,cappa)
c
      if(isw .eq. 4) goto 4
      if(isw .eq.13) goto 4
c
c.... stiffness / residual computation
c.... full or selective integration for shear
c
      if(igausb .eq. igauss) then
        kstre = 5              ! full integration
      else
        kstre = 5-1            ! selective integration for shear
      endif
c
c.... integration for bending and membran (no.gauss-points bending = igausb)
c.... also for shear if igauss = igausb
c
      call gauss08(igausb,lint,sg,wg)
c
      do l=1,lint
        call shap08(sg(l),xl,shp,dsh,ndm,nel,detj)
        radius=glin(xl,sg(l),ndm,nel,1)
        dv = 2*pi*detj*wg(l)*radius
        call bmat08(bmat,cbmat,shp,dsh,radius,sib,cob,nel,nst)
        call matmulf(cmat,bmat,5,5,nst,cb)

        do ievab=1,nst
          do jevab=ievab,nst
c....       s = int bt x cb x dv
            do istre=1,kstre
              s(ievab,jevab) = s(ievab,jevab) +
     1                     bmat(istre,ievab)*cb(istre,jevab)*dv
            end do ! istre 
c....       add elastic foundation term 
            s(ievab,jevab) = s(ievab,jevab) +
     1                     cbmat(ievab)*cbmat(jevab)*d(5)*dv
          end do ! jevab 
        end do ! ievab 
      end do ! l
c
      if(igausb .eq. igauss) goto 350
c.... integration for shear   (no.gauss-points shear = igauss)
c.... done only if   igauss .ne. ibausb
      call gauss08(igauss,lint,sg,wg)
      do l=1,lint
        call shap08(sg(l),xl,shp,dsh,ndm,nel,detj)
        radius = glin(xl,sg(l),ndm,nel,1)
        dv = 2*pi*detj*wg(l)*radius
        call bmat08(bmat,cbmat,shp,dsh,radius,sib,cob,nel,nst)
        call matmulf(cmat,bmat,5,5,nst,cb)

c....   s = int bt x cb x dv
        do ievab=1,nst
          do jevab=ievab,nst
            do istre=5,5
              s(ievab,jevab) = s(ievab,jevab) +
     1                         bmat(istre,ievab)*cb(istre,jevab)*dv
            end do ! istre 
          end do ! jevab 
        end do ! ievab 
      end do ! l
350   continue
c.... make stiffness symmetric
      do 330 i=1,nst
        do 330 j=i,nst
          s(j,i) = s(i,j)
330   continue
c
c.... transform stiffness to global
      call trans08 ( s,co,si,nst,ndf,1 )
c 
c.... form element-load vektor
      help=abs(d(9))+abs(d(10))+abs(d(11))+abs(d(12))+abs(d(14))
      if(help .ne. 0.d0) then
c...    local vector
        call pzero(pel1,12)
        call eload081(pel1,d,nel,xl,tl,ndm,nst,ndf,pi)
c....   transform element load vector to global
        call trans08(pel1,co,si,nst,ndf,2)
c...    global vector
        call pzero(pel2,12)
        call eload082(pel2,d,nel,xl,tl,ndm,nst,ndf,pi)
c...    store into p        
        do 340 inst=1,nst
          p(inst) = p(inst)+pel1(inst)+pel2(inst)
340     continue
      end if
c
cww      if(isw.eq.6) then  ! if not comment: no convergence and add per tang 1 lin.sol.
c....   compute a residual 
        do 370 j=1,nst
          do 370 k=1,nst
            p(j) = p(j) - s(j,k)*ul(k)
370     continue
cww      end if
c
      return
c-----------------------------------------------------------------
c     isw = 4  compute strain and stresses at gauss-points (shear)
c-----------------------------------------------------------------
4     call gauss08(igauss,lint,sg,wg)
c
c.... transform global disp to local
      call trans08(ul,co,si,nst,ndf,3)
c
      thick = d(3)
      do l=1,lint
        call shap08(sg(l),xl,shp,dsh,ndm,nel,detj)
        radius = glin(xl,sg(l),ndm,nel,1)
        call stre08(d,xl,ul,tl,shp,dsh,radius,ndm,sigr,eps,nst,nel,ndf)
        if(isw .eq. 4) then
c....    stresses from resultants
         sig(1) = sigr(1)/thick
         sig(2) = sigr(2)/thick
         sig(3) = sigr(3)*6/(thick*thick)
         sig(4) = sigr(4)*6/(thick*thick)
         sig(5) = sigr(5)/thick
c
c... output strain , (stress) and stress-resultants
c
        mct=mct-3
        hight=glin(xl,sg(l),ndm,nel,2)
        if(mct .le. 0) then
          write(iow,2001) o,head
          if(ior .lt. 0) write(*,2001) o,head
          mct = 50
        endif
c                     write(iow,2002) n,radius,eps,ma,hight,sig,sigr
c        if(ior.lt.0) write(*  ,2002) n,radius,eps,ma,hight,sig,sigr
                     write(iow,2002) n,ma,radius,hight,sigr
        if(ior.lt.0) write(*  ,2002) n,ma,radius,hight,sigr
      end if
      end do ! l
c
c.... plot p,v,m diagrams on frame  ! only 2 node element
      if(isw.eq.13) then
        if(nfp.lt.1.or.nfp.gt.5) return
        klayf = 1
        if(flfp) then
          ccfp  = sigr(nfp)
          xmaxf = max(xmaxf,ccfp)
          xminf = min(xminf,ccfp)
          ccfp  = abs(sigr(nfp))
          cfp   = max(cfp,ccfp)
        else
          call pzero(xll,12)
          xll(1,1) = xl(1,1)
          xll(2,1) = xl(2,1)
          xll(1,2) = xl(1,nel)
          xll(2,2) = xl(2,nel)
          xll(1,3) = xl(1,nel) + si*sigr(nfp)*cfp
          xll(2,3) = xl(2,nel) - co*sigr(nfp)*cfp
          xll(1,4) = xl(1,1)   + si*sigr(nfp)*cfp
          xll(2,4) = xl(2,1)   - co*sigr(nfp)*cfp
cwwc..... choose color: blue = +  red = -
cww       if(sigr(nfp).gt.0.0d0) then
cww         call pppcol(3)
cww       else
cww         call pppcol(2)
cww       endif
c
          call pppcolf(sigr(nfp))
c
          call plot9s(ixl,xll,3,4)
        endif
      endif
      return
c
c.... formats for input-output
c
2000  format(/,5x,a12,' linear elastic shell element of revolution',//
     1 10x,'E-modulus       E....',g16.5/
     2 10x,'poisson ratio   nu...',g16.5/
     6 10x,'density         rho..',g16.5/
     5 10x,'Modul C_B       cb...',g16.5/
     5 10x,'SCF             kappa',g16.5/
     7 10x,'gauss pts bend  igb..',i3/
     8 10x,'gauss pts shear igs..',i3/
     3 10x,'thickness       h....',f16.5,//
     + 10x,'dead-load       g....',g16.5/
     1 10x,'press-const     pc...',g16.5/
     2 10x,'press-slope     p0...',g16.5/  
     3 10x,'z-niveau        z0...',g16.5,' in ',i2,' direction'/
     9 10x,'snow-load       s....',g16.5/
     4 10x,'alpha_t         a_t..',g16.5/
cww  + 10x,'base-temp            ',g16.5/
     6 10x,'temp-bottom     tb...',g16.5/
     7 10x,'temp-top        tt...',g16.5)
c
c2001  format(a1,20a4,//,5x,'element strains,stresses and resultants',//,
c     +  ' elmt  radius',
c     +  2x,'ss-strain',2x,'tt-strain',2x,'s-curvat ',2x,'t-curvat ',
c     +  2x,'shear-str',/,' matl  hight ',2x,'ss-stress',
c     +  2x,'tt-stress',2x,'ms-stress',2x,'mt-stress'2x'shear-stress',/,
c     +  15x,'s-force  ',2x,'t-force  ',2x,'s-moment',2x,' t-moment',
c     +  3x,'shear-force',/,39(' -'))
c2002  format(i5,0p1f9.3,1p5e11.3/i5,0p1f9.3,1p5e11.3/
c     +       ,'              ',1p5e11.3/)
2001  format(a1,20a4,/,' stress resultants for axisymmmetric shell',/,
     +  ' elmt matl   x_1      x_2    ',
     +  '  ns-force ','  nt-force ','  ms-moment','  mt-moment',
     +  '  qs-force',/)
c
2002  format(i5,i5,0p2f9.3,1p5e11.3)
5001  format(' input: e, nu, thick,rho',/3x,'> ',$)
5002  format(' input: deadload, press-const, press-basis, z-nineau,'
     +       /3x,'>',$)
5003  format(' input: alpa,temp-bast,temp-bottom,temp-top,'/3x,'>',$)
cww5003  format(' input: elm-typ (1=2-node  2=3-node 3=4-node,'
cww     +        ,'g-pts bending, g-pts shear ',/3x,'> ',$)
      end
c
      subroutine cmat08(c,d,cappa)
c----------------------------------------------------------
c      compute material-matrix c  (size 5x5)
c----------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension c(5,5),d(*)
      e    = d(1)
      xnu  = d(2)
      h    = d(3)
      g    = e/(2.d0*(1.d0+xnu)) 

      c = 0.d0
      c(1,1) = e*h/(1.d0-xnu*xnu)
      c(2,2) = c(1,1)
      c(1,2) = c(1,1)*xnu
      c(2,1) = c(1,2)
      c(3,3) = (e*h*h*h)/(12.d0*(1.d0-xnu*xnu))
      c(4,4) = c(3,3)
      c(4,3) = c(3,3)*xnu
      c(3,4) = c(4,3)
      c(5,5) = cappa*g
      return
      end
c
      subroutine bmat08(bmat,cbmat,shp,dsh,r,si,co,nel,nst)
c----------------------------------------------------
c     set b-matrix
c----------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension bmat(5,nst),shp(*),dsh(*),cbmat(nst)
      bmat = 0.d0
c
c.... evaluate B-matrix
      j3 = 0
      do inode=1,nel
        j1=j3+1
        j2=j1+1
        j3=j1+2
        bmat(1,j1) =  dsh(inode)
        bmat(2,j1) =  shp(inode)*si/r
        bmat(2,j2) = -shp(inode)*co/r
        bmat(3,j3) = -dsh(inode)
        bmat(4,j3) = -shp(inode)*si/r
        bmat(5,j2) = -dsh(inode)
        bmat(5,j3) = -shp(inode)
      end do
c.... evaluate N-matrix
      j1=-1
      cbmat = 0 
      do inode=1,nel
        j1=j1+3
        cbmat(j1) = shp(inode)
      end do
      return
      end
c
      subroutine stre08(d,xl,ul,tl,shp,dsh,radius,ndm,sigr,eps,nst,nel,
     +                  ndf)
c----------------------------------------------------
c     compute strains, stresses and stress-resultants
c----------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension d(*),xl(ndm,*),ul(*),tl(*),shp(*),dsh(*)
      dimension sigr(5),eps(5),bmat(5,12),cbmat(12)

      common /ec08a08/ cmat(5,5),hyp,co,si,cob,sib
!$OMP THREADPRIVATE (/ec08a08/)  
c
c.....zero stress and strain vector
      sigr=0
      eps =0
c
c.... set up b-matrix
      call bmat08(bmat,cbmat,shp,dsh,radius,sib,cob,nel,nst)
c
c.... compute strain eps = b * u
      call mvmul(bmat,ul,5,nst,eps)
cwwfg>change curvatures!!
      eps(3) = -eps(3)
      eps(4) = -eps(4)
cwwfg<
c.... compute stress-resultants  sigr = c * eps
      call mvmul(cmat,eps,5,5,sigr)
c
c...  compute thermal stresses
c
      if(d(14) .gt. 1.0d-15) then
        ctemp = d(1)*d(3)*d(14)/(1.0-d(2))
        btemp = ctemp*(d(3)*d(3))/12.d0
        tg = 0.5d0*(d(16)+d(17))
        do 110 inod = 1,nel
          tg = tg + tl(inod)*shp(inod)
110     continue
        sigt1 = ctemp * tg
        tdiff = d(16)-d(17)
        sigt2 = btemp*tdiff/d(3)
c
c.....add stresses due to initial strain (e.g. thermal strain)
c
        sigr(1) = sigr(1) - sigt1
        sigr(2) = sigr(2) - sigt1
cwwfg>sigr(3) = sigr(3) - sigt2   change due to def. of m_s
c       sigr(4) = sigr(4) - sigt2
      sigr(3) = sigr(3) + sigt2
        sigr(4) = sigr(4) + sigt2
cwwfg<
      end if
      return
      end
c
      subroutine eload081(pel,d,nel,xl,tl,ndm,nst,ndf,pi)
c--------------------------------------------------------------------
c     compute local element-load vectors
c--------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension pel(*),d(*),xl(ndm,*),tl(*)
      dimension sg(2),wg(2),shp(4),dsh(4)
      common /ec08a08/ cmat(5,5),hyp,co,si,cob,sib
!$OMP THREADPRIVATE (/ec08a08/)  
c
c.... set integation-order
      k   = 2
c
      pconst = d(11)
      po     = d(12)
      zniv   = d(13)
      thick  = d(3)
      alp    = d(14)
      tb     = d(16)
      tt     = d(17)
      idir   = d(18)
      ctemp  = d(1)*d(3)*d(14)/(1.0-d(2))
      btemp  = ctemp*(thick*thick)/12.
c
      call gauss08(k,lint,sg,wg)
c
      do 100 l=1,lint
        call shap08(sg(l),xl,shp,dsh,ndm,nel,detj)
        radius = glin(xl,sg(l),ndm,nel,1)
        hight  = glin(xl,sg(l),ndm,nel,2)
        dv     = 2*pi*detj*wg(l)*radius
c
c.... pressure due to hydrostatic load
c
        if(idir.eq.1) dept   = zniv-radius
        if(idir.eq.2) dept   = zniv-hight
        if(dept .le. 1.0d-10) then
          phw = 0.0d0
        else
          phw = po*dept
        endif
c
c.... pressure due to thermal load
c
cww     tg = 0.5*(tt+tb)-t0
        tg = 0.5*(tt+tb)
        do 110 inod = 1,nel
          tg = tg + tl(inod)*shp(inod)
110     continue
        sigt1 = ctemp *tg
        tdiff = tb-tt
        sigt2 = btemp*tdiff/thick
c
c.... sum  pressure due to hydrostatic and const.pressure
c
        press2 = phw + pconst
c
        j = 1
        do 120 inod=1,nel
          pel(j+1) = pel(j+1)+shp(inod) * press2 * dv
c...      add thermal load
          if((abs(sigt1)+abs(sigt2)) .gt. 1.0d-15) then
            pel(j)   = pel(j)   + (dsh(inod)+shp(inod)*sib/radius)
     +                          * sigt1 * dv
            pel(j+1) = pel(j+1) + sigt1*(-shp(inod))*cob/radius * dv
            pel(j+2) = pel(j+2) + (-dsh(inod)-shp(inod)*cob/radius)
     +                          * sigt2 * dv
          endif
          j=j+3
120     continue
100   continue
      return
      end
c
      subroutine eload082(pel,d,nel,xl,tl,ndm,nst,ndf,pi)
c--------------------------------------------------------------------
c     compute global element-load vectors
c--------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension pel(*),d(*),xl(ndm,*),tl(*)
      dimension sg(2),wg(2),shp(4),dsh(4)
      common /ec08a08/ cmat(5,5),hyp,co,si,cob,sib
!$OMP THREADPRIVATE (/ec08a08/)  
c
c.... set integation-order
      k   = 2
      call gauss08(k,lint,sg,wg)
c
      do 100 l=1,lint
        call shap08(sg(l),xl,shp,dsh,ndm,nel,detj)
        radius = glin(xl,sg(l),ndm,nel,1)
        hight  = glin(xl,sg(l),ndm,nel,2)
        dv     = 2*pi*detj*wg(l)*radius
c
        press  = 0.d0
c....   pressure due to snow load
c
        sload  = d(9)
        if(abs(sload) .gt. 1.0d-10) then
          press = press - sload *dabs(co)
        endif
c
c.... pressure due to deadload
c
        dload  = d(10)
        if(dload .gt. 1.0d-10) then
          thick = d(3)
          press = press -dload * thick
        endif
c
c.... set load vector
c
        j = 1
        do 120 inod=1,nel
          pel(j+1) = pel(j+1)+shp(inod) * press * dv
          j=j+3
120     continue
100   continue
      return
      end
c
      subroutine shap08 (s,xl,shp,dsh,ndm,nel,detj)
c------------------------------------------------------------------
c     shape functions one-dimensional c-0 element
c     the 3rd node in the quadric elem. must be in the middle
c------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension xl(ndm,nel),sl(4),shp(nel),dsh(nel)
c
      dx   = xl(1,nel) - xl(1,1)
      dy   = xl(2,nel) - xl(2,1)
      eleng= sqrt(dx*dx + dy*dy)
      detj = 0.0d0
      go to (1,2,3) nel-1
c
c.... linear element  2-nodes
c
    1 shp(1) = (1.0-s)*0.5
      shp(2) = (1.0+s)*0.5
      detj   = 0.5*eleng
      dsh(1) =-0.5
      dsh(2) = 0.5
      goto 10
c
c.... quadratic element  3-nodes
c
    2 shp(1) = (s*s-s)*0.5
      shp(2) = 1.0-s*s
      shp(3) = (s*s+s)*0.5
      detj   = 0.5*eleng
      dsh(1) = (2.0*s-1.0)*0.5
      dsh(2) = -2.0*s
      dsh(3) = (2.0*s+1.0)*0.5
      goto 10
c
c.... cubic element  4-nodes
c
    3 continue
      do 100 i=1,2
        dx = xl(1,nel-i) - xl(1,1)
        dy = xl(2,nel-i) - xl(2,1)
        sx = sqrt(dx*dx + dy*dy)
100    sl(nel-i) = eleng -sx
      sl(1) = 0.0d0
      sl(4) = eleng

      s2 = s*s
      s3 = s2*s
      shp(1) = -9.0/16.0*(s3-s/9.0-s2+1.0/9.0)
      shp(2) = 27.0/16.0*(s3-s-s2/3.0+1.0/3.0)
      shp(3) =-27.0/16.0*(s3-s+s2/3.0-1.0/3.0)
      shp(4) =  9.0/16.0*(s3-s/9.0+s2-1.0/9.0)
      dsh(1) = -9.0/16.0*(3.0*s2-1.0/9.0-2.0*s)
      dsh(2) = 27.0/16.0*(3.0*s2-1.0-2.0/3.0*s)
      dsh(3) =-27.0/16.0*(3.0*s2-1.0+2.0/3.0*s)
      dsh(4) =  9.0/16.0*(3.0*s2-1.0/9.0+2.0*s)
c
c.....compute jacobian-det
      do 20 i=1,nel
        detj=detj + dsh(i)*sl(i)
  20    continue
        if(detj .lt. 1.0d-10) then
          write(*,*) '***error*** jacb-det lower/equal zero'
          stop
        endif
c
c.....calculate cartesian derivatives
c
10    continue
      do 30 inode = 1,nel
        dsh(inode) = dsh(inode)/detj
30    continue
      return
      end
c
c------------------------------------------------------------------
      subroutine gauss08(l,lint,s,w)
      implicit real*8 (a-h,o-z)
      dimension s(*),w(*)
c
      lint = l
      if (lint.gt.5) go to 6
      go to (1,2,3,4,5) lint
c.... 1 point integration
    1 s(1) = 0.d0
      w(1) = 2.d0
      return
c.... 2 point integration
    2 s(1) = -0.577350269189626d0
      s(2) = -s(1)
      w(1) = 1.d0
      w(2) = w(1)
      return
c.... 3 point integration
    3 s(1) = -0.774596669241483d0
      s(2) = 0.d0
      s(3) = -s(1)
      w(1) = 0.555555555555556d0
      w(2) = 0.888888888888889d0
      w(3) = w(1)
      return
c.... 4 point integration
    4 s(1) = -0.861136311594053d0
      s(2) = -0.339981043584856d0
      s(3) = -s(2)
      s(4) = -s(1)
      w(1) = 0.347854845137454d0
      w(2) = 0.652145154862546d0
      w(3) = w(2)
      w(4) = w(1)
      return
c.... 5 point integration
    5 s(1) = -0.906179845938664d0
      s(2) = -0.538469310105683d0
      s(3) =  0.d0
      s(4) = -s(2)
      s(5) = -s(1)
      w(1) = 0.236926885056189d0
      w(2) = 0.478628670499366d0
      w(3) = 0.568888888888889d0
      w(4) = w(2)
      w(5) = w(1)
      return

    6 write(6,2000) lint
      write(*,2000) lint
 2000 format(' ** error **  no',i3,' point integration allowed'/)
      stop
      end
c
c-------------------------------------------------------------------------
      subroutine trans08(s,cs,snt,nst,ndf,itype)
      implicit real*8 (a-h,o-z)
c
c.... itype:
c        1  transform matrix s(nst,nst) from local  to global
c        2  transform vector s(nst,1)   from local  to global
c        3  transform vector s(nst,1)   from global to local
c
      dimension s(nst,*)
c
      if(cs.gt.0.999999d0) return
      sn = snt
      go to (1,2,3) itype
1     do 12 i = 1,nst,ndf
        j = i + 1
        do 11 n = 1,nst
          t      = s(n,i)*cs - s(n,j)*sn
          s(n,j) = s(n,i)*sn + s(n,j)*cs
          s(n,i) = t
11      continue
12    continue
      do 14 i = 1,nst,ndf
        j = i + 1
        do 13 n = 1,nst
          t      = s(i,n)*cs - s(j,n)*sn
          s(j,n) = s(i,n)*sn + s(j,n)*cs
          s(i,n) = t
13      continue
14    continue
      return
2     sn = -sn
3     do 31 i=1,nst,ndf
        j = i + 1
        t      =  s(i,1)*cs + s(j,1)*sn
        s(j,1) = -s(i,1)*sn + s(j,1)*cs
        s(i,1) =  t
31    continue
      return
      end
c
c----------------------------------------------------------------------
      double precision function glin (xl,ps,ndm,nel,k)
      implicit real*8 (a-h,o-z)
c
c.... linear interpolation between end-nodes
c.... k=1 => x-direction   k=2 => y-direction
c
      dimension xl(ndm,*)
c
      glin = 0.d0
      an1  = 0.5d0 - 0.5d0*ps
      an2  = 0.5d0 + 0.5d0*ps
      glin = xl(k,1)*an1 + xl(k,nel)*an2
      return
      end
c
      subroutine shearfac08(d,hyp,cappa)
c-----------------------------------------------------------------------
c.... shear correction factor due to size
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  d(*)
      cappa = d(15)
      if(cappa.lt.0.d0) then
        xnu = d(2)   
        hs  = d(3)
        ae = hyp*hyp
        cappa = dabs(d(15))/(1.d0+ae/(2.d0*hs*hs*(1.d0+xnu))) ! Tessler
      else if(cappa.eq.0.d0) then
        cappa =5.0d0/6.0d0          
      end if
      return
      end  
 