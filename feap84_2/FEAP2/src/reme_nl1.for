      subroutine reme_nl(u,dr,
     +      iek,x,ike,ikz,ianp,iael,
     +      x0,iek0,ike0,ikz0,iael0,ianp0,mikno,
czr  +      erron0,erros,erroe,errof,erro,
     +      erron0,e_ome, erro,
     +      u0,hst0,
     +      nen1,ndm,ndf,lct,ct,fint2)
c-----------------------------------------------------------------------
c.... adaptive mesh refinement for
c.... 4-node shell-element
c.... adjustments for nonlinear analysis: 11/96 ziegler)
c-----------------------------------------------------------------------
      USE cdata
      USE comfil
      USE errin1
      USE errin2
      USE errnam
      USE fdata
      USE iofile
      USE ldata
cba      USE sdata
      implicit double precision (a-h,o-z)
      real*8  ct(3,*)
      real*8 u(ndf,*),dr(*)               ! --> SR formfe
      logical pcomp,l_four
      logical exst,hflgu
      logical fa                           ! false
      character*4  unif,adap,fixd,frac
      character*4  lct(*)
      character*229 fint2
cba   integer*2 l1(5),l2(5)
czr   logical lad,pcomp,lae
      common /adap1/ naiter,mnph,mnphdt,nhist,nhsw
      common /gener/ kk,km,ke
cba      common /gene0/ kk0,ke0
cba      common /gefehl/ nfehler
cba      common /adap1/  lad,naiter,lae
cba      common /adap3/  n9e,asteu,fsteu,fmax,nerro
cba      common /adap4/  nadap
cba      common /prio1/  npiter,nprior
cba      common /tola/  tolp
cba      common /poll /  npoll0
czr   dimension xx(6),nk(5),mikno(3,*),iek0(nen1,*)
czr   dimension erro(*),erron(*),erron0(*)
      dimension x(ndm,*),x0(ndm,*),
     +      iek(nen1,*),ikz(*),ike(*),iael(*),ianp(*),
     +      mikno(3,*),iek0(nen1,*),
     +      iael0(*),ianp0(*),ike0(*),ikz0(*),
     +      erron0(*),
     +      erro(*),e_ome(numel,2),                !......
     +      u0(*)
      data unif /'unif'/ , adap /'adap'/, fixd/'fixd'/, frac/'frac'/
      data fa /.false./
      call pzero (u0,5*numel*ndf)                  ! geom  nonlin
      call pzero (erro,numnp)                      !
c     call pzero (erros,numnp)                     !
c     call pzero (erroe,numnp)                     !
c     call pzero (errof,numnp)                     !
      call pzero (e_ome,2*numnp)                   !
      call pzeroi (iael,numel)                     !
      call pzeroi (ianp,numnp)                     !
c..................................................
C.... read restart file
      inquire(file = fint2,exist=exst)
      if( exst ) then
czr     call rmstrt (fint2,ndm,ndf,nen1,1,1,     ! ascii
        call rmstrt (fint2,ndm,ndf,nen1,1,0,     ! binary
     1        kk,ke,naiter,iek,iek0,iael,iael0,ianp,ianp0,
     2        ct)                     
      else
        call matcoi(ianp,numnp,1,ianp0)
        call matcoi(iael,numel,1,iael0)
        call matcoi(iek,nen1,numel,iek0)
      endif
c.... initialisation
      kf=1                             ! used in xymod
      ke0=numel
      kk0=numnp
      ke=numel
      kk=numnp
      naiter=naiter+1
      write(*,3002)   naiter
      write(iow,3002) naiter
c.... ->  iek0 = iek
      call matcoi(ike,numel*(nen  ),1,ike0)
      call matcoi(ikz,numnp+1,1,ikz0)
c.... (generieren des grundnetzes)
      continue
c...........................................check element for error +1/-1
      if (naiter.gt.1) then
        if (pcomp(lct(l),adap,4)) then
c.... check stress errors
          if(fl(11)) then
c.... with respect to 'eval' percent of energy (default eval = 5 perc.)
            eval = ct(1,l)
            if(eval.lt.0.0001) eval = 5.          
            do ierror = 1,numerr             
              e_om(ierror) = 0.0d0
              e_bar(ierror)  = eval/100.*sqrt(u_om(ierror)/numel)
            enddo
            iet(1) = 3                         ! switch set to refinement
c.... save information on file (ioerr = 1)
            ioerr = 0
            hflgu = .false.
c.... loop over elements
c            ietold=iet       ! avoid printing of error
c            iet = 4
            call formfe(u,dr,dr,dr,fa,fa,fa,fa,9,1,ke,1)
c            iet=ietold
c.... equivalent to  -xsi-  eq. 5.120  Vertieferarbeit Eidel
czr300398   do ierror = 1,numerr
czr300398              e_bar(ierror) = 
czr300398     +            eval/100.*sqrt((u_om(ierror)+e_om(ierror))/numel)
czr300398    enddo
c.... check type of error desired as refinement indicator
            ictnam = ct(2,l)
c
            ierror =  min(max(1,iet(2)),numerr )
            write(*,3003)   e_name(ierror),eval
            write(iow,3003) e_name(ierror),eval
            call matcop(e_ome(1,ierror),numel,1,erro)
          else
c.... Compute nodal stresses before remeshing
            write(*,3001)
            write(iow,3001)
            goto 998  
          endif
        elseif (pcomp(lct(l),fixd,4)) then
c         write(*,'(A)')' - using fixed error'
c	    iet    = 3
c.... save information on file (ioerr = 1)
c	    ioerr = 0
c  	    hflgu = .false.
c.... loop over elements
c         ietold=iet       ! avoid printing of error
c         iet = 4
c         call formfe(u,dr,dr,dr,fa,fa,fa,fa,9,1,ke,1)
c         iet=ietold
czr       call matcop(erroe,numel,1,erro)
c         call matcop(e_ome(1,1),numel,1,erro)
        elseif (pcomp(lct(l),'frac',4)) then
          if(fl(11)) then
c            write(*,'(A)')' - fracture as indicator'
c.... save information on file (io = 1)
c            ioerr = 0
c            hflgu = .false.
cc.... loop over elements
c            iet = 6
c            call formfe(u,dr,dr,dr,fa,fa,fa,fa,9,1,ke,1)
c            call matcop(e_ome(3,1),numel,1,erro)     
c            write(*,*)'not yet tested   -> reme_nl1.f' 
          else
c.... compute nodal stresses before remeshing
            write(*,3001)
            write(iow,3001)
            goto 998  
          endif
        endif
      endif
c......................................................set error +1/-1
c.... loop over elements  (sets error +1 / -1)
c     write(*,'(A)')' - set error'
      do 121 i=1,ke
c....
czr        if(nprior.eq.1)then                                  
czr          if(erro(i).gt.tolp)then
czr            erro(i)=1.d0
czr          else
czr           erro(i)=-1.d0
czr          endif
czr        elseif(naiter.eq.1) then
      if(naiter.eq.1) then
        erro(i)=1.d0
czr        elseif(npoll0.eq.1) then
czr           call polout(i,nen1,ndm,iek,x,nvx)
czr           if(nvx.eq.1) then
czr            erro(i)=1.d0
czr           elseif(nvx.eq.0)then
czr            erro(i)=-1.d0
czr           else
czr            stop
czr           endif
      elseif(pcomp(lct(l),unif,4)) then
        erro(i) = 1.d0
      elseif(pcomp(lct(l),adap,4)
     +   .or.pcomp(lct(l),fixd,4)
     +   .or.pcomp(lct(l),frac,4))then
c.... set relative error in element
c.... equivalent to  -xsi-  eq. 5.121  Vertieferarbeit Eidel
czr     erro(i) = sqrt(erro(i))/e_bartemp 
        if(erro(i).gt.1.d0)then
          erro(i)=1.d0
        else
          erro(i)=-1.d0
        endif
      endif
121   continue                            ! enddo first loop over elements
c.... tieed nodes for mesh generation
      if(ndf.eq.6.and.naiter.gt.1
     +      .and.(pcomp(lct(l),fixd,4).or.pcomp(lct(l),fixd,4)))then
         call tiegen (x,iek0,ikz0,ike0,nen1,ndm,erro)
         write(*,'(A)')'ndf.eq.6.and.naiter.gt.1    -> reme_nl1.f'
      endif
c....
      fmax=0.d0
C.... loop over elements..........................refine error elements
c     write(*,'(A)')' - refine error elements'
      i=0
122   i=i+1
c.... check if 2 or 4 node element
        l_four = .false.
        if (iek(nen-1,i).ne.0) l_four = .true.
c
czr     if(erro(i).gt.0.d0)then
        if(erro(i).gt.0.d0 .and. l_four)then
c.... convert y elements to x (K)
          if(iael0(i).lt.0)then            !... convert Y-elements to x (K)
            nst=1
            call xymod_nl
     +         (nst,i,u0,u,hst0,x0,x,iek0,ikz0,ike0,
     +          mikno,nen1,nhist,ndm,ndf,
     +          ianp,ianp0,iael,iael0,naiter, erro, kf )
          else                             !... generate x(k) element
c.... regular refinement (subdivide with a cross)
            iael0(i)=naiter
            call xgen_nl (0,i,u0,u,hst0,x0,x,iek0,ikz0,ike0,
     +            mikno,nen1,nhist,ndm,ndf,
     +            ianp,ianp0,iael,iael0,naiter,erro,    erron0)
          endif
        else
          call ogenes (0,i,iek0,nen1)
        endif
      continue                     ! enddo .. second loop over elements
      if(i.lt.ke0) goto122         ! end   counter modification
c.... 
czr   if (pcomp(flag,adaf)) ke0=ke
czr   if (pcomp(flag,adaf)) kk0=kk
c................................................generate compatibility
      i=0
131   continue              ! not possible with do (counter modified)
      i=i+1
czr     if(naiter.eq.1.and.i.lt.nadap)then    !nadap (n1 feap 2. Zeile) 
czr        iek0(nen1,i)=abs(iek0(nen1,i))*(-1.d0)
czr        do 117 il =ke0+(i-1)*3+1,ke0+(i-1)*3+3
czr117     iek0(nen1,il)=abs(iek0(nen1,il))*(-1.d0)
czr     endif
c.... check if element is marked as 'neighboring' element
        if( iek0(nen1,i).gt.0) goto 130
c.... switch 2/4 node element 
        if (iek(nen-1,i).ne.0) then    ! case of 4 node element
c
          k=0
          if(naiter.eq.1) goto 40
c.... check for marked nodes   (iek lt 0 -> 'new' edge node)
          do j=1,4
            if(ianp0(abs(iek0(j,i))).lt.0)k=k+1
          enddo
40        continue
c.... no marked node / x (regular) element
        if(k.eq.0.and.iael0(i).ge.0.and.naiter.gt.1)then
          call xgen_nl (1,i,u0,u,hst0,x0,x,iek0,ikz0,ike0,
     +         mikno,nen1,nhist,ndm,ndf,
     1         ianp,ianp0,iael,iael0,naiter,erro,    erron0)
czr       elseif( k.eq.1 ) then
        elseif( k.eq.1 .or. iael0(i).lt.0 ) then
          if(iael0(i).lt.0) then
c.... one marked node / y element
            nst=0
            call xymod_nl
     +        (nst,i,u0,u,hst0,x0,x,iek0,ikz0,ike0,
     +         mikno,nen1,nhist,ndm,ndf,
     +         ianp,ianp0,iael,iael0,naiter,erro, kf)
            elseif( iael0(i).ge.0 ) then   
c.... one marked node / x element
c.... generate y element
              call ygen_nl(0,i,u0,u,hst0,   x0,x,iek0,ikz0,ike0,
     +            mikno,nen1,nhist,ndm,ndf,
     +            ianp,ianp0,iael,iael0,naiter,erro)
            endif
c.... two marked nodes / no action required
          elseif (k.eq.2) then
            call ogenes (1,i,iek0,nen1)
          endif
        else                                   ! case of 2 node element
              nst = 0
c1              call l_gener
c1czr  +        (nst,i,u0,u,hst0,x0,x,iek0,ikz0,ike0,
c1     +        (nst,i,u0,u,hst0,x0,x,iek0,ikz0,ike0,ikz,ike,
c1     +         mikno,nen1,nhist,ndm,ndf,
c1     +         ianp,ianp0,iael,iael0,naiter,erro, kf)
          call l_gener 
     +         (1,i,u0,u,hst0,x0,x,iek0,ikz0,ike0,
     +         mikno,nen1,nhist,ndm,ndf,
     1         ianp,ianp0,iael,iael0,naiter,erro,    erron0)
        endif
130   continue
      if(i.lt.ke) goto 131       ! end   ...   third loop over elements
c....
czr      if(fsteu.eq.1.and.kf.lt.ct2.and.naiter.ne.1)then
czr        kf=kf+1
czr        write(*,*) 'adaf kf=',kf,'fmax=',fmax
czr       if(pcomp(flag,adaf)) ke0=ke
czr       if(pcomp(flag,adaf)) kk0=kk
czr        goto 100
czr      endif
998   continue
3001  format('**WARNING** Compute nodal stresses before remeshing')
3002  format(I2,1x,'iteration')
3003  format(2x, a15 ,f5.1,' %')
      end                                                  
c......................................................... end.SR.admess
