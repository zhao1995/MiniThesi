      subroutine mreme (d,id,ix,u,dr,lct,ct,ndm,ndf,nen1,llreme)
c --------------------------------------------------------------------
c.... macro instruction subprogram for adaptivity
c.... including switch linear nonlinear analysis
c
c     WW BS KIT 04/15 korrekte Einbindung in PMACR  
c     WW BS KIT 04/15 soweit möglich ALLOCATE 
c                     Rest geht nicht wg. diverser Feldüberschreitungen  
c     WW BS KIT 04/15 etwas aufgeräumt 
c
c--------------------------------------------------------------------
      USE adap
      USE cdata
      USE comfil
      USE errin1
      USE errin2
      USE errin3
      USE hdatam
      USE iofile
      USE ldata
      USE mdata
      USE plong
      USE psize
      USE rdata
      USE rfeap
      USE strnam
      USE doalloc
      implicit double precision (a-h,o-z)
      logical pcomp,lirm,fa,tr
      character*4 lct(*)
      character*229 fint,fint1,fint2
      character*80 yyy
      integer id(*)
      real*8  ct(3,*)
      real*8 u(*),dr(*)
      dimension d(*)
      common /adap1/ naiter,mnph,mnphdt,nhist,nhsw
      common /gener/ kk,km,ke
czr   common /error1/ nerrs,nerre,nerrf
      common /rmsh1/ nx,melem  ! ### delete when adap ok    

      common m(maxm)
      data fa,tr/.false.,.true./

      if (pcomp(lct(l),'dman',4) ) then
c....   manipulate d-field to 'assist' convergence
        ichange = ct(1,l)
        write(*,731) ichange, ct(2,l)
731     format('Change d: position',i2,'value',f8.3)
        d(ichange) = ct(2,l)
        return

      else if (pcomp(lct(l),'wdfm',4) ) then
c....   write deformed mesh to file
        call write_defm(coor,u,ndm)
        return

      else if (pcomp(lct(l),'rest',4) ) then
        call pseta(mnph,numnp*nhist,ipr,tflb,'REST-1')    ! --> hst field
        call pseta(mnph0,numnp*nhist,ipr,tflb,'REST-2')   ! --> hst field
        nhsw = 2                                 ! proj. np->gp
        call rest_adap(fres,u,m(mnph0),ix,
     +      ndm,ndf,nen1,nhist,1,1,1.0d0)
        call rearrange(m(mnph),m(mnph0),nhist,nhsw)
        hflgu  = .true.
        h3flgu = .true.
        call formfe(u,dr,dr,dr,fa,fa,fa,fa,23,1,numel,1)
        return
      end if

      kk=0
      km=0
      ke=0
      lirm = .false.
c.... transfer to correct process
czr   ne = 0
c.... remesh FE-mesh
c     write(iow,*) '[reme,(adap,unif,new,old),n1,n2]'

c.... open new input file
      fint1 = finp  ! nfile 
      fint2 = finp  ! nfile.mrm
      call dochar2(finp,ipos)
      if(finp(ipos:ipos).eq.'n') then
        lirm = .true.
        call addext(fint1,'tmp ')
        call addext(fint2,'mrm ')
      else
         call dochar1(fint1,'n',ipos) ! ifile->nfile
         fint2 = fint1
         call addext(fint2,'mrm ')
      end if
      inpf = abs(ior) + abs(iow) ! ww##### 31! to check in combination with other files
c....
      open(inpf,file=fint1,form='formatted',status='unknown') ! nfile

c.... read new input file (nfile)
      if (pcomp(lct(l),'new',3) ) then
        call dochar2(finp,ipos)
        if(finp(ipos:ipos).eq.'n') then
          close(inpf,status='delete')
        else
c....
          close(inpf)
          close(abs(ior))
          open(abs(ior),file=fint1,form='formatted',status='unknown') ! nfile
          finp = fint1
        end if
        rewind(iabs(ior))
c....   clear and start
c       kmax = 0 !??
        call deallocall
        llreme=-2  ! exit pmacr  
        irfeap=3   ! restart for REME
        return
      end if

c.... read old input file (nfile)
      if (pcomp(lct(l),'old',3) ) then
        call dochar2(finp,ipos)
        if(finp(ipos:ipos).eq.'n') then
          close(inpf,status='delete')
        end if
        rewind(iabs(ior))
c....   clear and start
c       kmax = 0
        call deallocall
        llreme=-2  ! exit pmacr  
        irfeap=3   ! restart for REME
        return
      end if

c.... REME, ADAP etc 
      naiter = 0
      nt = 2
cww   nfdel  = 4**nt*2*numnp
      mnumnp = (4**nt-1)*2*numel + numnp 
c....
cww   kmaxold = kmax                          ! reset m

      call pseta(   n9a,numel*(nen+1),   1,tflb,'RMESH-ike   ')   
      call pseta(   n9b,numel*(nen+1),   1,tflb,'RMESH-ikz   ')   
      call pseta(   n9c,numnp,           1,tflb,'RMESH-ianp  ')   
      call pseta(   n9d,numel*4,         1,tflb,'RMESH-iael  ')   
      call pseta( melem,mnumnp,          1,tflb,'RMESH-iek0  ')   
      call pseta(mnodel,mnumnp*nen1,     1,tflb,'RMESH-ike0  ')   
      call pseta(mnodez,mnumnp*nen1,     1,tflb,'RMESH-ikz0  ')   
      call pseta(   mel,mnumnp,          1,tflb,'RMESH-iael0 ')   
      call pseta(   mnp,mnumnp,          1,tflb,'RMESH-ianp0 ')   
      call pseta(  mmik,mnumnp,          1,tflb,'RMESH-mikno ')   
      call pseta( mcoor,numnp+1,       ipr,tflb,'RMESH-x0    ')   
      call pseta( merro,mnumnp*ndm,    ipr,tflb,'RMESH-erron0')   
                                                                 
cww   geht nicht wg. diverser Feldüberschreitungen 
cww   call ialloc(adaike,   numel*(nen+1),'REMESH-ike', tflb)         ! n9a       
cww   call ialloc(adaikz,   numel*(nen+1),'REMESH-ikz', tflb)         ! n9b
cww   call ialloc(adaianp,  numnp,        'REMESH-ianp',tflb)         ! n9c
cww   call ialloc(adaiael,  numel*4,      'REMESH-iael',tflb)         ! n9d
cww   call ialloc(adaiek0,  mnumnp,       'RMESH-iek0 melem',tflb)    ! melem
cww   call ialloc(adaike0,  mnumnp*nen1,  'RMESH-ike0 mnodel',tflb)   ! mnodel
cww   call ialloc(adaikz0,  mnumnp*nen1,  'RMESH-ikz0 mnodez',tflb)   ! mnodez
cww   call ialloc(adaiael0, mnumnp,       'RMESH-iael0',tflb)         ! mel
cww   call ialloc(adaianp0, mnumnp,       'RMESH-ianp0',tflb)         ! mnp 
cww   call ialloc(adamikno, mnumnp,       'RMESH-mikno',tflb)         ! mmik 
cww   call ralloc(adax0,    numnp+1,      'RMESH-x0 mcoor',tflb)      ! mcoor   
cww   call ralloc(adaerron0,mnumnp*ndm,   'RMESH-erron0',tflb)        ! merro   

c.... determine node-element connectivities for adaptivity and others
C.... felder m(ne) und m(nnum) werden nur temp. benoetigt  allocate local
      ne   = kmax
      nnum = ne + numnp + mod(numnp,ipr)
      call elnode (econ,m(n9a),m(n9b),nen1,numnp,numel,nmat,
     +             m(ne),m(nnum))    

c.... mesh refinement
      if  ((pcomp(lct(l),'unif',4)).or.(pcomp(lct(l),'adap',4))
     + .or.(pcomp(lct(l),'fixd',4)).or.(pcomp(lct(l),'frac',4)) ) then
        iet(1) = 3        ! switch in element for 'macro>reme,' mode
        iet(2) = ct(2,l)  ! switch in element for error type

cww     call pseta (nerrt,      numel, ipr,tflb,'REMESH-errot/tem')   ! --> erro/temp  
cww     call pseta( nerr,numerr*numel, ipr,tflb,'REMESH-erron/new')   ! --> erro/new!

        call ralloc(adaerrot,    numel,'REMESH-errot/tem',tflb)        ! --> erro/new!
        call ralloc(e_ome,numerr*numel,'REMESH-erron/new',tflb)       ! --> erro/new!

c       call pseta(nerrs,  numel, ipr,tflb,'REMESH-erros   ')   ! --> erros
c       call pseta(nerre,  numel, ipr,tflb,'REMESH-erroe   ')   ! --> erroe
c       call pseta(nerrf,  numel, ipr,tflb,'REMESH-frac_err')   ! --> frac_err

c.... linear/nonlinear analysis ww: only linear!
cww     if (linear) then
          call admess(u,dr,econ,coor,m(n9a),m(n9b),m(n9c),m(n9d),
     1      m(mcoor),m(melem),m(mnodel),m(mnodez),m(mel),m(mnp),m(mmik),
     2      m(merro),e_ome,adaerrot,nen1,ndm,ndf,lct,ct,fint2)

czr  2      adaerron0,m(nerrs),m(nerre)         ,adaerrot,
czr  2      adaerron0,m(nerrs),m(nerre),m(nerrf),adaerrot,

cww     else
cww       call pseta(mnoddisp,5*numel*ndf,ipr,tflb)    ! --> u0
cwwcpl>
cww       call pseta(mnph,numnp*nhist,ipr,tflb)         ! hist.var. at np's
cww       call pseta(mnphdt,numnp*nhist,ipr,tflb)       ! nodal averaging
cww       call pzero (m(mnph),numnp*nhist)
cww       call pzero (m(mnphdt),numnp*nhist)
cww       nhsw = 1                                      ! proj. gp->np
cww       hflgu  = .false.
cww       h3flgu = .false.
cww       call formfe(u,dr,dr,dr,
cww     +    .false.,.false.,.false.,.false.,19,1,numel,1)
cww       call haverage(m(mnphdt),m(mnph),numnp)
cww       call pseta(mnph0,9*numel*nhist,ipr,tflb)       ! hist.var. at np's
cww       call pzero (m(mnph0),9*numel*nhist)
cww       call rearrange(m(mnph),m(mnph0),nhist,1)
cwwcpl<
cww       call reme_nl(u,dr,
cww     +      m(n9),m(n8),adaike,adaikz,adaianp,adaiael,
cww     +      adax0,adaiek0,adaike0,adaikz0,adaiael0,adaianp0,adamikno,
cwwczr  +      adaerron0,m(nerrs),m(nerre),m(nerrf),adaerrot,
cww     +      adaerron0,e_ome, adaerrot,
cwwcnl>
cww     +      m(mnoddisp),m(mnph0),
cww     +      nen1,ndm,ndf,lct,ct,fint2)
cww     end if
      else
        write(yyy,*)'Specify REME options  [adap/unif/new,old],n1,n2'
        call drawmess(yyy,1,0)
        return
      end if

c.... calculate boundary conditions for new nodes
      call pseta(nbou , kk*ndf, 1, tflb,'REMESH-newbc')   ! --> new bc's
cww   call ialloc(adanewbc,kk*ndf,'REMESH-newbc',tflb)    ! --> new bc's

      call calboun(psid,m(mcoor),m(nbou),ndm,ndf)

c.... set pointers for new - connected - arrays
cww   call pseta(nx ,kk*ndm,  ipr,tflb,'REMESH-nx')       ! --> x
cww   call pseta(nf ,kk*ndf,  ipr,tflb,'REMESH-nf')       ! --> f
cww   call pseta(nu ,kk*ndf,  ipr,tflb,'REMESH-nu')       ! --> f
cww   call pzero(m(nf),kk*ndf)
cww   call pzero(m(nu),kk*ndf)

      call ralloc(adanx,kk*ndm,'REMESH-nx',tflb)          ! --> x
      call ralloc(adanf,kk*ndf,'REMESH-nf',tflb)          ! --> f
      call ralloc(adanu,kk*ndf,'REMESH-nu',tflb)          ! --> f

c.... connect arrays
      call conarray (coor,m(mcoor),adanx,numnp*ndm,(kk-numnp)*ndm)

c.... calculate new surface loads
      call calload(adanx,m(melem),m(mnodel),m(mnodez),adanf,ndm,ndf)

c.... write new FEAP input file (- feap,coor,elem,boun,forc)
      call winput(coor,m(mcoor),m(melem),m(nbou),adanf,
     +                ndm,nen1,ndf,inpf,finp)

c.... write restart (nodal/element markers)
      fint = fsav

c.....set parameter ascii/binary for restart-file: ibin=0=binary/1=ascii
      ibin=0
      call rmstrt(fint2,ndm,ndf,nen1,2,ibin,
     1   kk,ke,naiter,econ,m(melem),m(n9d),m(mel),m(n9c),m(mnp),ct)

c.... write restart 'rinput' (nodal displacements, hist.var.,...)
cww   if  (.not.linear) then
cww     call conarray (u,m(mnoddisp),adanu,numnp*ndf,(kk-numnp)*ndf)
cww     numeltmp = numel
cww     numnptmp = numnp
cww     numel = ke
cww     numnp = kk
cww     call rest_adap(fsav,adanu,m(mnph0),adaiek0,
cww     +      ndm,ndf,nen1,nhist,1,2,1.0d0)
cww     numel = numeltmp
cww     numnp = numnptmp
cww   end if

c.... info old/new nodes/elements

      write(*,1001) numel, ke
      write(*,1002) numnp, kk

c....
      if(lirm) then                          ! file with .irm exists
        rewind (inpf)
        rewind (iabs(ior))
20      continue                             ! do forever
          read(inpf,1003,err=25,end=25) yyy
          write(iabs(ior),1003) yyy
        goto 20
25      continue
        rewind (inpf)
        write(inpf,'(A)') '...temp...'
        close (inpf,status='delete')
        return
      end if

      close (inpf,status='keep')
      return
1001  format(' elements : was',I5,' now',I5)
1002  format(' nodes    : was',I5,' now',I5)
1003  format(A80)
      end
c
      subroutine elem_ener(a,b,dv,numsig)
c----------------------------------------------------------------------
c.... compute element energy (with stresses in gauss-points)
c----------------------------------------------------------------------
      USE errin1
      USE errin2
      USE iofile
      implicit double precision (a-h,o-z)
      dimension a(10),b(10)
      do i = 1,numsig
        u_om(1) = u_om(1) + a(i) * b(i)*dv  !.."energy", E-normed
        u_om(2) = u_om(2) + a(i) * a(i)*dv  !.."energy", L_2-normed
      enddo
      end
c
      subroutine winput(x,x0,iek0,nboue,f,
     +      ndm,nen1  ,ndf,inpf,finp)
c----------------------------------------------------------------------
c.... writes input file - input/output
c.... modified for the needs of adaptiv mesh refinement ( <-- restrt )
c----------------------------------------------------------------------
      USE cdat1
      USE cdata
      USE iofile
      implicit double precision (a-h,o-z)
      logical pcomp
      character*80 yyy
      character*229 finp
      dimension x(ndm,*),iek0(nen1,*),x0(ndm,*)
      dimension nboue(ndf,*)
      dimension  f(ndf,*)
      common /gener/  kk,km,ke
      write(iow,'(A)')'write new I N P U T - file  --> n[file] '
c....
      nior = abs(ior)

cww   rewind nior  ! old, now file has to be opened ww 20.01.07

      open(nior,file=finp,form='formatted',status='unknown')

10    continue              ! external loop over 'old' input file
        read(nior,1003,end=501) yyy
c....   find/write  feap
        if( pcomp(yyy,'feap',4) )then
          write(inpf,1003) yyy
c....     vgl. feaps3.f
          iortemp = ior
          ior = nior
          call pintio(yyy,8)
          ior = iortemp
          read(yyy,1006,err=600,end=700)
     +         n1,   n2,   n3,    n4, n5, n6, n7, n8
c....          numnp,numel,nummat,ndm,ndf,nen,nad,ndd
700       continue
          write(inpf,1002) kk,ke,nummat,ndm,ndf,nen,n7,ndd
          write(inpf,*)'    '
          read(nior,1003,end=501) yyy

c...    write coordinates to file
        else if( pcomp(yyy,'coor',4)) then
          write(inpf,1003) yyy
210       continue
          read(nior,1003,end=501) yyy
          if ( yyy(1:10) .eq. '          ') goto 220
          goto 210
220       continue

          do inx = 1,kk
            if (inx.le.numnp) then
              write(inpf,1000) inx,(x(i,inx),i=1,ndm)
            else
              write(inpf,1000) inx,(x0(i,inx-numnp),i=1,ndm)
            end if
          end do
          write(inpf,1003) yyy

c....   write new element-nodal-connectivities to file
        else if( pcomp(yyy,'elem',4) ) then
          write(inpf,1003) yyy

310       continue
          read(nior,1003,end=501) yyy
          if ( yyy(1:10) .eq. '          ') goto 320
          goto 310
320       continue

          do ien = 1,ke
            if (iek0(nen1,3).ne.0)then
              write(inpf,1001)
     +            ien,abs(iek0(nen1,ien)),(iek0(i,ien),i=1,nen)
            else                                  ! 2 node line element
              write(inpf,1001)
     +            ien,abs(iek0(nen1,ien)),(iek0(i,ien),i=1,2)
            end if
          end do
          write(inpf,1003) yyy

c....   ignore 'bloc'  (not used in new input file)
        else if( pcomp(yyy,'bloc',4) ) then
410       continue
          read(nior,1003,end=501) yyy
          if ( yyy(1:10) .eq. '          ') goto 420
          goto 410
420       continue

c....   write new boundary conditions
        else if( pcomp(yyy,'boun',4) ) then
          write(inpf,1003) yyy

510       continue
          read(nior,1003,end=501) yyy
          if ( yyy(1:10) .eq. '          ') goto 520
          goto 510
520       continue

          do ni=1,kk
            nbsum = 0
            do nf=1,ndf
              nbsum = nbsum+nboue(nf,ni)
            end do
            if (nbsum.ne.0) then
              write(inpf,1005)  ni , (nboue(in,ni),in=1,ndf)
            end if
          end do
          write(inpf,1003) yyy

c....   write forces for input file
c        else if( pcomp(yyy,'forc',4).or. pcomp(yyy,'load',4) ) then ! geht nicht?? load funktioniert aber, solange nicht mehr Knoten belegt werden  WW??
        else if( pcomp(yyy,'forc',4) ) then
          write(inpf,1003) yyy
c
610       continue
          read(nior,1003,end=501) yyy
          if ( yyy(1:10) .eq. '          ') goto 620
          goto 610
620       continue
c
          do ikn = 1,kk
            if(ndf.eq.2) then
              if(abs(f(1,ikn)).gt.0.d0.or.abs(f(2,ikn)).gt.0.0d0)then
                write(inpf,1004) ikn, (f(ifhg,ikn), ifhg=1,ndf)
              end if
            else if(ndf.eq.3) then
              if(abs(f(1,ikn)).gt.0.d0.or.abs(f(2,ikn)).gt.0.d0
     +                                .or.abs(f(3,ikn)).gt.0.0d0)then
                write(inpf,1004) ikn, (f(ifhg,ikn), ifhg=1,ndf)
              end if
            else if(ndf.eq.5) then
              if(abs(f(1,ikn)).gt.0.d0.or.abs(f(2,ikn)).gt.0.d0
     +           .or.abs(f(3,ikn)).gt.0.0d0.or.abs(f(4,ikn)).gt.0.0d0
     +                                .or.abs(f(5,ikn)).gt.0.0d0)then
                write(inpf,1004) ikn, (f(ifhg,ikn), ifhg=1,ndf)
              end if
            end if
          end do
          write(inpf,1003) yyy

c..  .. ignore 'pola','tie'   (singles)  (not used in new input file)
        else if( pcomp(yyy,'pola',4).or. pcomp(yyy,'tie',3)) then
          continue

        else
          write(inpf,1003) yyy
        end if
      goto 10                ! end 'external' loop
501   continue
      return
600   write(*,*) 'error in winput!'
1000  format(I6, ',0,' , 3(e14.7,','))
1001  format(I6,',',I3,',', 4(I6,',') )
1002  format(2(I6,','),6(I4,',') )
1003  format(A80)
1004  format ( I4, ',0,', 7(g12.6 ,','))
1005  format ( I4, ',0,', 6( I2 ,','))
1006  format(10i8)
      end
c....................................................... end.SR.wrinput
      subroutine write_defm(x,u,ndm)
c----------------------------------------------------------------------
c.... writes input file - input/output
c.... modified for the needs of adaptiv mesh refinement ( <-- restrt )
c----------------------------------------------------------------------
      USE cdat1
      USE cdata
      USE comfil
      USE iofile
      implicit double precision (a-h,o-z)
      logical pcomp
      character*80 yyy
      character*229 fint
      dimension x(ndm,*),u(ndm,*)
      common /gener/  kk,km,ke
c.... open new input file
      fint = finp
      call dochar2(finp,ipos)
      call dochar1(fint,'f',ipos)
      inpf = abs(ior) + abs(iow)
      open(inpf,file=fint,form='formatted',status='unknown')
      write(iow,'(A)')'write new I N P U T - file  --> f[file] '
c....
      nior = abs(ior)
      rewind nior
10    continue              ! external loop over 'old' input file
        read(nior,1003,end=501,err=600) yyy
c.2.. write coordinates to file
        if( pcomp(yyy,'coor',4)) then
          write(inpf,1003) yyy
c ......................................................
            inx = 0
210         continue
              inx = inx + 1
              read(nior,1003,end=501) yyy
c              write(inpf,1000) inx,(x(i,inx)+u(i,inx),i=1,ndm)
              if ( yyy(1:10) .eq. '          ') goto 220
            goto 210
220         continue
c
          do inx = 1,numnp
            write(inpf,1000) inx,(x(i,inx)+u(i,inx),i=1,ndm)
          enddo
          write(inpf,1003) yyy
c ......................................................
        else
          write(inpf,1003) yyy
        endif
      goto 10                ! end 'eternal' loop
501   continue
      close(inpf)
      return
600   write(*,*) 'error in winput!'
1000  format(I6, ',0,' , 3(e14.7,','))
1003  format(A80)
      end
c................................................... end.SR.write_defm
c............................................................ SR.rmstrt
      subroutine rmstrt(fres,ndm,ndf,nen1,isw,iasc,
     1      kk,ke,naiter,iek,iek0,iael,iael0,ianp,ianp0, ct)
c----------------------------------------------------------------------
c.... restart files - input/output  in ascii for iasc=1
c     modified for the needs of adaptiv mesh refinement (orig. SR restrt)
c     write nodal and element markings
c----------------------------------------------------------------------
      USE cdata
      USE errin1
      USE errin2
      USE fdata
      USE hdatam
      USE iodata
      USE iofile
      USE ldata
      implicit double precision (a-h,o-z)
      logical exst
      character fres*229,yorn*1,y*80
      real*8 ct(3,*)
      dimension iek(nen1,*),iek0(nen1,*),
     +          iael(*),iael0(*),ianp(*),ianp0(*)
C
1     inquire(file=fres,exist=exst)
      if(.not.exst.and.isw.eq.1) then
        write(iow,3002) fres
        if(ior.lt.0) then
          write(*,3002) fres
          write(*,3003)
10        read (*,1000,err=11,end=12) yorn
          goto  13
11        call  errclr ('RESTRT')
          goto  10
12        call  endclr ('RESTRT',yorn)
13        if(yorn.eq.'y' .or. yorn.eq.'Y') then
             write(*,3004)
20           read (*,1001,err=21,end=22) fres
             goto  1
21           call  errclr ('RESTRT')
             goto  20
22           call  endclr ('RESTRT',fres)
             go to 1
           end if
        end if
        return
      end if
c.... open file
      if(iasc.eq.0) then
        if(exst) then
          open(ios,file=fres,form='unformatted',status='old',err=23)
        else
          open(ios,file=fres,form='unformatted',status='new')
        end if
      else if(iasc.ne.0) then
        if(exst) then
          open(ios,file=fres,form='formatted'  ,status='old',err=23)
        else
          open(ios,file=fres,form='formatted'  ,status='new')
        end if
      end if
      goto 24
23    write(*,25)
      return
25    format(' Error, the restart/save file is not un-/formatted!')
24    rewind ios
c....
      if(iasc.eq.0) then ! unformatted files
c....   read restart files
        if(isw.eq.1) then
c....     general values, element-node connectivities
          read(ios) nnpo,nnlo,nnmo,ndmo,ndfo,naitero
          if((nnpo.eq.numnp).and.(nnlo.eq.numel).and.(nnmo.eq.nummat)
     1         .and.(ndmo.eq.ndm).and.(ndfo.eq.ndf)) then
            write(*,*) 'found matching markings'
            if (naiter.lt.naitero) naiter = naitero
            call matcoi(iek,nen1,numel,iek0)
 
c....       element identification
            call readfiu(iael0,1,numel,ios)
c....       nodal identification
            call readfiu(ianp0,1,numnp,ios)
C....       case restart not identical
          else
            call matcoi(ianp,numnp,1,ianp0)
            call matcoi(iael,numel,1,iael0)
            call matcoi(iek,nen1,numel,iek0)
          end if
        elseif(isw.eq.2) then
c....     save information for restart during mesh refinement
c....     general values  -> FEAP
          write(ios) kk,ke,nummat,ndm,ndf,naiter
c....     element - nodal connectivities
c....     element identification
          call writefiu(iael0,1,ke,ios)
c....     nodal identification
          call writefiu(ianp0,1,kk,ios)
        end if
C
c...................................................................
      elseif(iasc.ne.0) then ! formatted files
c....   read restart files
        if(isw.eq.1) then
c....     general values, element-node connectivities
          read(ios,4020) y
          read(ios,4002) nnpo,nnlo,nnmo,ndmo,ndfo,naitero
          if((nnpo.eq.numnp).and.(nnlo.eq.numel).and.(nnmo.eq.nummat)
     1         .and.(ndmo.eq.ndm).and.(ndfo.eq.ndf)) then
            write(*,*) 'found matching markings'
            if (naiter.lt.naitero) naiter = naitero
            call matcoi(iek,nen1,numel,iek0)
c....       element identification
            read(ios,4020) y
            call readfi(iael0,1,numel,ios)
c....       nodal identification
            read(ios,4020) y
            call readfi(ianp0,1,numnp,ios)
C....       case restart not identical
          else
            call matcoi(ianp,numnp,1,ianp0)
            call matcoi(iael,numel,1,iael0)
            call matcoi(iek,nen1,numel,iek0)
          end if
        elseif(isw.eq.2) then
c....     save information for restart during mesh refinement
c....     general values  -> FEAP
          write(ios,*) ' numnp,numel,nummat,ndm,ndf,naiter '
          write(ios,4002) kk,ke,nummat,ndm,ndf,naiter
c....     element - nodal connectivities
c....     element identification
          write(ios,4004)
          call writefi(iael0,1,ke,ios)  !iael: - new y element
c....     nodal identification
          write(ios,4005)
          call writefi(ianp0,1,kk,ios)  !ianp: - new edge / + new center node
        end if
      end if
c.... close the file
      close(ios)
      return
1000  format(a1)
cww1001  format(a17)
cww3002  format(' **ERROR** Restart file ',a17,' does not exist')
1001  format(a229)
3002  format(' **ERROR** Restart file ',a229,/,' does not exist')
3003  format(11x,'Specify new name for restart file? (y or n) >',$)
3004  format(11x,'New Restart File Name >',$)
4002  format(6i7,l5)
czr4003  format(' iek0 - field')
4004  format(' iael0 - field')
4005  format(' ianp0 - field')
4020  format(a80)
      end
c.......................................................... end.SR.rmstrt
c............................................................. SR.writefi
      subroutine writefi(inf,n1,n2,ios)
c----------------------------------------------------------------------
c.....write an integer field on save file
c----------------------------------------------------------------------
      integer inf(n1,n2)
      do i = 1,n2
        write(ios,1000) (inf(k,i),k=1,n1)
      end do
1000  format(1x,8I12,/,1x,8I12)
      return
      end
c......................................................... end.SR.writefi
c.............................................................. SR.readfi
      subroutine readfi(inf,n1,n2,ios)
c----------------------------------------------------------------------
c.....read an integer field from restart file
c----------------------------------------------------------------------
      integer inf(n1,n2)
      do i = 1,n2
        read(ios,1000) (inf(k,i),k=1,n1)
      end do
1000  format(1x,8I12,/,1x,8I12)
      return
      end
c.......................................................... end.SR.readfi
c............................................................. SR.writefiu
      subroutine writefiu(inf,n1,n2,ios)
c----------------------------------------------------------------------
c.....write an integer field on save file
c----------------------------------------------------------------------
      integer inf(n1,n2)
      do i = 1,n2
        write(ios) (inf(k,i),k=1,n1)
      end do
      end
c......................................................... end.SR.writefiu
c.............................................................. SR.readfiu
      subroutine readfiu(inf,n1,n2,ios)
c----------------------------------------------------------------------
c.....read an integer field from restart file
c----------------------------------------------------------------------
      integer inf(n1,n2)
      do i = 1,n2
        read(ios) (inf(k,i),k=1,n1)
      end do
      end
c.......................................................... end.SR.readfiu
c...............................................................SR.tiegen
      subroutine tiegen  (x,iek,ikz,ike,nen1,ndm,erro)
c----------------------------------------------------------------------
c.... unterteilt die Elemente an Schalenverschneidungen
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /gener/ kk,km,ke
      dimension x(ndm,*),iek(nen1,*)
      dimension erro(*),ikz(*),ike(*)
      dimension x1(3),x2(3),x3(3),x4(3)
      tol=0.0001d0
      do 100 i=1,ke
        if(erro(i).gt.0.d0) then
          do 200 j=1,4
            ik1=abs(iek(j,i))
            if(j.eq.4)then
              ik2=abs(iek(1,i))
            else
              ik2=abs(iek(j+1,i))
            end if
            call kaelem(ik1,ik2,i2,i,ikz,ike)
            if(i2.ne.0)goto200
c.... second loop elements to refine
            do 300 ii=1,ke
              if(erro(ii).lt.0.d0) then
                do 400 jj=1,4
                  iik1=abs(iek(jj,ii))
                  if(jj.eq.4)then
                    iik2=abs(iek(1,ii))
                  else
                    iik2=abs(iek(jj+1,ii))
                  end if
                  call kaelem(iik1,iik2,ii2,ii,ikz,ike)
                  if(ii2.ne.0)goto 400
c.... compare nodes
                  do 500 ndim =1,3
                    x1(ndim) =x(ndim,ik1)
                    x2(ndim) =x(ndim,ik2)
                    x3(ndim) =x(ndim,iik1)
                    x4(ndim) =x(ndim,iik2)
500               continue
           xl1=dsqrt((x1(1)-x3(1))**2+(x1(2)-x3(2))**2+(x1(3)-x3(3))**2)
           xl2=dsqrt((x1(1)-x4(1))**2+(x1(2)-x4(2))**2+(x1(3)-x4(3))**2)
           xl3=dsqrt((x2(1)-x3(1))**2+(x2(2)-x3(2))**2+(x2(3)-x3(3))**2)
           xl4=dsqrt((x2(1)-x4(1))**2+(x2(2)-x4(2))**2+(x2(3)-x4(3))**2)
                  if(xl1+xl4.lt.tol) then
                    erro(ii)=1.0
                    write(*,*) 'chance elmt',ii
                    goto 600
                  end if
                  if(xl2+xl3.lt.tol) then
                    erro(ii)=1.0
                    write(*,*) 'chance elmt',ii
                    goto 600
                  end if
400             continue
              end if
600         continue
300         continue
200       continue
        end if
100   continue
      return
      end
c......................................................... end.SR.tiegen
c...................................................................SR.admess
      subroutine admess(u,dr,iek,x,ike,ikz,ianp,iael,
     1      x0,iek0,ike0,ikz0,iael0,ianp0,mikno,
     2      erron0,e_ome,erro,
     3      nen1,ndm,ndf,lct,ct,fint2)
c ---------------------------------------------------------------------------
c.... adaptive mesh refinement for
c.... 4-node shell-element
c....       last modification :27.10.1992  baumann
c....       (adjustments to feap.ww : 08/95 ziegler)
c.... Verzweigung --> adaptiv/uniform/fixed
c ---------------------------------------------------------------------------
      USE cdata
      USE comfil
      USE errin1
      USE errin2
      USE errnam
      USE fdata
      USE iofile
      USE ldata
      implicit double precision (a-h,o-z)
      real*8  ct(3,*)
      real*8 u(*),dr(*)               ! --> SR formfe
      logical pcomp
      logical exst,hflgu
      logical fa                      ! false
      character*4  unif,adap,fixd,frac
      character*4 lct(*)
      character*229 fint2
      common /adap1/ naiter,mnph,mnphdt,nhist,nhsw
      common /gener/ kk,km,ke
      dimension x(ndm,*),x0(ndm,*),
     + iek(nen1,*),iek0(nen1,*),
     + iael(*),iael0(*),
     + ianp(*),ianp0(*),
     + erro(numel),e_ome(numel,2),
     + ikz(*),ikz0(*),
     + ike(*),ike0(*),
     + mikno(3,*),
     + erron0(*)

      data unif /'unif'/ , adap /'adap'/, fixd/'fixd'/, frac/'frac'/
      data fa /.false./
      
cww   call pzero (erro,numnp)                      ! org=??
      call pzero (erro,numel)                      !
      call pzero (e_ome,2*numel)                   !
      call pzeroi (iael,numel)                     !
      call pzeroi (ianp,numnp)                     !

c.....set parameter ascii/binary for restart-file: ibin=0=binary/1=ascii
      ibin=0

c.... read restart file
      inquire(file = fint2,exist=exst)
      if( exst ) then
        call rmstrt (fint2,ndm,ndf,nen1,1,ibin,
     1       kk,ke,naiter,iek,iek0,iael,iael0,ianp,ianp0,ct)
      else
        call matcoi(ianp,numnp,1,ianp0)
        call matcoi(iael,numel,1,iael0)
        call matcoi(iek,nen1,numel,iek0)
      end if
c.... initialisation
      kf=1
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
c.... generieren des grundnetzes
      continue
      if (naiter.gt.1) then
        if (pcomp(lct(l),adap,4)) then
c.... check stress errors
          if(fl(11)) then
c.... with respect to 'eval' percent of energy (default eval = 5 perc.)
            eval = ct(1,l)
            if(eval.lt.0.0001) eval = 5.
            do ierror = 1,numerr
              e_om(ierror) = 0.0
              e_bar(ierror)  = eval/100.*sqrt(u_om(ierror)/numel)
            end do
c....       save information on file (ioerr = 1)
            ioerr = 0
            hflgu = .false.
c....       loop over elements
            call formfe(u,dr,dr,dr,fa,fa,fa,fa,9,1,ke,1)
            ierror =  min(max(1,iet(2)),numerr )
            write(*,3003)   e_name(ierror),eval
            write(iow,3003) e_name(ierror),eval
            call matcop(e_ome(1,ierror),numel,1,erro)
c
          else
c....       Compute nodal stresses before remeshing
            write(*,3001)
            write(iow,3001)
            return
          end if
c        elseif (pcomp(lct(l),fixd,4)) then
ccfx>
cc          write(*,'(A)')' - using fixed error'
c     iet    = 3
cc.... save information on file (ioerr = 1)
c     ioerr = 0
c       hflgu = .false.
cc.... loop over elements
c          ietold=iet       ! avoid printing of error
c          iet = 4
c          call formfe(u,dr,dr,dr,fa,fa,fa,fa,9,1,ke,1)
c          iet=ietold
c          call matcop(erroe,numel,1,erro)
c        elseif (pcomp(lct(l),'frac',4)) then
ccfx>
c          if(fl(11)) then
cC           write(*,'(A)')' - fracture as indicator'
c.... save information on file (io = 1)
c            ioerr = 0
c            hflgu = .false.
cc.... loop over elements
c            iet = 6
c            call formfe(u,dr,dr,dr,fa,fa,fa,fa,9,1,ke,1)
c            call matcop(errof,numel,1,erro)
c          else
cc.... Compute nodal stresses before remeshing
c            write(*,3001)
c            write(iow,3001)
c            return
c          endif
        end if
      end if
c.... loop over elements  (sets error +1 / -1)
C      write(*,'(A)')' - set error'
      do 121 i=1,ke
        if(naiter.eq.1) then
          erro(i)=1.d0
        else if(pcomp(lct(l),unif,4)) then
          erro(i) = 1.d0
        else if(pcomp(lct(l),adap,4)
     +      .or.pcomp(lct(l),fixd,4)
     +      .or.pcomp(lct(l),frac,4))then
          if(erro(i).gt.1.d0)then
            erro(i)=1.d0
          else
            erro(i)=-1.d0
          end if
        end if
121   continue                              ! enddo first loop over elements
C
c.... tieed nodes for mesh generation
      if(ndf.eq.6.and.naiter.gt.1.and.pcomp(lct(l),fixd,4))then
        call tiegen (x,iek0,ikz0,ike0,nen1,ndm,erro)
        write(*,'(A)')' ndf.eq.6 -> not tested remesh1.f '
      end if
c....
      fmax=0.d0
C.... loop over elements  (refines elementes with error > 0)
c      write(*,'(A)')' - refine error elements'
czr      do 120 i=1,ke0                    ! counter will be modified
      i=0
122   i=i+1
        if(erro(i).gt.0.d0)then
c....     convert Y-elements to K
          if(iael0(i).lt.0)then
czr         write(2,*) 'kymods(error refinement)',i
            nst=1
            call kymods
     1         (nst,i,x0,x,iek0,ikz0,ike0  ,mikno,nen1,ndm,
     2           ianp,ianp0,iael,iael0,naiter,
     3           erro, kf )
          else
c....       regular refinement
            iael0(i)=naiter
czr         write(2,*) 'kgenes(error refinement)',i
            call kgenes (0,i,x0,x,iek0,ikz0,ike0,   mikno,nen1,ndm,
     1         ianp,ianp0,iael,iael0,naiter,erro,    erron0)
          end if
        else
          call ogenes (0,i,iek0,nen1)
        end if
czr120continue                     ! enddo .. second loop over elements
      continue                     ! enddo .. second loop over elements
      if(i.lt.ke0) goto122         ! end   counter modification
c....
czr   if (pcomp(flag,adaf)) ke0=ke
czr   if (pcomp(flag,adaf)) kk0=kk
c.... loop over elements (generate compatibility)
c      write(*,'(A)')' - generate compatibility'
      i=0
131   continue              ! not possible with do (counter modified)
      i=i+1
czru      do 131 i=1,ke
czr      if(naiter.eq.1.and.i.lt.nadap)then    !nadap (n1 feap 2. Zeile)
czr         iek0(nen1,i)=abs(iek0(nen1,i))*(-1.d0)
czr         do 117 il =ke0+(i-1)*3+1,ke0+(i-1)*3+3
czr117      iek0(nen1,il)=abs(iek0(nen1,il))*(-1.d0)
czr      endif
c... check if element is marked as 'neighboring' element
        if( iek0(nen1,i).gt.0) goto 130
        k=0
        if(naiter.eq.1) goto 40
c....
        do j=1,4
          if(ianp0(abs(iek0(j,i))).lt.0)k=k+1
        end do
40      continue
c.... no marked node / k-element
        if(k.eq.0.and.iael0(i).ge.0.and.naiter.gt.1)then
czr       write(2,*) 'kgenes ',i
          call kgenes (1,i,x0,x,iek0,ikz0,ike0,   mikno,nen1,ndm,
     1         ianp,ianp0,iael,iael0,naiter,erro,    erron0)
cz      else if( k.eq.1 ) then
        else if( k.eq.1 .or. iael0(i).lt.0 ) then
          if(iael0(i).lt.0) then
c.... one marked node / y-element
czr         write(2,*) 'kymods',i
            nst=0
            call kymods
     1         (nst,i,x0,x,iek0,ikz0,ike0   ,mikno,nen1,ndm,
     2           ianp,ianp0,iael,iael0,naiter,
     3           erro, kf)
          else if( iael0(i).ge.0 ) then
c.... one marked node / k-element
czr         write(2,*) 'ygenes',i
            call ygenes (0,i,x0,x,iek0,ikz0,ike0     ,mikno,nen1,ndm,
     1         ianp,ianp0,iael,iael0,naiter,erro)
          end if
c.... two marked nodes / no action required
        elseif(k.eq.2)then
czr       write(2,*) 'ogenes',i
          call ogenes (1,i,iek0,nen1)
        end if
130     continue
czru131  continue                   ! enddo ...   third loop over elements
         if(i.lt.ke) goto 131       ! end   ...   third loop over elements
c....
czr      if(fsteu.eq.1.and.kf.lt.ct2.and.naiter.ne.1)then
czr        kf=kf+1
czr        write(*,*) 'adaf kf=',kf,'fmax=',fmax
czr       if(pcomp(flag,adaf)) ke0=ke
czr       if(pcomp(flag,adaf)) kk0=kk
czr        goto 100
czr      endif
czrc.... write new FEAP input file (- feap,coor,elem)
czr      call wrinput(kk,ke,x,x0,iek0,ndm,nen1  ,ndf)
czrc.... write restart
czr      fint = fsav
czr      call rmstrt(fint,ndm,ndf,nen1,2,1,
czr     1        kk,ke,naiter,iek,iek0,iael,iael0,ianp,ianp0,
czr     2        ct)
c...
      return
3001  format('**WARNING** Compute nodal stresses before remeshing')
3002  format(I2,1x,'iteration')
3003  format(2x, a15 ,f5.1,' %')
      end
c...............................................................end.SR.admess
c------------------------------------------------------------------------ende
