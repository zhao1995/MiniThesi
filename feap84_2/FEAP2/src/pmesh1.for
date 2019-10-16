      subroutine pmesh1(numnp,numel,numat,ndm,ndf,nen)
      USE pathn
      USE pdata2
      implicit double precision (a-h,o-z)
      character hfile*60,hfile1*229,hfile2*60
      logical lexst1,lexst2,lexst3,lexst4
#ifndef _NOLIC_
      if(idev.le.2) then
cww     hfile=' '
cww     call getenv('HOME',hfile)   !  getenv is unix command
cww     do  i=60,1,-1
cww       if (hfile(i:i).ne.' ') goto 100
cww     end do
cww100  continue
cww     hfile(i+1:i+9)='/.cal.dll'
        hfile = 'c:\windows\system\cal.dll'
      end if

      if(idev.eq.3.or.idev.eq.4) hfile  = 'c:\windows\system\cal.dll'
      if(idev.eq.3.or.idev.eq.4) hfile2 = 'c:\windows\system32\cal.dll'
      inquire(file=hfile,exist=lexst1)
      inquire(file=hfile2,exist=lexst2)
      lexst3=.true.
c     hfile1 = file(5)
c     inquire(file=hfile1,exist=lexst3)
      lexst4=.true.
      call pmesh3(itd,itm,ity)

      call tend(lexst4,itd,itm,ity)
      if(.not.lexst4) then
        il=license()
        stop
      end if
c
      write(16,*) lexst1,lexst2,lexst3,lexst4
      if((lexst1 .or. lexst2) .and. lexst3) then
        return
      else
        numnp=numel+1
        ndm=ndf+2
        nen=ndf-1
        ndf=numel
      end if
#else

#endif
      return
      end
c
      subroutine pmesh2(n1,n2,n3,n4,n5,n6,n7,n8,n9,
     +                  ie,numnp,numel,nummat,nie)
      implicit double precision (a-h,o-z)
      dimension ie(nie,*)
      return
      end
c
      subroutine pmesh3(itd,itm,ity)
      USE feapprog
      USE iodata
      USE iofile
      implicit double precision (a-h,o-z)
      logical lreg
      character key*72, citd*1,citm*1,
     +          cd4*4,cd5*5,cd7*7,cd12*12,cd13*13,cd21*21
      character*4 city,city1*1,city2*1,city3*1,city4*1
      character*225 feapreg
      character regfile*25
      regfile = 'feap_registration_key.txt'
      np1 = ipos(fpath,229)
      feapreg = ' '
      feapreg(1:np1)=fpath(1:np1)
      feapreg(np1+1:np1+26)=regfile
      inquire(file=feapreg,exist=lreg)
      if(lreg) then
        open(ios,file=feapreg,status='unknown')
        read(ios,1000) key
        close(ios)
        read(key,'(a5,a1,a7,a1,a4,a1,a12,a1,a13,a1,a4,a1,a21)')
     +  cd5,city4,cd7,city1,cd4,citm,cd12,city3,cd13,citd,cd4,city2,cd21
        itd = ichar(citd)-64
        itm = ichar(citm)-64
        ity1= (ichar(city1)-48)*1000
        ity2= (ichar(city2)-48)*100
        ity3= (ichar(city3)-48)*10
        ity4= (ichar(city4)-64)
        if(ity4.eq.10) ity4=0
        ity = ity1+ity2+ity3+ity4
      else
        itd=01
        itm=01
        ity=2000
      end if
      if(ior.lt.0) write(  *,1001) itd,itm,ity
                   write(iow,1001) itd,itm,ity
      return
1000  format(a72)
1001  format(/,'  Expiration date: ',i2,'.',i2,'.',i4,/)
      end