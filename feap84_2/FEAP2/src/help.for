c-----------------------------------------------------------------------
c
      subroutine pman(isw,name1)
c-----------------------------------------------------------------------
c
c.... Purpose: open pdf-manual for chosen macro command 
c
c     Input:
c       isw  - 1=mesh,2=macro,3=plot,4=titl 
c       name - name of macro    
c
c     Output:
c       open pdf-manual
c
c     Comments:  (start) program.exe "file docu"
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------

      implicit integer*4 (i-n)
      character name1*4,name*4,ext*3,extname*20
 
      name(1:4)=name1(1:4) 

c.... set extension
      if(isw.eq.1) ext='me_'        
      if(isw.eq.2) ext='ma_'        
      if(isw.eq.3) ext='pl_'        
      if(isw.eq.4) ext='   '        

      extname=' '
      
      if(ext.eq.'   ') then 
        extname(1:4)=name(1:4)
      else
        extname(1:3)= ext(1:3)
        extname(4:7)=name(1:4)
      end if

      call start_manual(extname)

      return
      end
      
c-----------------------------------------------------------------------
      subroutine start_manual(name)
c-----------------------------------------------------------------------
c
c.... Purpose: start manual with special entry 
c
c     Inputs:
c       fadobe - acrobat reader
c       fpath  - path to pdf manual
c       name   - position to go    
c       fcis1 - not used
c
c     W. Wagner BS KIT 11/12
c-----------------------------------------------------------------------
      USE feapprog
      character fcis*229,fcis1*229,psfile*20,name*20
      psfile='feap-man.pdf'

c.... Acrobat _/A "nameddest=ext_name"  "path\feap-man.pdf"     

c.... set Program
      fcis = ' '
      np1=0
        
      np2 = ipos(fadobe,229)
      fcis(np1+1:np1+np2) =  fadobe(1:np2)
      np1=np1+np2

c     ityp=1 acrobat.exe ityp=2 gotoacro.exe
      ityp=1

      fcis1 = ' '
      np1=0

      if(ityp.eq.1) then
c....   [ /A "nameddest=]
        fcis1(np1+1:np1+15) = ' /A nameddest=' 
        np1=np1+14 

c....   name
        np2 = ipos(name,20)
        fcis1(np1+1:np1+np2) = name(1:np2)
        np1 = np1+np2      
      
c....   [_"path\feap-man.pdf"]
        fcis1(np1+1:np1+3) = ' "' 
        np1 = np1+2      
      
        np2 = ipos(fpath,229)
        fcis1(np1+1:np1+np2) = fpath(1:np2)
        np1 = np1+np2      
           
        np2 = ipos(psfile,20)
        fcis1(np1+1:np1+np2) = psfile(1:np2)
        np1 = np1+np2      

        fcis1(np1+1:np1+1) = '"' 

      else if(ityp.eq.2) then
c....   ["path\feap-man.pdf"]
        fcis1(np1+1:np1+1) = '"' 
        np1 = np1+1      
      
        np2 = ipos(fpath,229)
        fcis1(np1+1:np1+np2) = fpath(1:np2)
        np1 = np1+np2      
           
        np2 = ipos(psfile,20)
        fcis1(np1+1:np1+np2) = psfile(1:np2)
        np1 = np1+np2      

        fcis1(np1+1:np1+2) = '" ' 
        np1 = np1+2      

c....   name
        np2 = ipos(name,20)
        fcis1(np1+1:np1+np2) = name(1:np2)
      
      end if

c      call mess_win(fcis,-3)
c      call mess_win(fcis1,-3)



c.... start process 
      call start_pprocess(fcis,fcis1)

      return 
      end      
c
c-----------------------------------------------------------------------
c
            subroutine phelpm(wd,nwd,clab2)
c-----------------------------------------------------------------------
c
c.... Purpose: open list of available macro commands 
c
c     Input:
c       wd     - list of commands 
c       nwd    - length of list    
c       clab2 
c
c     Output:
c        
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata2
      character*4 wd(nwd),wdd(150),yyy*80,name*4,clab2
      integer ii
        
      nwdd = nwd + 6
      do  i = 1,150
          wdd(i) = '.'
      end do

c.... copy
      do  i = 1,nwd
          wdd(i) = wd(i)
      end do

c.... add some macros
      wdd(nwd+1) = 'end '
      wdd(nwd+2) = 'exit'
      wdd(nwd+3) = 'quit'
      wdd(nwd+4) = 'hist'
      wdd(nwd+5) = 'help'
      wdd(nwd+6) = 'proc'

c.... sort list of macros  
      call psortm(wdd,nwdd,4)

      if(idev.eq.1.or.idev.eq.2.or.idev.eq.3) then
c....   print in console window

10      write(*,2000)
        ii = int((nwdd)/9)+1
        do 11 i = 1,ii
          iwd = (i-1)*9
          ir  = nwdd - iwd
          ia  = iwd+1
          if(ir.gt.9) ie  = ia + 8
          if(ir.le.9) ie  = nwdd
          write(yyy,'(9(3x,a4))') (wdd(k),k=ia,ie)
          write(*,1000) yyy
11      continue
        write(*,2002)  
        yyy=' '
        read(*,1000) yyy
        name = yyy(1:4)
        if(name.eq.' ') return

c....   show macro
        call pman(2,name)

        goto 10

      else if(idev.eq.4) then
c....   show windows help menue
40      name=' '
c....   find macro
        call pwinhelp(name,wdd,nwdd,'Macro')
        if(name.eq.' ') return

c....   show macro
        call pman(2,name)
        goto 40 

      endif
      return
1000  format(a80)
2000  format(/' The following MACRO commands are available:',/,
     +        ' Overview on MACRO commands:   Index',/,
     +        ' Overview on possible actions: Action')
2002  format(' EXIT with <CR> or  enter <name> for more help:',$)
      end
c

c-----------------------------------------------------------------------
c
      subroutine phelpmp(wd,nwd,wrd,c2)
c-----------------------------------------------------------------------
c
c.... Purpose: open list of available mesh/plot commands 
c
c     Input:
c       wd     - list of commands 
c       nwd    - length of list    
c       wrd    -
c       c2 
c
c     Output:
c        
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata2
      integer  ii
      character*4 wd(nwd),wdd(150),wrd,yyy*80,name*4,c2,wrd5*5
      character*35 title

      write(title,2005) wrd
      do  i = 1,150
        wdd(i) = '.'
      end do

c.... copy
      do  i = 1,nwd
        wdd(i) = wd(i)
      end do

c.... sort
      call psortm(wdd,nwd,4)

      if(idev.eq.1.or.idev.eq.2.or.idev.eq.3) then
c....   print in console window

10      write(*,2000) wrd,wrd
        ii = int((nwd)/9)+1
        do 11 i = 1,ii
          iwd = (i-1)*9
          ir  = nwd - iwd
          ia  = iwd+1
          if(ir.gt.9) ie  = ia + 8
          if(ir.le.9) ie  = nwd
          write(yyy,'(9(3x,a4))') (wdd(k),k=ia,ie)
          write(*,1000) yyy
11      continue
        write(*,2002)  
        yyy=' '
        read(*,1000) yyy
        name = yyy(1:4)
        if(name.eq.' ') return
    
        if(wrd.eq.'MESH') iwrd=1
        if(wrd.eq.'PLOT') iwrd=3

c....   show macro
        call pman(iwrd,name)

        goto 10

      elseif(idev.eq.4) then
c....   show windows help menue

        if(wrd.eq.'MESH') then
        iwrd=1
        wrd5='MESH '

        else if(wrd.eq.'PLOT') then
        iwrd=3
        wrd5='PLOT '
        end if

40      name=' '
c....   find macro
        call pwinhelp(name,wdd,nwd,wrd5)
        if(name.eq.' ') return

c....   show macro
        call pman(iwrd,name)

        goto 40 
      end if

      return
1000  format(a80)
2000  format(/' The following ',a4,' commands are available:',/,
     +        ' Overview on   ',a4,' commands: Index',/,
     +        ' Overview on possible actions:  Action')
2002  format(' EXIT with <CR> or  enter <name> for more help:',$)
2005  format('FEAP ',a4,' Commands')
      end
c

c-----------------------------------------------------------------------
c
      subroutine psortm(wdd,nwdd,nc)
c-----------------------------------------------------------------------
c
c.... Purpose: open pdf-manual for chosen macro command 
c
c     Input:
c       wdd  - list of macros unsorted
c       nwdd - number of inputs in list    
c       nc   -
c
c     Output:
c       wdd  - list of macros sorted
c
c     W. Wagner BS UKA 04/09
c--------------------------------------------------------------------
      implicit integer(i-n)
      character*1 wdd(nc,nwdd),a(80)
      character*1 a1,a2
      logical lgt
10    m = 0
      do i = 1,nwdd-1
        a1 = wdd(1,i)
        a2 = wdd(1,i+1)
c.....first character
        if(lgt(a1,a2)) then
          do j = 1,nc
            a(j)=wdd(j,i)
            wdd(j,i)=wdd(j,i+1)
            wdd(j,i+1)=a(j)
            m = 1
          enddo
c.....second character
        elseif(a1.eq.a2) then
          a1 = wdd(2,i)
          a2 = wdd(2,i+1)
          if(lgt(a1,a2)) then
            do j = 1,nc
              a(j)=wdd(j,i)
              wdd(j,i)=wdd(j,i+1)
              wdd(j,i+1)=a(j)
              m = 1
            enddo
c.....third character
          elseif(a1.eq.a2) then
            a1 = wdd(3,i)
            a2 = wdd(3,i+1)
            if(lgt(a1,a2)) then
              do j = 1,nc
                a(j)=wdd(j,i)
                wdd(j,i)=wdd(j,i+1)
                wdd(j,i+1)=a(j)
                m = 1
              enddo
c.....fourth character
            elseif(a1.eq.a2) then
              a1 = wdd(4,i)
              a2 = wdd(4,i+1)
              if(lgt(a1,a2)) then
                do j = 1,nc
                  a(j)=wdd(j,i)
                  wdd(j,i)=wdd(j,i+1)
                  wdd(j,i+1)=a(j)
                  m = 1
                enddo
              endif
            endif
          endif
        endif
      enddo
      if(m.eq.1) goto 10
      return
      end
c
