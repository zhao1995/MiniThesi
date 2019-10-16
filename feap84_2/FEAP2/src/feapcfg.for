      subroutine feapcfg
c-----------------------------------------------------------------------
c
c....   Purpose: configuration for FEAP under windows:
c                read/set data from feapwin.ini
c
c      Input:
c
c      Output:
c
c     WW BS UKA  7/99                                                  |
c     WW IBS KIT 01/14  add search once for Adobe Reader
c-----------------------------------------------------------------------
      USE feapprog
      USE pathn
      USE plodfb
      USE vdata
      USE wincfg
      logical lfil,pcomp
      character*229, allocatable, dimension(:) :: data_ini
      integer*4 iocfg,nperf,ic,im
      real*8    td(10)
      character*40 title
      character yyy*229
      character*229 feapini,ppath
      character*6 wd(17),cc
      character*11 fini
      common /winperf/nperf
c
      data wd/'manual','wininp','wingra','winhlp','editor','psview',
     1        'perfor','feapsa','feapin','feaped','fnege ','fproce',
     2        'intnet','tplopt','fadobe','fgmesh','fcylt'/
      data list/17/
c
      iocfg = 32
c
c.... default-werte
c...1   [manual] manual path
        yyy      = ' '
        file(1)  = 'c:\feap\man\mesh\'
        file(2)  = 'c:\feap\man\macro\'
        file(3)  = 'c:\feap\man\plot\'
        file(4)  = 'c:\feap\man\elmt\'
        file(5)  = 'c:\programme\feap\add\feap.hlp'
c
c...2   [wininp] input window
        title = ' '
        title1 ='F E A P    Dialog   Window'
        fpx1 = 0.60d0
        fpy1 = 0.00d0
        fwx1 = 0.40d0
        fwy1 = 0.95d0

c
c...3   [wingra] graphic window
        title2 ='F E A P    Graphic   Window'
        fpx2 = 0.0d0
        fpy2 = 0.0d0
        fwx2 = 0.595d0
        fwy2 = 0.72d0

c
c...4   [winhlp] user manual window
        title3 ='F E A P    Help   Window'
        fpx3   = 0.0
        fpy3   = 0.725
        fwx3   = 0.60
        fwy3   = 0.23

c
c...5   [editor] user defined editor
         editor ='c:\program files (x86)\UltraEdit\uedit32.exe'

c
c...6   [psview] postscript viewer
         psview ='c:\program files (x86)\gs\gsview\gsview32.EXE'

c
c...7   [perfor] no.eq. to plot performance bar
        nperf = 2500

c
c...8   [feapsa] feap SALFORD
        feapsal ='c:\program files (x86)\feap\feap.exe'

c
c...9   [feapin] feap INTEL
        feapint ='c:\program files (x86)\feap\feap.exe'

c
c..10   [feaped] feap_ed
        feaped ='c:\program files (x86)\feap\add\feap_ed.exe'

c
c..11   [fnege ] netzgenerierer
        fnege  ='c:\program files (x86)\feap\add\negef.exe'

c
c..12   [fproce] procedure editor
        fproce ='c:\program files (x86)\feap\add\f_proced.exe'

c
c..13   [intnet] internet explorer
        fintnet ='c:\program files (x86)\Internet Explorer\iexplore.exe'

c
c..14   [tplopt] no.of points for tplo
        nplo  = 1000

c
c..15   [adobe] acrobat reader
        fadobe =
     +  '"c:\program files (x86)\Adobe\Reader 11.0\Reader\AcroRd32.exe"'
        ifad=1

c
c..16   [fgmesh] netzgenerierer
        fgmesh ='c:\program files (x86)\feap\gmesh.exe'

c..17   [fcylt ] ylt-netzgenerierer
        fcylt ='start c:\program files (x86)\feap\cylt.exe'

c...  end default definitions
c
c
c.... read values from ini-file
c.... get path to prog dir (PCFG)
      ppath = ' '
      call getparap(ppath)
      np = ipos(ppath,229)
c.... get path to ini file (FCFG)
      feapini = ' '
      call getparaf(fpath,fini)
      nf = ipos(fpath,229)
      feapini(1:nf)=fpath(1:nf)
      feapini(nf+1:nf+12)=fini

c.... look to see if ini.file exist
      inquire(file=feapini,exist=lfil)
      if(.not.lfil) then
        write(*,*) 'FEAPwin.ini: Default values are used.'
        write(*,*) 'FEAPwin.ini at ',fpath,feapini
        return
      end if
      if(lfil) then
c....   read values from file
        call opencf(iocfg,feapini)
c
        ic=0 ! count number of input lines
103     read(iocfg,1001) cc
        ic=ic+1
        if(pcomp(cc,'end',3)) then ! end of file reached
          close(iocfg)
c....     modify feapwin.ini if fadobe=' '
          if(ifad.eq.0) then
            ic=ic-1
            allocate(data_ini(ic))
            data_ini=' '
c....       read feapwin.ini
            call opencf(iocfg,feapini)
            do i = 1,ic
              read(iocfg,1000) yyy
              data_ini(i)=yyy
              if(pcomp(yyy,'fadobe',6)) im=i+1
            end do
c....       set path to Adobe Reader
            data_ini(im)=fadobe
            close(iocfg)
c....       write feapwin.ini
            call opencf(iocfg,feapini)
            rewind iocfg
            do i = 1,ic
              write(iocfg,1000) data_ini(i)
            end do
            close(iocfg)
            deallocate(data_ini)
          end if
          return
        end if
c
c....   read data if in list
        do i = 1,list
          if(pcomp(cc,wd(i),6)) go to 104
        end do
        go to 103
104     continue
        ic=ic+1
c
c----------------------------------------------------------------------
c              m  w  w  w  e  p  p  f  f  f  f  f  i  t  f  f  f |
c              a  i  i  i  d  s  e  e  e  e  n  p  n  p  a  g  c |
c              n  n  n  n  i  v  r  a  a  a  e  r  t  l  d  m  y |
c              u  i  g  h  t  i  f  p  p  p  g  o  n  o  o  e  l |
c              a  n  r  l  o  e  o  s  i  e  e  c  e  p  b  s  t |
c              l  p  a  p  r  w  r  a  n  d     e  t  t  e  h    |
      go to ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17) i
c----------------------------------------------------------------------
c
c....   [manual] manual path
1       read(iocfg,1000) yyy
        if(yyy.ne.' ') call cfg_path(file(5),yyy,fpath,nf,ppath,np)
c       call mess_win(file(5),-3)
        goto 103
c
c....   [wininp] input window
2       read(iocfg,1000) title
        if(title.ne.' ') title1 = title
        call dinput2(td,4)
        ic=ic+1
        fsum = dot(td,td,4)
        if(fsum.eq.0.d0) goto 103
        fpx1 = td(1)
        fpy1 = td(2)
        fwx1 = td(3)
        fwy1 = td(4)
c
        goto 103
c
c....   [wingra] graphic window
3       read(iocfg,1000) title
        if(title.ne.' ') title2 = title
        call dinput2(td,4)
        ic=ic+1
        fsum = dot(td,td,4)
        if(fsum.eq.0.d0) goto 103
        fpx2 = td(1)
        fpy2 = td(2)
        fwx2 = td(3)
        fwy2 = td(4)
c
        ic=ic+1
        goto 103
c
c....   [winhlp] user manual window
4       read(iocfg,1000) title
        if(title.ne.' ') title3 = title
        call dinput2(td,4)
        ic=ic+1
        fsum = dot(td,td,4)
        if(fsum.eq.0.d0) goto 103
        fpx3 = td(1)
        fpy3 = td(2)
        fwx3 = td(3)
        fwy3 = td(4)
c
        goto 103
c
c....   [editor] user defined editor
5       read(iocfg,1000)  yyy
        if(yyy.ne.' ') call cfg_path(editor,yyy,fpath,nf,ppath,np)
        goto 103
c
c....   [psview] postscript viewer
6       read(iocfg,1000)  yyy
        if(yyy.ne.' ') call cfg_path(psview,yyy,fpath,nf,ppath,np)
        goto 103
c
c....   [perfor] no.eq. to plot performance bar
7       call dinput2(td,1)
        if(td(1).ne.0.d0) nperf = td(1)
c
        goto 103
c
c....   [feapsa] feap-programm SALFORD
8       read(iocfg,1000)  yyy
        if(yyy.ne.' ') call cfg_path(feapsal,yyy,fpath,nf,ppath,np)
        goto 103
c
c....   [feapin] feap-programm INTEL
9       read(iocfg,1000)  yyy
        if(yyy.ne.' ') call cfg_path(feapint,yyy,fpath,nf,ppath,np)
        goto 103
c
c....   [feaped] feap editor
10      read(iocfg,1000)  yyy
        if(yyy.ne.' ') call cfg_path(feaped,yyy,fpath,nf,ppath,np)
        goto 103
c
c....   [fnege ] netzgenerierer
11      read(iocfg,1000)  yyy
        if(yyy.ne.' ') call cfg_path(fnege,yyy,fpath,nf,ppath,np)
        goto 103
c
c....   [fproce] procedure editor
12      read(iocfg,1000)  yyy
        if(yyy.ne.' ') call cfg_path(fproce,yyy,fpath,nf,ppath,np)
        goto 103
c
c....   [intnet] internet explorer
13      read(iocfg,1000)  yyy
        if(yyy(229:229).ne.' ') goto 900
        if(yyy.ne.' ') fintnet = yyy
c
        goto 103
c
c....   [tplopt] no.of points for tplo
14      call dinput2(td,1)
        if(td(1).ne.0.d0) nplo  = td(1)
c
        goto 103
c
c....   [adobe] Acrobat Reader
15      read(iocfg,1000)  yyy
        if(yyy.eq.' ') then  ! search for position
          call regcheck(fadobe,nf)
          ifad=0
        else                 ! read from ini-file
          call cfg_path(fadobe,yyy,fpath,nf,ppath,np)
          ifad=1
        end if
        goto 103
c
c....   [fgmesh] netzgenerierer
16      read(iocfg,1000)  yyy
        if(yyy.ne.' ') call cfg_path(fgmesh,yyy,fpath,nf,ppath,np)
        goto 103

c....   [fcylt] ylt-netzgenerierer
17      read(iocfg,1000)  yyy
        if(yyy.ne.' ') call cfg_path(fcylt,yyy,fpath,nf,ppath,np)
        goto 103

c
      end if ! from lfil
      write(*,*) 'FEAPwin.ini: Default values are used.'
      write(*,*) 'FEAPwin.ini at ',fpath,feapini

      return

900   write(*,*) 'Path or Filename in FEAPwin.ini is to long(<229!)'
      write(*,*) yyy
      stop
c
1000  format(a)
1001  format(a6)
      end
c
      subroutine cfg_path(yout,y,fpath,nf,ppath,np)
c-----------------------------------------------------------------------
c
c....   Purpose: set path in feapwin.ini variables for FCFG/PCFG
c
c      Input:   y     - original String with FCFG/PCFG
c               fpath - value of FCFG
c               nf    - length of fpath
c               ppath - value of PCFG
c               np    - length of ppath
c
c      Output:  yout  - final    string with PATH from FCFG/PCFG
c               ny    - length of ny
c
c
c     WW IBS KIT  11/12
c     WW IBS KIT  01/14 ppath added for prog directory                                  |
c-----------------------------------------------------------------------
      integer*4 nf,np,ny
      character y*229,fpath*229,ppath*229,yout*229

      yout=' '
      if (y(1:4).eq.'FCFG') then
         ny = ipos(y,229)
         yout(1:nf)=fpath(1:nf)
         yout(nf+1:nf+ny-5)=y(6:ny)  ! fpath with \
      else if (y(1:4).eq.'PCFG') then
         ny = ipos(y,229)
         yout(1:np)=ppath(1:np)
         yout(np+1:np+ny-5)=y(6:ny)  ! ppath with \
      else
         yout=y
      end if
      if(yout(229:229).ne.' ') then
         write(*,*) 'Path/Filename in FEAPwin.ini is to long(<229!)'
         write(*,*)  yout
      end if

c     write(*,*)  yout

      return
      end
c
      subroutine regcheck(path,np)
c-----------------------------------------------------------------------
c
c....   Purpose: find Acrobat reader or similar in registry
c
c      Input:
c
c      Output:  path  - position of Acrobat reader in registry
c               np    - length of path
c
c      Examples:
c        c:\Program Files (x86)\Adobe\Acrobat 10.0\Acrobat\Acrobat.exe
c        c:\Program Files (x86)\Adobe\Reader 10.0\Reader\AcroRd32.exe
c
c     DH TUD  12/13
c-----------------------------------------------------------------------

c....... Import WinAPI routines

#ifdef __INTEL_
      use advapi32
#endif
      implicit none
#ifdef __SALFORD_
      stdcall RegOpenKeyEx 'RegOpenKeyExA' (VAL,STRING,VAL,VAL,REF):
     +    integer
      stdcall RegEnumKeyEx 'RegEnumKeyExA' (VAL,VAL,REF,REF,REF,
     +    REF,REF,REF):integer
      stdcall RegQueryInfoKey 'RegQueryInfoKeyA'(VAL,REF,REF,REF,REF,
     +    REF,REF,REF,REF,REF,REF,REF):integer
      stdcall RegQueryValueEx 'RegQueryValueExA'(VAL,STRING,REF,REF,
     +    REF,REF):integer
      stdcall RegCloseKey 'RegCloseKey' (VAL):integer
#endif
      character(len=80):: subkey, subkey2, enum_subkey
      character(len=229):: cur_path
      character(len=229), intent(out) :: path

      integer :: i,j, len_subkey, reader_found, look_for_pro
      integer, intent(out) :: np

#ifdef __INTEL_
      integer(DWORD) :: buf_size, num_subkeys
      integer(LONG) :: ret
      integer(LPHANDLE) :: hkcu_handle, key_handle, subkey_handle
#endif
c........ FTN95 uses differently sized handles and doesn't know constants
#ifdef __SALFORD_
      integer :: buf_size, num_subkeys
      integer*4 :: ret
      integer   :: hkcu_handle, key_handle, subkey_handle
      integer, parameter :: ERROR_SUCCESS = 0, KEY_READ = Z'20019',
     +    HKEY_CURRENT_USER = Z'80000001'
#endif

      real :: best_version, cur_version

      buf_size = 80
      best_version = 0.e0

c......... clear path first
      do i = 1,len(cur_path)
          cur_path(i:i) = ' '
      end do

c......... try to find Acrobat Reader first
c......... if unsuccessful, find Acrobat Pro

#ifdef __INTEL_
      subkey = "SOFTWARE\\Adobe\\Acrobat Reader"//CHAR(0)
#endif
#ifdef __SALFORD_
      subkey = "SOFTWARE\Adobe\Acrobat Reader"
#endif
      reader_found = 0
      look_for_pro = 0

      write (*,*) "FEAP CONFIGURATION with feapwin.ini"
      write (*,*) "Fetching registry key for Acrobat Reader path..."

#ifdef __INTEL_
c......... open HKCU
      ret = RegOpenCurrentUser(KEY_READ, hkcu_handle)

      if(ret .ne. ERROR_SUCCESS) then
          write (*,'(A,1X,I4)')
     +    "Regcheck error at RegOpenCurrentUser | "//
     +    "Couldn't open HKEY_CURRENT_USER. Windows error code:", ret
          stop
      end if
#endif
#ifdef __SALFORD_
      hkcu_handle = HKEY_CURRENT_USER
#endif

c......... open Software\Adobe\Acrobat Reader
100   ret = RegOpenKeyEx(hkcu_handle, subkey, 0, KEY_READ,
     +            key_handle)
      if(ret .eq. 2) then
          write (*,*)
     +    "Regcheck error at RegOpenKeyEx | "//
     +    "Key = "//subkey(1:len_trim(subkey)-1)
          write(*,*)
     +    "Couldn't open key, error 2: File not found."
          write (*,*) "Program not installed?"
      else if(ret .ne. ERROR_SUCCESS) then
          write (*,'(A,1X,I4)')
     +    "Regcheck error at RegOpenKeyEx | "//
     +    "Key = "//subkey(1:len_trim(subkey)-1)
          write(*,*)
     +    "Couldn't open key. Windows error code:", ret
      else !if ret.eq. ERROR_SUCCESS

c......... number of subkeys?
#ifdef __INTEL_
          ret = RegQueryInfoKey(key_handle, 0, 0, 0, LOC(num_subkeys),
     +            0, 0, 0, 0, 0, 0, 0)
#endif
#ifdef __SALFORD_
          ret = RegQueryInfoKey(key_handle, CORE4(0), CORE4(0),
     +            CORE4(0), num_subkeys, CORE4(0), CORE4(0),
     +            CORE4(0), CORE4(0), CORE4(0), CORE4(0), CORE4(0))
#endif
          if(ret .ne. ERROR_SUCCESS) then
              write (*,'(A,1X,I4)')
     +        "Regcheck error at RegQueryInfoKey | "//
     +        "Couldn't query info for key. Windows error code:", ret
              stop
          end if

c      write (*,'(A,1X,A,1X,I3,1X,A)')
c     + subkey(1:35), "has", num_subkeys, "subkey(s):"


c......... query install path for each subkey, try to use the last
          do i = 0, num_subkeys-1
c......... identify subkey
#ifdef __INTEL_
              ret = RegEnumKeyEx(key_handle, i, enum_subkey,
     +                LOC(buf_size), 0, 0, 0, 0)
#endif
#ifdef __SALFORD_
              ret = RegEnumKeyEx(key_handle, i, enum_subkey,
     +                buf_size, CORE4(0), CORE4(0), CORE4(0), CORE4(0))
#endif
              ! reset buffer size because WinApi overrides it
              buf_size = 80
              if(ret .ne. ERROR_SUCCESS) then
                  write (*,'(A,1X,I4)')
     +            "Regcheck error at RegEnumKeyEx | "//
     +            "Couldn't enumerate subkey. Windows error code:", ret
              else
                  len_subkey = 0
                  j = 1
                  do while (enum_subkey(j:j).ne.CHAR(0))
                      len_subkey = len_subkey+1
                      j = j+1
                  end do
#ifdef __INTEL_
                  subkey2 = subkey(1:len_trim(subkey)-1) // "\\" //
     +            enum_subkey(1:len_subkey) // "\\InstallPath"//CHAR(0)
#endif
#ifdef __SALFORD_
                  subkey2 = subkey(1:len_trim(subkey))//"\"//
     +            enum_subkey(1:len_subkey)//"\InstallPath"
#endif
c              write (*,*) "  * "//subkey2(1:len_trim(subkey2))
              end if

c......... open subkey
              ret = RegOpenKeyEx(hkcu_handle, subkey2, 0,
     +                KEY_READ, subkey_handle)
              if(ret .ne. ERROR_SUCCESS) then
                  write (*,'(A,1X,I4)')
     +        " Regcheck error at RegOpenKeyEx (subkey) | "//
     +        "Subkey = "//subkey2(1:len_trim(subkey2))
                  write (*,*)
     +        "Couldn't open subkey. Windows error code:", ret
c          else
c              write (*,'(A,1X,A,1X,A)') "Subkey",subkey2,"opened."
              else
c......... query install path for subkey
#ifdef __INTEL_
                  ret = RegQueryValueEx(subkey_handle, "", 0, 0,
     +                    LOC(cur_path),  LOC(buf_size))
#endif
#ifdef __SALFORD_
                  ret = RegQueryValueEx(subkey_handle, "", CORE4(0),
     +                     CORE4(0), cur_path,  buf_size)
#endif
                  ! reset buffer size because WinApi overrides it
                  buf_size = 80
                  if(ret .eq. 234) then
                   write (*,*) "RegQueryValueEx error: Data too long!"//
     +             "Provide larger buffer!"
                  else if(ret .ne. ERROR_SUCCESS) then
                      write (*,'(A,1X,I4)')
     +                "Regcheck error at RegQueryValueEx | "//
     +                "Couldn't query value. Windows error code:", ret
                  else
                      reader_found = 1

                      read(enum_subkey(1:len_subkey),*) cur_version
                      write (*,'(1X,A,1X,F4.1)')
     +                "Found version", cur_version

                      if(cur_version.gt.best_version) then
                          write (*,*) "It's the new best version."
                          best_version = cur_version
                          path = cur_path
                      end if
                  end if
              end if

          end do

      end if

      if (reader_found.eq.0.and.look_for_pro.eq.0) then
          write (*,*) " "
          write (*,*) "Fetching for Acrobat Pro instead..."
#ifdef __INTEL_
          subkey = "SOFTWARE\\Adobe\\Adobe Acrobat"//CHAR(0)
#endif
#ifdef __SALFORD_
          subkey = "SOFTWARE\Adobe\Adobe Acrobat"
#endif
          look_for_pro = 1
          goto 100
      else if (reader_found.eq.0.and.look_for_pro.eq.1)then
          write (*,*) "No Acrobat installation found!"
          path = ""
          np = 0
          return
      else if (reader_found.eq.1.and.look_for_pro.eq.0)then
c......... Acrobat Reader found!
          np = len_trim(path)-1  ! trailing '\0' still in string
          path = path(1:np)//"\AcroRd32.exe"
          np = len_trim(path)
      else if (reader_found.eq.1.and.look_for_pro.eq.1)then
c......... Acrobat Pro found!
          np = len_trim(path)-1  ! trailing '\0' still in string
          path = path(1:np)//"\Acrobat.exe"
          np = len_trim(path)
      end if

      write (*,*) "ADOBE Path found: ", path(1:np)
c      write (*,'(1X,A,I3)') "Length of path: ", np


      write (*,*) " "
      write (*,*) "Feapwin.ini will be modified with this path"
      write (*,*) "Continuation of Program with input <CR>!"

      read (*,*)

#ifdef __SALFORD_
      call close_output_window
#endif

      return
      end