c
c     FEAP_SHMEM
c
c
c     Feap interface to handle FE^2 macro/micro communication
c     via Windows-API shared memory file mappings.
c
c     Feap must be executed with administrator privileges
c     to be able to use these file mappings.
c
c     Only used in Intel version, FTN95 does not support accessing
c     shared memory segments with several process instances of the same
c     executable.
c
c     D.Heller TUD Apr2011-Jun2012
c
c-----------------------------------------------------------------------
c     Compiler Options INTEL
c     Fortran>Optimization>Parallelization: /Qparallel
c     Fortran>Language>Process OpenMP Directives> Disable  in  DEBUG
c     Fortran>Language>Process OpenMP Directives> Disable  in  DEBUG-S
c     Fortran>Libraries>Use Intel Math Kernel Library: /Qmkl:parallel
c     WW KIT 05/11
c----------------------------------------------------------------------

! c wrapper funktionen fuer memory mapped file zeugs
! unter linux: mmap shm_open munmap shm_unlink
! memcpy gibts es unter windows und linux
! ACHTUNG fortran pointer != c pointer
      module memorymap

	! O_CREAT | O_RDWR
	integer, parameter :: ocreate = 66
	! O_RDONLY
	integer, parameter :: oread = 0
	! PROT_READ | PROT_WRITE
	integer, parameter :: preadwrite = 3
	! PROT_READ
	integer, parameter :: pread = 1
	! MAP_SHARED
	integer, parameter :: mshare = 1
	! S_IRUSR | S_IWUSR
	integer, parameter :: fperm = 384

      interface
      subroutine  memcpy(dest, src, n) bind(C,name='memcpy')
      use iso_c_binding
      INTEGER(c_intptr_t), value:: dest
      INTEGER(c_intptr_t), value:: src
      integer(c_size_t), value :: n
      end subroutine memcpy
      end interface

#ifdef linux
	interface
      INTEGER(c_intptr_t) function mmap(addr,len,prot, \
	flags,fildes,off) result(result) bind(c,name='mmap') 
      use iso_c_binding 
      INTEGER(c_intptr_t), value :: addr 
      integer(c_size_t), value :: len 
      integer(c_int), value :: prot 
      integer(c_int), value :: flags 
      integer(c_int), value :: fildes 
      integer(c_size_t), value :: off 
      end function mmap 
      end interface
      
      interface
      integer(c_int) function munmap(addr, len) \
	bind(c,name='munmap')
      use iso_c_binding 
      INTEGER(c_intptr_t), value :: addr 
      integer(c_size_t), value :: len 
      end function munmap
      end interface

      interface
      integer(c_int) function ftruncate(fd, len) \
	bind(c,name='ftruncate')
      use iso_c_binding 
      integer(c_int), value :: fd
      integer(c_size_t), value :: len
      end function ftruncate
      end interface
	  
	interface
      integer(c_int) function close(fd) \
	bind(c,name='close')
      use iso_c_binding 
      integer(c_int), value :: fd
      end function close
      end interface
  
      interface
      integer(c_int) function shm_open(name,oflag,mode) \
	bind(c,name='shm_open')
      use iso_c_binding 
      character(kind=c_char) :: name(*)
      integer(c_int), value :: oflag
      integer(c_int16_t), value :: mode
      end function shm_open
      end interface
  
      interface
      integer(c_int) function shm_unlink(name) \
      bind(c,name='shm_unlink')
      use iso_c_binding
      character(kind=c_char) :: name(*)
      end function shm_unlink
      end interface
#endif
    
      end module memorymap

#ifdef __INTEL_

#define SIZEOF_INT  4
#define SIZEOF_REAL 8
#define SIZEOF_CHAR 1


#define LEN_OF_MAPNAME 25

      SUBROUTINE feap_shmem_ParseRestartMappingName(filename,mapname)
c ----------------------------------------------------------------------
c.... Purpose: Convert string with path and filename of a restart file
c....   (e.g. 'C:\FEAP\exe\RVE\4\rrve4_01.0001') to a mapping name
c....   by splitting the string and only returning the filename without
c....   path. (here: 'rrve4_01.0001')
c ----------------------------------------------------------------------
c      USE KERNEL32
      use memorymap
      use iso_c_binding
      IMPLICIT NONE
      CHARACTER(len=229), INTENT(IN) :: filename
      CHARACTER(len=24), INTENT(OUT) :: mapname
      INTEGER :: ifirstchar, ilastchar, idx

c  find last '\' character and first blank
      ifirstchar = 0
      ilastchar = LEN(filename)

      DO idx = 1, 229
        IF (filename(idx:idx).eq.'/'.or.filename(idx:idx).eq.'\') THEN
            ifirstchar = idx+1
        END IF
        IF (filename(idx:idx).eq.' '.and.ifirstchar.ne.0) THEN
            ilastchar = idx-1
            EXIT
        END IF
      END DO

      mapname = 'Global\\rrve4_XX.XXXX'C

c     this command creates a blank char at index (10+ilastchar-ifirstchar), so...
      mapname(8:9+ilastchar-ifirstchar) = filename(ifirstchar:ilastchar)

c     terminate string
      idx = 0 ! dummy 0 in memory
      CALL memcpy ( LOC(mapname)+9+ilastchar-ifirstchar, LOC(idx),1)

      END SUBROUTINE

      SUBROUTINE feap_shmem_CreateRestartMapping(rmap,name,irecl)
c ----------------------------------------------------------------------
c.... Purpose: Create a file-mapping in memory for storing restart
c....   values in non-linear problems, according to the routine
c....   restrt(.) in feaps5.for
c....  The mapping size is determined by irecl (number of entries),
c....   multiplied with 4 (bytes per entry).
c ----------------------------------------------------------------------
#ifndef _LINUX_
      USE KERNEL32
#else
      use ifport
#endif
      use, intrinsic :: iso_c_binding

      use memorymap
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: irecl
      INTEGER(c_intptr_t), INTENT(OUT) :: rmap
      INTEGER :: len
      INTEGER(c_intptr_t) :: view
      CHARACTER(len=28) :: name
      INTEGER :: i

      len = irecl * 4
#ifdef _LINUX_
      rmap = shm_open(name, ocreate, fperm)
      if (rmap == -1) then
          stop 'shm_open error'
      end if
      i = ftruncate(rmap, len)
#else
      rmap = CreateFileMapping( INVALID_HANDLE_VALUE, 0,
     &                              PAGE_READWRITE, 0,
     &                              len,
     &                              name )
#ifdef __SHMDEBUG_
      WRITE (*,*) "Created Restart Mapping '",name
      WRITE (*,*) "'."
#endif
#endif
      IF (rmap <= 0) THEN
        WRITE (*,*) "Warning: Failed to create restart file mappin
     &g. Please run with administrator privileges. (", name, "  )"
        WRITE (*,*) "   Error code:", GetLastError()
        READ (*,*)
        STOP
      END IF

      END SUBROUTINE

      SUBROUTINE feap_shmem_OpenEMapping(emap,basename)

#ifndef _LINUX_
      USE KERNEL32
#endif
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE
      INTEGER(c_intptr_t), INTENT(OUT) :: emap
      CHARACTER(len=LEN_OF_MAPNAME), INTENT(IN) :: basename
      CHARACTER(len=LEN_OF_MAPNAME) :: name

c...  basename is 'Global\FEAPmNNNTTTEEEEEE'C where NNN will be
c...  the lgp number and TTT is replaced with EMA, CSM or NSS.
c...  EEEEEE = ngp i.e. element number

      CALL memcpy( LOC(name), LOC(basename),
     &  LEN_OF_MAPNAME*SIZEOF_CHAR )
      name(16:18) = 'EMA'
      CALL feap_shmem_OpenMapping(emap,name)

      END SUBROUTINE

      SUBROUTINE feap_shmem_OpenCSMapping(csmap,basename)

#ifndef _LINUX_
      USE KERNEL32
#endif
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE
      INTEGER(c_intptr_t), INTENT(OUT) :: csmap
      CHARACTER(len=LEN_OF_MAPNAME), INTENT(IN) :: basename
      CHARACTER(len=LEN_OF_MAPNAME) :: name

c...  basename is 'Global\FEAPmNNNTTT'C where NNN will be
c...  with lgp number and TTT is replaced with EMA, CSM or NSS.

      CALL memcpy( LOC(name), LOC(basename),
     &      LEN_OF_MAPNAME*SIZEOF_CHAR )
      name(16:18) = 'CSM'
      CALL feap_shmem_OpenMapping(csmap,name)

      END SUBROUTINE

      SUBROUTINE feap_shmem_OpenNSSMapping(nssmap,basename)

#ifndef _LINUX_
      USE KERNEL32
#endif
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE
      INTEGER(c_intptr_t), INTENT(OUT) :: nssmap
      CHARACTER(len=LEN_OF_MAPNAME), INTENT(IN) :: basename
      CHARACTER(len=LEN_OF_MAPNAME) :: name

c...  basename is 'Global\FEAPmNNNTTT'C where NNN will be
c...  with lgp number and TTT is replaced with EMA, CSM or NSS.

      CALL memcpy( LOC(name), LOC(basename),
     &      LEN_OF_MAPNAME*SIZEOF_CHAR )
      name(16:18) = 'NSS'
      CALL feap_shmem_OpenMapping(nssmap,name)

      END SUBROUTINE

      SUBROUTINE feap_shmem_OpenMapping(hndl,name)
c ----------------------------------------------------------------------
c.... Purpose: Open File Mappings only (no views of file) for e and c,s
c....       which have already been created and return their handles.
c ----------------------------------------------------------------------
#ifndef _LINUX_
      USE KERNEL32
#else
      use ifport
#endif
      use, intrinsic :: iso_c_binding
      use memorymap

      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(OUT) :: hndl
      CHARACTER(len=LEN_OF_MAPNAME), INTENT(IN) :: name
#ifdef _LINUX_
      hndl = shm_open(name, oread, fperm)
      IF (hndl == -1) THEN
        WRITE (*,*) "Warning: Could not open mapping '",name,"'!"
        WRITE (*,*) "   Error code:", GetLastError()
        READ (*,*)
        STOP
      END IF
#else
      hndl = OpenFileMapping( FILE_MAP_ALL_ACCESS, FALSE,
     &                              name )
      IF (hndl.eq.0) THEN
        WRITE (*,*) "Warning: Could not open mapping '",name,"'!"
        WRITE (*,*) "   Error code:", GetLastError()
        READ (*,*)
        STOP
      END IF
#endif
      END SUBROUTINE

      SUBROUTINE feap_shmem_CreateMappings(emap,csmap,nssmap,nss,lgp,
     &ngp)
c ----------------------------------------------------------------------
c.... Purpose: Create File Mappings only (no views of file) for e and c,s
c....       and return their handles.
c....    * nss defines the size of arrays E, C, S
c....    * lgp is the current gauss point number and is used to identify
c....      the distinct shared memory segments. Must be 3 digits or shorter!
c....    * ngp is the current element number and will be needed for parallel
c....      element loop calculations as well as using different materials.
c....      Must be 6 digits or shorter! (ie. only 10^6 - 1 elements supported)
c....
c ----------------------------------------------------------------------
#ifndef _LINUX_
      USE KERNEL32
#endif
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER, INTENT(IN):: nss,lgp,ngp
      CHARACTER(len=25) :: name
      CHARACTER(len=3) :: slgp
      CHARACTER(len=6) :: sngp

      INTEGER(c_intptr_t), INTENT(OUT) :: emap, csmap, nssmap
      integer :: i
      
      WRITE(slgp,'(i3.3)') lgp
      WRITE(sngp,'(i6.6)') ngp

#ifdef _LINUX_
      name = 'Global\\FEAPmNNNTTTEEEEEE'C

      name(13:15) = slgp
      name(19:24) = sngp
      name(16:18) = 'EMA'
c      write (*,*) name
      emap = shm_open(name, ocreate, fperm)
      if (emap == -1) then
          stop 'shm_open error'
      end if
      i = ftruncate(emap, (nss+4)*SIZEOF_REAL + 2*SIZEOF_INT)
#ifdef __SHMDEBUG_
      WRITE (*,*) "Created Mapping '",name
      WRITE (*,*) "'."
#endif

      name(16:18) = 'CSM'
c      write (*,*) name
      csmap = shm_open(name, ocreate, fperm)
      if (csmap == -1) then
          stop 'shm_open error'
      end if
      i = ftruncate(csmap, (nss*nss+nss)*SIZEOF_REAL)
#ifdef __SHMDEBUG_
      WRITE (*,*) "Created Mapping '",name
      WRITE (*,*) "'."
#endif

      name(16:18) = 'NSS'
c      write (*,*) name
      nssmap =shm_open(name, ocreate, fperm)
      if (nssmap == -1) then
          stop 'shm_open error'
      end if
      i = ftruncate(nssmap,1*SIZEOF_INT)
#else
c     * NNN is a placeholder for the gauss point number LGP
c     * TTT is a placeholder for the type specification: EMA, CSM or NSS
c       for the E-mapping, CS-mapping and NSS-mapping respectively
c     * EEEEEE is a placeholder for the element number NGP
      name = 'Global\\FEAPmNNNTTTEEEEEE'C

      name(13:15) = slgp
      name(19:24) = sngp
      name(16:18) = 'EMA'
c      write (*,*) name
      emap = CreateFileMapping( INVALID_HANDLE_VALUE, 0,
     &                              PAGE_READWRITE, 0,
     &                              (nss+4)*SIZEOF_REAL + 2*SIZEOF_INT,
     &                              name )
#ifdef __SHMDEBUG_
      WRITE (*,*) "Created Mapping '",name
      WRITE (*,*) "'."
#endif

      name(16:18) = 'CSM'
c      write (*,*) name
      csmap = CreateFileMapping(INVALID_HANDLE_VALUE, 0,
     &                              PAGE_READWRITE, 0,
     &                              (nss*nss+nss)*SIZEOF_REAL,
     &                              name )
#ifdef __SHMDEBUG_
      WRITE (*,*) "Created Mapping '",name
      WRITE (*,*) "'."
#endif

      name(16:18) = 'NSS'
c      write (*,*) name
      nssmap = CreateFileMapping(INVALID_HANDLE_VALUE, 0,
     &                              PAGE_READWRITE, 0,
     &                              1*SIZEOF_INT,
     &                              name )
#ifdef __SHMDEBUG_
      WRITE (*,*) "Created Mapping '",name
      WRITE (*,*) "'."
#endif
#endif
      END SUBROUTINE

      SUBROUTINE feap_shmem_CloseHandle(hndl)
c ----------------------------------------------------------------------
c.... Purpose: Close file mapping handle opened before
c ----------------------------------------------------------------------
#ifndef _LINUX_
      USE KERNEL32
#endif
      use memorymap
      use, intrinsic :: iso_c_binding
      IMPLICIT NONE

      INTEGER(c_intptr_t) :: hndl
      INTEGER :: ret
#ifdef _LINUX_
      if(hndl.gt.0) THEN
          close(hndl)
      end if
#else
      IF (hndl .ne. 0) THEN
        ret = CloseHandle( hndl )
      END IF
#endif
      END SUBROUTINE

      SUBROUTINE feap_shmem_OpenView(mapping,size,view)
c ----------------------------------------------------------------------
c.... Purpose: Open and return a map view to a previously opened mapping.
c....       If size is zero, the whole mapping will be available in the
c....       opened view.
c ----------------------------------------------------------------------
#ifndef _LINUX_
      USE KERNEL32
#endif
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: size
      INTEGER(c_intptr_t), INTENT(IN) :: mapping
      INTEGER(c_intptr_t), INTENT(OUT) :: view
#ifdef _LINUX_
      view = mmap(loc(0), size, preadwrite, mshare, mapping, 0)
#else
      view = MapViewOfFile( mapping, FILE_MAP_ALL_ACCESS, 0, 0, size )
#endif
      END SUBROUTINE

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccc ab hier nur map/unmap/memcpy aufrufe cccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE feap_shmem_writeNSS(nssmap,nss)
c ----------------------------------------------------------------------
c.... Purpose: Write nss into given file mapping
c....       nssmap.
c ----------------------------------------------------------------------
#ifndef _LINUX_
      USE KERNEL32
#endif
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nss

      INTEGER(c_intptr_t), INTENT(IN) :: nssmap
      INTEGER(c_intptr_t) :: view
      INTEGER :: ret
#ifdef _LINUX_
      view = mmap(loc(0), 1*SIZEOF_INT, preadwrite, mshare, nssmap, 0)
      CALL memcpy( view, LOC(nss), 1*SIZEOF_INT )
      ret = munmap(view, 1*SIZEOF_INT )
#else
      view = MapViewOfFile( nssmap, FILE_MAP_ALL_ACCESS, 0, 0,
     &                      1*SIZEOF_INT )

      CALL memcpy( view, LOC(nss), 1*SIZEOF_INT )
      ret = UnmapViewOfFile(nssmap)
#endif

      END SUBROUTINE

      SUBROUTINE feap_shmem_readNSS(nssmap,nss)
c ----------------------------------------------------------------------
c.... Purpose: Read nss from given file mapping
c....       nssmap.
c ----------------------------------------------------------------------
#ifndef _LINUX_
      USE KERNEL32
#endif
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: nss

      INTEGER(c_intptr_t), INTENT(IN) :: nssmap
      INTEGER(c_intptr_t) :: view
      INTEGER :: ret
#ifdef _LINUX_
      view = mmap(loc(0), 1*SIZEOF_INT, preadwrite, mshare, nssmap, 0)
      CALL memcpy( LOC(nss), view, 1*SIZEOF_INT )
      ret = munmap(view, 1*SIZEOF_INT )
#else
      view = MapViewOfFile( nssmap, FILE_MAP_ALL_ACCESS, 0, 0,
     &                      1*SIZEOF_INT )

      CALL memcpy( LOC(nss), view, 1*SIZEOF_INT )
      ret = UnmapViewOfFile(nssmap)
#endif

      END SUBROUTINE

      SUBROUTINE feap_shmem_writeE(emap,pe,pdt,pttim,dvp,skfy,skfz,iswm,
     &                             nss,ntyp)
c ----------------------------------------------------------------------
c.... Purpose: Write e(1,...,nss),dt,ttim,iswm into given file mapping
c....       emap.
c.... Enhanced: also write skfy, skfz (double), ntyp (int)
c.... further enhanced: write dvp (double)
c....
c.... Make sure to update mapping length when adding further variables!
c ----------------------------------------------------------------------
#ifndef _LINUX_
      USE KERNEL32
#endif
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: pdt, pttim, dvp, skfy, skfz
      REAL*8, DIMENSION(nss), INTENT(IN):: pe
      INTEGER, INTENT(IN) :: iswm, nss, ntyp
      INTEGER :: offset

      INTEGER(c_intptr_t), INTENT(IN) :: emap
      INTEGER(c_intptr_t) :: view
      INTEGER :: ret

#ifdef _LINUX_
      view = mmap(loc(0), (nss+5)*SIZEOF_REAL + 2*SIZEOF_INT, 
     &preadwrite, mshare, emap, 0)
#else
      view = MapViewOfFile( emap, FILE_MAP_ALL_ACCESS, 0, 0,
     &                      (nss+5)*SIZEOF_REAL + 2*SIZEOF_INT )
#endif

      offset = 0
c...    copy data to memory

      ! E
      CALL memcpy( view+offset,       ! c-pointer to destination
     &                 LOC(pe),           ! pointer to source
     &                 nss*SIZEOF_REAL    ! bytes to copy
     &                )
      offset = offset + nss*SIZEOF_REAL

      ! dt
      CALL memcpy( view+offset,       ! c-pointer to destination
     &                 LOC(pdt),          ! pointer to source
     &                 1*SIZEOF_REAL      ! bytes to copy
     &                )
      offset = offset + 1*SIZEOF_REAL

      ! ttim
      CALL memcpy( view+offset,       ! c-pointer to destination
     &                 LOC(pttim),        ! pointer to source
     &                 1*SIZEOF_REAL      ! bytes to copy
     &                )
      offset = offset + 1*SIZEOF_REAL

      ! iswm
      CALL memcpy( view+offset,       ! c-pointer to destination
     &                 LOC(iswm),         ! pointer to source
     &                 1*SIZEOF_INT       ! bytes to copy
     &                )
      offset = offset + 1*SIZEOF_INT

      ! ntyp
      CALL memcpy( view+offset,       ! c-pointer to destination
     &                 LOC(ntyp),         ! pointer to source
     &                 1*SIZEOF_INT       ! bytes to copy
     &                )
      offset = offset + 1*SIZEOF_INT

      ! dvp
      CALL memcpy( view+offset,       ! c-pointer to destination
     &                 LOC(dvp),          ! pointer to source
     &                 1*SIZEOF_REAL      ! bytes to copy
     &                )
      offset = offset + 1*SIZEOF_REAL

      ! skfy
      CALL memcpy( view+offset,       ! c-pointer to destination
     &                 LOC(skfy),         ! pointer to source
     &                 1*SIZEOF_REAL      ! bytes to copy
     &                )
      offset = offset + 1*SIZEOF_REAL

      ! skfz
      CALL memcpy( view+offset,       ! c-pointer to destination
     &                 LOC(skfz),         ! pointer to source
     &                 1*SIZEOF_REAL      ! bytes to copy
     &                )
      offset = offset + 1*SIZEOF_REAL

#ifdef _LINUX_
      ret = munmap(view, (nss+5)*SIZEOF_REAL + 2*SIZEOF_INT )
#else
      ret = UnmapViewOfFile(emap)
#endif

      END SUBROUTINE

      SUBROUTINE feap_shmem_readE(emap,pe,pdt,pttim,dvp,skfy,skfz,iswm,
     &                            nss,ntyp)
c ----------------------------------------------------------------------
c.... Purpose: Read e,dt,ttim,iswm from given file mapping
c....       emap.
c....   additional variables to read: dvp, skfy, skfz (all double), ntyp(int)
c ----------------------------------------------------------------------
#ifndef _LINUX_
      USE KERNEL32
#endif
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      REAL*8, INTENT(OUT) :: pdt, pttim, dvp, skfy, skfz
      REAL*8, DIMENSION(nss), INTENT(OUT):: pe
      INTEGER, INTENT(OUT) :: iswm, ntyp
      INTEGER, INTENT(IN)  :: nss
      INTEGER :: offset

      INTEGER(c_intptr_t), INTENT(IN) :: emap
      INTEGER(c_intptr_t) :: view
      INTEGER :: ret

#ifdef _LINUX_
      view = mmap(loc(0), (nss+5)*SIZEOF_REAL + 2*SIZEOF_INT, 
     &preadwrite, mshare, emap, 0)
#else
      view = MapViewOfFile( emap, FILE_MAP_ALL_ACCESS, 0, 0,
     &                      (nss+5)*SIZEOF_REAL + 2*SIZEOF_INT )
#endif

c...    copy from memory
      offset = 0

      ! E
      CALL memcpy( LOC(pe),           ! c-pointer to destination
     &                 view+offset,       ! pointer to source
     &                 nss*SIZEOF_REAL    ! bytes to copy
     &                )
      offset = offset + nss*SIZEOF_REAL

      ! dt
      CALL memcpy( LOC(pdt),          ! c-pointer to destination
     &                 view+offset,       ! pointer to source
     &                 1*SIZEOF_REAL      ! bytes to copy
     &                )
      offset = offset + 1*SIZEOF_REAL

      ! ttim
      CALL memcpy( LOC(pttim),        ! c-pointer to destination
     &                 view+offset,       ! pointer to source
     &                 1*SIZEOF_REAL      ! bytes to copy
     &                )
      offset = offset + 1*SIZEOF_REAL

      ! iswm
      CALL memcpy( LOC(iswm),         ! c-pointer to destination
     &                 view+offset,       ! pointer to source
     &                 1*SIZEOF_INT       ! bytes to copy
     &                )
      offset = offset + 1*SIZEOF_INT

      ! ntyp
      CALL memcpy( LOC(ntyp),         ! c-pointer to destination
     &                 view+offset,       ! pointer to source
     &                 1*SIZEOF_INT       ! bytes to copy
     &                )
      offset = offset + 1*SIZEOF_INT

      ! dvp
      CALL memcpy( LOC(dvp),          ! c-pointer to destination
     &                 view+offset,       ! pointer to source
     &                 1*SIZEOF_REAL      ! bytes to copy
     &                )
      offset = offset + 1*SIZEOF_REAL

      ! skfy
      CALL memcpy( LOC(skfy),         ! c-pointer to destination
     &                 view+offset,       ! pointer to source
     &                 1*SIZEOF_REAL      ! bytes to copy
     &                )
      offset = offset + 1*SIZEOF_REAL

      ! skfz
      CALL memcpy( LOC(skfz),         ! c-pointer to destination
     &                 view+offset,       ! pointer to source
     &                 1*SIZEOF_REAL      ! bytes to copy
     &                )
      offset = offset + 1*SIZEOF_REAL

#ifdef _LINUX_
      ret = munmap(view, (nss+5)*SIZEOF_REAL + 2*SIZEOF_INT )
#else
      ret = UnmapViewOfFile(emap)
#endif

      END SUBROUTINE

      SUBROUTINE feap_shmem_writeCS(csmap,psig,pcmat,nss)
c ----------------------------------------------------------------------
c.... Purpose: Write sig, cmat into given file mapping
c....       csmap.
c ----------------------------------------------------------------------
#ifndef _LINUX_
      USE KERNEL32
#endif
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nss

      REAL*8, DIMENSION(nss), INTENT(IN):: psig
      REAL*8, DIMENSION(nss,nss), INTENT(IN) :: pcmat
      REAL*8, DIMENSION(nss*nss) :: p1d_real_buf
      INTEGER :: offset,i,j

      INTEGER(c_intptr_t), INTENT(IN) :: csmap
      INTEGER(c_intptr_t) :: view
      INTEGER :: ret

#ifdef _LINUX_
      view = mmap(loc(0), (nss*nss+nss)*SIZEOF_REAL, 
     &preadwrite, mshare, csmap, 0)
#else
      view = MapViewOfFile( csmap, FILE_MAP_ALL_ACCESS, 0, 0,
     &                      (nss*nss+nss)*SIZEOF_REAL )
#endif

c...    copy C to memory-continous 1d array
      DO i=0,nss-1
          DO j=1,nss
              p1d_real_buf(i*nss+j) = pcmat(i+1,j)
          END DO
      END DO

      offset = 0
      CALL memcpy( view+offset,       ! c-pointer to destination
     &                 LOC(psig),         ! pointer to source
     &                 nss*SIZEOF_REAL    ! bytes to copy
     &                )
      offset = offset + nss*SIZEOF_REAL

      CALL memcpy( view+offset,       ! c-pointer to destination
     &                 LOC(p1d_real_buf), ! pointer to source
     &                 (nss*nss)*SIZEOF_REAL  ! bytes to copy
     &                )
      offset = offset + (nss*nss)*SIZEOF_REAL

#ifdef _LINUX_
      ret = munmap(view, (nss*nss+nss)*SIZEOF_REAL )
#else
      ret = UnmapViewOfFile(csmap)
#endif

      END SUBROUTINE

      SUBROUTINE feap_shmem_readCS(csmap,psig,pcmat,nss)
c ----------------------------------------------------------------------
c.... Purpose: Read sig, cmat from given file mapping
c....       csmap.
c ----------------------------------------------------------------------
#ifndef _LINUX_
      USE KERNEL32
#endif
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nss

      REAL*8, DIMENSION(nss), INTENT(OUT):: psig
      REAL*8, DIMENSION(nss,nss), INTENT(OUT) :: pcmat
      REAL*8, DIMENSION(nss*nss) :: p1d_real_buf
      INTEGER :: offset,i,j

      INTEGER(c_intptr_t), INTENT(IN) :: csmap
      INTEGER(c_intptr_t) :: view
      INTEGER :: ret

#ifdef _LINUX_
      view = mmap(loc(0), (nss*nss+nss)*SIZEOF_REAL, 
     &preadwrite, mshare, csmap, 0)
#else
      view = MapViewOfFile( csmap, FILE_MAP_ALL_ACCESS, 0, 0,
     &                      (nss*nss+nss)*SIZEOF_REAL )
#endif

      offset = 0
      CALL memcpy( LOC(psig),         ! c-pointer to destination
     &                 view+offset,       ! pointer to source
     &                 nss*SIZEOF_REAL      ! bytes to copy
     &                )
      offset = offset + nss*SIZEOF_REAL

      CALL memcpy( LOC(p1d_real_buf), ! c-pointer to destination
     &                 view+offset,       ! pointer to source
     &                 (nss*nss)*SIZEOF_REAL  ! bytes to copy
     &                )
      offset = offset + (nss*nss)*SIZEOF_REAL

c...    copy from 1d to 2d array
      DO i=0,nss-1
          DO j=1,nss
              pcmat(i+1,j) = p1d_real_buf(nss*i+j)
          END DO
      END DO

#ifdef _LINUX_
      ret = munmap(view, (nss*nss+nss)*SIZEOF_REAL )
#else
      ret = UnmapViewOfFile(csmap)
#endif

      END SUBROUTINE


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccc ab hier nur memcpy aufrufe  ccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE feap_shmem_writeGenValsDispl(view,offset,lo,numnp,
     &      numel,nummat,ndm,ndf,nrt,fl9,ttim,b,nneq,iu)
c ----------------------------------------------------------------------
c.... Purpose: Write general values and displacements to restart mapping
c....   with offset 'offset'.
c....
c....  Written data in this order:
c....   numnp, numel, nummat, ndm, ndf, nrt, fl9, ttim, b(1,nneq*iu)
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      INTEGER, INTENT(IN):: numnp,numel,nummat,ndm,ndf,nrt
      INTEGER, INTENT(IN):: nneq,iu
      REAL*8, INTENT(IN) :: ttim
      REAL*8, DIMENSION(nneq*iu), INTENT(IN) :: b
      LOGICAL, INTENT(IN):: fl9
      INTEGER :: i, len

      lo = 0

      i = 1
      len = (nneq*iu - i + 1) * SIZEOF_REAL

      call memcpy(view+offset+lo, LOC(numnp), SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(view+offset+lo, LOC(numel), SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(view+offset+lo, LOC(nummat), SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(view+offset+lo, LOC(ndm), SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(view+offset+lo, LOC(ndf), SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(view+offset+lo, LOC(nrt), SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(view+offset+lo, LOC(fl9), SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(view+offset+lo, LOC(ttim), SIZEOF_REAL)
      lo = lo + SIZEOF_REAL

      call memcpy(view+offset+lo, LOC(b(i)), len)
      lo = lo + len

      END SUBROUTINE

      SUBROUTINE feap_shmem_readGenValsDispl(view,offset,lo,nnpo,
     &      nnlo,nnmo,ndmo,ndfo,nrt,fl9,ttim,b,nneq,iu)
c ----------------------------------------------------------------------
c.... Purpose: Read general values and displacements from restart mapping
c....   with offset 'offset'.
c....
c....  Read data in this order:
c....   nnpo,nnlo,nnmo,ndmo,ndfo,nrt,fl9,ttim,b(1,nneq*iu)
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      INTEGER, INTENT(OUT):: nnpo,nnlo,nnmo,ndmo,ndfo,nrt
      INTEGER, INTENT(IN):: nneq,iu
      REAL*8, INTENT(OUT) :: ttim
      REAL*8, DIMENSION(nneq*iu), INTENT(OUT) :: b
      LOGICAL, INTENT(OUT):: fl9
      INTEGER :: i, len

      lo = 0

      i = 1
      len = (nneq*iu - i + 1) * SIZEOF_REAL

      call memcpy(LOC(nnpo), view+offset+lo, SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(LOC(nnlo), view+offset+lo, SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(LOC(nnmo), view+offset+lo, SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(LOC(ndmo), view+offset+lo, SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(LOC(ndfo), view+offset+lo, SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(LOC(nrt), view+offset+lo, SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(LOC(fl9), view+offset+lo, SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(LOC(ttim), view+offset+lo, SIZEOF_INT)
      lo = lo + SIZEOF_REAL

      call memcpy(LOC(b(i)), view+offset+lo, len)
      lo = lo + len

      END SUBROUTINE

      SUBROUTINE feap_shmem_readArcLenVals(view,offset,lo,prop,rlnew,
     &          c0,cs1,cs2,ds0,r,det0,xn)
c ----------------------------------------------------------------------
c.... Purpose: Read general values and displacements from restart mapping
c....   with offset 'offset'.
c....
c....  Read data in this order:
c....   prop, rlnew, c0, cs1, cs2, ds0, r, det0, xn
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      REAL*8, INTENT(OUT) :: prop,rlnew,c0,cs1,cs2,ds0,r,det0,xn

      lo = 0

      call memcpy(LOC(prop), view+offset+lo, SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(LOC(rlnew), view+offset+lo, SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(LOC(c0), view+offset+lo, SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(LOC(cs1), view+offset+lo, SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(LOC(cs2), view+offset+lo, SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(LOC(ds0), view+offset+lo, SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(LOC(r), view+offset+lo, SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(LOC(det0), view+offset+lo, SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(LOC(xn), view+offset+lo, SIZEOF_REAL)
      lo = lo + SIZEOF_REAL

      END SUBROUTINE

      SUBROUTINE feap_shmem_writeArcLenVals(view,offset,lo,prop,rlnew,
     &          c0,cs1,cs2,ds0,r,det0,xn)
c ----------------------------------------------------------------------
c.... Purpose: Write general values and displacements to restart mapping
c....   with offset 'offset'.
c....
c....  Written data in this order:
c....   prop, rlnew, c0, cs1, cs2, ds0, r, det0, xn
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      REAL*8, INTENT(IN) :: prop,rlnew,c0,cs1,cs2,ds0,r,det0,xn

      lo = 0

      call memcpy(view+offset+lo, LOC(prop), SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(view+offset+lo, LOC(rlnew), SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(view+offset+lo, LOC(c0), SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(view+offset+lo, LOC(cs1), SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(view+offset+lo, LOC(cs2), SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(view+offset+lo, LOC(ds0), SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(view+offset+lo, LOC(r), SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(view+offset+lo, LOC(det0), SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(view+offset+lo, LOC(xn), SIZEOF_REAL)
      lo = lo + SIZEOF_REAL

      END SUBROUTINE

      SUBROUTINE feap_shmem_readTransFields(view,offset,lo,m,nv,nrt,
     &              nneq,ipr)
c ----------------------------------------------------------------------
c.... Purpose: Read transient fields from restart mapping
c....   with offset 'offset'.
c....
c....  Read data in this order:
c....   m(nv, nv+nrt*nneq*ipr)
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      INTEGER, INTENT(IN) :: nv,nrt,nneq,ipr
      REAL*8, DIMENSION(*), INTENT(OUT):: m
      INTEGER :: i, len

      lo = 0

      i = 1 ! start offset
      len = nrt*nneq * SIZEOF_REAL ! length of field

      call memcpy(LOC(m(i)), view+offset+lo, len)
      lo = lo + len

      END SUBROUTINE

      SUBROUTINE feap_shmem_writeTransFields(view,offset,lo,m,nv,nrt,
     &              nneq,ipr)
c ----------------------------------------------------------------------
c.... Purpose: Write transient fields to restart mapping
c....   with offset 'offset'.
c....
c....  Written data in this order:
c....   m(nv, nv+nrt*nneq*ipr)
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      INTEGER, INTENT(IN) :: nv,nrt,nneq,ipr
      REAL*8, DIMENSION(*), INTENT(INOUT):: m
      INTEGER :: i, len

      lo = 0

      i = 1 ! start offset
      len = nrt*nneq * SIZEOF_REAL ! length of field

      call memcpy(view+offset+lo, LOC(m(i)), len)
      lo = lo + len

      END SUBROUTINE

      SUBROUTINE feap_shmem_readHistFields(view,offset,lo,ttim)
c ----------------------------------------------------------------------
c.... Purpose: Read history fields from restart mapping
c....   with offset 'offset'.
c....
c....  Read data in this order:
c....   ttim, lgh1,lgh3, gh1,gh2,gh3
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      USE hdata
      USE doalloc
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      REAL*8, INTENT(OUT) :: ttim
      INTEGER :: sgh1, sgh3
      INTEGER :: i, len
      logical ldummy

      lo = 0


      call memcpy(LOC(ttim), view+offset+lo, SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(LOC(sgh1), view+offset+lo, SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(LOC(sgh3), view+offset+lo, SIZEOF_INT)
      lo = lo + SIZEOF_INT

      call ralloc(gh1,sgh1,'History-NH1',ldummy)
      call ralloc(gh2,sgh1,'History-NH2',ldummy)
      call ralloc(gh3,sgh3,'History-NH3',ldummy)

      len = sgh1*SIZEOF_REAL
      call memcpy(LOC(gh1), view+offset+lo, len)
      lo = lo + len
      len = sgh3*SIZEOF_REAL
      call memcpy(LOC(gh3), view+offset+lo, len)
      lo = lo + len

      END SUBROUTINE

      SUBROUTINE feap_shmem_writeHistFields(view,offset,lo,ttim)
c ----------------------------------------------------------------------
c.... Purpose: Write history fields to restart mapping
c....   with offset 'offset'.
c....
c....  Written data in this order:
c....   ttim, lgh1,lgh3, gh1,gh2,gh3
c ----------------------------------------------------------------------
c     USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      USE hdata
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      REAL*8, INTENT(IN) :: ttim
      INTEGER :: sgh1,sgh3
      INTEGER :: i, len, j

      lo = 0
c      write(*,*) 'test'
c      pause
      sgh1 = size(gh1)
      sgh3 = size(gh3)


      call memcpy(view+offset+lo, LOC(ttim), SIZEOF_REAL)
      lo = lo + SIZEOF_REAL
      call memcpy(view+offset+lo, LOC(sgh1), SIZEOF_INT)
      lo = lo + SIZEOF_INT
      call memcpy(view+offset+lo, LOC(sgh3), SIZEOF_INT)
      lo = lo + SIZEOF_INT
      len = sgh1 * SIZEOF_REAL
      call memcpy(view+offset+lo, LOC(gh1), len)
      lo = lo + len
      len = sgh3 * SIZEOF_REAL
      call memcpy(view+offset+lo, LOC(gh3), len)
      lo = lo + len

      END SUBROUTINE

      SUBROUTINE feap_shmem_readEdgeData(view,offset,lo,ne5,ii,m)
c ----------------------------------------------------------------------
c.... Purpose: Read edge data from restart mapping
c....   with offset 'offset'.
c....
c....  Read data in this order:
c....   m(ne5, ne5+ii)
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      INTEGER, INTENT(IN) :: ne5, ii
      INTEGER, DIMENSION(*), INTENT(INOUT):: m
      INTEGER :: i, len

      lo = 0

      i = ne5 ! start offset
      len = (ii-ne5+1) * SIZEOF_INT ! length of field

      call memcpy(LOC(m(i)), view+offset+lo, len)
      lo = lo + len

      END SUBROUTINE

      SUBROUTINE feap_shmem_writeEdgeData(view,offset,lo,ne5,ii,m)
c ----------------------------------------------------------------------
c.... Purpose: Write edge data to restart mapping
c....   with offset 'offset'.
c....
c....  Written data in this order:
c....   m(ne5, ne5+ii)
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      INTEGER, INTENT(IN) :: ne5, ii
      INTEGER, DIMENSION(*), INTENT(IN):: m
      INTEGER :: i, len

      lo = 0

      i = ne5 ! start offset
      len = (ii-ne5+1) * SIZEOF_INT ! length of field

      call memcpy(view+offset+lo, LOC(m(i)), len)
      lo = lo + len

      END SUBROUTINE

      SUBROUTINE feap_shmem_readSingleInt(view,offset,lo,iint)
c ----------------------------------------------------------------------
c.... Purpose: Read single integer value from restart mapping
c....   with offset 'offset'.
c....
c....  Read data in this order:
c....   iint
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      INTEGER, INTENT(OUT) :: iint

      lo = 0

      call memcpy(LOC(iint), view+offset+lo, SIZEOF_INT)
      lo = lo + SIZEOF_INT

      END SUBROUTINE

      SUBROUTINE feap_shmem_writeSingleInt(view,offset,lo,iint)
c ----------------------------------------------------------------------
c.... Purpose: Write iint to restart mapping
c....   with offset 'offset'.
c....
c....  Written data in this order:
c....   iint
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      INTEGER, INTENT(IN) :: iint

      lo = 0

      call memcpy(view+offset+lo, LOC(iint), SIZEOF_INT)
      lo = lo + SIZEOF_INT

      END SUBROUTINE

      SUBROUTINE feap_shmem_readCrackValues(view,offset,lo,nc1,
     &              ncmax,m)
c ----------------------------------------------------------------------
c.... Purpose: Read crack values from restart mapping
c....   with offset 'offset'.
c....
c....  Read data in this order:
c....   m(nc1, ncmax)
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      INTEGER, INTENT(IN) :: nc1, ncmax
      INTEGER, DIMENSION(*), INTENT(INOUT):: m
      INTEGER :: i, len

      lo = 0

      i = nc1
      len = (ncmax - nc1 + 1)*SIZEOF_INT

      call memcpy(LOC(m(i)), view+offset+lo, len)
      lo = lo + len

      END SUBROUTINE

      SUBROUTINE feap_shmem_writeCrackValues(view,offset,lo,nc1,
     &              ncmax,m)
c ----------------------------------------------------------------------
c.... Purpose: Write crack values to restart mapping
c....   with offset 'offset'.
c....
c....  Written data in this order:
c....   m(nc1, ncmax)
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      INTEGER, INTENT(IN) :: nc1, ncmax
      INTEGER, DIMENSION(*), INTENT(IN):: m
      INTEGER :: i, len

      lo = 0

      i = nc1
      len = (ncmax - nc1 + 1)*SIZEOF_INT

      call memcpy(view+offset+lo, LOC(m(i)), len)
      lo = lo + len

      END SUBROUTINE

      SUBROUTINE feap_shmem_readDirectorVals(view,offset,lo,mdir,
     &              mdirmax,m)
c ----------------------------------------------------------------------
c.... Purpose: Read director values from restart mapping
c....   with offset 'offset'.
c....
c....  Read data in this order:
c....   m(mdir, mdirmax)
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      INTEGER, INTENT(IN) :: mdir, mdirmax
      INTEGER, DIMENSION(*), INTENT(INOUT):: m
      INTEGER :: i, len

      lo = 0

      i = mdir
      len = (mdirmax - mdir + 1)*SIZEOF_INT

      call memcpy(LOC(m(i)), view+offset+lo, len)
      lo = lo + len

      END SUBROUTINE

      SUBROUTINE feap_shmem_writeDirectorVals(view,offset,lo,mdir,
     &              mdirmax,m)
c ----------------------------------------------------------------------
c.... Purpose: Write director values to restart mapping
c....   with offset 'offset'.
c....
c....  Written data in this order:
c....   m(mdir, mdirmax)
c ----------------------------------------------------------------------
c      USE KERNEL32
      use, intrinsic :: iso_c_binding
      use memorymap
      IMPLICIT NONE

      INTEGER(c_intptr_t), INTENT(IN)::view
      INTEGER, INTENT(IN) :: offset
      INTEGER, INTENT(OUT) :: lo ! local offset
      INTEGER, INTENT(IN) :: mdir, mdirmax
      INTEGER, DIMENSION(*), INTENT(IN):: m
      INTEGER :: i, len

      lo = 0

      i = mdir
      len = (mdirmax - mdir + 1)*SIZEOF_INT

      call memcpy(view+offset+lo, LOC(m(i)), len)
      lo = lo + len

      END SUBROUTINE


#endif