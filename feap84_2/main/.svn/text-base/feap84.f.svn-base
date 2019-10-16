c$Id:$
      program feap

c-----[--.----+----.----+----.-----------------------------------------]

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c     The University of California  does not guaranteed this program to
c     be error free. Every effort has been made to ensure proper coding
c     and declarations, but testing of all program options has not been
c     performed.

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    11/11/2008
c       1.  Version 8.4.1a                                  23/06/2013
c       2.  Version 8.4.1b                                  01/07/2013
c       3.  Version 8.4.1c                                  01/12/2013
c       3.  Version 8.4.1d                                  01/01/2014
c       4.  Version 8.4.1e                                  02/07/2014
c-----[--.----+----.----+----.-----------------------------------------]

c     Finite Element Analysis Program - (FEAP)  for solution of general
c     problem classes using the finite element method.  Problem size is
c     controlled by the amount of main memory in each computer.

c     Programmed by:
c                R. L. Taylor
c                Department of Civil and Environmental Engineering
c                University of California at Berkeley
c                Berkeley, California 94720-1710

c     E-mail:
c                feap@berkeley.edu

c     Web address:
c                http://www.ce.berkeley.edu/feap

c     User forum:
c                http://feap.berkeley.edu
c-----[--.----+----.----+----.-----------------------------------------]

c     Notes:

c     1. Precision is controlled by ipr:

c        Set ipr = 1 for 8-byte integers; = 2 for 4-byte integers.

c     2. Pointers used to address location of main arrays in memory.
c        Generally, for 32 bit processors compilation should be made
c        using common blocks defined in: ./include/integer4
c        For 64 bit processors compilation should be made using those
c        in: ./include/integer8
c        N.B. It is also necessary to always use those in: ./include

c     3. User written subprograms should include type specification
c        for all variables in each subprogram.

c        e.g.    real*8    a(12)
c                integer   name
c                character word*6
c                logical   flag
c                etc.

c     4. FEAP may create temporary input files during use.
c        Users should periodically check and delete files
c        which are no longer needed.  File names are normally
c        either the name of the data input file with an extender
c        or the name of the plot save file with an extender.

c     5. System dependent routines are included in:

c           a. File for subroutine doargs.

c           b. Plot routines include special characters which
c              must have case preserved (i.e., upper and lower
c              case letters in formats, etc.).


c     6. Input/Output is performed to files during execution of FEAP.
c        In general, the following files are used during executions:

c           a.  User logical unit numbers should be from 1 to 9.
c           b.  ilg = 10 : Used for write of log file.
c           c.  iop = 11 : Used for read/write delayed inputs.
c           d.  ios = 12 : Used for read/write scratch files.
c           e.  ird = 13 : Used to read results data from disk.
c           f.  iwd = 14 : Used to write results data to disk.
c           g.  ior = 15 : Use to read from the input data file.
c                          (specified when a problem is initiated).
c           h.  iow = 16 : Use to write output result data to file.
c                          (specified when a problem is initiated).
c                          (N.B. Multi-problems use iow = 8).
c           i.  lun = 17 : For PostScript file outputs.
c           j.  icl = 18+: Used for include file inputs.
c                          (additional include files may be opened).
c           k.  24+      : Used to save time history data.

c     End of Notes
c-----[--.----+----.----+----.-----------------------------------------]

c     Set variable types

      implicit none

      include 'cdata.h'
      include 'codat.h'
      include 'hlpdat.h'
      include 'iodata.h'
      include 'iofile.h'
      include 'prmptd.h'
      include 'psize.h'
      include 'pathn.h'
      include 'setups.h'
      include 'vdata.h'

c-----[--.----+----.----+----.-----------------------------------------]
c     Set version header for output to file and screen

      versn(1) = 'Release 8.4.1d'
      versn(2) = '01 January 2014'

c-----[--.----+----.----+----.-----------------------------------------]
c     Set ratio for real to integer variables: Set ipr = 1 or 2
                          ! ipr = 1 for equal  length real to integers
      ipr = 2             ! ipr = 2 for double length real to integers

c-----[--.----+----.----+----.-----------------------------------------]
c     Set default logical unit numbers for files

      ilg = 10
      iop = 11
      ios = 12
      ird = 13
      iwd = 14
      ior = 15
      iow = 16
      lun = 17
      icl = 18

c-----[--.----+----.----+----.-----------------------------------------]
c     Set data input parsing flag

      coflg = .true.  ! Parse all input as expressions   (slower mode)
                      ! N.B. Use of 'parse' and 'noparse' mesh commands
                      !      permit change of this default setting.
c     ciflg = .true.  ! Get input file from menu window if .true.
      ciflg = .false. ! or from keyboard entry if .false.
                      ! N.B. Used for Windows version only

c-----[--.----+----.----+----.-----------------------------------------]
c     Set graphics default options

      defalt = .true. ! Graphics runs with default contour interval
      prompt = .true. ! Prompt for graphics inputs when defalt = .false.

c-----[--.----+----.----+----.-----------------------------------------]
c     Set PostScript default mode

      pscolr = .true. ! PostScript outputs are in color
      psrevs = .false.! Color order is normal

c-----[--.----+----.----+----.-----------------------------------------]
c     Set help display level:

      hlplev = 0      ! Basic

c-----[--.----+----.----+----.-----------------------------------------]
c     Set solver flag: Program - solver= .true.; User - solver =.false.
c                      soltyp    = Profile solver type: 1 = by 1 column
c                                                       2 = by 2 column
      soltyp    = 2
      solver    = .true.

c-----START SOLUTION---------------------------------------------------]

c     Initialize solution system

      call pstart()

c     Display messages

      call pmessage()

c     Solve problem

      call pcontr()

c-----END SOLUTION-----------------------------------------------------]

c     N.B. No additional calls permitted after PCONTR

      end ! FEAP8.4
