c$Id:$
      subroutine doargs(inp,outp,res,sav,plt,narg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Replace MSFLIB by DFLIB                          08/12/2006
c       2. Add 'reflg' to perform restart on command line   12/10/2007
c          and error output.
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Use command line arguments to set file names for
c               input/output, restarts, and plots.

c      Inputs:
c         None

c      Outputs:
c         inp     - Name of input file
c         outp    - Name of output file
c         res     - Name of restart read file
c         sav     - Name of restart save file
c         plt     - Name of plot file
c         narg    - Number of items input from command line
c         reflg   - Restart file on command line
c-----[--.----+----.----+----.-----------------------------------------]
      use       DFLIB

      implicit  none

      include  'prmptd.h'

      character inp*128,outp*128,res*128,sav*128,plt*128,argv*130
      integer*2 i2, nchar
      integer   i , narg, nchars

      save

      narg  = nargs() - 1

      if(narg.gt.0) then

c       Set files to blank

        inp = ' '
        outp= ' '
        plt = ' '
        res = ' '
        sav = ' '

c       Check arguments set on command line

        do i = 1, narg

          i2 = i
          call getarg(i2,argv,nchar)
          nchars = nchar

c         Device specification

          if (argv(1:1) .eq. '-') then

c           Input file specification

            if      (argv(2:2).eq.'i') then
              inp = argv(3:nchars)

c           Output file specification

            else if (argv(2:2).eq.'o') then
              outp = argv(3:nchars)

c           Plot file specification

            else if (argv(2:2).eq.'p') then
              plt = argv(3:nchars)

c           Restart read file specification

            else if (argv(2:2).eq.'r') then
              res = argv(3:nchars)
              reflg = .true.

c           Restart save file specification

            else if (argv(2:2).eq.'s') then
              sav = argv(3:nchars)

c           Error on command line

            else
              write( *, 2000) argv(2:nchars)
              call plstop()
            endif
          else
            write( *, 2001)  argv(1:nchars)
            call plstop()
          endif
        end do ! i

c     Check that files are correct if narg > 0

        if(inp.ne.' ') then
          call dosets(inp,outp,res,sav,plt)
        else
          write( *, 2002)
          call plstop()
        endif
      endif

2000  format(' Unknown command line option   -> ',a)
2001  format(' Unknown command line argument -> ',a)
2002  format(' *ERROR*  Command line arguments must include an',
     &       ' INPUT filename')

      end
