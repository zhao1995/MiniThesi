
                               ! CONTACT PARAMETERS

      integer    c_nr0,c_ncc,c_ncs,c_ncel,c_lp1,c_lp3,c_lmv

      parameter (c_nr0 = 18)   ! # of available rows in all tables
      parameter (c_ncc = 12)   ! # of available contact commands
      parameter (c_ncs = 200)  ! # of available command strings
      parameter (c_ncel= 26)   ! # of available contact elements
      parameter (c_lp1 = 200)  ! # of available history variables defs
                               !      for vectors CH1 & CH2
      parameter (c_lp3 = 100)  ! # of available history variables defs
                               !      for vectors CH3
      parameter (c_lmv = 50)   ! # of available material variables
