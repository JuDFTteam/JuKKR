    Real *8 Function rnuctab(z)
      Use mod_datatypes, Only: dp
      Integer, Intent (Out) :: z

      z = 0
      rnuctab = 0E0_dp
      Stop ' < RNUCTAB > : NUCLEUS <> 0 not implemented '
    End Function
