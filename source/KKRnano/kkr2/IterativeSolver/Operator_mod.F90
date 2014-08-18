module Operator_mod
  implicit none

  type, abstract :: Operator
    contains
    procedure(apply), deferred :: apply
  end type

  interface
    subroutine apply(self)
      import Operator
      class(Operator) :: self
    end subroutine
  end interface

end module
