module Precond_mod
  implicit none

  type, abstract :: Precond 
    contains
    procedure(apply), deferred :: apply
  end type

  interface
    subroutine apply(self)
      import Precond 
      class(Precond) :: self
    end subroutine
  end interface

end module
