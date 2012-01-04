module scsystem
    use double
    use environment
    use scsolve
    implicit none
    type :: SelfConsistentSystem
        real(kind=DP), dimension(:), allocatable :: tolerance
        type(SelfConsistentEq), dimension(:), allocatable :: equations
    end type
contains
    function scSystemSolve(system, env) result(error)
        type(SelfConsistentSystem), intent(in) :: system
        type(Environ), intent(inout) :: env
        integer :: error
    end function
end module scsystem
