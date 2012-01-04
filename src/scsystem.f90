module scsystem
    use double
    use environment
    use scsolve
    implicit none
    type :: SelfConsistentSystem
        real(kind=DP), dimension(:), allocatable :: tolerances
        type(SelfConsistentEq), dimension(:), allocatable :: equations
    contains
        procedure :: isSolved => scSystemIsSolved
        procedure :: solvedUpTo => scSystemSolvedUpTo
    end type
contains
    function scSystemSolvedUpTo(this, env, maxIndex) result(isSolved)
        class(SelfConsistentSystem), intent(in) :: this
        type(Environ), intent(in) :: env
        integer, intent(inout) :: maxIndex
        logical :: isSolved
        integer :: n, numEquations
        numEquations = size(this%equations)
        if (maxIndex > numEquations) then
            maxIndex = numEquations
        end if
        do n=1,maxIndex
            if (abs(this%equations(n)%absError(env)) > this%tolerances(n)) then
                isSolved = .false.
                return
            end if
        end do
        isSolved = .true.
    end function

    function scSystemIsSolved(this, env) result(isSolved)
        class(SelfConsistentSystem), intent(in) :: this
        type(Environ), intent(in) :: env
        logical :: isSolved
        integer :: numEquations
        numEquations = size(this%equations)
        if (numEquations == 0) then
            isSolved = .true.
        else
            isSolved = this%solvedUpTo(env, numEquations)
        end if
    end function

    function scSystemSolve(system, env) result(error)
        type(SelfConsistentSystem), intent(in) :: system
        type(Environ), intent(inout) :: env
        integer :: error, i
        i = 0
        do
            if (system%isSolved(env)) then
                exit
            end if
            error = scEqSolve(system%equations(i), env)
            if (error > 0) then
                return
            end if
            if (.not. system%solvedUpTo(env, i)) then
                i = 0
            else
                i = i + 1
            end if
        end do
    end function
end module scsystem
