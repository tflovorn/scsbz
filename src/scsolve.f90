module scsolve
    use double
    use environment
    implicit none
    type :: SelfConsistentEq
        real(kind=DP) :: argMin, argMax
        procedure(setArgInterface), pointer, nopass :: setArg
        procedure(absErrorInterface), pointer, nopass :: absError
    end type

    abstract interface
        function setArgInterface(env, x)
            use double
            use environment
            type(Environ), intent(in) :: env
            real(kind=DP), intent(in) :: x
            real(kind=DP) :: setArg
        end function
        function absErrorInterface(env)
            use double
            use environment
            type(Environ), intent(in) :: env
            real(kind=DP) :: absError
        end function
    end interface
contains

end module scsolve
