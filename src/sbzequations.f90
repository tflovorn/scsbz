module sbzequations
    use double
    use environment
    use brillouin
    use scsolve
    use scsystem
    implicit none
contains
    ! Fermi distribution function.
    function fermi(env, energy)
        type(Environ), intent(in) :: env
        real(kind=DP), intent(in) :: energy
        real(kind=DP) :: fermi
        fermi = 1.0_DP / (exp(env%beta * energy) + 1.0_DP)
    end function

    ! Bose distribution function.
    function bose(env, energy)
        type(Environ), intent(in) :: env
        real(kind=DP), intent(in) :: energy
        real(kind=DP) :: bose
        bose = 1.0_DP / (exp(env%beta * energy) - 1.0_DP)
    end function
end module
