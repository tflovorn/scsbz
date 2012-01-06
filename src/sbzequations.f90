module sbzequations
    use double
    use environment
    use brillouin
    use scsolve
    use scsystem
    implicit none
contains
    function sbzEqMuB(env)
        type(Environ), intent(in) :: env
        type(SelfConsistentEq) :: sbzEqMuB
        ! (Mostly) arbitrary bounds (muB must be negative and nonzero).
        sbzEqMuB%argMin = -10.0 * abs(env%t)
        sbzEqMuB%argMax = -1e-9
        sbzEqMuB%setArg => setMuB
        sbzEqMuB%absError => absErrorMuB
    end function

    function absErrorMuB(env)
        type(Environ), intent(in) :: env
        real(kind=DP) :: absErrorMuB, rhs
        rhs = BrilSum(env, sumFuncMuB) / (env%zoneLength ** 3)
        absErrorMuB = env%x - rhs
    end function

    function sumFuncMuB(env, k)
        type(Environ), intent(in) :: env
        real(kind=DP), intent(in) :: k(1:3)
        real(kind=DP) :: sumFuncMuB
        sumFuncMuB = bose(env, xiB(env, k))
    end function

    function setMuB(env, muB)
        type(Environ), intent(inout) :: env
        real(kind=DP), intent(in) :: muB
        real(kind=DP) :: setMuB
        setMuB = env%muB
        env%muB = muB
    end function

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
