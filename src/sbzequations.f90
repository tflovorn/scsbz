module sbzequations
    use double
    use environment
    use brillouin
    use scsolve
    use scsystem
    implicit none
contains
    ! --- D equation ---
    function sbzEqD(env)
        type(Environ), intent(in) :: env
        type(SelfConsistentEq) :: sbzEqD
        ! (Mostly) arbitrary bounds (should D be positive?).
        sbzEqD%argMin = 1e-9
        sbzEqD%argMax = 10.0 * abs(env%t)
        sbzEqD%setArg => setD
        sbzEqD%absError => absErrorD
    end function

    function absErrorD(env)
        type(Environ), intent(in) :: env
        real(kind=DP) :: absErrorD, rhs
        rhs = BrilSum(env, sumFuncD) / (env%zoneLength ** 3)
        absErrorD = env%D - rhs
    end function

    function sumFuncD(env, k)
        type(Environ), intent(in) :: env
        real(kind=DP), intent(in) :: k(1:3)
        real(kind=DP) :: sumFuncD
        sumFuncD = cos(k(1)) * bose(env, xiB(env, k))
    end function

    function setD(env, D)
        type(Environ), intent(inout) :: env
        real(kind=DP), intent(in) :: D
        real(kind=DP) :: setD
        setD = env%D
        env%D = D
    end function

    ! --- B equation ---
    function sbzEqB(env)
        type(Environ), intent(in) :: env
        type(SelfConsistentEq) :: sbzEqB
        ! (Mostly) arbitrary bounds (should B be positive?).
        sbzEqB%argMin = 1e-9
        sbzEqB%argMax = 10.0 * abs(env%t)
        sbzEqB%setArg => setB
        sbzEqB%absError => absErrorB
    end function

    function absErrorB(env)
        type(Environ), intent(in) :: env
        real(kind=DP) :: absErrorB, rhs
        rhs = BrilSum(env, sumFuncB) / (env%zoneLength ** 3)
        absErrorB = env%B - rhs
    end function

    function sumFuncB(env, k)
        type(Environ), intent(in) :: env
        real(kind=DP), intent(in) :: k(1:3)
        real(kind=DP) :: sumFuncB
        sumFuncB = 0.5_DP * cos(k(1)) * ((xiF(env, k) / bogoEnergy(env, k)) * (2.0_DP * fermi(env, xiF(env, k)) - 1.0_DP) + 1.0_DP)
    end function

    function setB(env, B)
        type(Environ), intent(inout) :: env
        real(kind=DP), intent(in) :: B
        real(kind=DP) :: setB
        setB = env%B
        env%B = B
    end function

    ! --- A equation ---
    function sbzEqA(env)
        type(Environ), intent(in) :: env
        type(SelfConsistentEq) :: sbzEqA
        ! (Mostly) arbitrary bounds (should A be positive?).
        sbzEqA%argMin = 1e-9
        sbzEqA%argMax = 10.0 * abs(env%t)
        sbzEqA%setArg => setA
        sbzEqA%absError => absErrorA
    end function

    function absErrorA(env)
        type(Environ), intent(in) :: env
        real(kind=DP) :: absErrorA, rhs
        rhs = BrilSum(env, sumFuncA) / (env%zoneLength ** 3)
        absErrorA = env%A - rhs
    end function

    function sumFuncA(env, k)
        type(Environ), intent(in) :: env
        real(kind=DP), intent(in) :: k(1:3)
        real(kind=DP) :: sumFuncA
        sumFuncA = 0.5_DP * cos(k(1)) * (Delta(env, k) / bogoEnergy(env, k)) * (2.0_DP * fermi(env, xiF(env, k)) - 1.0_DP)
    end function

    function setA(env, A)
        type(Environ), intent(inout) :: env
        real(kind=DP), intent(in) :: A
        real(kind=DP) :: setA
        setA = env%A
        env%A = A
    end function

    ! --- muB equation ---
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

    ! --- muF equation ---
    function sbzEqMuF(env)
        type(Environ), intent(in) :: env
        type(SelfConsistentEq) :: sbzEqMuF
        ! (Mostly) arbitrary bounds (muF must be positive and nonzero).
        sbzEqMuF%argMin = 10.0 * abs(env%t)
        sbzEqMuF%argMax = 1e-9
        sbzEqMuF%setArg => setMuF
        sbzEqMuF%absError => absErrorMuF
    end function

    function absErrorMuF(env)
        type(Environ), intent(in) :: env
        real(kind=DP) :: absErrorMuF, rhs
        rhs = BrilSum(env, sumFuncMuF) / (env%zoneLength ** 3)
        absErrorMuF = 1.0_DP - env%x - rhs
    end function

    function sumFuncMuF(env, k)
        type(Environ), intent(in) :: env
        real(kind=DP), intent(in) :: k(1:3)
        real(kind=DP) :: sumFuncMuF
        sumFuncMuF = (xiF(env, k) / bogoEnergy(env, k)) * (2.0_DP * fermi(env, xiF(env, k)) - 1.0_DP) + 1.0_DP
    end function

    function setMuF(env, muF)
        type(Environ), intent(inout) :: env
        real(kind=DP), intent(in) :: muF
        real(kind=DP) :: setMuF
        setMuF = env%muF
        env%muF = muF
    end function

    ! --- Distribution functions. ---
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
