test_suite scsystem

test linear_system_solve
    use double
    use environment
    use scsolve
    implicit none
    type(Environ) :: env
    type(SelfConsistentEq) :: linearEqPlus, linearEqMinus
    real(kind=DP) :: t, tc, x, expectedD, expectedA, tolerance
    real(kind=DP), dimension(:), allocatable :: tolerances
    type(SelfConsistentEq), dimension(:), allocatable :: equations
    type(SelfConsistentSystem) :: system
    integer :: error
    t = 1.0_DP
    tc = 0.1_DP
    x = 0.1_DP
    expectedD = -(1.0_DP / 13.0_DP) * x / t
    expectedA = -(8.0_DP / 13.0_DP) * x / tc
    tolerance = 1e-10
    env = NewEnv(64, t, tc, 1.0_DP, x, 0.1_DP)
    linearEqPlus%argMin = -5.0_DP * abs(expectedD)
    linearEqPlus%argMax = 5.0_DP * abs(expectedD)
    linearEqPlus%setArg => setD
    linearEqPlus%absError => linearFuncPlus
    linearEqMinus%argMin = -5.0_DP * abs(expectedA)
    linearEqMinus%argMax = 5.0_DP * abs(expectedA)
    linearEqMinus%setArg => setA
    linearEqMinus%absError => linearFuncMinus
    allocate(tolerances(1:2))
    allocate(equations(1:2))
    tolerances(1) = tolerance
    tolerances(2) = tolerance
    equations(1) = linearEqPlus
    equations(2) = linearEqMinus
    system%tolerances = tolerances
    system%equations = equations
    error = scSystemSolve(system, env)
    assert_equal(0, error)
    assert_real_equal(env%D, expectedD)
    assert_real_equal(env%A, expectedA)
end test

function setD(env, D)
    use double
    use environment
    implicit none
    type(Environ), intent(out) :: env
    real(kind=DP), intent(in) :: D
    real(kind=DP) :: setD
    setD = env%D
    env%D = D
end function

function setA(env, A)
    use double
    use environment
    implicit none
    type(Environ), intent(out) :: env
    real(kind=DP), intent(in) :: A
    real(kind=DP) :: setA
    setA = env%A
    env%A = A
end function

function linearFuncPlus(env)
    use double
    use environment
    implicit none
    type(Environ), intent(in) :: env
    real(kind=DP) :: linearFuncPlus
    linearFuncPlus = 5.0_DP * env%t * env%D + env%tc * env%A + env%x
end function

function linearFuncMinus(env)
    use double
    use environment
    implicit none
    type(Environ), intent(in) :: env
    real(kind=DP) :: linearFuncMinus
    linearFuncMinus = env%t * env%D - 5.0_DP * env%tc * env%A - 3.0_DP * env%x
end function

end test_suite
