test_suite scsolve

test linear_solve
    use double
    use environment
    implicit none
    type(Environ) :: env
    type(SelfConsistentEq) :: linearEq
    real(kind=DP) :: t, tc
    integer :: error
    t = 1.0_DP
    tc = 0.1_DP
    env = NewEnv(64, t, tc, 1.0_DP, 0.1_DP, 0.1_DP)
    linearEq%argMin = -5.0_DP * tc / t
    linearEq%argMax = 5.0_DP * tc / t
    linearEq%setArg => setD
    linearEq%absError => linearFuncD
    error = scEqSolve(linearEq, env)
    assert_equal(0, error)
    assert_real_equal(-tc / t, env%D)
end test

function linearFuncD(env)
    use double
    use environment
    implicit none
    type(Environ), intent(in) :: env
    real(kind=DP) :: linearFuncD
    linearFuncD = env%t * env%D + env%tc
end function

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

end test_suite
