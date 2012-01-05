test_suite brillouin

test simple_sums
    use double
    use environment
    implicit none
    type(Environ) :: env
    real :: sum
    integer :: zoneLength
    zoneLength = 64
    env = NewEnv(zoneLength, 1.0_DP, 0.1_DP, 1.0_DP, 0.1_DP, 0.1_DP)
    sum = BrilSum(env, constantValue)
    assert_real_equal(real(zoneLength ** 3), real(BrilSum(env, constantValue)))
    assert_equal_within(0.0, real(BrilSum(env, sinXValue)), 1e-10)
    assert_equal_within(0.0, real(BrilSum(env, sinYValue)), 1e-10)
    assert_equal_within(0.0, real(BrilSum(env, sinZValue)), 1e-10)
end test

function constantValue(env, k)
    use double
    use environment
    implicit none
    type(Environ), intent(in) :: env
    real(kind=DP), intent(in) :: k(1:3)
    real(kind=DP) :: constantValue
    constantValue = 1.0
end function

function sinXValue(env, k)
    use double
    use environment
    implicit none
    type(Environ), intent(in) :: env
    real(kind=DP), intent(in) :: k(1:3)
    real(kind=DP) :: sinXValue
    sinXValue = sin(k(1))
end function

function sinYValue(env, k)
    use double
    use environment
    implicit none
    type(Environ), intent(in) :: env
    real(kind=DP), intent(in) :: k(1:3)
    real(kind=DP) :: sinYValue
    sinYValue = sin(k(2))
end function

function sinZValue(env, k)
    use double
    use environment
    implicit none
    type(Environ), intent(in) :: env
    real(kind=DP), intent(in) :: k(1:3)
    real(kind=DP) :: sinZValue
    sinZValue = sin(k(3))
end function

end test_suite
