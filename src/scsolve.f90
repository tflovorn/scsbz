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
        ! Set the variable in env which corresponds to this equation to x.
        ! Return the value of the equation's argument in env before it is set.
        function setArgInterface(env, x)
            use double
            use environment
            implicit none
            type(Environ), intent(out) :: env
            real(kind=DP), intent(in) :: x
            real(kind=DP) :: setArgInterface
        end function
        ! Return the (signed) absolute error of the equation using env.
        function absErrorInterface(env)
            use double
            use environment
            implicit none
            type(Environ), intent(in) :: env
            real(kind=DP) :: absErrorInterface
        end function
    end interface
contains
    ! Return the (signed) absolute error of this equation, using env for all
    ! arguments except the one corresponding to this equation. Set that
    ! argument to env, evaluate the error, then replace the argument with its
    ! original value.
    function absErrorWith(scEquation, env, x)
        type(SelfConsistentEq), intent(in) :: scEquation
        type(Environ), intent(inout) :: env
        real(kind=DP), intent(in) :: x
        real(kind=DP) :: absErrorWith, orig
        orig = scEquation%setArg(env, x)
        absErrorWith = scEquation%absError(env)
        ! Shouldn't really be assigning here, but need to for compilation.
        ! Assignment doesn't affect result.
        orig = scEquation%setArg(env, orig)
    end function

    ! Modify env to minimize the magnitude of absError(env), solving
    ! scEquation. Returns 1 if no solution could be found, 0 otherwise.
    function scEqSolve(scEquation, env) result(error)
        type(SelfConsistentEq), intent(in) :: scEquation
        type(Environ), intent(inout) :: env
        integer :: error
        real(kind=DP) :: leftBracket, rightBracket, solution
        call BracketRoot(scEquation, env, leftBracket, rightBracket, error)
        if (error > 0) then
            return
        end if
        solution = BisectionRoot(scEquation, env, leftBracket, rightBracket,&
                                 error)
        if (error > 0) then
            return
        end if
        ! Shouldn't really assign here, same as in scEqAbsErrorWith.
        solution = scEquation%setArg(env, solution)
    end function

    ! Return a bracket for the leftmost root of scEquation between argMin
    ! and argMax. If error = 1, no root was found.
    subroutine BracketRoot(scEquation, env, leftBracket, rightBracket, error)
        real(kind=DP), intent(out) :: leftBracket, rightBracket
        type(SelfConsistentEq), intent(in) :: scEquation
        type(Environ) :: env
        integer, intent(out) :: error
        integer :: minSteps, maxSteps, numSteps, n
        real(kind=DP) :: stepLength, left, right, a, b, leftBound, rightBound
        if (scEquation%argMax < scEquation%ArgMin) then
            leftBound = scEquation%argMax
            rightBound = scEquation%argMin
        else
            leftBound = scEquation%argMin
            rightBound = scEquation%argMax
        end if
        ! minSteps and maxSteps are arbitary
        ! (could reduce number of fn calls by memoizing already-seen values)
        minSteps = 1
        maxSteps = minSteps * 256
        numSteps = minSteps
        varyStepSize: do
            stepLength = (rightBound - leftBound) / numSteps
            left = leftBound
            right = leftBound + stepLength
            varyFunctionArgs: do n=1, numSteps
                a = absErrorWith(scEquation, env, left)
                b = absErrorWith(scEquation, env, right)
                if (a == 0 .or. b == 0 .or. (a > 0 .and. b < 0) .or. &
                  (a < 0 .and. b > 0)) then
                    error = 0
                    leftBracket = left
                    rightBracket = right
                    return
                end if
                left = right
                right = right + stepLength
            end do varyFunctionArgs
            numSteps = numSteps * 2
            if (numSteps > maxSteps) then
                exit
            end if
        end do varyStepSize
        ! if we get here, overflowed step size limit
        error = 1
    end subroutine

    ! Find a root of targetFunc between leftBound and rightBound by bisection.
    ! Assumes that leftBracket and rightBracket bracket a root. If error = 1, 
    ! no root is bracketed. Assumes that targetFunc is continuous.
    ! Looks for a zero to an accuracy equal to the machine precision.
    function BisectionRoot(scEquation, env, leftBracket, rightBracket, error)
        real(kind=DP), intent(in) :: leftBracket, rightBracket
        integer, intent(out) :: error
        type(SelfConsistentEq), intent(in) :: scEquation
        type(Environ) :: env
        real(kind=DP) :: BisectionRoot, a, b, low, high, mid
        a = absErrorWith(scEquation, env, leftBracket)
        b = absErrorWith(scEquation, env, rightBracket)
        if (.not. ((a >= 0 .and. b <= 0) .or. (a <= 0 .and. b >= 0))) then
            error = 1
            return
        end if
        if (a <= 0) then
            low = leftBracket
            high = rightBracket
        else
            low = rightBracket
            high = leftBracket
        end if
        mid = low + (high - low) / 2.0_DP
        do
            if ((mid == low) .or. (mid == high)) then
                error = 0
                BisectionRoot = mid
                return
            end if
            if (absErrorWith(scEquation, env, mid) <= 0) then
                low = mid
            else
                high = mid
            end if
            mid = low + (high - low) / 2.0_DP
        end do
    end function
end module scsolve
