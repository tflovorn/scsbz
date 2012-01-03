module rootfinder
    use double
contains
    ! Return a bracket for the leftmost root of targetFunc between leftBound
    ! and rightBound. If error = 1, no root was found. Must have rightBound
    ! greater than leftBound; if this is violated, error = 2.
    subroutine BracketRoot(targetFunc, leftBound, rightBound, leftBracket,&
                           rightBracket, error)
        ! minSteps and maxSteps are arbitary
        ! (could reduce number of fn calls by memoizing already-seen values)
        real(kind=DP), intent(in) :: leftBound, rightBound
        real(kind=DP), intent(out) :: leftBracket, rightBracket
        interface
            function targetFunc(x)
                use double
                real(kind=DP), intent(in) :: x
                real(kind=DP) :: targetFunc
            end function
        end interface
        integer, intent(out) :: error
        integer :: minSteps, maxSteps, numSteps, n
        real(kind=DP) :: stepLength, left, right, a, b
        if (rightBound < leftBound) then
            error = 2
            return
        end if
        minSteps = 1
        maxSteps = minSteps * 256
        numSteps = minSteps
        varyStepSize: do
            stepLength = (rightBound - leftBound) / numSteps
            left = leftBound
            right = leftBound + stepLength
            varyFunctionArgs: do n=1, numSteps
                a = targetFunc(left)
                b = targetFunc(right)
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
    function BisectionRoot(targetFunc, leftBracket, rightBracket, error)
        real(kind=DP), intent(in) :: leftBracket, rightBracket
        integer, intent(out) :: error
        interface
            function targetFunc(x)
                use double
                real(kind=DP), intent(in) :: x
                real(kind=DP) :: targetFunc
            end function
        end interface
        real(kind=DP) :: BisectionRoot, a, b, low, high, mid
        a = targetFunc(leftBracket)
        b = targetFunc(rightBracket)
        if (.not. ((a >= 0 .and. b <= 0) .or. (a <= 0 .and. b >= 0))) then
            error = 1
            return
        end if
        if (a <= 0) then
            low = a
            high = b
        else
            low = b
            high = a
        end if
        mid = low + (high - low) / 2.0_DP
        do
            if ((mid == low) .or. (mid == high)) then
                error = 0
                BisectionRoot = mid
                return
            end if
            if (targetFunc(mid) <= 0) then
                low = mid
            else
                high = mid
            end if
            mid = low + (high - low) / 2.0_DP
        end do
    end function
end module rootfinder
