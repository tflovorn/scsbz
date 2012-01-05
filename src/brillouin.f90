module brillouin
    use double
    use environment
contains
    function BrilSum(env, sumFunc)
        ! uses Kahan summation algorithm to improve accuracy
        implicit none
        type(Environ), intent(in) :: env
        interface
            function sumFunc(env, k)
                use double
                use environment
                implicit none
                type(Environ), intent(in) :: env
                real(kind=DP), intent(in) :: k(1:3)
                real(kind=DP) :: sumFunc
            end function sumFunc
        end interface
        real(kind=DP) :: BrilSum, start, finish, step, k(1:3), y, t, c
        integer :: nx, ny, nz
        start = -ACos(-1.0_DP) ! start = -pi
        finish = -start
        step = (finish - start) / env%zoneLength
        BrilSum = 0.0_DP
        c = 0.0_DP
        k = start
        xLoop: do nx=1, env%zoneLength
            yLoop: do ny=1, env%zoneLength
                zLoop: do nz=1, env%zoneLength
                    y = sumFunc(env, k) - c
                    t = BrilSum + y
                    c = (t - BrilSum) - y
                    BrilSum = t
                    k(3) = k(3) + step
                end do zLoop
                k(2) = k(2) + step
                k(3) = start
            end do yLoop
            k(1) = k(1) + step
            k(2) = start
        end do xLoop
    end function BrilSum
end module brillouin       
