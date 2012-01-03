module environment
    use double
    implicit none
    type :: Environ
        ! static parameters
        integer :: zoneLength   ! number of points on a side in the Brill. zone
        real(kind=DP) :: t,&    ! in-plane hopping energy
                         tc,&   ! c-direction hopping energy
                         beta,& ! inverse temperature
                         x,&    ! hole doping
                         J      ! in-plane singlet energy
        ! self-consistent parameters (expectation values)
        real(kind=DP) :: D,&    ! in-plane holon hopping
                         Dc,&   ! c-direction holon hopping
                         B,&    ! in-plane spinon hopping
                         Bc,&   ! c-direction spinon hopping
                         A,&    ! in-plane singlet formation
                         Ac     ! c-direction singlet formation
    end type Environ
contains
        ! dependent parameters
        function Jc(env)
            implicit none
            type(Environ), intent(in) :: env
            real(kind=DP) :: Jc
            Jc = (env%tc / env%t) ** 2.0_DP
        end function Jc

        function Eh(env)
            implicit none
            type(Environ), intent(in) :: env
            real(kind=DP) :: Eh
            Eh = env%t * env%D + 0.5_DP * env%J * env%B
        end function Eh

        function Ehc(env)
            implicit none
            type(Environ), intent(in) :: env
            real(kind=DP) :: Ehc
            Ehc = env%tc * env%Dc + 0.5_DP * Jc(env) * env%Bc
        end function Ehc

        function epsilonMin(env)
            implicit none
            type(Environ), intent(in) :: env
            real(kind=DP) :: epsilonMin
            epsilonMin = -2.0_DP * (2.0_DP * Eh(env) + Ehc(env))
        end function epsilonMin

        ! constructor
        function NewEnv(zoneLength, t, tc, beta, x, J)
            implicit none
            integer, intent(in) :: zoneLength
            real(kind=DP), intent(in) :: t, tc, beta, x, J
            type(Environ) :: NewEnv
            ! copy parameters
            NewEnv%zoneLength = zoneLength
            NewEnv%t = t
            NewEnv%tc = tc
            NewEnv%beta = beta
            NewEnv%x = x
            NewEnv%J = J
            ! default values (TODO: change to reasonable guesses)
            NewEnv%D = 0.1
            NewEnv%Dc = 0.1
            NewEnv%B = 0.1
            NewEnv%Bc = 0.1
            NewEnv%A = 0.1
            NewEnv%Ac = 0.1
        end function NewEnv
end module environment
