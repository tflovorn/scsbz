module environment
    use double
    implicit none
    type :: Environ
        ! static parameters
        integer :: zoneLength   ! number of points on a side in the Brill. zone
        real(kind=DP) :: t,&    ! bare in-plane hopping energy
                         tc,&   ! bare c-direction hopping energy
                         beta,& ! inverse temperature
                         x,&    ! hole doping
                         J      ! in-plane singlet energy
        ! self-consistent parameters (expectation values)
        real(kind=DP) :: D,&    ! in-plane holon hopping
                         Dc,&   ! c-direction holon hopping
                         B,&    ! in-plane spinon hopping
                         Bc,&   ! c-direction spinon hopping
                         A,&    ! in-plane singlet formation
                         Ac,&   ! c-direction singlet formation
                         muf,&  ! spinon chemical potential
                         mub    ! holon chemical potential   
    end type Environ
contains
    ! dependent parameters
    ! c-direction singlet energy
    function Jc(env)
        type(Environ), intent(in) :: env
        real(kind=DP) :: Jc
        Jc = (env%tc / env%t) ** 2.0_DP
    end function Jc

    ! effective in-plane spinon hopping energy
    function Eh(env)
        type(Environ), intent(in) :: env
        real(kind=DP) :: Eh
        Eh = env%t * env%D + 0.5_DP * env%J * env%B
    end function Eh

    ! effective c-direction spinon hopping energy
    function Ehc(env)
        type(Environ), intent(in) :: env
        real(kind=DP) :: Ehc
        Ehc = env%tc * env%Dc + 0.5_DP * Jc(env) * env%Bc
    end function Ehc

    ! minimum of the spinon spectrum
    function epsilonMin(env)
        type(Environ), intent(in) :: env
        real(kind=DP) :: epsilonMin
        epsilonMin = -2.0_DP * (2.0_DP * Eh(env) + Ehc(env))
    end function epsilonMin

    ! Momentum-dependent functions.
    ! Spinon spectrum.
    function epsilonBar(env, k)
        type(Environ), intent(in) :: env
        real(kind=DP), dimension(1:3) :: k
        real(kind=DP) :: epsilonBar
        epsilonBar = -2.0_DP * (Eh(env) * (cos(k(1)) + cos(k(2))) +&
                                Ehc(env) * cos(k(3)))
    end function epsilonBar

    ! Spinon (diagonal-part) hopping energy, including chemical potential.
    ! Minimum is -mu_f.
    function xi(env, k)
        type(Environ), intent(in) :: env
        real(kind=DP), dimension(1:3) :: k
        real(kind=DP) :: xi
        xi = epsilonBar(env, k) - epsilonMin(env) - env%muf
    end function xi

    ! Gap function (spinon k-up, minus-k-down coupled part).
    function Delta(env, k)
        type(Environ), intent(in) :: env
        real(kind=DP), dimension(1:3) :: k
        real(kind=DP) :: Delta
        Delta = -2.0_DP * (env%J * env%A * (cos(k(1)) - cos(k(2))) +&
                           Jc(env) * env%Ac * cos(k(3)))
    end function Delta

    ! Energy for spinon-based pseudoparticles obtained after doing Bogoliubvov
    ! transformation.
    function bogoEnergy(env, k)
        type(Environ), intent(in) :: env
        real(kind=DP), dimension(1:3) :: k
        real(kind=DP) :: bogoEnergy
        bogoEnergy = sqrt(Delta(env, k) ** 2.0_DP + xi(env, k) ** 2.0_DP)
    end function bogoEnergy

    ! Constructor for Environ.
    function NewEnv(zoneLength, t, tc, beta, x, J)
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
        NewEnv%muf = 0.1    ! should be positive
        NewEnv%mub = -0.1   ! should be negative
    end function NewEnv
end module environment
