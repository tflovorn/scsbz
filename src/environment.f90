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
                         muF,&  ! spinon chemical potential
                         muB    ! holon chemical potential   
    end type Environ
contains
    ! --- Dependent parameters. ---
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

    ! Minimum of the spinon spectrum.
    ! Should work whether Eh and Ehc are positive or negative.
    function epsilonMinF(env)
        type(Environ), intent(in) :: env
        real(kind=DP) :: epsilonMinF, abPart, cPart
        abPart = abs(2.0_DP * Eh(env))
        cPart = abs(Ehc(env))
        epsilonMinF = -2.0_DP * (abPart + cPart)
    end function epsilonMinF

    ! Minimum of the holon spectrum.
    ! Should work whether tB and tcBc are positive or negative.
    function epsilonMinB(env)
        type(Environ), intent(in) :: env
        real(kind=DP) :: epsilonMinB, abPart, cPart
        abPart = abs(2.0_DP * env%t * env%B)
        cPart = abs(env%tc * env%Bc)
        epsilonMinB = -4.0_DP * (abPart + cPart)
    end function epsilonMinB

    ! --- Momentum-dependent functions. ---
    ! Spinon spectrum.
    function epsilonBarF(env, k)
        type(Environ), intent(in) :: env
        real(kind=DP), dimension(1:3) :: k
        real(kind=DP) :: epsilonBarF
        epsilonBarF = -2.0_DP * (Eh(env) * (cos(k(1)) + cos(k(2))) +&
                                 Ehc(env) * cos(k(3)))
    end function epsilonBarF

    ! Spinon (diagonal-part) hopping energy, including chemical potential.
    ! Minimum is -muF.
    function xiF(env, k)
        type(Environ), intent(in) :: env
        real(kind=DP), dimension(1:3) :: k
        real(kind=DP) :: xiF
        xiF = epsilonBarF(env, k) - epsilonMinF(env) - env%muF
    end function xiF

    ! Holon spectrum.
    function epsilonBarB(env, k)
        type(Environ), intent(in) :: env
        real(kind=DP), dimension(1:3) :: k
        real(kind=DP) :: epsilonBarB
        epsilonBarB = -4.0_DP * (env%t * env%B * (cos(k(1)) + cos(k(2))) +&
                                 env%tc * env%Bc * cos(k(3)))
    end function epsilonBarB

    ! Holon hopping energy, including chemical potential.
    ! Minimum is -muB.
    function xiB(env, k)
        type(Environ), intent(in) :: env
        real(kind=DP), dimension(1:3) :: k
        real(kind=DP) :: xiB
        xiB = epsilonBarB(env, k) - epsilonMinB(env) - env%muB
    end function xiB

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
        bogoEnergy = sqrt(Delta(env, k) ** 2.0_DP + xiF(env, k) ** 2.0_DP)
    end function bogoEnergy

    ! --- Constructor for Environ. ---
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
        NewEnv%D = 0.1_DP
        NewEnv%Dc = 0.1_DP
        NewEnv%B = 0.1_DP
        NewEnv%Bc = 0.1_DP
        NewEnv%A = 0.1_DP
        NewEnv%Ac = 0.1_DP
        NewEnv%muF = 0.1_DP    ! should be positive
        NewEnv%muB = -0.1_DP   ! should be negative
    end function NewEnv
end module environment
