module driver
    use double
    use environment
    use sbzequations
    implicit none
    integer, parameter :: fileNameMaxLen = 256
contains
    ! Read data in specified input file, calculate self-consistent solution,
    ! and write solution to output file. Returns 0 if no errors encountered.
    function sbzRun(inputFileName, outputFileName)
        integer, parameter :: inputUnit = 23, outputUnit = 24, recLen = 250
        real(kind=DP), parameter :: tolerance = 1e-6
        character(len=fileNameMaxLen), intent(in) :: inputFileName,&
                                                     outputFileName
        integer :: sbzRun, error, zoneLength, runCount, i
        type(Environ) :: env
        type(SelfConsistentSystem) :: system
        real(kind=DP), dimension(:), allocatable :: toleranceValues
        real(kind=DP) :: t, tc, beta, x, J
        allocate(toleranceValues(1:sbzNumEquations))
        toleranceValues(1:sbzNumEquations) = tolerance
        open(inputUnit, file=inputFileName, recl=recLen)
        open(outputUnit, file=outputFileName, recl=recLen)
        read (inputUnit, *) runCount
        write (outputUnit, *) runCount
        sbzRun = 0
        do i=1,runCount
            read (inputUnit, *) zoneLength, t, tc, beta, x, J
            ! Start with a system with no c-direction hopping to get reasonable
            ! initial values for system with tc > 0.
            env = NewEnv(zoneLength, t, 0.0_DP, beta, x, J)
            system = sbzSystemNoC(env, toleranceValues)
            error = scSystemSolve(system, env)
            if (error > 0) then
                write (outputUnit, *) "!tc=0 run failed; error = ", error
                sbzRun = error
                cycle
            end if
            ! Done with tc = 0 approximation. Move on to real system.
            env%tc = tc
            system = sbzSystem(env, toleranceValues)
            error = scSystemSolve(system, env)
            if (error > 0) then
                write (outputUnit, *) "!run failed; error = ", error
                sbzRun = error
                cycle
            end if
            write (outputUnit, *) env%muF, env%muB, env%D, env%B, env%A
            write (outputUnit, *) env%Dc, env%Bc, env%Ac
            flush (outputUnit)
        end do
        close (inputUnit)
        close (outputUnit)
    end function
end module
