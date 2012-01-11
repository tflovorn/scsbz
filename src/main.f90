program main
    use driver
    implicit none
    call runMain()
contains
    subroutine runMain()
        integer :: length, stat, error
        character(len=fileNameMaxLen) :: inputFileName, outputFileName
        ! TODO: should probably check number of arguments here
        call get_command_argument(1, inputFileName, length, stat)
        if (stat /= 0) then
            print *, "Failed to read input file from command line"
            return
        end if
        call get_command_argument(2, outputFileName, length, stat)
        if (stat /= 0) then
            print *, "Failed to read output file from command line"
            return
        end if
        error = sbzRun(inputFileName, outputFileName)
        if (error /= 0) then
            print *, "sbzRun quit with errors"
        end if
    end subroutine
end program
