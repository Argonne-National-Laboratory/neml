      program main
            use iso_c_binding
            implicit none
            include "neml_interface.f"

            type(c_ptr) :: model
            integer :: ier, nstore, i
            double precision, allocatable, dimension(:) :: h_n, h_np1
            double precision :: s_n(6), s_np1(6), e_n(6), e_np1(6)
            double precision :: ee_np1(6)
            double precision :: A_np1(6,6)
            double precision :: temp_np1, temp_n, time_np1, time_n,
     &            u_np1, u_n, p_np1, p_n

            character(len=64) :: fname_arg, mname_arg, e_arg, time_arg
            character(len=64) :: nsteps_arg, temp_arg
            character(len=65,kind=c_char) :: fname, mname
            
            double precision :: temp, time, e
            integer :: nsteps
            
            ! Setup arguments
            if (COMMAND_ARGUMENT_COUNT() .ne. 6) then
                  write(*,*) "Expected 6 arguments"
                  write(*,*) "   XML file, model name, max strain"
                  write(*,*) "   max time, nsteps, temperature"
                  stop
            endif

            call getarg(1, fname_arg)
            call getarg(2, mname_arg)
            call getarg(3, e_arg)
            call getarg(4, time_arg)
            call getarg(5, nsteps_arg)
            call getarg(6, temp_arg)

            fname = trim(fname_arg)//C_NULL_CHAR
            mname = trim(mname_arg)//C_NULL_CHAR

            read(e_arg,*) e
            read(time_arg,*) time
            read(nsteps_arg,*) nsteps
            read(temp_arg,*) temp

            model = create_nemlmodel(fname, mname, ier)
            if (ier .ne. 0) then
                  write(*,*) "Loading the model failed"
                  stop
            end if

            ! Allocate and step
            nstore = nstore_nemlmodel(model)
            allocate(h_np1(nstore))
            allocate(h_n(nstore))

            call init_store_nemlmodel(model, h_n, ier)
            if (ier .ne. 0) then
                  write(*,*) "Problem setting up history"
                  deallocate(h_np1)
                  deallocate(h_n)
                  stop
            end if

            e_n = 0.0
            s_n = 0.0

            time_n = 0.0
            temp_np1 = temp
            temp_n = temp

            u_n = 0.0
            p_n = 0.0

            ! Loop and update
            do i=1,nsteps
                  time_np1 = i/DBLE(nsteps) * time
                  e_np1 = 0.0
                  e_np1(1) = i/DBLE(nsteps) * e

                  call update_sd_nemlmodel(model, e_np1, e_n, temp_np1,
     &                  temp_n, time_np1, time_n, s_np1, s_n, h_np1, 
     &                  h_n, A_np1, u_np1, u_n, p_np1, p_n, ier)
                  if (ier .ne. 0) then
                        write(*,*) "Error in model update"
                        deallocate(h_np1)
                        deallocate(h_n)
                        stop
                  end if

                  call elastic_strains_nemlmodel(model, s_np1, temp_np1,
     &                  ee_np1, ier)

                  ! write(*,*) s_np1(1)

                  s_n = s_np1
                  e_n = e_np1
                  h_n = h_np1
                  time_n = time_np1
                  u_n = u_np1
                  p_n = p_np1
            end do

            ! Deallocate
            deallocate(h_np1)
            deallocate(h_n)

            call destroy_nemlmodel(model, ier)
            if (ier .ne. 0) then
                  write(*,*) "Destroying the model failed"
                  stop
            end if

      end program
