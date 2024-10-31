      subroutine convertv(input, indexes, multipliers, res)
            implicit none

            double precision, dimension(6), intent(in) :: input
            integer, dimension(6), intent(in) :: indexes
            double precision, dimension(6), intent(in) :: multipliers

            double precision, dimension(6), intent(out) :: res

            integer :: i

            do i=1,6
                  res(i) = input(indexes(i)) * multipliers(i)
            end do

            return
      end

      subroutine bconvertv(input, indexes, multipliers, res)
            implicit none

            double precision, dimension(6), intent(in) :: input
            integer, dimension(6), intent(in) :: indexes
            double precision, dimension(6), intent(in) :: multipliers

            double precision, dimension(6), intent(out) :: res

            integer :: i

            do i=1,6
                  res(indexes(i)) = input(i) / multipliers(indexes(i))
            end do

            return
      end

      subroutine bconvertt(input, indexes, multi, multj, res)
            implicit none

            double precision, dimension(6,6), intent(in) :: input
            integer, dimension(6), intent(in) :: indexes
            double precision, dimension(6), intent(in) :: multi, multj

            double precision, dimension(6,6), intent(out) :: res

            integer :: i, j

            do i=1,6
            do j=1,6
             res(indexes(i),indexes(j)) = input(i,j) / (
     1            multi(indexes(i)) * multj(indexes(j)))
            end do
            end do

            return
      end

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      use, intrinsic :: iso_c_binding
      INCLUDE 'ABA_PARAM.INC'
      include 'neml_interface.f'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
c
c           Hard-coded model names
c     
      character(len=64) :: fname_hc, mname_hc
      parameter(fname_hc='/home/messner/Documents/Projects/appendix-z/
     &abaqus/neml.xml')
      parameter(mname_hc='abaqus')
c
c           Used for NEML call
c
      character(len=65,kind=c_char) :: fname, mname
      type(c_ptr) :: model
      integer :: ier
      double precision, dimension(6) :: e_np1, e_n, s_np1, s_n
      integer, dimension(6) :: imap
      double precision, dimension(6) :: smult, emult
      double precision, dimension(6,6) :: A_np1
      double precision, dimension(NSTATV) :: h_np1, h_n
      double precision :: temp_np1, temp_n, time_np1, time_n,
     1 u_np1, u_n, p_np1, p_n
c
c           Setup the names as C strings
c
      fname = trim(fname_hc)//C_NULL_CHAR
      mname = trim(mname_hc)//C_NULL_CHAR
c
c           Setup the maps
c           These go from NEML -> ABAQUS
c
c           Indexes     
      imap(1) = 1
      imap(2) = 2
      imap(3) = 3
      imap(4) = 6
      imap(5) = 5
      imap(6) = 4
c
c           Stress
c
      smult(1) = 1.0
      smult(2) = 1.0
      smult(3) = 1.0
      smult(4) = sqrt(2.0)
      smult(5) = sqrt(2.0)
      smult(6) = sqrt(2.0)
c
c           Strain
c
      emult(1) = 1.0
      emult(2) = 1.0
      emult(3) = 1.0
      emult(4) = sqrt(2.0) / 2.0
      emult(5) = sqrt(2.0) / 2.0
      emult(6) = sqrt(2.0) / 2.0
c
c           Load the model
c
      model = create_nemlmodel(fname, mname, ier)
      if (ier .ne. 0) then
            write(*,*) "ERROR: Could not load NEML model!"
            stop
      end if
c
c           Map over quantities
c                 NEML        ABAQUS
c                 e_np1       STRAN + DSTRAN
c                 e_n         STRAN
c                 temp_np1    TEMP + DTEMP
c                 temp_n      TEMP
c                 time_np1    TIME(2) + DTIME
c                 time_n      TIME(2)
c                 s_np1       STRESS
c                 s_n         STRESS
c                 h_np1       STATEV
c                 h_n         STATEV
c                 A_np1       DDSDDE
c                 u_np1       n/a
c                 u_n         SSE + SPD
c                 p_np1       SPD
c                 p_n         SPD
c
c           Can get SSE as u_np1 - p_np1
c           Ignore SCD
c
      call convertv(STRAN + DSTRAN, imap, emult, e_np1)
      call convertv(STRAN, imap, emult, e_n)
      temp_np1 = TEMP + DTEMP
      temp_n = TEMP
      time_np1 = TIME(2) + DTIME
      time_n = TIME(2)
      call convertv(STRESS, imap, smult, s_n)
      h_n = STATEV
      u_n = SSE + SPD
      p_n = SPD
c
      call update_sd_nemlmodel(model, e_np1, e_n, temp_np1, temp_n,
     1 time_np1, time_n, s_np1, s_n, h_np1, h_n, A_np1, u_np1, u_n,
     2 p_np1, p_n, ier)
c
c           Only thing to do is reduce the timestep
c
      if (ier .ne. 0) then
            write(*,*) "WARNING: Model requested step reduction!"
            PNEWDT = 0.5
            return
      end if
c
c     Translate back
c
      call bconvertt(transpose(A_np1), imap, smult, emult, DDSDDE)
      call bconvertv(s_np1, imap, smult, STRESS)
      STATEV = h_np1
      SSE = u_np1 - p_np1
      SPD = p_np1
      SCD = 0.0
c
c           Destroy the model
c
      call destroy_nemlmodel(model, ier)
      if (ier .ne. 0) then
            write(*,*) "WARNING: Destroying the NEML model failed!"
            ! No error
      end if
c
      return

      END
