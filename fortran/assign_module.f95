module assign_module

	implicit none
	public :: read_control_file

	contains

	subroutine read_control_file(file_control, n_table, file_sequence, &
		file_initial_group, file_out_data, file_out_tables,  &
		file_connection_tab, file_spectra_name, group_size, pool_size, N_step, &
		N_try, num_free, mutate_rate, mutate_ad_rate, cross_rate, p_null, iseed)

		! Input Varables
		character*128, intent(in):: file_control
		! Output Varables
		integer, intent(out) :: n_table
		character*128, intent(out) :: file_sequence
		character*128, intent(out) :: file_initial_group
		character*128, intent(out) :: file_out_data
		character*128, intent(out) :: file_out_tables
		character*128, intent(out) :: file_connection_tab
		character*128, allocatable, intent(out) :: file_spectra_name(:)
		integer, intent(out) :: group_size, pool_size, N_step, N_try, num_free
		real, intent(out) :: mutate_rate, mutate_ad_rate, cross_rate, p_null
		integer, intent(out) :: iseed(12)
		! Temporary varables
		integer error, k
		character*128 :: temp_c
		logical	alive

		inquire(file = file_control, exist = alive)
		if (.not. alive)	then
			write(*,*)	'Cannot find file: ', file_control
			stop
		end if

		!	Read input file (including the files' locations and parameters)
		open(10, file = file_control, status = 'old', action = 'read', &
			iostat = error)

		! The first line of the file should begin with NSGA2_ASSIGN
		read(10, *)	temp_c
		temp_c = trim(adjustL(temp_c))
		if (temp_c(1:12) .ne. 'NSGA2_ASSIGN')	then
			write(*,*)	'The control file does not begin with NSGA2_ASSIGN.'
			stop
		end if

		! Brun a line
		read(10, *)

		! Read the name of the sequence file
		read(10, "(A128)", iostat=error) temp_c
		temp_c = trim(adjustL(temp_c))
		file_sequence = temp_c(1:index(temp_c,' ')-1)

		! Read the number of tables
		read(10, *, iostat=error)	n_table

		! Read the file spectra names
		allocate(file_spectra_name(n_table))
		do k = 1, n_table
			read(10, "(A128)", iostat=error) temp_c
			temp_c = trim(adjustL(temp_c))
			file_spectra_name(k) = temp_c(1:index(temp_c,' ')-1)
		end do

		! Read the connection table file name
		read(10, "(A128)", iostat=error) temp_c
		temp_c = trim(adjustL(temp_c))
		file_connection_tab = temp_c(1:index(temp_c,' ')-1)

		! Read the initial group file name
		read(10, "(A128)", iostat=error) temp_c
		temp_c = trim(adjustL(temp_c))
		file_initial_group = temp_c(1:index(temp_c,' ')-1)

		! Read the output file name
		read(10, "(A128)", iostat=error) temp_c
		temp_c = trim(adjustL(temp_c))
		file_out_data = temp_c(1:index(temp_c,' ')-1)

		! Read the output tab directory name
		read(10, "(A128)", iostat=error) temp_c
		temp_c = trim(adjustL(temp_c))
		file_out_tables = temp_c(1:index(temp_c,' ')-1)

		! Burn another line
		read(10, *)

		! Read the rest of the NSGA2 paramaters
		read(10, *, iostat=error) 	group_size
		read(10, *, iostat=error) 	pool_size
		read(10, *, iostat=error) 	N_step
		read(10, *, iostat=error) 	N_try
		read(10, *, iostat=error) 	num_free
		read(10, *, iostat=error) 	mutate_rate
		read(10, *, iostat=error) 	mutate_ad_rate
		read(10, *, iostat=error) 	cross_rate
		read(10, *, iostat=error) 	p_null
		read(10, *, iostat=error) 	iseed

		close(10)

	end subroutine

end module
