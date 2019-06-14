module assign

	implicit none
	public :: read_control_file, read_sequence_file, read_peak_files, &
		read_connection_file, initialize, evaluate_peak_rsd, evaluate_idv, &
		evaluate_group, delete_z_i, write_output, write_tables

	type info_spectra
	!
	! Spectra Information
	! ===================

	! peak numbers of the spectrums
	integer :: num_peak = 0
	! numbers of the kinds of chemical shifts
	integer :: num_freq = 0
	! chemical shifts
	real, allocatable :: ch_shifts(:,:)
	! resolution (half width)
	real, allocatable :: e_width(:,:)
	! degeneracy number
	integer, allocatable ::	degeneracy(:)
	! possible residue corresponding to each peak
	character*20, allocatable :: poss_rsd(:)
	! the used times of each peak
	character*100, allocatable :: poss_rsd_str(:)
	! number of time the peak was used.
	integer, allocatable ::	num_used(:)
	! the possible peaks for each residue
	integer, allocatable ::	peak_seq(:,:)
	! number of possible peaks for each residue
	integer, allocatable ::	num_poss_peak(:)
	! the possible residues for each peak
	integer, allocatable ::	prsd_peak(:,:)
	! number of possible residues for each peak
	integer, allocatable ::	num_poss_rsd(:)

	end type


	type idv
	!
	! Individual Information
	! ======================
	integer, allocatable:: 	rsd_pk(:, :)
	integer :: 	n_good = 0
	integer ::	n_bad = 0
	integer ::	n_edge = 0
	integer ::	n_used = 0

	end type

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

	end subroutine read_control_file


	subroutine read_sequence_file(file_sequence, N_seq, residue_seq)
		!
		! Sequence must be less then 500 residues.

		implicit none

		! Input Varables.
		character(len=128),intent(in) :: file_sequence
		! Output Varables.
		integer, intent(out) :: N_seq
		character(len=1), allocatable, intent(out) :: residue_seq(:)

		! Temporary varables.
		integer error, i, k, k1
		character*128 :: temp_c
		character*500	seq_temp
		logical	alive
		character*1, allocatable ::	rsd_temp(:)

		open(10, file = file_sequence, status = "old", iostat = error)
		if(error /= 0) then
			write(*,*) "Fail to open sequence file: ", file_sequence
			stop
		endif

		allocate (rsd_temp(500))
		N_seq = 0
		do while (.true.)
			read(10, *, iostat= error) seq_temp
			if (error /= 0) exit
			k = lnblnk(seq_temp)
			do k1 = 1, k
				N_seq = N_seq + 1
				rsd_temp(N_seq) = seq_temp(k1:k1)
			end do
		end do
		allocate (residue_seq(N_seq))

		do i = 1, N_seq
			residue_seq(i) = rsd_temp(i)
		end do

		deallocate (rsd_temp)

		close(10)
	end subroutine read_sequence_file


	subroutine read_peak_files(n_table, file_spectra_name, N_seq, info_spectrum)

		implicit none

		integer, intent(in) :: n_table, N_seq
		character(len=128), intent(in) :: file_spectra_name(n_table)
		type(info_spectra), intent(out) :: info_spectrum(n_table)

		integer :: error, i, k, k1, k2, p1
		integer :: lnblnk
		integer :: num_peak_temp, num_freq_temp
		character(len=100) :: prsd_temp, prsd_str
		integer	:: lock_n, bg_n

		! Read peak files
		do k = 1, n_table
			open(10, file = file_spectra_name(k), status = "old", iostat = error)
			if(error /= 0) then
				write(*,*) "Fail to open spectrum file ",k,": ", &
					file_spectra_name(k)
				stop
			endif

			read(10, *) num_peak_temp, num_freq_temp

			info_spectrum(k)%num_peak = num_peak_temp
			info_spectrum(k)%num_freq = num_freq_temp
			!
			allocate (info_spectrum(k)%ch_shifts(num_peak_temp,num_freq_temp))
			allocate (info_spectrum(k)%e_width(num_peak_temp,num_freq_temp))
			allocate (info_spectrum(k)%degeneracy(num_peak_temp))
			allocate (info_spectrum(k)%poss_rsd(num_peak_temp))
			allocate (info_spectrum(k)%poss_rsd_str(num_peak_temp))
			allocate (info_spectrum(k)%num_used(num_peak_temp))
			allocate (info_spectrum(k)%num_poss_peak(N_seq))
			allocate (info_spectrum(k)%num_poss_rsd(num_peak_temp))

			do i = 1, num_peak_temp
				read(10,*) (info_spectrum(k)%ch_shifts(i, k1), k1 = 1, num_freq_temp), &
					&	(info_spectrum(k)%e_width(i, k1), k1 = 1, num_freq_temp), &
					&	info_spectrum(k)%degeneracy(i), info_spectrum(k)%poss_rsd_str(i)
				info_spectrum(k)%num_used(i) = 0
				! store the possible residue type into info_spectrum(k)%poss_rsd(i)
				prsd_str = info_spectrum(k)%poss_rsd_str(i)
				prsd_temp = ''
				lock_n = 0
				p1 = 0
				do k2 = 1, lnblnk(prsd_str)
					if (prsd_str(k2:k2) .eq. '(') then
						lock_n = 1
						cycle
					elseif (prsd_str(k2:k2) .eq. ')') then
						lock_n = 0
						cycle
					end if
					if (lock_n .eq. 0) then
						p1 = p1+1
						prsd_temp(p1:p1) = prsd_str(k2:k2)
					end if
				end do
				info_spectrum(k)%poss_rsd(i) = prsd_temp(1:p1)
			end do

			close(10)
		end do

	end subroutine read_peak_files


	subroutine read_connection_file(file_connection_tab, cnn_row, &
		connection_table)

		character(len=128), intent(in) :: file_connection_tab
		integer, intent(out) :: cnn_row
		integer, allocatable, intent(out)::	connection_table(:,:)

		integer :: error, i, k

		open(10, file = file_connection_tab, status = "old", iostat = error)
		if(error /= 0) then
			write(*,*) "Fail to open connection table file: ", file_connection_tab
			stop
		endif

		! read the nuber of rows.s
		read(10, *) cnn_row

		! read the connection table.
		allocate (connection_table(cnn_row, 6))
		do k = 1, cnn_row
			read(10, *) (connection_table(k,i), i = 1, 6)
		end do
		close (10)

	end subroutine read_connection_file


	subroutine initialize(group_size, n_table, N_seq, residue_seq, &
		info_spectrum, parents, offsprings, temp_p)

		implicit none

		integer, intent(in)	:: group_size, n_table, N_seq
		character(len=1), intent(in) :: residue_seq(N_seq)
		type(info_spectra), intent(inout) :: info_spectrum(n_table)
		type(idv), intent(out) :: parents(group_size)
		type(idv), intent(out) :: offsprings(group_size)
		type(idv), intent(out) :: temp_p(group_size)

		! Temp varables.
		! peak_seq(:,k) contains the peak possiblites for the kth residue
		integer, allocatable ::	peak_seq_temp(:,:)
		integer :: i, k, k1, k2, k3, k4, p1, p2
		integer :: num_peak_temp
		character(len=20) :: idx_str
		character(len=100) :: prsd_temp, prsd_str
		integer :: new_asg, old_asg
		integer :: idx_table
		real :: rand_1
		integer :: residue_peak(N_seq, n_table)

		! initialize
		do k = 1, group_size
			allocate (parents(k)%rsd_pk(N_seq, n_table))
			allocate (temp_p(k)%rsd_pk(N_seq, n_table))
			allocate (offsprings(k)%rsd_pk(N_seq, n_table))
			parents(k)%rsd_pk = 0
			temp_p(k)%rsd_pk = 0
			offsprings(k)%rsd_pk = 0
		end do

		!**************************************************************
		!	generate peak-residue corresponding matrix
		!	peak_seq(:,k) contains the possible peaks corresponding to the kth residue
		!
		do k = 1, n_table
			allocate (peak_seq_temp(info_spectrum(k)%num_peak, N_seq))
			peak_seq_temp = 0
			info_spectrum(k)%num_poss_peak = 0
			!
			do i = 1, N_seq
				k4 = 0
				do k1 = 1, info_spectrum(k)%num_peak
					prsd_temp = info_spectrum(k)%poss_rsd(k1)
					do k2 = 1, lnblnk(prsd_temp)
						if (residue_seq(i) .eq. prsd_temp(k2:k2)) then
							k4 = k4 + 1
							peak_seq_temp(k4,i) = k1
							exit
						endif
					end do
				end do
				info_spectrum(k)%num_poss_peak(i) = k4
			end do
			!
			p1 = maxval(info_spectrum(k)%num_poss_peak(:))
			allocate(info_spectrum(k)%peak_seq(p1, N_seq))
			info_spectrum(k)%peak_seq(1:p1, 1:N_seq) = peak_seq_temp(1:p1, 1:N_seq)
			deallocate (peak_seq_temp)
		end do
		!
		do k = 1, n_table
			info_spectrum(k)%num_poss_rsd(:) = 0
			do k1 = 1, N_seq
				do k2 = 1, info_spectrum(k)%num_poss_peak(k1)
					i = info_spectrum(k)%peak_seq(k2,k1)
					info_spectrum(k)%num_poss_rsd(i) = info_spectrum(k)%num_poss_rsd(i)+1
				end do
			end do
			num_peak_temp = info_spectrum(k)%num_peak
			i = maxval(info_spectrum(k)%num_poss_rsd(:))
			allocate(info_spectrum(k)%prsd_peak(i,num_peak_temp))
		end do
		do k = 1, n_table
			info_spectrum(k)%prsd_peak(:,:) = 0
			num_peak_temp = info_spectrum(k)%num_peak
			do k1 = 1, num_peak_temp
				p1 = 0
				do k2 = 1, N_seq
					if (any(info_spectrum(k)%peak_seq(:,k2) .eq. k1)) then
						p1 = p1+1
						info_spectrum(k)%prsd_peak(p1, k1) = k2
					end if
				end do
			end do
		end do
		!
		do k = 1, n_table
			do k1 = 1, info_spectrum(k)%num_peak
				i = 0
				prsd_temp = info_spectrum(k)%poss_rsd(k1)
				do k2 = 1, lnblnk(prsd_temp)
					k4 = ichar(prsd_temp(k2:k2))
					if((k4.gt.47).and.(k4.lt.58))	then
						i = i + 1
						idx_str(i:i) = prsd_temp(k2:k2)
					endif
				end do
				if (i .gt. 0) then
					read(idx_str(1:i), '(i4)') k3
					info_spectrum(k)%num_poss_peak(k3) = -1
					write(*,*) 'assigned peaks:'
					write(*,*) 'table: ', k, ' residue: ', k3, ' peaks: ', k1
					do k2 = 1, group_size
						parents(k2)%rsd_pk(k3, k) = k1
					end do
				endif
			end do
		end do
		! Initialize the group members
		! if (file_initial_group .eq. 'NULL' &
		! 	.or. file_initial_group .eq. 'Null' &
		! 	.or. file_initial_group .eq. 'null') then
		! 	write(*,*)	'Initialize the group RANDOMLY'
		do k1 =1, group_size
			residue_peak = 0
			do idx_table = 1, n_table
				do k2 = 1, N_seq
					if (info_spectrum(idx_table)%num_poss_peak(k2) .eq. -1) then
						new_asg = parents(k1)%rsd_pk(k2, idx_table)
						old_asg = 0
						info_spectrum(idx_table)%num_used(new_asg) = &
						&info_spectrum(idx_table)%num_used(new_asg)+1
						parents(k1)%n_used = parents(k1)%n_used+1
						!
						residue_peak(k2, idx_table) = new_asg
						else
							k = info_spectrum(idx_table)%num_poss_peak(k2)
							call random_number(rand_1)
							k3 = int(rand_1*(k+1))+1
							if (k3 .le. k) then
							new_asg = info_spectrum(idx_table)%peak_seq(k3, k2)
							p1 = 0
							do p2 = 1, N_seq
								if (residue_peak(p2, idx_table) .eq. new_asg) p1 = p1+1
							end do
							if (p1 .eq. info_spectrum(idx_table)%degeneracy(new_asg)) cycle
							old_asg = 0
							parents(k1)%n_used = parents(k1)%n_used+1
							!
							residue_peak(k2, idx_table) = new_asg
						end if
					end if
				end do
			end do
			parents(k1)%rsd_pk = residue_peak
		end do
		! else
		! 	write(*,*)	'Initialize the group using the input file: ', file_initial_group
		! 	open(10, file = file_initial_group, status = "old", iostat = error)
		! 	if(error /= 0) then
		! 		write(*,*) 'Fail to open the initializing file: ', file_initial_group
		! 		stop
		! 	endif
		! 	!
		! 	read(10,*)	temp_c
		! 	if (trim(adjustL(temp_c)) .ne. 'Output_data')	then
		! 		write(*,*)	'The title of the input initial data file does not begin with "Output_data", please check.'
		! 		stop
		! 	end if
		! 	!
		! 	do k1 = 1, group_size
		! 		do idx_table = 1, n_table
		! 			read(10, *) (parents(k1)%rsd_pk(k2, idx_table), k2 = 1, N_seq)
		! 		end do
		! 	end do
		! 	close(10)
		! end if

	end subroutine initialize

	subroutine evaluate_peak_rsd (residue_peak, info_spectrum, &
		connection_table, n_table, N_seq, cnn_row, table_type, ind_rsd, asg, &
		ngood, nbad, nedge)
	!
	! Evaluate peak to residue matach subroutine
	! ------------------------------------------
	!
	! Evaluates the residue-to-peak connections.
	!
	implicit none

	integer	n_table, N_seq, cnn_row, ind_rsd, asg, table_type
	integer ngood, nbad, nedge
	integer ::	residue_peak(N_seq, n_table)
	type(info_spectra) :: info_spectrum(n_table)
	integer ::	connection_table(cnn_row, 6)

	integer k1, k2, peak_1, peak_2
	real	delta_1, delta_2, e_1, e_2, ch_1, ch_2
	integer	tab1, tab2, col1, col2, shift_ind, shift_ind_1, shift_ind_2
	integer ind_der, ind_tab
	integer	::	ngood_tab(n_table, 3), nbad_tab(n_table, 3)
	integer	::	nedge_tab(n_table, 3)

	ngood = 0
	nbad = 0
	nedge = 0
	do k1 = 1, n_table
		do k2 = 1, 3
			ngood_tab(k1, k2) = 0
			nbad_tab(k1, k2) = 0
			nedge_tab(k1, k2) = 0
		end do
	end do

	! evaluate the assignment 'asg' in the 'ind_rsd' residue
	do k1 = 1, cnn_row
		tab1 = connection_table(k1, 1)
		tab2 = connection_table(k1, 2)
		col1 = connection_table(k1, 3)
		col2 = connection_table(k1, 4)
		shift_ind_1 = connection_table(k1, 5)
		shift_ind_2 = connection_table(k1, 6)
		shift_ind = shift_ind_1-shift_ind_2
		if (tab1 .ne. table_type) then
			if (tab2 .ne. table_type)	cycle
		endif
		if (tab1 .eq. table_type) then
			peak_1 = asg
			ind_tab = tab2
			ind_der = shift_ind+2
			if (((ind_rsd+shift_ind) .gt. 0) &
			    .and. ((ind_rsd+shift_ind) .le. N_seq))then
				peak_2 = residue_peak(ind_rsd+shift_ind, tab2)
			else
				if ((shift_ind_1+ind_rsd) .gt. N_seq &
				   .or. (shift_ind_1+ind_rsd) .le. 0) then
					if ((info_spectrum(tab1)%ch_shifts(peak_1, col1) .gt. 0) &
					.and.&
					(info_spectrum(tab1)%ch_shifts(peak_1, col1) .lt. 1000)) then
						nbad_tab(ind_tab, ind_der) = nbad_tab(ind_tab, ind_der) + 1
					end if
				end if
				cycle
			endif
		else
			peak_2 = asg
			ind_tab = tab1
			ind_der = -shift_ind + 2
			if (((ind_rsd-shift_ind) .gt. 0) &
			   .and. ((ind_rsd-shift_ind) .le. N_seq)) then
				peak_1 = residue_peak(ind_rsd-shift_ind, tab1)
			else
				if ((ind_rsd+shift_ind_2) .gt. N_seq &
				    .or. (ind_rsd+shift_ind_2) .le. 0) then
					if ((info_spectrum(tab2)%ch_shifts(peak_2, col2) .gt. 0) &
					.and. &
					& (info_spectrum(tab2)%ch_shifts(peak_2, col2) .lt. 1000))	then
						nbad_tab(ind_tab, ind_der) = nbad_tab(ind_tab, ind_der) + 1
					end if
				end if
				cycle
			endif
		endif
		if (peak_1 .gt. 0) then
			e_1	 = info_spectrum(tab1)%e_width(peak_1, col1)
			ch_1 = info_spectrum(tab1)%ch_shifts(peak_1, col1)
			if (ch_1 .gt. 1000)		cycle
		endif
		if (peak_2 .gt. 0) then
			e_2	 = info_spectrum(tab2)%e_width(peak_2, col2)
			ch_2 = info_spectrum(tab2)%ch_shifts(peak_2, col2)
			if (ch_2 .gt. 1000)		cycle
		endif
		if (peak_1*peak_2 .eq. 0) then
			if ((peak_1+peak_2) .ne. 0)  then
				! edge assignment
				nedge_tab(ind_tab, ind_der) = nedge_tab(ind_tab, ind_der) + 1
			endif
		else
			delta_1 = e_1*e_1 + e_2*e_2
			delta_2 = ch_1 - ch_2
			delta_2 = delta_2*delta_2
			if (delta_1 .gt. delta_2) then
				! good assignment
				ngood_tab(ind_tab, ind_der) = ngood_tab(ind_tab, ind_der) + 1
			else
				! bad assignment
				nbad_tab(ind_tab, ind_der) = nbad_tab(ind_tab, ind_der) + 1
			endif
		endif
	end do
	do k1 = 1, n_table
		do k2 = 1, 3
			if (ngood_tab(k1, k2) .gt. 0) ngood = ngood + 1
			if (nbad_tab(k1, k2) .gt. 0) nbad = nbad + 1
			if (nedge_tab(k1, k2) .gt. 0) nedge = nedge + 1
		end do
	end do

	end subroutine evaluate_peak_rsd


	subroutine evaluate_idv(residue_peak, info_spectrum, connection_table, &
		n_table, N_seq, cnn_row, ngood, nbad, nedge, nused)
	!
	! Evaluate individuals subroutine
	! ===============================
	!

	implicit none

	integer, intent(in) :: n_table, N_seq, cnn_row
	integer, intent(in) ::	residue_peak(N_seq, n_table)
	type(info_spectra), intent(in)::	info_spectrum(n_table)
	integer, intent(in) ::	connection_table(cnn_row, 6)
	integer, intent(out) ::	ngood, nbad, nedge, nused

	integer asg, k1, k2
	integer :: rsd_peak_temp(N_seq, n_table)
	integer	ngood_new, ngood_old, nbad_new, nbad_old, nedge_new, nedge_old

	ngood = 0
	nbad = 0
	nedge = 0
	nused = 0
	rsd_peak_temp = 0
	do k1 = 1, N_seq
		do k2 = 1, n_table
			asg = residue_peak(k1, k2)

			call evaluate_peak_rsd(rsd_peak_temp, info_spectrum,  &
			connection_table, n_table, N_seq, cnn_row, k2, k1, asg, ngood_new, &
	 		nbad_new, nedge_new)

			call evaluate_peak_rsd(rsd_peak_temp, info_spectrum, &
			connection_table, n_table, N_seq, cnn_row, k2, k1, 0, ngood_old, &
			nbad_old, nedge_old)

			ngood = ngood+ngood_new-ngood_old
			nbad = nbad+nbad_new-nbad_old
			nedge = nedge+nedge_new-nedge_old
			!
			if (asg .ne. 0) nused = nused+1
			rsd_peak_temp(k1,k2) = asg
		end do
	end do
	end subroutine evaluate_idv


	subroutine evaluate_group(group_size, n_table, N_seq, cnn_row, &
		info_spectrum, connection_table, parents)
		!
		! Evaluate individuals in group subroutine
		! ============================================
		!

		implicit none

		integer, intent(in) :: group_size, n_table, N_seq, cnn_row
		type(info_spectra), intent(in) :: info_spectrum(n_table)
		integer, intent(in) :: connection_table(cnn_row, 6)
		type(idv), intent(inout) :: parents(group_size)

		integer :: residue_peak(N_seq, n_table)
		integer	::	n_good_new, n_edge_new, n_bad_new, n_used_new
		integer :: k1

		!	evaluate the individuals in the group
		do k1 = 1, group_size
			residue_peak = parents(k1)%rsd_pk

			call evaluate_idv(residue_peak, info_spectrum, connection_table, &
				n_table, N_seq, cnn_row, n_good_new, n_bad_new, n_edge_new, &
				n_used_new)

			parents(k1)%n_good = n_good_new
			parents(k1)%n_bad = n_bad_new
			parents(k1)%n_edge = n_edge_new
			parents(k1)%n_used = n_used_new
		end do

	end subroutine evaluate_group

	subroutine delete_z_i(residue_peak, parents, offsprings, ind_new_group, &
		N_seq, n_table, size_group, del_note)
	!
	! Delete Individuals with zero
	! ============================
	!
	! Delete the all-zero individuals and identical ones in the group.
	!
	! Input
	! -----
	! N_seq
	! n_table
	! size_group
	! residue_peak(N_seq, n_table)
	! parents(size_group)
	! offsprings(size_group)
	!
	! Output
	! ------
	! parents(size_group), offsprings(size_group)
	! del_note
	!
	implicit none

	integer, intent(in) ::	N_seq, n_table, size_group, ind_new_group
	integer, intent(in) ::	residue_peak(N_seq, n_table)
	type(idv), intent(in) :: parents(size_group), offsprings(size_group)
	integer, intent(out) ::	del_note

	integer ::	diff_residue_peak(N_seq)
	integer	k, i

	! delete the zero ones & identical one in the group
	if (all(residue_peak .eq. 0)) then
		del_note = 1
		return
	end if
	do k = 1, size_group
		del_note = 1
		do i = 1, n_table
			diff_residue_peak = parents(k)%rsd_pk(:, i) - residue_peak(:, i)
			if (all(diff_residue_peak .eq. 0)) cycle
			del_note = 0
			exit
		end do
		if (del_note .eq. 0) cycle
		exit
	end do
	if (del_note .eq. 1) then
		return
	end if

	if (ind_new_group > 1) then
		do k = 1, ind_new_group
			del_note = 1
			do i = 1, n_table
				diff_residue_peak = offsprings(k)%rsd_pk(:, i)-residue_peak(:, i)
				if (all(diff_residue_peak .eq. 0)) then
					cycle
				end if
				del_note = 0
				exit
			end do
			if (del_note .eq. 0) cycle
			exit
		end do
		if (del_note .eq. 1) then
			return
		end if
	end if
	end subroutine delete_z_i


subroutine write_output(file_out_data, group_size, n_table, parents)

	integer, intent(in) :: group_size, n_table
	character(len=128), intent(in) :: file_out_data
	type(idv), intent(in) :: parents(n_table)

	integer :: i, k1

	open(10, file = file_out_data)
	write(10,*)	 'Output_data'
	do i = 1, group_size
		do k1 = 1, n_table
			write(10,*) parents(i)%rsd_pk(:, k1)
		end do
	end do
	close(10)
end subroutine


subroutine write_tables(group_size, N_seq, n_table, file_out_tables, parents, &
	Pareto_order, info_spectrum, residue_seq)

	integer, intent(in) :: group_size, N_seq, n_table
	character(len=128), intent(in) :: file_out_tables
	type(idv), intent(in) :: parents(group_size)
	integer, intent(in) :: Pareto_order(group_size)
	type(info_spectra), intent(in) :: info_spectrum(n_table)
	character(len=1), intent(in):: residue_seq(N_seq)

	real:: freq(n_table, info_spectrum(1)%num_freq)
	integer :: i, k1, k2, k3, k4, k5, p1, p2, p3
	integer	::	idx_table
	integer	lock_n, bg_n, k_peak
	character(len=100) :: prsd_str, prsd_temp
	character(len=20) :: idx_str
	character(len=80) :: fmt_s

	i = 0
	do k1 = 1, group_size
		i = i+1
		p1 = i/100
		p2 = (i-p1*100)/10
		p3 = i-p1*100-p2*10
		open(11, file = trim(adjustL(file_out_tables))//achar(48+p1)//achar(48+p2)//achar(48+p3)//'.txt')
		write(11,"('Ng, Nb, Ne, Nu and Pareto order: 'I4,2x,I4,2x,I4,2x,I4,2x,I3)") &
			& parents(i)%n_good, parents(i)%n_bad, parents(i)%n_edge, parents(i)%n_used, Pareto_order(k1)
		do k2 = 1, N_seq
			prsd_temp = ''
			k4 = 0
			do idx_table = 1, n_table
				if (parents(k1)%rsd_pk(k2, idx_table) .ne. 0) then
					do k3 = 1, info_spectrum(idx_table)%num_freq
						freq(idx_table,k3) = info_spectrum(idx_table)%ch_shifts(parents(k1)%rsd_pk(k2, idx_table),k3)
						if (freq(idx_table, k3) .ge. 1000)	freq(idx_table, k3) = 0.0
					end do
				else
					freq(idx_table,:) = 0
				end if
				k_peak = parents(k1)%rsd_pk(k2, idx_table)
				prsd_str = info_spectrum(idx_table)%poss_rsd_str(k_peak)
				lock_n = 0
				bg_n = 0
				do k3 = 1, lnblnk(prsd_str)
					if (prsd_str(k3:k3) .eq. '(') then
						lock_n = 1
						cycle;
					elseif (prsd_str(k3:k3) .eq. ')') then
						lock_n = 0
						bg_n = 0
						cycle;
					end if
					if (prsd_str(k3:k3) .eq. residue_seq(k2)) then
						if (lock_n == 0) then
							bg_n = 1
							cycle;
						end if
					else
						if (lock_n == 0) then
							bg_n = 0
							cycle;
						end if
					end if
					if (lock_n .eq. 1 .and. bg_n .eq. 1) then
						k4 = k4+1
						prsd_temp(k4:k4) = prsd_str(k3:k3)
					end if
				end do
			end do
			prsd_str = prsd_temp
			do k3 = lnblnk(prsd_temp), 2, -1
				do k5 = 1, k3-1
					if (prsd_temp(k3:k3) .eq. prsd_temp(k5:k5)) then
						prsd_str = prsd_str(1:k3-1)//prsd_str(k3+1:k4)
						k4 = k4-1
						exit
					end if
				end do
			end do
			write(11,*) k2, residue_seq(k2), (parents(k1)%rsd_pk(k2, idx_table), freq(idx_table,:), idx_table = 1, n_table), &
				&prsd_str(1:k4)

		end do
		close(11)
	end do

end subroutine write_tables

end module
