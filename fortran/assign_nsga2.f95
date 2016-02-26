program assign_nsga2

	use qsort_c_module
	use assign_module
	use nsga2_module
	implicit none

	! parameters of typical example from input file.
	character*128 :: file_control
	integer	::	n_table
	character*128 :: file_sequence
	character*128 :: file_initial_group
	character*128 :: file_out_data
	character*128 :: file_out_tables
	character*128 :: file_connection_tab
	character*128, allocatable :: file_spectra_name(:)
	integer	:: group_size, pool_size, N_step, N_try, num_free
	real	:: mutate_rate, mutate_ad_rate, cross_rate, p_null
	integer :: iseed(12)

	! residue sequence length
	integer	::	N_seq
	character*1, allocatable::	residue_seq(:)

	! the possible peaks corresponding to each residue
	integer, allocatable::	residue_peak(:, :)
	! store the infomation of each table (spectra)
	type(info_spectra), allocatable	::	info_spectrum(:)
	! the residue sequence
	integer	::	cnn_row
	! the connections in different spectrum
	integer, allocatable::	connection_table(:,:)
	!*****************************************************
	type(idv), allocatable	::	parents(:), offsprings(:), temp_p(:)
	integer, allocatable	::	n_good_all(:), n_edge_all(:), n_bad_all(:), n_used_all(:)
	integer, allocatable	::	num_set(:), Pareto_set(:), Pareto_order(:), ind_pool(:)
	real, allocatable		::	dist_group(:)
	!
	integer	rank
	integer	::	idx_table, ind_new_group
	integer	:: 	err_sign
	integer	::	n_good_new, n_edge_new, n_bad_new, n_used_new
	real	::	lucky_ratio
	!*****************************************************
	! temperary parameters
	character*80	fmt_s
	integer	ind_ns, ind_nt
	type(idv)	::	parents2(2), child1, children2(2)
	integer i, error, k, k1, k2, k3, k4, k5, p1, p2, p3
	integer lnblnk
	integer num_peak_temp, num_freq_temp
	character*100	prsd_temp, prsd_str
	integer	lock_n, bg_n, k_peak
	character*20	::	idx_str
	integer 	new_asg, old_asg
	real		rand_1, rand_2, rand_3
	integer	::	timearray(3), new_timearray(3)
	integer	num_new_group, size_all
	integer	::	ind_k
	integer	del_note, note_rstr, num_obj
	integer, allocatable	::	ind_temp(:), ind_temp2(:)
	integer,allocatable	::	diff_residue_peak(:)
	real, allocatable	::	obj_4(:,:), dist_s(:), freq(:,:)
	integer, allocatable ::	peak_seq_temp(:,:)

	character*128 :: temp_c
	logical	alive

    write(*,*)  'Enter control file name:'
    read(*, "(A)") file_control

	call read_control_file(file_control, n_table, file_sequence, &
		file_initial_group, file_out_data, file_out_tables,  &
		file_connection_tab, file_spectra_name, group_size, pool_size, N_step, &
		N_try, num_free, mutate_rate, mutate_ad_rate, cross_rate, p_null, iseed)


	! check if inputs are right
	! write(*,*)	'Check the following input: '
	! write(*,*) 	'Sequence file name is: ', trim(file_sequence)
	! write(*,*)	'Number of tables (spectra) is: ', n_table
	! do k = 1, n_table
	! 	write(*,*) 	'Spectrum ',k,' file name is: ', trim(file_spectra_name(k))
	! end do
	! write(*,*) 	'Connection table file name is: ', trim(file_connection_tab)
	! write(*,*) 	'Initial data file name is: ', trim(file_initial_group)
	! write(*,*) 	'Output file (data) name is: ', trim(file_out_data)
	! write(*,*) 	'Output file (tables) name is: ', trim(file_out_tables)
	! write(*,*)	'Input Parameter: '
	! write(*,*)	'group size is: ', group_size
	! write(*,*)	'gene pool size is: ', pool_size
	! write(*,*)	'number of steps is: ', N_step
	! write(*,*)	'number of attempts is: ', N_try
	! write(*,*)	'number of free steps is: ', num_free
	! write(*,*)	'mutation rate is: ', mutate_rate
	! write(*,*)	'additional mutation rate is: ', mutate_ad_rate
	! write(*,*)	'crossover rate is: ', cross_rate
	! write(*,*)	'null probability is: ', p_null


	if (pool_size>group_size) then
		write(*,*) 'Wrong setting for pool_size!'
		stop
	end if

	allocate(info_spectrum(n_table))
	allocate(parents(group_size))
	allocate(offsprings(group_size))
	allocate(temp_p(group_size))
	allocate(n_good_all(2 * group_size))
	allocate(n_bad_all(2 * group_size))
	allocate(n_edge_all(2 * group_size))
	allocate(n_used_all(2 * group_size))
	allocate(num_set(2 * group_size))
	allocate(Pareto_set(2 * group_size))
	allocate(Pareto_order(group_size))
	allocate(ind_pool(pool_size))
	allocate(dist_group(group_size))
	allocate(ind_temp(2 * group_size))
	allocate(ind_temp2(2 * group_size))
	allocate(dist_s(2 * group_size))
	allocate(obj_4(2 * group_size, 4))


	! Read the infomation from files
	call read_sequence_file(file_sequence, N_seq, residue_seq)

	! write(*,*)'Input sequence has ', N_seq, ' residues: ', &
	!	(residue_seq(k), k=1,N_seq)

	allocate (residue_peak(N_seq, n_table))
	allocate (diff_residue_peak(N_seq))

	call read_peak_files(n_table, file_spectra_name, N_seq, info_spectrum)

	call read_connection_file(file_connection_tab, cnn_row, &
		connection_table)

	! Main Program
	call random_seed(put = abs(iseed))

	call initialize(group_size, n_table, N_seq, residue_seq, info_spectrum, &
		parents, offsprings, temp_p)

	!	evaluate the individuals in the group
	call evaluate_group(group_size, n_table, N_seq, cnn_row, &
		info_spectrum, connection_table, parents)

	!	calculate the pareto_order
	Pareto_set = 0
	num_set = 0
	obj_4 = 0.0
	do k1 = 1, group_size
		obj_4(k1, 1) = parents(k1)%n_good*1.0
		obj_4(k1, 2) = parents(k1)%n_used*1.0
		obj_4(k1, 3) = -parents(k1)%n_bad*1.0
		obj_4(k1, 4) = 10.0*parents(k1)%n_good-20.0*parents(k1)%n_bad-3.0*parents(k1)%n_edge&
		&+parents(k1)%n_used
	end do
	!
	num_obj = 4
	call pareto(obj_4(1:group_size, 1:num_obj), num_obj, group_size, num_set(1:group_size), Pareto_set(1:group_size), rank)
	!
	ind_k = 0
	do k1 = 1, rank
		ind_temp(1: num_set(k1)) = ind_k +(/(i, i=1, num_set(k1))/)
		ind_temp2(1: num_set(k1)) = Pareto_set(ind_temp)
		call crowd_dist(obj_4(ind_temp2(1:num_set(k1)), 1:num_obj), num_obj, num_set(k1), dist_s(1:num_set(k1)))
		ind_k = ind_k+num_set(k1)
		do k2 = 1, num_set(k1)
			Pareto_order(ind_temp2(k2)) = k1
			dist_group(ind_temp2(k2)) = dist_s(k2)
		end do
	end do
	!
	call itime(timearray)
	write(*,*) 'Begin the evolution at: ', timearray(1), ':', timearray(2), ':', timearray(3)
	!
	write(*,*) '--------------------------'
	write(*,*) 'After initialization: '
	write(*,*) 'number of good connections range: [', minval(parents(:)%n_good),&
	& ', ', maxval(parents(:)%n_good), ']'
	write(*,*) 'number of bad connections range: [', minval(parents(:)%n_bad),&
	& ', ', maxval(parents(:)%n_bad), ']'
	write(*,*) 'number of edges range: [', minval(parents(:)%n_edge),&
	& ', ', maxval(parents(:)%n_edge), ']'
	write(*,*) 'number of used peaks range: [', minval(parents(:)%n_used),&
	& ', ', maxval(parents(:)%n_used), ']'
	! ********************* evolution *****************************
	do ind_ns = 1, N_step
		lucky_ratio = exp(-(ind_ns-N_step/2+2)**2/10.0)
		if (lucky_ratio .lt. 1e-5)	lucky_ratio = 0
		do ind_nt = 1, N_try
			! initiate the evaluation values of the parents + offsprings
			n_good_all = 0
			n_bad_all = 0
			n_edge_all = 0
			n_used_all = 0
			do k1 = 1, group_size
				n_good_all(k1) = parents(k1)%n_good
				n_bad_all(k1) = parents(k1)%n_bad
				n_edge_all(k1) = parents(k1)%n_edge
				n_used_all(k1) = parents(k1)%n_used
			end do
			!
			! choose N individuals for mutation and crossover (N = new_group_size)
			if (pool_size .eq. group_size) then
				ind_pool = (/(k1, k1 = 1, pool_size)/)
			else
				k1 = 0
				ind_pool = 0
				do while (k1 < pool_size)
					call random_number(rand_1)
					call random_number(rand_2)
					p1 = int(rand_1*group_size)+1
					p2 = int(rand_2*group_size)+1
					do while (p1 .eq. p2)
						call random_number(rand_2)
						p2 = int(rand_2*group_size)+1
					end do
					k1 = k1+1
					if (Pareto_order(p1)>Pareto_order(p2)) then
						ind_pool(k1) = p2
					else if (Pareto_order(p1)<Pareto_order(p2)) then
						ind_pool(k1) = p1
					else if (dist_group(p1) < dist_group(p2)) then
						ind_pool(k1) = p2
					else
						ind_pool(k1) = p1
					end if
					if (all(ind_pool(1: k1-1) .ne. ind_pool(k1))) cycle
					ind_pool(k1) = 0
					k1 = k1-1
				end do
			end if
			! ******** generate new group ********
			ind_new_group = 0
			do while (ind_new_group < group_size)
				call random_number(rand_1)
				if (rand_1 < cross_rate) then
					! crossover
					call random_number(rand_2)
					call random_number(rand_3)
					p1 = int(rand_2*pool_size)+1
					p2 = int(rand_3*pool_size)+1
					do while (p1 .eq. p2)
						call random_number(rand_3)
						p2 = int(rand_3*pool_size)+1
					end do
					diff_residue_peak = 0
					do idx_table = 1, n_table
						diff_residue_peak = diff_residue_peak+abs(parents(ind_pool(p1))%rsd_pk(:, idx_table)&
						& - parents(ind_pool(p2))%rsd_pk(:, idx_table))
					end do
					if (any(diff_residue_peak .ne. 0)) then
						parents2(1) = parents(ind_pool(p1))
						parents2(2) = parents(ind_pool(p2))
						call crossover(parents2, info_spectrum, connection_table, N_seq, n_table, cnn_row, children2, err_sign)
						if (err_sign .eq. 0) then
							!
							do k2 = 1, 2
								residue_peak = children2(k2)%rsd_pk
								if (ind_new_group .eq. group_size) exit
								! delete the zero ones & identical one in the group
								call delete_z_i(residue_peak, parents, offsprings, ind_new_group, N_seq, n_table,&
								& group_size, del_note)
								if  (del_note .eq. 1) cycle
								ind_new_group = ind_new_group + 1
								offsprings(ind_new_group) = children2(k2)
								!
								n_good_all(group_size+ind_new_group) = children2(k2)%n_good
								n_bad_all(group_size+ind_new_group) = children2(k2)%n_bad
								n_edge_all(group_size+ind_new_group) = children2(k2)%n_edge
								n_used_all(group_size+ind_new_group) = children2(k2)%n_used
							end do
						end if
					end if
				end if
				if (ind_new_group .eq. group_size) exit
				call random_number(rand_1)
				if (rand_1 < mutate_rate) then
					! mutation
					call random_number(rand_2)
					p1 = int(rand_2*pool_size)+1
					call mutation(parents(ind_pool(p1)), info_spectrum, connection_table, N_seq, &
					& n_table, cnn_row, p_null, child1)
					!
					! delete the zero ones & identical one in the group
					call delete_z_i(child1%rsd_pk, parents, offsprings, ind_new_group, N_seq, n_table, &
					& group_size, del_note)
					if (del_note .eq. 0) then
						ind_new_group = ind_new_group + 1
						offsprings(ind_new_group) = child1
						!
						n_good_all(group_size+ind_new_group) = child1%n_good
						n_bad_all(group_size+ind_new_group) = child1%n_bad
						n_edge_all(group_size+ind_new_group) = child1%n_edge
						n_used_all(group_size+ind_new_group) = child1%n_used
					end if
				end if
				if (ind_new_group .eq. group_size) exit
				call random_number(rand_1)
				if (rand_1 < mutate_ad_rate) then
					! additional mutation
					call random_number(rand_2)
					p1 = int(rand_2*pool_size)+1
					if (all(parents(ind_pool(p1))%rsd_pk .eq. 0)) then
						err_sign = 1
					else
						call random_number(rand_3)
						if (rand_3 .ge. 0.5)	then
							call mutation_ad_2(parents(ind_pool(p1)), info_spectrum, connection_table, &
							& N_seq, n_table, cnn_row, child1, err_sign)
						else
							call mutation_ad_1(parents(ind_pool(p1)), info_spectrum, connection_table, &
							& N_seq, n_table, cnn_row, child1, err_sign)
						end if
					end if
					if (err_sign .eq. 1) cycle
					!
					! delete the zero ones & identical one in the group
					call delete_z_i(child1%rsd_pk, parents, offsprings, ind_new_group, N_seq, n_table, &
					& group_size, del_note)
					if (del_note .eq. 1) cycle
					ind_new_group = ind_new_group + 1
					offsprings(ind_new_group) = child1
					!
					n_good_all(group_size+ind_new_group) = child1%n_good
					n_bad_all(group_size+ind_new_group) = child1%n_bad
					n_edge_all(group_size+ind_new_group) = child1%n_edge
					n_used_all(group_size+ind_new_group) = child1%n_used
				end if
			end do

			! ******** choose the next group based on the fitness of the new group + old group*********
			num_new_group = ind_new_group
			size_all = num_new_group+group_size
			!
			note_rstr = ind_ns-num_free
			!
			call elitism(parents, offsprings, n_good_all, n_bad_all, n_edge_all, n_used_all, group_size, &
			& N_seq, n_table, note_rstr, temp_p, pareto_order, dist_group,lucky_ratio)
			parents = temp_p
		end do
		write(*,*) '--------------------------'
		write(*,*) 'Finish step: ', ind_ns
		write(*,*) 'Lucky ratio is: ', lucky_ratio
		write(*,*) 'number of good connections range: [', minval(parents(:)%n_good),&
		& ', ', maxval(parents(:)%n_good), ']'
		write(*,*) 'number of bad connections range: [', minval(parents(:)%n_bad),&
		& ', ', maxval(parents(:)%n_bad), ']'
		write(*,*) 'number of edges range: [', minval(parents(:)%n_edge),&
		& ', ', maxval(parents(:)%n_edge), ']'
		write(*,*) 'number of used peaks range: [', minval(parents(:)%n_used),&
		& ', ', maxval(parents(:)%n_used), ']'
	end do
	! **************************************** output results *****************************************
	n_good_all = 0
	n_bad_all = 0
	n_edge_all = 0
	n_used_all = 0
	do k1 = 1, group_size
		obj_4(k1,1) = parents(k1)%n_good*1.0
		obj_4(k1,2) = parents(k1)%n_good*10.0-parents(k1)%n_bad*20.0-parents(k1)%n_edge*3.0+parents(k1)%n_used*1.0
		obj_4(k1,3) = parents(k1)%n_used*1.0
		obj_4(k1,4) = -parents(k1)%n_bad*1.0
	end do
	!
	! sort the group based on the number of bad connections
	call pareto(obj_4(1:group_size, :), 4, group_size, num_set(1:group_size), Pareto_set(1:group_size), rank)
	ind_k = 0
	obj_4(:, 4) = -obj_4(:, 4)
	do k2 = 1, rank
		ind_temp(1:num_set(k2)) = (/(i, i=1,num_set(k2))/)
		call QsortC(obj_4(ind_k+1: ind_k+num_set(k2), 4), ind_temp(1:num_set(k2)))
		do k1 = 1, num_set(k2)
			temp_p(k1+ind_k) = parents(Pareto_set(ind_k+ind_temp(k1)))
			n_good_all(k1+ind_k) = parents(Pareto_set(ind_k+ind_temp(k1)))%n_good
			n_bad_all(k1+ind_k) = parents(Pareto_set(ind_k+ind_temp(k1)))%n_bad
			n_edge_all(k1+ind_k) = parents(Pareto_set(ind_k+ind_temp(k1)))%n_edge
			n_used_all(k1+ind_k) = parents(Pareto_set(ind_k+ind_temp(k1)))%n_used
			Pareto_order(k1+ind_k) = k2
		end do
		ind_k = ind_k+num_set(k2)
	end do
	parents = temp_p
	!
	write(*,*) 	'--------- final results ------------'
	write(*,*)	'Ng,    Nb,    Ne,    Nu,    Pareto order'
	do k = 1, group_size
		write(*, '(i4,2x,i4,3x,i4,3x,i4,8x,i3)')	parents(k)%n_good, parents(k)%n_bad, parents(k)%n_edge, &
		&parents(k)%n_used, Pareto_order(k)
	end do
	!
	call itime(new_timearray)
	write(*,*)	'Stop at: ', new_timearray(1), ':', new_timearray(2), ':', new_timearray(3)
	k = new_timearray(3)+new_timearray(2)*60+new_timearray(1)*3600&
	&-timearray(3)-timearray(2)*60-timearray(1)*3600
	p1 = k/3600
	p2 = (k-p1*3600)/60
	p3 = k-p1*3600-p2*60
	write(*,*) 	'Cost time: ', p1, ' hours ', p2, ' minutes ', p3, ' seconds.'
	!
	! Output results, generate data files
	open(10, file = file_out_data)
	write(10,*)	 'Output_data'
	do i = 1, group_size
		do k1 = 1, n_table
			write(10,*) parents(i)%rsd_pk(:, k1)
		end do
	end do
	close(10)
	!
	allocate(freq(n_table, info_spectrum(1)%num_freq))
	write(fmt_s,'(ai1ai1aa)') '(i4, a2,(',n_table,'(i4,',info_spectrum(1)%num_freq,'f8.1)),a6)'
	i = 0
	do k1 = 1, group_size
		i = i+1
		p1 = i/100
		p2 = (i-p1*100)/10
		p3 = i-p1*100-p2*10
		open(11, file = trim(adjustL(file_out_tables))//achar(48+p1)//achar(48+p2)//achar(48+p3)//'.txt')
		write(11,"('Ng, Nb, Ne, Nu and Pareto order: 'I4,2x,I4,2x,I4,2x,I4,2x,I3)") parents(i)%n_good, parents(i)%n_bad, &
		&parents(i)%n_edge, parents(i)%n_used, Pareto_order(k1)
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
			write(11,fmt_s) k2, residue_seq(k2), (parents(k1)%rsd_pk(k2, idx_table), freq(idx_table,:), idx_table = 1, n_table), &
			&prsd_str(1:k4)
		end do
		close(11)
	end do
	!
	deallocate(freq)
end program
