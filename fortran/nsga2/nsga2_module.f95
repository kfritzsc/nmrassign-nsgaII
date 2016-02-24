module nsga2_module

	public :: crossover, mutation, mutation_ad_1, mutation_ad_2, &
		evaluate_peak_rsd, pareto, evaluate_idv, crowd_dist, delete_z_i

	private ::	shuffle_y


	type	info_spectra
	!
	! Spectra Information
	! ===================
	!

	! peak numbers of the spectrums
	integer ::	num_peak = 0
	! numbers of the kinds of chemical shifts
	integer ::	num_freq = 0
	! chemical shifts
	real, allocatable ::	ch_shifts(:,:)
	! resolution (half width)
	real, allocatable :: e_width(:,:)
	! degeneracy number
	integer, allocatable ::	degeneracy(:)
	! possible residue corresponding to each peak
	character*20, allocatable	::	poss_rsd(:)
	! the used times of each peak
	character*100, allocatable	::	poss_rsd_str(:)
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
	!
	integer, allocatable:: 	rsd_pk(:, :)
	integer :: 	n_good = 0
	integer ::	n_bad = 0
	integer ::	n_edge = 0
	integer ::	n_used = 0

	end type


	contains


	! Private subroutines


	subroutine shuffle_y(a, len_a)
	!
	! Shuffle Subroutine (private)
	! ==================
	! Random shuffle a.
	!
	! Input
	! -----
	! a
	! len_a
	!
	! Output
	! ------
	! a
	!

	implicit none

	integer, intent(in):: len_a
	integer, intent(inout):: a(len_a)

	integer	k1, k2, k3, len_k
	real		rand_1

	len_k = len_a
	do k1 = 1, len_a
		call random_number(rand_1)
		k2 = int(rand_1*len_k)+1
		k3 = a(k2+k1-1)
		a(k2+k1-1) = a(k1)
		a(k1) = k3
		len_k = len_k-1
	end do
	end subroutine shuffle_y


	! Public subroutines


	subroutine crossover(parents2, info_spectrum, connection_table, N_seq, &
		n_table, cnn_row, children2, err_sign)
	!
	! cross-over subroutine
	! =====================
	!
	! Choose a random residue p and exchanges the rd_pk(p, ind_table) in parent1
	! with rsd_pk(p, ind_table) in parent2. The subroutine also evaluates
	! the two children.
	!
	! Input
	! ----------
	! parent2: two individuals in the generation group
	!
	! Output
	! -------
	! children2: two children
	!
	implicit none

	type(idv), intent(in)	::	parents2(2)
	integer, intent(in)	::	N_seq, n_table, cnn_row
	type(info_spectra), intent(in) ::	info_spectrum(n_table)
	integer, intent(in) ::	connection_table(cnn_row, 6)
	type(idv), intent(out)	::	children2(2)
	integer, intent(out)	::	err_sign

	integer	::	rsd_peak_2(N_seq, n_table, 2), new_asg_2(2)
	integer	new_asg, old_asg
	integer	ind_rsd, new_idx, k, p1, p2, k2, num_pos, idx_table
	integer	ngood_new, ngood_old
	integer nbad_new, nbad_old
	integer nedge_new, nedge_old
	integer nused_new, nused_old
	integer	stop_sign, i
	integer	:: 	diff_rsd_abs(N_seq)
	integer	::	ind_a(N_seq)
	real	rand_1

	stop_sign = 0
	err_sign = 0

	rsd_peak_2(:, :, 1) = parents2(1)%rsd_pk
	rsd_peak_2(:, :, 2) = parents2(2)%rsd_pk

	call random_number(rand_1)
	idx_table = int(n_table*rand_1) + 1
	diff_rsd_abs = abs(rsd_peak_2(:, idx_table, 1)-rsd_peak_2(:, idx_table, 2))
	num_pos = count(diff_rsd_abs .ne. 0)
	if (num_pos .lt. 2) then
	 	err_sign = 1
		return
	end if

	p1 = 0
	do k = 1, N_seq
		if (diff_rsd_abs(k) .eq. 0) 	cycle
		p1 = p1+1
		ind_a(p1) = k
	end do

	call shuffle_y(ind_a(1:num_pos), num_pos)
	do i = 1, num_pos
		ind_rsd = ind_a(i)
		if (info_spectrum(idx_table)%num_poss_peak(ind_rsd) .lt. 1) then
			new_asg_2 = rsd_peak_2(ind_rsd, idx_table, :)
			cycle
		end if

		new_asg_2(2) = rsd_peak_2(ind_rsd, idx_table,1)
		new_asg_2(1) = rsd_peak_2(ind_rsd, idx_table, 2)
		new_idx = ind_rsd

		do k = 1, 2
			if (new_asg_2(k) .ne. 0 &
			    .and. new_asg_2(k) .ne. rsd_peak_2(ind_rsd, idx_table, k)) then
				p1 = 0
				do p2 = 1, N_seq
					if (rsd_peak_2(p2, idx_table, k) .eq. new_asg_2(k)) then
						p1 = p1+1
					end if
				end do
				if (p1 .eq. info_spectrum(idx_table)%degeneracy(new_asg_2(k))) then
					new_asg_2(k) = rsd_peak_2(ind_rsd, idx_table, k)
				end if
			end if
		end do

		if (new_asg_2(1) .eq. rsd_peak_2(ind_rsd, idx_table, 1) &
		    .or. new_asg_2(2) .eq. rsd_peak_2(ind_rsd, idx_table, 2)) then
			stop_sign = 0
		else
			stop_sign = 1
			exit
		end if
	end do

	if (stop_sign .eq. 1) then
		err_sign = 0
	else
		err_sign = 1
	end if

	if (err_sign .eq. 1)	return

	!	Evaluate child1 and child2
	do k2 = 1, 2
		! evaluate new group member
		children2(k2) = parents2(k2)

		old_asg = parents2(k2)%rsd_pk(new_idx, idx_table)
		new_asg = new_asg_2(k2)

		call evaluate_peak_rsd (children2(k2)%rsd_pk, info_spectrum, &
			connection_table, n_table, N_seq, cnn_row, idx_table, new_idx, &
			new_asg, ngood_new, nbad_new, nedge_new)

		call evaluate_peak_rsd (children2(k2)%rsd_pk, info_spectrum, &
			connection_table, n_table, N_seq, cnn_row, idx_table, new_idx, &
			old_asg, ngood_old, nbad_old, nedge_old)

		if (old_asg .ne. 0) then
			nused_old = 1
		else
			nused_old = 0
		end if
		if (new_asg .ne. 0) then
			nused_new = 1
		else
			nused_new = 0
		end if

		children2(k2)%n_good = children2(k2)%n_good-ngood_old+ngood_new
		children2(k2)%n_bad = children2(k2)%n_bad-nbad_old+nbad_new
		children2(k2)%n_edge = children2(k2)%n_edge-nedge_old+nedge_new
		children2(k2)%n_used = children2(k2)%n_used+nused_new-nused_old
		children2(k2)%rsd_pk(new_idx, idx_table) = new_asg
	end do

	end subroutine crossover


	subroutine mutation(parent1, info_spectrum, connection_table, N_seq, &
		n_table, cnn_row, p_null, child1)
	!
	! Mutation #0 Subroutine
	! ===================
	! Choose a random residue p and change the rsd_peak_k(p, ind_table) value
	! to make a new assignment where rsd_peak_k is the kth individual of the
	! generation group
	!
	! Input
	! -----
	! N_seq, n_table, cnn_row
	!
	! Output
	! ------
	! new_asg
	! new_idx

	implicit none

	integer, intent(in) ::	N_seq, n_table, cnn_row
	real, intent(in) ::	p_null
	type(idv), intent(in) ::	parent1
	type(info_spectra), intent(in) ::	info_spectrum(n_table)
	integer, intent(in) ::	connection_table(cnn_row, 6)
	type(idv), intent(out) ::	child1

	integer	::	rsd_peak_k(N_seq, n_table)
	real	rand_1, rand_2
	integer	ind_rsd, ind_peak, idx_table,new_idx, new_asg, old_asg
	integer	ngood_new, ngood_old
	integer nbad_new, nbad_old
	integer nedge_new, nedge_old
	integer nused_new, nused_old
	integer	num_peak, p1, p2, stop_sign

	stop_sign = 0

	rsd_peak_k = parent1%rsd_pk
	do while (stop_sign .eq. 0)
		call random_number(rand_1)
		ind_rsd = int(rand_1*N_seq) + 1
		!
		call random_number(rand_1)
		idx_table = int(n_table*rand_1) + 1
		new_asg = rsd_peak_k(ind_rsd, idx_table)
		if (info_spectrum(idx_table)%num_poss_peak(ind_rsd) .lt. 1) cycle
		num_peak = info_spectrum(idx_table)%num_poss_peak(ind_rsd)
		do while (new_asg .eq. rsd_peak_k(ind_rsd, idx_table))
			call random_number(rand_2)
			if (rand_2 < p_null) then
				new_asg = 0
			else
				call random_number(rand_2)
				ind_peak = int(num_peak*rand_2) + 1
				new_asg = info_spectrum(idx_table)%peak_seq(ind_peak, ind_rsd)
			end if
		end do
		new_idx = ind_rsd
		if (new_asg .ne. 0 &
		    .and. new_asg .ne. rsd_peak_k(ind_rsd, idx_table)) then
			p1 = 0
			do p2 = 1, N_seq
				if (rsd_peak_k(p2, idx_table) .eq. new_asg) p1 = p1+1
			end do
			if (p1 .eq. info_spectrum(idx_table)%degeneracy(new_asg)) then
				new_asg = rsd_peak_k(ind_rsd, idx_table)
			end if
		end if
		if (new_asg .eq. rsd_peak_k(ind_rsd, idx_table)) then
			stop_sign = 0
		else
			stop_sign = 1
		end if
	end do

	!	evaluate the new child
	child1 = parent1
	old_asg = rsd_peak_k(new_idx, idx_table)
	call evaluate_peak_rsd(rsd_peak_k, info_spectrum, connection_table, &
		n_table, N_seq, cnn_row, idx_table, new_idx, new_asg, ngood_new, &
		nbad_new, nedge_new)
	call evaluate_peak_rsd(rsd_peak_k, info_spectrum, connection_table, &
		n_table, N_seq, cnn_row, idx_table, new_idx, old_asg, ngood_old ,&
		nbad_old, nedge_old)

	if (old_asg .ne. 0) then
		nused_old = 1
	else
		nused_old = 0
	end if
	if (new_asg .ne. 0) then
		nused_new = 1
	else
		nused_new = 0
	end if

	child1%n_good = child1%n_good-ngood_old+ngood_new
	child1%n_bad = child1%n_bad-nbad_old+nbad_new
	child1%n_edge = child1%n_edge-nedge_old+nedge_new
	child1%n_used = child1%n_used-nused_old+nused_new
	rsd_peak_k(new_idx, idx_table) = new_asg
	child1%rsd_pk = rsd_peak_k

	end subroutine mutation


	subroutine mutation_ad_1(parent1, info_spectrum, connection_table, N_seq, &
		n_table, cnn_row, child1, err_sign)
	!
	! Mutation #1 subroutine
	! ==============================
	!
	! Choose a random residue p1 and another random residue p2 then exchange
	! rsd_peak_k(p1, ind_table) and rsd_peak_k(p2, ind_table) where rsd_peak_k
	! is the kth individual of the generation group
	!
	! Output
	! ------
	! new_asg
	! new_idx
	!

	implicit none

	integer, intent(in) ::	N_seq, n_table, cnn_row
	type(idv), intent(in) ::	parent1
	type(info_spectra), intent(in) ::	info_spectrum(n_table)
	integer, intent(in) ::	connection_table(cnn_row, 6)
	type(idv), intent(out) ::	child1
	integer, intent(out) ::	err_sign

	integer	::	new_asg_2(2), new_idx_2(2), old_asg_2(2)
	integer	::	new_asg, old_asg
	integer	::	ngood_new, ngood_old
	integer nbad_new, nbad_old
	integer	nedge_new, nedge_old
	integer	nused_new, nused_old
	integer	::	rsd_peak_k(N_seq, n_table)
	real	rand_1
	integer	ind_rsd, idx_table
	integer	p1, p2, p3
	integer	num_pos, k, p, k2
	integer	::	ind_a(N_seq)

	err_sign = 0
	new_idx_2 = 0
	new_asg_2 = 0
	old_asg_2 = 0

	rsd_peak_k = parent1%rsd_pk
	call random_number(rand_1)
	idx_table = int(n_table*rand_1)+1
	num_pos = count(rsd_peak_k(:, idx_table) .ne. 0)
	if (num_pos .eq. 0) then
		err_sign = 1
		return
	end if
	p = 0
	do k = 1, N_seq
		if (rsd_peak_k(k, idx_table) .ne. 0) then
			p = p+1
			ind_a(p) = k
		end if
		if (p .eq. num_pos) exit
	end do
	call shuffle_y(ind_a(1:num_pos), num_pos)

	do k = 1, num_pos
		ind_rsd = ind_a(k)
		new_idx_2(1) = ind_rsd
		old_asg_2(1) = rsd_peak_k(ind_rsd, idx_table)
		if (info_spectrum(idx_table)%num_poss_peak(ind_rsd) .lt. 1) then
			new_asg_2(1) = old_asg_2(1)
			cycle
		end if
		p1 = info_spectrum(idx_table)%num_poss_rsd(old_asg_2(1))
		if (p1 < 2) then
			new_asg_2(1) = old_asg_2(1)
			cycle
		end if
		call random_number(rand_1)
		p2 = int(p1*rand_1) + 1
		p3 = info_spectrum(idx_table)%prsd_peak(p2, old_asg_2(1))
		do while (p3 .eq. ind_rsd)
			call random_number(rand_1)
			p2 = int(p1*rand_1) + 1
			p3 = info_spectrum(idx_table)%prsd_peak(p2, old_asg_2(1))
		end do
		new_idx_2(2) = p3
		old_asg_2(2) = rsd_peak_k(new_idx_2(2), idx_table)
		if (info_spectrum(idx_table)%num_poss_peak(new_idx_2(2)) .lt. 1) then
			new_asg_2(2) = old_asg_2(2)
			cycle
		end if
		new_asg_2(1) = old_asg_2(2)
		new_asg_2(2) = old_asg_2(1)
		if (new_asg_2(1) .ne. 0) then
			if (all(info_spectrum(idx_table)%prsd_peak(:, new_asg_2(1)) &
			    .ne. new_idx_2(1))) 	then
				new_asg_2(1) = 0
			end if
		end if
		exit
	end do
	if (all((new_asg_2-old_asg_2) .eq. 0) .or. new_idx_2(2) .eq. 0) then
		err_sign = 1
	else
		err_sign = 0
	end if

	if (err_sign .eq. 1)	return

	!	evaluate child1
	child1 = parent1
	do k2 = 1, 2
		old_asg = rsd_peak_k(new_idx_2(k2), idx_table)
		new_asg = new_asg_2(k2)
		if (new_asg .eq. old_asg) cycle
		call evaluate_peak_rsd(rsd_peak_k, info_spectrum, connection_table, &
			n_table, N_seq, cnn_row, idx_table, new_idx_2(k2), new_asg, &
			ngood_new, nbad_new, nedge_new)

		call evaluate_peak_rsd(rsd_peak_k, info_spectrum, connection_table, &
		n_table, N_seq, cnn_row, idx_table, new_idx_2(k2), old_asg, ngood_old,&
		nbad_old, nedge_old)

		if (old_asg .ne. 0) then
			nused_old = 1
		else
			nused_old = 0
		end if
		if (new_asg .ne. 0) then
			nused_new = 1
		else
			nused_new = 0
		end if

		rsd_peak_k(new_idx_2(k2), idx_table) = new_asg
		child1%n_good = child1%n_good-ngood_old+ngood_new
		child1%n_bad = child1%n_bad-nbad_old+nbad_new
		child1%n_edge = child1%n_edge-nedge_old+nedge_new
		child1%n_used = child1%n_used-nused_old+nused_new
	end do
	child1%rsd_pk = rsd_peak_k

	end subroutine mutation_ad_1


	subroutine mutation_ad_2(parent1, info_spectrum, connection_table, N_seq, &
		n_table, cnn_row, child1, err_sign)
	!
	! Mutation #2 subroutine
	! ======================
	!
	! Choose a random residue p1 and another random residue p2, exchange
	! rsd_peak_k(p1, ind_table) and rsd_peak_k(p2, ind_table) where rsd_peak_k
	! is the kth individual of the generation group
	!
	! Input
	! ----
	! parent1
	! N_seq
	! n_table
	! cnn_row
	!
	! Output
	! ------
	! child1
	! err_sign
	!

	implicit none
	type(idv), intent(in) ::	parent1
	type(info_spectra), intent(in) ::	info_spectrum(n_table)
	integer, intent(in) ::	N_seq, n_table, cnn_row
	integer, intent(in) ::	connection_table(cnn_row, 6)
	type(idv), intent(out) ::	child1
	integer, intent(out) ::	err_sign

	integer	::	new_asg_2n(n_table, 2), new_idx_2(2), old_asg_2n(n_table, 2)
	integer	::	new_asg, old_asg
	integer	::	ngood_new, ngood_old
	integer	::	nbad_new, nbad_old
	integer	::	nedge_new, nedge_old
	integer	::	nused_new, nused_old
	integer	::	rsd_peak_k(N_seq, n_table)
	real	rand_1
	integer	ind_rsd, idx_table
	integer	p1, p2, p3
	integer	num_pos, k, p, k2, i
	integer	::	ind_a(N_seq)
	integer	::	ind_t(n_table)

	err_sign = 0
	new_idx_2 = 0
	new_asg_2n = 0
	old_asg_2n = 0
	!
	rsd_peak_k = parent1%rsd_pk
	ind_t = (/(i, i = 1, n_table)/)
	call shuffle_y(ind_t, n_table)
	!
	idx_table = ind_t(1)
	num_pos = count(rsd_peak_k(:, idx_table) .ne. 0)
	if (num_pos .eq. 0) then
		err_sign = 1
		return
	end if
	p = 0
	do k = 1, N_seq
		if (rsd_peak_k(k, idx_table) .ne. 0) then
			p = p+1
			ind_a(p) = k
		end if
		if (p .eq. num_pos) exit
	end do
	call shuffle_y(ind_a(1:num_pos), num_pos)
	do k = 1, num_pos
		ind_rsd = ind_a(k)
		new_idx_2(1) = ind_rsd
		old_asg_2n(idx_table,1) = rsd_peak_k(ind_rsd, idx_table)
		if (info_spectrum(idx_table)%num_poss_peak(ind_rsd) .lt. 1) then
			new_asg_2n(idx_table,1) = old_asg_2n(idx_table,1)
			cycle
		end if
		p1 = info_spectrum(idx_table)%num_poss_rsd(old_asg_2n(idx_table,1))
		if (p1 < 2) then
			new_asg_2n(idx_table,1) = old_asg_2n(idx_table,1)
			cycle
		end if
		call random_number(rand_1)
		p2 = int(p1*rand_1) + 1
		p3 = info_spectrum(idx_table)%prsd_peak(p2, old_asg_2n(idx_table,1))
		do while (p3 .eq. ind_rsd)
			call random_number(rand_1)
			p2 = int(p1*rand_1) + 1
			p3 = info_spectrum(idx_table)%prsd_peak(p2, old_asg_2n(idx_table,1))
		end do
		new_idx_2(2) = p3
		old_asg_2n(idx_table, 2) = rsd_peak_k(new_idx_2(2), idx_table)
		if (info_spectrum(idx_table)%num_poss_peak(new_idx_2(2)) .lt. 1) then
			new_asg_2n(idx_table,2) = old_asg_2n(idx_table,2)
			cycle
		end if
		new_asg_2n(idx_table, 2) = old_asg_2n(idx_table, 1)
		new_asg_2n(idx_table, 1) = old_asg_2n(idx_table, 2)
		if (new_asg_2n(idx_table, 1) .ne. 0) then
			if (all(info_spectrum(idx_table)%prsd_peak(:, &
			 	new_asg_2n(idx_table, 1)).ne. new_idx_2(1))) 	then
				new_asg_2n(idx_table, 1) = 0
			end if
		end if
		exit
	end do
	if (all((new_asg_2n(idx_table,:)-old_asg_2n(idx_table,:)) .eq. 0) &
	    .or. new_idx_2(2) .eq. 0) then
		err_sign = 1
		return
	else
		err_sign = 0
	end if
	do k = 2, n_table
		p = ind_t(k)
		old_asg_2n(p, 1) = rsd_peak_k(new_idx_2(1), p)
		old_asg_2n(p, 2) = rsd_peak_k(new_idx_2(2), p)
		new_asg_2n(p, 1) = old_asg_2n(p, 2)
		new_asg_2n(p, 2) = old_asg_2n(p, 1)
		if (new_asg_2n(p, 1) .ne. 0) then
			if (all(info_spectrum(p)%prsd_peak(:, new_asg_2n(p,1)) &
			    .ne. new_idx_2(1))) 	then
				new_asg_2n(p, 1) = 0
			end if
		end if
		if (new_asg_2n(p, 2) .ne. 0) then
			if (all(info_spectrum(p)%prsd_peak(:, new_asg_2n(p,2)) &
			    .ne. new_idx_2(2))) 	then
				new_asg_2n(p, 2) = 0
			end if
		end if
		if (all(new_asg_2n(p, :) .eq. 0))	then
			new_asg_2n(p, :) = old_asg_2n(p, :)
		end if
	end do
	!
	!	evaluate child1
	child1 = parent1
	do idx_table = 1, n_table
		do k2 = 1, 2
			old_asg = rsd_peak_k(new_idx_2(k2), idx_table)
			new_asg = new_asg_2n(idx_table, k2)
			if (new_asg .eq. old_asg) cycle
			call evaluate_peak_rsd(rsd_peak_k, info_spectrum, &
				connection_table, n_table, N_seq, cnn_row, idx_table, &
				new_idx_2(k2), new_asg, ngood_new, nbad_new, nedge_new)

			call evaluate_peak_rsd(rsd_peak_k, info_spectrum, &
				connection_table, n_table, N_seq, cnn_row, idx_table, &
				new_idx_2(k2), old_asg, ngood_old, nbad_old, nedge_old)

			if (old_asg .ne. 0) then
				nused_old = 1
			else
				nused_old = 0
			end if
			if (new_asg .ne. 0) then
				nused_new = 1
			else
				nused_new = 0
			end if

			rsd_peak_k(new_idx_2(k2), idx_table) = new_asg
			child1%n_good = child1%n_good-ngood_old+ngood_new
			child1%n_bad = child1%n_bad-nbad_old+nbad_new
			child1%n_edge = child1%n_edge-nedge_old+nedge_new
			child1%n_used = child1%n_used-nused_old+nused_new
		end do
	end do
	child1%rsd_pk = rsd_peak_k

	end subroutine mutation_ad_2


	subroutine evaluate_peak_rsd (residue_peak, info_spectrum, &
		connection_table, n_table, N_seq, cnn_row, table_type, ind_rsd, asg, &
		ngood, nbad, nedge)
	!
	! evaluate_peak_rsd subroutine
	! ---------------------------
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

	end subroutine


	subroutine pareto(obj, num_obj, size_all, num_set, P_set, rank)

	!
	! pareto subroutine
	! =================
	!
	! Calculate the pareto order of the whole group(new generation +
	! old generation)
	!
	! Input
	! -----
	! num_obj
	! size_all
	! obj(size_all, num_obj)
	!
	! Output
	! ------
	! num_set(size_all)
	! P_set(size_all)
	! rank
	!

	implicit none

	integer, intent(in)	::	num_obj, size_all
	real, intent(in)	::	obj(size_all, num_obj)
	integer, intent(out)	::	num_set(size_all), P_set(size_all)
	integer, intent(out)	::	rank

	integer	::	np(size_all), Sp(size_all, size_all)
	integer	k1, k2, ind_1, ind_2, ind_g, num_ss, ind_temp
	real	::	a(num_obj)
	integer	e1

	np = 0

	ind_g = 0
	num_set = 0
	Sp = 0
	P_set = 0
	do k1 = 1, size_all
		ind_1 = 1
		do k2 = 1, size_all
			if (k2 .eq. k1) cycle
			! if k2 dominated k1, all the objectives of k2 are larger than k1
			! if k1 is dominated by k2, then np(k1) = np(k1)+1
			! if k2 is dominated by k1, add k2 into the Sp(k1,:) array
			a = obj(k2, :) - obj(k1, :)
			if (all(a .eq. 0)) cycle
	 		if (all(a .ge. 0))	then
				np(k1) = np(k1)+1
			end if
	 		if (all(a .le. 0)) then
				Sp(k1, ind_1) = k2
				ind_1 = ind_1+1
				cycle
			end if
		end do
		if (np(k1) .eq. 0)	then
			ind_g = ind_g+1
			P_set(ind_g) = k1
			num_set(1) = num_set(1)+1
		end if
	end do
	!
	rank = 1
	ind_temp = 0
	e1=0
	do while (ind_g .lt. size_all)
		e1= e1+1
		if (e1>size_all) then
			write(*,*) 'error pareto'
			stop
			return
		end if
		num_ss = 0
		do k1 = 1, num_set(rank)
			ind_1 = P_set(ind_temp+k1)
			do k2 = 1, size_all
				if (Sp(ind_1, k2) .eq. 0)  exit
				ind_2 = Sp(ind_1, k2)
				np(ind_2) = np(ind_2)-1
				if (np(ind_2) .eq. 0) then
					ind_g = ind_g+1
					P_set(ind_g) = ind_2
					num_ss = num_ss+1
				end if
			end do
		end do
		ind_temp = ind_temp+num_set(rank)
		rank = rank+1
		num_set(rank) = num_ss
	end do

	end subroutine pareto


	subroutine pareto_rstr(obj, rstr, num_obj, size_all, num_set, P_set, rank)
	!
	! Pareto with restraints subroutine
	! =================================
	!
	! Calculates the pareto order of the whole group (new generation +
	! old generation) with restraints.
	!
	! Input
	! -----
	! num_obj, size_all, obj(size_all, num_obj), rstr(size_all)
	!
	! Output
	! ------
	! num_set(size_all), P_set(size_all), rank
	!

	implicit none

	integer, intent(in)	::	num_obj, size_all
	real, intent(in)	::	obj(size_all, num_obj)
	integer, intent(in)	::	rstr(size_all)
	integer, intent(out)	::	num_set(size_all), P_set(size_all)
	integer, intent(out)	::	rank

	integer	::	np(size_all), Sp(size_all, size_all)
	integer	k1, k2, ind_1, ind_2, ind_g, num_ss, ind_temp
	real	::	a(num_obj)
	integer	e1

	np = 0

	ind_g = 0
	num_set = 0
	Sp = 0
	P_set = 0
	do k1 = 1, size_all
		ind_1 = 1
		do k2 = 1, size_all
			if (k2 .eq. k1) cycle
			! if k2 dominated k1, all the objectives of k2 are larger than k1
			! if k1 is dominated by k2, then np(k1) = np(k1)+1
			! if k2 is dominated by k1, add k2 into the Sp(k1,:) array
			if (rstr(k1) > rstr(k2)) then
				np(k1) = np(k1)+1
				cycle
			elseif (rstr(k1) < rstr(k2)) then
				Sp(k1, ind_1) = k2
				ind_1 = ind_1+1
				cycle
			end if
			a = obj(k2, :) - obj(k1, :)
			if (all(a .eq. 0)) cycle
	 		if (all(a .ge. 0))	then
				np(k1) = np(k1)+1
				cycle
			end if
	 		if (all(a .le. 0)) then
				Sp(k1, ind_1) = k2
				ind_1 = ind_1+1
				cycle
			end if
		end do
		if (np(k1) .eq. 0)	then
			ind_g = ind_g+1
			P_set(ind_g) = k1
			num_set(1) = num_set(1)+1
		end if
	end do
	!
	rank = 1
	ind_temp = 0
	e1=0
	do while (ind_g .lt. size_all)
		e1= e1+1
		if (e1>size_all) then
			write(*,*) 'error pareto'
			stop
			return
		end if
		num_ss = 0
		do k1 = 1, num_set(rank)
			ind_1 = P_set(ind_temp+k1)
			do k2 = 1, size_all
				if (Sp(ind_1, k2) .eq. 0)  exit
				ind_2 = Sp(ind_1, k2)
				np(ind_2) = np(ind_2)-1
				if (np(ind_2) .eq. 0) then
					ind_g = ind_g+1
					P_set(ind_g) = ind_2
					num_ss = num_ss+1
				end if
			end do
		end do
		ind_temp = ind_temp+num_set(rank)
		rank = rank+1
		num_set(rank) = num_ss
	end do

	end subroutine pareto_rstr


	subroutine crowd_dist(obj, num_obj, size_s, dist_s)
	!
	! Crowding-distance subroutine
	! ============================
	!
	! Calculates the crowding-distance of the chosen set.
	!
	! Uses
	! ----
	! qsort_c_module
	!
	! Input
	! -----
	! num_obj, size_s, obj(size_s, num_obj)
	!
	! Output
	! ------
	! num_set(size_all), P_set(size_all), rank
	!

	use qsort_c_module
	implicit none

	integer,intent(in) ::	num_obj, size_s
	real,intent(in) ::	obj(size_s, num_obj)
	real,intent(out) ::	dist_s(size_s)

	real :: dist_temp(size_s)
	real :: a(size_s)
	integer	k, k1
	integer :: 	ind_array(size_s), ind_array_temp(size_s)
	real diff
	!
	dist_s = 0
	do k = 1, num_obj
		a = obj(:, k)
		ind_array = (/(k, k=1, size_s)/)

		dist_temp = 0
		ind_array_temp = ind_array
		call QsortC(a, ind_array_temp)
		dist_temp(ind_array_temp(1)) = 1.0
		dist_temp(ind_array_temp(size_s)) = 1.0
		if (size_s .le. 2)	then
			dist_s = dist_temp
			return
		end if
		diff = (a(size_s)-a(1))*1.0
		if (diff .eq. 0.0) then
			dist_temp = 0.0
		else
			do k1 = 2, size_s-1
				dist_temp(ind_array_temp(k1)) = dist_temp(ind_array_temp(k1)) &
					+(a(k1+1)-a(k1-1))/diff
			end do
		end if
		dist_s = dist_temp+dist_s
	end do

	end subroutine crowd_dist


	subroutine elitism(group_1, group_2, n_good_all, n_bad_all, n_edge_all, &
		n_used_all, size_group, N_seq, n_table, note_rstr, group_new, &
		pareto_order, dist_group, lucky_ratio)
	!
	! Elitism subroutine
	! ============================
	!
	! Select N individuals from N+N individuals
	!
	! Uses
	! ----
	! qsort_c_module
	!
	! Input
	! -----
	! size_group, N_seq, n_table, note_rstr, group_1(size_group),
	! group_2(size_group), n_good_all(2*size_group), n_bad_all(2*size_group)
	! n_edge_all(2*size_group), n_used_all(2*size_group), lucky_ratio
	!
	! Output
	! ------
	! group_new(size_group), pareto_order(size_group), dist_group(size_group)
	!

	use qsort_c_module
	implicit none

	integer, intent(in) :: size_group, N_seq, n_table, note_rstr
	type(idv), intent(in) :: group_1(size_group), group_2(size_group)
	integer, intent(in) :: n_good_all(2*size_group)
	integer, intent(in) :: n_bad_all(2*size_group)
	integer, intent(in) :: n_edge_all(2*size_group)
	integer, intent(in) :: n_used_all(2*size_group)
	real, intent(in) :: lucky_ratio
	type(idv), intent(out) :: group_new(size_group)
	integer, intent(out) ::	pareto_order(size_group)
	real, intent(out) :: dist_group(size_group)

	integer	i, ind_k, k1, k2, k, rank
	integer	size_all
	integer :: ind_temp(2*size_group)
	integer :: ind_temp2(2*size_group)
	integer :: ind_dist(2*size_group)
	integer :: Pareto_set(2*size_group), num_set(2*size_group)
	real ::	dist_s(2*size_group), dist_temp(2*size_group)
	real ::	obj(2*size_group, 5)
	integer :: temp_one(2*size_group)
	integer	num_obj

	pareto_order = 0
	dist_group = 0

	size_all = 2*size_group

	temp_one = 10*n_good_all-20*n_bad_all-3*n_edge_all+n_used_all
	if (note_rstr .ge. 0) then
		obj(:, 1) = n_good_all*1.0
		obj(:, 2) = n_used_all*1.0
		call random_number(obj(:,3))
		obj(:, 3) = obj(:,3) * lucky_ratio &
			+ (1-lucky_ratio)*(temp_one/maxval(abs(temp_one)))
		num_obj = 3
	else
		obj(:, 1) = n_good_all*1.0
		obj(:, 2) = n_used_all*1.0
		obj(:, 3) = -n_bad_all*1.0
		call random_number(obj(:,4))
		obj(:, 4) = obj(:,4) * lucky_ratio &
			+ (1-lucky_ratio)*(temp_one/maxval(abs(temp_one)))
		num_obj = 4
	end if

	! pareto order
	if (note_rstr .lt. 0) then
		call pareto(obj(:, 1:num_obj), num_obj, size_all, num_set, Pareto_set, &
			rank)
	else
		call pareto_rstr(obj(:, 1:num_obj), n_bad_all, num_obj, size_all, &
			num_set, Pareto_set, rank)
	end if
	ind_k = 0
	do k1 = 1, rank
		! calculate the crowd-distance
		ind_temp(1: num_set(k1)) = ind_k +(/(i, i=1, num_set(k1))/)
		ind_temp2(1: num_set(k1)) = Pareto_set(ind_temp(1:num_set(k1)))
		call crowd_dist(obj(ind_temp2(1: num_set(k1)), 1:num_obj), num_obj, &
			num_set(k1), dist_s(1:num_set(k1)))

		if (size_group .ge. (num_set(k1)+ind_k)) then
			do k2 = 1, num_set(k1)
				k = Pareto_set(ind_k + k2)
				if (k .le. size_group) then
					group_new(ind_k+k2) = group_1(k)
				else
					group_new(ind_k+k2) = group_2(k-size_group)
				end if
				pareto_order(ind_k+k2) = k1
				dist_group(ind_k+k2) = dist_s(k2)
			end do
			ind_k = ind_k+num_set(k1)
		else if (size_group .eq. ind_k) then
			exit
		else
			ind_dist(1:num_set(k1)) = (/(i, i=1,num_set(k1))/)
			dist_temp(1:num_set(k1)) = maxval(dist_s(1:num_set(k1))) &
				- dist_s(1:num_set(k1))
			call QsortC(dist_temp(1:num_set(k1)), ind_dist(1:num_set(k1)))
			do k2 = 1, size_group-ind_k
				k = ind_temp2((ind_dist(k2)))
				if (k .le. size_group) then
					group_new(ind_k+k2) = group_1(k)
				else
					group_new(ind_k+k2) = group_2(k-size_group)
				end if
				pareto_order(ind_k+k2) = k1
				dist_group(ind_k+k2) = dist_s(ind_dist(k2))
			end do
			ind_k = size_group
			exit
		end if
	end do
	end subroutine elitism


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

	integer, intent(in) ::	N_seq, n_table, size_group
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


	subroutine evaluate_idv(residue_peak, info_spectrum, connection_table, &
		n_table, N_seq, cnn_row, ngood, nbad, nedge, nused)
	!
	! Evaluate individuals subroutine
	! ===============================
	!
	! Input
	! -----
	! n_table
	! N_seq
	! cnn_row
	! residue_peak(N_seq, n_table)
	! info_spectrum(n_table)
	! connection_table(cnn_row, 6)
	!
	! Output
	! ------
	! ngood, nbad, nedge, nused
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

end module nsga2_module
