	program mc_hms_single

C CHANGES FOR OPTICS TESTING:
C 1. Remove target multiple scattering/energy loss

C+______________________________________________________________________________
!
! Monte-Carlo of HMS spectrometer using uniform illumination.
!   This version uses TRANSPORT right-handed coordinate system.
!
! Author: David Potterveld, March-1993
!
! Modification History:
!
!  11-Aug-1993	(D. Potterveld) Modified to use new transformation scheme in
!		which each transformation begins at the pivot.
!
!  19-AUG-1993  (D. Potterveld) Modified to use COSY INFINITY transformations.
C-______________________________________________________________________________

	implicit none

	include 'hms_simc/struct_hms.inc'
	include 'spectrometers.inc'
	include 'constants.inc'

C HBOOK/NTUPLE common block and parameters.
	integer*4	pawc_size
	parameter	(pawc_size = 80000)
	common		/pawc/ hbdata(pawc_size)
	integer*4	hbdata
	character*8	hut_nt_names(13)/
     >			'hsxfp', 'hsyfp', 'hsxpfp', 'hsypfp',
     >			'hsytari', 'hsdeltai', 'hsyptari', 'hsxptari',
     >			'hsytar', 'hsdelta', 'hsyptar', 'hsxptar','fry'/
	real*4		hut(13)

C Local declarations.
	integer*4	i,
     >			chanin	/1/,
     >			chanout	/2/,
     >			n_trials,trial,
     >			tmp_int

	logical*4	iss

C Event limits, topdrawer limits, physics quantities
	real*8 gen_lim(8)			!M.C. phase space limits.

	real*8 gen_lim_up(3)
	real*8 gen_lim_down(3)

	real*8 cut_dpp,cut_dth,cut_dph,cut_z	!cuts on reconstructed quantities
	real*8 xoff,yoff,zoff                   !Beam offsets (z target offset)
	real*8 spec_xoff,spec_yoff,spec_zoff    !Spectrometer offsets
	real*8 spec_xpoff, spec_ypoff           !Spectrometer angle offsets
	real*8 cos_ts,sin_ts			!cos and sin of spectrometer angle
	real*8 th_ev,cos_ev,sin_ev		!cos and sin of event angle
	real*8 x,y,z,dpp,dxdz,dydz,t1,t2,t3,t4	!temporaries
	real*8 musc_targ_len			!target length for multiple scattering
	real*8 m2				!particle mass squared.
	real*8 rad_len_cm			!conversion r.l. to cm for target
	real*8 pathlen				!path length through spectrometer.
C DJG Variables used for tuna can calcs
	real*8 t,atmp,btmp,ctmp,z_can
	real*8 side_path,costmp,th_can,s_Al
	logical*4 ok_spec			!indicates whether event makes it in MC
c	integer*4 hit_calo                      !flag for hitting the calorimeter

C Initial and reconstructed track quantities.
	real*8 dpp_init,dth_init,dph_init,xtar_init,ytar_init,ztar_init
	real*8 dpp_recon,dth_recon,dph_recon,ytar_recon,ztar_recon
	real*8 x_fp,y_fp,dx_fp,dy_fp		!at focal plane
	real*8 fry,fr1,fr2
	real*8 p_spec,th_spec			!spectrometer setting
	real*8 resmult

C Control flags (from input file)
	integer*4 p_flag			!particle identification
	logical*4 ms_flag
	logical*4 wcs_flag

C Hardwired control flags.
	logical*4 hut_ntuple	/.true./
	logical*4 decay_flag	/.false./

	real*8	dpp_var(2),dth_var(2),dph_var(2),ztg_var(2)
	real*8	stime,etime

	character*132	str_line

C Function definitions.

	integer*4	last_char
	logical*4	rd_int,rd_real
	real*8          grnd,gauss1

        character*80 rawname, filename
	real*4  secnds

        integer iquest
        common/quest/iquest(100)

	save		!Remember it all!

C ================================ Executable Code =============================

C Initialize
C xiaochao:
C using SIMC unstructured version
C
	hSTOP_trials	= 0
	hSTOP_slit_hor	= 0
	hSTOP_slit_vert	= 0
	hSTOP_slit_oct	= 0
	hSTOP_Q1_in	= 0
	hSTOP_Q1_mid	= 0
	hSTOP_Q1_out	= 0
	hSTOP_Q2_in	= 0
	hSTOP_Q2_mid	= 0
	hSTOP_Q2_out	= 0
	hSTOP_Q3_in	= 0
	hSTOP_Q3_mid	= 0
	hSTOP_Q3_out	= 0
	hSTOP_D1_in	= 0
	hSTOP_D1_out	= 0
	hSTOP_hut	= 0
	hSTOP_dc1	= 0
	hSTOP_dc2	= 0
	hSTOP_scin	= 0
	hSTOP_cal	= 0
	hSTOP_successes	= 0

C Open setup file.

	write(*,*)'Enter input filename (assumed to be in infiles dir)'
	read(*,1968) rawname
 1968	format(a)
	print *,'open infiles'
	filename = 'infiles/'//rawname(1:last_char(rawname))//'.inp'
	write(6,*) filename,'--- opened'
	open(unit=chanin,status='old',name=filename,readonly)

C Initialize HBOOK/NTUPLE if used.
	if (hut_ntuple) then
	  call hlimit(pawc_size)
	  filename = 'worksim/'//rawname(1:last_char(rawname))//'.rzdat'
!	  call hropen(30,'HUT',filename,'N',1024,i)
	  iquest(10) = 256000
	  iquest(10) = 510000
! see for example
!   http://wwwasd.web.cern.ch/wwwasd/cgi-bin/listpawfaqs.pl/7
! the file size is limited to ~260M no matter how I change iquest !
	  call hropen(30,'HUT',filename,'NQ',4096,i) !CERNLIB
 
	  if (i.ne.0) then
	    type *,'HROPEN error: istat = ',i
	    stop
	  endif
	  call hbookn(1,'HUT NTUPLE',13,'HUT',10000,hut_nt_names)
	endif	   

C Open Output file.
	print *,'open outfile'
	filename = 'outfiles/'//rawname(1:last_char(rawname))//'.hist'
	write(6,*) filename,'--- opened'
	open (unit=chanout,status='unknown',name=filename)

C Read in real*8's from setup file

	str_line = '!'

C Strip off header

	do while (str_line(1:1).eq.'!')
	  type *,str_line(1:last_char(str_line))
	  read (chanin,1001) str_line
	enddo

! Read data lines.

! N_TRIALS:
c	read (chanin,1001) str_line
	type *,str_line(1:last_char(str_line))
	iss = rd_int(str_line,n_trials)
	if (.not.iss) stop 'ERROR (ntrials) in setup!'

! Spectrometer momentum:
	read (chanin,1001) str_line
	type *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,p_spec)
	if (.not.iss) stop 'ERROR (Spec momentum) in setup!'

! Spectrometer angle:
	read (chanin,1001) str_line
	type *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,th_spec)
	if (.not.iss) stop 'ERROR (Spec theta) in setup!'
	th_spec = abs(th_spec) / degrad
	cos_ts = cos(th_spec)
	sin_ts = sin(th_spec)

! M.C. limits (half width's for dp,th,ph, full width's for x,y,z)
	do i=1,3
	  read (chanin,1001) str_line
	  type *,str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	  gen_lim_down(i) = gen_lim(i)
	  read (chanin,1001) str_line
	  type *,str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	  gen_lim_up(i) = gen_lim(i)
	enddo
!	do i=1,3
!	   write(*,*)'gen_lim_up/down = ',gen_lim_up(i),' ',gen_lim_down(i)
!	enddo

	do i = 4,6
	  read (chanin,1001) str_line
	  type *,str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	enddo

! Raster size

	do i=7,8
	   read (chanin,1001) str_line
	   type *,str_line(1:last_char(str_line))
	   iss = rd_real(str_line,gen_lim(i))
	   if (.not.iss) stop 'ERROR (Fast Raster) in setup'
	enddo

! Cuts on reconstructed quantities
	read (chanin,1001) str_line
	type *,str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dpp)) stop 'ERROR (CUT_DPP) in setup!'

	read (chanin,1001) str_line
	type *,str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dth)) stop 'ERROR (CUT_DTH) in setup!'

	read (chanin,1001) str_line
	type *,str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dph)) stop 'ERROR (CUT_DPH) in setup!'

	read (chanin,1001) str_line
	type *,str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_z)) stop 'ERROR (CUT_Z) in setup!'

! Read in radiation length of target material in cm
	read (chanin,1001) str_line
	type *,str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,rad_len_cm)) stop 'ERROR (RAD_LEN_CM) in setup!'

! Beam and target offsets
	read (chanin, 1001) str_line
	type *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,xoff)
	if(.not.iss) stop 'ERROR (xoff) in setup!'

	read (chanin, 1001) str_line
	type *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,yoff)
	if(.not.iss) stop 'ERROR (yoff) in setup!'

	read (chanin, 1001) str_line
	type *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,zoff)
	if(.not.iss) stop 'ERROR (zoff) in setup!'

! Spectrometer offsets
	read (chanin, 1001) str_line
	type *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_xoff)
	if(.not.iss) stop 'ERROR (spect. xoff) in setup!'

	read (chanin, 1001) str_line
	type *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_yoff)
	if(.not.iss) stop 'ERROR (spect. yoff) in setup!'

	read (chanin, 1001) str_line
	type *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_zoff)
	if(.not.iss) stop 'ERROR (spect. zoff) in setup!'

	read (chanin, 1001) str_line
	type *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_xpoff)
	if(.not.iss) stop 'ERROR (spect. xpoff) in setup!'

	read (chanin, 1001) str_line
	type *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_ypoff)
	if(.not.iss) stop 'ERROR (spect. ypoff) in setup!'

! read in flag for particle type.
	read (chanin,1001) str_line
	type *,str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,p_flag)) stop 'ERROR: p_flag in setup file!'

! Read in flag for aerogel usage in the HMS.
c	read (chanin,1001) str_line
c	type *,str_line(1:last_char(str_line))
c	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: use_aer in setup file!'
c	if (tmp_int.eq.1) use_aer = .true.

! Read in flag for multiple scattering.
	read (chanin,1001) str_line
	type *,str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: ms_flag in setup file!'
	if (tmp_int.eq.1) ms_flag = .true.

! Read in flag for wire chamber smearing.
	read (chanin,1001) str_line
	type *,str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: wcs_flag in setup file!'
	if (tmp_int.eq.1) wcs_flag = .true.

! Read in flag for dumping ALL events into the HUT NTUPLE (note that
! the ...recon quantities will be ill defined but the FAIL_ID could be
! used to tell...
!	read (chanin,1001) str_line
!	type *,str_line(1:last_char(str_line))
!	if (.not.rd_int(str_line,tmp_int)) stop
!	1    'ERROR:dump_all_in_ntuple in setup file!'
!	if (tmp_int.eq.1) dump_all_in_ntuple = .true.

C Set particle masses.
	m2 = me2			!default to electron
	if(p_flag.eq.0) then
	  m2 = me2
	else if(p_flag.eq.1) then
	  m2 = mp2
	else if(p_flag.eq.2) then
	  m2 = md2
	else if(p_flag.eq.3) then
	  m2 = mpi2
	else if(p_flag.eq.4) then
	  m2 = mk2
	endif

C------------------------------------------------------------------------------C
C                           Top of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C

	stime = secnds(0.0)
c	print *,'Enter total number of trials'
c	read *,n_trials
	do trial = 1,n_trials

	  if(mod(trial,5000).eq.0) write(*,*)'event #: ',trial,'       successes: ',hSTOP_successes

C Pick starting point within target. Z is picked uniformly, X and Y are
C chosen as truncated gaussians, using numbers picked above.
C Units are cm.

	  x = gauss1(3.0) * gen_lim(4) / 6.0			!beam width
	  y = gauss1(3.0) * gen_lim(5) / 6.0			!beam height
	  z = (grnd() - 0.5) * gen_lim(6)                       !along target

C DJG Assume flat raster
	  fr1 = (grnd() - 0.5) * gen_lim(7)   !raster x
	  fr2 = (grnd() - 0.5) * gen_lim(8)   !raster y

	  fry = -fr2  !+y = up, but fry needs to be positive when pointing down

	  x = x + fr1
	  y = y + fr2

	  x = x + xoff
	  y = y + yoff
	  z = z + zoff

C Pick scattering angles and DPP from independent, uniform distributions.
C dxdz and dydz in HMS TRANSPORT coordinates.

	  dpp  = grnd()*(gen_lim_up(1)-gen_lim_down(1))
     &             + gen_lim_down(1)
	  dydz = grnd()*(gen_lim_up(2)-gen_lim_down(2))
     &          /1000.   + gen_lim_down(2)/1000.
	  dxdz = grnd()*(gen_lim_up(3)-gen_lim_down(3))
     &          /1000.   + gen_lim_down(3)/1000.

C Transform from target to HMS (TRANSPORT) coordinates.
C Note that this assumes that HMS is on the right-hand side of the beam
C line (looking downstream).
	  xs    = -y
	  ys    = x * cos_ts + z * sin_ts
	  zs    = z * cos_ts - x * sin_ts

C DJG Apply spectrometer offsets
C DJG If the spectrometer if too low (positive x offset) a particle
C DJG at "x=0" will appear in the spectrometer to be a little high
C DJG so we should subtract the offset
	  xs = xs - spec_xoff
	  ys = ys - spec_yoff
	  zs = zs - spec_zoff

C Version for spectrometer on the left-hand side:
!	  xs    = -y
!	  ys    = x * cos_ts - z * sin_ts
!	  zs    = z * cos_ts + x * sin_ts

	  dpps  = dpp
	  dxdzs = dxdz
	  dydzs = dydz

C DJG Apply spectrometer angle offsets
	  dxdzs = dxdzs - spec_xpoff/1000.0
	  dydzs = dydzs - spec_ypoff/1000.0


C Save init values for later.
	  xtar_init = xs
	  ytar_init = ys
	  ztar_init = z
	  dpp_init = dpp
	  dth_init = dydzs*1000.		!mr
	  dph_init = dxdzs*1000.		!mr

C Drift back to zs = 0, the plane through the target center
	  xs = xs - zs * dxdzs
	  ys = ys - zs * dydzs
	  zs = 0.0

	  cos_ev = (cos_ts+dydzs*sin_ts)/sqrt(1+dydzs**2+dxdzs**2)
	  th_ev = acos(cos_ev)
	  sin_ev = sin(th_ev)

C Calculate multiple scattering length of target

C Case 1 : extended target:  TUNA CAN ONLY!!!!!, BEER can is old news.
C   cryo LH2 (4.0 cm diamater = 2.0 radius)
C   Liquid + 2 x 0.13mm Al (X0=8.89cm) tuna can
C
	  if (abs(gen_lim(6)).gt.3.) then   

C DJG Here I'm just copying stuff from SIMC - hope it works
C JRA this is ugly.  Solve for z position where particle intersects can.  The
C JRA pathlength is then (z_intersect - z_scatter)/cos(theta)
C JRA Angle from center to z_intersect is acos(z_intersect/R).  Therefore the
C JRA angle between the particle and wall is pi/2 - (theta - theta_intersect)
	     t=tan(th_ev)**2
	     atmp=1+t
	     btmp=-2*z*t
	     ctmp=z**2*t-(gen_lim(6)/2.)**2
	     z_can=(-btmp+sqrt(btmp**2-4.*atmp*ctmp))/2./atmp
	     side_path = (z_can - z)/abs(cos_ev)
	     costmp=z_can/(gen_lim(6)/2.)
	     if (abs(costmp).le.1) then
		th_can=acos(z_can/(gen_lim(6)/2.))
	     else if (abs(costmp-1.).le.0.000001) then
		th_can=0.	!extreme_trip_thru_target can give z/R SLIGHTLY>1.0
	     else
		stop 'z_can > can radius in target.f !!!'
	     endif
	     s_Al =  0.0050*2.54/abs(sin(pi/2 - (th_ev - th_can)))
	     musc_targ_len = side_path/rad_len_cm + s_Al/8.89
	  else
C Case 2 solid target
	     musc_targ_len = abs(gen_lim(6)/2. - z)/rad_len_cm/cos_ev
	  endif

C Scattering before magnets:  Approximate all scattering as occuring AT TARGET.
C  16 mil Al scattering chamber window (X0=8.89cm)
C  15(?) cm air (X0=30420cm)
C spectrometer entrance window
C  15 mil Kevlar (X0=74.6cm)
C  5 mil Mylar  (X0=28.7cm)             Total of 0.60% rad. length.

	  musc_targ_len = musc_targ_len + .016*2.54/8.89 +
     >          15/30420 + .015*2.54/74.6 + .005*2.54/28.7

C Begin transporting particle.  

	  if (ms_flag) call musc(m2,p_spec*(1.+dpps/100.),musc_targ_len,dydzs,dxdzs)

	  call mc_hms(p_spec, th_spec, dpps, xs, ys, zs, dxdzs, dydzs,
     >		x_fp, dx_fp, y_fp, dy_fp, m2,
     >		ms_flag, wcs_flag, decay_flag, resmult, fry, ok_spec, pathlen)

! Option for dumping all events is not implemented properly.  Skip it for now
!
!	  if (ok_spec.or.dump_all_in_ntuple) then !Success, increment arrays

	  if (ok_spec) then !Success, increment arrays
	    dpp_recon = dpps
	    dth_recon = dydzs*1000.			!mr
	    dph_recon = dxdzs*1000.			!mr
	    ytar_recon = + ys

C Output NTUPLE entry.

	    if (hut_ntuple) then
	      hut(1) = x_fp
	      hut(2) = y_fp
	      hut(3) = dx_fp
	      hut(4) = dy_fp
	      hut(5) = ytar_init
	      hut(6) = dpp_init
	      hut(7) = dth_init/1000.
	      hut(8) = dph_init/1000.
	      hut(9) = ytar_recon
	      hut(10)= dpp_recon
	      hut(11)= dth_recon/1000.
	      hut(12)= dph_recon/1000.
	      hut(13) = fry
!	      hut(13)= hit_calo 
	      call hfn(1,hut)
	    endif

C Cut on reconstructed quantities.
	    if ((abs(dpp_recon).gt.cut_dpp) .or.
     >		(abs(dth_recon).gt.cut_dth) .or.
     >		(abs(dph_recon).gt.cut_dph) .or.
     >		(abs(ztar_recon).gt.cut_z)) then
	      goto 500		!quit if failed
	    endif

C Compute sums for calculating reconstruction variances.
	    dpp_var(1) = dpp_var(1) + (dpp_recon - dpp_init)
	    dth_var(1) = dth_var(1) + (dth_recon - dth_init)
	    dph_var(1) = dph_var(1) + (dph_recon - dph_init)
	    ztg_var(1) = ztg_var(1) + (ztar_recon - ztar_init)

	    dpp_var(2) = dpp_var(2) + (dpp_recon - dpp_init)**2
	    dth_var(2) = dth_var(2) + (dth_recon - dth_init)**2
	    dph_var(2) = dph_var(2) + (dph_recon - dph_init)**2
	    ztg_var(2) = ztg_var(2) + (ztar_recon - ztar_init)**2
	  endif			!Incremented the arrays

C We are done with this event, whether GOOD or BAD.
C Loop for remainder of trials.

500	  continue
	enddo				!End of M.C. loop

C------------------------------------------------------------------------------C
C                           End of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C

	etime = secnds(stime)
	type *,'Elapsed time = ',etime,' seconds'
	type *,' '

C Close NTUPLE file.

	call hrout(1,i,' ')
	call hrend('HUT')

	write (chanout,1002)
	write (chanout,1003) p_spec,th_spec*degrad
        write (chanout,1004) (gen_lim(i),i=1,6)

	write (chanout,1005) n_trials

C Indicate where particles are lost in spectrometer.

	write (chanout,1015)
     >	hSTOP_slit_hor,hSTOP_slit_vert,hSTOP_slit_oct,
     >	hSTOP_Q1_in,hSTOP_Q1_mid,hSTOP_Q1_out,
     >	hSTOP_Q2_in,hSTOP_Q2_mid,hSTOP_Q2_out,
     >	hSTOP_Q3_in,hSTOP_Q3_mid,hSTOP_Q3_out,
     >	hSTOP_D1_in,hSTOP_D1_out

	write (chanout,1006)
     >	hSTOP_trials,hSTOP_hut,hSTOP_dc1,hSTOP_dc2,hSTOP_scin,hSTOP_cal,
     >  hSTOP_successes,hSTOP_successes

C Compute reconstruction resolutions.

	if (hSTOP_successes.eq.0) hSTOP_successes=1
	t1 = sqrt(max(0.,dpp_var(2)/hSTOP_successes - (dpp_var(1)/hSTOP_successes)**2))
	t2 = sqrt(max(0.,dth_var(2)/hSTOP_successes - (dth_var(1)/hSTOP_successes)**2))
	t3 = sqrt(max(0.,dph_var(2)/hSTOP_successes - (dph_var(1)/hSTOP_successes)**2))
	t4 = sqrt(max(0.,ztg_var(2)/hSTOP_successes - (ztg_var(1)/hSTOP_successes)**2))

	write (chanout,1011) dpp_var(1)/hSTOP_successes,t1,dth_var(1)/hSTOP_successes,
     >		t2,dph_var(1)/hSTOP_successes,t3,ztg_var(1)/hSTOP_successes,t4

	write(6,*) hSTOP_trials,' Trials',hSTOP_successes,' Successes'
	write (6,1011) dpp_var(1)/hSTOP_successes,t1,dth_var(1)/hSTOP_successes,
     >		t2,dph_var(1)/hSTOP_successes,t3,ztg_var(1)/hSTOP_successes,t4

C ALL done!

	stop ' '

C =============================== Format Statements ============================

1001	format(a)
1002	format('!',/,'! Uniform illumination Monte-Carlo results')
1003	format('!',/'! Spectrometer setting:',/,'!',/,
     >	g18.8,' =  P  spect (MeV)',/,
     >	g18.8,' =  TH spect (deg)')

1004	format('!',/'! Monte-Carlo limits:',/,'!',/,
     >	g18.8,' =  GEN_LIM(1) - DP/P                    (half width,% )',/,
     >	g18.8,' =  GEN_LIM(2) - Theta                   (half width,mr)',/,
     >	g18.8,' =  GEN_LIM(3) - Phi                     (half width,mr)',/,
     >	g18.8,' =  GEN_LIM(4) - HORIZ (full width of 3 sigma cutoff,cm)',/,
     >	g18.8,' =  GEN_LIM(5) - VERT  (full width of 3 sigma cutoff,cm)',/,
     >	g18.8,' =  GEN_LIM(6) - Z                       (Full width,cm)')

!inp     >	,/,
!inp     >	g18.8,' =  Hor. 1/2 gap size (cm)',/,
!inp     >	g18.8,' =  Vert. 1/2 gap size (cm)')

1005	format('!',/,'! Summary:',/,'!',/,
     >	i,' Monte-Carlo trials:')

1006	format(i,' Initial Trials',/
     >	i,' Trials made it to the hut',/
     >	i,' Trial cut in dc1',/
     >	i,' Trial cut in dc2',/
     >	i,' Trial cut in scin',/
     >	i,' Trial cut in cal',/
     >	i,' Trials made it thru the detectors and were reconstructed',/
     >	i,' Trials passed all cuts and were histogrammed.',/
     >	)

1008	format(8i)
1009	format(1x,i4,g,i)
1010	format(a,i)
1011	format(
     >	'DPP ave error, resolution = ',2g18.8,' %',/,
     >	'DTH ave error, resolution = ',2g18.8,' mr',/,
     >	'DPH ave error, resolution = ',2g18.8,' mr',/,
     >	'ZTG ave error, resolution = ',2g18.8,' cm')

1012	format(1x,16i4)

1015	format(/,
     >	i,' stopped in the FIXED SLIT HOR',/
     >	i,' stopped in the FIXED SLIT VERT',/
     >	i,' stopped in the FIXED SLIT OCTAGON',/
     >	i,' stopped in Q1 ENTRANCE',/
     >	i,' stopped in Q1 MIDPLANE',/
     >	i,' stopped in Q1 EXIT',/
     >	i,' stopped in Q2 ENTRANCE',/
     >	i,' stopped in Q2 MIDPLANE',/
     >	i,' stopped in Q2 EXIT',/
     >	i,' stopped in Q3 ENTRANCE',/
     >	i,' stopped in Q3 MIDPLANE',/
     >	i,' stopped in Q3 EXIT',/
     >	i,' stopped in D1 ENTRANCE',/
     >	i,' stopped in D1 EXIT',/
     >	)

1100	format('!',79('-'),/,'! ',a,/,'!')
1200	format(/,'! ',a,' Coefficients',/,/,
     >	(5(g,','))
     >	)
1300	format(/,'! ',a,' Coefficient uncertainties',/,/,
     >	(5(g,','))
     >	)

	end
