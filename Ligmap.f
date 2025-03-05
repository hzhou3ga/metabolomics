! ncludes maximal sequence identity of the native ligand
! pv<0.003 ni=6
! pid<=0.6
!outputs to pocket_MET_production to pdb_heme_docked_homola
! fixed misassignment of nh; now correctly outputs pdb files
!
! really includes homology cutoff
! includes homology cutoff
! fixed ps and pv cutoff
! now sequence identity cutoff
!4/3/2024/ pid>0.3 excluded   now includes full pocket
! for hem template library
! adds touching constraint  full pocket 
! ni>4 and np>1
! like  pocket_rings_to_pdb_sel_pos_fixed2_401_all_b_dock_b2.f but selects particular user defined pockets
! also outputs just the coordinates
!io>2 is rejected
!analyzes trajectory
! generalized to arbitrary file
! has backbone excluded volumen check ni>=4 allowed, io<nh allowed
! only docks adjacent ligands to known ligand binding sites
! like pocket_rings_to_pdb_sel_pos_fixed2_401_all_b.f but outputs the docked ligands
! pv>0.01x
! ni>=4 allowed
! uses positions of the ligsnds
! includes pocket information
! includes parent rings
! based on pocket_matchbottompid3_pdb_rings_pocket_3.f
! much match contacting ligands exactly
!like pocket_matchbottompid3_pdb_rings_adjacent2_2.f but now considers ligands in the pocket
! ni>=5
! updated to conhet3
! also no directly focuses on the individual contacting rings
! like pocket_matchbottompid3.f but uses pdb files numberging
! 3/27/20/
! how often can you predict true pockets?
! 	pocket_match.f
!used to construct a set of conditional probabilities of interactions.

	PARAMETER(NMAXP=1200)
	PARAMETER(NMAX=1200)

	PARAMETER(NCASES=120000)
	PARAMETER(NCASES_TAR=120000)
!
	CHARACTER*1 A1,CHAIN
	CHARACTER*3 AA(0:19),LIG3,LIG_LIB(500,200000)

	CHARACTER*3 MET(1000),LIG_TYPE,LTYPE(1000)
	CHARACTER*4 HET_PAR(100,NCASES),het_par_tem(100,1000,1000)

	CHARACTER*4 NAME4
	CHARACTER*5 NAME5,NAMEPDB(NCASES_TAR),NAME5TEM,NAMELIG(200000),NAMESEQ(NCASES)
	CHARACTER*5 ARGV,name5h,name52
	CHARACTER*6 ATOM
	!CHARACTER*7 ECNUMBER(NCASES)
	CHARACTER*12 NAME12  
	CHARACTER*12 NAMEPARENT(NCASES),name12max
	CHARACTER*255 PDBLOCATION,DATA6,ALIGNMENT,PDBLOCATION23!,TITLE(NCASES)

	LOGICAL*1 TLIG(1000),TASS(1000)!,TMET(1000)
	

	DIMENSION :: X(6)
	DIMENSION :: CA_PAR(3,NMAX,NCASES),CA(3,NMAX),XHET(3,100,NCASES)
	DIMENSION :: XROT(3,100),RMSD(500,1000)
	DIMENSION :: PID(NCASES),PMAX(NCASES)
	REAL*4 XHET_PAR(3,100,500,1000)

	INTEGER :: MM1A(NMAXP),MM2A(NMAXP)
	INTEGER :: JCODE(NMAXP),JLIG(1000)
	INTEGER :: ICODE(NMAX,NCASES)
	INTEGER :: NLIG(1000),IMAP_MET(1000)
	INTEGER :: NHET(500,1000)
	INTEGER :: NH_PAR(NCASES),NTOT(1000),NGLIG(1000)
	INTEGER :: NRES(NCASES),LIG_TOT(200000),MET_TYPE(200000)

	REAL*8 PV
	DOUBLE PRECISION R_1(3,NMAX),R_2(3,NMAX),DRMS,W(NMAX)
	DOUBLE PRECISION U(3,3),T(3)

	DATA AA/ 'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO','MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG','HIS','PHE','TYR','TRP'/
!	*********************************************************************************************


	data6='../pdb_ligands/'
	alignment='alignmentbottom_MET_enz/'!alignmentbottom_poc_MET_v2_pocket_enzyme2/'
	pdblocation='/local/images/pdb-202109a/str/'
	pdblocation23='../../2023/pdbdirectory/pdb_library/'


	do i=1,nmax
	w(i)=1.d0
	end do
	
	open(unit=1,action='read',file='LIST.metabolites')
	rewind(1)
	read(1,*)NMET 
	do im=1,NMET
		read(1,*)MET(im)
	end do
	close(1)
!	-----------------------
!      Read in template structures and their associated ligands:

	open(unit=20,action='read',file='LIST.MET_2021_pockets_good') ! list of metabolites being screened
	rewind(20)
	read(20,*)npar
	do ik=1,npar
		read(20,*)nameparent(ik)!,id1,id2,id3
		name12=nameparent(ik)
		LIG_TYPE=name12(7:9)
!	write(6,*)ik,name12
		Do IM=1,NMET
			if(LIG_TYPE==MET(IM))then
				MET_TYPE(IK)=IM
				exit
			end if
		end do
 	    open(unit=54, action='read',file='../pdb_ligands/'//name12(2:3)//'/'//trim(adjustl(name12))//'.hetcord') ! 3d coordinates of metabolites
		rewind(54)
        read(54,*)nh_par(ik)
	    do id=1,nh_par(ik)
                  read(54,*)idr,(xhet(j,id,ik),j=1,3),het_par(id,ik)
		end do
		close(54)
		open(unit=45,action='read',file='../pdb_ligands/'//name12(2:3)//'/'//name12(1:5)//'.cacbfiles') ! ca and cb coordinatesx and aminoacid code of residues
			rewind(45)
			read(45,*)nres(ik)
			do i=1,Nres(ik)
				read(45,*)id,(ca_par(j,i,ik),j=1,3),(X(J),j=1,3),icode(i,ik)
			end do
			close(45)
!	go to 157
!156	name5=name12(1:5)
!	write(66,66)name5
!	write(6,66)name5
!66	format(a5)
157	end do
	close(20)
	open(unit=43,action='read',file='LIST.pockets_good_enz')
!	open(unit=43,action='read',file='LIST.
	rewind(43)
	read(43,*)NPDB
	do IPDB=1,NPDB
		read(43,*)NAMEPDB(IPDB)
	end do
	close(43)
! what are the ligands that bind to the given native structure?
!	open(unit=53,action='read',file='LIST.enzymes_ligands.summary') ! these are the native ligands if you want to check them
	open(unit=53,action='read',file='LIST.metabolites_pockets')
	rewind(53)
	read(53,*)MLIG
	do ilig=1,MLIG
		read(53,*)namelig(ilig),n,(lig_lib(j,ilig),j=1,n)
		lig_tot(ilig)=n
	end do
	close(53)
        I=1
	 call get_command_argument(i, argv)
	    read(argv,"(I)")input
	do 2999 IK=(INPUT-1)*50+1,min(input*50,NPDB)    
!	do 2999 IK=1,NPDB
		name5=namepdb(ik)
		write(6,*)IK,name5
		DO IM=1,NMET
			TLIG(IM)=.false. ! is the ligand found in template matching?
			TASS(IM)=.FALSE. ! is this ia native ligand
		END DO
	
		do ilig=1,MLIG	
			if(name5==namelig(ilig))then
				do j=1,lig_tot(ilig)
!				write(6,*)j,lig_lib(j,ilig)
						DO IM=1,NMET
							if(lig_lib(j,ilig)== MET(IM))then
								TASS(IM)=.true.
!								write(6,*)name5,MET(IM),' name5 met'
								EXIT
							end if
						end do
				end do
			end if
		end do

		! list of unique native ligands that are bound
		KMET=0
			DO IM=1,NMET	
				IF(TASS(IM))then
					KMET=KMET+1
					pmax(kmet)=0.
					IMAP_MET(KMET)=IM
				end if
			end do
!50		CONTINUE

		write(6,222)ik,namepdb(ik),kmet,(met(imap_met(j)),j=1,kmet)
222	format(i5,1x,a5,1x,i4,1x,200(a3,1x))
		open(unit=45,action='read',file='../pdb_ligands/'//name5(2:3)//'/'//name5//'.cacbfiles')
		rewind(45)
		read(45,*,END=2999,ERR=2999)MRES
		do i=1,MRES
			read(45,*)jd,(ca(j,i),j=1,3),(x(j),j=1,3),jcode(i)
		end do
		close(45)
		open(unit=55,action='read',file='../pdb2024/'//name5(2:3)//'/'//name5//'.homol_human')
		rewind(55)
		read(55,*,END=2999,ERR=2999)nhom
		do ihom=1,nhom
			read(55,*)name5h,pid(ihom)
			nameseq(ihom)=name5h
		end do
		close(55)
551		CONTINUE

		NG=0
		DO IM=1,NMET
            NTOT(IM)=0
		    NGLIG(IM)=0
			NLIG(IM)=0 ! now generalized for multiple ligands which ones bind to the target structure?
		end do
	
		open(unit=591,action='read',file=trim(adjustl(alignment))//namePDB(IK)//'.aln1')
		rewind(591)
!		read(591,*,END=2999,ERR=2999)ntem
	read(591,*)ntem
		psel=0
		do 200 lk=1,ntem
     		read(591,592)name12,naln,nr2,nr1,ps,ncv,tm,pidtem,psim,pv
592     	format(a12,1x,3(i4,1x),1f7.3,1x,i4,1x,3(1f7.3,1x),1pd9.3)
      		read(591,593)(mm1a(i),mm2a(i),i=1,naln)
593			format(20(i5,1x,i5,1x))
		if(name52==name5)go to 200
!            if(ps.lt.0.35)go to 200 ! abort this will not satisfy pv<0.05
            if(pv>0.003D0)go to 200
			JKP=0
			do ipar=1,NPAR
				if(trim(adjustl(name12))== trim(adjustl(nameparent(ipar))))then
					JKP=ipar
					ILIG_TYPE=MET_TYPE(IPAR)! what type of metabolite is the template ligand?
					exit
				end if
			end do
			IF(JKP==0) go to 200
			name5tem=name12(1:5)
			do ih=1,nhom
				if(name5tem==nameseq(ih))then	
					if(pid(ih)>0.3)then
					go to 200
					end if
					psel=pid(ih)
				end if
			end do
552			continue

			ni=0
			do ia=1,naln
	   ! remember order is flipped
				i=mm2a(ia)
				ires=jcode(i)
				j=mm1a(ia)
				jres=icode(j,JKP)
				if(ires==jres)ni=ni+1
	    	end do
		if(ni <6) go to 200 ! must match at least 4 residues and be part of the pocket
	
		ll=naln
! rotate template to target pocket
			do ia=1,naln
				i=mm2a(ia)
				j=mm1a(ia)
				do jj=1,3
!	            r_1(jj,ia)=ca(jj,i)
!		      r_2(jj,ia)=ca_par(jj,j,JKP)
	            	r_2(jj,ia)=ca(jj,i)
		      		r_1(jj,ia)=ca_par(jj,j,JKP)
				end do
			enddo
  	 		call u3b(w,r_1,r_2,ll,1,drms,u,t,ier) !u rotate r_1 to r_2
	     	armsd=dsqrt(drms/ll)
			if(armsd>3.5)go to 200
!

      		nh=nh_par(JKP)
       		do ih=1,nh
       			do j=1,3
                    x(j)=xhet(j,ih,JKP)
            	end do
            	xrot(1,ih)=t(1)+u(1,1)*x(1)+u(1,2)*x(2)+u(1,3)*x(3)
          		xrot(2,ih)=t(2)+u(2,1)*x(1)+u(2,2)*x(2)+u(2,3)*x(3)
            	xrot(3,ih)=t(3)+u(3,1)*x(1)+u(3,2)*x(2)+u(3,3)*x(3)
       		end do

			io=0
			itouch=0
			do i=1,mres
			do ih=1,nh
				dd=0.
				do j=1,3
					dd=dd+(xrot(j,ih)-ca(j,i))**2 ! refers to target protein
				enddo
!				if(dd.lt.25)then	! cannot have too many excluded volume overlaps
				if(dd.lt.9)then	! cannot have too many excluded volume overlaps
					io=io+1 
				end if
				if(dd.lt.30)then ! ligand must touch the target proten
					itouch=itouch+1
				end if		
			end do
			end do

			if(itouch <3) go to 200! rotated ligand must touch the protein
			if(io.gt.2*nh)then !cannot have too many ligand-protein overlaps
				go to 200
			end if
!	     write(6,*)name12,armsd

			NG=NG+1
!	write(6,*)name12,ni,MET(ILIG_TYPE)
			IF(NG>500)GO TO 200
			TLIG(ILIG_TYPE)=.true.
!			RMSD(NG)=armsd
			NLIG(ILIG_TYPE)=NLIG(ILIG_TYPE)+1
			NN=NLIG(ILIG_TYPE)
			RMSD(NN,ILIG_TYPE)=Armsd
			NHET(NN,ILIG_type)=NH
			DO Imm=1,KMET
				if(imap_met(imm)==ILIG_TYPE)then
					if(pmax(imm)<psel)then
						pmax(imm)=psel
					end if
				end if
			end do
!				----------------:
!ligand is now accepted
! now rotate into the crystal structure' reference system
			!NHET(NG,ILIG_TYPE)=NH this ignores ligand type
			do ih=1,nh
				het_par_tem(ih,NN,ILIG_TYPE)=het_par(ih,JKP)
				do j=1,3
					xhet_par(j,ih,NN,ILIG_TYPE)=xrot(j,ih)
				end do
			end do
200		continue
		close(591)

		DO IM=1,NMET
			IF(TASS(IM))then
				NTOT(IM)=NTOT(IM)+1
				IF(TLIG(IM))NGLIG(IM)=NGLIG(IM)+1
!			write(6,*)MET(IM),NTOT(IM),NGLIG(IM)
			end if
		end do
        open(unit=19,action='write',file='results_fixed_0.003_6_v2_enz/native/'//namepdb(IK)//'.native')
! list of matched native metabolitesx
		write(19,19)namepdb(ik),kmet,(met(imap_met(j)),nlig(imap_met(j)),j=1,kmet)
19		format(a5,1x,i4,1x,100(a3,1x,i3))
		ntt=0
		do j=1,kmet
		if(nlig(imap_met(j))>0)ntt=ntt+1
		end do
		write(19,*)ntt
		do j=1,kmet
		if(nlig(imap_met(j))>0)then
		write(19,191)met(imap_met(j)),nlig(imap_met(j)),pmax(j)
		end if
		end do
191		format(a3,1x,i4,1x,1f7.3)
		
        close(19)

	open(unit=59,Action='write',file='results_fixed_0.003_6_v2_enz/summary/'//namepdb(IK)//'.summary')
! summary of predicted metabolites
	NCS=0
	DO IM=1,NMET
		IF(TLIG(IM))then
			NCS=NCS+1
			JLIG(NCS)=NLIG(IM)
			LTYPE(NCS)=MET(IM)
		END IF
	END DO
	write(59,9)NCS,(LTYPE(J),JLIG(J),J=1,NCS) !total number of ligands\
	!write(9,9)namepdb(ik),NCS,(LTYPE(J),JLIG(J),J=1,NCS)
	!write(6,9)namepdb(ik),NCS,(LTYPE(J),JLIG(J),J=1,NCS)
9	format(i4,1x,300(a3,1x,i3,1x))
	IF(NCS==0)GO TO 2999! don't output the pdb file
	DO IM=1,NMET
		IF(TLIG(IM))then
			write(59,91)nlig(im),MET(IM)
91	format(i6,1x,a3)
		END IF
	END DO
            close(59)
	A1=namepdb(ik)(5:5)
	name4=namepdb(ik)(1:4)
	open(UNIT=10,action='read',file=trim(pdblocation23)//name4//'.pdb')
	rewind(10)
	read(10,*,ERR=192,END=192)
	rewind(10)
	go to 193
192	open(UNIT=10,action='read',file=trim(pdblocation)//name4//'.pdb')	
	rewind(10)
193	CONTINUE
	open(unit=11,action='write',file='results_fixed_0.003_6_v2_enz/pdb/'//trim(adjustl(namepdb(ik)))//'.pdb')
	rewind(11)

! coordinates are rotated in the coordinate system of the parent template
	!write(6,*)'entering pdbside_target'
		        call pdbside_target(A1)
	!write(6,*)'exiting pdbside_target'
			ATOM='HETATM'
			ib=mres
		write(11,3)'TER_'
3	format(a3)
		write(11,*)NCS
	DO IM=1,NMET
		IF(TLIG(IM))then
			write(11,11)'REMARK',nlig(IM),MET(IM)
			lig3=MET(IM)
11	format(a6,1x,i4,1X,A3)
			CHAIN='Q'
			do ics=1,nlig(im)
				write(11,111)'REMARK',rmsd(ics,IM),NHET(ICS,IM)
111	format(a6,1x,1f7.3,1x,i4)
					ib=0
				do ih=1,nhet(ICS,IM)
					ib=ib+1
!					write(11,10)atom,ib,het_par_tem(ih,ics,im),lig3,chain,ih,(xhet_par(j,ih,ics,im),j=1,3)
					write(11,10)atom,ib,het_par_tem(ih,ics,im),lig3,chain,ih,(xhet_par(j,ih,ics,im),j=1,3)
10     				format(a6,i5,1x,a4,1x,a3,1x,a1,i4,1x,3x,3f8.3)
					!write(12,12)ih,(xhet_par(j,ih,ics,ipar),j=1,3),het_par(ih,1)
12 	format(i3,1x,3f8.3,1x,a4)
				end do
				write(11,112)'TER'
112	format(a3)
 			END DO
		END IF
	END DO
	close(11)
	close(12)
	close(10)

2999	continue
	END

      subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
      integer ip(9), ip2312(4), i, j, k, l, m1, m, ier, n, mode
      PARAMETER(nmax=1200)
      double precision w(nmax), x(3, nmax), y(3, nmax), u(3, 3), t(3), rms, sigma
      double precision r(3, 3), xc(3), yc(3), wc, a(3, 3), b(3, 3), e0, 
     &e(3), e1, e2, e3, d, spur, det, cof, h, g, cth, sth, sqrth, p, tol
     &, rr(6), rr1, rr2, rr3, rr4, rr5, rr6, ss(6), ss1, ss2, ss3, ss4, 
     &ss5, ss6, zero, one, two, three, sqrt3
      equivalence (rr(1), rr1), (rr(2), rr2), (rr(3), rr3), (rr(4), rr4)
     &, (rr(5), rr5), (rr(6), rr6), (ss(1), ss1), (ss(2), ss2), (ss(3), 
     &ss3), (ss(4), ss4), (ss(5), ss5), (ss(6), ss6), (e(1), e1), (e(2)
     &, e2), (e(3), e3)
      data sqrt3 / 1.73205080756888d+00 /
      data tol / 1.0d-2 /
      data zero / 0.0d+00 /
      data one / 1.0d+00 /
      data two / 2.0d+00 /
      data three / 3.0d+00 /
      data ip / 1, 2, 4, 2, 3, 5, 4, 5, 6 /
      data ip2312 / 2, 3, 1, 2 /
c 156 "rms.for"
      wc = zero
	!!write(6,*)'n',n
      rms = 0.0
      e0 = zero
      do 1 i = 1, 3
      xc(i) = zero
      yc(i) = zero
      t(i) = 0.0
      do 1 j = 1, 3
      d = zero
      if (i .eq. j) d = one
      u(i,j) = d
      a(i,j) = d
    1 r(i,j) = zero
      ier = -1
c**** DETERMINE CENTROIDS OF BOTH VECTOR SETS X AND Y
c 170 "rms.for"
      if (n .lt. 1) return 
c 172 "rms.for"
      ier = -2
      do 2 m = 1, n
      if (w(m) .lt. 0.0) return 
      wc = wc + w(m)
      do 2 i = 1, 3
      xc(i) = xc(i) + (w(m) * x(i,m))
    2 yc(i) = yc(i) + (w(m) * y(i,m))
      if (wc .le. zero) return 
      do 3 i = 1, 3
      xc(i) = xc(i) / wc
c**** DETERMINE CORRELATION MATRIX R BETWEEN VECTOR SETS Y AND X
c 182 "rms.for"
    3 yc(i) = yc(i) / wc
c 184 "rms.for"
      do 4 m = 1, n
      do 4 i = 1, 3
      e0 = e0 + (w(m) * (((x(i,m) - xc(i)) ** 2) + ((y(i,m) - yc(i)) ** 
     &2)))
c 187 "rms.for"
      d = w(m) * (y(i,m) - yc(i))
      do 4 j = 1, 3
c**** CALCULATE DETERMINANT OF R(I,J)
c 189 "rms.for"
    4 r(i,j) = r(i,j) + (d * (x(j,m) - xc(j)))
c 191 "rms.for"
      det = ((r(1,1) * ((r(2,2) * r(3,3)) - (r(2,3) * r(3,2)))) - (r(1,2
     &) * ((r(2,1) * r(3,3)) - (r(2,3) * r(3,1))))) + (r(1,3) * ((r(2,1)
     & * r(3,2)) - (r(2,2) * r(3,1))))
c**** FORM UPPER TRIANGLE OF TRANSPOSED(R)*R
c 194 "rms.for"
      sigma = det
c 196 "rms.for"
      m = 0
      do 5 j = 1, 3
      do 5 i = 1, j
      m = m + 1
c***************** EIGENVALUES *****************************************
c**** FORM CHARACTERISTIC CUBIC  X**3-3*SPUR*X**2+3*COF*X-DET=0
c 200 "rms.for"
    5 rr(m) = ((r(1,i) * r(1,j)) + (r(2,i) * r(2,j))) + (r(3,i) * r(3,j)
     &)
c 203 "rms.for"
      spur = ((rr1 + rr3) + rr6) / three
      cof = ((((((rr3 * rr6) - (rr5 * rr5)) + (rr1 * rr6)) - (rr4 * rr4)
     &) + (rr1 * rr3)) - (rr2 * rr2)) / three
c 205 "rms.for"
      det = det * det
      do 6 i = 1, 3
    6 e(i) = spur
c**** REDUCE CUBIC TO STANDARD FORM Y**3-3HY+2G=0 BY PUTTING X=Y+SPUR
c 208 "rms.for"
      if (spur .le. zero) goto 40
c 210 "rms.for"
      d = spur * spur
      h = d - cof
c**** SOLVE CUBIC. ROOTS ARE E1,E2,E3 IN DECREASING ORDER
c 212 "rms.for"
      g = (((spur * cof) - det) / two) - (spur * h)
c 214 "rms.for"
      if (h .le. zero) goto 8
      sqrth = dsqrt(h)
      d = ((h * h) * h) - (g * g)
      if (d .lt. zero) d = zero
      d = datan2(dsqrt(d),- g) / three
      cth = sqrth * dcos(d)
      sth = (sqrth * sqrt3) * dsin(d)
      e1 = (spur + cth) + cth
      e2 = (spur - cth) + sth
      e3 = (spur - cth) - sth
c.....HANDLE SPECIAL CASE OF 3 IDENTICAL ROOTS
c 224 "rms.for"
      if (mode) 10, 50, 10
c**************** EIGENVECTORS *****************************************
c 226 "rms.for"
    8 if (mode) 30, 50, 30
c 228 "rms.for"
   10 do 15 l = 1, 3, 2
      d = e(l)
      ss1 = ((d - rr3) * (d - rr6)) - (rr5 * rr5)
      ss2 = ((d - rr6) * rr2) + (rr4 * rr5)
      ss3 = ((d - rr1) * (d - rr6)) - (rr4 * rr4)
      ss4 = ((d - rr3) * rr4) + (rr2 * rr5)
      ss5 = ((d - rr1) * rr5) + (rr2 * rr4)
      ss6 = ((d - rr1) * (d - rr3)) - (rr2 * rr2)
      j = 1
      if (dabs(ss1) .ge. dabs(ss3)) goto 12
      j = 2
      if (dabs(ss3) .ge. dabs(ss6)) goto 13
   11 j = 3
      goto 13
   12 if (dabs(ss1) .lt. dabs(ss6)) goto 11
   13 d = zero
      j = 3 * (j - 1)
      do 14 i = 1, 3
      k = ip(i + j)
      a(i,l) = ss(k)
   14 d = d + (ss(k) * ss(k))
      if (d .gt. zero) d = one / dsqrt(d)
      do 15 i = 1, 3
   15 a(i,l) = a(i,l) * d
      d = ((a(1,1) * a(1,3)) + (a(2,1) * a(2,3))) + (a(3,1) * a(3,3))
      m1 = 3
      m = 1
      if ((e1 - e2) .gt. (e2 - e3)) goto 16
      m1 = 1
      m = 3
   16 p = zero
      do 17 i = 1, 3
      a(i,m1) = a(i,m1) - (d * a(i,m))
   17 p = p + (a(i,m1) ** 2)
      if (p .le. tol) goto 19
      p = one / dsqrt(p)
      do 18 i = 1, 3
   18 a(i,m1) = a(i,m1) * p
      goto 21
   19 p = one
      do 20 i = 1, 3
      if (p .lt. dabs(a(i,m))) goto 20
      p = dabs(a(i,m))
      j = i
   20 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((a(k,m) ** 2) + (a(l,m) ** 2))
      if (p .le. tol) goto 40
      a(j,m1) = zero
      a(k,m1) = - (a(l,m) / p)
      a(l,m1) = a(k,m) / p
   21 a(1,2) = (a(2,3) * a(3,1)) - (a(2,1) * a(3,3))
      a(2,2) = (a(3,3) * a(1,1)) - (a(3,1) * a(1,3))
c****************** ROTATION MATRIX ************************************
c 282 "rms.for"
      a(3,2) = (a(1,3) * a(2,1)) - (a(1,1) * a(2,3))
c 284 "rms.for"
   30 do 32 l = 1, 2
      d = zero
      do 31 i = 1, 3
      b(i,l) = ((r(i,1) * a(1,l)) + (r(i,2) * a(2,l))) + (r(i,3) * a(3,l
     &))
c 288 "rms.for"
   31 d = d + (b(i,l) ** 2)
      if (d .gt. zero) d = one / dsqrt(d)
      do 32 i = 1, 3
   32 b(i,l) = b(i,l) * d
      d = ((b(1,1) * b(1,2)) + (b(2,1) * b(2,2))) + (b(3,1) * b(3,2))
      p = zero
      do 33 i = 1, 3
      b(i,2) = b(i,2) - (d * b(i,1))
   33 p = p + (b(i,2) ** 2)
      if (p .le. tol) goto 35
      p = one / dsqrt(p)
      do 34 i = 1, 3
   34 b(i,2) = b(i,2) * p
      goto 37
   35 p = one
      do 36 i = 1, 3
      if (p .lt. dabs(b(i,1))) goto 36
      p = dabs(b(i,1))
      j = i
   36 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((b(k,1) ** 2) + (b(l,1) ** 2))
      if (p .le. tol) goto 40
      b(j,2) = zero
      b(k,2) = - (b(l,1) / p)
      b(l,2) = b(k,1) / p
   37 b(1,3) = (b(2,1) * b(3,2)) - (b(2,2) * b(3,1))
      b(2,3) = (b(3,1) * b(1,2)) - (b(3,2) * b(1,1))
      b(3,3) = (b(1,1) * b(2,2)) - (b(1,2) * b(2,1))
      do 39 i = 1, 3
      do 39 j = 1, 3
c****************** TRANSLATION VECTOR *********************************
c 320 "rms.for"
   39 u(i,j) = ((b(i,1) * a(j,1)) + (b(i,2) * a(j,2))) + (b(i,3) * a(j,3
     &))
   40 do 41 i = 1, 3
c********************** RMS ERROR **************************************
c 323 "rms.for"
   41 t(i) = ((yc(i) - (u(i,1) * xc(1))) - (u(i,2) * xc(2))) - (u(i,3)
     & * xc(3))
   50 do 51 i = 1, 3
      if (e(i) .lt. zero) e(i) = zero
   51 e(i) = dsqrt(e(i))
      ier = 0
      if (e2 .le. (e1 * 1.0d-05)) ier = -1
      d = e3
      if (sigma .ge. 0.0) goto 52
      d = - d
      if ((e2 - e3) .le. (e1 * 1.0d-05)) ier = -1
   52 d = (d + e2) + e1
      rms = (e0 - d) - d
      if (rms .lt. 0.0) rms = 0.0
      return 
c.....END U3B...........................................................
c----------------------------------------------------------
c                       THE END
c----------------------------------------------------------
c 338 "rms.for"
      end

	SUBROUTINE PDBSIDE_target(A1)
	PARAMETER (NMAX=1200)
	PARAMETER (NMAXP=1200)
	characTER*1 A1,CHAIN
	characTER*3 aa(0:19)
	CHARACTER*4 ITYPE
	CHARACTER*6 ATOM
	CharacTER*3 ARES
	!DIMENSION  SCM(3)
	DIMENSION X(3)
	!DOUBLE PRECISION U(3,3),T(3,3)
        data aa/'GLY','ALA','SER','CYS',
     &     'VAL','THR','ILE','PRO','MET',
     &	   'ASP','ASN','LEU','LYS','GLU',
     &	   'GLN','ARG','HIS','PHE','TYR',
     &	   'TRP'/
	!write(6,*)'entering pdb files'
9	CONTINUE
      READ(10,10,END=22,ERR=9)ATOM,ID,ITYPE,ARES,CHAIN,INUM,X(1),X(2),X(3)
	!write(6,10)ATOM,ID,ITYPE,ARES,CHAIN,INUM,X(1),X(2),X(3)
10	FORMAT(A6,I5,1X,A4,1x,a3,1X,a1,I4,1x,3X,3F8.3)
	IF(CHAIN /= A1)go to 9
	IF(ATOM(1:3) == 'TER')GO TO 20
	IF(ATOM(1:5) .eq.'MODEL')then
	go to 9
	ELSEIF(ATOM .eq. 'ENDMDL')then
	go to 22
	end if
	IF(ATOM(1:4).ne.'ATOM') go to 9
	write(11,10)ATOM,ID,ITYPE,ARES,CHAIN,INUM,X(1),X(2),X(3)
!	write(6,10)ATOM,ID,ITYPE,ARES,CHAIN,INUM_R,X(1),X(2),X(3)
	GO TO 9
20	continue
	write(11,21)'TER'
21	format(a3)
	!write(6,*)'exiting'
22	RETURN
	END


