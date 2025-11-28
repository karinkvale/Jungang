#if defined O_plastic
!=======================================================================
!     read and interpolate plastic distribution data
!=======================================================================

      subroutine getplastic

      implicit none

      character(120) :: fname, name, vname, new_file_name, text

      integer i, iou, j, n, ln, ib(10), ic(10)

      logical first_time, intrp, exists, inqvardef

      real data_time, wt3, wt1, yrl(3), iyr(3), p1e10

      real, allocatable :: time(:)

      save time, ln, yrl, first_time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "atm.h"
      include "cembm.h"
      include "levind.h"
      include "tmngr.h"

      real dmsk(imt,jmt), tmp, tmpij(imtm2,jmtm2)

      name = "O_plastic.nc"

      vname = "O_PLASTIC"

# if defined O_plastic_conc
! KK incoming units are particles per km^2 but I think this is misleading
! neuston nets sample the upper meter or so, but this isn't strict
! Kooi et al. (2018?) models plastic distributed throughout upper 100 m
! consider input particles per km^2*m^-1, average over k=1, convert to
! particles per cm^3
      p1e10 = 1.e-10
# else
! incoming units are particles per m^2 per second. These input fluxes are based on an assumption of 100% total historical plastic production (Geyer et al. 2017) is microplastic. Tons were converted to particles in the input file using a microplastic mass conversion from Eriksen et al (2014) and scaled against the F_co2dist.nc file.
! Now assume here 100% of the total plastic production enters the sea as mis-managed waste and 100% of this is microplastic
      p1e10 = 1.e-4 
# endif


      if (.not. allocated (time)) then
        fname = new_file_name (name)
        inquire (file=trim(fname), exist=exists)
        if (.not. exists) then
          print*, "Warning => ", trim(fname), " does not exist."
          ln = 3
          allocate ( time(ln) )
          time(:) = year0
          plastic(:,:,:) = 1.
          first_time = .false.
        else
          call openfile (fname, iou)
          call getdimlen ('TIME', iou, ln)
          allocate ( time(ln) )
          exists = inqvardef(trim(vname), iou)
          if (.not. exists) then
            print*, "==>  Warning: F_plastic data does not exist."
          endif
          ib(:) = 1
          ic(:) = ln
          call getvara ('TIME', iou, ln, ib, ic, time, c1, c0)
          text = 'years'
          call getatttext (iou, 'TIME', 'units', text)
          if (trim(text) .eq. "days since 1-1-1")
     &      time(:) = time(:)/yrlen - 1.
          if (trim(text) .eq. "days since 0-1-1")
     &       time(:) = time(:)/yrlen
          if (trim(text) .eq. "years since 1-1-1")
     &      time(:) = time(:) - 1.
          first_time = .true.
        endif
        iyr(:) = 0
        yrl(:) = 0.
      else
        first_time = .false.
      endif

      intrp = .false.

      data_time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      yrl(2) = min(time(ln), max(time(1), data_time))

      if (yrl(2) .gt. time(1) .and. yrl(2) .lt. time(ln)) intrp = .true.

      if (first_time .or. yrl(2) .gt. yrl(3)) then
!       read data
        ib(:) = 1
        ic(:) = 1
        ic(1) = imtm2
        ic(2) = jmtm2
        fname = new_file_name (name)
        inquire (file=trim(fname), exist=exists)
        if (exists) then
          dmsk(:,:) = 1.
          if (intrp) then
            do n=2,ln
              if (time(n-1) .le. yrl(2) .and. time(n) .ge. yrl(2)) then
                yrl(1) = time(n-1)
                iyr(1) = n-1
                yrl(3) = time(n)
                iyr(3) = n
              endif
            enddo
            call openfile (fname, iou)
            ib(3) = iyr(1)
            print*, "=> reading plastic data for year:",yrl(1)
            call getvara (trim(vname), iou, imtm2*jmtm2, ib, ic
     &,       tmpij, p1e10, c0)
            plastic(2:imtm1,2:jmtm1,1) = tmpij(1:imtm2,1:jmtm2)
           do j=1,jmt
              do i=1,imt
                 if (plastic(i,j,1).lt.0) then
                    plastic(i,j,1)=0.
                 endif
              enddo
           enddo
            call embmbc (plastic(:,:,1))
            call areatot (plastic(:,:,1), dmsk, tmp)
            if (tmp .ne. 0) plastic(:,:,1) = plastic(:,:,1)!*atmsa/tmp
            ib(3) = iyr(3)
            print*, "=> reading plastic data for year:",yrl(3)
            call getvara (trim(vname), iou, imtm2*jmtm2, ib, ic
     &,       tmpij, p1e10, c0)
            plastic(2:imtm1,2:jmtm1,3) = tmpij(1:imtm2,1:jmtm2)
           do j=1,jmt
              do i=1,imt
                 if (plastic(i,j,3).lt.0) then
                    plastic(i,j,3)=0.
                 endif
              enddo
           enddo
            call embmbc (plastic(:,:,3))
            call areatot (plastic(:,:,3), dmsk, tmp)
            if (tmp .ne. 0) plastic(:,:,3) = plastic(:,:,3)!*atmsa/tmp
          else
            if (yrl(2) .le. time(1)) then
              n = 1
              yrl(1) = time(1)
              yrl(3) = time(1)
              iyr(n) = 1
            else
              n = 3
              yrl(1) = time(ln)
              yrl(3) = time(ln)
              iyr(n) = ln
            endif
            call openfile (fname, iou)
            ib(3) = iyr(n)
            print*, "=> reading plastic data for year:",yrl(2)
            call getvara (trim(vname), iou, imtm2*jmtm2, ib, ic
     &,       tmpij, p1e10, c0)
            plastic(2:imtm1,2:jmtm1,2) = tmpij(1:imtm2,1:jmtm2)
           do j=1,jmt
              do i=1,imt
                 if (plastic(i,j,2).lt.0) then
                    plastic(i,j,2)=0.
                 endif
              enddo
           enddo
            call embmbc (plastic(:,:,2))
            call areatot (plastic(:,:,2), dmsk, tmp)
            if (tmp .ne. 0) plastic(:,:,2) = plastic(:,:,2)!*atmsa/tmp
            plastic(:,:,1) = plastic(:,:,2)
            plastic(:,:,3) = plastic(:,:,2)
          endif
        endif
      endif

      if (intrp) then
!       interpolate data
        wt1 = 1.
        if (yrl(3) .ne. yrl(1)) wt1 = (yrl(3)-yrl(2))/(yrl(3)-yrl(1))
        wt1 = max(0., min(1., wt1))
        wt3 = 1. - wt1
        do j=1,jmt
          do i=1,imt
            plastic(i,j,2) = plastic(i,j,1)*wt1 + plastic(i,j,3)*wt3
          enddo
        enddo
      elseif (yrl(2) .le. time(1)) then
        plastic(:,:,2) = plastic(:,:,1)
      elseif (yrl(2) .ge. time(ln)) then
        plastic(:,:,2) = plastic(:,:,3)
      endif
      call embmbc (plastic(1,1,2))
!      print*,'KK2',plastic(:,:,2), atmsa, tmp

      return
      end

#endif
