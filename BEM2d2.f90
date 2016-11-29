! ********************  B E M 2 d  ***********************
!
! A 2D incompressible, unsteady panel code using constant
! strength source/doublet panels for the perturbation
! velocity potential.
!
! AUTHOR: Trevor Wood,  Boston University
! DATE:   c. 2000 
! Under the direction of  Professor Sheryl Grace, Boston University
!
! ********************************************************


module BEM2d_main

use TypeKinds
use BEM2d_typ
use Globals
use IO
use Geometry
use BEM2d_mod
implicit none
private

public  :: Initialize, Output
private :: InitOldRun, InitNewRun

contains

subroutine Initialize(W)
   type(Wing)             , intent(inout) :: W

   integer, dimension(:,:), pointer :: IdxWake
   real(kind=dp) :: AOA, Tf
   integer       :: i,ii, Nx, Nt, fid, iargc
   character(len=4)                 :: str
   
   nullify(W%Xc,W%B,W%W)
   if(iargc() < 1) then
      Filetag = " "
   else
      call getarg(1,Filetag)
      i = len_trim(Filetag)+1
      Filetag(i:i) = "."
   end if

   ! read input data
   open(unit=1,file="BEM."//trim(Filetag)//"in", action="read",status="old")
   read(unit=1,fmt=*) GeomPath
   read(unit=1,fmt=*) DoField, DoVort, FreeVort, Fixed_Wake, Steady, WriteDataFile
   read(unit=1,fmt=*) Mach,      dn,       delta,  AOA
   read(unit=1,fmt=*) Tf,        Nt
   read(unit=1,fmt=*) TEcond,   BCmode,    BCfreq
   read(unit=1,fmt=*) Npieces

   Beta  = sqrt(1.0_dp - Mach**2)

   print *, "GeomPath: ", trim(GeomPath)
   print *, "Filetag:  ", trim(Filetag)
   if(Steady) then
      print *, "Performing steady-state calculation"
   else
      print *, "Performing time-dependent calculation",FreeVort
   end if
 
   ! Initialize Time vector
   nullify(Time)
   allocate(Time(Nt))
   Time = Tf * real((/ (ii, ii = 1, Nt) /), kind=dp) / real(Nt,kind=dp)

   if(Npieces == 0)then
      call InitOldRun(W)
   else
      call InitNewRun(W,AOA,Nt)
   end if
   close(unit=1)

   call getCollPts(W)

   call openmatfile("BEM_data."//trim(Filetag),"w", fid, WriteDataFile)
   call savematfile(fid,"Time", Time, WriteDataFile)
   if(WriteDataFile)then
      write(unit=fid,fmt=*) "Npieces = ",Npieces," ;"
   else
      write(unit=fid) Npieces
   end if
   do i = 1, Npieces
      write(unit=str,fmt="( i4)")i
      call savematfile(fid,"Xb{"//str//"}", W%B(i,1)%X, WriteDataFile)
   end do
   call savematfile(fid,"Xc"  , W%Xc, WriteDataFile)
   call closematfile(fid)
end subroutine Initialize

subroutine InitNewRun(W,AOA,Nt)
   type(Wing)   , intent(inout) :: W
   real(kind=dp), intent(in)    :: AOA
   integer      , intent(in)    :: Nt

   real(kind=dp), dimension(:,:), pointer :: Xbtmp
   
   integer, dimension(:,:), pointer :: IdxWake
   character(len=max_str_len)       :: AFtype, AFspacing
   character(len=4)                 :: str
   real(kind=dp), dimension(2)      :: trans,scales
   integer                          :: i, Nx

   ! Initialize Airfoil shape
   nullify(W%B,W%W,W%idx, IdxWake,Xbtmp)
   allocate(W%B(Npieces,Nt),W%W(Npieces,Nt),W%idx(Npieces))

   do i = 1, Npieces
      read(unit=1,fmt=*) AFtype, Nx, AFspacing
      read(unit=1,fmt=*) trans,  scales
      call afgeom(AFtype,Nx,Xbtmp,IdxWake,AFspacing)
         Nx = size(Xbtmp,1)
         nullify(W%B(i,1)%X, W%W(i,1)%X, W%B(i,1)%phi, W%W(i,1)%phi)
         allocate(W%B(i,1)%X(Nx,2)  , W%W(i,1)%X(2,2))
         allocate(W%B(i,1)%phi(Nx-1), W%W(i,1)%phi(1))
         W%B(i,1)%X = Xbtmp(Nx:1:-1,(/1,3/))   ! Clockwise orient_n of AF geom !
         deallocate(Xbtmp,IdxWake)
      call aforient(W%B(i,1)%X,trans,scales,-AOA*pi/180.0_dp)

      call WakeTE(W%B(i,1),W%W(i,1))
   end do ! Npieces
end subroutine InitNewRun

subroutine InitOldRun(W)
   type(Wing), intent(inout)   :: W

   real(kind=dp), dimension(:), pointer :: t
   
   character(len=max_str_len)  :: ftag
   character(len=11)           :: tmp
   character(len=13)           :: cell
   character(len=4)            :: str
   logical                     :: ftag_asc
   integer                     :: fid, Nt, i,j, N

   read(unit=1,fmt=*) ftag
   read(unit=1,fmt=*) ftag_asc

   print *, "InitOldRun: Reading data from prior run"

   call openmatfile("BEM_data."//trim(ftag)//".","r", fid, ftag_asc)
   call readmatfile(fid,"Time", t, ftag_asc)
   Nt = size(Time)
	ot = .true.
	oldtime = t(2) - t(1)
   if(ftag_asc)then
      read(unit=fid,fmt="(a11,i)") tmp,Npieces
   else
      read(unit=fid) Npieces
   end if

   ! Initialize Airfoil shape
   nullify(W%B,W%W,W%idx)
   allocate(W%B(Npieces,Nt),W%W(Npieces,Nt),W%idx(Npieces))
   do i = 1, Npieces
      write(unit=str,fmt="( i4)")i
      call readmatfile(fid,"Xb{"//str//"}",W%B(i,1)%X,ftag_asc,(/i/))
   end do
   call readmatfile(fid,"Xc",W%Xc,ftag_asc,(/i/))

   ! Read through data file of prior run, storing only final timestep
   !!! NOTE: only valid for incompressible flow !!!
   do j = 1, size(t)
   do i = 1, Npieces
      write(unit=cell,fmt="(""{"",i5,"","",i5,""}"")") i,Timeiter
      call readmatfile(fid,"Phi" //cell, W%B(i,1)%phi, ftag_asc )
      call readmatfile(fid,"Xw"  //cell, W%W(i,1)%X  , ftag_asc )
      call readmatfile(fid,"Phiw"//cell, W%W(i,1)%phi, ftag_asc )
   end do
   end do
   call closematfile(fid)
   deallocate(t)

   Timeiter = 1
!   call EvolveWake(W,Time(1))
!   do i = 1, Npieces
!      deallocate(W%W(i,1)%X,W%W(i,1)%phi)
!      nullify(W%W(i,1)%X,W%W(i,1)%phi)
!      N = size(W%W(i,2)%X,1)
!      allocate(W%W(i,1)%X(N,2),W%W(i,1)%phi(N-1))
!      W%W(i,1)%X   = W%W(i,2)%X
!      W%W(i,1)%phi = W%W(i,2)%phi
!      deallocate(W%W(i,2)%X,W%W(i,2)%phi,W%B(i,2)%phi)
!   end do
   print *, "InitOldRun: Done"
end subroutine InitOldRun

subroutine Output(W)
   type(Wing)       , intent(in)     :: W
   integer                           :: fid, i
   character(len=13)                 :: cell

   call openmatfile("BEM_data."//trim(Filetag),"u", fid, WriteDataFile)
   do i = 1, Npieces
      write(unit=cell,fmt="(""{"",i5,"","",i5,""}"")") i,Timeiter
      call savematfile(fid,"Phi" //cell, W%B(i,Timeiter)%phi, WriteDataFile )
      call savematfile(fid,"Xw"  //cell, W%W(i,Timeiter)%X  , WriteDataFile )
      call savematfile(fid,"Phiw"//cell, W%W(i,Timeiter)%phi, WriteDataFile )
   end do
   call closematfile(fid)
end subroutine Output

end module BEM2d_main



program BEM2d
   use TypeKinds
   use BEM2d_typ
   use Globals
   use Vortex_mod
   use BEM2d_mod
   use Farfield
   use BEM2d_main
   implicit none

   type(Wing)                             :: W
   real(kind=dp)                          :: dt
   integer                                :: i

   call Initialize(W)
   if(DoField)then
      call InitFField()
   end if
   if(DoVort)then
      call InitVortex()

   end if

   Timeiter = 1
   dt = Time(1)

   call Solve(W, .false.)
   if(DoField)then
      call CalcFField(W, .false.)
   end if
   
   call Output(W)
   print "("" iteration: "",i4,"", time: "",f10.5)", 1, Time(1)
   
   do i = 2, size(Time)
      dt = Time(i) - Time(i-1)
      call EvolveWake(W,dt)
      if(DoVort)then
         call EvolVort(W,dt)
      end if
      Timeiter = i

      call Solve(W, .false.)
      if(DoField)then
         call CalcFField(W, .false.)
      end if

      call Output(W)
      print "("" iteration: "",i4,"", time: "",f10.5)", i, Time(i)
   end do

   if(DoField)then
      call WriteFField()
   end if
   if(DoVort)then
      call OutputVort()
   !else
	!   call OutputPhi(W,size(W%B(1,1)%phi+1))
   end if

	!if (Timeiter .eq. size(Time)) then
   !   print *, 'Calling Velfield'
	!	call Velfield(W,dt)
	!endif

   stop
end program BEM2d