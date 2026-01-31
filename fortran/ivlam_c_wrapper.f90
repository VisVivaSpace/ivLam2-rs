! ivlam_c_wrapper.f90 — Thin bind(C) wrapper for ivLam2 Fortran library.
! Provides C-callable entry points for Rust FFI.
!
! NOTE: ivLam assumes mu=1. Callers must scale inputs/outputs accordingly.

module ivlam_c_wrapper
  use iso_c_binding
  implicit none

contains

  subroutine ivlam_init_c(nmax, path, pathlen, info) bind(C, name="ivlam_init_c")
    integer(c_int), intent(in), value :: nmax
    type(c_ptr), intent(in), value :: path
    integer(c_int), intent(in), value :: pathlen
    integer(c_int), intent(out) :: info

    character(len=:), allocatable :: fpath
    character(len=1), pointer :: path_chars(:)
    integer :: i, info_internal

    call c_f_pointer(path, path_chars, [pathlen])
    allocate(character(len=pathlen) :: fpath)
    do i = 1, pathlen
      fpath(i:i) = path_chars(i)
    end do

    call ivLam_initialize(nmax, fpath, info_internal)
    info = int(info_internal, c_int)

    deallocate(fpath)
  end subroutine

  subroutine ivlam_zero_rev_c(r1, r2, tof, direction, v1, v2, info, halfrev) &
      bind(C, name="ivlam_zero_rev_c")
    real(c_double), intent(in) :: r1(3), r2(3), tof
    integer(c_int), intent(in), value :: direction
    real(c_double), intent(out) :: v1(3), v2(3)
    integer(c_int), intent(out) :: info, halfrev

    integer :: dir_internal, info_internal, halfrev_internal
    real(8) :: r1_f(3), r2_f(3), tof_f, v1_f(3), v2_f(3)

    r1_f = r1
    r2_f = r2
    tof_f = tof
    dir_internal = direction

    call ivLam_zeroRev(r1_f, r2_f, tof_f, dir_internal, v1_f, v2_f, &
                       info_internal, halfrev_internal)

    v1 = v1_f
    v2 = v2_f
    info = int(info_internal, c_int)
    halfrev = int(halfrev_internal, c_int)
  end subroutine

  subroutine ivlam_ntilde_with_derivs_c(r1, r2, tof, direction, ntilde, &
      v1, v2, info, halfrev, include_second, dzdy_t, d2zdy_t) &
      bind(C, name="ivlam_ntilde_with_derivs_c")
    real(c_double), intent(in) :: r1(3), r2(3), tof
    integer(c_int), intent(in), value :: direction, ntilde
    real(c_double), intent(out) :: v1(3), v2(3)
    integer(c_int), intent(out) :: info, halfrev
    integer(c_int), intent(in), value :: include_second
    real(c_double), intent(out) :: dzdy_t(7, 6)
    real(c_double), intent(out) :: d2zdy_t(7, 7, 6)

    integer :: dir_i, ntilde_i, info_i, halfrev_i
    logical :: inc_second
    real(8) :: r1_f(3), r2_f(3), tof_f, v1_f(3), v2_f(3)
    real(8) :: dzdy_f(7, 6), d2zdy_f(7, 7, 6)

    r1_f = r1
    r2_f = r2
    tof_f = tof
    dir_i = direction
    ntilde_i = ntilde
    inc_second = (include_second /= 0)

    call ivLam_NtildeWithDerivs(r1_f, r2_f, tof_f, dir_i, ntilde_i, &
                                v1_f, v2_f, info_i, halfrev_i, &
                                inc_second, dzdy_f, d2zdy_f)

    v1 = v1_f
    v2 = v2_f
    info = int(info_i, c_int)
    halfrev = int(halfrev_i, c_int)
    dzdy_t = dzdy_f
    d2zdy_t = d2zdy_f
  end subroutine

end module ivlam_c_wrapper
