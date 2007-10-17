module precd

  type precdata
     character(len=10)  :: lv1, lv2    ! First and second level prec type
     integer            :: nlev
     integer            :: novr        ! number of overlapping levels
     integer            :: restr       ! restriction  over application of as
     integer            :: prol        ! prolongation over application of as
     integer            :: ftype1      ! Factorization type: ILU, SuperLU, UMFPACK. 
     integer            :: fill1       ! Fill-in for factorization 1
     real(kind(1.d0))   :: thr1        ! Threshold for fact. 1 ILU(T)
     integer            :: mltype      ! additive or multiplicative 2nd level prec
     integer            :: aggr        ! local or global aggregation
     integer            :: smthkind    ! smoothing type
     integer            :: cmat        ! coarse mat
     integer            :: smthpos     ! pre, post, both smoothing
     integer            :: glbsmth     ! global smoothing
     integer            :: ftype2      ! Factorization type: ILU, SuperLU, UMFPACK. 
     integer            :: fill2       ! Fill-in for factorization 1
     real(kind(1.d0))   :: thr2        ! Threshold for fact. 1 ILU(T)
     integer            :: jswp        ! Jacobi sweeps
     real(kind(1.d0))   :: omega       ! smoother omega
     character(len=40)  :: descr       ! verbose description of the prec
  end type precdata

end module precd
