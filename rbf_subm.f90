
!----------------------
!*** END MLS
!*** START RBF MESHLESS
!----------------------

!---------------------------
!*** radius and derivatives
!*** determined
!*** chk0,1,2
!---------------------------
  SUBROUTINE rbf_radius(ndi, x, xi, r, drdx, d2rdx2)
    IMPLICIT REAL(8) (a - h, o - z)
    REAL(8), PARAMETER::small = 1.0d-10
    REAL(8), DIMENSION(ndi)::x, xi
    REAL(8), DIMENSION(ndi)::drdx
    REAL(8), DIMENSION(ndi, ndi)::d2rdx2
    r = rnorm2(ndi, x - xi)
    drdx = 0.0d00
    d2rdx2 = 0.0d00
    IF (r .GT. small) THEN
      DO id = 1, ndi
        drdx(id) = (x(id) - xi(id))/r
        DO jd = 1, ndi
          d2rdx2(id, jd) = deltak(id, jd)/r - (1.0d00/(r**3.0d00))*(x(id) - xi(id))*(x(jd) - xi(jd))
        END DO
      END DO
    END IF
  END SUBROUTINE rbf_radius

!--------------------------
!*** radial basis function
!*** basis function
!*** chk0,1
!--------------------------
  SUBROUTINE rbf_c2(rmax, r, a, dadr, d2adr2)
    IMPLICIT REAL(8) (a - h, o - z)
    REAL(8)::rmax, r, a, dadr, d2adr2
    REAL(8), PARAMETER::small = 1.0d-10
    a = 0.0d00
    dadr = 0.0d00
    d2adr2 = 0.0d00
    IF (r .LE. 0.99999d00*rmax) THEN
      a = ((1.0d00 - (r/rmax))**4)*(1.0d00 + 4.0d00*(r/rmax))
      dadr = (20.0d00*r*(r - rmax)**3.0d00)/(rmax**5.0d00)
      d2adr2 = 20.0d00*((r - rmax)**2.0d00)*(4.0d00*r - rmax)/(rmax**5.0d00)
    END IF
  END SUBROUTINE rbf_c2

!------------------------------------
!*** radial and derivatives
!*** xi is the "node" coordinate
!*** x is the gauss point coordinate
!*** chk0,1
!------------------------------------
  SUBROUTINE rbf_one(rmax, ndi, x, xi, a, dadx, d2adx2)
    IMPLICIT REAL(8) (a - h, o - z)
    REAL(8)::rmax, r, a, dadr, d2adr2
    REAL(8), DIMENSION(ndi)::x, xi
    REAL(8), DIMENSION(ndi)::dadx, drdx
    REAL(8), DIMENSION(ndi, ndi)::d2adx2, d2rdx2
    CALL rbf_radius(ndi, x, xi, r, drdx, d2rdx2)
    CALL rbf_c2(rmax, r, a, dadr, d2adr2)
    DO id = 1, ndi
      dadx(id) = dadr*drdx(id)
      DO jd = 1, ndi
        d2adx2(id, jd) = d2adr2*drdx(id)*drdx(jd) + dadr*d2rdx2(id, jd)
      END DO
    END DO
  END SUBROUTINE rbf_one

!-------------------------------------------
!*** obtain shape functions and derivatives
!*** from a given list of nodes
!*** at one point with coordinates x
!*** chk0
!-------------------------------------------
  SUBROUTINE rbf_shapefunctions(imeshless, rmax, ndi, n, xn, x, xc, ff, dff, dff2)
    IMPLICIT REAL(8) (a - h, o - z)
    INTEGER, PARAMETER::mpol = 20
    REAL(8)::rmax
    REAL(8), DIMENSION(n, n)::amatrix, invmatrix, gmatrix
    REAL(8), DIMENSION(n, mpol)::bmatrix
    REAL(8), DIMENSION(mpol, n)::hmatrix
    REAL(8), DIMENSION(ndi, ndi, n)::dff2
    REAL(8), DIMENSION(ndi, n)::dff
    REAL(8), DIMENSION(n)::ff, a
    REAL(8), DIMENSION(ndi, n)::xn
    REAL(8), DIMENSION(ndi)::x, xc
    REAL(8), DIMENSION(mpol)::b
    REAL(8), DIMENSION(ndi, mpol)::dbdx
    REAL(8), DIMENSION(ndi, ndi, mpol)::dbdx2
    REAL(8), DIMENSION(ndi, ndi, mpol)::d2bdx2
    REAL(8), DIMENSION(mpol, mpol)::b1b, invb1b
    REAL(8), DIMENSION(ndi, n)::dadx
    REAL(8), DIMENSION(ndi, ndi, n)::d2adx2
    REAL(8), DIMENSION(mpol)::polyn
    REAL(8), DIMENSION(ndi, mpol)::dpolyn
    REAL(8), DIMENSION(ndi, ndi, mpol)::dpolyn2
!------------------
!*** forms amatrix
!------------------
    DO jn=1,n
      DO in=1,n
!-----------------------------------
!*** dadx and d2adx2 are trash here
!-----------------------------------
        CALL rbf_one(rmax, ndi, xn(1:ndi, in), xn(1:ndi, jn), amatrix(in, jn), dadx(1:ndi, in), d2adx2(1:ndi, 1:ndi, in))
      END DO
    END DO
!--------------------------
!*** invmatrix=amatrix^-1
!--------------------------
    CALL invmat(n, deta, amatrix, invmatrix)
!----------------------------------
!*** forms bmatrix
!*** checks effect of dtemp and xc
!----------------------------------
    mp=0
    IF(imeshless.NE.0)THEN
      DO in = 1, n
        mp = 0
        SELECT CASE (imeshless)
        CASE (1)
          CALL mls_polynlinbase(ndi, 1.0d00, xn(1:ndi, in), xc(1:ndi), polyn(1:mpol), dpolyn(1:ndi, 1:mpol), dpolyn2(1:ndi, 1:ndi, 1:mpol), mp)
        CASE (2)
          CALL mls_polynqbase(ndi, 1.0d00, xn(1:ndi, in), xc(1:ndi), polyn(1:mpol), dpolyn(1:ndi, 1:mpol), dpolyn2(1:ndi, 1:ndi, 1:mpol), mp)
        CASE (3)
          CALL mls_polyncubbase(ndi, 1.0d00, xn(1:ndi, in), xc(1:ndi), polyn(1:mpol), dpolyn(1:ndi, 1:mpol), dpolyn2(1:ndi, 1:ndi, 1:mpol), mp)
        END SELECT
        DO ip = 1, mp
          bmatrix(in, ip) = polyn(ip)
        END DO
      END DO
    END IF
!-----------------------
!*** forms big matrices
!-----------------------
    IF(imeshless.NE.0)THEN
      DO jd = 1, mp
        DO id = 1, mp
          b1b(id, jd) = 0.0d00
          DO in = 1, n
            DO jn = 1, n
              b1b(id, jd) = b1b(id, jd) + bmatrix(in, id)*invmatrix(in, jn)*bmatrix(jn, jd)
            END DO
          END DO
        END DO
      END DO
    END IF
!-----------------------
!*** !!!above seems ok
!-----------------------
!----------------
!*** invert b1b
!*** into invb1b
!----------------
    IF(mp.GT.0)THEN
      CALL invmat(mp, det, b1b(1:mp, 1:mp), invb1b(1:mp, 1:mp))
!------------
!*** hmatrix
!------------
      DO in = 1, n
        DO id = 1, mp
          hmatrix(id, in) = 0.0d00
          DO jn = 1, n
            DO jd = 1, mp
              hmatrix(id, in) = hmatrix(id, in) + invb1b(id, jd)*bmatrix(jn, jd)*invmatrix(jn, in)
            END DO
          END DO
        END DO
      END DO
      DO jn = 1, n
        DO in = 1, n
          gmatrix(in, jn) = invmatrix(in, jn)
          DO id = 1, mp
            DO kn = 1, n
              gmatrix(in, jn) = gmatrix(in, jn) - invmatrix(in, kn)*bmatrix(kn, id)*hmatrix(id, jn)
            END DO
          END DO
        END DO
      END DO
    END IF
!----------------------------------
!*** forms vectors a and b
!*** checks effect of dtemp and xc
!----------------------------------
    DO in = 1, n
      CALL rbf_one(rmax, ndi, x(1:ndi), xn(1:ndi, in), a(in), dadx(1:ndi, in), d2adx2(1:ndi, 1:ndi, in))
      SELECT CASE (imeshless)
      CASE (1)
        CALL mls_polynlinbase(ndi, 1.0d00, x(1:ndi), xc(1:ndi), b, dbdx(1:ndi, 1:mp), dbdx2(1:ndi, 1:ndi, 1:mp), mp)
      CASE (2)
        CALL mls_polynqbase(ndi, 1.0d00, x(1:ndi), xc(1:ndi), b, dbdx(1:ndi, 1:mp), dbdx2(1:ndi, 1:ndi, 1:mp), mp)
      CASE (3)
        CALL mls_polyncubbase(ndi, 1.0d00, x(1:ndi), xc(1:ndi), b, dbdx(1:ndi, 1:mp), dbdx2(1:ndi, 1:ndi, 1:mp), mp)
      END SELECT
    END DO
!------------------------------------------------------
!*** now determines rb shape functions and derivatives
!------------------------------------------------------
    IF (imeshless .EQ. 0) THEN
      DO in = 1, n
        ff(in) = dotprod(n, a(1:n), invmatrix(1:n, in))
        DO id = 1, ndi
          dff(id, in) = dotprod(n, dadx(id, 1:n), invmatrix(1:n, in))
          DO jd = 1, ndi
            dff2(id, jd, in) = dotprod(n, d2adx2(id, jd, 1:n), invmatrix(1:n, in))
          END DO
        END DO
      END DO
    ELSE
      DO in = 1, n
        ff(in) = dotprod(n, a(1:n), gmatrix(1:n, in)) + dotprod(mp, b(1:mp), hmatrix(1:mp, in))
        DO id = 1, ndi
          dff(id, in) = dotprod(n, dadx(id, 1:n), gmatrix(1:n, in)) + dotprod(mp, dbdx(id, 1:mp), hmatrix(1:mp, in))
          DO jd = 1, ndi
            dff2(id, jd, in) = dotprod(n, d2adx2(id, jd, 1:n), gmatrix(1:n, in)) + dotprod(mp, dbdx2(id, jd, 1:mp), hmatrix(1:mp, in))
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE rbf_shapefunctions

!-------------------------------------------------
!*** from an element support estimate, determines
!*** gauss point support radius
!-------------------------------------------------
  SUBROUTINE rbf_gausslevelshapefunctions(imeshless, rcentroid, ndi, xg, xc, ntot, xtot, ff, dff, dff2)
    IMPLICIT REAL(8) (a - h, o - z)
    INTEGER::imeshless
    REAL(8)::rcentroid, rgauss
    REAL(8), DIMENSION(ndi)::xg, xc
    REAL(8), DIMENSION(ntot)::ff
    REAL(8), DIMENSION(ndi, ntot)::dff
    REAL(8), DIMENSION(ndi, ndi, ntot)::dff2
    REAL(8), DIMENSION(ndi, ntot)::xtot
    rgauss = rcentroid-rnorm2(ndi,xg-xc)
    CALL rbf_testingallrequired(imeshless, rgauss, ndi, ntot, xtot, xg, xc, ff, dff, dff2)
  END SUBROUTINE rbf_gausslevelshapefunctions

!---------------------------------------------------------------------
!*** for representation purposes only:
!*** get all nodes, makes a selection
!*** and returns all quantities for all nodes, not only the selection
!---------------------------------------------------------------------
  SUBROUTINE rbf_testingallrequired(imeshless, rmax, ndi, ntot, xtot, x, xc, ff, dff, dff2)
    IMPLICIT REAL(8) (a - h, o - z)
    INTEGER::imeshless
    INTEGER, PARAMETER::mpol = 20
    INTEGER, DIMENSION(:), ALLOCATABLE::listn
    REAL(8)::rmax
    REAL(8), DIMENSION(ntot)::ff, ffloc
    REAL(8), DIMENSION(ndi, ntot)::dff, dffloc
    REAL(8), DIMENSION(ndi, ndi, ntot)::dff2, dff2loc
    REAL(8), DIMENSION(ndi, ntot)::xtot, xn
    REAL(8), DIMENSION(ndi)::x, xc
    REAL(8), DIMENSION(mpol)::polyn
    REAL(8), DIMENSION(ndi, mpol)::dpolyn
    REAL(8), DIMENSION(ndi, ndi, mpol)::dpolyn2
    REAL(8), DIMENSION(mpol, ntot)::u2
!-----------------------
!*** gets closest nodes
!*** to x with corected
!*** radius above
!-----------------------
    CALL mls_getsnodes(rmax, x, ndi, ntot, xtot, n, listn)
    DO i = 1, n
      DO id = 1, ndi
        xn(id, i) = xtot(id, listn(i))
      END DO
    END DO
!-------------------------------
!*** shape function derivatives
!-------------------------------
    IF (imeshless .GE. 0) THEN
      CALL rbf_shapefunctions(imeshless, rmax, ndi, n, xn, x, xc, ffloc, dffloc, dff2loc)
    ELSE
!---------
!**** mls
!---------
      tol = 1.0d-2
      CALL mls_u2atacoordinate(imeshless, rmax, tol, ndi, n, xc, xn, x, u2(1:mpol, 1:n))
      SELECT CASE (imeshless)
      CASE (-1)
        CALL mls_polynlinbase(ndi, rmax, x(1:ndi), xc(1:ndi), polyn(1:mpol), dpolyn(1:ndi, 1:mpol), dpolyn2(1:ndi, 1:ndi, 1:mpol), m)
      CASE (-2)
        CALL mls_polynqbase(ndi, rmax, x(1:ndi), xc(1:ndi), polyn(1:mpol), dpolyn(1:ndi, 1:mpol), dpolyn2(1:ndi, 1:ndi, 1:mpol), m)
      CASE (-3)
        CALL mls_polyncubbase(ndi, rmax, x(1:ndi), xc(1:ndi), polyn(1:mpol), dpolyn(1:ndi, 1:mpol), dpolyn2(1:ndi, 1:ndi, 1:mpol), m)
      END SELECT
      CALL mls_sfder(ndi, n, m, polyn(1:m), dpolyn(1:ndi, 1:m), dpolyn2(1:ndi, 1:ndi, 1:m), u2(1:m, 1:n), ffloc, dffloc, dff2loc)
    END IF
    ff = 0.0d00
    dff = 0.0d00
    dff2 = 0.0d00
    DO i = 1, n
      ff(listn(i)) = ffloc(i)
      dff(1:ndi, listn(i)) = dffloc(1:ndi, i)
      dff2(1:ndi, 1:ndi, listn(i)) = dff2loc(1:ndi, 1:ndi, i)
    END DO
    DEALLOCATE (listn)
  END SUBROUTINE rbf_testingallrequired
