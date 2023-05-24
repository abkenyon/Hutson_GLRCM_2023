      subroutine accadia(native, values, dest, dx, dy, neighbor_locs,
     &                   interped, tp, nn, np, two)
        integer np
        integer tp
        integer nn
        integer two
        integer neighbor_locs(tp, nn)
        real native(np, two)
        real values(np)
        real dest(tp, two)
        real dx(tp)
        real dy(tp)
        real interped(tp)

Cf2py intent(in) native
Cf2py intent(in) values
Cf2py intent(in) dest
Cf2py intent(in) dx
Cf2py intent(in) dy
Cf2py intent(in) neighbor_locs
Cf2py intent(in) np
Cf2py intent(in) tp
Cf2py intent(in) nn
Cf2py intent(in) two
Cf2py depend(tp) interped
Cf2py intent(out) interped

        integer p, n, nneighbors, neighbor_min_dist
        real neighbor_x, neighbor_y, an_x, an_y, analysis
        real subgrid_dx, subgrid_dy, dist, distmin
        real subgrid_x(5, 5)
        real subgrid_y(5, 5)
        real subgrid_factor(5)

        subgrid_factor(1) = -2.0
        subgrid_factor(2) = -1.0
        subgrid_factor(3) = 0.0
        subgrid_factor(4) = 1.0
        subgrid_factor(5) = 2.0

        interped(:) = 0.0

        do p = 1, tp
          nneighbors = 0
          subgrid_x(:, :) = 0.0
          subgrid_y(:, :) = 0.0

          an_x = dest(p, 1)
          an_y = dest(p, 2)

          do n = 1, nn
            if ( neighbor_locs(p, n) .gt. 0 ) then
              nneighbors = nneighbors + 1
            endif
          enddo

          if ( nneighbors .eq. 0) cycle

          subgrid_dx = dx(p) / 5.0
          subgrid_dy = dy(p) / 5.0

          do i = 1, 5
            do j = 1, 5
              subgrid_x(i, j) = an_x + subgrid_factor(j) * subgrid_dx
              subgrid_y(i, j) = an_y + subgrid_factor(i) * subgrid_dy
            enddo
          enddo

          analysis = 0.0
          do i = 1, 5
            do j = 1, 5
              distmin = subgrid_dx * 10.0
              neighbor_min_dist = 1
              do n = 1, nneighbors
                neighbor_x = native(neighbor_locs(p, n), 1)
                neighbor_y = native(neighbor_locs(p, n), 2)
                dist = ((neighbor_x - subgrid_x(i, j)) ** 2 +
     &                  (neighbor_y - subgrid_y(i, j)) ** 2) ** 0.5
                if ( dist .lt. distmin ) then
                  neighbor_min_dist = n
                endif
              enddo
              analysis = analysis +
     &                   values(neighbor_locs(p, neighbor_min_dist))
            enddo
          enddo
          interped(p) = analysis / 25.0
        enddo
      end
