#![forbid(unsafe_code)]

/// Radiation/spreading from point sources on ternary grids.

/// Set cells within Manhattan distance `radius` of (source_x, source_y) to `intensity`.
pub fn irradiate(grid: &mut [i8], width: usize, source_x: usize, source_y: usize, intensity: i8, radius: usize) {
    let height = grid.len() / width;
    for y in 0..height {
        for x in 0..width {
            let dist = (x as isize - source_x as isize).unsigned_abs()
                + (y as isize - source_y as isize).unsigned_abs();
            if dist <= radius {
                let idx = y * width + x;
                grid[idx] = intensity;
            }
        }
    }
}

/// Exponential decay of all nonzero cells toward 0 by half_life factor.
pub fn decay(grid: &mut [i8], _width: usize, half_life: f64) {
    let factor = 0.5f64.powf(1.0 / half_life);
    for cell in grid.iter_mut() {
        if *cell != 0 {
            let v = *cell as f64 * factor;
            let rounded = v.round() as i8;
            // Preserve sign but don't cross zero
            if (*cell > 0 && rounded <= 0) || (*cell < 0 && rounded >= 0) {
                *cell = 0;
            } else {
                *cell = rounded;
            }
        }
    }
}

/// Diffuse values to 4-connected neighbors by `rate` fraction.
pub fn spread(grid: &mut [i8], width: usize, rate: f64) {
    let len = grid.len();
    let mut delta = vec![0i16; len]; // use i16 to avoid overflow during accumulation
    for i in 0..len {
        if grid[i] == 0 { continue; }
        let contribution = (grid[i] as f64 * rate) as i8;
        if contribution == 0 { continue; }
        let ns = _neighbors4(i, width, len);
        let share = contribution / ns.len() as i8;
        for &n in &ns {
            delta[n] += share as i16;
        }
        // Subtract from source to conserve mass
        delta[i] -= contribution as i16;
    }
    for i in 0..len {
        let new_val = grid[i] as i16 + delta[i];
        grid[i] = new_val.clamp(-128, 127) as i8;
    }
}

fn _neighbors4(idx: usize, width: usize, len: usize) -> Vec<usize> {
    let mut ns = Vec::new();
    if idx >= width { ns.push(idx - width); }
    if idx % width > 0 { ns.push(idx - 1); }
    if idx % width + 1 < width { ns.push(idx + 1); }
    if idx + width < len { ns.push(idx + width); }
    ns
}

/// Compute cumulative dose map from multiple sources. Each source: (x, y, intensity).
/// Uses inverse-Manhattan-distance weighting.
pub fn dose_map(sources: &[(usize, usize, i8)], width: usize, height: usize) -> Vec<f64> {
    let mut dose = vec![0.0f64; width * height];
    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            for &(sx, sy, intensity) in sources {
                let dist = (x as isize - sx as isize).unsigned_abs()
                    + (y as isize - sy as isize).unsigned_abs();
                if dist == 0 {
                    dose[idx] += intensity as f64;
                } else {
                    dose[idx] += intensity as f64 / (dist as f64 * dist as f64);
                }
            }
        }
    }
    dose
}

/// Compute dose with shielding: cells with value -1 in `grid` block radiation.
/// A -1 cell absorbs all radiation passing through it (simple model: ray-based Manhattan check).
pub fn shielding(grid: &[i8], width: usize, source_x: usize, source_y: usize) -> Vec<f64> {
    let height = grid.len() / width;
    let mut dose = vec![0.0f64; width * height];
    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            // Check if line from source to (x,y) passes through any -1 cell
            let blocked = _line_blocked(grid, width, source_x, source_y, x, y);
            if !blocked {
                let dist = (x as isize - source_x as isize).unsigned_abs()
                    + (y as isize - source_y as isize).unsigned_abs();
                dose[idx] = if dist == 0 { 1.0 } else { 1.0 / dist as f64 };
            }
        }
    }
    dose
}

/// Simple Bresenham-like check for shields along Manhattan path.
fn _line_blocked(grid: &[i8], width: usize, x0: usize, y0: usize, x1: usize, y1: usize) -> bool {
    if x0 == x1 && y0 == y1 { return false; }
    // Walk along x then y (L-shaped path)
    let mut cx = x0 as isize;
    let mut cy = y0 as isize;
    let dx = (x1 as isize - x0 as isize).signum();
    let dy = (y1 as isize - y0 as isize).signum();

    // Walk x
    while cx != x1 as isize {
        cx += dx;
        if cx >= 0 && cy >= 0 {
            let ux = cx as usize;
            let uy = cy as usize;
            if ux < width && uy < grid.len() / width {
                let idx = uy * width + ux;
                if grid[idx] == -1 && !(cx == x1 as isize && cy == y1 as isize) {
                    return true;
                }
            }
        }
    }
    // Walk y
    while cy != y1 as isize {
        cy += dy;
        if cx >= 0 && cy >= 0 {
            let ux = cx as usize;
            let uy = cy as usize;
            if ux < width && uy < grid.len() / width {
                let idx = uy * width + ux;
                if grid[idx] == -1 && !(cx == x1 as isize && cy == y1 as isize) {
                    return true;
                }
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_irradiate_center() {
        let mut grid = vec![0i8; 9];
        irradiate(&mut grid, 3, 1, 1, 5, 1);
        assert_eq!(grid[4], 5); // center
        assert_eq!(grid[1], 5); // above
        assert_eq!(grid[7], 5); // below
        assert_eq!(grid[3], 5); // left
        assert_eq!(grid[5], 5); // right
        assert_eq!(grid[0], 0); // diagonal, out of Manhattan radius 1
    }

    #[test]
    fn test_irradiate_radius_zero() {
        let mut grid = vec![0i8; 9];
        irradiate(&mut grid, 3, 1, 1, 7, 0);
        assert_eq!(grid[4], 7); // only center
        assert_eq!(grid[0], 0);
    }

    #[test]
    fn test_irradiate_large_radius() {
        let mut grid = vec![0i8; 9];
        irradiate(&mut grid, 3, 1, 1, 1, 10);
        assert!(grid.iter().all(|&v| v == 1)); // all cells hit
    }

    #[test]
    fn test_decay_simple() {
        let mut grid = vec![10i8, -10];
        decay(&mut grid, 2, 1.0); // half_life=1
        assert!(grid[0] < 10 && grid[0] > 0);
        assert!(grid[1] > -10 && grid[1] < 0);
    }

    #[test]
    fn test_decay_zeros_unchanged() {
        let mut grid = vec![0i8; 4];
        decay(&mut grid, 2, 5.0);
        assert!(grid.iter().all(|&v| v == 0));
    }

    #[test]
    fn test_decay_reduces_magnitude() {
        let mut grid = vec![100i8];
        decay(&mut grid, 1, 2.0);
        assert!(grid[0] < 100);
    }

    #[test]
    fn test_spread_basic() {
        let mut grid = vec![0i8; 9];
        grid[4] = 10; // center
        spread(&mut grid, 3, 0.5);
        // Neighbors should get some value
        assert!(grid[1] > 0 || grid[3] > 0 || grid[5] > 0 || grid[7] > 0);
    }

    #[test]
    fn test_spread_preserves_total_approx() {
        let mut grid = vec![0i8; 9];
        grid[4] = 8;
        spread(&mut grid, 3, 0.5);
        let total: i32 = grid.iter().map(|&v| v as i32).sum();
        assert!(total >= 6 && total <= 10, "total was {total}"); // approximately conserved
    }

    #[test]
    fn test_dose_map_single_source() {
        let dose = dose_map(&[(2, 2, 10)], 5, 5);
        assert_eq!(dose[12], 10.0); // at source
        // Nearby cell should have positive dose
        assert!(dose[11] > 0.0);
    }

    #[test]
    fn test_dose_map_two_sources() {
        let dose = dose_map(&[(0, 0, 5), (4, 4, 5)], 5, 5);
        assert!((dose[0] - 5.0).abs() < 1.0); // ≈5 at first source
        assert!((dose[24] - 5.0).abs() < 1.0); // ≈5 at second source
    }

    #[test]
    fn test_shielding_no_shield() {
        let grid = vec![0i8; 9];
        let dose = shielding(&grid, 3, 1, 1);
        assert!(dose[4] > 0.0); // source itself
        assert!(dose[0] > 0.0); // no shield, should get dose
    }

    #[test]
    fn test_shielding_with_shield() {
        let mut grid = vec![0i8; 15]; // 5x3
        grid[7] = -1; // shield at (2,1)
        let dose = shielding(&grid, 5, 4, 1); // source at (4,1)
        // source itself should have dose
        assert!(dose[9] > 0.0); // (4,1) = idx 9
        // Cells behind shield get zero or less dose
        assert_eq!(dose[5], 0.0); // (0,1) behind shield
    }

    #[test]
    fn test_dose_map_empty_sources() {
        let dose = dose_map(&[], 3, 3);
        assert!(dose.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_irradiate_corner() {
        let mut grid = vec![0i8; 9];
        irradiate(&mut grid, 3, 0, 0, 3, 1);
        assert_eq!(grid[0], 3);
        assert_eq!(grid[1], 3); // right
        assert_eq!(grid[3], 3); // below
        assert_eq!(grid[8], 0); // far corner
    }
}
