//! Radiation and energy propagation on ternary grids {-1, 0, +1}.
//!
//! Provides point sources, inverse-square irradiance, diffusion,
//! shadow casting, half-life decay, cascade events, and shielding.

#![forbid(unsafe_code)]
#![no_std]

extern crate alloc;

use alloc::vec::Vec;

// ── Helpers ──────────────────────────────────────────────────────────

fn idx(x: usize, y: usize, w: usize) -> usize {
    y * w + x
}

fn in_bounds(x: i32, y: i32, w: usize, h: usize) -> bool {
    x >= 0 && y >= 0 && (x as usize) < w && (y as usize) < h
}

fn snap(v: f64, lo: f64, hi: f64) -> i8 {
    if v < lo { -1 }
    else if v > hi { 1 }
    else { 0 }
}

fn neighbors4(x: usize, y: usize, w: usize, h: usize) -> Vec<(usize, usize)> {
    let mut ns = Vec::new();
    for (dx, dy) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
        let nx = x as i32 + dx;
        let ny = y as i32 + dy;
        if in_bounds(nx, ny, w, h) {
            ns.push((nx as usize, ny as usize));
        }
    }
    ns
}

// ── RadiationSource ──────────────────────────────────────────────────

/// A point radiation source with intensity and decay rate.
#[derive(Clone, Debug)]
pub struct RadiationSource {
    pub x: usize,
    pub y: usize,
    pub intensity: f64,
    pub decay_rate: f64,
}

impl RadiationSource {
    pub fn new(x: usize, y: usize, intensity: f64, decay_rate: f64) -> Self {
        Self { x, y, intensity, decay_rate }
    }

    /// Compute irradiance at (px, py) using inverse-square law with decay.
    pub fn irradiance_at(&self, px: usize, py: usize) -> f64 {
        let dx = px as f64 - self.x as f64;
        let dy = py as f64 - self.y as f64;
        let dist_sq = dx * dx + dy * dy;
        if dist_sq < 1.0 {
            self.intensity
        } else {
            self.intensity / (1.0 + self.decay_rate * dist_sq)
        }
    }
}

// ── IrradianceField ──────────────────────────────────────────────────

/// Compute irradiance field from multiple sources on a grid.
pub fn compute_irradiance(
    sources: &[RadiationSource],
    w: usize,
    h: usize,
) -> Vec<f64> {
    let mut field = alloc::vec![0.0f64; w * h];
    for y in 0..h {
        for x in 0..w {
            let mut total = 0.0;
            for src in sources {
                total += src.irradiance_at(x, y);
            }
            field[idx(x, y, w)] = total;
        }
    }
    field
}

// ── TernaryIrradianceMap ─────────────────────────────────────────────

/// Convert continuous irradiance field to ternary grid.
pub fn to_ternary(field: &[f64], w: usize, h: usize, lo: f64, hi: f64) -> Vec<i8> {
    field.iter().map(|&v| snap(v, lo, hi)).collect()
}

// ── Diffusion ────────────────────────────────────────────────────────

/// Ternary diffusion: each cell averages its 4-neighbors, then snaps to {-1, 0, +1}.
pub fn diffuse(grid: &[i8], w: usize, h: usize, steps: usize) -> Vec<i8> {
    let mut current = grid.to_vec();
    let mut next = current.clone();

    for _ in 0..steps {
        for y in 0..h {
            for x in 0..w {
                let ns = neighbors4(x, y, w, h);
                let mut sum = current[idx(x, y, w)] as f64;
                let mut count = 1.0;
                for &(nx, ny) in &ns {
                    sum += current[idx(nx, ny, w)] as f64;
                    count += 1.0;
                }
                next[idx(x, y, w)] = snap(sum / count, -0.33, 0.33);
            }
        }
        core::mem::swap(&mut current, &mut next);
    }
    current
}

// ── Absorption ───────────────────────────────────────────────────────

/// Beam attenuation through a ternary medium along a row.
/// +1 cells absorb (reduce intensity), -1 cells amplify, 0 cells are neutral.
pub fn absorption(grid: &[i8], w: usize, h: usize, row: usize, initial: f64, coeff: f64) -> Vec<f64> {
    let mut intensity = initial;
    let mut result = Vec::new();
    if row >= h { return result; }
    for x in 0..w {
        let cell = grid[idx(x, row, w)];
        let factor = match cell {
            1 => (-coeff).exp(),   // absorb
            -1 => coeff.exp(),     // amplify
            _ => 1.0,              // neutral
        };
        intensity *= factor;
        result.push(intensity);
    }
    result
}

// ── ShadowCast ───────────────────────────────────────────────────────

/// Directional shadow casting. +1 cells cast shadows in direction (dx, dy).
/// Shadow cells are marked -1, lit cells +1, shadow boundaries 0.
pub fn shadow_cast(grid: &[i8], w: usize, h: usize, dx: i32, dy: i32) -> Vec<i8> {
    let mut shadow = alloc::vec![1i8; w * h];

    for y in 0..h {
        for x in 0..w {
            if grid[idx(x, y, w)] == 1 {
                // Cast shadow in direction
                let mut sx = x as i32 + dx;
                let mut sy = y as i32 + dy;
                while in_bounds(sx, sy, w, h) {
                    let si = idx(sx as usize, sy as usize, w);
                    if shadow[si] == 1 {
                        shadow[si] = -1;
                    }
                    sx += dx;
                    sy += dy;
                }
            }
        }
    }

    // Mark cells that are shadow boundaries (0: adjacent to both +1 and -1)
    let mut result = shadow.clone();
    for y in 0..h {
        for x in 0..w {
            if shadow[idx(x, y, w)] == -1 {
                let ns = neighbors4(x, y, w, h);
                for &(nx, ny) in &ns {
                    if shadow[idx(nx, ny, w)] == 1 {
                        result[idx(x, y, w)] = 0;
                        break;
                    }
                }
            }
        }
    }
    result
}

// ── HalfLife ─────────────────────────────────────────────────────────

/// Simulate radioactive decay: +1 → 0 → -1 with given probability per step.
pub fn half_life(grid: &[i8], w: usize, h: usize, steps: usize, decay_prob: f64) -> Vec<i8> {
    let mut current = grid.to_vec();
    let mut next = current.clone();

    for _ in 0..steps {
        for y in 0..h {
            for x in 0..w {
                let i = idx(x, y, w);
                next[i] = match current[i] {
                    1 => if pseudo_random(x, y, steps) < decay_prob { 0 } else { 1 },
                    0 => if pseudo_random(x + 100, y + 100, steps) < decay_prob { -1 } else { 0 },
                    _ => current[i],
                };
            }
        }
        core::mem::swap(&mut current, &mut next);
    }
    current
}

/// Simple deterministic pseudo-random for reproducible tests (not crypto).
fn pseudo_random(x: usize, y: usize, seed: usize) -> f64 {
    let v = ((x.wrapping_mul(2654435761)).wrapping_add(y.wrapping_mul(2246822519)).wrapping_add(seed.wrapping_mul(3266489917))) as u64;
    ((v ^ (v >> 16)) & 0xFFFF) as f64 / 65535.0
}

// ── CascadeEvent ─────────────────────────────────────────────────────

/// A cascade: one cell flips, neighbors may flip with given probability.
/// Returns the set of flipped cell indices and final grid.
pub fn cascade(grid: &[i8], w: usize, h: usize, start_x: usize, start_y: usize, flip_prob: f64) -> (Vec<usize>, Vec<i8>) {
    let mut current = grid.to_vec();
    let mut flipped = Vec::new();
    let mut frontier = alloc::vec![(start_x, start_y)];

    // Flip the start cell
    let si = idx(start_x, start_y, w);
    current[si] = flip_val(current[si]);
    flipped.push(si);

    let mut step = 0usize;
    while !frontier.is_empty() {
        let mut next_frontier = Vec::new();
        for &(fx, fy) in &frontier {
            let ns = neighbors4(fx, fy, w, h);
            for &(nx, ny) in &ns {
                let ni = idx(nx, ny, w);
                if !flipped.contains(&ni) && pseudo_random(nx, ny, step) < flip_prob {
                    current[ni] = flip_val(current[ni]);
                    flipped.push(ni);
                    next_frontier.push((nx, ny));
                }
            }
        }
        frontier = next_frontier;
        step += 1;
    }

    (flipped, current)
}

fn flip_val(v: i8) -> i8 {
    match v {
        1 => -1,
        -1 => 1,
        _ => 0,
    }
}

// ── IrradianceProfile ────────────────────────────────────────────────

/// 1D irradiance profile across a row or column.
pub fn profile_row(field: &[f64], w: usize, h: usize, row: usize) -> Vec<f64> {
    if row >= h { return Vec::new(); }
    (0..w).map(|x| field[idx(x, row, w)]).collect()
}

pub fn profile_col(field: &[f64], w: usize, h: usize, col: usize) -> Vec<f64> {
    if col >= w { return Vec::new(); }
    (0..h).map(|y| field[idx(col, y, w)]).collect()
}

// ── Shielding ────────────────────────────────────────────────────────

/// Compute how many 0 cells (shields) are between two points on the grid.
/// Uses Bresenham-style line traversal.
pub fn shielding(grid: &[i8], w: usize, h: usize, x0: usize, y0: usize, x1: usize, y1: usize) -> usize {
    let dx = (x1 as i64 - x0 as i64).abs();
    let dy = (y1 as i64 - y0 as i64).abs();
    let sx: i64 = if x0 < x1 { 1 } else { -1 };
    let sy: i64 = if y0 < y1 { 1 } else { -1 };

    let mut x = x0 as i64;
    let mut y = y0 as i64;
    let mut err = dx - dy;
    let mut count = 0;

    loop {
        if (x as usize) < w && (y as usize) < h {
            if grid[idx(x as usize, y as usize, w)] == 0 {
                count += 1;
            }
        }
        if x == x1 as i64 && y == y1 as i64 { break; }
        let e2 = 2 * err;
        if e2 > -dy { err -= dy; x += sx; }
        if e2 < dx { err += dx; y += sy; }
    }
    count
}

/// Check if there is a clear line of sight (no 0 shields) between two points.
pub fn has_line_of_sight(grid: &[i8], w: usize, h: usize, x0: usize, y0: usize, x1: usize, y1: usize) -> bool {
    shielding(grid, w, h, x0, y0, x1, y1) == 0
}

// ── Tests ────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn grid3x3(vals: &[i8]) -> (Vec<i8>, usize, usize) {
        (vals.to_vec(), 3, 3)
    }

    #[test]
    fn source_irradiance_at_center() {
        let src = RadiationSource::new(1, 1, 10.0, 1.0);
        let v = src.irradiance_at(1, 1);
        assert_eq!(v, 10.0); // distance 0
    }

    #[test]
    fn source_irradiance_decay() {
        let src = RadiationSource::new(0, 0, 10.0, 1.0);
        let v0 = src.irradiance_at(0, 0);
        let v1 = src.irradiance_at(1, 0);
        let v2 = src.irradiance_at(2, 0);
        assert!(v0 > v1);
        assert!(v1 > v2);
    }

    #[test]
    fn compute_irradiance_basic() {
        let sources = vec![RadiationSource::new(0, 0, 5.0, 0.0)];
        let field = compute_irradiance(&sources, 3, 3);
        assert!(field[idx(0, 0, 3)] > field[idx(2, 2, 3)]);
    }

    #[test]
    fn to_ternary_thresholds() {
        let field = vec![-2.0, -0.5, 0.0, 0.5, 2.0];
        let t = to_ternary(&field, 5, 1, -0.33, 0.33);
        assert_eq!(t, vec![-1, -1, 0, 1, 1]);
    }

    #[test]
    fn diffuse_uniform() {
        let grid = vec![1, 1, 1, 1, 1, 1, 1, 1, 1];
        let result = diffuse(&grid, 3, 3, 1);
        assert!(result.iter().all(|&v| v == 1)); // uniform stays uniform
    }

    #[test]
    fn diffuse_spreads() {
        let mut grid = vec![0; 25];
        grid[12] = 1; // center of 5x5
        let result = diffuse(&grid, 5, 5, 5);
        // After diffusion, neighbors should have some non-zero values
        let non_zero = result.iter().filter(|&&v| v != 0).count();
        assert!(non_zero > 1); // spread beyond center
    }

    #[test]
    fn absorption_row() {
        let grid = vec![1, 0, -1, 0, 1];
        let result = absorption(&grid, 5, 1, 0, 100.0, 0.5);
        assert_eq!(result.len(), 5);
        assert!(result[0] < 100.0); // +1 absorbs
        assert_eq!(result[1], result[0]); // 0 neutral
        assert!(result[2] > result[1]); // -1 amplifies
    }

    #[test]
    fn shadow_cast_basic() {
        // 3x3 with a wall at (1,1)
        let grid = vec![0, 0, 0, 0, 1, 0, 0, 0, 0];
        let shadow = shadow_cast(&grid, 3, 3, 1, 1); // direction: down-right
        assert_eq!(shadow[idx(2, 2, 3)], -1); // shadow behind wall
    }

    #[test]
    fn shadow_cast_origin_lit() {
        let grid = vec![0, 0, 0, 0, 1, 0, 0, 0, 0];
        let shadow = shadow_cast(&grid, 3, 3, 0, 1); // direction: down
        assert_eq!(shadow[idx(1, 1, 3)], 1); // the wall itself is lit
    }

    #[test]
    fn half_life_decay() {
        let grid = vec![1; 25];
        let result = half_life(&grid, 5, 5, 10, 1.0); // guaranteed decay
        // With prob=1.0 all should decay
        assert!(result.iter().any(|&v| v != 1));
    }

    #[test]
    fn half_life_no_decay() {
        let grid = vec![1, 0, -1, 1, 0];
        let result = half_life(&grid, 5, 1, 5, 0.0); // prob=0, no decay
        assert_eq!(result, grid);
    }

    #[test]
    fn cascade_single_flip() {
        let grid = vec![0; 25]; // all zeros, nothing to cascade to
        let (flipped, _) = cascade(&grid, 5, 5, 2, 2, 0.0);
        assert_eq!(flipped.len(), 1); // only the start cell
    }

    #[test]
    fn cascade_propagation() {
        let grid = vec![1; 25]; // all +1
        let (flipped, result) = cascade(&grid, 5, 5, 2, 2, 1.0); // prob=1
        assert!(flipped.len() > 1);
        assert!(result.iter().any(|&v| v == -1)); // flipped to -1
    }

    #[test]
    fn profile_row_basic() {
        let field = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let p = profile_row(&field, 3, 2, 0);
        assert_eq!(p, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn profile_col_basic() {
        let field = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let p = profile_col(&field, 3, 2, 1);
        assert_eq!(p, vec![2.0, 5.0]);
    }

    #[test]
    fn shielding_line_of_sight() {
        let grid = vec![1, 0, 1, 0, 0, 0, 1, 0, 1]; // 3x3
        let s = shielding(&grid, 3, 3, 0, 0, 2, 2);
        assert!(s > 0); // has shields in between
    }

    #[test]
    fn has_los_clear() {
        let grid = vec![1, 1, 1, 1, 1, 1, 1, 1, 1];
        assert!(has_line_of_sight(&grid, 3, 3, 0, 0, 2, 2));
    }

    #[test]
    fn has_los_blocked() {
        let mut grid = vec![1; 9];
        grid[idx(1, 1, 3)] = 0; // shield in center
        assert!(!has_line_of_sight(&grid, 3, 3, 0, 0, 2, 2));
    }

    #[test]
    fn multiple_sources_additive() {
        let sources = vec![
            RadiationSource::new(0, 0, 5.0, 0.5),
            RadiationSource::new(4, 4, 5.0, 0.5),
        ];
        let field = compute_irradiance(&sources, 5, 5);
        // Center should have contributions from both
        let center = field[idx(2, 2, 5)];
        assert!(center > 0.0);
    }
}
