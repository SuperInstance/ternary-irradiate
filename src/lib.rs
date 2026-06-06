//! Ternary irradiation: radiation damage, cascade simulation, annealing, defect tracking.

use std::collections::HashMap;

/// Ternary lattice for radiation simulation
#[derive(Clone, Debug)]
pub struct TernaryLattice {
    size: usize,
    /// -1 = damaged, 0 = vacant, +1 = intact
    cells: Vec<i8>,
    defect_log: Vec<DefectEvent>,
}

#[derive(Clone, Debug)]
pub struct DefectEvent {
    pub tick: usize,
    pub position: usize,
    pub old_state: i8,
    pub new_state: i8,
    pub cause: DefectCause,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DefectCause {
    PrimaryKnock,    // Direct radiation hit
    Cascade,         // Secondary cascade
    ThermalRecovery, // Annealing
    Recombination,   // Defect-vacancy recombination
    Spontaneous,     // Spontaneous recovery
}

impl TernaryLattice {
    pub fn new(size: usize) -> Self {
        Self { size, cells: vec![1; size], defect_log: Vec::new() }
    }

    pub fn with_defects(size: usize, intact_ratio: f64, rng_vals: &[f64]) -> Self {
        let mut cells = vec![1i8; size];
        let defect_count = (size as f64 * (1.0 - intact_ratio)) as usize;
        for i in 0..defect_count.min(rng_vals.len()) {
            let idx = (rng_vals[i] * size as f64) as usize % size;
            cells[idx] = -1;
        }
        Self { size, cells, defect_log: Vec::new() }
    }

    pub fn get(&self, idx: usize) -> i8 { self.cells[idx] }
    pub fn set(&mut self, idx: usize, v: i8) { self.cells[idx] = v; }
    pub fn len(&self) -> usize { self.size }

    /// Simulate a radiation hit: primary knock-on + cascade
    pub fn irradiate(&mut self, hit_pos: usize, energy: f64, tick: usize, rng: &[f64]) -> usize {
        let mut defects = 0;
        // Primary knock
        if self.cells[hit_pos] == 1 {
            self.defect_log.push(DefectEvent { tick, position: hit_pos, old_state: 1, new_state: -1, cause: DefectCause::PrimaryKnock });
            self.cells[hit_pos] = -1;
            defects += 1;
        }

        // Cascade: each unit of energy can create additional defects
        let cascade_size = (energy * 3.0) as usize;
        let mut rng_idx = 0;
        for _ in 0..cascade_size {
            if rng_idx + 1 >= rng.len() { break; }
            // Random walk from hit position
            let offset = (rng[rng_idx] * 10.0) as i32 - 5;
            rng_idx += 1;
            let target = ((hit_pos as i32 + offset).rem_euclid(self.size as i32)) as usize;
            if self.cells[target] == 1 {
                self.defect_log.push(DefectEvent { tick, position: target, old_state: 1, new_state: -1, cause: DefectCause::Cascade });
                self.cells[target] = -1;
                defects += 1;
            }
        }
        defects
    }

    /// Thermal annealing: damaged cells spontaneously recover with given probability
    pub fn anneal(&mut self, temperature: f64, tick: usize, rng: &[f64]) -> usize {
        let recovery_prob = (temperature / 10.0).min(0.95);
        let mut recovered = 0;
        for i in 0..self.size {
            if self.cells[i] == -1 {
                let r = rng[i % rng.len()];
                if r < recovery_prob {
                    self.defect_log.push(DefectEvent { tick, position: i, old_state: -1, new_state: 1, cause: DefectCause::ThermalRecovery });
                    self.cells[i] = 1;
                    recovered += 1;
                }
            }
        }
        recovered
    }

    /// Recombination: damaged cell next to vacancy can fill it
    pub fn recombine(&mut self, tick: usize) -> usize {
        let mut recombined = 0;
        let mut changes = Vec::new();
        for i in 0..self.size {
            if self.cells[i] == 0 {
                // Check neighbors for damaged cells
                for &neighbor in &[i.saturating_sub(1), (i+1).min(self.size-1)] {
                    if neighbor < self.size && self.cells[neighbor] == -1 {
                        changes.push((i, 1)); // vacancy filled
                        changes.push((neighbor, 0)); // damaged becomes vacant (moved)
                        break;
                    }
                }
            }
        }
        for (pos, new_state) in changes {
            let old = self.cells[pos];
            if old != new_state {
                self.defect_log.push(DefectEvent { tick, position: pos, old_state: old, new_state, cause: DefectCause::Recombination });
                self.cells[pos] = new_state;
                recombined += 1;
            }
        }
        recombined
    }

    /// Integrity: fraction of intact cells
    pub fn integrity(&self) -> f64 {
        self.cells.iter().filter(|&&v| v == 1).count() as f64 / self.size as f64
    }

    /// Defect density
    pub fn defect_density(&self) -> f64 {
        self.cells.iter().filter(|&&v| v == -1).count() as f64 / self.size as f64
    }

    /// Vacancy density
    pub fn vacancy_density(&self) -> f64 {
        self.cells.iter().filter(|&&v| v == 0).count() as f64 / self.size as f64
    }

    /// Statistics
    pub fn stats(&self) -> LatticeStats {
        let intact = self.cells.iter().filter(|&&v| v == 1).count();
        let damaged = self.cells.iter().filter(|&&v| v == -1).count();
        let vacant = self.cells.iter().filter(|&&v| v == 0).count();
        LatticeStats { intact, damaged, vacant, total: self.size }
    }
}

#[derive(Debug)]
pub struct LatticeStats {
    pub intact: usize,
    pub damaged: usize,
    pub vacant: usize,
    pub total: usize,
}

/// Full irradiation simulation
pub struct IrradiationSim {
    pub lattice: TernaryLattice,
    pub tick: usize,
    pub history: Vec<LatticeStats>,
}

impl IrradiationSim {
    pub fn new(lattice_size: usize) -> Self {
        Self { lattice: TernaryLattice::new(lattice_size), tick: 0, history: Vec::new() }
    }

    /// Run one tick: irradiate + anneal + recombine
    pub fn step(&mut self, dose_rate: f64, temperature: f64, rng: &[f64]) -> (usize, usize, usize) {
        let n_hits = (dose_rate * self.lattice.size as f64) as usize;
        let mut total_defects = 0;
        let mut offset = 0;
        for _ in 0..n_hits {
            if offset + 10 >= rng.len() { break; }
            let pos = (rng[offset] * self.lattice.size as f64) as usize;
            let energy = rng[offset + 1] * 5.0 + 1.0;
            offset += 2;
            total_defects += self.lattice.irradiate(pos, energy, self.tick, &rng[offset..]);
            offset += (energy * 3.0) as usize;
        }
        let recovered = self.lattice.anneal(temperature, self.tick, rng);
        let recombined = self.lattice.recombine(self.tick);
        self.history.push(self.lattice.stats());
        self.tick += 1;
        (total_defects, recovered, recombined)
    }

    /// Run simulation for n ticks
    pub fn run(&mut self, ticks: usize, dose_rate: f64, temperature: f64, rng: &[f64]) {
        let mut rng_offset = 0;
        for _ in 0..ticks {
            let needed = (dose_rate * self.lattice.size as f64 * 15.0) as usize;
            if rng_offset + needed >= rng.len() { break; }
            self.step(dose_rate, temperature, &rng[rng_offset..]);
            rng_offset += needed;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_rng(n: usize) -> Vec<f64> {
        (0..n).map(|i| ((i as u64 * 6364136223846793005 + 1442695040888963407) as f64 / u64::MAX as f64)).collect()
    }

    #[test]
    fn test_fresh_lattice_full_integrity() {
        let l = TernaryLattice::new(100);
        assert!((l.integrity() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_irradiate_creates_defects() {
        let mut l = TernaryLattice::new(100);
        let rng = make_rng(200);
        let defects = l.irradiate(50, 2.0, 0, &rng);
        assert!(defects > 0);
        assert!(l.integrity() < 1.0);
    }

    #[test]
    fn test_anneal_recovers() {
        let mut l = TernaryLattice::new(100);
        let rng = make_rng(200);
        l.irradiate(50, 3.0, 0, &rng);
        let damaged_before = l.defect_density();
        l.anneal(5.0, 1, &rng);
        let damaged_after = l.defect_density();
        assert!(damaged_after <= damaged_before);
    }

    #[test]
    fn test_recombination() {
        let mut l = TernaryLattice::new(20);
        l.set(5, -1); // damaged
        l.set(6, 0);  // vacancy
        let r = l.recombine(0);
        assert!(r > 0);
    }

    #[test]
    fn test_simulation_runs() {
        let rng = make_rng(50000);
        let mut sim = IrradiationSim::new(200);
        sim.run(50, 0.01, 2.0, &rng);
        assert_eq!(sim.tick, 50);
        assert!(!sim.history.is_empty());
    }

    #[test]
    fn test_equilibrium_dose_recovery() {
        let rng = make_rng(100000);
        let mut sim = IrradiationSim::new(500);
        sim.run(200, 0.005, 8.0, &rng); // high temperature, low dose
        // Should recover most defects
        assert!(sim.lattice.integrity() > 0.5);
    }

    #[test]
    fn test_defect_log_tracks_events() {
        let mut l = TernaryLattice::new(50);
        let rng = make_rng(200);
        l.irradiate(25, 1.0, 0, &rng);
        assert!(!l.defect_log.is_empty());
        assert_eq!(l.defect_log[0].cause, DefectCause::PrimaryKnock);
    }

    #[test]
    fn test_stats_consistency() {
        let mut l = TernaryLattice::new(100);
        let rng = make_rng(200);
        l.irradiate(50, 2.0, 0, &rng);
        let s = l.stats();
        assert_eq!(s.intact + s.damaged + s.vacant, s.total);
    }
}
