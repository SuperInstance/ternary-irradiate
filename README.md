# ternary-irradiate

**Ternary radiation damage simulation: primary knock-on events, cascade propagation, thermal annealing, defect-vacancy recombination, and full dose-response modeling.**

## Background

Radiation damage in materials is a cascade process: a high-energy particle (neutron, ion, electron) strikes an atom, displacing it from its lattice site. The displaced atom — a primary knock-on atom (PKA) — carries kinetic energy and can displace further atoms, creating a damage cascade. The resulting lattice has three types of sites: intact (+1), vacant (0, where an atom was ejected), and damaged/interstitial (−1, where an atom came to rest in a non-lattice position).

Over time, thermal energy allows some defects to anneal (spontaneous recovery) and interstitials to recombine with vacancies (Frenkel pair recombination). The competition between damage creation and recovery determines the material's steady-state defect concentration.

`ternary-irradiate` models this process using a ternary lattice where each cell is in one of three states:

| State | Value | Physical meaning |
|-------|-------|------------------|
| Intact | +1 | Normal lattice site |
| Vacant | 0 | Empty lattice site (vacancy) |
| Damaged | −1 | Interstitial or displaced atom |

The simulation tracks the full lifecycle: radiation hits → cascade damage → annealing → recombination → equilibrium.

## How It Works

### Ternary Lattice

```rust
let lattice = TernaryLattice::new(10_000);  // All intact (+1)
let lattice = TernaryLattice::with_defects(10_000, 0.95, &rng); // 95% intact
```

### Radiation Events

`irradiate(hit_pos, energy, tick, rng)` simulates a single radiation event:

1. **Primary knock-on**: The atom at `hit_pos` is displaced (intact → damaged).
2. **Cascade**: `energy × 3` additional atoms near the hit site are randomly displaced using a random-walk offset from the impact point. Each displacement is logged as `DefectCause::Cascade`.

The energy parameter controls cascade size. Higher energy → more secondary defects.

### Recovery Mechanisms

- **Thermal annealing** (`anneal(temperature, tick, rng)`): Each damaged site recovers with probability `min(temperature / 10, 0.95)`. Higher temperature → faster recovery. Logged as `DefectCause::ThermalRecovery`.
- **Recombination** (`recombine(tick)`): A vacancy (0) adjacent to a damaged site (−1) can "fill" the vacancy — the damaged atom moves into the vacancy. The damaged site becomes vacant; the vacancy becomes intact. Logged as `DefectCause::Recombination`.

### Full Simulation

`IrradiationSim` runs a time-stepped simulation:

```rust
let mut sim = IrradiationSim::new(lattice_size);
sim.run(ticks, dose_rate, temperature, &rng);
```

Each tick:
1. Apply `n = dose_rate × lattice_size` radiation hits with random positions and energies.
2. Run thermal annealing at the given temperature.
3. Run defect-vacancy recombination.
4. Record lattice statistics (intact, damaged, vacant counts).

### Defect Logging

Every state change is recorded as a `DefectEvent`:
```rust
pub struct DefectEvent {
    tick: usize,
    position: usize,
    old_state: i8,
    new_state: i8,
    cause: DefectCause,  // PrimaryKnock | Cascade | ThermalRecovery | Recombination | Spontaneous
}
```

### Statistics

```rust
lattice.integrity()       // Fraction of intact cells
lattice.defect_density()  // Fraction of damaged cells
lattice.vacancy_density() // Fraction of vacant cells
lattice.stats()           // LatticeStats { intact, damaged, vacant, total }
```

## Experimental Results

The test suite verifies (using a deterministic LCG random number generator):
- **Fresh lattice**: 100% integrity.
- **Irradiation creates defects**: Integrity drops below 1.0 after irradiation.
- **Annealing recovers defects**: Defect density decreases after annealing.
- **Recombination works**: Adjacent damage-vacancy pairs are resolved.
- **Full simulation**: 50-tick simulation runs to completion with recorded history.
- **Equilibrium**: At high temperature (8.0) and low dose rate (0.005), integrity recovers above 50%.
- **Defect logging**: Events are recorded with correct causes (PrimaryKnock for first hit, Cascade for secondary).
- **Conservation**: `intact + damaged + vacant = total` at all times.

## Impact

This crate demonstrates that ternary state systems are a natural fit for radiation damage modeling. The three states map directly to the physical picture, and the ternary morphological operations from `ternary-morph` can be applied to the resulting lattices for spatial analysis. The defect logging system provides full provenance tracking — every state change has a timestamp, position, and cause.

## Use Cases

1. **Nuclear Materials Research** — Simulate radiation damage in reactor materials (steel, silicon carbide, tungsten) under various dose rates and temperatures. Tune parameters to match experimental data.
2. **Space Electronics** — Model single-event upsets (SEUs) in ternary logic circuits exposed to cosmic radiation. The cascade model captures multi-bit upsets from a single particle.
3. **Ion Implantation** — Simulate focused ion beam processing where controlled damage is desired (e.g., semiconductor doping). The cascade model predicts damage spread.
4. **Defect Engineering** — Explore annealing schedules (temperature profiles over time) to minimize defect density while controlling vacancy concentration.

## Open Questions

1. **3D extension** — The current model is 1D (linear lattice). Real crystal lattices are 3D. Would the same cascade model extend to `TernaryGrid` (2D) or a 3D lattice?
2. **Cascade energy distribution** — The current model uses uniform random energy. Real PKAs follow a distribution (e.g., Kinchin-Pease model). Would a more realistic distribution change the equilibrium behavior?
3. **Defect clustering** — Do damaged sites cluster spatially under high dose rates? Morphological analysis with `ternary-morph` could quantify this.

## Connection to Oxide Stack

`ternary-irradiate` is the physics simulation layer of the ternary fleet. It consumes `TernaryGrid` concepts from `ternary-core`, produces lattices suitable for analysis with `ternary-morph` (spatial damage patterns) and `ternary-signals` (spectral analysis of damage propagation), and uses `ternary-walk` concepts for cascade random walks. The simulation demonstrates the practical value of ternary state systems: when your domain has three natural states, ternary is not an approximation — it's the correct representation.
