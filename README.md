# ternary-irradiate

**Radiation and energy propagation on ternary grids. Point sources, inverse-square falloff, diffusion, shadows, and half-life decay.**

How does energy spread from a point source across a grid? Inverse-square law says the intensity falls off as 1/r². But on a discrete ternary grid, the continuous function is quantized: at each cell, the irradiance is snapped to {-1, 0, +1} based on configurable thresholds. Energy propagates outward, diffuses into neighbors, decays over time with configurable half-life, and is blocked by shadow-casting obstacles.

This crate models all of that: point sources, field computation, diffusion, shadow casting, cascade events (chain reactions where irradiated cells become sources), and shielding effectiveness.

## What's Inside

- **`irradiance_field(sources, width, height)`** — compute irradiance at every cell from point sources using inverse-square law, snapped to ternary
- **`diffuse(grid, width, rate)`** — one step of diffusion: each cell blends with its neighbors
- **`shadow_cast(sources, obstacles, width, height)`** — compute shadow regions behind obstacles
- **`half_life_decay(grid, half_life)`** — exponential decay: intensity halves every `half_life` ticks
- **`cascade(grid, width, threshold)`** — cells above threshold become new sources (chain reaction)
- **`shield_effectiveness(grid, width, shield_value)`** — how well does a shield value block propagation?

## Quick Example

```rust
use ternary_irradiate::*;

// Single point source at center of 20x20 grid
let sources = vec![(10, 10)];
let field = irradiance_field(&sources, 20, 20);
// Center is +1, nearby cells are 0, distant cells are -1 (below threshold)

// Add obstacles that cast shadows
let obstacles = vec![(10, 5), (10, 6), (10, 7)]; // wall above source
let shadow = shadow_cast(&sources, &obstacles, 20, 20);
// Cells behind the wall (rows 0-4) are in shadow

// Cascade: irradiated cells become sources
let cascade_field = cascade(&field, 20, 1);
// Any cell at +1 becomes a source, irradiating its neighbors

// Decay over time
let decayed = half_life_decay(&field, 10);
// After 10 ticks, intensity has halved
```

## The Deeper Truth

**Ternary irradiance is quantized energy.** The continuous irradiance field (0 to ∞) is compressed into three states: below threshold (-1), within range (0), above threshold (+1). This quantization loses the *magnitude* of the energy but preserves the *topology* of the field — where the peaks are, where the valleys are, where the shadow boundaries fall. For many applications (coverage planning, safety zones, sensor placement), the topology is what matters, not the exact numbers.

The cascade function is where ternary irradiance gets genuinely interesting: a cell that receives enough energy (+1) becomes a new source, which can trigger its neighbors, which trigger *their* neighbors. This is a chain reaction — and the ternary threshold determines whether it spreads or dies out. Set the threshold too high and nothing cascades. Set it too low and everything lights up. The sweet spot creates branching, tree-like propagation patterns that look exactly like lightning, neural firing, and epidemic spreading.

## `#![no_std]`

This crate runs without an allocator. Use it on microcontrollers, embedded systems, or anywhere you need radiation modeling without an OS.

**Use cases:**
- **Wireless signal propagation** — model coverage areas with ternary signal strength
- **Environmental monitoring** — radiation or pollution spreading with ternary severity levels
- **Game AI** — light/shadow/line-of-sight on ternary terrain
- **Sensor networks** — coverage planning and gap detection
- **Epidemic modeling** — cascade dynamics for ternary infection states

## See Also

- **ternary-field** — static field analysis (no propagation dynamics)
- **ternary-fire** — fire spreading (a cascade with specific rules)
- **ternary-shield** — containment analysis (how effective are obstacles at blocking?)
- **ternary-diffusion** — (if it exists) continuous diffusion without the ternary snap
- **ternary-sandpile** — another cascade model (toppling vs. irradiation)

## Install

```bash
cargo add ternary-irradiate
```

## License

MIT
