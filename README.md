# Fredholm Integral Equation Solver

**A Rust library for solving Fredholm integral equations of the second kind** using the Nyström method with trapezoidal quadrature and Gaussian elimination.

## Why It Matters

Fredholm integral equations of the second kind, `φ(x) = f(x) + λ∫K(x,t)φ(t)dt`, appear throughout physics and engineering: quantum scattering theory, electrostatics boundary-value problems, radiative heat transfer, and signal deconvolution. The Nyström method discretizes the integral using numerical quadrature, converting the integral equation into a standard linear system `(I - λDK)φ = f` that can be solved with Gaussian elimination. This library implements the full pipeline — from kernel and source function closures to solution vector — in pure Rust with no dependencies.

## How It Works

The solver takes `n` equally-spaced quadrature nodes on `[a, b]` with spacing `h = (b-a)/(n-1)`. Using the trapezoidal rule, the quadrature weights are `h/2` at the endpoints and `h` at interior points. This converts the integral equation into the matrix equation:

```
(I - λ·D·K) φ = f
```

where `D = diag(w₀, w₁, ..., wₙ₋₁)` contains the quadrature weights and `Kᵢⱼ = kernel(xᵢ, xⱼ)`. The resulting `n×n` dense linear system is solved via partial-pivot Gaussian elimination with **O(n³)** complexity. The implementation handles the full pipeline: node generation, matrix assembly, pivoting, forward elimination, and back substitution.

## Quick Start

```rust
use fredholm_equation::FredholmSolver;

fn main() {
    // Solve φ(x) = f(x) + 0.5·∫₀¹ K(x,t)·φ(t) dt
    // with K(x,t) = x·t and f(x) = 1
    let solver = FredholmSolver::new(51, 0.0, 1.0);

    let phi = solver.solve(
        0.5,                            // λ
        |x, t| x * t,                   // kernel K(x,t)
        |_x| 1.0,                       // f(x)
    );

    let nodes = solver.nodes();
    for (i, &x) in nodes.iter().enumerate() {
        println!("φ({:.3}) = {:.6}", x, phi[i]);
    }
}
```

## API

| Type / Method | Complexity | Description |
|---|---|---|
| `FredholmSolver::new(n, a, b)` | **O(1)** | Create solver with `n` quadrature nodes on `[a, b]` |
| `solve(lambda, kernel, f)` | **O(n³)** | Solve via Nyström + Gaussian elimination |
| `nodes()` | **O(n)** | Return the quadrature node positions |

## Architecture Notes

Part of the SuperInstance numerical methods collection. Companion crates include `gauss-markov` (stochastic processes) and `hermite-polynomial` (quadrature). See the [Architecture Guide](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md).

## License

MIT
