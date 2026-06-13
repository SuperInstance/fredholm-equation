# fredholm-equation: Nyström Method Solver for Fredholm Integral Equations

Solves **Fredholm integral equations of the second kind**:

$$\varphi(x) = f(x) + \lambda \int_a^b K(x,t)\,\varphi(t)\,dt$$

using the **Nyström method** with trapezoidal quadrature and Gaussian elimination. Given a kernel K(x,t), forcing function f(x), and eigenvalue λ, it computes the solution φ(x) at N quadrature nodes.

## Why It Matters

Fredholm equations arise in quantum scattering theory, radiative transfer, signal restoration, electrostatics, and continuum mechanics. They are the canonical form for any linear inverse problem where the relationship between input and output is a convolution-like integral. The Nyström method is the simplest numerical approach — it discretizes the integral into a linear system and solves exactly, making it ideal for prototyping and teaching before moving to Galerkin or collocation methods.

## How It Works

### Discretization (Trapezoidal Rule)

Partition [a, b] into N equally spaced nodes with spacing h = (b−a)/(N−1). The trapezoidal rule gives quadrature weights:

```
w₀ = w_{N-1} = h/2
w_j = h       for j = 1, ..., N-2
```

Substituting into the integral equation at each node xᵢ:

```
φ(xᵢ) = f(xᵢ) + λ Σⱼ wⱼ · K(xᵢ, xⱼ) · φ(xⱼ)
```

### Linear System

This yields the N×N system **(I − λ·D·K)φ = f** where:

- **I** is the identity matrix
- **D = diag(w₀, w₁, ..., w_{N-1})** is the weight matrix
- **K** is the kernel matrix with K_{ij} = K(xᵢ, xⱼ)
- **φ** is the unknown solution vector
- **f** is the forcing vector

### Gaussian Elimination with Partial Pivoting

The solver uses partial pivoting (selecting the row with maximum |pivot|) for numerical stability:

```
for col in 0..N:
    find pivot row = argmax|aug[row][col]|
    swap rows
    eliminate below
back-substitute
```

**Complexity**: O(N³) for the elimination, O(N²) for back-substitution. Memory: O(N²) for the augmented matrix.

### Correctness Check

When λ = 0, the integral term vanishes and φ = f exactly. The test suite verifies this with f(x) = x² on [0, 1] to 10⁻⁶ precision.

### Error Analysis

The trapezoidal Nyström method has error O(h²) for smooth kernels and solutions. For N = 21 nodes on [0, 1]: h ≈ 0.05, giving error ~2.5 × 10⁻³. Increasing N improves this quadratically.

| Component | Time | Space |
|-----------|------|-------|
| Matrix assembly | O(N²) | O(N²) |
| Gaussian elimination | O(N³) | O(N²) |
| Back substitution | O(N²) | O(N) |
| Total | O(N³) | O(N²) |

## Quick Start

```rust
use fredholm_equation::FredholmSolver;

let solver = FredholmSolver::new(21, 0.0, 1.0);

// Solve φ(x) = f(x) + λ∫K(x,t)φ(t)dt
let phi = solver.solve(
    0.5,                           // λ (eigenvalue)
    |_x, t| (x * t).exp(),        // K(x,t) = e^{xt}
    |x| x * x,                     // f(x) = x²
);

let nodes = solver.nodes();
for (i, &x) in nodes.iter().enumerate() {
    println!("φ({:.2}) = {:.6}", x, phi[i]);
}
```

## API

### `FredholmSolver`

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(n_points: usize, a: f64, b: f64) -> Self` | Create solver on [a, b] with N nodes |
| `solve` | `<F, G>(&self, lambda: f64, kernel: F, f: G) -> Vec<f64>` | Solve for φ, given K and f as closures |
| `nodes` | `(&self) -> Vec<f64>` | Return the quadrature node positions |

The closures `F: Fn(f64, f64) -> f64` (kernel) and `G: Fn(f64) -> f64` (forcing) make the solver generic over any problem.

## Architecture Notes

This is a **γ (gamma)** module — pure mathematics, deterministic, no I/O. In the γ + η = C framework, it provides the numerical analysis foundation. An **η** layer could wrap it with adaptive quadrature, multi-resolution refinement, or eigenvalue sweeping (finding all λ for which the equation is solvable — the Fredholm alternative).

## References

- Atkinson, K. E. (1997). *The Numerical Solution of Integral Equations of the Second Kind*. Cambridge University Press.
- Delves, L. M. & Mohamed, J. L. (1985). *Computational Methods for Integral Equations*. Cambridge.
- Kress, R. (2014). *Linear Integral Equations* (3rd ed.). Springer.

## License

MIT
