/// Solve Fredholm integral equation of the second kind:
/// φ(x) = f(x) + λ ∫ K(x,t)φ(t) dt  via Nyström method (trapezoidal rule)
pub struct FredholmSolver {
    pub n_points: usize,
    pub a: f64,
    pub b: f64,
}

impl FredholmSolver {
    pub fn new(n_points: usize, a: f64, b: f64) -> Self {
        Self { n_points, a, b }
    }

    /// Solve using trapezoidal quadrature
    /// kernel(x, t) and f(x) provided as closures
    pub fn solve<F, G>(&self, lambda: f64, kernel: F, f: G) -> Vec<f64>
    where
        F: Fn(f64, f64) -> f64,
        G: Fn(f64) -> f64,
    {
        let n = self.n_points;
        let h = (self.b - self.a) / (n - 1) as f64;
        let nodes: Vec<f64> = (0..n).map(|i| self.a + i as f64 * h).collect();

        // Build the system (I - λ*D*K) φ = f
        let mut mat = vec![vec![0.0; n]; n];
        let mut rhs = vec![0.0; n];
        for i in 0..n {
            rhs[i] = f(nodes[i]);
            mat[i][i] = 1.0;
            for j in 0..n {
                let w = if j == 0 || j == n - 1 { h / 2.0 } else { h };
                mat[i][j] -= lambda * w * kernel(nodes[i], nodes[j]);
            }
        }

        // Gaussian elimination
        let mut aug = vec![vec![0.0; n + 1]; n];
        for i in 0..n {
            for j in 0..n { aug[i][j] = mat[i][j]; }
            aug[i][n] = rhs[i];
        }
        for col in 0..n {
            let pivot_row = (col..n).max_by(|&a, &b| aug[a][col].abs().partial_cmp(&aug[b][col].abs()).unwrap()).unwrap();
            aug.swap(col, pivot_row);
            if aug[col][col].abs() < 1e-14 { continue; }
            for row in (col + 1)..n {
                let factor = aug[row][col] / aug[col][col];
                for j in col..=n { aug[row][j] -= factor * aug[col][j]; }
            }
        }
        // Back substitution
        let mut phi = vec![0.0; n];
        for i in (0..n).rev() {
            if aug[i][i].abs() < 1e-14 { continue; }
            phi[i] = aug[i][n];
            for j in (i + 1)..n { phi[i] -= aug[i][j] * phi[j]; }
            phi[i] /= aug[i][i];
        }
        phi
    }

    pub fn nodes(&self) -> Vec<f64> {
        let h = (self.b - self.a) / (self.n_points - 1) as f64;
        (0..self.n_points).map(|i| self.a + i as f64 * h).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_identity_kernel() {
        let solver = FredholmSolver::new(21, 0.0, 1.0);
        let phi = solver.solve(0.0, |_x, _t| 0.0, |x| x * x);
        let nodes = solver.nodes();
        for (i, &n) in nodes.iter().enumerate() {
            assert!((phi[i] - n * n).abs() < 1e-6);
        }
    }
}
