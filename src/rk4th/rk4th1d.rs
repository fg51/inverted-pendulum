pub trait RungeKutta4thOrder2nd {
    fn cal_ddy(t: f64, v: f64, x: f64) -> f64;
    fn cal_dy(t: f64, v: f64, x: f64) -> f64;
    fn delta_time(&self) -> f64;
    fn current_value(&self) -> (f64, f64);
    fn update(&mut self, vs: &[f64]);

    fn solve(&mut self, t: f64) {
        let dt = self.delta_time();
        let (v, x) = self.current_value();

        let (k1v, k1x) = self.cal_slope(dt, (t, &[v, x]));
        let (k2v, k2x) = self.cal_slope(dt, (t + dt / 2., &[v + k1v / 2., x + k1x / 2.]));
        let (k3v, k3x) = self.cal_slope(dt, (t + dt / 2., &[v + k2v / 2., x + k2x / 2.]));
        let (k4v, k4x) = self.cal_slope(dt, (t + dt, &[v + k3v, x + k3x]));

        let dv = (k1v + 2. * k2v + 2. * k3v + k4v) / 6.;
        let dx = (k1x + 2. * k2x + 2. * k3x + k4x) / 6.;
        self.update(&[dv, dx]);
    }

    fn cal_slope(&self, dt: f64, (t, xs): (f64, &[f64])) -> (f64, f64) {
        let (v, x) = (xs[0], xs[1]);
        let kv1 = Self::cal_ddy(t, v, x) * dt;
        let kx1 = Self::cal_dy(t, v, x) * dt;
        (kv1, kx1)
    }
}

pub trait RungeKutta4th {
    fn cal_slope(&self, dt: f64, t: f64, vs: &[f64]) -> Vec<f64>;
    fn delta_time(&self) -> f64;
    fn current_value(&self) -> Vec<f64>;
    fn update(&mut self, vs: &[f64]);

    fn solve(&mut self, t: f64) {
        let dt = self.delta_time();
        let vs = self.current_value();

        let k1s = self.cal_slope(dt, t, &vs);
        let k2s = self.cal_slope(dt, t + dt / 2., &Self::next_half_point(&vs, &k1s));
        let k3s = self.cal_slope(dt, t + dt / 2., &Self::next_half_point(&vs, &k2s));
        let k4s = self.cal_slope(dt, t + dt, &Self::next_point(&vs, &k3s));

        let mut dvs = vec![0.; vs.len()];
        for i in 0..vs.len() {
            dvs[i] = (k1s[i] + 2. * k2s[i] + 2. * k3s[i] + k4s[i]) / 6.;
        }

        self.update(&dvs);
    }

    fn next_half_point(vs: &[f64], ks: &[f64]) -> Vec<f64> {
        vs.iter().zip(ks.iter()).map(|(v, k)| v + k / 2.).collect()
    }

    fn next_point(vs: &[f64], ks: &[f64]) -> Vec<f64> {
        vs.iter().zip(ks.iter()).map(|(v, k)| v + k).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::RungeKutta4th;
    use super::RungeKutta4thOrder2nd;

    const GRAVITY: f64 = 9.8;

    struct FreeFall1 {
        delta_time: f64,
        velocity: f64,
        position: f64,
    }

    impl FreeFall1 {
        fn new(delta_time: f64) -> Self {
            Self {
                delta_time,
                velocity: 0.,
                position: 0.,
            }
        }
    }

    impl RungeKutta4thOrder2nd for FreeFall1 {
        fn cal_ddy(_t: f64, _v: f64, _x: f64) -> f64 {
            -GRAVITY
        }

        fn cal_dy(_t: f64, v: f64, _x: f64) -> f64 {
            v
        }

        fn delta_time(&self) -> f64 {
            self.delta_time
        }

        fn current_value(&self) -> (f64, f64) {
            (self.velocity, self.position)
        }

        fn update(&mut self, vs: &[f64]) {
            self.velocity += vs[0];
            self.position += vs[1];
        }
    }

    struct FreeFall2 {
        delta_time: f64,
        velocity: f64,
        position: f64,
    }

    impl FreeFall2 {
        fn new(delta_time: f64) -> Self {
            Self {
                delta_time,
                velocity: 0.,
                position: 0.,
            }
        }

        fn ddy(_t: f64, _v: f64, _x: f64) -> f64 {
            -GRAVITY
        }

        fn dy(_t: f64, v: f64, _x: f64) -> f64 {
            v
        }
    }

    impl RungeKutta4th for FreeFall2 {
        fn delta_time(&self) -> f64 {
            self.delta_time
        }

        fn current_value(&self) -> Vec<f64> {
            vec![self.velocity, self.position]
        }

        fn cal_slope(&self, dt: f64, t: f64, vs: &[f64]) -> Vec<f64> {
            let dy = Self::ddy(t, vs[0], vs[1]) * dt;
            let y = Self::dy(t, vs[0], vs[1]) * dt;
            vec![dy, y]
        }

        fn update(&mut self, vs: &[f64]) {
            self.velocity += vs[0];
            self.position += vs[1];
        }
    }

    #[test]
    fn test_rk4th2order() {
        let error = 1E-6;
        let dt = 0.01;

        let mut t = 0.;
        let mut x = FreeFall1::new(dt);

        for _ in 0..4 {
            t += dt;
            x.solve(t);
            assert!((x.velocity - (-GRAVITY * t)).abs() < error);
            assert!((x.position - (1. / 2. * (-GRAVITY) * t * t)).abs() < error);
        }
    }

    #[test]
    fn test_rk4th() {
        let error = 1E-6;
        let dt = 0.01;

        let mut t = 0.;
        let mut x = FreeFall2::new(dt);

        for _ in 0..4 {
            t += dt;
            x.solve(t);
            assert!((x.velocity - (-GRAVITY * t)).abs() < error);
            assert!((x.position - (1. / 2. * (-GRAVITY) * t * t)).abs() < error);
        }
    }
}
