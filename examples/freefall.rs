use inverted_pendulum as lib;

use lib::vector::{Point3, Vector3};

use lib::rk4th::rk4th3d::RungeKutta4thOrder2nd;

const GRAVITY: f64 = 9.8;
const STEP: f64 = 0.01;

pub fn main() {
    let mut p = Particle::new(initial_position());

    for i in 0..100 {
        let t = i as f64 * STEP;
        p.solve(t);
        // println!("{},{},{}", t, p.velocity[2], p.position[2]);
        let actual = p.position[2];
        let theory = -GRAVITY * (t + STEP) * (t + STEP) / 2.;
        assert!((actual - theory) < 1E-4);
    }
}

pub fn initial_position() -> Point3 {
    [0., 0., 0.].into()
}

#[derive(Debug)]
pub struct Particle {
    pub position: Point3,
    pub velocity: Vector3,
    delta_time: f64,
}

impl Particle {
    pub fn new(position: Point3) -> Self {
        Self {
            position,
            velocity: Vector3::zeros(),
            delta_time: STEP,
        }
    }
}

impl RungeKutta4thOrder2nd for Particle {
    fn delta_time(&self) -> f64 {
        self.delta_time
    }

    fn cal_ddy(_t: f64, _v: &Vector3, _x: &Vector3) -> Vector3 {
        [0., 0., -GRAVITY].into()
    }

    fn cal_dy(_t: f64, v: &Vector3, _x: &Vector3) -> Vector3 {
        v.clone()
    }

    fn current_value(&self) -> (&Vector3, &Vector3) {
        (&self.velocity, &self.position.coords)
    }

    fn update(&mut self, vs: &[Vector3]) {
        self.velocity += vs[0];
        self.position += vs[1];
    }
}
