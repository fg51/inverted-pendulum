use inverted_pendulum as lib;

use lib::Particle;

use lib::rk4th::rk4th3d::RungeKutta4thOrder2nd;

const BALL_MASS: f64 = 1.;
//const BAR_LENGTH: f64 = 1.;
//const BALL_POSITION: (f64, f64, f64) = (0., 0., 1.);

pub fn main() {
    //let p = Point3::new(0., 0., 1.0);
    //let mut p = Particle::new(BALL_MASS, point3(0., 0., 0.));
    let mut p = Particle::new(BALL_MASS, [0., 0., 0.].into());

    for i in 0..100 {
        let t = i as f64 * 0.01;
        p.solve(t);
        println!("{},{},{}", t, p.velocity[2], p.position[2]);
    }

    // rk
}
