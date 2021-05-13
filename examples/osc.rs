// \frac{d^2\theta}{dt^2} = - \frac{g}{l}sin\theta

use core::f64::consts::PI;

const LENGTH: f64 = 1.; // [m]
const GRAVITY: f64 = 9.8; // [m/s^2]

fn f(radian: f64, _v: f64, _t: f64) -> f64 {
    -GRAVITY / LENGTH * radian.sin()
}

pub fn main() {
    // let m = 1.0; // mass: 1 [Kg]
    let t0 = 0.0;
    let t1 = 10.0;

    let num = 50;

    let dt = (t1 - t0) / num as f64; // time grid

    let ts = ts(dt, num);
    let mut thetas = vec![];
    let mut vs = vec![];

    // initial condition
    //let theta0 = 90. - 5.; // [degree]
    let theta0 = 180. - 5.; // [degree]
    let theta0 = PI * theta0 / 180.;
    let v0 = 0.0;

    let (mut theta, mut v) = (theta0, v0);

    for t in ts {
        thetas.push(theta);
        vs.push(v);
        let k1v = f(theta, v, t) * dt;
        let k1x = v * dt;

        let k2v = f(theta + k1x / 2., v + k1v / 2., t + dt / 2.) * dt;
        let k2x = (v + k1v / 2.) * dt;

        let k3v = f(theta + k2x / 2., v + k2v / 2., t + dt / 2.) * dt;
        let k3x = (v + k2v / 2.) * dt;

        let k4v = f(theta + k3x, v + k3v, t + dt) * dt;
        let k4x = (v + k3v) * dt;

        v += (k1v + 2. * k2v + 2. * k3v + k4v) / 6.;
        theta += (k1x + 2. * k2x + 2. * k3x + k4x) / 6.;
    }

    //let mut fout = BufWriter::new(File::create("out.csv").unwrap());

    println!("t[sec],dtheta[rad],theta[rad]");
    for i in 0..num {
        println!("{},{},{}", i as f64 * dt, vs[i], thetas[i]);
    }
}

fn ts(dt: f64, num: usize) -> Vec<f64> {
    let mut ts = vec![];
    for i in 0..num {
        ts.push(dt * i as f64);
    }
    ts
}
