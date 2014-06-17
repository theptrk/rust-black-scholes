/* The Black and Scholes (1973) Stock option formula */

use std::f64;

fn black_scholes(putCallFlag:&str, s:f64, x:f64, t:f64, r:f64, v:f64) -> f64 {

  let d1 = ( (s / x).ln() + (r + v * v / 2.0) * t) / ( v * t.sqrt() );
  let d2 = d1 - v * t.sqrt();


  if putCallFlag == "c" {
    return s * cnd(d1)-x * (-r * t).exp() * cnd(d2);
  } else {
    return x * (-r * t).exp() * cnd(-d2) - s * cnd(-d1);
  }

}

/* The cummulative Normal distribution function: */

fn cnd(x:f64) -> f64{

  let a1 = 0.31938153; 

  if x < 0.0 {
    return 1.0 - cnd(-x);
  } else {
    let k = 1.0 / (1.0 + 0.2316419 * x);
    return 1.0 - (-x * x / 2.0).exp()/ (2.0*pi()).sqrt() * k * (a1 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + k * 1.330274429)))) ;

  }

}

fn pi() -> f64{
    f64::consts::PI
}


fn main() {
    println!("{}",black_scholes("c", 31.55, 22.75, 3.5, 0.05, 0.50))
}
