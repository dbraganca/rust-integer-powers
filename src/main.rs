#[cfg(target_pointer_width = "64")]
#[allow(unused_imports)]
#[allow(unused_variables)]
use rug::{Complex, Float, Integer, float::Constant, ops::Pow};
use std::collections::HashMap;
use std::time::Instant;

// mod ltrian;
mod jfunc;


const PREC: u32 = 190;

fn main() {
    let mut cache = jfunc::ltrian::TopologyCache {
        // tadpole_cache: HashMap::new(),
        bubble_cache: HashMap::new(),
        trian_cache: HashMap::new(),
        ltrian_cache: HashMap::new(),
    };

    // let m1 = &Complex::with_val(PREC, (0.001, -0.0082));
    // let m2 = &Complex::with_val(PREC, (0.001, -0.0082));
    let m3 = &Complex::with_val(PREC, (0.000076,-0.0013));

    let kval1 = &Float::with_val(PREC, 1.5);
    let kval2 = &Float::with_val(PREC, 1.4);
    let kval3 = &Float::with_val(PREC, 1.3);
    let t1 = Instant::now();

    println!("ltrian: {}", jfunc::ltrian::ltrian(-1,0,-1,0,0,1, 
                            kval1, 
                            kval2, 
                            kval3, 
                            m3, m3, m3, 1, 2, 3, &mut cache));
    println!("computeJ: {}", jfunc::compute_j(1,1,1,1,1,1,kval1,kval2,kval3, &mut cache));
    println!(
        "Time taken bubble: {:#?}",
        Instant::now().duration_since(t1)
    );

}