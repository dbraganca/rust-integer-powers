#[cfg(target_pointer_width = "64")]
#[allow(unused_imports)]
#[allow(unused_variables)]
use rug::{Complex, Float, Integer, float::Constant, ops::Pow};
use num_integer::binomial;
use std::collections::HashMap;
use std::time::Instant;

const PREC: u32 = 190;
const PI: f64 = 3.14159265358979323846264338327950288f64;
const SQRT_PI: f64 = 1.77245385090551602729816748334114518279754945612238f64;
const CHOP_TOL: f64 = 1e-30;

struct TopologyCache {
    // tadpole_cache: HashMap<(isize, isize), Complex>,
    bubble_cache: HashMap<(isize, isize, 
                            usize, usize), Complex>,
    trian_cache: HashMap<(isize, isize, isize, 
                            usize, usize, usize), Complex>,
    // Possibly add here 
    // ltrian_cache: HashMap<(isize, isize, isize, isize, isize, isize, usize, usize, usize), Complex>,
    // Also possibly add in bubble_cache and trian_cache a index corresponding to the k - 1, 2 or 3
}

fn main() {
    let mut cache = TopologyCache {
        // tadpole_cache: HashMap::new(),
        bubble_cache: HashMap::new(),
        trian_cache: HashMap::new(),
    };

    let m1 = &Complex::with_val(PREC, (0.034, 0.069));
    let m2 = &Complex::with_val(PREC, (0.001, -0.0082));
    // let m3 = &Complex::with_val(PREC, (0.001, 0.0082));
    let m3 = &Complex::with_val(PREC, (0.000076, -0.0013));

    let kval1 = &Float::with_val(PREC, 1.5);
    let kval2 = &Float::with_val(PREC, 1.4);
    let kval3 = &Float::with_val(PREC, 1.3);
    
    let t1 = Instant::now();

    println!("ltrian: {}", ltrian(0, 1, 0, 1, 0, 2, 
                            kval1, 
                            kval2, 
                            kval3, 
                            m3, m3, m3, 1, 2, 3, &mut cache));
    println!(
        "Time taken trian 1: {:#?}",
        Instant::now().duration_since(t1)
    );

    // let t1 = Instant::now();

    // println!("ltrian: {}", ltrian(0, 10, 0, 10, 0, 10, 
    //                         kval, 
    //                         kval, 
    //                         kval, 
    //                         m1, m2, m3, 1, 2, 3, &mut cache));
    // println!(
    //     "Time taken trian 2: {:#?}",
    //     Instant::now().duration_since(t1)
    // );
}

fn k1dotk2(k12: &Float, k22: &Float, ksum2: &Float) -> Float {
    return (ksum2.clone() - k12 - k22)/2
}

fn trinomial(n: usize, k: usize, i: usize) -> usize {
    return binomial(n, k + i) * binomial(k + i, i)
}

fn get_coef_simple(n1: usize, kmq2exp: usize) -> usize {
    return binomial(n1, kmq2exp)
}

fn get_coef(n1: usize, k2exp: usize, q2exp: usize, kmq: bool) -> isize {
    let kqexp: u32 = (n1 - k2exp - q2exp) as u32;
    let mut sign: isize = 1;

    if kmq & (kqexp%2 != 0) {
        sign = -1;
    }
    return sign * (trinomial(n1, k2exp, q2exp) as isize) * 2isize.pow(kqexp)
}

fn num_terms(n1: usize, kmq: bool) -> (Vec<isize>, Vec<[usize;3]>) {
    // expands terms of the type (k-q)^(2n1) if kmq = True
	// expands terms of the type (k+q)^(2n1) if kmq = False
    let list_length: usize = ((1 + n1)*(2 + n1)/2) as usize;
    let mut term_list: Vec<isize> = vec![0; list_length];
    let mut exp_list = vec![[0,0,0]; list_length];
    let mut i: usize = 0;

    for k2exp in 0..n1+1 {
        for q2exp in 0..n1-k2exp+1{
            term_list[i] = get_coef(n1,k2exp,q2exp,kmq);
			exp_list[i][0] = k2exp;
			exp_list[i][1] = q2exp;
			exp_list[i][2] = n1-k2exp-q2exp;
            i += 1;
        }
    }
    return (term_list, exp_list)
} 

fn dim_gen(expnum: isize, expden: isize, m: &Complex) -> Complex {
    let complex_zero = Complex::new(PREC);

    if m == &0. {
        return complex_zero
    }

    let expnum_float: Float = Float::with_val(PREC,expnum);
    let expden_float: Float = Float::with_val(PREC,expden);

    return (2./SQRT_PI) * m.clone().pow(expnum as i64 - expden as i64 + 1) * m.clone().sqrt()
        * (expnum_float.clone() + 1.5f64).gamma() 
        * (expden_float.clone() - expnum_float - 1.5f64).gamma() / expden_float.gamma()        
}

fn expand_massive_num(n1: usize) -> (Vec<usize>, Vec<[usize;2]>) {
	//expand term of the form ((kmq)^2+mNum)^n1
	//returns:
	// - term_list: list of coefficients
	// - exp_list: list of exponents [expkmq2,mNumexp]

    let list_length: usize = n1+1;
    let mut term_list: Vec<usize> = vec![0; list_length];
    let mut exp_list = vec![[0,0]; list_length];

	for kmq2exp in 0..n1+1 {
		term_list[kmq2exp] = get_coef_simple(n1, kmq2exp);
		exp_list[kmq2exp][0] = kmq2exp;
		exp_list[kmq2exp][1] = n1-kmq2exp;
    }

    return (term_list, exp_list)
}

fn dim_result(expkmq2: usize, expden: usize, k2: &Float, m_den: &Complex, kmq: bool) -> Complex {
    // computes integrals of the type ((k-q)^2)^expnum / (q^2 + m_den)^expden
    // if kmq is TRUE, consider (k-q)
    // if kmq is FALSE, consider (k+q)

    let list_length: usize = (1 + expkmq2)*(2 + expkmq2)/2;
    let (term_list, exp_list) = num_terms(expkmq2,kmq);

    let mut res: Complex = Complex::new(PREC);

    for i in 0..list_length {
        let k2exp = exp_list[i][0];
        let q2exp = exp_list[i][1];
        let kqexp = exp_list[i][2];

        if kqexp%2 == 0 {
            if kqexp != 0 {
                res += term_list[i] as i64
                    * dim_gen((q2exp+kqexp/2) as isize, expden as isize, m_den) 
                    * k2.clone().pow((k2exp + kqexp/2) as u64)/((1+kqexp) as u64) ;
            } 
            else {
                res += term_list[i] as i64 
                    * dim_gen(q2exp as isize, expden as isize, m_den) 
                    * k2.clone().pow(k2exp as u64);
            }
        }
    }    
    return res
}

fn compute_massive_num(exp_num: usize, exp_den: usize, k2: &Float, 
                    m_num: &Complex, m_den: &Complex, kmq: bool) -> Complex {
    // computes integrals of the type ((k-q)^2+mNum)^expnum / (q^2 + mDen)^expden if kmq = true 
    // or ((k+q)^2+mNum)^expnum / (q^2 + mDen)^expden if kmq = false
    // by summing over many dim_result

    let complex_zero = Complex::new(PREC);
    
    if m_num.clone() == complex_zero {
        return dim_result(exp_num, exp_den, k2, m_den, kmq);
    }

    let list_length: usize = exp_num + 1;
    let mut kmq2_exp: usize;
    let mut m_num_exp: usize;
    let (term_list, exp_list) = expand_massive_num(exp_num);
        
    let mut res= Complex::new(PREC);

    for i in 0..list_length {
        kmq2_exp = exp_list[i][0];
		m_num_exp = exp_list[i][1];

		res += term_list[i] as u64
                *(m_num.clone().pow(m_num_exp as u64))*dim_result(kmq2_exp, exp_den, k2, m_den, kmq);
    }

    return res
}


fn tadn(n: isize, m: &Complex) -> Complex {
    let complex_zero = Complex::new(PREC);
    
    if m == &0. {
        return complex_zero
    }
    
    if n == 0 {
        return complex_zero;
    } 
    if n == 1 {
        return (-2) * m.clone().sqrt() * SQRT_PI;
    }
    if n < 0 {
        return complex_zero
    }

    // Code to use cache for TadN:
    // if let Some(c) = cache.tadpole_cache.get(&(n, pdg)) {
    //     return c.clone();
    // }
    // let dim = Float::with_val(PREC, 3);
    // let nu = n - 1;
    // let c0 = 1 / m.clone() * (1 - dim / Float::with_val(PREC, 2 * nu));

    // let r: Complex = c0 * tadn(nu, m, pdg, cache);
    // cache.tadpole_cache.insert((n, pdg), r.clone());
    return dim_gen(0, n, m)
}

fn bubmaster_zero_masses(k2: &Float) -> Complex {
    return Complex::with_val(PREC, (PI * SQRT_PI / k2.clone().sqrt(),0))
}

fn bubmaster(k2: &Float, m1in: &Complex, m2in: &Complex) -> Complex {

    if m1in == &0. && m2in == &0. {
        return bubmaster_zero_masses(k2)
    }

    let m1: Complex = m1in.clone()/k2;
    let m2: Complex = m2in.clone()/k2;
    let mut sign = 0.;

    let arg_log_1: Complex = (m1.clone() - &m2 - 1i8).mul_i(false) + 2 * m1.clone().sqrt();
    let arg_log_2: Complex = (m1.clone() - &m2 + 1i8).mul_i(false) + 2 * m2.clone().sqrt();

    if arg_log_1.imag() > &0. && arg_log_2.imag() < &0. 
    {
        sign = 1.;
    }
    
    let bubmaster = (SQRT_PI/k2.clone().sqrt()) * ((arg_log_1.ln() - arg_log_2.ln()).mul_i(false) + 2.* PI * sign);
    
    return bubmaster
}

fn bubn(
    n1: isize,
    n2: isize,
    k2: &Float,
    m1: &Complex,
    m2: &Complex,
    pdg1: usize,
    pdg2: usize,
    cache: &mut TopologyCache,
) -> Complex {

    if let Some(r) = cache.bubble_cache.get(&(n1, n2, pdg1, pdg2)) {
        return r.clone();
    }

    if n1 == 0 {
        return tadn(n2, m2);
    } 
    if n2 == 0 {
        return tadn(n1, m1);
    } 
    if n1 == 1 && n2 == 1 {
        return bubmaster(k2, m1, m2);
    }

    // Code to deal with numerators
    if n1 < 0 || n2 < 0 {
        if m1 == &0. && m2 == &0. {
            return Complex::new(PREC)
        } 
        if n1 < 0 && n2 > 0 {
        // m1 is the mass in the numerator
        // m2 is the mass in the denominator 
            return compute_massive_num((-n1) as usize, n2 as usize, k2, m1, m2, true)
        }
        if n2 < 0 && n1 > 0 {
        // m2 is the mass in the numerator
        // m1 is the mass in the denominator
            return compute_massive_num((-n2) as usize, n1 as usize, k2, m2, m1, true)
        }
    }

    let k1s = k2.clone() + m1.clone() + m2;
    let jac = &k1s * &k1s - 4 * m1.clone() * m2;
    let dim: isize = 3;
    let (nu1, nu2, cpm0, cmp0, c000) = if n1 > 1 {
        let nu1 = n1 - 1 ;
        let nu2 = n2;
        let ndim = (dim - nu1 - nu2) as i64;

        (
            nu1,
            nu2,
            k1s.clone(),
            -2 * m2.clone() / nu1 as i64 * nu2 as i64,
            (2 * m2.clone() - &k1s) / nu1 as i64 * ndim - 2 * m2.clone() + k1s * nu2 as i64 / nu1 as i64,
        )
    } else {
        let nu1 = n1;
        let nu2 = n2 - 1;
        let ndim = (dim - nu1 - nu2) as i64;

        (
            nu1,
            nu2,
            -2 * m1.clone() / nu2 as i64 * nu1 as i64,
            k1s.clone(),
            (2 * m1.clone() - &k1s) / nu2 as i64 * ndim + k1s / nu2 as i64 * nu1 as i64 - 2 * m1.clone(),
        )
    };

    let c000 = c000 / &jac;
    let cmp0 = cmp0 / &jac;
    let cpm0 = cpm0 / &jac;

    let r: Complex = c000 * bubn(nu1, nu2, k2, m1, m2, pdg1, pdg2, cache)
        + cpm0 * bubn(nu1 + 1, nu2 - 1, k2, m1, m2, pdg1, pdg2, cache)
        + cmp0 * bubn(nu1 - 1, nu2 + 1, k2, m1, m2, pdg1, pdg2, cache);

    cache.bubble_cache.insert((n1, n2, pdg1, pdg2), r.clone());

    return r
}

fn trian_kinem(k21: &Float , k22: &Float, k23: &Float , 
                m1: &Complex, m2: &Complex, m3: &Complex) -> [Complex; 7] {

    let k1s =  k21.clone() + m1.clone() + m2;
    let k2s =  k22.clone() + m2.clone() + m3;
    let k3s =  k23.clone() + m3.clone() + m1;

    let jac = -4*m1.clone() * m2 * m3 + k1s.clone() * &k1s * m3 
                + k2s.clone() * &k2s * m1 + k3s.clone() * &k3s * m2 
                - k1s.clone() * &k2s * &k3s;
    let jac = 2*jac;

    let ks11 = (- 4 * m1.clone() * m2 + &k1s * &k1s) / &jac;
	let ks12 = (-2 * k3s.clone() * m2 + &k1s * &k2s) / &jac;
	let ks22 = (-4* m2.clone() * m3 + &k2s * &k2s) / &jac;
    let ks23 = (-2* k1s.clone() * m3 + &k2s * &k3s) / &jac;
    let ks31 = (-2* k2s.clone() * m1 + &k1s * &k3s)/&jac;
	let ks33 = (-4 * m1.clone() * m3 + &k3s * &k3s)/&jac;

    return [jac, ks11, ks22, ks33, ks12, ks23, ks31]
}

fn tria_master_zero_masses(k21: &Float, k22: &Float, k23: &Float) -> Complex {
    return Complex::with_val(PREC,(PI * SQRT_PI/k21.clone().sqrt()/k22.clone().sqrt()/k23.clone().sqrt(),0))
}

fn tria_master(k21: &Float, k22: &Float, k23: &Float,
    m1: &Complex,
    m2: &Complex,
    m3: &Complex) -> Complex {
    // PLACEHOLDER
    if m1 == &0. && m2 == &0. && m3 == &0. {
        return tria_master_zero_masses(k21, k22, k23)
    }
    return trmxy(1, k21, k22, k23, m1, m2, m3) - trmxy(0, k21, k22, k23, m1, m2, m3)
}

fn trmxy(y: i64, k21: &Float, k22: &Float, k23: &Float,
         m1: &Complex, m2: &Complex, m3: &Complex) -> Complex {
    let num1: Float = 4*k22.clone()*y + 2*k21.clone() - 2*k22.clone() - 2*k23.clone();
    let num0: Complex = -4*k22.clone()*y+2*m2.clone()-2*m3.clone()+2*k22.clone();
	let deltar2: Float = -k21.clone()*y + k23.clone()*y - k23;
	let deltar1: Complex = -m2.clone()*y+m3.clone()*y+k21.clone()*y-k23.clone()*y+ m1 - m3 + k23;
	let deltar0: Complex = m2.clone()*y-m3.clone()*y+m3;
	let deltas2: Float = (-1)*k21.clone().pow(2) + 2*k21.clone()*k22 + 2*k21.clone()*k23 - k22.clone().pow(2) +2*k22.clone()*k23-k23.clone().pow(2);
	let deltas1: Complex = -4*m1.clone()*k22 - 2*m2.clone()*k21 + 2*m2.clone()*k22
                           +2*m2.clone()*k23 + 2*m3.clone()*k21 + 2*m3.clone()*k22
                           -2*m3.clone()*k23 - 2*k21.clone()*k22 + 2*k22.clone().pow(2)-2*k22.clone()*k23;
	let deltas0: Complex = (-1)*(m2.clone().pow(2)) + 2*m2.clone()*m3 - 2*m2.clone()*k22
                           -m3.clone().pow(2)-2*m3.clone()*k22-k22.clone().pow(2);

    let diakrs: Complex = Complex::from(deltas1.clone().pow(2) - 4 * deltas2.clone() * deltas0).sqrt();
    let sols1: Complex = (-deltas1.clone() + &diakrs)/2/&deltas2;
    let sols2: Complex = (-deltas1 - &diakrs)/2/deltas2;

    let cf2: Complex = -(sols2.clone() * &num1  + &num0)/&diakrs;
    let cf1: Complex = (sols1.clone() * num1 + num0)/&diakrs;

    let diakrr: Complex = Complex::from(deltar1.clone().pow(2) - 4 * deltar2.clone() * deltar0).sqrt();
    let solr1: Complex = (-deltar1.clone() + &diakrr)/2/&deltar2;
    let solr2: Complex = (-deltar1 - diakrr)/2/&deltar2;
    // println!("cf1: {}, cf2: {}", cf1, cf2);
	// println!("Fint cf2 = {}",f_int(&deltar2, &solr1, &solr2, &sols2));
	// println!("Fint cf1 = {}",f_int(&deltar2, &solr1, &solr2, &sols1));
    if (cf1.clone()).abs().real() < &CHOP_TOL {
        return cf2 * f_int(&deltar2, &solr1, &solr2, &sols2)
    }
    if (cf2.clone()).abs().real() < &CHOP_TOL {
        return cf1 * f_int(&deltar2, &solr1, &solr2, &sols1)
    }
    return cf2 * f_int(&deltar2, &solr1, &solr2, &sols2) + cf1 * f_int(&deltar2, &solr1, &solr2, &sols1)
}

fn f_int(aa: &Float, y1in: &Complex, y2in: &Complex, x0in: &Complex) -> Complex {
    
    let y1: Complex = if (y1in.imag().clone()).abs() < CHOP_TOL {
        Complex::from(y1in.real().clone())
    } else {y1in.clone()};

    let y2: Complex = if (y2in.imag().clone()).abs() < CHOP_TOL {
        Complex::from(y2in.real().clone())
    } else {y2in.clone()};

    let x0: Complex = if (x0in.imag().clone()).abs() < CHOP_TOL {
        Complex::from(x0in.real().clone())
    } else {x0in.clone()};

    let rey1 = y1.real();
    let imy1 = y1.imag();
    let rey2 = y2.real();
    let imy2 = y2.imag();
    let rex0 = x0.real();
    let imx0 = x0.imag();

    let c: Float = imy1.clone().pow(2) * imy2 * rex0 - imy1 * imy2.clone().pow(2)*rex0
          - imx0.clone().pow(2) * imy2 * rey1 + imx0 * imy2.clone().pow(2) * rey1
          - imy2 * rex0.clone().pow(2) * rey1 + imy2 * (rex0 * rey1.clone().pow(2)) 
          + imx0.clone().pow(2) * imy1 * rey2 - imx0 * imy1.clone().pow(2) * rey2
          + imy1 * rex0.clone().pow(2) * rey2 - imx0 * rey1.clone().pow(2) * rey2
          - imy1 * (rex0 * rey2.clone().pow(2)) + imx0 * (rey1 * rey2.clone().pow(2));
    
    let a: Float = imy1.clone() * rex0 - imy2 * rex0 
                 - imx0 * rey1 + imy2 * rey1
                 + imx0 * rey2 - imy1 * rey2;
    
    let b: Float = (-1)*imx0.clone().pow(2) * imy1 + imx0 * imy1.clone().pow(2)
          + imx0.clone().pow(2) * imy2 - imy1.clone().pow(2) * imy2 
          - imx0 * imy2.clone().pow(2) + imy1 * imy2.clone().pow(2)
          - imy1 * rex0.clone().pow(2) + imy2 * rex0.clone().pow(2)
          + imx0 * rey1.clone().pow(2) - imy2 * rey1.clone().pow(2)
          - imx0 * rey2.clone().pow(2) + imy1 * rey2.clone().pow(2);
    
    // println!("a: {}, b: {}, c: {}", a, b, c);
    // if x0 is real there will always be a crossing through i or -i, which gives a cut of pi/2 instead of pi
    let mut cutx0 = Complex::new(PREC);

    if &0 < rex0 && rex0 < &1 && imx0.clone().abs() < CHOP_TOL {
		let derivcritx0: Complex = (y1.clone() - &y2)/2/(rex0 - y2.clone())/(-(-y1.clone() + rex0).pow(2u8)).sqrt();
		let signx0: f64;
        if derivcritx0.real() < &0 {
			signx0 = 1.
        } else {
			signx0 = -1.	
        }
		cutx0 = signx0 * PI/((y1.clone() - rex0).sqrt()*(rex0 - y2.clone()).sqrt())
    } 
    
    // find remaining crossings of the imaginary axis
    // let xsol: Vec<Float> = Vec::with_capacity(2);
    let xsol: Vec<Float> =  
        if a.clone().abs() < CHOP_TOL {
            if b != 0 {
                vec![-c/b]
            } else {
                vec![]
            }
        } else {
            let delta: Float = b.clone().pow(2)  - 4 * c * &a;
            if delta > 0 {
                vec![(delta.clone().sqrt() - &b)/2/&a, (-b - delta.sqrt())/2/a] 
            } else {
                // case where there is no intersection of the real axis (includes double zero)
                vec![]
            }
        };

    let mut xsolnew: Vec<Float> = Vec::with_capacity(2);
    for i in 0..xsol.len() {
        if xsol[i] > CHOP_TOL && xsol[i] < Float::with_val(PREC,1) - CHOP_TOL 
           && (&xsol[i] - x0.clone()).abs().real() > &CHOP_TOL {
            xsolnew.push(xsol[i].clone());
        }
    }

    let numxsol = xsolnew.len();
    let mut num_branch_points = 0;
    let mut xbranch = Vec::with_capacity(2);
    let mut cut = Complex::new(PREC);
    // WE CAN SIMPLIFY THIS ALGORITHM A LOT - NO NEED FOR 2 FOR LOOPS
    if  numxsol > 0 {
        let mut atanarglist: Vec<Complex> = Vec::with_capacity(2);
        let mut imcrit: Vec<Float> = Vec::with_capacity(2);
        let mut recrit: Vec<Float> = Vec::with_capacity(2);
        for i in 0..numxsol {
            let atanarg = (&xsolnew[i] - y1.clone()).sqrt() * (x0.clone() - &y2).sqrt()
                     / (-x0.clone() + &y1).sqrt() / (&xsolnew[i] - y2.clone()).sqrt();
            atanarglist.push(atanarg.clone());
            imcrit.push(atanarg.imag().clone());
            recrit.push(atanarg.real().clone());
        }

        for i in 0..numxsol {
            if imcrit[i] > 1 && recrit[i].clone().abs() < CHOP_TOL {
                num_branch_points += 1;
                xbranch.push(xsolnew[i].clone());
            }
        }

        if num_branch_points == 1 {
            let sign: f64;
            let xcrit = &xbranch[0];
            let derivcrit: Complex = (x0.clone() - &y2).sqrt() * (&y1 - y2.clone())/2/(xcrit - y1.clone()).sqrt()/(-x0.clone() + &y1).sqrt()/ (xcrit - y2.clone())/ (xcrit - y2.clone()).sqrt();
            if derivcrit.real() < &0 {
                sign = 1.
            } else {
                sign = -1.
            }
            cut = sign * PI * 2./(-x0.clone() + &y1).sqrt()
                  /(x0.clone() - &y2).sqrt();        
        }
    }
    
    let prefac0 = prefactor(aa, &y1, &y2);
    // println!("aa: {}, y1: {}, y2: {}, x0: {}", aa, y1, y2, x0);
    // println!("prefac0: {}, cut: {}, cutx0: {}", prefac0, cut, cutx0);
    // println!("anti_deriv1: {}, anti_deriv0: {}", anti_deriv(1, &y1, &y2, &x0), anti_deriv(0, &y1, &y2, &x0));
    let result = prefac0 * (SQRT_PI/2.) * (cut + cutx0 + anti_deriv(1, &y1, &y2, &x0) - anti_deriv(0, &y1, &y2, &x0));

    return result
}

fn prefactor(a: &Float, y1: &Complex, y2: &Complex) -> Complex {
    // We can simplify this function a lot more
    let y2re = y2.real();
    let y1re = y1.real();
    let y2im = y2.imag();
    let y1im = y1.imag();

    if y2im.clone().abs() < CHOP_TOL && y1im.clone().abs() < CHOP_TOL {
		if y1re.clone().abs() >= CHOP_TOL && y2re.clone().abs() >= CHOP_TOL { 
			return Complex::from((-y1re.clone()).sqrt()*(-y2re.clone()).sqrt()/(a.clone() * y1re * y2re).sqrt())
        }
		else if y1re.clone().abs() < CHOP_TOL && y2re.clone().abs() >= CHOP_TOL {
			return Complex::from((-y2re.clone()).sqrt()/(-a.clone()*y2re).sqrt())
	    }	
        else if y1re.clone().abs() >= CHOP_TOL && y2re.clone().abs() < CHOP_TOL {
			return Complex::from((-y1re.clone()).sqrt()/(-a.clone()*y1re).sqrt())
        }
		else {
        // if y1re.clone().abs() < CHOP_TOL && y2re.clone().abs() < CHOP_TOL
			return Complex::from(1/a.clone().sqrt())
        }

	} else if y2im.clone().abs() >= CHOP_TOL && y1im.clone().abs() < CHOP_TOL {
		if y1re.clone().abs() >= CHOP_TOL { 
			return (-y1re.clone()).sqrt()*(-y2.clone()).sqrt()/(y2.clone() * a * y1re).sqrt()
        }
		else {
			return (-y2.clone()).sqrt()/(-y2.clone()*a).sqrt()
        }

	} else if y2im.clone().abs() < CHOP_TOL && y1im.clone().abs() >= CHOP_TOL {
		if y2re.clone().abs() >= CHOP_TOL { 
			return (-y1.clone()).sqrt()*(-y2re.clone()).sqrt()/(a*y1.clone()*y2re).sqrt()
        }
		else {
			return (-y1.clone()).sqrt()/(-y1.clone() * a).sqrt()
        }
	} else {
		// case where abs(y2.imag) >= CHOP_TOL && abs(y1.imag) >= CHOP_TOL
		return (-y1.clone()).sqrt()*(-y2.clone()).sqrt()/(a*y1.clone()*y2).sqrt()
    }
}

fn anti_deriv(x: i64, y1: &Complex, y2: &Complex, x0: &Complex) -> Complex {
    if (x0.clone() - y2).abs().real() < &CHOP_TOL {
        if (x - y2.clone()).abs().real() < &CHOP_TOL {
            return Complex::new(PREC)
        }
        else {
            return 2. * (x - y1.clone()).sqrt()/(-x0.clone() + y1)/(x - y2.clone()).sqrt()
        }
    }

    if (x0.clone() - y1).abs().real() < &CHOP_TOL {
        println!("WARNING: switching var in anti_deriv");
        return anti_deriv(x, y2, y1, x0)
    }

    let prefac: Complex = 2./(-x0.clone() + y1).sqrt()/(x0.clone() - y2).sqrt();
    let temp: Complex = (x - y1.clone()).sqrt() * (x0.clone() - y2).sqrt() / (-x0.clone() + y1).sqrt();

    if x == 1 && (y2.clone() - 1f64).abs().real() < &CHOP_TOL {
        let lim_arctan: Complex = (-temp.clone().pow(2u8)).sqrt().mul_i(false) * PI / 2 / &temp;
        return lim_arctan * prefac
    }
    if x == 0 && y2.clone().abs().real() < &CHOP_TOL {
        let lim_arctan: Complex = (temp.clone().pow(2u8)).sqrt() * PI / 2 / temp;
        return lim_arctan * prefac
    }
    
    return prefac * (temp/(x-y2.clone()).sqrt()).atan()
}

fn trian(
    n1: isize,
    n2: isize,
    n3: isize,
    k21: &Float,
    k22: &Float,
    k23: &Float,
    m1: &Complex,
    m2: &Complex,
    m3: &Complex,
    pdg1: usize,
    pdg2: usize,
    pdg3: usize,
    cache: &mut TopologyCache,
) -> Complex {

    if let Some(r) = cache.trian_cache.get(&(n1, n2, n3, pdg1, pdg2, pdg3)) {
        return r.clone();
    }

    if n1 == 0 {
        return bubn(n2, n3, k22, m2, m3, pdg2, pdg3, cache)
    }
    if n2 == 0 {
        return bubn(n3, n1, k23, m3, m1, pdg3, pdg1, cache)
    }
    if n3 == 0 {
        return bubn(n1, n2, k21, m1, m2, pdg1, pdg2, cache)
    }

    if n1 == 1 && n2 == 1 && n3 == 1 {
        return tria_master(k21, k22, k23, m1, m2, m3)
    }

    if n1 < 0 || n2 < 0 || n3 < 0 {
        if n1 < -4 || n2 < -4 || n3 < -4 {
            println!("ERROR: case not considered -  n1 = {}, n2 = {}, n3 = {}", n1, n2, n3);
        }
        if n1 < 0 {
            if n2 > 0 && n3 > 0 {
                // in this case m1 is necessarily 0
                return tri_dim(-n1, n2, n3, k21, k22, k23, m2, m3, pdg2, pdg3, cache)
            } else if n2 < 0 {
                // in this case m1 and m2 are necessarily 0
                return tri_dim_two(-n2, -n1, n3, k22, k23, k21, m3, pdg3, cache)
            } else if n3 < 0 {
                // in this case m1 and m3 are necessarily 0
                return tri_dim_two(-n1, -n3, n2, k21, k22, k23, m2, pdg2, cache)
            }
        } else if n2 < 0 {
            if n3 > 0 {
                return tri_dim(-n2,n1,n3,k21,k23,k22,m1,m3, pdg1, pdg3, cache)
            } else if n3 < 0 {
                return tri_dim_two(-n3,-n2,n1,k23,k21,k22,m1, pdg1, cache)
            }
        } else if n3 < 0 {
            if n1 > 0 && n2 > 0 {
                return tri_dim(-n3,n1,n2,k23,k21,k22,m1,m2, pdg2, pdg3, cache)
            } else {
                println!("ERROR: case not considered");
            }
        }
    }  

    let kinem = trian_kinem(k21, k22, k23, m1, m2, m3);
    let [_jac, ks11, ks22, ks33, ks12, ks23, ks31] = kinem;

    let dim: i64 = 3;
    let (nu1, nu2, nu3, cpm0, cmp0, cm0p, cp0m, c0pm, c0mp, c000) = if n1 > 1 {
        let nu1: i64 = n1 as i64 - 1;
        let nu2: i64 = n2 as i64;
        let nu3: i64 = n3 as i64;

        let ndim: i64 = dim - nu1 - nu2 - nu3;
        let cpm0 = -ks23.clone();
        let cmp0 = (ks22.clone() * nu2)/nu1;
        let cm0p = (ks22.clone() * nu3)/nu1;
        let cp0m = -ks12.clone();
        let c0pm = -(ks12.clone() * nu2)/nu1;
        let c0mp = -(ks23.clone() * nu3)/nu1;
        let c000 = (-nu3 + ndim)*ks12/nu1 - (-nu1 + ndim)*ks22/nu1 + (-nu2 + ndim)*ks23/nu1;

        (
            nu1,
            nu2,
            nu3, 
            cpm0,
            cmp0,
            cm0p,
            cp0m,
            c0pm,
            c0mp,
            c000
        )
    } else if n2 > 1 {
        let nu1: i64 = n1 as i64;
        let nu2: i64 = n2 as i64 - 1;
        let nu3: i64 = n3 as i64;
        let ndim: i64 = dim - nu1 - nu2 - nu3;

        let cpm0 = (ks33.clone()*nu1)/nu2;
        let cmp0 = -ks23.clone();
        let cm0p = -(ks23.clone()*nu3)/nu2;
        let cp0m = -(ks31.clone()*nu1)/nu2;
        let c0pm = -ks31.clone();
        let c0mp = (ks33.clone()*nu3)/nu2;
        let c000 = (-nu1 + ndim)*ks23/nu2 + (-nu3 + ndim)*ks31/nu2 - (-nu2 + ndim)*ks33/nu2;

        (
            nu1,
            nu2,
            nu3, 
            cpm0,
            cmp0,
            cm0p,
            cp0m,
            c0pm,
            c0mp,
            c000
        )
    } else {
        // case where n3 > 1
        let nu1: i64 = n1 as i64;
        let nu2: i64 = n2 as i64;
        let nu3: i64 = n3 as i64 - 1;
        let ndim: i64 = dim - nu1 - nu2 - nu3;

        let cpm0 =  -(ks31.clone()*nu1)/nu3;
        let cmp0 = -(ks12.clone()*nu2)/nu3;
        let cm0p = -ks12.clone();
        let cp0m =  (ks11.clone()*nu1)/nu3;
        let c0pm = (ks11.clone()*nu2)/nu3;
        let c0mp = -ks31.clone();
        let c000 = -(-nu3 + ndim)*ks11/nu3 + (-nu1 + ndim)*ks12/nu3 + (-nu2 + ndim)*ks31/nu3;

        (
            nu1,
            nu2,
            nu3, 
            cpm0,
            cmp0,
            cm0p,
            cp0m,
            c0pm,
            c0mp,
            c000
        )
    };

    let nu1 = nu1 as isize;
    let nu2 = nu2 as isize;
    let nu3 = nu3 as isize;

    let res: Complex = c000*trian(nu1, nu2, nu3, k21, k22, k23, m1, m2, m3, pdg1, pdg2, pdg3, cache) 
                    + c0mp*trian(nu1, nu2-1, nu3+1, k21, k22, k23, m1, m2, m3, pdg1, pdg2, pdg3, cache) 
                    + c0pm*trian(nu1, nu2+1, nu3-1, k21, k22, k23, m1, m2, m3, pdg1, pdg2, pdg3, cache)
                    + cm0p*trian(nu1-1, nu2, nu3+1, k21, k22, k23, m1, m2, m3, pdg1, pdg2, pdg3, cache)
                    + cp0m*trian(nu1+1, nu2, nu3-1, k21, k22, k23, m1, m2, m3, pdg1, pdg2, pdg3, cache)
                    + cmp0*trian(nu1-1, nu2+1, nu3, k21, k22, k23, m1, m2, m3, pdg1, pdg2, pdg3, cache)
                    + cpm0*trian(nu1+1, nu2-1, nu3, k21, k22, k23, m1, m2, m3, pdg1, pdg2, pdg3, cache);

    cache.trian_cache.insert((n1, n2, n3, pdg1, pdg2, pdg3), res.clone());
    return res
}


fn tri_dim(n1: isize, d1: isize, d2: isize, numk2: &Float, denk2: &Float, ksum2: &Float,
            m1: &Complex, m2: &Complex, pdg1: usize, pdg2: usize, cache: &mut TopologyCache) -> Complex {

    // integral of (numk-q)^2n1/(q^2+m1)^2d1/((denk+q)^2+m2)^2d2
    // numerator (numk-q)^2n1 is massless
    // m1 is mass of d1 propagator, which is (q^2+m1)^2d1,
    // m2 is mass of d2 propagator, which is ((denk+q)^2+m2)^2d2

    let list_length = (1 + n1 as usize)*(2 + n1 as usize)/2;
    let (term_list, exp_list) = num_terms(n1 as usize, true);
    let mut res = Complex::new(PREC);

    for i in 0..list_length {
        let mut term = Complex::new(PREC);
        let k2exp = exp_list[i][0];
		let q2exp = exp_list[i][1];
		let kqexp = exp_list[i][2];
        if kqexp == 0 {
            if q2exp == 0 {
                term = bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache);
            } else if q2exp == 1 {
                term = bubn(d1-1, d2, denk2, m1, m2, pdg1, pdg2, cache)-m1.clone()*bubn(d1,d2,denk2,m1,m2, pdg1, pdg2, cache);
            } else if q2exp == 2 {
                term = bubn(d1-2,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                    - 2*m1.clone()*bubn(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                    + m1.clone().pow(2)*bubn(d1,d2,denk2,m1,m2, pdg1, pdg2, cache);
            } else if q2exp == 3 {
                term = bubn(d1-3,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                - 3*m1.clone()*bubn(d1-2,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                + 3*m1.clone().pow(2)*bubn(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                - m1.clone().pow(3)*bubn(d1,d2,denk2,m1,m2, pdg1, pdg2, cache);
            } else if q2exp == 4 {
                term = bubn(d1-4,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                - 4*m1.clone()*bubn(d1-3,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                + 6*m1.clone().pow(2)*bubn(d1-2,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                - 4*m1.clone().pow(3)*bubn(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                + m1.clone().pow(4)*bubn(d1,d2,denk2,m1,m2, pdg1, pdg2, cache);
            } else {
                println!("ERROR: exceeded calculable power in tri_dim");
            }
        } else if kqexp == 1 {
            if q2exp == 0 {
                term = num_one_pow(d1, d2, denk2, m1, m2, pdg1, pdg2, cache)*k1dotk2(numk2,denk2,ksum2);
            } else if q2exp == 1 {
                term = (num_one_pow(d1-1, d2, denk2, m1, m2, pdg1, pdg2, cache)-m1.clone()*num_one_pow(d1,d2,denk2,m1,m2, pdg1, pdg2, cache))*k1dotk2(numk2,denk2,ksum2);
            } else if q2exp == 2 {
                term = (num_one_pow(d1-2,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                    - 2*m1.clone()*num_one_pow(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                    + m1.clone().pow(2)*num_one_pow(d1,d2,denk2,m1,m2, pdg1, pdg2, cache))*k1dotk2(numk2,denk2,ksum2);
            } else if q2exp == 3 {
                term = (num_one_pow(d1-3,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                - 3*m1.clone()*num_one_pow(d1-2,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                + 3*m1.clone().pow(2)*num_one_pow(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                - m1.clone().pow(3)*num_one_pow(d1,d2,denk2,m1,m2, pdg1, pdg2, cache))*k1dotk2(numk2,denk2,ksum2);
            }
        } else if kqexp == 2 {
            let (delta_coef, dkcoef) = num_two_pow(d1,d2,denk2,m1,m2, pdg1, pdg2, cache);
			if q2exp == 0 {
				term = numk2*delta_coef + k1dotk2(numk2,denk2,ksum2).pow(2)/denk2*dkcoef;
            } else if q2exp == 1 {
				let (delta_coef2, dkcoef2) = num_two_pow(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache);
				term = -m1.clone()*(numk2.clone()*delta_coef + k1dotk2(numk2,denk2,ksum2).pow(2)/(denk2.clone())*dkcoef);
				term += numk2.clone()*delta_coef2 + k1dotk2(numk2,denk2,ksum2).pow(2)/(denk2.clone())*dkcoef2;
            } else if q2exp == 2 {
				let (delta_coef2, dkcoef2) = num_two_pow(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache);
				let (delta_coef3, dkcoef3) = num_two_pow(d1-2,d2,denk2,m1,m2, pdg1, pdg2, cache);
				term = numk2.clone()*delta_coef3 + k1dotk2(numk2,denk2,ksum2).pow(2)/(denk2.clone())*dkcoef3;
				term += -2*m1.clone()*(numk2.clone()*delta_coef2 + k1dotk2(numk2,denk2,ksum2).pow(2)/(denk2.clone())*dkcoef2);
				term += m1.clone().pow(2)*(numk2.clone()*delta_coef + k1dotk2(numk2,denk2,ksum2).pow(2)/denk2.clone()*dkcoef);
            } else {
                println!("ERROR: exceeded calculable power in tri_dim");
            }
        } else if kqexp == 3 {
            let (delta_coef, dkcoef) = num_three_pow(d1, d2, denk2, m1, m2, pdg1, pdg2, cache);
            if q2exp == 0 {
                term = numk2.clone()*delta_coef*k1dotk2(numk2,denk2,ksum2)/(denk2.clone().sqrt()) 
                     + dkcoef*k1dotk2(numk2,denk2,ksum2).pow(3)/(denk2.clone()*denk2.clone().sqrt());
            } else if q2exp == 1 {
                let (delta_coef2, dkcoef2) = num_three_pow(d1 - 1, d2, denk2, m1, m2, pdg1, pdg2, cache);
                term = numk2.clone()*delta_coef2*k1dotk2(numk2,denk2,ksum2)/(denk2.clone().sqrt()) 
                     + dkcoef2*k1dotk2(numk2,denk2,ksum2).pow(3)/(denk2.clone()*denk2.clone().sqrt());
                term += -m1.clone()*(numk2.clone()*delta_coef*k1dotk2(numk2,denk2,ksum2)/(denk2.clone().sqrt()) 
                        + dkcoef*k1dotk2(numk2,denk2,ksum2).pow(3)/(denk2.clone()*denk2.clone().sqrt()));
            } else {
                println!("ERROR: exceeded calculable power in tri_dim");
            }
        } else if kqexp == 4 {
            if q2exp == 0 {
                let (coef1, coef2, coef3) = num_four_pow(d1, d2, denk2, m1, m2, pdg1, pdg2, cache);
                term = coef1*numk2.clone().pow(2) + numk2.clone()*k1dotk2(numk2,denk2,ksum2).pow(2)*coef2/denk2.clone() 
                     + coef3*k1dotk2(numk2,denk2,ksum2).pow(4)/denk2.clone().pow(2)
            } else {
                println!("ERROR: exceeded calculable power in tri_dim");
            }
        } else {
            println!("ERROR: exceeded calculable power in tri_dim");
        }
        res += term * (term_list[i] as i64) * numk2.clone().pow(k2exp as u64);
        // println!("{:?}",res);
    }
    return res
}

fn num_one_pow(d1: isize, d2: isize, denk2: &Float, m1: &Complex, m2:&Complex, pdg1: usize, pdg2: usize, cache: &mut TopologyCache) -> Complex {
    return (bubn(d1,d2-1,denk2,m1,m2, pdg1, pdg2, cache) 
            - bubn(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache) 
            - (denk2.clone() + m2.clone() - m1)*bubn(d1,d2,denk2,m1,m2, pdg1, pdg2, cache))/2/denk2
}

fn num_two_pow(d1: isize, d2: isize, denk2: &Float, m1: &Complex, m2:&Complex, pdg1: usize, pdg2: usize, cache: &mut TopologyCache) -> (Complex, Complex) {
    let coef1: Complex = (bubn(d1,d2-2,denk2,m1,m2, pdg1, pdg2, cache) 
                         - 2*(denk2.clone() + m2.clone() - m1.clone())*bubn(d1,d2-1,denk2,m1,m2, pdg1, pdg2, cache) 
                        + (denk2.clone() + m2.clone() - m1.clone()).pow(2)*bubn(d1,d2,denk2,m1,m2, pdg1, pdg2, cache)
                        -2*bubn(d1-1,d2-1,denk2,m1,m2, pdg1, pdg2, cache) + 2*(denk2.clone() + m2.clone() - m1.clone())*bubn(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache) 
                        + bubn(d1-2,d2,denk2,m1,m2, pdg1, pdg2, cache))*(-1)/(8*denk2.clone()) 
                        + bubn(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache)/2 - m1.clone()*bubn(d1,d2,denk2,m1,m2, pdg1, pdg2, cache)/2;
    let coef2: Complex = bubn(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache) - m1.clone()*bubn(d1,d2,denk2,m1,m2, pdg1, pdg2, cache) - 3 * coef1.clone();
return (coef1, coef2)
}

fn num_three_pow(d1: isize, d2: isize, denk2: &Float, m1: &Complex, m2:&Complex, pdg1: usize, pdg2: usize, cache: &mut TopologyCache) -> (Complex, Complex) {
    let coef1: Complex =  3./(16 * denk2.clone() * denk2.clone().sqrt())*(bubn(d1-3,d2,denk2,m1,m2, pdg1, pdg2, cache) - 3*bubn(d1-2,d2-1,denk2,m1,m2, pdg1, pdg2, cache) - 4 * denk2.clone() * bubn(d1-2,d2,denk2,m1,m2, pdg1, pdg2, cache)
                        + 3*(denk2.clone() + m2.clone() - m1.clone())*bubn(d1-2,d2,denk2,m1,m2, pdg1, pdg2, cache) + 3*bubn(d1-1,d2-2,denk2,m1,m2, pdg1, pdg2, cache)
                        + 4*denk2.clone()*bubn(d1-1,d2-1,denk2,m1,m2, pdg1, pdg2, cache) - 6*(denk2.clone() + m2.clone() - m1.clone())*bubn(d1-1,d2-1,denk2,m1,m2, pdg1, pdg2, cache)
                        - 4*denk2.clone()*(denk2.clone() + m2.clone() - m1.clone())*bubn(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache) + 3*(denk2.clone() + m2.clone() - m1.clone()).pow(2)*bubn(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache)
                        + 4* denk2.clone()*m1.clone()*bubn(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache) - bubn(d1,d2-3,denk2,m1,m2, pdg1, pdg2, cache) + 3*(denk2.clone() + m2.clone() - m1.clone())*bubn(d1,d2-2,denk2,m1,m2, pdg1, pdg2, cache)
                        - 3*(denk2.clone() + m2.clone() - m1.clone()).pow(2)*bubn(d1,d2-1,denk2,m1,m2, pdg1, pdg2, cache) - 4*denk2.clone()*m1.clone()*bubn(d1,d2-1,denk2,m1,m2, pdg1, pdg2, cache)
                        + (denk2.clone() + m2.clone() - m1.clone()).pow(3)*bubn(d1,d2,denk2,m1,m2, pdg1, pdg2, cache) + 4*denk2.clone()*(denk2.clone() + m2.clone() - m1.clone())*m1.clone()*bubn(d1,d2,denk2,m1,m2, pdg1, pdg2, cache));

    let coef2: Complex = 1./(2*denk2.clone().sqrt())*(bubn(d1-1,d2-1,denk2,m1,m2, pdg1, pdg2, cache) - bubn(d1-2,d2,denk2,m1,m2, pdg1, pdg2, cache)
    -(denk2.clone() + m2.clone() - m1)*bubn(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache) + m1.clone()*bubn(d1-1,d2,denk2,m1,m2, pdg1, pdg2, cache) - m1.clone()*bubn(d1,d2-1,denk2,m1,m2, pdg1, pdg2, cache)
    +(denk2.clone() + m2.clone() - m1)*m1.clone()*bubn(d1,d2,denk2,m1,m2, pdg1, pdg2, cache))-5*coef1.clone()/3;
    
    return (coef1, coef2)
}

fn num_four_pow(d1: isize, d2: isize, denk2: &Float, m1: &Complex, m2:&Complex, pdg1: usize, pdg2: usize, cache: &mut TopologyCache) -> (Complex, Complex, Complex) {
     let coef1: Complex = 3*(bubn(-4 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*bubn(-3 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*denk2.clone()*bubn(-3 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*m1.clone()*bubn(-3 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                4*m2.clone()*bubn(-3 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                6*bubn(-2 + d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                4*denk2.clone()*bubn(-2 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                12*m1.clone()*bubn(-2 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                12*m2.clone()*bubn(-2 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                6*denk2.clone().pow(2)*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                12*denk2.clone()*m1.clone()*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                6*m1.clone().pow(2)*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*denk2.clone()*m2.clone()*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                12*m1.clone()*m2.clone()*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                6*m2.clone().pow(2)*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*bubn(-1 + d1, -3 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                4*denk2.clone()*bubn(-1 + d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                12*m1.clone()*bubn(-1 + d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                12*m2.clone()*bubn(-1 + d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                4*denk2.clone().pow(2)*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                8*denk2.clone()*m1.clone()*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                12*m1.clone().pow(2)*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                8*denk2.clone()*m2.clone()*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                24*m1.clone()*m2.clone()*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                12*m2.clone().pow(2)*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*denk2.clone().pow(3)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                12*denk2.clone().pow(2)*m1.clone()*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                12*denk2.clone()*m1.clone().pow(2)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*m1.clone().pow(3)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*denk2.clone().pow(2)*m2.clone()*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                8*denk2.clone()*m1.clone()*m2.clone()*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                12*m1.clone().pow(2)*m2.clone()*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                4*denk2.clone()*m2.clone().pow(2)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                12*m1.clone()*m2.clone().pow(2)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                4*m2.clone().pow(3)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                bubn(d1, -4 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*denk2.clone()*bubn(d1, -3 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                4*m1.clone()*bubn(d1, -3 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*m2.clone()*bubn(d1, -3 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                6*denk2.clone().pow(2)*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*denk2.clone()*m1.clone()*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                6*m1.clone().pow(2)*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                12*denk2.clone()*m2.clone()*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                12*m1.clone()*m2.clone()*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                6*m2.clone().pow(2)*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*denk2.clone().pow(3)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*denk2.clone().pow(2)*m1.clone()*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                4*denk2.clone()*m1.clone().pow(2)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                4*m1.clone().pow(3)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                12*denk2.clone().pow(2)*m2.clone()*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                8*denk2.clone()*m1.clone()*m2.clone()*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                12*m1.clone().pow(2)*m2.clone()*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                12*denk2.clone()*m2.clone().pow(2)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                12*m1.clone()*m2.clone().pow(2)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*m2.clone().pow(3)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                denk2.clone().pow(4)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                4*denk2.clone().pow(3)*m1.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                6*denk2.clone().pow(2)*m1.clone().pow(2)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                4*denk2.clone()*m1.clone().pow(3)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                m1.clone().pow(4)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                4*denk2.clone().pow(3)*m2.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                4*denk2.clone().pow(2)*m1.clone()*m2.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*denk2.clone()*m1.clone().pow(2)*m2.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*m1.clone().pow(3)*m2.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                6*denk2.clone().pow(2)*m2.clone().pow(2)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*denk2.clone()*m1.clone()*m2.clone().pow(2)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                6*m1.clone().pow(2)*m2.clone().pow(2)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                4*denk2.clone()*m2.clone().pow(3)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                4*m1.clone()*m2.clone().pow(3)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                m2.clone().pow(4)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache))/(128*denk2.clone().pow(2));
    let coef2: Complex = -3*(5*bubn(-4 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            20*bubn(-3 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            4*denk2.clone()*bubn(-3 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            20*m1.clone()*bubn(-3 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            20*m2.clone()*bubn(-3 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            30*bubn(-2 + d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            12*denk2.clone()*bubn(-2 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            60*m1.clone()*bubn(-2 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            60*m2.clone()*bubn(-2 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            2*denk2.clone().pow(2)*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            12*denk2.clone()*m1.clone()*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            30*m1.clone().pow(2)*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            12*denk2.clone()*m2.clone()*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            60*m1.clone()*m2.clone()*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            30*m2.clone().pow(2)*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            20*bubn(-1 + d1, -3 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            36*denk2.clone()*bubn(-1 + d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            60*m1.clone()*bubn(-1 + d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            60*m2.clone()*bubn(-1 + d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            12*denk2.clone().pow(2)*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            24*denk2.clone()*m1.clone()*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            60*m1.clone().pow(2)*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            72*denk2.clone()*m2.clone()*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            120*m1.clone()*m2.clone()*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            60*m2.clone().pow(2)*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            4*denk2.clone().pow(3)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            4*denk2.clone().pow(2)*m1.clone()*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            12*denk2.clone()*m1.clone().pow(2)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            20*m1.clone().pow(3)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            12*denk2.clone().pow(2)*m2.clone()*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            24*denk2.clone()*m1.clone()*m2.clone()*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            60*m1.clone().pow(2)*m2.clone()*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            36*denk2.clone()*m2.clone().pow(2)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            60*m1.clone()*m2.clone().pow(2)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            20*m2.clone().pow(3)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            5*bubn(d1, -4 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            20*denk2.clone()*bubn(d1, -3 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            20*m1.clone()*bubn(d1, -3 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            20*m2.clone()*bubn(d1, -3 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            30*denk2.clone().pow(2)*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            36*denk2.clone()*m1.clone()*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            30*m1.clone().pow(2)*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            60*denk2.clone()*m2.clone()*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            60*m1.clone()*m2.clone()*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            30*m2.clone().pow(2)*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            20*denk2.clone().pow(3)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            12*denk2.clone().pow(2)*m1.clone()*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            12*denk2.clone()*m1.clone().pow(2)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            20*m1.clone().pow(3)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            60*denk2.clone().pow(2)*m2.clone()*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            72*denk2.clone()*m1.clone()*m2.clone()*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            60*m1.clone().pow(2)*m2.clone()*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            60*denk2.clone()*m2.clone().pow(2)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            60*m1.clone()*m2.clone().pow(2)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            20*m2.clone().pow(3)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            5*denk2.clone().pow(4)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            4*denk2.clone().pow(3)*m1.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            2*denk2.clone().pow(2)*m1.clone().pow(2)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            4*denk2.clone()*m1.clone().pow(3)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            5*m1.clone().pow(4)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            20*denk2.clone().pow(3)*m2.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            12*denk2.clone().pow(2)*m1.clone()*m2.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            12*denk2.clone()*m1.clone().pow(2)*m2.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            20*m1.clone().pow(3)*m2.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            30*denk2.clone().pow(2)*m2.clone().pow(2)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            36*denk2.clone()*m1.clone()*m2.clone().pow(2)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            30*m1.clone().pow(2)*m2.clone().pow(2)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            20*denk2.clone()*m2.clone().pow(3)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
            20*m1.clone()*m2.clone().pow(3)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
            5*m2.clone().pow(4)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache))/(64*denk2.clone().pow(2));
    let coef3: Complex = (-1)*(-35*bubn(-4 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                140*bubn(-3 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                20*denk2.clone()*bubn(-3 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                140*m1.clone()*bubn(-3 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                140*m2.clone()*bubn(-3 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                210*bubn(-2 + d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                180*denk2.clone()*bubn(-2 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                420*m1.clone()*bubn(-2 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                420*m2.clone()*bubn(-2 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                18*denk2.clone().pow(2)*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                60*denk2.clone()*m1.clone()*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                210*m1.clone().pow(2)*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                180*denk2.clone()*m2.clone()*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                420*m1.clone()*m2.clone()*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                210*m2.clone().pow(2)*bubn(-2 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                140*bubn(-1 + d1, -3 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                300*denk2.clone()*bubn(-1 + d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                420*m1.clone()*bubn(-1 + d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                420*m2.clone()*bubn(-1 + d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                180*denk2.clone().pow(2)*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                360*denk2.clone()*m1.clone()*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                420*m1.clone().pow(2)*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                600*denk2.clone()*m2.clone()*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                840*m1.clone()*m2.clone()*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                420*m2.clone().pow(2)*bubn(-1 + d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                20*denk2.clone().pow(3)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                36*denk2.clone().pow(2)*m1.clone()*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                60*denk2.clone()*m1.clone().pow(2)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                140*m1.clone().pow(3)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                180*denk2.clone().pow(2)*m2.clone()*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                360*denk2.clone()*m1.clone()*m2.clone()*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                420*m1.clone().pow(2)*m2.clone()*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                300*denk2.clone()*m2.clone().pow(2)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                420*m1.clone()*m2.clone().pow(2)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                140*m2.clone().pow(3)*bubn(-1 + d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                35*bubn(d1, -4 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                140*denk2.clone()*bubn(d1, -3 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                140*m1.clone()*bubn(d1, -3 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                140*m2.clone()*bubn(d1, -3 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                210*denk2.clone().pow(2)*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                300*denk2.clone()*m1.clone()*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                210*m1.clone().pow(2)*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                420*denk2.clone()*m2.clone()*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                420*m1.clone()*m2.clone()*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                210*m2.clone().pow(2)*bubn(d1, -2 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                140*denk2.clone().pow(3)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                180*denk2.clone().pow(2)*m1.clone()*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                180*denk2.clone()*m1.clone().pow(2)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                140*m1.clone().pow(3)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                420*denk2.clone().pow(2)*m2.clone()*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                600*denk2.clone()*m1.clone()*m2.clone()*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                420*m1.clone().pow(2)*m2.clone()*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                420*denk2.clone()*m2.clone().pow(2)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                420*m1.clone()*m2.clone().pow(2)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                140*m2.clone().pow(3)*bubn(d1, -1 + d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                35*denk2.clone().pow(4)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                20*denk2.clone().pow(3)*m1.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                18*denk2.clone().pow(2)*m1.clone().pow(2)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                20*denk2.clone()*m1.clone().pow(3)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                35*m1.clone().pow(4)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                140*denk2.clone().pow(3)*m2.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                180*denk2.clone().pow(2)*m1.clone()*m2.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                180*denk2.clone()*m1.clone().pow(2)*m2.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                140*m1.clone().pow(3)*m2.clone()*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                210*denk2.clone().pow(2)*m2.clone().pow(2)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                300*denk2.clone()*m1.clone()*m2.clone().pow(2)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                210*m1.clone().pow(2)*m2.clone().pow(2)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                140*denk2.clone()*m2.clone().pow(3)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) + 
                140*m1.clone()*m2.clone().pow(3)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache) - 
                35*m2.clone().pow(4)*bubn(d1, d2, denk2, m1, m2, pdg1, pdg2, cache))/(128*denk2.clone().pow(2));

return (coef1, coef2, coef3)
}

fn tri_dim_two(n1: isize, n2: isize, d1: isize, numk21: &Float, numk22: &Float, ksum2: &Float,
    dm: &Complex, _pdg1: usize, _cache: &mut TopologyCache) -> Complex {
    // integral of (k1 - q)^2^n1 (k2 + q)^2^n2/(q2+dm)^d1

	// term_list1 are the coefficients of (k1 - q)^2^n1 corresponding to the exponents in exp_list   
	// exp_list1 are the exponents of (k1 - q)^2^n1 of the form k1^2^k2exp1*q^2^q2exp1*(k.q)^kqexp1, 
	// written as (k2exp1, q2exp1, kqexp1) 

	// term_list2 are the coefficients of (k2 + q)^2^n2 corresponding to the exponents in exp_list   
	// exp_list2 are the exponents of (k2 + q)^2^n2 of the form k2^2^k2exp2*q^2^q2exp2*(k.q)^kqexp2, 
	// written as (k2exp2, q2exp2, kqexp2) 
    let list_length_n1: usize = ((1 + n1)*(2 + n1)/2) as usize;
    let list_length_n2: usize = ((1 + n2)*(2 + n2)/2) as usize;

    let (term_list1, exp_list1) = num_terms(n1 as usize, true);
    let (term_list2, exp_list2) = num_terms(n2 as usize, false);

    let mut res = Complex::new(PREC);

    for i1 in 0..list_length_n1 {
        for i2 in 0..list_length_n2 {
            let k2exp1 = exp_list1[i1][0] as i64;
			let k2exp2 = exp_list2[i2][0] as i64;
			let q2exp = (exp_list1[i1][1] + exp_list2[i2][1]) as isize;
			let kqexp1 = exp_list1[i1][2] as i64;
			let kqexp2 = exp_list2[i2][2] as i64;
			let kqexp = kqexp1 + kqexp2;
            if kqexp%2 == 0 {
                let mut term = Complex::new(PREC);
                if kqexp != 0 {
                    // cases where kqexp == 2
                    if kqexp1 == 2 && kqexp2 == 0 {
                        term = dim_gen(q2exp+1,d1,dm)*numk21.clone()/3;
                    } else if kqexp1 == 0 && kqexp2 == 2 {
                        term = dim_gen(q2exp+1,d1,dm)*numk22.clone()/3;
                    } else if kqexp1 == 1 && kqexp2 == 1 {
                        term = dim_gen(q2exp+1,d1,dm)*(k1dotk2(numk21,numk22,ksum2))/3;
                    }    
                    // cases where kqexp == 4
                    else if kqexp1 == 0 && kqexp2 == 4 {
						term = dim_gen(q2exp+2,d1,dm)*(numk22.clone().pow(2))/5
					} else if kqexp1 == 4 && kqexp2 == 0 {
						term = dim_gen(q2exp+2,d1,dm)*(numk21.clone().pow(2))/5
					} else if kqexp1 == 1 && kqexp2 == 3 {
						term = dim_gen(q2exp+2,d1,dm)*(k1dotk2(numk21,numk22,ksum2)*numk22.clone())/5
					} else if kqexp1 == 3 && kqexp2 == 1 {
						term = dim_gen(q2exp+2,d1,dm)*(k1dotk2(numk21,numk22,ksum2)*numk21.clone())/5
					} else if kqexp1 == 2 && kqexp2 == 2 {
						term = dim_gen(q2exp+2,d1,dm)*(numk21*numk22 + 2*k1dotk2(numk21,numk22,ksum2).pow(2))/15
                    }
                    // cases where kqexp == 6
                    else if kqexp1 == 6 && kqexp2 == 0 {
						term = dim_gen(q2exp + 3, d1, dm)*numk21.clone().pow(3)/7
					} else if kqexp1 == 0 && kqexp2 == 6 {
						term = dim_gen(q2exp + 3, d1, dm)*numk22.clone().pow(3)/7
					} else if kqexp1 == 5 && kqexp2 == 1 {
						term = dim_gen(q2exp + 3, d1, dm)*numk21.clone().pow(2)*k1dotk2(numk21,numk22,ksum2)/7
					} else if kqexp1 == 1 && kqexp2 == 5 {
						term = dim_gen(q2exp + 3, d1, dm)*numk22.clone().pow(2)*k1dotk2(numk21,numk22,ksum2)/7
					} else if kqexp1 == 4 && kqexp2 == 2 {
						term = dim_gen(q2exp + 3, d1, dm)*(numk21.clone().pow(2)*numk22.clone() + 4*k1dotk2(numk21,numk22,ksum2).pow(2)*numk21)/35
					} else if kqexp1 == 3 && kqexp2 == 3 {
						term = dim_gen(q2exp + 3, d1, dm)*(3*numk21.clone()*numk22.clone()*k1dotk2(numk21,numk22,ksum2) + 2*k1dotk2(numk21,numk22,ksum2).pow(3))/35
					} else if kqexp1 == 2 && kqexp2 == 4 {
						term = dim_gen(q2exp + 3,d1,dm)*(numk22.clone().pow(2)*numk21.clone() + 4*k1dotk2(numk21,numk22,ksum2).pow(2)*numk22.clone())/35
                    }
                    // cases where kqexp == 8
                    else if kqexp1 == 4 && kqexp2 == 4 {
                        term = dim_gen(q2exp + 4, d1, dm)*(3*numk21.clone().pow(2)*numk22.clone().pow(2) + 24*numk21.clone()*numk22.clone()*k1dotk2(numk21,numk22,ksum2).pow(2) + 8*k1dotk2(numk21,numk22,ksum2).pow(4))/315
                    } else {
                        println!("ERROR: case not considered in tri_dim_two");
                    }
                } else {
                    // case where kqexp == 0
                    term = dim_gen(q2exp,d1,dm);
                }                
                res += term * (term_list1[i1] * term_list2[i2]) as i64 * numk21.clone().pow(k2exp1) * numk22.clone().pow(k2exp2)
            }       
        }
    }
   
    return res
}

fn ltrian(n1: isize, d1: isize, 
          n2: isize, d2: isize, 
          n3: isize, d3: isize,
          k21: &Float, k22: &Float, k23: &Float,
          m1: &Complex, m2: &Complex, m3: &Complex,
          pdg1: usize, pdg2: usize, pdg3: usize, 
          cache: &mut TopologyCache) -> Complex {

    if n1 == 0 && n2 == 0 && n3 == 0 {
        return trian(d1,d2,d3,k21,k22,k23,m1,m2,m3, pdg1, pdg2, pdg3, cache)
    }
    if d1 == 0 && n1 != 0 {
        return ltrian(0,-n1,n2,d2,n3,d3,k21,k22,k23, &Complex::new(PREC), m2, m3, 0, pdg2, pdg3, cache)
    }
    if d2 == 0 && n2 != 0 {
        return ltrian(n1,d1,0,-n2,n3,d3,k21,k22,k23, m1, &Complex::new(PREC), m3, pdg1, 0, pdg3, cache)
    }
    if d3 == 0 && n3 != 0 {
        return ltrian(n1,d1,n2,d2,0,-n3,k21,k22,k23, m1, m2, &Complex::new(PREC), pdg1, pdg2, 0, cache)
    }
    if n1 > 0 {
        return ltrian(n1-1,d1-1,n2,d2,n3,d3,k21,k22,k23,m1,m2,m3, pdg1, pdg2, pdg3, cache) 
                - m1*ltrian(n1-1,d1,n2,d2,n3,d3,k21,k22,k23,m1,m2,m3, pdg1, pdg2, pdg3, cache)
    }
    if n2 > 0 {
        return ltrian(n1,d1,n2-1,d2-1,n3,d3,k21,k22,k23,m1,m2,m3, pdg1, pdg2, pdg3, cache) 
                - m2*ltrian(n1,d1,n2-1,d2,n3,d3,k21,k22,k23,m1,m2,m3, pdg1, pdg2, pdg3, cache)
    }
    if n3 > 0 {
        return ltrian(n1,d1,n2,d2,n3-1,d3-1,k21,k22,k23,m1,m2,m3, pdg1, pdg2, pdg3, cache) 
                - m3*ltrian(n1,d1,n2,d2,n3-1,d3,k21,k22,k23,m1,m2,m3, pdg1, pdg2, pdg3, cache)
    }
    if n1 < 0 {
        return (ltrian(n1,d1-1,n2,d2,n3,d3,k21,k22,k23,m1,m2,m3, pdg1, pdg2, pdg3, cache) 
                - ltrian(n1+1,d1,n2,d2,n3,d3,k21,k22,k23,m1,m2,m3, pdg1, pdg2, pdg3, cache))/m1
    }
    if n2 < 0 {
        return (ltrian(n1,d1,n2,d2-1,n3,d3,k21,k22,k23,m1,m2,m3, pdg1, pdg2, pdg3, cache) 
                - ltrian(n1,d1,n2+1,d2,n3,d3,k21,k22,k23,m1,m2,m3, pdg1, pdg2, pdg3, cache))/m2
    }
    if n3 < 0 {
        return (ltrian(n1,d1,n2,d2,n3,d3-1,k21,k22,k23,m1,m2,m3, pdg1, pdg2, pdg3, cache) 
                - ltrian(n1,d1,n2,d2,n3+1,d3,k21,k22,k23,m1,m2,m3, pdg1, pdg2, pdg3, cache))/m3
    }
    println!("Error: case not considered in Ltrian");
    return Complex::new(PREC)
}