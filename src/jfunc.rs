use lazy_static::lazy_static;
use rug::{Complex, Float, ops::Pow};
use ndarray::{Array1,s};

pub mod ltrian;

const PREC: u32 = 190;
const PI: f64 = 3.14159265358979323846264338327950288f64;
const SQRT_PI: f64 = 1.77245385090551602729816748334114518279754945612238f64;


lazy_static! {

    static ref KUV1: Complex = Complex::with_val(PREC,(0.,0.069));
    static ref KUV2: Complex = Complex::with_val(PREC,(0.,0.0082));
    static ref KUV3: Complex = Complex::with_val(PREC,(0.,0.0013));
    static ref KUV4: Complex = Complex::with_val(PREC,(0.,0.0000135));

    // static ref kpeak1: Complex = Complex::with_val(PREC,(0.034,0.));
    // static ref kpeak2: Complex = Complex::with_val(PREC,(0.001,0.));
    // static ref kpeak3: Complex = Complex::with_val(PREC,(0.000076,0.));
    // static ref kpeak4: Complex = Complex::with_val(PREC,(0.0000156,0.));

    static ref M1: Complex = Complex::with_val(PREC,(0.034,-0.069));
    static ref M1CONJ: Complex = Complex::with_val(PREC,(0.034,0.069));
    static ref M2: Complex = Complex::with_val(PREC,(0.001,-0.0082));
    static ref M2CONJ: Complex = Complex::with_val(PREC,(0.001,0.0082));
    static ref M3: Complex = Complex::with_val(PREC,(0.000076,-0.0013));
    static ref M3CONJ: Complex = Complex::with_val(PREC,(0.000076,0.0013));
    static ref M4: Complex = Complex::with_val(PREC,(0.0000156,-0.0000135));
    static ref M4CONJ: Complex = Complex::with_val(PREC,(0.0000156,0.0000135));

    static ref FBABISMASSES: Vec<Complex> = vec![Complex::new(PREC),
        M2.clone(),M2CONJ.clone(),
        M3.clone(),M3CONJ.clone(),M3.clone(),M3CONJ.clone(),
        M4.clone(),M4CONJ.clone(),
        M1.clone(),M1CONJ.clone(),
        M2.clone(),M2CONJ.clone(),
        M3.clone(),M3CONJ.clone(),
        M4.clone(),M4CONJ.clone(),
        M1.clone(),M1CONJ.clone(),
        M2.clone(),M2CONJ.clone(),
        M3.clone(),M3CONJ.clone(),
        M4.clone(),M4CONJ.clone(),
        M1.clone(),M1CONJ.clone(),
        M2.clone(),M2CONJ.clone(),
        M3.clone(),M3CONJ.clone(),
        M4.clone(),M4CONJ.clone()];

    static ref FBABISMASSIND: [usize;33] = [0,2,6,3,7,3,7,4,8,1,5,2,6,3,7,4,8,1,5,2,6,3,7,4,8,1,5,2,6,3,7,4,8];

    static ref FBABISEXP: Vec<[isize; 2]> = vec![[0,0],[1,1],[1,1],
    [0,1],[0,1],[0,2],[0,2],[0,1],[0,1],[0,1],[0,1],[1,2],[1,2],[0,3],[0,3],
    [0,2],[0,2],[0,2],[0,2],[1,3],[1,3],[0,4],[0,4],[0,3],[0,3],[0,3],[0,3],
    [1,4],[1,4],[0,5],[0,5],[0,4],[0,4]];

    static ref DBASIS1: Array1<isize> = Array1::from_vec(vec![0]);
    static ref DBASIS2: Array1<isize> = Array1::from_vec(vec![1,2]);
    static ref DBASIS3: Array1<isize> = Array1::from_vec(vec![3,4,5,6]);
    static ref DBASIS4: Array1<isize> = Array1::from_vec(vec![7,8]);
    static ref DBASIS5: Array1<isize> = Array1::from_vec(vec![9,10]);
    static ref DBASIS6: Array1<isize> = Array1::from_vec(vec![1,2,11,12]);   
    static ref DBASIS7: Array1<isize> = Array1::from_vec(vec![3,4,5,6,13,14]);
    static ref DBASIS8: Array1<isize> = Array1::from_vec(vec![7,8,15,16]);
    static ref DBASIS9: Array1<isize> = Array1::from_vec(vec![9,10,17,18]);
    static ref DBASIS10: Array1<isize> = Array1::from_vec(vec![1,2,11,12,19,20]);
    static ref DBASIS11: Array1<isize> = Array1::from_vec(vec![3,4,5,6,13,14,21,22]);
    static ref DBASIS12: Array1<isize> = Array1::from_vec(vec![7,8,15,16,23,24]);
    static ref DBASIS13: Array1<isize> = Array1::from_vec(vec![9,10,17,18,25,26]);
    static ref DBASIS14: Array1<isize> = Array1::from_vec(vec![1,2,11,12,19,20,27,28]);
    static ref DBASIS15: Array1<isize> = Array1::from_vec(vec![3,4,5,6,13,14,21,22,29,30]);
    static ref DBASIS16: Array1<isize> = Array1::from_vec(vec![7,8,15,16,23,24,31,32]);

    static ref MATCOEF1: Array1<Complex> = Array1::from_vec(vec![Complex::with_val(PREC,1.0)]);
    static ref MATCOEF2: Array1<Complex> = Array1::from_vec(vec![-200.*KUV2.clone(), 200.*KUV2.clone()]);
    static ref MATCOEF3: Array1<Complex> = Array1::from_vec(vec![-KUV3.clone()/4., KUV3.clone()/4.,KUV3.clone().pow(2i32)/4.,KUV3.clone().pow(2i32)/4.]);
    static ref MATCOEF4: Array1<Complex> = Array1::from_vec(vec![-KUV4.clone()/2., KUV4.clone()/2.]);
    static ref MATCOEF5: Array1<Complex> = Array1::from_vec(vec![-KUV1.clone()/2., KUV1.clone()/2.]);
    static ref MATCOEF6: Array1<Complex> = Array1::from_vec(vec![-100.*KUV2.clone(),100.*KUV2.clone(),100.*KUV2.clone().pow(2i32),100.*KUV2.clone().pow(2i32)]);
    static ref MATCOEF7: Array1<Complex> = Array1::from_vec(vec![-KUV3.clone()*3./16.,KUV3.clone()*3./16.,KUV3.clone().pow(2i32)*3./16.,KUV3.clone().pow(2i32)*3./16.,-KUV3.clone().pow(3 as i32)/8.,KUV3.clone().pow(3 as i32)/8.]);
    static ref MATCOEF8: Array1<Complex> = Array1::from_vec(vec![-KUV4.clone()/4.,KUV4.clone()/4.,KUV4.clone().pow(2i32)/4.,KUV4.clone().pow(2i32)/4.]);
    static ref MATCOEF9: Array1<Complex> = Array1::from_vec(vec![-KUV1.clone()/4.,KUV1.clone()/4.,KUV1.clone().pow(2 as i32)/4.,KUV1.clone().pow(2 as i32)/4.]);
    static ref MATCOEF10: Array1<Complex> = Array1::from_vec(vec![-75.*KUV2.clone(),75.*KUV2.clone(),75.*KUV2.clone().pow(2 as i32),75.*KUV2.clone().pow(2 as i32),-50.*KUV2.clone().pow(3 as i32),50.*KUV2.clone().pow(3 as i32)]);
    static ref MATCOEF11: Array1<Complex> = Array1::from_vec(vec![-KUV3.clone()*5./32.,KUV3.clone()*5./32.,KUV3.clone().pow(2 as i32)*5./32.,KUV3.clone().pow(2 as i32)*5./32.,-KUV3.clone().pow(3 as i32)/8.,KUV3.clone().pow(3 as i32)/8.,KUV3.clone().pow(4 as i32)/16.,KUV3.clone().pow(4 as i32)/16.]);
    static ref MATCOEF12: Array1<Complex> = Array1::from_vec(vec![-KUV4.clone()*3./16.,KUV4.clone()*3./16.,KUV4.clone().pow(2 as i32)*3./16.,KUV4.clone().pow(2 as i32)*3./16.,-KUV4.clone().pow(3 as i32)/8.,KUV4.clone().pow(3 as i32)/8.]);
    static ref MATCOEF13: Array1<Complex> = Array1::from_vec(vec![-KUV1.clone()*3./16.,KUV1.clone()*3./16.,KUV1.clone().pow(2 as i32)*3./16.,KUV1.clone().pow(2 as i32)*3./16.,-KUV1.clone().pow(3 as i32)/8.,KUV1.clone().pow(3 as i32)/8.]);
    static ref MATCOEF14: Array1<Complex> = Array1::from_vec(vec![-KUV2.clone()*125./2.,KUV2.clone()*125./2.,KUV2.clone().pow(2 as i32)*125./2.,KUV2.clone().pow(2 as i32)*125./2.,-50.*KUV2.clone().pow(3 as i32),50.*KUV2.clone().pow(3 as i32),25.*KUV2.clone().pow(4 as i32),25.*KUV2.clone().pow(4 as i32)]);
    static ref MATCOEF15: Array1<Complex> = Array1::from_vec(vec![-KUV3.clone()*35./256.,KUV3.clone()*35./256.,KUV3.clone().pow(2 as i32)*35./256.,KUV3.clone().pow(2 as i32)*35./256.,-KUV3.clone().pow(3 as i32)*15./128.,KUV3.clone().pow(3 as i32)*15./128.,KUV3.clone().pow(4 as i32)*5./64.,KUV3.clone().pow(4 as i32)*5./64.,-KUV3.clone().pow(5 as i32)/32.,KUV3.clone().pow(5 as i32)/32.]);
    static ref MATCOEF16: Array1<Complex> = Array1::from_vec(vec![-KUV4.clone()*5./32.,KUV4.clone()*5./32.,KUV4.clone().pow(2 as i32)*5./32.,KUV4.clone().pow(2 as i32)*5./32.,-KUV4.clone().pow(3 as i32)/8.,KUV4.clone().pow(3 as i32)/8.,KUV4.clone().pow(4 as i32)/16.,KUV4.clone().pow(4 as i32)/16.]);
    


    // static ref cache = ltrian::ltrian::TopologyCache {
    //     // tadpole_cache: HashMap::new(),
    //     bubble_cache: HashMap::new(),
    //     trian_cache: HashMap::new(),
    // };

}

fn dbabis_coef(d: isize) -> (Array1<isize>, Array1<Complex>) {
    if d == 0 {
        return (DBASIS1.clone(), MATCOEF1.clone());
    } else if d == 1 {
        return (DBASIS2.clone(), MATCOEF2.clone());
    } else if d == 2 {
        return (DBASIS3.clone(), MATCOEF3.clone());
    } else if d == 3 {
        return (DBASIS4.clone(), MATCOEF4.clone());
    } else if d == 4 {
        return (DBASIS5.clone(), MATCOEF5.clone());
    } else if d == 5 {
        return (DBASIS6.clone(), MATCOEF6.clone());
    } else if d == 6 {
        return (DBASIS7.clone(), MATCOEF7.clone());
    } else if d == 7 {
        return (DBASIS8.clone(), MATCOEF8.clone());
    } else if d == 8 {
        return (DBASIS9.clone(), MATCOEF9.clone());
    } else if d == 9 {
        return (DBASIS10.clone(), MATCOEF10.clone());
    } else if d == 10 {
        return (DBASIS11.clone(), MATCOEF11.clone());
    } else if d == 11 {
        return (DBASIS12.clone(), MATCOEF12.clone());
    } else if d == 12 {
        return (DBASIS13.clone(), MATCOEF13.clone());
    } else if d == 13 {
        return (DBASIS14.clone(), MATCOEF14.clone());
    } else if d == 14 {
        return (DBASIS15.clone(), MATCOEF15.clone());
    } else {
        return (DBASIS16.clone(), MATCOEF16.clone());
    };
}


fn compute_l(d1new: Vec<isize>, d2basis: Array1<isize>, d3basis: Array1<isize>, n1: isize, n2: isize, n3: isize, 
    d1: isize, d2: isize, d3: isize, k1sq: &Float, k2sq: &Float, k3sq: &Float, cache: &mut ltrian::TopologyCache) -> Float {

    let lend1 = d1new.len();
    let lend2 = d2basis.len();
    let lend3 = d3basis.len();

    // let gil = pyo3::Python::acquire_gil();
    // let mut matmul = PyArray3::zeros(gil.python(), [5,5,5], false);
    // let mut matmul: Array3<Complex64> = Array3::zeros((lend1,lend2,lend3));
    // let mut matmul: Vec<Vec<Vec<Complex>>> = vec![vec![vec![Complex::new(PREC);lend1];lend2];lend3];
    let mut res: Complex = Complex::new(PREC);
    let coef1 = dbabis_coef(d1).1;
    let coef2 = dbabis_coef(d2).1;
    let coef3 = dbabis_coef(d3).1;
    // coef1.assign(&coef1.clone().slice(s![..;2]));
    // let mut coef1new = coef1.clone().slice(s![..;2]);
    // println!("{}", lend2);
    // println!("{}", lend3);
    // println!("{:?}", coef1);
    // println!("{:?}", coef2);
    // println!("{:?}", coef3);

    for indx1 in 0..lend1 {
        for indx2 in 0..lend2 {
            for indx3 in 0..lend3 {
                let i: usize = d1new[indx1] as usize;
                let j: usize = d2basis[indx2] as usize;
                let l: usize = d3basis[indx3] as usize;
                // println!("{},{},{},{},{},{}",-n2 + FBABISEXP[j][0], FBABISEXP[j][1], -n1 + 
                //     FBABISEXP[i][0], FBABISEXP[i][1], -n3 + FBABISEXP[l][0], FBABISEXP[l][1]);

                res += coef1[2*indx1].clone()*coef2[indx2].clone()*coef3[indx3].clone()*ltrian::ltrian(-n2 + FBABISEXP[j][0], FBABISEXP[j][1], -n1 + 
                    FBABISEXP[i][0], FBABISEXP[i][1], -n3 + FBABISEXP[l][0], FBABISEXP[l][1], k1sq, k2sq, k3sq, 
                    &FBABISMASSES[j], &FBABISMASSES[i], &FBABISMASSES[l],FBABISMASSIND[j], FBABISMASSIND[i]+9, FBABISMASSIND[l]+18, cache);
                // println!("{}", res);
                // matmul[[indx1, indx2, indx3]] = Complex64::new(res.real().to_f64(),res.imag().to_f64());
                // println!("{}",-n2 as i32 + FBABISEXP[j][0] as i32);
                // println!("{}", coef1[2*indx1].clone()*coef2[indx2].clone()*coef3[indx3].clone()*ltrian::ltrian(-n2 as isize + FBABISEXP[j][0] as isize, 4, 5, 3, -n3 +FBABISEXP[i][0], FBABISEXP[l][1], k1sq, k2sq, k3sq, 
                //     &FBABISMASSES[j], &FBABISMASSES[i], &FBABISMASSES[l],1, 2, 3, cache));
            }
        }
    
    };
    // let c = einsum("ijk, i, j, k -> i", &[&matmul, &Array::from_vec(coef1.clone()).slice(s![..;2]), &Array::from_vec(coef2.clone()), &Array::from_vec(coef3.clone())]).unwrap();
    let fin: &Float = res.real();
    return fin/(4. * Float::with_val(PREC,PI * SQRT_PI));
    // return 2*np.real(np.einsum('ijk, i, j, k', matmul, MATCOEFs[d1][::2], MATCOEFs[d2], MATCOEFs[d3]))/(8 * PI * SQRT_PI)   
}

fn compute_l2(d2new: Vec<isize>, d3basis: Array1<isize>, n1: isize, n2: isize, n3: isize, 
    d2: isize, d3: isize, k1sq: &Float, k2sq: &Float, k3sq: &Float, cache: &mut ltrian::TopologyCache) -> Float {

    let lend2 = d2new.len();
    let lend3 = d3basis.len();
    let mut res: Complex = Complex::new(PREC);
    let coef2 = dbabis_coef(d2).1;
    let coef3 = dbabis_coef(d3).1;
    // coef2.assign(&coef2.clone().slice(s![..;2]));

    for indx2 in 0..lend2 {
        for indx3 in 0..lend3 {
            let i: usize = d2new[indx2] as usize;
            let j: usize = d3basis[indx3] as usize;
            // println!("{}, {}, {}, {}, {}, {}", -n2 + FBABISEXP[i][0], FBABISEXP[i][1], -n1, 0, 
            //     -n3 + FBABISEXP[j][0], FBABISEXP[j][1]);
            res += coef2[2*indx2].clone()*coef3[indx3].clone()*ltrian::ltrian(-n2 + FBABISEXP[i][0], FBABISEXP[i][1], -n1, 0, 
                    -n3 + FBABISEXP[j][0], FBABISEXP[j][1], k1sq, k2sq, k3sq, 
                    &FBABISMASSES[i], &Complex::new(PREC), &FBABISMASSES[j],FBABISMASSIND[i], 9usize, FBABISMASSIND[j]+18, cache);

        }
    };

    // let fin: &Float = res.real();
    return res.real().clone()/(4. * PI * SQRT_PI);
}

pub fn compute_j(n1: isize, n2: isize, n3: isize, 
    d1: isize, d2: isize, d3: isize, k1sq: &Float, k2sq: &Float, k3sq: &Float, cache: &mut ltrian::TopologyCache) -> Float {
    // # n1, n2, n3 are the exponents of q, k1pq, k2mq in the denominator
    let mut res: Complex = Complex::new(PREC);
    let mut d1basis = dbabis_coef(d1).0;
    let mut d2basis = dbabis_coef(d2).0;
    let (mut d3basis, coef3) = dbabis_coef(d3);
    let d3new = d3basis.slice_mut(s![..;2]).to_vec();
    let lend3 = d3new.len();
    // d3basis.assign(&d3basis.clone().slice(s![..;2]));

    if d1 == 0 {
        if d2 == 0 {
            if d3 == 0 {
                // let fin: &Float = ltrian::ltrian(-n2, 0, -n1, 0, -n3, 0, k1sq, k2sq, k3sq, 
                //     &Complex::new(PREC), &Complex::new(PREC), &Complex::new(PREC), 1, 2, 3, cache).real()/(8. *Float::with_val(PREC, PI * SQRT_PI));
                return ltrian::ltrian(-n2, 0, -n1, 0, -n3, 0, k1sq, k2sq, k3sq, 
                    &Complex::new(PREC), &Complex::new(PREC), &Complex::new(PREC), 0, 9, 18, cache).real().clone()/(8.  *PI * SQRT_PI);
            } else {


                for indx in 0..lend3 {
                    let i: usize = d3new[indx] as usize;
                    res += coef3[2*indx].clone()*ltrian::ltrian(-n2, 0, -n1, 0, -n3 + 
                    FBABISEXP[i][0], FBABISEXP[i][1], k1sq, k2sq, k3sq, 
                    &Complex::new(PREC), &Complex::new(PREC), &FBABISMASSES[i], 0, 9, FBABISMASSIND[i]+18, cache);
                };

                return res.real().clone()/(4. * PI * SQRT_PI);
            }             
        } else {
            let d2new = d2basis.slice_mut(s![..;2]).to_vec();
            return compute_l2(d2new, d3basis, n1, n2, n3, d2, d3, k1sq, k2sq, k3sq, cache);
        }
    } else {
        let d1new = d1basis.slice_mut(s![..;2]).to_vec();
        // println!("{:?}", d1new);
        return compute_l(d1new, d2basis, d3basis, n1, n2, n3, d1, d2, d3, k1sq, k2sq, k3sq, cache);
    };
}

