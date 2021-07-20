const GMP_VERSION: c_int = 6;
const GMP_VERSION_MINOR: c_int = 2;
const GMP_VERSION_PATCHLEVEL: c_int = 1;
const GMP_LIMB_BITS: c_int = 64;
const GMP_NAIL_BITS: c_int = 0;
type GMP_LIMB_T = c_ulong;
const GMP_CC: *const c_char = b"/home/diogo/miniconda3/bin/x86_64-conda-linux-gnu-cc\0".as_ptr() as _;
const GMP_CFLAGS: *const c_char = b"-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /home/diogo/miniconda3/include\0".as_ptr() as _;
