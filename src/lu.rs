pub struct LU {
    pub(crate) lu_size: i32,
    pub(crate) lu_nz: Vec<f64>,
    pub(crate) lu_rowind: Vec<i32>,
    pub(crate) l_colptr: Vec<i32>,
    pub(crate) u_colptr: Vec<i32>,

    pub(crate) row_perm: Vec<i32>,
    pub(crate) col_perm: Vec<i32>,
}

impl LU {
    pub(crate) fn new(nrow: i32, ncol: i32, lu_size: i32) -> Self {
        Self {
            lu_size,
            lu_nz: vec![0.0; lu_size as usize],
            lu_rowind: vec![0; lu_size as usize],
            l_colptr: vec![0; ncol as usize],
            u_colptr: vec![0; (ncol + 1) as usize],
            row_perm: vec![0; nrow as usize],
            col_perm: vec![0; ncol as usize],
        }
    }
}

#[derive(Clone)]
pub struct GP {
    /// 0=none, 1=partial, 2=threshold
    pub pivot_policy: i32,
    pub pivot_threshold: f64,
    pub drop_threshold: f64,
    pub col_fill_ratio: f64,
    pub fill_ratio: f64,
    pub expand_ratio: f64,
    pub col_perm: Option<Vec<i32>>,
    // col_perm_length: i32,
    pub col_perm_base: i32,
    // statistics_reporter_t reporter_func;
    // void*                 reporter_ctxt;
}

impl Default for GP {
    fn default() -> Self {
        Self {
            pivot_policy: 1,
            pivot_threshold: 1.0,
            drop_threshold: 0.0,
            col_fill_ratio: -1.0,
            fill_ratio: 4.0,
            expand_ratio: 1.2,
            col_perm: None,
            // col_perm_length = 0;
            col_perm_base: 0, // reporter_func   = NULL;
        }
    }
}

pub struct CSC {
    pub m: i32,
    pub n: i32,
    pub nnz: i32,
    pub base: i32,
    pub colptr: Vec<i32>,
    pub rowind: Vec<i32>,
}
