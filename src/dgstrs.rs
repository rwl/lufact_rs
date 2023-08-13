use crate::lu::*;
use crate::lufact::{lsolve, ltsolve, usolve, utsolve};

/// Solve `Ax=b` for one or more right-hand-sides given the numeric
/// factorization of A from [dgstrf()].
pub fn dgstrs(
    _gp: &GP,
    trans: char,
    n: i32,
    nrhs: i32,
    lu: &mut LU,
    ia: i32,
    ja: i32,
    b: &mut [f64],
    ib: i32,
    jb: i32,
    // desc_b: &CSC,
) -> i32 {
    // scalar_t*    rwork  = NULL;

    // int b_desc_type, b_ld;

    // statistics_reporter_t reporter_func;
    // void*                 reporter_ctxt;

    // let mut flops = 0.0;

    // (*info) = 0;

    /* extract data from gp object */

    // if ((gp) == NULL) {
    //   (*info) = -1;
    //   goto free_and_exit;
    // }

    // gp_get_statistics_reporter_c ((gp),&reporter_func, &reporter_ctxt);

    if ia != 1 {
        return -5;
    }
    if ja != 1 {
        return -6;
    }
    if ib != 1 {
        return -8;
    }
    if jb != 1 {
        return -9;
    }

    if nrhs != 1 {
        return -3;
    }

    /* we do not need this now since we assume a single dense rhs */
    /*
    b_desc_type = desc_b[0];
    if (b_desc_type != DESC_TYPE_DENSE) { (*info) = -10; goto free_and_exit; }
    b_ld = desc_b[1];
    */

    // if ((rwork  =
    //        (scalar_t*) malloc( n * sizeof(scalar_t) )) == NULL)
    //   { (*info) = -999; goto free_and_exit; }
    let mut rwork = vec![0.0; n as usize];

    if trans.to_ascii_uppercase() == 'N' {
        let _info = lsolve(
            n,
            &lu.lu_nz,
            &lu.lu_rowind,
            &lu.l_colptr,
            &lu.u_colptr,
            &lu.row_perm,
            &lu.col_perm,
            b,
            &mut rwork,
        );

        let _info = usolve(
            n,
            &lu.lu_nz,
            &lu.lu_rowind,
            &lu.l_colptr,
            &lu.u_colptr,
            &lu.row_perm,
            &lu.col_perm,
            &rwork,
            b,
        );
    } else if trans.to_ascii_uppercase() == 'T' {
        let _info = utsolve(
            n,
            &lu.lu_nz,
            &lu.lu_rowind,
            &lu.l_colptr,
            &lu.u_colptr,
            &lu.row_perm,
            &lu.col_perm,
            b,
            &mut rwork,
        );

        let _info = ltsolve(
            n,
            &lu.lu_nz,
            &lu.lu_rowind,
            &lu.l_colptr,
            &lu.u_colptr,
            &lu.row_perm,
            &lu.col_perm,
            &rwork,
            b,
        );
    } else {
        return -1;
    }

    let flops = (2 * (lu.u_colptr[n as usize] - 1)) as f64;

    // free_and_exit:
    //   if (rwork)  free(rwork);

    // if (reporter_func)
    //   (*reporter_func)(reporter_ctxt,"FLOPS",&flops);
    log::info!("flops: {}", flops);

    0
}
