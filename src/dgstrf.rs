use crate::lu::*;
use crate::lufact::{lucomp, lucopy, ludfs, maxmatch};

pub fn dgstrf(gp: &GP, mrows: i32, ncols: i32, a_nz: &[f64], desc_a: &mut CSC) -> Result<LU, i32> {
    // let a_rowind: Vec<i32>;
    // let a_colptr: Vec<i32>;

    /* work arrays */

    // let rwork: Vec<f64>;
    // let twork: Vec<f64>;

    // let found: Vec<i32>;
    // let parent: Vec<i32>;
    // let child: Vec<i32>;

    // let pattern: Vec<i32>;

    // let cmatch: Vec<i32>;
    // let rmatch: Vec<i32>;

    /* copies of object parameters */

    // let pivot_policy: i32;
    // let  pivot_threshold: f64;
    // let  drop_threshold: f64;
    // let  col_fill_ratio: f64;
    // let  fill_ratio: f64;
    // let  expand_ratio: f64;

    /* local variables */

    let nrow: i32 = mrows;
    let ncol: i32 = ncols;

    // let a_desc_type: i32;
    // let a_m: i32;
    // let a_n: i32;
    // let a_nnz: i32;
    // let a_base: i32;

    // let jcol: i32;
    // let i: i32;
    // let lasta: i32;
    // let lastlu: i32;
    let mut zpivot: i32 = 0;

    // let local_pivot_policy: i32;
    let mut nz_count_limit: i32;

    // let user_col_perm: Vec<i32>;
    // let user_col_perm_length: i32;
    // let user_col_perm_base: i32;

    // statistics_reporter_t reporter_func;
    // void*                 reporter_ctxt;

    let mut flops: f64 = 0.0;

    // let ujj: f64;
    // let minujj: f64;

    // int       out_of_mem = FALSE;
    // int       eline = -1;

    let mut pivt_row: i32;
    let mut orig_row: i32;
    let mut this_col: i32;
    let mut othr_col: i32;

    /* constants */

    // const izero: i32 = 0;
    // const zero: f64 = 0.0; /* this is not good for complex !!! replace with macro */
    /* extract data from gp object */

    // if ((gp) == NULL) {
    //   (*info) = -1;
    //   goto free_and_exit;
    // }
    // gp_get_pivot_policy_c        ((gp),&pivot_policy);
    // gp_get_pivot_threshold_c     ((gp),&pivot_threshold);
    // gp_get_drop_threshold_c      ((gp),&drop_threshold);
    // gp_get_col_fill_ratio_c      ((gp),&col_fill_ratio);
    // gp_get_fill_ratio_c          ((gp),&fill_ratio);
    // gp_get_expand_ratio_c        ((gp),&expand_ratio);
    // gp_get_statistics_reporter_c ((gp),&reporter_func,
    //                                                   &reporter_ctxt);
    // gp_get_col_perm_c            ((gp),&user_col_perm,
    //                                                   &user_col_perm_length,
    //                                                   &user_col_perm_base);

    let pivot_policy: i32 = gp.pivot_policy;
    let pivot_threshold: f64 = gp.pivot_threshold;
    let drop_threshold: f64 = gp.drop_threshold;
    let col_fill_ratio: f64 = gp.col_fill_ratio;
    let fill_ratio: f64 = gp.fill_ratio;
    let expand_ratio: f64 = gp.expand_ratio;

    let user_col_perm: &Option<Vec<i32>> = &gp.col_perm;
    // let user_col_perm_length: i32 = gp.;
    let user_col_perm_base: i32 = gp.col_perm_base;

    // println!(
    //     "piv pol={} piv_thr={} drop_thr={} col_fill_rt={}\n",
    //     pivot_policy, pivot_threshold, drop_threshold, col_fill_ratio
    // );

    // pivot_threshold = 0.001;
    // if pivot_threshold == 0.0 {
    //     pivot_policy = 0; // no pivoting
    // }
    // pivot_policy = 0; // no pivoting

    // if a column permutation is specified,
    // it must be a length ncol permutation.
    if let Some(col_perm) = user_col_perm {
        if col_perm.len() != ncol as usize {
            return Err(-1);
        }
    }

    // extract data from a's array descriptor

    // a_desc_type = desc_a[CSC_DESC_TYPE];
    // if (a_desc_type != DESC_TYPE_CSC) {
    //   (*info) = -5;
    //   goto free_and_exit;
    // }
    let _a_m = desc_a.m;
    let a_n = desc_a.n;
    let a_nnz = desc_a.nnz;
    let a_base = desc_a.base;
    let a_colptr = &mut desc_a.colptr;
    let a_rowind = &mut desc_a.rowind;

    // convert the descriptor to 1-base if necessary
    if a_base == 0 {
        for jcol in 0..a_n + 1 {
            a_colptr[jcol as usize] += 1;
        }
        for jcol in 0..a_nnz {
            a_rowind[jcol as usize] += 1;
        }
        desc_a.base = 1;
        // a_base = 1;
    }

    // Allocate work arrays.
    let mut rwork = vec![0.0; nrow as usize];
    let mut twork = vec![0.0; nrow as usize];

    let mut found = vec![0; nrow as usize];
    let mut child = vec![0; nrow as usize];
    let mut parent = vec![0; nrow as usize];

    let mut pattern = vec![0; nrow as usize];

    // Create lu structure
    let mut lu = LU::new(nrow, ncol, (a_nnz as f64 * fill_ratio) as i32);

    // Compute max matching. We use elements of the lu structure
    // for all the temporary arrays needed.
    let mut cmatch = vec![0; ncol as usize];
    let mut rmatch = vec![0; nrow as usize];
    {
        lu.l_colptr.fill(0);
        lu.u_colptr.fill(0);
        lu.col_perm.fill(0);
        lu.lu_rowind.fill(0);
        rmatch.fill(0);
        cmatch.fill(0);
        lu.row_perm.fill(0);
    }

    maxmatch(
        nrow,              // in.
        ncol,              // in.
        a_colptr,          // in.
        a_rowind,          // in.
        &mut lu.l_colptr,  // work. prevcl(cols)
        &mut lu.u_colptr,  // work. prevrw(cols)
        &mut lu.row_perm,  // work. marker(rows)
        &mut lu.col_perm,  // work. tryrow(cols)
        &mut lu.lu_rowind, // work. nxtchp(cols)
        &mut rmatch,       // out.  rowset(rows)
        &mut cmatch,       // out.  colset(cols)
    );

    for jcol in 0..ncol {
        if cmatch[jcol as usize] == 0 {
            println!("Warning: Perfect matching not found");
            break;
        }
    }

    /*
      for (jcol = 0; jcol < ncol; jcol++)
        cmatch[jcol] = rmatch[jcol] = jcol+1;
    */
    /* Initialize useful values and zero out the dense vectors.
    If we are threshold pivoting, get row counts. */

    let mut lastlu = 0;

    let mut local_pivot_policy = pivot_policy;
    // (*info) = 0;
    let _lasta = a_colptr[ncol as usize] - 1;
    lu.u_colptr[0] = 1;

    // ifill_ (pattern, &nrow, &izero);
    // ifill_ (found, &nrow, &izero);
    // rfill_ (rwork, &nrow, &zero);
    // ifill_ ((*lu)->row_perm, &nrow, &izero);
    pattern.fill(0);
    found.fill(0);
    rwork.fill(0.0);
    lu.row_perm.fill(0);

    if let Some(user_col_perm) = user_col_perm {
        println!("user_col_perm_base = {}", user_col_perm_base);
        for jcol in 0..ncol {
            lu.col_perm[jcol as usize] = user_col_perm[jcol as usize] + (1 - user_col_perm_base);
        }
    } else {
        for jcol in 0..ncol {
            lu.col_perm[jcol as usize] = jcol + 1;
        }
    }

    // compute one column at a time
    for jcol in 1..=ncol {
        // mark pointer to new column, ensure it is large enough

        if lastlu + nrow >= lu.lu_size {
            let new_size: i32 = (lu.lu_size as f64 * expand_ratio) as i32;

            // eprintln!("expanding to %d nonzeros...",new_size);

            // if ((lu.lu_nz =
            //       (scalar_t*) realloc( lu.lu_nz,
            //                            (new_size * sizeof(scalar_t)) )) == NULL)
            //   { out_of_mem = TRUE; eline = __LINE__; goto free_and_exit; }
            //
            // if ((lu.lu_rowind =
            //      (int*) realloc( lu.lu_rowind,
            //                            (new_size * sizeof(int)) )) == NULL)
            //   { out_of_mem = TRUE; eline = __LINE__; goto free_and_exit; }

            let mut lu_nz = vec![0.0; new_size as usize];
            lu_nz[..lu.lu_size as usize].copy_from_slice(&lu.lu_nz[..]);
            lu.lu_nz = lu_nz;

            let mut lu_rowind = vec![0; new_size as usize];
            lu_rowind[..lu.lu_size as usize].copy_from_slice(&lu.lu_rowind[..]);
            lu.lu_rowind = lu_rowind;

            lu.lu_size = new_size;
        }

        /* Set up nonzero pattern */

        {
            // int jjj;

            let jjj = lu.col_perm[(jcol - 1) as usize];
            for i in a_colptr[(jjj - 1) as usize]..a_colptr[jjj as usize] {
                pattern[(a_rowind[(i - 1) as usize] - 1) as usize] = 1;
            }

            this_col = lu.col_perm[(jcol - 1) as usize];
            orig_row = cmatch[(this_col - 1) as usize];

            pattern[(orig_row - 1) as usize] = 2;

            if lu.row_perm[(orig_row - 1) as usize] != 0 {
                println!("ERROR: PIVOT ROW FROM MAX-MATCHING ALREADY USED.");
                // exit(1);
                return Err(1);
            }

            // pattern[ this_col - 1 ] = 2;
        };

        // Depth-first search from each above-diagonal nonzero of column
        // jcol of A, allocating storage for column jcol of U in
        // topological order and also for the non-fill part of column
        // jcol of L.

        let info = ludfs(
            jcol,
            a_nz,
            a_rowind,
            a_colptr,
            &mut lastlu,
            &mut lu.lu_rowind,
            &mut lu.l_colptr,
            &mut lu.u_colptr,
            &mut lu.row_perm,
            &mut lu.col_perm,
            &mut rwork,
            &mut found,
            &mut parent,
            &mut child,
        );

        if info != 0 {
            // info = -100;
            // goto free_and_exit;
            // return info;
            return Err(-100);
        }

        // Compute the values of column jcol of L and U in the dense
        // vector, allocating storage for fill in L as necessary.

        flops = lucomp(
            jcol,
            &mut lastlu,
            &mut lu.lu_nz,
            &mut lu.lu_rowind,
            &mut lu.l_colptr,
            &mut lu.u_colptr,
            &lu.row_perm,
            &lu.col_perm,
            &mut rwork,
            &mut found,
            &mut pattern,
        );

        if rwork[(orig_row - 1) as usize] == 0.0 {
            println!("WARNING: MATCHING TO A ZERO");

            for i in a_colptr[(jcol - 1) as usize]..a_colptr[jcol as usize] {
                print!(
                    "({},{}) ",
                    a_rowind[(i - 1) as usize],
                    a_nz[(i - 1) as usize]
                );
            }
            println!(". orig_row={}", orig_row);
        }

        // Copy the dense vector into the sparse data structure, find the
        // diagonal element (pivoting if specified), and divide the
        // column of L by it.

        nz_count_limit = (col_fill_ratio
            * ((a_colptr[this_col as usize] - a_colptr[(this_col - 1) as usize] + 1) as f64))
            as i32;

        lucopy(
            local_pivot_policy,
            pivot_threshold,
            drop_threshold,
            nz_count_limit,
            jcol,
            ncol,
            &mut lastlu,
            &mut lu.lu_nz,
            &mut lu.lu_rowind,
            &mut lu.l_colptr,
            &mut lu.u_colptr,
            &mut lu.row_perm,
            &mut lu.col_perm,
            &mut rwork,
            &mut pattern,
            &mut twork,
            &mut flops,
            &mut zpivot,
        );

        if zpivot == -1 {
            // info = jcol;
            // goto free_and_exit;
            // return info;
            return Err(jcol);
        }

        {
            // int jjj;

            let jjj = lu.col_perm[(jcol - 1) as usize];
            for i in a_colptr[(jjj - 1) as usize]..a_colptr[jjj as usize] {
                pattern[(a_rowind[(i - 1) as usize] - 1) as usize] = 0;
            }

            pattern[(orig_row - 1) as usize] = 0;

            pivt_row = zpivot;
            othr_col = rmatch[(pivt_row - 1) as usize];

            cmatch[(this_col - 1) as usize] = pivt_row;
            cmatch[(othr_col - 1) as usize] = orig_row;
            rmatch[(orig_row - 1) as usize] = othr_col;
            rmatch[(pivt_row - 1) as usize] = this_col;

            // pattern[ this_col - 1 ] = 0;
        }

        // If there are no diagonal elements after this column, change
        // the pivot mode.

        if jcol == nrow {
            local_pivot_policy = -1;
        }
    } // end of jcol loop

    // Fill in the zero entries of the permutation vector, and renumber the
    // rows so the data structure represents L and U, not PtL and PtU.

    let mut jcol = ncol + 1;
    for i in 0..nrow {
        if lu.row_perm[i as usize] == 0 {
            lu.row_perm[i as usize] = jcol;
            jcol = jcol + 1;
        }
    }

    for i in 0..lastlu {
        lu.lu_rowind[i as usize] = lu.row_perm[(lu.lu_rowind[i as usize] - 1) as usize];
    }
    /* Return */

    // free_and_exit:

    // println!("rperm: {}", lu.row_perm);
    // println!("cperm: {}", lu.col_perm);

    // if (out_of_mem) {
    //   fprintf(stderr,
    //           "Out of space in gstrf_gp. Limit of maxlu=%d exceeded at column %d line %d\n",
    //           (*lu)->lu_size,jcol,eline);
    //   (*info) = -999;
    // }

    // if (rmatch) free(rmatch);
    // if (cmatch) free(cmatch);

    // if (pattern) free(pattern);

    // if (parent) free(parent);
    // if (child)  free(child);
    // if (found)  free(found);
    // if (twork)  free(rwork);
    // if (rwork)  free(rwork);

    // if ((*info) != 0) {
    //   if (*lu) {
    //
    //     if ((*lu)->row_perm)  free((*lu)->col_perm);
    //     if ((*lu)->row_perm)  free((*lu)->row_perm);
    //     if ((*lu)->u_colptr)  free((*lu)->u_colptr);
    //     if ((*lu)->l_colptr)  free((*lu)->l_colptr);
    //     if ((*lu)->lu_rowind) free((*lu)->lu_rowind);
    //     if ((*lu)->lu_nz)     free((*lu)->lu_nz);
    //
    //     free (*lu);
    //     *lu = NULL;
    //   }
    // } else {
    let mut minujj = f64::INFINITY; // 1.0 / 0.0;

    for jcol in 1..=ncol {
        let ujj = f64::abs(lu.lu_nz[(lu.l_colptr[(jcol - 1) as usize] - 2) as usize]);
        if ujj < minujj {
            minujj = ujj;
        }
    }

    // println!(">>> last = {}, min = {}",ujj,minujj);
    // }

    // if (reporter_func) {
    //   (*reporter_func)(reporter_ctxt,"FLOPS",&flops);
    //   flops = (double) lastlu;
    //   (*reporter_func)(reporter_ctxt,"NONZEROS",&flops);
    // }
    println!("FLOPS: {}", flops);
    println!("NONZEROS: {}", lastlu);

    Ok(lu)
}
