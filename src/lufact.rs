pub(crate) fn maxmatch(
    nrows: i32,
    ncols: i32,
    colstr: &[i32],
    rowind: &[i32],
    prevcl: &mut [i32],
    prevrw: &mut [i32],
    marker: &mut [i32],
    tryrow: &mut [i32],
    nxtchp: &mut [i32],
    rowset: &mut [i32],
    colset: &mut [i32],
) {
    unsafe {
        lufact_sys::maxmatch_(
            &nrows,
            &ncols,
            colstr.as_ptr(),
            rowind.as_ptr(),
            prevcl.as_mut_ptr(),
            prevrw.as_mut_ptr(),
            marker.as_mut_ptr(),
            tryrow.as_mut_ptr(),
            nxtchp.as_mut_ptr(),
            rowset.as_mut_ptr(),
            colset.as_mut_ptr(),
        );
    }
}

pub(crate) fn ludfs(
    jcol: i32,
    a: &[f64],
    arow: &[i32],
    acolst: &[i32],
    lastlu: &mut i32,
    lurow: &mut [i32],
    lcolst: &mut [i32],
    ucolst: &mut [i32],
    rperm: &mut [i32],
    cperm: &mut [i32],
    dense: &mut [f64],
    found: &mut [i32],
    parent: &mut [i32],
    child: &mut [i32],
) -> i32 {
    let mut error: i32 = 0;
    unsafe {
        lufact_sys::ludfs_(
            &jcol,
            a.as_ptr(),
            arow.as_ptr(),
            acolst.as_ptr(),
            lastlu,
            lurow.as_mut_ptr(),
            lcolst.as_mut_ptr(),
            ucolst.as_mut_ptr(),
            rperm.as_mut_ptr(),
            cperm.as_mut_ptr(),
            dense.as_mut_ptr(),
            found.as_mut_ptr(),
            parent.as_mut_ptr(),
            child.as_mut_ptr(),
            &mut error,
        )
    }
    error
}

pub(crate) fn lucomp(
    jcol: i32,
    lastlu: &mut i32,
    lu: &mut [f64],
    lurow: &mut [i32],
    lcolst: &mut [i32],
    ucolst: &mut [i32],
    rperm: &[i32],
    cperm: &[i32],
    dense: &mut [f64],
    found: &mut [i32],
    pattern: &mut [i32],
) -> f64 {
    let mut flops: f64 = 0.0;
    unsafe {
        lufact_sys::lucomp_(
            &jcol,
            lastlu,
            lu.as_mut_ptr(),
            lurow.as_mut_ptr(),
            lcolst.as_mut_ptr(),
            ucolst.as_mut_ptr(),
            rperm.as_ptr(),
            cperm.as_ptr(),
            dense.as_mut_ptr(),
            found.as_mut_ptr(),
            pattern.as_mut_ptr(),
            &mut flops,
        )
    }
    flops
}

pub(crate) fn lucopy(
    pivot: i32,
    pthresh: f64,
    mut dthresh: f64,
    mut nzcount: i32,
    jcol: i32,
    ncol: i32,
    lastlu: &mut i32,
    lu: &mut [f64],
    lurow: &mut [i32],
    lcolst: &mut [i32],
    ucolst: &mut [i32],
    rperm: &mut [i32],
    cperm: &mut [i32],
    dense: &mut [f64],
    pattern: &mut [i32],
    twork: &mut [f64],
    flops: &mut f64,
    zpivot: &mut i32,
) {
    unsafe {
        lufact_sys::lucopy_(
            &pivot,
            &pthresh,
            &mut dthresh,
            &mut nzcount,
            &jcol,
            &ncol,
            lastlu,
            lu.as_mut_ptr(),
            lurow.as_mut_ptr(),
            lcolst.as_mut_ptr(),
            ucolst.as_mut_ptr(),
            rperm.as_mut_ptr(),
            cperm.as_mut_ptr(),
            dense.as_mut_ptr(),
            pattern.as_mut_ptr(),
            twork.as_mut_ptr(),
            flops,
            zpivot,
        );
    }
}

pub(crate) fn lsolve(
    n: i32,
    lu: &[f64],
    lurow: &[i32],
    lcolst: &[i32],
    ucolst: &[i32],
    rperm: &[i32],
    cperm: &[i32],
    b: &[f64],
    x: &mut [f64],
) -> i32 {
    let mut error: i32 = 0;
    unsafe {
        lufact_sys::lsolve_(
            &n,
            lu.as_ptr(),
            lurow.as_ptr(),
            lcolst.as_ptr(),
            ucolst.as_ptr(),
            rperm.as_ptr(),
            cperm.as_ptr(),
            b.as_ptr(),
            x.as_mut_ptr(),
            &mut error,
        );
    }
    error
}

pub(crate) fn usolve(
    n: i32,
    lu: &[f64],
    lurow: &[i32],
    lcolst: &[i32],
    ucolst: &[i32],
    rperm: &[i32],
    cperm: &[i32],
    b: &[f64],
    x: &mut [f64],
) -> i32 {
    let mut error: i32 = 0;
    unsafe {
        lufact_sys::usolve_(
            &n,
            lu.as_ptr(),
            lurow.as_ptr(),
            lcolst.as_ptr(),
            ucolst.as_ptr(),
            rperm.as_ptr(),
            cperm.as_ptr(),
            b.as_ptr(),
            x.as_mut_ptr(),
            &mut error,
        );
    }
    error
}

pub(crate) fn ltsolve(
    n: i32,
    lu: &[f64],
    lurow: &[i32],
    lcolst: &[i32],
    ucolst: &[i32],
    rperm: &[i32],
    cperm: &[i32],
    b: &[f64],
    x: &mut [f64],
) -> i32 {
    let mut error: i32 = 0;
    unsafe {
        lufact_sys::ltsolve_(
            &n,
            lu.as_ptr(),
            lurow.as_ptr(),
            lcolst.as_ptr(),
            ucolst.as_ptr(),
            rperm.as_ptr(),
            cperm.as_ptr(),
            b.as_ptr(),
            x.as_mut_ptr(),
            &mut error,
        );
    }
    error
}

pub(crate) fn utsolve(
    n: i32,
    lu: &[f64],
    lurow: &[i32],
    lcolst: &[i32],
    ucolst: &[i32],
    rperm: &[i32],
    cperm: &[i32],
    b: &[f64],
    x: &mut [f64],
) -> i32 {
    let mut error: i32 = 0;
    unsafe {
        lufact_sys::utsolve_(
            &n,
            lu.as_ptr(),
            lurow.as_ptr(),
            lcolst.as_ptr(),
            ucolst.as_ptr(),
            rperm.as_ptr(),
            cperm.as_ptr(),
            b.as_ptr(),
            x.as_mut_ptr(),
            &mut error,
        );
    }
    error
}
