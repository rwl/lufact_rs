# lufact

Sparse LU factorization with partial pivoting as described in
"Sparse Partial Pivoting in Time Proportional to Arithmetic
Operations" by John R. Gilbert and Tim Peierls.

```bibtex
@article{Gilbert1988,
  doi = {10.1137/0909058},
  url = {https://doi.org/10.1137/0909058},
  year  = {1988},
  month = {sep},
  publisher = {Society for Industrial {\&} Applied Mathematics ({SIAM})},
  volume = {9},
  number = {5},
  pages = {862--874},
  author = {John R. Gilbert and Tim Peierls},
  title = {Sparse Partial Pivoting in Time Proportional to Arithmetic Operations},
  journal = {{SIAM} Journal on Scientific and Statistical Computing}
}
```

The original FORTRAN source is distributed in Sivan Toledo's work on
incomplete-factorization, from PARC in the early 1990s, and can be
found in the [ILU](http://www.netlib.org/linalg/ilu.tgz) package on Netlib.

## License

Licensed, with permission from John Gilbert and Tim Peierls, under either the

* Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or https://www.apache.org/licenses/LICENSE-2.0) or
* MIT license ([LICENSE-MIT](LICENSE-MIT) or https://opensource.org/licenses/MIT)

at your option.
