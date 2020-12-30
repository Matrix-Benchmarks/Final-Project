
# Low Rank Matrix Completion Benchmarks

This repository contains a collection of low rank matrix completion algorithms, a framework to test those algorithms on many different datasets, and a website for displaying those results, partially inspired by what [ann-benchmarks](https://github.com/erikbern/ann-benchmarks) does for near neighbour search. This repository was originally forked from the [the MatrixIRLS repo](https://github.com/ckuemmerle/MatrixIRLS), a repository the authors of [Escaping Saddle Points in Ill-Conditioned Matrix Completion with a Scalable Second Order Method](https://arxiv.org/pdf/2009.02905.pdf) wrote to test their algorithm against other low rank matrix completion algorithms. Testing of algorithms can be seen at the following [webpage](https://matrix-benchmarks.github.io/Low-Rank-Completion).

## Run instructions
* Open MATLAB and run `setup` or, alternatively, add subfolders manually to path. Then run some or all of the tests in the compare folder, and then run the python file compare.py in the output folder to generate the graphs. This will also allow local generation of the website using the generated plots.

The algorithms tested can be found in subfolders of `algorithms`. Some of the algorithms use methods for sparse access of factorized matrices implemented in C compiled as `.mex` files, and other auxiliary files contained in the `tools` folder.

## List of algorithms
We gratefully acknowledge the authors of the following matrix completion algorithms. For re-use of the algorithms, please refer to the provided links and contact the authors if the respective terms of usage are unclear.

* 'MatrixILRS' (_Matrix Iteratively Reweighted Least Squares_) by Christian Kümmerle and Claudio M. Verdun [[KV20]](https://arxiv.org/pdf/2009.02905.pdf),
[available](https://github.com/ckuemmerle/MatrixIRLS) as the parent github repository.
* `ASD` (_Alternating Steepest Descent_) and `ScaledASD` (_Scaled Alternating Steepest Descent_) by Jared Tanner and Ke Wei [[TW16]](https://doi.org/10.1016/j.acha.2015.08.003), [available](http://www.sdspeople.fudan.edu.cn/weike/code/mc20140528.tar) at [Ke Wei's website](http://www.sdspeople.fudan.edu.cn/weike/publications.html).
* `NIHT` (_Normalized Iterative Hard Thresholding_ [[TW12]](https://doi.org/10.1137/120876459)) and `CGIHT` (_Conjugate Gradient Iterative Hard Thresholding_ [[BTW15]](https://doi.org/10.1093/imaiai/iav011)) by Jeffrey Blanchard, Jared Tanner and Ke Wei, [available](http://www.sdspeople.fudan.edu.cn/weike/code/mc20140528.tar) at [Ke Wei's website](http://www.sdspeople.fudan.edu.cn/weike/publications.html).
* `LMaFit` (Low-rank Matrix Fitting algorithm [[WYZ12]](https://doi.org/10.1007/s12532-012-0044-1)) by Zaiwen Wen, Wotao Yin, and Yin Zhang, available [here](http://lmafit.blogs.rice.edu).
* `LRGeomCG` (Low-rank matrix completion by Geometric CG [[V13]](https://doi.org/10.1137/110845768)) by Bart Vandereycken, [available]((http://www.unige.ch/math/vandereycken/matrix_completion.html)) at [Bart Vandereycken's website](http://www.unige.ch/math/vandereycken/research.php).
* `MatrixIRLS` (Matrix Iteratively Reweighted Least Squares). This paper/repository.
* `HM-IRLS` (Harmonic Mean Iteratively Reweighted Least Squares [[KS18]](http://www.jmlr.org/beta/papers/v19/17-244.html)) and `AM-IRLS` by Christian Kümmerle and Juliane Sigl, available [here](https://github.com/ckuemmerle/hm_irls).
* `IRLS-p` and `sIRLS-p`: IRLS with weight operator acting on row space only, solving linear systems by gradient descent [[MF12]](http://www.jmlr.org/beta/papers/v13/mohan12a.html), by Karthik Mohan and Maryam Fazel, [available](https://faculty.washington.edu/mfazel/IRLS_final.zip) at [Maryam Fazel's website](https://faculty.washington.edu/mfazel/).
* `IRucLq` and `tIRucLq` ((truncated) Iterative Reweighted unconstrained Lq for low-rank matrix recovery [[LXY13]](https://epubs.siam.org/doi/abs/10.1137/110840364)) by Zaiwen Wen, Wotao Yin and Yin Zhang, [available](https://xu-yangyang.github.io/codes/IRucLq.zip) at [Yangyang Xu's website](https://xu-yangyang.github.io/papers.html). 
* `IRLS-col` and `IRLS-row`: IRLS with weight operators that act on the column or row space, respectively, and thus very similar to algorithms of [[FRW11]](https://epubs.siam.org/doi/abs/10.1137/100811404) and [[MF12]](http://www.jmlr.org/beta/papers/v13/mohan12a.html). Main purpose: illustrate the influence of the choice of weight operator. 
* `R2RILS` (Rank 2r Iterative Least Squares [[BN20]](https://arxiv.org/abs/2002.01849)) by Jonathan Bauch and Boaz Nadler, available [here](https://github.com/Jonathan-WIS/R2RILS), see also [here](http://www.wisdom.weizmann.ac.il/~nadler/Projects/R2RILS/R2RILS.html). 
* `R3MC` (Riemannian three-factor algorithm for low-rank matrix completion [[MS14]](https://doi.org/10.1109/CDC.2014.7039534)) by Bamdev Mishra and Rodolphe Sepulchre, [available](https://dl.dropboxusercontent.com/s/qzxgax0bg3s8oe2/R3MC_17feb_2017.zip) at [Bamdev Mishra's website](https://bamdevmishra.in/codes/r3mc/). We also included `R3MC-rankupd`, a variant of `R3MC` which optimizes on fixed-rank manifolds with increasing rank (see also [[MS14]](https://doi.org/10.1109/CDC.2014.7039534)).
* `RTRMC` (Riemannian trust-region method for low-rank matrix completion [[BA15]](https://doi.org/10.1016/j.laa.2015.02.027)) by Nicholas Boumal and Pierre-Antoine Absil, available [here](http://web.math.princeton.edu/~nboumal/RTRMC/index.html). The version provided in this repository uses [Manopt 6.0](https://www.manopt.org).
* `ScaledGD` (Scaled Gradient Descent [[TMC20]](https://arxiv.org/abs/2005.08898)) by Tian Tong, Cong Ma and Yuejie Chi. Available [here](https://github.com/Titan-Tong/ScaledGD).

## About this repository
##### Rice Comp 414 Final Project: 
* Joshua Engels
* Richard Morse
##### Authors of original MatrixICLR repository:
* Christian Kümmerle 
* Claudio M. Verdun

## References
 -  [[KV20]](https://arxiv.org/pdf/2009.02905.pdf) Christian Kümmerle and Claudio M. Verdun, [**Escaping Saddle Points in Ill-Conditioned Matrix Completion with a Scalable Second Order Method**](https://arxiv.org/pdf/2009.02905.pdf)'MatrixILRS' (_Matrix Iteratively Reweighted Least Squares_) by Christian Kümmerle and Claudio M. Verdun [[KV20]](https://arxiv.org/pdf/2009.02905.pdf), _arXiv preprint_, arXiv:2009.02905, 2020.
[available](https://github.com/ckuemmerle/MatrixIRLS) as the parent github repository
 - [[KS18]](http://www.jmlr.org/beta/papers/v19/17-244.html) Christian Kümmerle and Juliane Sigl, [**Harmonic Mean Iteratively Reweighted Least Squares for Low-Rank Matrix Recovery**](http://www.jmlr.org/beta/papers/v19/17-244.html). _J. Mach. Learn. Res._, 19(47):1–49, 2018.
- [[FRW11]](https://epubs.siam.org/doi/abs/10.1137/100811404) Massimo Fornasier, Holger Rauhut and Rachel Ward, [**Low-rank matrix recovery via iteratively reweighted least squares minimization**](https://epubs.siam.org/doi/abs/10.1137/100811404). _SIAM J. Optim._, 21(4):1614–1640, 2011.
- [[MF12]](http://www.jmlr.org/beta/papers/v13/mohan12a.html) Karthik Mohan and Maryam Fazel, [**Iterative reweighted algorithms for matrix rank minimization**](http://www.jmlr.org/beta/papers/v13/mohan12a.html), _J. Mach. Learn. Res._, 13 (1):3441–3473, 2012.
- [[LXY13]](https://epubs.siam.org/doi/abs/10.1137/110840364) Ming-Jun Lai, Yangyang Xu and Wotao Yin, [**Improved iteratively reweighted least squares for unconstrained smoothed $\ell_q$ minimization**](https://epubs.siam.org/doi/abs/10.1137/110840364), _SIAM J. Numer. Anal._, 51(2):927-957, 2013.
- [[TW16]](https://doi.org/10.1016/j.acha.2015.08.003) Jared Tanner and Ke Wei, [**Low rank matrix completion by alternating steepest descent methods**](https://doi.org/10.1016/j.acha.2015.08.003). _Appl. Comput. Harmon. Anal._, 40(2):417–429, 2016.
- [[TW12]](https://doi.org/10.1137/120876459) Jared Tanner and Ke Wei, [**Normalized Iterative Hard Thresholding for Matrix Completion**](https://doi.org/10.1137/120876459)), _SIAM J. Sci. Comput._, 35(5):S104–S125, 2012.
- [[BTW15]](https://doi.org/10.1093/imaiai/iav011) Jeffrey D. Blanchard, Jared Tanner and Ke Wei, [**CGIHT: conjugate gradient iterative hard thresholding for compressed sensing and matrix completion**](https://doi.org/10.1093/imaiai/iav011), _Inf. Inference_, 4(4):289-327, 2015.
- [[WYZ12]](https://doi.org/10.1007/s12532-012-0044-1) Zaiwen Wen, Wotao Yin and Yin Zhang, [**Solving a low-rank factorization model for matrix completion by a nonlinear successive over-relaxation algorithm**](https://doi.org/10.1007/s12532-012-0044-1), _Math. Prog. Comp._ 4(4):333–361, 2012.
-  [[V13]](https://doi.org/10.1137/110845768) Bart Vandereycken, [**Low-Rank Matrix Completion by Riemannian Optimization**](https://doi.org/10.1137/110845768), _SIAM J. Optim._, 23(2):1214-1236, 2013.
- [[BN20]](https://arxiv.org/abs/2002.01849), Jonathan Bauch and Boaz Nadler, [**Rank 2r iterative least squares: efficient recovery of ill-conditioned low rank matrices from few entries**](https://arxiv.org/abs/2002.01849), _arXiv preprint_, arXiv:2002.01849, 2020.
- [[MS14]](https://doi.org/10.1109/CDC.2014.7039534) Bamdev Mishra and Rodolphe Sepulchre, [**R3MC: A Riemannian three-factor algorithm for low-rank matrix completion**](https://doi.org/10.1109/CDC.2014.7039534), In _53rd IEEE Conference on Decision and Control_, 1137-1142. IEEE, 2014.
- [[BA15]](https://doi.org/10.1016/j.laa.2015.02.027) Nicholas Boumal and Pierre-Antoine Absil, [**Low-rank matrix completion via preconditioned optimization on the Grassmann manifold**](https://doi.org/10.1016/j.laa.2015.02.027), _Linear Algebra Appl._, 15(475):200–239, 2015.
- [[TMC20]](https://arxiv.org/abs/2005.08898) Tian Tong, Cong Ma and Yuejie Chi, [**Accelerating Ill-Conditioned Low-Rank Matrix Estimation via Scaled Gradient Descent**](https://arxiv.org/abs/2005.08898), _arXiv preprint_, arXiv:2005.08898, 2020.
