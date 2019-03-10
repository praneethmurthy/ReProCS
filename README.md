This folder contains the code accompanying pre-print.

[1] "Provable Dynamic Robust PCA or Robust Subspace Tracking", Praneeth Narayanamurthy and Namrata Vaswani, IEEE Trans. Info. Theory, 2019 (available at arXiv:1705.08948).

If you use this code please also cite the following papers

[2] "An online algorithm for separating sparse  and low-dimensional signal sequences from their sum", Han Guo, Chenlu Qiu, and Namrata Vaswani, IEEE Trans. Sig. Proc., 2014.

[3] "Recursive Robust PCA or Recursive Sparse Recovery in Large but Structure Noise", Chenlu Qiu, Namrata Vaswani, Brain Lois, and Leslie Hogben, IEEE Trans. Info. Theory., 2014.

[4] "Real-time Robust Principal Components' Pursuit", Chenlu Qiu, and Namrata Vaswani, Allerton, 2010.

List of main files:
1. DemoFB.m - Wrapper containing the real video foreground background separation. The sparse recovery step here uses ell-1 minimization. 
2. DemoDynRPCA.m - Wrapper containing the simulated data experiments. The sparse recovery step here uses CoSAMP. We can use ell-1 too if necessary.
3. AutoReProCS - main function which implements the ReProCS algororithm using CoSaMP.
4. ReProCS_real - main function which implements the ReProCS algorithm using ell-1.

Folders:
Data folder (120 MB) contains .mat files for videos. Pushing this too, because the original source web-page (http://perception.i2r.a-star.edu.sg/bk_model/bk_index.html) is down.
YALL1 - folder containing files to implement ell-1 minimization.


Other files:

ncrpca -- code implemented Non-convex Robust PCA, NIPS 14 downloaded from authors' website and its accompaniments
cgls -- fast method to implement least squares


For any further questions/suggestions please contact me @ pkurpadn iastate edu (insert obvious symbols)

