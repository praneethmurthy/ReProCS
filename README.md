This folder contains the code accompanying the following paper.

	[1] "Provable Dynamic Robust PCA or Robust Subspace Tracking", Praneeth Narayanamurthy and Namrata Vaswani, IEEE Trans. Info. Theory, 2019 (available at arXiv:1705.08948).

List of main files:

	1. DemoFB.m - Wrapper containing the real video foreground background separation. The sparse recovery step here uses ell-1 minimization. 
	1. DemoDynRPCA.m - Wrapper containing the simulated data experiments. The sparse recovery step here uses CoSAMP. We can use ell-1 too if necessary.
	1. AutoReProCS - main function which implements the ReProCS algororithm using CoSaMP.
	1. ReProCS_real - main function which implements the ReProCS algorithm using ell-1.

Folders:

	1. Data folder (120 MB) contains .mat files for videos. Pushing this too, because the original source web-page (http://perception.i2r.a-star.edu.sg/bk_model/bk_index.html) is down.
	1. YALL1 - folder containing files to implement ell-1 minimization.


Helper files:

	1. ncrpca -- code implemented Non-convex Robust PCA, NIPS 14 downloaded from authors' website and its accompaniments
	1. cgls -- fast method to implement least squares


For any further questions/suggestions please contact me @ pkurpadn iastate edu (insert obvious symbols)

