This folder contains the code accompanying the paper 

[1] "New Results for Provable Dynamic Robust PCA", Praneeth Narayanamurthy and Namrata Vaswani, arXiv:1705.08948, 2017.

List of main files:
1. DemoFB.m - Wrapper containing the real video foreground background separation. The sparse recovery step here uses ell-1 minimization. 
2. DemoDynRPCA.m - Wrapper containing the simulated data experiments. The sparse recovery step here uses CoSAMP. We can use ell-1 too if necessary.
3. AutoReProCS - main function which implements the ReProCS algororithm using CoSaMP.
4. ReProCS_real - main function which implements the ReProCS algorithm using ell-1.

Folders:
Data folder (120 MB) contains .mat files for videos. Pushing this too, because the original source web-page is down.
YALL1 - folder containing files to implement ell-1 minimization.


Other files:

ncrpca -- code implemented Non-convex Robust PCA, NIPS 14 downloaded from authors' website and its accompaniments
cgls -- fast method to implement least squares


For any further questions/suggestions please contact me @ pkurpadn iastate edu (insert obvious symbols)

