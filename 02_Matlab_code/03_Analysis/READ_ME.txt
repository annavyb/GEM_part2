###################

ALGORITHM ( version of the 2018/06/05 —> implementing version 2018/06/05 see algoV1)

(1) Set the GFT bandwidth (the number of the GFT frequencies to consider), W
(2) Set the size of the ROI (the number of the electrodes around the electrode with the highest spike amplitude), Ns
(3) Calculate the Shannon number, i.e. K = Ns*W/N, where N is the total number of the electrodes (204 in our case)
(4) Knowing the Shannon number, we have the information about number of the Slepians with the optimal energy concentration within he region of interest, therefore, we can construct the matrix Sopt (K x N).
(5) Project the resting state data matrix X (NxNt), where Nt is the number of time points), Xsl = SoptX of dimensions (K x Nt) 
(6) Compute the covariance matrix of the matrix Xsl, Csl = Xsl(Xsl)^t
(7) The eigenvectors of the Csl = QM(Q)^(-1) describe the axes of the hyperellipsoid within the data
(8) Obtain the distribution of the Malahanobis distance of each resting state data time point 
(9) Perform the step (5) for the spike data Xspike, then project the data into the space spanned by the eigenvectors of the covariance matrix Csl, Q. Xspike_sl = QXspike (Ns x Nt_spike)
(10) Calculate the Malahanobis distance of each timepoint in the spike data
(11) Define the alpha level and identify which are the time points of the transformed spike data falling within the defined alpha level of the distribution obtained with the resting state data 

#################

