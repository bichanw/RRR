import numpy as np
# not tested!!!


def alignment_output(X, Y, W):
    """
    Compute the alignment index of the output population activities with communication subspace.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_input_features)
        Input activities.
    Y : ndarray of shape (n_samples, n_output_features)
        Output activities.
    W : ndarray of shape (n_input_features, n_output_features)
        Communication weights.

    Returns
    -------
    outputalignmentidx : float
        Normalized alignment index (between -1 and 1).
    """

    # ----- calculate raw alignment index -----
    # Do PCA on population covariance
    upop, spopvec, _ = np.linalg.svd(np.cov(Y, rowvar=False))
    spopcum = np.cumsum(spopvec)
    spopvec_nrm = spopvec / np.sum(spopvec)
    spopcum_nrm = np.cumsum(spopvec_nrm)

    # Project communication covariance onto PCs
    cov_predicted = np.cov(X @ W, rowvar=False)
    scomvec = np.diag(upop.T @ cov_predicted @ upop)
    # scomvec_nrm = scomvec / np.sum(scomvec)
    scomcum_nrm = np.cumsum(scomvec) / np.sum(scomvec)

    # Compute alignment index
    muscom = np.mean(scomcum_nrm)
    muspop = np.mean(spopcum_nrm)
    alignment_raw = np.dot(scomvec, spopvec) # dot product of scomvec and spopvec

    # Compute communication fraction
    commfrac = np.sum(scomvec) / np.sum(spopvec)

    # ----- normalize alignment index -----
    totcom = np.sum(scomvec)

    # Max possible alignment
    ii = np.argmax(spopcum > totcom + 1e-10)
    scommax = spopvec.copy()
    scommax[ii:] = 0
    scommax[ii] = totcom - np.sum(scommax[:ii])
    scommax_cum = np.cumsum(scommax) / np.sum(scommax)
    a_max = np.dot(scommax, spopvec)

    # Min possible alignment (flip order of eigenvalues)
    spopvec_rev = np.flipud(spopvec)
    spopcum_rev = np.cumsum(spopvec_rev)
    ii = np.argmax(spopcum_rev > totcom + 1e-10)
    scommin = spopvec_rev.copy()
    scommin[ii:] = 0
    scommin[ii] = totcom - np.sum(scommin[:ii])
    scommin = np.flipud(scommin)  # flip back
    scommin_cum = np.cumsum(scommin) / np.sum(scommin)
    a_min = np.dot(scommin, spopvec)

    # Rescale alignment
    outputalignmentidx = (alignment_raw - a_min) / (a_max - a_min)

    return outputalignmentidx, commfrac


def alignment_input(X, W, r=None, C=None):
    """
    Calculate how much the communication weights W align with 
    the principal components of the input X.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_input_neurons)
        Input matrix (stimuli).
    W : ndarray of shape (n_input_neurons, n_output_neurons)
        Communication weights.
    r : int, optional
        Rank of the communication weights W (if not provided,
        estimate from W).
    C : ndarray, optional
        Covariance matrix of X. If None, will compute cov(X).

    Returns
    -------
    aa : float
        Alignment index (0-1), where 1 is maximally aligned.
    p : None
        Placeholder (for compatibility with MATLAB signature).
    aa_rand : None
        Placeholder (for compatibility with MATLAB signature).
    """
    # Covariance of inputs
    if C is None:
        C = np.cov(X, rowvar=False)

    # PCA of covariance matrix
    _, Spcavec, _ = np.linalg.svd(C)

    # SVD of weights
    _, Swvec, _ = np.linalg.svd(W)
    # Pad singular values to length n_input_neurons
    Swvec_padded = np.concatenate([Swvec, np.zeros(W.shape[0] - len(Swvec))])

    # Compute alignment index
    amax = Spcavec @ (Swvec_padded ** 2)             # maximal value
    amin = Spcavec @ (np.flipud(Swvec_padded ** 2))  # minimal value
    araw = np.trace(W.T @ C @ W)                     # test statistic

    aa = (araw - amin) / (amax - amin)

    # Match MATLAB output signature
    return aa, None, None
