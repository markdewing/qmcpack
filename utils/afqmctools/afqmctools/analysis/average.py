"""Simple extraction of afqmc rdms."""
import h5py
import numpy
import scipy.stats
from afqmctools.analysis.extraction import (
        get_metadata,
        extract_observable
        )


# enumish
WALKER_TYPE = ['undefined', 'closed', 'collinear', 'non_collinear']


def average_one_rdm(filename, estimator='back_propagated', eqlb=1, skip=1, ix=None):
    """Average AFQMC 1RDM.

    Returns P_{sij} = <c_{is}^+ c_{js}^> as a (nspin, M, M) dimensional array.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    eqlb : int
        Number of blocks for equilibration. Default 1.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).

    Returns
    -------
    one_rdm : :class:`numpy.ndarray`
        Averaged 1RDM.
    one_rdm_err : :class:`numpy.ndarray`
        Error bars for 1RDM elements.
    """
    md = get_metadata(filename)
    mean, err = average_observable(filename, 'one_rdm', eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    nbasis = md['NMO']
    wt = md['WalkerType']
    try:
        walker = WALKER_TYPE[wt]
    except IndexError:
        print('Unknown walker type {}'.format(wt))

    if walker == 'closed':
        return 2*mean.reshape((1,nbasis,nbasis)), err.reshape((1,nbasis, nbasis))
    elif walker == 'collinear':
        return mean.reshape((2,nbasis,nbasis)), err.reshape((2, nbasis, nbasis))
    elif walker == 'non_collinear':
        return mean.reshape((1,2*nbasis,2*nbasis)), err.reshape((1,2*nbasis, 2*nbasis))
    else:
        print('Unknown walker type.')
        return None

def average_diag_two_rdm(filename, estimator='back_propagated', eqlb=1, skip=1, ix=None):
    """Average diagonal part of 2RDM.

    Returns <c_{is}^+ c_{jt}^+ c_{jt} c_{is}> as a (2M,2M) dimensional array.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    eqlb : int
        Number of blocks for equilibration. Default 1.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).

    Returns
    -------
    two_rdm : :class:`numpy.ndarray`
        Averaged diagonal 2RDM.
    two_rdm_err : :class:`numpy.ndarray`
        Error bars for diagonal 2RDM elements.
    """
    md = get_metadata(filename)
    mean, err = average_observable(filename, 'diag_two_rdm', eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    nbasis = md['NMO']
    wt = md['WalkerType']
    try:
        walker = WALKER_TYPE[wt]
    except IndexError:
        print('Unknown walker type {}'.format(wt))

    if walker == 'closed':
        dm_size = nbasis*(2*nbasis-1) - nbasis*(nbasis-1) // 2
        assert mean.shape == dm_size
        two_rdm = numpy.zeros((2*nbasis, 2*nbasis), dtype=mean.dtype)
        two_rdm_err = numpy.zeros((2*nbasis, 2*nbasis), dtype=mean.dtype)
        ij = 0
        for i in range(nbasis):
            for j in range(i+1, 2*nbasis):
                two_rdm[i,j] = mean[ij]
                two_rdm_err[i,j] = err[ij]
                ij += 1
        two_rdm[nbasis:,nbasis:] = two_rdm[:nbasis,:nbasis].copy()
        two_rdm_err[nbasis:,nbasis:] = two_rdm_err[:nbasis,:nbasis].copy()
    elif walker == 'collinear':
        dm_size = nbasis*(2*nbasis-1)
        two_rdm = numpy.zeros((2*nbasis, 2*nbasis), dtype=mean.dtype)
        two_rdm_err = numpy.zeros((2*nbasis, 2*nbasis), dtype=mean.dtype)
        ij = 0
        for i in range(2*nbasis):
            for j in range(i+1, 2*nbasis):
                two_rdm[i,j] = mean[ij]
                two_rdm_err[i,j] = err[ij]
                ij += 1
    elif walker == 'non_collinear':
        print("Non-collinear wavefunction not supported.")
        return None
    else:
        print('Unknown walker type.')
        return None
    # Diagonal is zero
    two_rdm = 0.5 * (two_rdm + two_rdm.conj().T)
    two_rdm_err = 0.5 * (two_rdm_err + two_rdm_err.T)
    return two_rdm, two_rdm_err

def average_on_top_pdm(filename, estimator='back_propagated', eqlb=1, skip=1, ix=None):
    """Average on-top pair density matrix.

    Returns n(r,r) for a given real space grid.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    eqlb : int
        Number of blocks for equilibration. Default 1.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).

    Returns
    -------
    opdm : :class:`numpy.ndarray`
        Averaged diagonal on-top pair density matrix.
    opdm_err : :class:`numpy.ndarray`
        Error bars for diagonal on-top pair density matrix elements.
    """
    md = get_metadata(filename)
    mean, err = average_observable(filename, 'on_top_pdm', eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    # TODO: Update appropriately.
    return mean, err

def average_observable(filename, name, eqlb=1, estimator='back_propagated',
                       ix=None, skip=1):
    """Compute mean and error bar for AFQMC HDF5 observable.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    name : string
        Name of observable (see estimates.py for list).
    eqlb : int
        Number of blocks for equilibration. Default 1.
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).

    Returns
    -------
    mean : :class:`numpy.ndarray`
        Averaged quantity.
    err : :class:`numpy.ndarray`
        Error bars for quantity.
    """
    md = get_metadata(filename)
    free_proj = md['FreeProjection']
    if free_proj:
        mean = None
        err = None
        print("# Error analysis for free projection not implemented.")
    else:
        data = extract_observable(filename, name=name, estimator=estimator, ix=ix)
        mean = numpy.mean(data[eqlb:len(data):skip], axis=0)
        err = scipy.stats.sem(data[eqlb:len(data):skip].real, axis=0)
    return mean, err

def average_gen_fock(filename, fock_type='plus', estimator='back_propagated',
                     eqlb=1, skip=1, ix=None):
    """Average AFQMC genralised Fock matrix.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    fock_type : string
        Which generalised Fock matrix to extract. Optional (plus/minus).
        Default: plus.
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    eqlb : int
        Number of blocks for equilibration. Default 1.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).

    Returns
    -------
    gfock : :class:`numpy.ndarray`
        Averaged 1RDM.
    gfock_err : :class:`numpy.ndarray`
        Error bars for 1RDM elements.
    """
    md = get_metadata(filename)
    name = 'gen_fock_' + fock_type
    mean, err = average_observable(filename, name, eqlb=eqlb, skip=skip,
                                   estimator=estimator, ix=ix)
    nbasis = md['NMO']
    wt = md['WalkerType']
    try:
        walker = WALKER_TYPE[wt]
    except IndexError:
        print('Unknown walker type {}'.format(wt))

    if walker == 'closed':
        return 2*mean.reshape(1,nbasis,nbasis), err.reshape(1,nbasis, nbasis)
    elif walker == 'collinear':
        return mean.reshape((2,nbasis,nbasis)), err.reshape((2, nbasis, nbasis))
    elif walker == 'non_collinear':
        return mean.reshape((1,2*nbasis,2*nbasis)), err.reshape((1,2*nbasis, 2*nbasis))
    else:
        print('Unknown walker type.')
        return None

def get_noons(filename, estimator='back_propagated', eqlb=1, skip=1, ix=None,
              nsamp=20, screen_factor=1, cutoff=1e-14):
    """Get NOONs from averaged AFQMC RDM.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    eqlb : int
        Number of blocks for equilibration. Default 1.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).
    nsamp : int
        Number of perturbed RDMs to construct to estimate of error bar. Optional.
        Default: 20.
    screen_factor : int
        Zero RDM elements with abs(P[i,j]) < screen_factor*Perr[i,j]. Optional
        Default: 1.

    Returns
    -------
    noons : :class:`numpy.ndarray`
        NOONS.
    noons_err : :class:`numpy.ndarray`
        Estimate of error bar on NOONs.
    """
    P, Perr = average_one_rdm(filename, estimator='back_propagated', eqlb=1,
                              skip=1, ix=ix)
    if Perr.shape[0] == 2:
        # Collinear
        Perr = numpy.sqrt((Perr[0]**2 + Perr[1]**2))
    else:
        # Non-collinear / Closed
        Perr = Perr[0]
    # Sum over spin.
    P = numpy.sum(P, axis=0)
    P = 0.5 * (P + P.conj().T)
    Perr = 0.5 * (Perr + Perr.T)
    P[numpy.abs(P) < screen_factor*Perr] = 0.0
    noons = numpy.zeros((nsamp, P.shape[-1]))
    for s in range(nsamp):
        PT, X = regularised_ortho(P, cutoff=cutoff)
        Ppert = gen_sample_matrix(PT, Perr)
        e, ev = numpy.linalg.eigh(Ppert)
        noons[s] = e[::-1]

    PT, X = regularised_ortho(P, cutoff=cutoff)
    e, ev = numpy.linalg.eigh(PT)
    return e[::-1], numpy.std(noons, axis=0, ddof=1)

def analyse_ekt(filename, estimator='back_propagated', eqlb=1, skip=1, ix=None,
                nsamp=20, screen_factor=1, cutoff=1e-14):
    """Perform EKT analysis.

    Parameters
    ----------
    filename : string
        QMCPACK output containing density matrix (*.h5 file).
    estimator : string
        Estimator type to analyse. Options: back_propagated or mixed.
        Default: back_propagated.
    eqlb : int
        Number of blocks for equilibration. Default 1.
    skip : int
        Number of blocks to skip in between measurements equilibration.
        Default 1 (use all data).
    ix : int
        Back propagation path length to average. Optional.
        Default: None (chooses longest path).
    nsamp : int
        Number of perturbed RDMs to construct to estimate of error bar. Optional.
        Default: 20.
    screen_factor : int
        Zero RDM elements with abs(P[i,j]) < screen_factor*Perr[i,j]. Optional
        Default: 1.

    Returns
    -------
    [eip, eip_err] : list
        Ionisation potentials and estimates of their errors.
    [eea, eea_err] : list
        Electron affinities and estimates of their errors.
    """
    P, Perr = average_one_rdm(filename, estimator='back_propagated', eqlb=eqlb,
                              skip=skip, ix=ix)
    Fp, Fperr = average_gen_fock(filename, fock_type='plus',
                                 estimator='back_propagated', eqlb=eqlb,
                                 skip=skip, ix=ix)
    Fm, Fmerr = average_gen_fock(filename, fock_type='minus',
                                 estimator='back_propagated', eqlb=eqlb,
                                 skip=skip, ix=ix)
    P[numpy.abs(P) < screen_factor*Perr] = 0.0
    Fm[numpy.abs(Fm) < screen_factor*Fmerr] = 0.0
    Fp[numpy.abs(Fp) < screen_factor*Fperr] = 0.0
    # TODO : Not quite sure if this is the correct way to treat spin.
    # Ionisation potential
    if P.shape[0] == 2:
        P = P[0] + P[1]
        Perr = 0.5*Perr[0] + Perr[1]
        Fp = Fp[0] + Fp[1]
        Fm = Fm[0] + Fm[1]
        Fperr = Fperr[0] + Fperr[1]
        Fmerr = Fmerr[0] + Fmerr[1]
    else:
        P = P[0]
        Perr = Perr[0]
        Fp = Fp[0]
        Fm = Fm[0]
        Fperr = Fperr[0]
        Fmerr = Fmerr[0]
    P = 0.5*(P+P.conj().T)
    I = numpy.eye(P.shape[-1])
    gamma = P.real
    gamma_ = gamma.copy()
    # Diagonalise gamma and discard problematic singular values.
    eip, eip_vec = solve_gen_eig(gamma, Fm)
    eip_err = estimate_error_eig(gamma_, Perr, Fm, Fmerr, nsamp=nsamp, cutoff=cutoff)
    # Electron affinity.
    gamma = 2*I - P.T
    gamma_ = gamma.copy()
    eea, eaa_vec = solve_gen_eig(gamma, Fp)
    eea_err = estimate_error_eig(gamma_, Perr, Fp, Fperr, nsamp=nsamp, cutoff=cutoff)

    return (eip, eip_err), (eea, eea_err)

def solve_gen_eig(gamma, fock, cutoff=1e-14):
    gamma, X = regularised_ortho(gamma, cutoff=cutoff)
    fock_trans = numpy.dot(X.conj().T, numpy.dot(fock, X))
    return numpy.linalg.eigh(fock_trans)

def estimate_error_eig(gamma, gamma_err, fock, fock_err, nsamp=20, cutoff=1e-14):
    """Bootstrap estimate of error in eigenvalues."""
    eigs_tot = numpy.zeros((nsamp, gamma.shape[-1]))
    # TODO FIX THIS
    for s in range(nsamp):
        gamma_p = gen_sample_matrix(gamma, gamma_err)
        fock_p = gen_sample_matrix(fock, fock_err, herm=False)
        eigs, eigv = solve_gen_eig(gamma_p, fock_p)
        eigs_tot[s,:len(eigs)] = eigs
    return numpy.std(eigs_tot, axis=0, ddof=1)

def regularised_ortho(S, cutoff=1e-14):
    """Get orthogonalisation matrix."""
    sdiag, Us = numpy.linalg.eigh(S)
    sdiag[sdiag<cutoff] = 0.0
    X = Us[:,sdiag>cutoff] / numpy.sqrt(sdiag[sdiag>cutoff])
    Smod = numpy.dot(Us[:,sdiag>cutoff], numpy.diag(sdiag[sdiag>cutoff]))
    Smod = numpy.dot(Smod, Us[:,sdiag>cutoff].T.conj())
    return Smod, X

def gen_sample_matrix(mat, err, herm=True):
    """Perturb matrix by error bar."""
    a = numpy.zeros_like(mat)
    nbasis = mat.shape[0]
    assert a.shape[1] == nbasis
    assert a.shape == mat.shape
    # Dangerous. Discards imaginary part.
    for i in range(nbasis):
        for j in range(nbasis):
            a[i,j] = numpy.random.normal(loc=mat[i,j], scale=err[i,j], size=1)
    if herm:
        return 0.5*(a+a.conj().T)
    else:
        return a
