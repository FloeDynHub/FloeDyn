import siconos.numerics as sn
import scipy.io as sio
import scipy.linalg as sla


lcpmat = sio.loadmat('./mat_jKtlx_high_requirement_new_perturbs_old_solvers.mat')
tol = 1e-12

def run_case(id, solver_id, solver_func):
    MM =lcpmat['M']
    qq = lcpmat['Q']
    M = MM[id,0]
    q = qq[id,0][:,0]
    # construction du LCP
    lcp = sn.LCP(M, q)
    # Definition des options
    so = sn.SolverOptions(lcp, solver_id)
    # Vecteurs vides
    z = np.zeros((lcp.size,), dtype=np.float64)
    w = np.zeros((lcp.size,), dtype=np.float64)
    # Appel du solveur
    info = solver_func(lcp, z,w,so)
    
    wold = w.copy()

    # Calcul d'erreur 
    error = sn.lcp_compute_error(lcp, z, w, tol)
    error2 = sla.norm((w-wold))
    return error, error2


## Exemple :
# run_case(1160, sn.SICONOS_LCP_LEMKE, sn.lcp_lexicolemke)

def run_all_solver(id, solver_id):
    MM =lcpmat['M']
    qq = lcpmat['Q']
    M = MM[id,0]
    q = qq[id,0][:,0]
    # construction du LCP
    lcp = sn.LCP(M, q)
    # Definition des options
    so = sn.SolverOptions(lcp, solver_id)
    so.filterOn = 0 # pour ne pas appeler lcp_compute_error
    # Vecteurs vides
    z = np.zeros((lcp.size,), dtype=np.float64)
    w = np.zeros((lcp.size,), dtype=np.float64)
    # Appel du solveur
    info = sn.lcp_driver_DenseMatrix(lcp, z, w, so)

    wold = w.copy()

    # Calcul d'erreur 
    error = sn.lcp_compute_error(lcp, z, w, tol)
    error2 = sla.norm((w-wold))
    return error, error2



for i in range(2100):
    error, error2 = run_all_solver(i, sn.SICONOS_LCP_LEMKE)
    if error > 1e-6:
        print(i, error, error2)

