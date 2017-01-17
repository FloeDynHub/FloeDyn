/* Siconos-Numerics, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*!
 * \file floe/lcp/solver/lexicolemke.hpp
 * \brief LCP Solver using Lemke algorithm with lexicographical ordering.
 * \author Roland Denis
 */

#include "floe/lcp/solver/lexicolemke.hpp"

#include <cmath>
#include <cassert>
#include <cstdio>
#include <cstdlib>

namespace floe { namespace lcp { namespace solver
{

//! Frontend specialization for double fundamental type.
template <>
bool lexicolemke<double>( floe::lcp::LCP<double>& lcp)
{
    int info;
    lcp_lexicolemke( 
        lcp.dim, 
        lcp.A.data().begin(), 
        lcp.q.data().begin(), 
        lcp.z.data().begin(),
        lcp.w.data().begin(), 
        &info 
    );
    return info == 0;
}

void lcp_lexicolemke(int dim, const double * M, const double * q, double *zlem , double *wlem , int *info)
{
    double tol = 0;

    /* matrix M of the lcp */
    //double * M = problem->M->matrix0;
    assert(M);
    /* size of the LCP */
    //int dim = problem->size;
    assert(dim>0);
    int dim2 = 2 * (dim + 1);

    int i, drive, block, Ifound;
    int ic, jc;
    int ITER;
    int nobasis;
    //int itermax = options->iparam[0];
    int itermax = std::min(40*dim, 2000);

    i=0;
    int n = dim;
    //double *q = problem->q;

    while ((i < (n - 1)) && (q[i] >= 0.)) 
        i++;

    if ((i == (n - 1)) && (q[n - 1] >= 0.))
    {

        /* TRIVIAL CASE : q >= 0
         * z = 0 and w = q is solution of LCP(q,M)
         */
        for (int j = 0 ; j < n; j++)
        {
            zlem[j] = 0.0;
            wlem[j] = q[j];
        }
        *info = 0;
        //options->iparam[1] = 0;   /* Number of iterations done */
        //options->dparam[1] = 0.0; /* Error */
        //if (verbose > 0)
        //printf("lcp_lexicolemke: found trivial solution for the LCP (positive vector q => z = 0 and w = q). \n");
        return ;
    }

    double z0, zb, dblock;
    double pivot, tovip;
    double tmp;
    int *basis;
    double** A;

    /*output*/

    //options->iparam[1] = 0;

    /* Allocation */

    basis = (int *)malloc(dim * sizeof(int));

    /*
    A = (double **)malloc(dim * sizeof(double*));
    for (ic = 0 ; ic < dim; ++ic)
    A[ic] = (double *)malloc(dim2 * sizeof(double));
    */

    // Better allocation
    A = (double **) malloc( dim * sizeof(double*) );
    A[0] = (double *) malloc( dim*dim2 * sizeof(double) );
    for (ic = 1; ic < dim; ++ic)
        A[ic] = A[ic-1] + dim2;


    /* construction of A matrix such that
     * A = [ q | Id | -d | -M ] with d = (1,...1)
     */

    /* We need to init only the part corresponding to Id */
    for (ic = 0 ; ic < dim; ++ic)
        for (jc = 1 ; jc <= dim; ++jc)
            A[ic][jc] = 0.0;

    for (ic = 0 ; ic < dim; ++ic)
        for (jc = 0 ; jc < dim; ++jc)
            A[ic][jc + dim + 2] = -M[dim * ic + jc];

    assert(q);

    for (ic = 0 ; ic < dim; ++ic) A[ic][0] = q[ic];

    for (ic = 0 ; ic < dim; ++ic) A[ic][ic + 1 ] =  1.0;
    for (ic = 0 ; ic < dim; ++ic) A[ic][dim + 1] = -1.0;

    /* End of construction of A */

    Ifound = 0;


    for (ic = 0 ; ic < dim  ; ++ic) basis[ic] = ic + 1;

    drive = dim + 1;
    block = 0;
    z0 = A[block][0];
    ITER = 0;

    /* Start research of argmin lexico */
    /* With this first step the covering vector enter in the basis */

    for (ic = 1 ; ic < dim ; ++ic)
    {
        zb = A[ic][0];
        if (zb < z0)
        {
            z0    = zb;
            block = ic;
        }
        else if (zb == z0)
        {
            for (jc = 0 ; jc < dim ; ++jc)
            {
                dblock = A[block][1 + jc] - A[ic][1 + jc];
                if (dblock < 0)
                {
                    break;
                }
                else if (dblock > 0)
                {
                    block = ic;
                    break;
                }
            }
        }
    }

    /* Stop research of argmin lexico */

    pivot = A[block][drive];
    tovip = 1.0 / pivot;

    /* Pivot < block , drive > */

    A[block][drive] = 1;
    for (ic = 0       ; ic < drive ; ++ic) A[block][ic] = A[block][ic] * tovip;
    for (ic = drive + 1 ; ic < dim2  ; ++ic) A[block][ic] = A[block][ic] * tovip;

    /* */

    for (ic = 0 ; ic < block ; ++ic)
    {
        tmp = A[ic][drive];
        for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
    }
    for (ic = block + 1 ; ic < dim ; ++ic)
    {
        tmp = A[ic][drive];
        for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
    }

    nobasis = basis[block];
    basis[block] = drive;

    while (ITER < itermax && !Ifound)
    {

        ++ITER;

        if (nobasis < dim + 1)      drive = nobasis + (dim + 1);
        else if (nobasis > dim + 1) drive = nobasis - (dim + 1);

        /* Start research of argmin lexico for minimum ratio test */
        pivot = 1e20;
        block = -1;

        for (ic = 0 ; ic < dim ; ++ic)
        {
            zb = A[ic][drive];
            if (zb > 0.0)
            {
                z0 = A[ic][0] / zb;
                if (z0 > pivot+tol) continue;
                if (z0 < pivot-tol)
                {
                    pivot = z0;
                    block = ic;
                }
                else
                {
                    for (jc = 1 ; jc < dim + 1 ; ++jc)
                    {
                        assert(block >=0 && "lcp_lexicolemke: block <0");
                        dblock = A[block][jc] / pivot - A[ic][jc] / zb;
                        if (dblock < 0.0-tol) break;
                        else if (dblock > 0.0+tol)
                        {
                            block = ic;
                            break;
                        }
                    }
                }
            }
        }
        if (block == -1)
        {
            Ifound = 0;
            /*
            DEBUG_PRINT("The pivot column is nonpositive !\n"
                    "It either means that the algorithm failed or that the LCP is infeasible\n"
                    "Check the class of the M matrix to find out the meaning of this\n");
            */
            break;
        }

        if (basis[block] == dim + 1) Ifound = 1;

        /* Pivot < block , drive > */

        pivot = A[block][drive];
        tovip = 1.0 / pivot;
        A[block][drive] = 1;

        for (ic = 0       ; ic < drive ; ++ic) A[block][ic] = A[block][ic] * tovip;
        for (ic = drive + 1 ; ic < dim2  ; ++ic) A[block][ic] = A[block][ic] * tovip;

        /* */

        for (ic = 0 ; ic < block ; ++ic)
        {
            tmp = A[ic][drive];
            for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
            //for (jc = 0; jc < dim2; ++jc) A[0][ic*dim2+jc] -= tmp * A[0][block*dim2+jc];
            
        }
        for (ic = block + 1 ; ic < dim ; ++ic)
        {
            tmp = A[ic][drive];
            for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
            //for (jc = 0; jc < dim2; ++jc) A[0][ic*dim2+jc] -= tmp * A[0][block*dim2+jc];
        }

        nobasis = basis[block];
        basis[block] = drive;

    } /* end while*/

    for (ic = 0 ; ic < dim; ++ic)
    {
        drive = basis[ic];
        if (drive < dim + 1)
        {
            zlem[drive - 1] = 0.0;
            wlem[drive - 1] = A[ic][0];
        }
        else if (drive > dim + 1)
        {
            zlem[drive - dim - 2] = A[ic][0];
            wlem[drive - dim - 2] = 0.0;
        }
    }

    //options->iparam[1] = ITER;

    if (Ifound) *info = 0;
    else *info = 1;

    free(basis);

    free(A[0]);
    free(A);
}

/*
   int linearComplementarity_lexicolemke_setDefaultSolverOptions(SolverOptions* options)
   {
   if (verbose > 0)
   {
   printf("Set the Default SolverOptions for the Lemke Solver\n");
   }

   options->solverId = SICONOS_LCP_LEMKE;
   options->numberOfInternalSolvers = 0;
   options->isSet = 1;
   options->filterOn = 1;
   options->iSize = 5;
   options->dSize = 5;
   options->iparam = (int *)calloc(options->iSize, sizeof(int));
   options->dparam = (double *)calloc(options->dSize, sizeof(double));
   options->dWork = NULL;
   options->iWork = NULL;   options->callback = NULL; options->numericsOptions = NULL;
   options->dparam[0] = 1e-6;
   options->iparam[0] = 10000;
   return 0;
   }
   */

}}} // namespace floe::lcp::solver

