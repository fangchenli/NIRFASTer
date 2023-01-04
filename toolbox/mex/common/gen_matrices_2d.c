#include <stdio.h>
#include <math.h>

#include "mex.h"
#include "matrix.h"

/* Calculates local matrices for FEM process to make overall matrix
*
* Part of NIRFAST package
* H Dehghani
* Last Updated - 7/13/09 M Jermyn
*/

void gg_and_det(double g[3][2], double gg[3][3], double *det);
void gradphi(double gg[3][3], double kappa[3], double det, double val[3][3]);
void phidotphi(double g[3][3], double mua[3], double det,double val[3][3]);
void boundary_int(double g[3][2], mxUint8 bnd_index[3],	\
                  double ksir[3], double valr[2][2]);
void bound2(double g[2][2], int il, int im, double ksir[2], double *val);

/* -------- Heart of the mex file----------- */
void mainloop(double *nodes,
              mxUint32 *elements,
              mxUint8 *bndvtx,
              double *mua,
              double *kappa,
              double *ksir,
              double *c,
              double omega,
              mwSize nodem,
              mwSize noden,
              mwSize elemm,
              mwSize elemn,
              mxUint32 *I,
              mxUint32 *J,
              mxComplexDouble *Sc)
{
    mwIndex ele, i, j, k, index[3];
    double tri_vtx[3][2];
    double kappa_index[3], mua_index[3], c_index[3];
    mxUint8 bnd_index[3];
    double k_val[3][3], c_val[3][3], b_val[3][3], f_valr[2][2], GG[3][3];
    double ksir_index[3];
    double det;

    k = 0;
    for (ele=0; ele<elemm; ++ele){
        for (i=0; i<elemn; ++i){
            index[i] = *(elements+(ele+(i*elemm)));;
        }

        for (i=0; i<elemn; ++i){
            kappa_index[i] = *(kappa+(index[i]-1));
            mua_index[i] = *(mua+(index[i]-1));
            c_index[i] = omega / *(c+(index[i]-1)) * -1;
            bnd_index[i] = *(bndvtx+(index[i]-1));
            ksir_index[i] = *(ksir+(index[i]-1));
            for (j=0; j<noden; ++j){
                tri_vtx[i][j] = *(nodes+(index[i]-1+(j*nodem)));
            }
        }
        gg_and_det(tri_vtx, GG, &det);
        gradphi(GG, kappa_index, det, k_val);
        phidotphi(GG, mua_index, det, c_val);
        phidotphi(GG, c_index, det, b_val);

        for (i=0; i<3; ++i){
            for (j=0; j<3; ++j){
                I[k] = index[i];
                J[k] = index[j];
                Sc[k].real = k_val[j][i]+c_val[j][i];
                Sc[k].imag = b_val[j][i];
                ++k;
            }
        }


        if (bnd_index[0]+bnd_index[1]+bnd_index[2] != 0){
            boundary_int(tri_vtx, bnd_index, ksir_index, f_valr);
            if (bnd_index[0]==1 && bnd_index[2]==1){
                int index1[2];
                index1[0] = index[0];
                index1[1] = index[2];
                for (i=0; i<2; ++i){
                    for (j=0; j<=i; ++j){
                        I[k] = index1[i];
                        J[k] = index1[j];
                        Sc[k].real = f_valr[i][j];
                        Sc[k].imag = 0.0;
                        ++k;
                    }
                }
            }
            else if (bnd_index[0]==1 && bnd_index[1]==1){
                int index1[2];
                index1[0] = index[0];
                index1[1] = index[1];
                for (i=0; i<2; ++i){
                    for (j=0; j<=i; ++j){
                        I[k] = index1[i];
                        J[k] = index1[j];
                        Sc[k].real = f_valr[i][j];
                        Sc[k].imag = 0.0;
                        ++k;
                    }
                }
            }
            else if (bnd_index[1]==1 && bnd_index[2]==1){
                int index1[2];
                index1[0] = index[1];
                index1[1] = index[2];
                for (i=0; i<2; ++i){
                    for (j=0; j<=i; ++j){
                        I[k] = index1[i];
                        J[k] = index1[j];
                        Sc[k].real = f_valr[i][j];
                        Sc[k].imag = 0.0;
                        ++k;
                    }
                }
            }
        }


    }
    return;
}

void gg_and_det(double g[3][2], double gg[3][3], double *det)
{
    static int L[2][3] = {{-1.0,1.0,0.0}, {-1.0,0.0,1.0}};
    double Jt[2][2], dJt, iJt[2][2], S[3], G[2][3], GG[3][3];

    Jt[0][0] = L[0][0]*g[0][0] + L[0][1]*g[1][0] + L[0][2]*g[2][0];
    Jt[0][1] = L[0][0]*g[0][1] + L[0][1]*g[1][1] + L[0][2]*g[2][1];
    Jt[1][0] = L[1][0]*g[0][0] + L[1][1]*g[1][0] + L[1][2]*g[2][0];
    Jt[1][1] = L[1][0]*g[0][1] + L[1][1]*g[1][1] + L[1][2]*g[2][1];

    dJt = (Jt[0][0]*Jt[1][1]) - (Jt[0][1]*Jt[1][0]);

    iJt[0][0] = +(Jt[1][1])/dJt;
    iJt[0][1] = -(Jt[0][1])/dJt;
    iJt[1][0] = -(Jt[1][0])/dJt;
    iJt[1][1] = +(Jt[0][0])/dJt;

    dJt = sqrt(dJt*dJt);

    G[0][0] = iJt[0][0]*L[0][0] + iJt[0][1]*L[1][0];
    G[0][1] = iJt[0][0]*L[0][1] + iJt[0][1]*L[1][1];
    G[0][2] = iJt[0][0]*L[0][2] + iJt[0][1]*L[1][2];
    G[1][0] = iJt[1][0]*L[0][0] + iJt[1][1]*L[1][0];
    G[1][1] = iJt[1][0]*L[0][1] + iJt[1][1]*L[1][1];
    G[1][2] = iJt[1][0]*L[0][2] + iJt[1][1]*L[1][2];

    GG[0][0] = G[0][0]*G[0][0] + G[1][0]*G[1][0];
    GG[0][1] = G[0][0]*G[0][1] + G[1][0]*G[1][1];
    GG[0][2] = G[0][0]*G[0][2] + G[1][0]*G[1][2];
    GG[1][0] = G[0][1]*G[0][0] + G[1][1]*G[1][0];
    GG[1][1] = G[0][1]*G[0][1] + G[1][1]*G[1][1];
    GG[1][2] = G[0][1]*G[0][2] + G[1][1]*G[1][2];
    GG[2][0] = G[0][2]*G[0][0] + G[1][2]*G[1][0];
    GG[2][1] = G[0][2]*G[0][1] + G[1][2]*G[1][1];
    GG[2][2] = G[0][2]*G[0][2] + G[1][2]*G[1][2];

    *det = dJt;
    return;
}

void gradphi(double GG[3][3], double kappa[3], double dJt, double val[3][3])
{
    int ii, jj;
    static double ip[3][2] = {{1.0/2.0,0.0},{1.0/2.0,1.0/2.0},{0.0,1.0/2.0}};
    static double w = 1 / 6;
    double S[3], tmp;

    for (ii=0; ii<3; ++ii){
        S[0] = 1-ip[ii][0]-ip[ii][1];
        S[1] = ip[ii][0];
        S[2] = ip[ii][1];
        tmp = kappa[0]*S[0] + kappa[1]*S[1] + kappa[2]*S[2];
        tmp *= w * dJt;
        val[0][0] += GG[0][0]*tmp;
        val[0][1] += GG[1][0]*tmp;
        val[0][2] += GG[2][0]*tmp;
        val[1][0] += GG[0][1]*tmp;
        val[1][1] += GG[1][1]*tmp;
        val[1][2] += GG[2][1]*tmp;
        val[2][0] += GG[0][2]*tmp;
        val[2][1] += GG[1][2]*tmp;
        val[2][2] += GG[2][2]*tmp;
    }
    return;
}

void phidotphi(double GG[3][3], double mua[3], double dJt, double val[3][3])
{
    int ii, jj;
    static double ip[3][2] = {{1.0/2.0,0.0},{1.0/2.0,1.0/2.0},{0.0,1.0/2.0}};
    static double w = 1 / 6;
    double S[3], tmp;

    for (ii=0; ii<3; ++ii){
        S[0] = 1-ip[ii][0]-ip[ii][1];
        S[1] = ip[ii][0];
        S[2] = ip[ii][1];
        tmp = mua[0]*S[0] + mua[1]*S[1] + mua[2]*S[2];
        tmp *= w * dJt;
        val[0][0] += S[0]*S[0]*tmp;
        val[0][1] += S[1]*S[0]*tmp;
        val[0][2] += S[2]*S[0]*tmp;
        val[1][0] += S[0]*S[1]*tmp;
        val[1][1] += S[1]*S[1]*tmp;
        val[1][2] += S[2]*S[1]*tmp;
        val[2][0] += S[0]*S[2]*tmp;
        val[2][1] += S[1]*S[2]*tmp;
        val[2][2] += S[2]*S[2]*tmp;
    }
    return;
}


void boundary_int(double gg[3][2], mxUint8 bnd_index[3], double ksir[3], \
                  double valr[2][2])
{
    double ksirtmp[2];
    int i, j;
    double g[2][2], val;

    for (i=0; i<2; ++i){
        for (j=0; j<2; ++j){
            valr[i][j]=0;
        }
    }

    if (bnd_index[0]==1 && bnd_index[2]==1){
        ksirtmp[0] = ksir[0];
        ksirtmp[1] = ksir[2];
        g[0][0] = gg[0][0];
        g[0][1] = gg[0][1];
        g[1][0] = gg[2][0];
        g[1][1] = gg[2][1];
        for (i=0; i<2; ++i){
            for (j=0; j<=i; ++j){
                bound2(g, i, j, ksirtmp, &val);
                valr[i][j] = val;
            }
        }
    }
    else if (bnd_index[0]==1 && bnd_index[1]==1){
        ksirtmp[0] = ksir[0];
        ksirtmp[1] = ksir[1];
        g[0][0] = gg[0][0];
        g[0][1] = gg[0][1];
        g[1][0] = gg[1][0];
        g[1][1] = gg[1][1];

        for (i=0; i<2; ++i){
            for (j=0; j<=i; ++j){
                bound2(g, i, j, ksirtmp, &val);
                valr[i][j] = val;
            }
        }
    }
    else if (bnd_index[1]==1 && bnd_index[2]==1){
        ksirtmp[0] = ksir[1];
        ksirtmp[1] = ksir[2];
        g[0][0] = gg[1][0];
        g[0][1] = gg[1][1];
        g[1][0] = gg[2][0];
        g[1][1] = gg[2][1];

        for (i=0; i<2; ++i){
            for (j=0; j<=i; ++j){
                bound2(g, i, j, ksirtmp, &val);
                valr[i][j] = val;
            }
        }
    }
    return;
}

void bound2(double g[2][2], int il, int im, double ksir[2], \
            double *val)
{
    static double w[2] = {1.0/2.0,1.0/2.0};
    static double ip[2] = {0.21132486540519,0.78867513459481};
    double dJt, S[2], tmp, tmpval;
    int ii;
    double x0 = g[0][0]; double x1 = g[1][0];
    double y0 = g[0][1]; double y1 = g[1][1];
    double k0 = ksir[0]; double k1 = ksir[1];

    dJt = sqrt(((g[1][0]-g[0][0])*(g[1][0]-g[0][0])) +	\
               ((g[1][1]-g[0][1])*(g[1][1]-g[0][1])));


    tmpval = 0.0;
    for (ii=0; ii<2; ++ii){
        S[0] = 1-ip[ii];
        S[1] = ip[ii];
        tmp = k0*S[0] + k1*S[1];
        tmp = tmp*w[ii];
        tmpval += S[il]*S[im]*tmp;
    }
    *val = tmpval*dJt;
    return;
}
/* -------- Gate-way to matlab  ------------ */

void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])

{
    double *nodes, *mua, *kappa, *ksir, *c, omega;
    mxUint8 *bndvtx;
    mxUint32 *elements;
    mwSize nodem, noden, elemm, elemn, nzmax;
    mxUint32 *I, *J;
    mxComplexDouble *S;

    /* Error checking  */

    if (nrhs < 8 )
        mexErrMsgTxt(" There is not enough input arguments");

    if(nlhs!=3)
        mexErrMsgTxt("This routine requires Three ouput arguments");

    nodes = mxGetDoubles(prhs[0]);          /*  nodes of mesh */
    elements = mxGetUint32s(prhs[1]);       /*  elements of mesh */
    bndvtx = mxGetUint8s(prhs[2]);          /*  boundary nodes */
    mua = mxGetDoubles(prhs[3]);            /*  absorption, nosal */
    kappa = mxGetDoubles(prhs[4]);          /*  diffusion, nodeal */
    ksir = mxGetDoubles(prhs[5]);           /*  reflection parameter, nodal  */
    c = mxGetDoubles(prhs[6]);              /*  speed of light in tissue, nodal */
    omega = mxGetScalar(prhs[7]);           /*  frequency of excitation */

    nodem=mxGetM(prhs[0]);          /*  Number of of nodes */
    noden=mxGetN(prhs[0]);          /*  Number of node freedom */
    elemm=mxGetM(prhs[1]);          /*  Number of elements */
    elemn=mxGetN(prhs[1]);          /*  Number of nodes per elements */

    nzmax = elemm*(elemn*(elemn+1));
    /*nzmax = nodem*nodem*0.1;*/
    /*mexPrintf("%d\n",nzmax);*/

    plhs[0] = mxCreateDoubleMatrix(nzmax, 1, mxREAL);    /* i */
    plhs[1] = mxCreateDoubleMatrix(nzmax, 1, mxREAL);    /* j */
    plhs[2] = mxCreateDoubleMatrix(nzmax, 1, mxCOMPLEX); /* s */
    I = mxGetUint32s(plhs[0]);
    J = mxGetUint32s(plhs[1]);
    S = mxGetComplexDoubles(plhs[2]);
    mainloop(nodes,elements,bndvtx,mua,kappa,ksir,c,omega,nodem,noden,elemm,elemn, I, J, S);
    return;
}
