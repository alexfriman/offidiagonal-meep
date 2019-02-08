/* Copyright (C) 2005-2014 Massachusetts Institute of Technology
%(С) 2015-2018 Lebedev Physical Institute of  the Russian Academy of Sciences (I really do not understand those copyright laws. Hope this line does not break any law. If it does, email me: friman_a@sci.lebedev.ru)
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/
#include <stdlib.h>
#include <string.h>
#include "meep.hpp"
//#include "vec.hpp" //for is_magnetic, is_electric functions
#include <complex>
#include <cmath>
#include "meep_internals.hpp"

using namespace std;

namespace meep {
    extern void abort(const char *, ...); 
    extern "C" { //to avoid using CLAPACK, we use Fortran LAPACK since it is already required in ./configure
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
    
   
}

    void matrices_add(realnum* A, realnum alpha, realnum* B,  realnum beta, int k, int m, realnum* C){
        for (int i=0; i<k;i++){
            for (int j=0; j<m;j++){
                C[i+j*k]=alpha*A[i+j*k]+beta*B[i+j*k];
            }
        }
    }
    /*A simple function for row major matrix multiplication. It runs only 6 times in the constructor, 
     * and the matrices are 3*3. It is not an algorithm requiring a lot of optimization.*/
    void matrix_mult_rowM(realnum* A, realnum* B, realnum* C, int k, int m){
         for (int i=0; i<k;i++){  
             for (int j=0; j<m;j++){
                 C[i+j*k]=0;
                 for(int l=0; l<k;l++){
                     C[i+j*k]+=A[i+l*k]*B[l+j*k];
                 }
             }
         } 
    }
    
    void rowMx9_to_2D3x3(realnum rM[9], realnum twoD[3][3]){
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                twoD[i][j]=rM[i*3+j];
            }
        }
    }

    
     gyroelectric_susceptibility::gyroelectric_susceptibility(double in_nu_col, double in_nu_cycl, double b_dirx, double b_diry, double b_dirz)  {
         nu_col=in_nu_col;
         nu_cycl=in_nu_cycl;
         realnum b_norm=pow(b_dirx*b_dirx+b_diry*b_diry+b_dirz*b_dirz,0.5);
           if (b_norm!=1.0){
               b_dirx=b_dirx/b_norm;
               b_diry=b_diry/b_norm;
               b_dirz=b_dirz/b_norm;
           }
        realnum sinus=pow((b_diry*b_diry + b_dirx*b_dirx),0.5); // sine of angle between the new  B_z and the Z
        realnum cosine=b_dirz; //cosine of the same angle
        //projection matrices 
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++ ){
                Pr_diag_cart[i][j]=0;
                Pr_offdiag_cart[i][j]=0;
                Pr_z_cart[i][j]=0;
            }
        }
           
        //end of projection matrices
        //calculating the rotation matrix
        realnum ident_rowM[]={1, 0, 0, 0, 1,0,0,0,1};
        realnum R_rowM[]={0, 0, 0, 0, 0,0,0,0,0};
        realnum R_inv_rowM[9];
        //v=bvector x e_z=[by,-1*bx,0]; vmat=[[0,-v_z,v_y],[v_z,0,-v_x],[-v_y,v_x,0]]
        realnum vmat_rowM[]={0, 0,b_dirx,0, 0, 1.0*b_diry, -1.0*b_dirx, -1.0*b_diry, 0};
        realnum vmatsq_rowM[9];

        matrix_mult_rowM(vmat_rowM,vmat_rowM,vmatsq_rowM,3,3);

        if (b_dirx==0 and b_diry==0 and b_dirz==-1){
            R_rowM[0]=1;
            R_rowM[1+1*3]=-1;
            R_rowM[2+2*3]=-1;//{-1,0,0, 0,-1,0, 0,0,-1};
        }else if(b_dirx==0 and b_diry==0 and b_dirz==1){
            R_rowM[0]=1;
            R_rowM[1+1*3]=1;
            R_rowM[2+2*3]=1;//{1,0,0, 0,1,0, 0,0,1};
        }else{ //the general case
            matrices_add(vmatsq_rowM,(1-cosine)/(sinus*sinus),                       vmat_rowM,1,3,3,R_rowM);
            matrices_add(R_rowM,1,                       ident_rowM,1,3,3,R_rowM);
        }
        int pivotArray[4]; //since our matrix has three rows
        int errorHandler[1];
        realnum buf9[9]; // temporary buffer
        int N3=3;
        int N9=9;
        for (int i=0; i<9; i++) R_inv_rowM[i]=R_rowM[i];
        dgetrf_(&N3,&N3,R_inv_rowM,&N3,pivotArray,errorHandler);//calculate LU
        dgetri_(&N3,R_inv_rowM,&N3,pivotArray,buf9,&N9,errorHandler);//calculate inverse R matrix
        //calculation of R*eps_{xx}*(R^-1)
        realnum espxx_rowM[]={1, 0, 0, 0, 1,0,0,0,0};
          
        realnum Pr_diag_rowM[9]; // temporary buffer
        matrix_mult_rowM(espxx_rowM,R_inv_rowM,buf9,3,3);
        matrix_mult_rowM(R_rowM,buf9,Pr_diag_rowM,3,3);
        rowMx9_to_2D3x3(Pr_diag_rowM, Pr_diag_cart);
        //calculation of R*eps_{xy}*(R^-1)
        realnum espxy_rowM[]={0, 1, 0, -1, 0,0,0,0,0};
        realnum Pr_offdiag_rowM[9]; // temporary buffer
        matrix_mult_rowM(espxy_rowM,R_inv_rowM,buf9,3,3);
        matrix_mult_rowM(R_rowM,buf9,Pr_offdiag_rowM,3,3);
        rowMx9_to_2D3x3(Pr_offdiag_rowM, Pr_offdiag_cart);
        //calculation of R*eps_{zz}*(R^-1)
        realnum espzz_rowM[]={0, 0, 0, 0, 0,0,0,0,1};
        realnum Pr_z_rowM[9]; // temporary buffer
        matrix_mult_rowM(espzz_rowM,R_inv_rowM,buf9,3,3);
        matrix_mult_rowM(R_rowM,buf9,Pr_z_rowM,3,3);
        rowMx9_to_2D3x3(Pr_z_rowM, Pr_z_cart);

    }
     
    void gyroelectric_susceptibility::update_P(realnum *W[NUM_FIELD_COMPONENTS][2], 
			realnum *W_prev[NUM_FIELD_COMPONENTS][2], 
			double dt, const grid_volume &gv,
			void *P_internal_data) const{
        #define SWAP(t,a,b) { t SWAP_temp = a; a = b; b = SWAP_temp; }

     // stable averaging of offdiagonal components
    #define OFFDIAG(u,g,sx,s) (0.25 * ((g[i]+g[i-sx])*u[i]		\
		   	         + (g[i+s]+g[(i+s)-sx])*u[i+s]))
        
        gyroelectric_susceptibility_data *d_int = (gyroelectric_susceptibility_data *) P_internal_data;
        //realnum dt2=0.5*dt;
        complex<double> j(0,1);
        complex<double> exp_decay = d_int->exp_decay;
        complex<double> exp_decay_zz = d_int->exp_decay_zz;
        FOR_COMPONENTS(c) DOCMP2 if (d_int->P[c][cmp]) { // the main loop
            const realnum *w = W[c][cmp], *s = sigma[c][component_direction(c)];//field and sigma
            //const realnum *w_prev=W_prev[c][cmp];
            
            if (w && s){
                realnum *p = d_int->P[c][cmp], *pp = d_int->P_prev[c][cmp];// polarisation and previous value
                complex<double>* p_c = new complex<double>[ntot];//dynamic buffer for complex polarization data. TODO: remove it and sum only real parts
                LOOP_OVER_VOL_OWNED(gv, c, i) {
                    p_c[i]=0;
                }
                 FOR_COMPONENTS(c1){
                     if (d_int->Pr_diag[c][c1]!=0 or d_int->Pr_offdiag[c][c1]!=0 or d_int->Pr_z[c][c1]!=0){
                        const direction d = component_direction(c);//CHANGE OF NAME d !
                        const direction d1 = component_direction(c1);
                        const int is = gv.stride(d) * (is_magnetic(c) ? -1 : +1);
                        int is1 = gv.stride(d1) * (is_magnetic(c) ? -1 : +1);
                        const realnum *w1 = W[c1][cmp];//field
                        //const realnum *w1_prev = W_prev[c1][cmp];
                        if (d_int->Pr_diag[c][c1]!=0 ){ //projections of diagonal elemets
                            complex<double> *conv_diag1 = d_int->conv_diag1[c][c1][cmp];
                            complex<double> *conv_diag2 = d_int->conv_diag2[c][c1][cmp];
                            LOOP_OVER_VOL_OWNED(gv, c, i) {
                                if (s[i]!=0.0){
                                    conv_diag1[i]=conv_diag1[i]+d_int->chi0_diag1[i]*w1[i];
                                    conv_diag2[i]=exp_decay*conv_diag2[i]+d_int->chi0_diag2[i]*w1[i]; 
                                    p_c[i]+=d_int->Pr_diag[c][c1]*(conv_diag1[i]+conv_diag2[i]);
                                }
                                
                            }
                        }
                        if (d_int->Pr_offdiag[c][c1]!=0 ){ // of off-diagonal elements
                            complex<double> *conv_offdiag1 = d_int->conv_offdiag1[c][c1][cmp];
                            complex<double> *conv_offdiag2 = d_int->conv_offdiag2[c][c1][cmp];
                            LOOP_OVER_VOL_OWNED(gv, c, i) {
                                if (s[i]!=0.0){
                                    conv_offdiag1[i]=conv_offdiag1[i]+ d_int->chi0_offdiag1[i]*w1[i];
                                    conv_offdiag2[i]=exp_decay*conv_offdiag2[i]+ d_int->chi0_offdiag2[i]*w1[i]; 
                                    p_c[i]+=d_int->Pr_offdiag[c][c1]*(conv_offdiag1[i]+conv_offdiag2[i]);
                                }
                                
                            }
                        }
                         if (d_int->Pr_z[c][c1]!=0 ){ // of off-diagonal elements
                            complex<double> *conv_z1 = d_int->conv_z1[c][c1][cmp];
                            complex<double> *conv_z2 = d_int->conv_z2[c][c1][cmp];
                            LOOP_OVER_VOL_OWNED(gv, c, i) {
                                if (s[i]!=0.0){
                                    conv_z1[i]=conv_z1[i]+d_int->chi0_zz1[i]*w1[i];
                                    conv_z2[i]=exp_decay_zz*conv_z2[i]+d_int->chi0_zz2[i]*w1[i];
                                    p_c[i]+=d_int->Pr_z[c][c1]*(conv_z1[i]+conv_z2[i]);
                                }
                                
                            }
                        }
                       
                }
            }
            //we are out of c1 iteration, time to sum up the c component of polarization 
           LOOP_OVER_VOL_OWNED(gv, c, i) {
               realnum pcur = p[i];
               p[i]=p_c[i].real();
               pp[i] = pcur;
          }
           delete [] p_c;
        }
        }

        
    };
 
 void gyroelectric_susceptibility::subtract_P(field_type ft,
			  realnum *f_minus_p[NUM_FIELD_COMPONENTS][2], 
			  void *P_internal_data) const {
  gyroelectric_susceptibility_data *d = (gyroelectric_susceptibility_data *) P_internal_data;
  field_type ft2 = ft == E_stuff ? D_stuff : B_stuff; // for sources etc.
  int ntot = d->ntot;
  FOR_FT_COMPONENTS(ft, ec) DOCMP2 if (d->P[ec][cmp]) {
    component dc = field_type_component(ft2, ec);
    if (f_minus_p[dc][cmp]) {
      realnum *p = d->P[ec][cmp];
      realnum *fmp = f_minus_p[dc][cmp];
      for (int i = 0; i < ntot; ++i)
      { 
          fmp[i] -= p[i]; 
      }
    }
  }
 };
 
 
    void *gyroelectric_susceptibility::new_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2],
				  const grid_volume &gv) const {
      //first determine projections
      int num = 0;
      int conv_num=0;
      FOR_COMPONENTS(c) DOCMP2 if (needs_P(c, cmp, W)){ 
        num += 2 * gv.ntot(); //check that 2*
        
      }
      /*It's important. At this stage we do not care if
       * we have magnetic or electric polarization. All we
       * need to know is the number of non-zero projection.
       * Field assignation will be later
       */
      for (int i=0; i<3; i++){
            for (int j=0; j<3; j++ ){
                if (Pr_diag_cart[i][j]!=0) conv_num++;
                if (Pr_offdiag_cart[i][j]!=0) conv_num++;
                if (Pr_z_cart[i][j]!=0) conv_num++;
            }
        }
      size_t sz = sizeof(gyroelectric_susceptibility_data) + sizeof(complex<double>)*(3*gv.ntot())+sizeof(realnum) * (4*num + 1) +sizeof(complex<double>)*(2*conv_num* gv.ntot());   //chi0's       
              //+fields
              //+convolutions
      gyroelectric_susceptibility_data *d = (gyroelectric_susceptibility_data *) malloc(sz);
      d->sz_data = sz;
      return (void*) d;
    };
    
    void gyroelectric_susceptibility::init_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2],
			  double dt, const grid_volume &gv, void *data) const {
        gyroelectric_susceptibility_data *d = (gyroelectric_susceptibility_data *) data;
        complex<double> j(0,1);
        size_t sz_data = d->sz_data;
        memset(d, 0, sz_data);
        d->sz_data = sz_data;
        int ntot = d->ntot = gv.ntot();
        
        /*here we assign field components for our projections. 
        * if we have an electric field polarization Pr_diag_cart[x][x] goes to d->Pr_diag[Ex][Ex], etc
        */ 
        bool elect=0;
        bool magn=0;
        component comp;
        component c1;
        component c2;
        FOR_COMPONENTS(c) DOCMP2 if (needs_P(c, cmp, W)) {
            if (is_electric(c)) elect=1;
            if (is_magnetic(c)) magn=1;
            
        }
        if (elect and magn) abort("We can not have electric and magnetic components in the same polarization");
        if (!elect and !magn) abort("We can not have not electric nor magnetic components in the same polarization");
        if (elect){
            comp=Ex;
        }else{
            comp=Hx;
        }
        for (int n=0; n<3; n++){
            for (int m=0; m<3; m++ ){
                c1=direction_component(comp,(direction)n);
                c2=direction_component(comp,(direction)m);
                d->Pr_diag[c1][c2]=Pr_diag_cart[n][m]; //here I should pray that prof Steven Jonson will not change enum direction in future versions. Fingers crossed.  
                d->Pr_offdiag[c1][c2]=Pr_offdiag_cart[n][m];
                d->Pr_z[c1][c2]=Pr_z_cart[n][m];
            }
        }
        
        realnum *P = d->data;
        realnum *P_prev = d->data + ntot;
        realnum *D_prev = d->data + 2*ntot;
        
       
          
        FOR_COMPONENTS(c) DOCMP2 if (needs_P(c, cmp, W)) {
          d->P[c][cmp] = P;
          d->P_prev[c][cmp] = P_prev;
          //another field type component
          field_type ft2;
          if (field_type(c)==E_stuff){
               ft2=D_stuff;
          }
          component dc = field_type_component(ft2, c);
          //end
          d->D_prev[dc][cmp] = D_prev;
          LOOP_OVER_VOL_OWNED(gv, dc, i) {
              D_prev[i]=0;
          }

          P += 3*ntot;
          P_prev += 3*ntot;
          D_prev += 3*ntot;
        }
        complex<double> *chi0_diag1= reinterpret_cast < complex<double>* >(D_prev += 1);
        complex<double> *chi0_offdiag1=chi0_diag1 + 1*ntot;
        complex<double> *chi0_diag2=chi0_diag1 + 2*ntot;
        complex<double> *chi0_offdiag2=chi0_diag1 + 3*ntot;
        complex<double> *chi0_zz1=chi0_diag1 + 4*ntot;
        complex<double> *chi0_zz2=chi0_diag1 + 5*ntot;
        d->chi0_diag1=chi0_diag1;
        d->chi0_diag2=chi0_diag2;
        d->chi0_offdiag1=chi0_offdiag1;
        d->chi0_offdiag2=chi0_offdiag2;
        d->chi0_zz1=chi0_zz1;
        d->chi0_zz2=chi0_zz2;
        //filling chi0's
        const realnum *s = sigma[comp][0];//sigma. Here we have 0 as some kind of magic number, but physically there is no way to have plasma frequency anisothropic 
            if (s){      
                for (int i = 0; i < ntot; ++i){
                                  
                        chi0_diag1[i]=dt*pow(s[i],2)*(nu_col+j*nu_cycl)/(pow(nu_col,2)+pow(nu_cycl,2));
                        chi0_offdiag1[i]=dt*pow(s[i],2)*(nu_cycl-j*nu_col)/(pow(nu_col,2)+pow(nu_cycl,2));
                        chi0_diag2[i]=(pow(s[i],2)*(nu_col+j*nu_cycl)/(pow(nu_col,2)+pow(nu_cycl,2)))*(1.0-exp(dt*(j*nu_cycl-nu_col)))/(j*nu_cycl-nu_col);
                        chi0_offdiag2[i]=(pow(s[i],2)*(nu_cycl-j*nu_col)/(pow(nu_col,2)+pow(nu_cycl,2)))*(1.0-exp(dt*(j*nu_cycl-nu_col)))/(j*nu_cycl-nu_col);
                        chi0_zz1[i]=dt*pow(s[i],2)/nu_col;
                        chi0_zz2[i]=(pow(s[i],2)/(pow(nu_col,2))*(exp(-1.0*nu_col*dt)-1.0));
                }
            }
        //end of filling chi0's
        complex<double> *dynamic_pointer = chi0_diag1 + 6*ntot; //here we use the dynamic_pointer since some steps might be omitted as the according projection is 0. 
        FOR_COMPONENTS(c) DOCMP2 if (needs_P(c, cmp, W)) {
            FOR_COMPONENTS(c1){
                if (d->Pr_diag[c][c1]!=0){
                    d->conv_diag1[c][c1][cmp]=dynamic_pointer;
                    dynamic_pointer+=ntot;
                    d->conv_diag2[c][c1][cmp]=dynamic_pointer;
                    dynamic_pointer+=ntot;
                    for (int i = 0; i < ntot; ++i){
                        d->conv_diag1[c][c1][cmp][i]=0;
                        d->conv_diag2[c][c1][cmp][i]=0;
                    }
                }
                if (d->Pr_offdiag[c][c1]!=0){
                    d->conv_offdiag1[c][c1][cmp]=dynamic_pointer;
                    dynamic_pointer+=ntot;
                    d->conv_offdiag2[c][c1][cmp]=dynamic_pointer;
                    dynamic_pointer+=ntot;
                    for (int i = 0; i < ntot; ++i){
                        d->conv_offdiag1[c][c1][cmp][i]=0;
                        d->conv_offdiag2[c][c1][cmp][i]=0;
                    }
                }
                if (d->Pr_z[c][c1]!=0){
                    d->conv_z1[c][c1][cmp]=dynamic_pointer;
                    dynamic_pointer+=ntot;
                    d->conv_z2[c][c1][cmp]=dynamic_pointer;
                    dynamic_pointer+=ntot;
                    for (int i = 0; i < ntot; ++i){
                        d->conv_z1[c][c1][cmp][i]=0;
                        d->conv_z2[c][c1][cmp][i]=0;
                    }
                }
            } 
                       
        }
       
        d->exp_decay=exp((j*nu_cycl-nu_col)*dt);
        d->exp_decay_zz=exp(-1.0*nu_col*dt);
  
    };
    

    
    void *gyroelectric_susceptibility::copy_internal_data(void *data) const { 
        gyroelectric_susceptibility_data *d = (gyroelectric_susceptibility_data *) data;
        if (!d) return 0;
        gyroelectric_susceptibility_data *dnew = (gyroelectric_susceptibility_data *) malloc(d->sz_data);
        memcpy(dnew, d, d->sz_data);
        int ntot = d->ntot;
        FOR_COMPONENTS(c){
            FOR_COMPONENTS(c1){
                dnew->Pr_diag[c][c1]=d->Pr_diag[c][c1];
                dnew->Pr_offdiag[c][c1]=d->Pr_offdiag[c][c1];
                dnew->Pr_z[c][c1]=d->Pr_z[c][c1];
            }
        }
        realnum *P = dnew->data;
        realnum *P_prev = dnew->data + ntot;
        realnum *D_prev = dnew->data + 2*ntot;
        
        FOR_COMPONENTS(c) DOCMP2 if (d->P[c][cmp]) {
          dnew->P[c][cmp] = P;
          dnew->P_prev[c][cmp] = P_prev;
          dnew->D_prev[c][cmp] = D_prev;
          
          P += 3*ntot;
          P_prev += 3*ntot;
          D_prev += 3*ntot;
          
        }
        //chi0 stuff here
        complex<double> *chi0_diag1= reinterpret_cast < complex<double>* >(D_prev += 1);
        complex<double> *chi0_offdiag1=chi0_diag1 + 1*ntot;
        complex<double> *chi0_diag2=chi0_diag1 + 2*ntot;
        complex<double> *chi0_offdiag2=chi0_diag1 + 3*ntot;
        complex<double> *chi0_zz1=chi0_diag1 + 4*ntot;
        complex<double> *chi0_zz2=chi0_diag1 + 5*ntot;
        dnew->chi0_diag1=chi0_diag1;
        dnew->chi0_diag2=chi0_diag2;
        dnew->chi0_offdiag1=chi0_offdiag1;
        dnew->chi0_offdiag2=chi0_offdiag2;
        dnew->chi0_zz1=chi0_zz1;
        dnew->chi0_zz2=chi0_zz2;
        //end of chi0 stuff
        complex<double> *dynamic_pointer = chi0_diag1 + 6*ntot;
     FOR_COMPONENTS(c) DOCMP2 if (d->P[c][cmp]) {
            FOR_COMPONENTS(c1){
                if (d->Pr_diag[c][c1]!=0){
                    d->conv_diag1[c][c1][cmp]=dynamic_pointer;
                    dynamic_pointer+=ntot;
                    d->conv_diag2[c][c1][cmp]=dynamic_pointer;
                    dynamic_pointer+=ntot;
                    
                }
                if (d->Pr_offdiag[c][c1]!=0){
                    d->conv_offdiag1[c][c1][cmp]=dynamic_pointer;
                    dynamic_pointer+=ntot;
                    d->conv_offdiag2[c][c1][cmp]=dynamic_pointer;
                    dynamic_pointer+=ntot;
                    
                }
                if (d->Pr_z[c][c1]!=0){
                    d->conv_z1[c][c1][cmp]=dynamic_pointer;
                    dynamic_pointer+=ntot;
                    d->conv_z2[c][c1][cmp]=dynamic_pointer;
                    dynamic_pointer+=ntot;
                    
                }
            } 
        
            
        }
        return (void*) dnew;
     };

    int gyroelectric_susceptibility::num_cinternal_notowned_needed(component c,
					    void *P_internal_data) const {
     gyroelectric_susceptibility_data *d = (gyroelectric_susceptibility_data *) P_internal_data;
     return d->P[c][0] ? 1 : 0;
    };
    realnum *gyroelectric_susceptibility::cinternal_notowned_ptr(int inotowned, component c, int cmp, 
					  int n, 
					  void *P_internal_data) const {
        gyroelectric_susceptibility_data *d = (gyroelectric_susceptibility_data *) P_internal_data;
        (void) inotowned; // always = 0
        if (!d || !d->P[c][cmp]) return NULL;
        return d->P[c][cmp] + n;
    };


}
