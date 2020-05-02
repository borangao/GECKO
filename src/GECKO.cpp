#include <fstream>
#include <RcppArmadillo.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;
vec WLS(mat x, vec y, vec w);
cube LDSC_Initial(int n1,int n2, int num_snp, mat Zscore, mat LDscore, mat Annotation_label,bool Fix_Ve, bool Weight);
void Max_Search_Exhaustive_Grid( double (*Log_L)(double,int,mat,mat,mat,mat,mat,bool),
                                double* a, double* fa,double* b, double* fb, int num_subintervals, double tolerance,int elem_index,mat Vg,mat Ve,mat num_weight_mat,mat ldscore,mat zscore,bool weight);
static int Stopping_Rule(double x0, double x1, double tolerance);
double Log_L(double grid_value, int elem_index,mat Vg,mat Ve,mat num_weight_mat,mat ldscore,mat zscore,bool weight);
double Grid_Search(int ain, int bin,int elem_indexin,int n1in,int n2in,mat LDscorein,mat Zscorein,bool Weightin,mat Vgin,mat Vein);
arma::cube EM(int n1,int n2,int num_snp,arma::cube Vg,arma::mat Ve,arma::mat ldscore,arma::mat zscore,arma::rowvec snp_num_ann,bool weight,bool FixVe);
arma::cube Newton_Raphson(int n1,int n2,int num_snp,arma::cube Vg,arma::mat Ve,arma::mat ldscore,arma::mat zscore,arma::rowvec snp_num_ann,bool weight,bool FixVe);
// [[Rcpp::depends(RcppArmadillo)]]

//' A matrix-based RcppArmadillo implementation for computing x'Ay //'
//' @param x A size n vector
//' @param y A size m vector
//' @param A A size (n x m) matrix
//' @return A scalar value evaluating x'Ay 
//' @export
// [[Rcpp::export]]
arma::cube GECKO(SEXP n1in,SEXP n2in,SEXP nsin,SEXP LDscorein,SEXP Zscorein,SEXP Annotation_labelin,SEXP Weightin,SEXP Fix_Vein,SEXP Test){
    
    const int n1 = Rcpp::as<int>(n1in);
    
    const int n2 = Rcpp::as<int>(n2in);
    
    const int ns = Rcpp::as<int>(nsin);
    
    const mat num_weight_mat ={{n1,sqrt(n1)*sqrt(n2)},{sqrt(n1)*sqrt(n2),n2}};
    
    mat LDscore = as<mat>(LDscorein);
    
    mat Zscore = as<mat>(Zscorein);
    
    mat Annotation_label = as<mat>(Annotation_labelin);
    
    const int num_snp = LDscore.n_rows;
    
    const int num_ann = Annotation_label.n_cols;
    
    rowvec snp_num_ann = sum(Annotation_label,0);
    
    bool weight = as<bool>(Weightin);
    
    bool Fix_ve = as<bool>(Fix_Vein);
    
    bool Test_Indicator = as<bool>(Test);
    
    int max_iter;
    if(Test_Indicator == TRUE){
        max_iter = 0;
    }else{
        max_iter =200;
    }
    // cout<<"Max iter"<<max_iter;
    cube Starting_Point = LDSC_Initial(n1, n2, num_snp, Zscore, LDscore, Annotation_label,Fix_ve, weight);
    
 //   cout<<"Starting_Point is "<<Starting_Point;
    mat grid_search_vg = Starting_Point.slice(1);
    mat grid_search_ve = Starting_Point.slice(0);
    
    cube search_vg_tmp(2,2,5);
    cube search_ve_tmp(2,2,5);
    
    vec LLK(5);
    if(Test_Indicator == FALSE){
        
        
        
        if(num_ann == 1){
            if(Fix_ve == TRUE){
                
                for(int i=0;i<5;i++){
                    grid_search_vg(0,0) = Grid_Search(0, 1,1,n1,n2, LDscore,Zscore,weight,grid_search_vg,grid_search_ve);
                    grid_search_vg(1,1) = Grid_Search(0, 1,2,n1,n2, LDscore,Zscore,weight,grid_search_vg,grid_search_ve);
                    double Vg_01 = Grid_Search(-1,1,3,n1,n2, LDscore,Zscore,weight,grid_search_vg,grid_search_ve);
                    //   cout<<grid_search_ve;
                    grid_search_vg(0,1) = Vg_01;
                    grid_search_vg(1,0) = Vg_01;
                    LLK(i) = Log_L(0,5,grid_search_vg,grid_search_ve,num_weight_mat,LDscore,Zscore,weight);
                    search_vg_tmp.slice(i) = grid_search_vg;
                }
                
            }else{
                
                for(int i=0;i<5;i++){
                    grid_search_vg(0,0) = Grid_Search(0, 1,1,n1,n2, LDscore,Zscore,weight,grid_search_vg,grid_search_ve);
                    grid_search_vg(1,1) = Grid_Search(0, 1,2,n1,n2, LDscore,Zscore,weight,grid_search_vg,grid_search_ve);
                    double Vg_01 = Grid_Search(-1,1,3,n1,n2, LDscore,Zscore,weight,grid_search_vg,grid_search_ve);
                    grid_search_vg(0,1) = Vg_01;
                    grid_search_vg(1,0) = Vg_01;
                    
                    double Ve_01 = Grid_Search(-1,1,4,n1,n2, LDscore,Zscore,weight,grid_search_vg,grid_search_ve);
                    grid_search_ve(0,1) = Ve_01;
                    grid_search_ve(1,0) = Ve_01;
                    LLK(i) = Log_L(0,5,grid_search_vg,grid_search_ve,num_weight_mat,LDscore,Zscore,weight);
                    search_vg_tmp.slice(i) = grid_search_vg;
                    search_ve_tmp.slice(i) = grid_search_ve;
                    
                }
                
                
                
            }
            
            uword i = index_max(LLK);
            grid_search_vg = search_vg_tmp.slice(i);
            grid_search_ve = search_ve_tmp.slice(i);
        }else{
            
            
        }
    }else{
        
        
    }
    cube Vg_initial(2,2,num_ann);
    mat Ve_initial(2,2);
    
    if(Test_Indicator == FALSE){
        if(num_ann == 1){
            Vg_initial.slice(0) = grid_search_vg;
        }else{
            Vg_initial = Starting_Point.slices(1,num_ann);
        }
        Ve_initial = grid_search_ve;
        
    }else{
        if(num_ann == 1){
            Vg_initial.slice(0) = grid_search_vg;
        }else{
            Vg_initial = Starting_Point.slices(1,num_ann);
        }
        Ve_initial = grid_search_ve;
        
    }
    //  cout<<"Ve_initial "<<Ve_initial;
    // Grid_Search(ain, bin,elem_indexin,n1in,n2in, LDscorein,Zscorein,Weightin,Vgin,Vein)
    // cout<<"Log Likelihood is "<<LLK;
    // cube result_output = join_slices(grid_search_vg,grid_search_ve);
    
    double LogL_Old = 0, LogL_New = 0, em_prec = pow(10,-6);
    int num_iter = 0;
    for(int t = 0;t < max_iter; t++){
        
        
        //     cout<<"Old LogLikelihood is "<<LogL_Old<<endl;
        cube out = EM(n1,n2,num_snp,Vg_initial,Ve_initial,LDscore,Zscore,snp_num_ann,weight,Fix_ve);
        /*   mat Ve_inter(2,2);
         cube Vg_inter(2,2,num_ann);
         Ve_inter = out.slice(0);
         Ve_inter(0,0) = 1;
         Ve_inter(1,1) = 1;
         for(int i=0;i<num_ann;i++){
         Vg_inter.slice(i) = out.slice(i+1);
         }
         */
        LogL_New = out(0,0,num_ann+1);
        //      cout<<"New LogLikelihood is "<<LogL_New<<endl;
        if(t!=0 && abs(LogL_New - LogL_Old)<em_prec){
            break;
        }
        LogL_Old = LogL_New;
        Ve_initial = out.slice(0);
        Ve_initial(0,0) =1;
         Ve_initial(1,1) =1;
        for(int i=0;i<num_ann;i++){
            Vg_initial.slice(i) = out.slice(i+1); //First Slice is Likelihood, Second Slice is Ve.
        }
        
 //       cout <<"Vg is "<<Vg_initial<<endl;
//	cout<<"Ve is "<<Ve_initial<<endl;
       /*      if(LogL_New > LogL_Old){
         Ve_initial = out.slice(0);
         Ve_initial(0,0) =1;
         Ve_initial(1,1) =1;
         for(int i=0;i<num_ann;i++){
         Vg_initial.slice(i) = out.slice(i+1); //First Slice is Likelihood, Second Slice is Ve.
         }
         }
         */
        //    cout<<"Vg Cov Est is "<<Vg_initial(0,1,0)<<"Vg H2 Est is "<<Vg_initial(0,0,0)<<" Ve Cov Est is "<<Ve_initial(0,1);
        num_iter = num_iter + 1;
        // cout<<"The iteration is "<<t<<" Vg 0 0 estimate is "<< Vg_initial(0,0,0)<<" Vg 1 1 estimate is "<< Vg_initial(1,1,0)<<" Vg 0 1 estimate is "<< Vg_initial(0,1,0)<<"The Likelihood is "<<LogL_Old<<endl;
    }
    // cout <<"EM iteration is "<<num_iter<<endl;
    // cout<<"Vg Cov Est is "<<Vg_initial(0,1,0)<<"Vg H2 Est is "<<Vg_initial(0,0,0)<<" Ve Cov Est is "<<Ve_initial(0,1);
    
    /* LogL_Old = EM_Initial(n1,n2,num_snp,Vg_initial,Ve_initial,LDscore,Zscore,snp_num_ann,weight);
     LogL_New = LogL_Old + 0.00002;
     num_iter = 0;*/
    
    LogL_New = 0;
    LogL_Old = 0;
    num_iter = 0;
    for(int t = 0;t < max_iter; t++){
        
        //  LogL_Old = LogL_New;
        //  cout<<"Old LogLikelihood is "<<LogL_Old<<endl;
        cube out = Newton_Raphson(n1,n2,num_snp,Vg_initial,Ve_initial,LDscore,Zscore,snp_num_ann,weight,Fix_ve);
        /*    mat Ve_inter(2,2);
         cube Vg_inter(2,2,num_ann);
         Ve_inter = out.slice(0);
         Ve_inter(0,0) = 1;
         Ve_inter(1,1) = 1;
         for(int i=0;i<num_ann;i++){
         Vg_inter.slice(i) = out.slice(i+1);
         }
         */
        LogL_New = out(0,0,num_ann+1);
        
        if(t!=0 && LogL_New - LogL_Old<em_prec){
            break;
        }
        
        LogL_Old = LogL_New;
        Ve_initial = out.slice(0);


        for(int i=0;i<num_ann;i++){
            Vg_initial.slice(i) = out.slice(i+1); //First Slice is Likelihood, Second Slice is Ve.
        }
        num_iter = num_iter + 1;
	// cout<<"NR Vg is "<<Vg_initial<<endl;
	// cout<<"NR Ve is "<<Ve_initial<<endl;
        //    cout<<"The iteration is "<<t<<" Vg 0 0 estimate is "<< Vg_initial(0,0,0)<<" Vg 1 1 estimate is "<< Vg_initial(1,1,0)<<" Vg 0 1 estimate is "<< Vg_initial(0,1,0)<<endl;
    }
    //  cout <<"NR iteration is "<<num_iter<<endl;
    cube out(2,2,num_ann+1);
    out.slice(0) = Ve_initial;
    out(0,1,0) =  out(0,1,0) /ns *sqrt(n1)*sqrt(n2);
    out(1,0,0) =  out(1,0,0) /ns *sqrt(n1)*sqrt(n2);
    if(Fix_ve == FALSE){
        for(int i=0;i<num_ann;i++){
            out.slice(i+1) = Vg_initial.slice(i);
            out.slice(0) = out.slice(0) - Vg_initial.slice(i);
        }
    }else{
        for(int i=0;i<num_ann;i++){
            out.slice(i+1) = Vg_initial.slice(i);
            out.slice(0) = out.slice(0) - Vg_initial.slice(i);
        }
        out(0,1,0) = 0;
        out(1,0,0) = 0;
        
    }
    
    //cube out = NR_ANN(n1,n2,num_snp,Vg_initial,Ve_initial,LDscore,Zscore,snp_num_ann,weight);
    return(out);
    
}
vec WLS(mat x, vec y, vec w){
    vec intercept(x.n_rows);
    intercept.ones();
    mat X_intercept = join_rows(intercept,x);
    //  cout<<"X_intercept col "<<X_intercept.n_cols<<"X_intercept row "<<X_intercept.n_rows<<"w dimesion is "<<w.n_elem;
    X_intercept.each_col()%= sqrt(w);
    vec y_weighted = y % sqrt(w);
    vec result = inv(X_intercept.t() * X_intercept) * X_intercept.t() * y_weighted;
    return(result);
}

cube LDSC_Initial(int n1,int n2, int num_snp, mat Zscore, mat LDscore, mat Annotation_label,bool Fix_Ve, bool Weight){
    vec Z1_square = square(Zscore.col(0));
    
    vec Z2_square = square(Zscore.col(1));
    
    vec LD_all_data(num_snp);
    
    
    int num_ann = size(Annotation_label)(1);
    
    rowvec snp_num_ann = sum(Annotation_label,0);
    
    LD_all_data.zeros();
    for(int i=0;i < num_ann;i++){
        LD_all_data += LDscore.col(i);
        
    }
    
    // Weighted Linear Regression for Heritability and Genetic Covariance
    double tau1 = mean(Z1_square-1)/(n1*mean(LD_all_data));
    double tau2 = mean(Z2_square-1)/(n2*mean(LD_all_data));
    
    vec w1 = 1 /(LD_all_data % (1 + n1 * tau1 * LD_all_data) % (1 + n1 * tau1 * LD_all_data));
    vec w2 = 1 /(LD_all_data % (1 + n2 * tau2 * LD_all_data) % (1 + n2 * tau2 * LD_all_data));
    
    // vec trait_1_lin = WLS (LDscore.col(0), Z1_square, w1);
    // vec trait_2_lin = WLS (LDscore.col(0), Z2_square, w2);
    
    
    vec trait_1_lin = WLS (LDscore, Z1_square, w1);
    vec trait_2_lin = WLS (LDscore, Z2_square, w2);
    
  //  cout << " trait_1_lin is "<<trait_1_lin;
 //   cout << " trait_2_lin is "<<trait_2_lin;
    //Weight for Genetic Covariance
    vec weight_1 = 1+(mean(Z1_square)-1)/mean(LD_all_data)*LD_all_data/num_snp;
    vec weight_2 = 1+(mean(Z2_square)-1)/mean(LD_all_data)*LD_all_data/num_snp;
    vec weight_3 = mean(Zscore.col(0)%Zscore.col(1))*LD_all_data;
    vec w = 1/(weight_1 % weight_2 + weight_3 % weight_3);
    
    vec trait_1_2_lin = WLS (LDscore, Zscore.col(0) % Zscore.col(1), w);
 //   cout << " trait_1_2_lin is "<<trait_1_2_lin;
    double corr_pheno;
    if(Fix_Ve == TRUE){
        corr_pheno = 0;
    }else{
        corr_pheno = trait_1_2_lin(0);
    }
    //  cout<<"Ve off diagonal is "<<corr_pheno<<endl;
    /*Initial For Covariance Estimate of Vg*/
    mat M_mat(num_ann,num_ann);
    vec V_vec(num_ann);
    for(int i=0;i<num_ann;i++) {
        V_vec(i)=sum(Zscore.col(0) % Annotation_label.col(i) % Zscore.col(1) % Annotation_label.col(i))/snp_num_ann(i)/sqrt(n1)/sqrt(n2);
        
        
        for(int j=0;j<num_ann;j++){
            M_mat(i,j)=sum(LDscore.col(j)%Annotation_label.col(i))/snp_num_ann(i)/snp_num_ann(j);
        }
    }
    
    vec rho_hat = inv(M_mat) * (V_vec-corr_pheno/sqrt(n1)/sqrt(n2));
    // cout<<"rho_hat is "<<rho_hat;
    cube Vg_search(2,2,num_ann);
    for(int i = 0; i < num_ann; i++){
        Vg_search.slice(i) = {{trait_1_lin(i+1) * snp_num_ann(i)/n1,rho_hat(i)},{rho_hat(i),trait_2_lin(i+1) * snp_num_ann(i)/n2}};
    }
    //  cout<<"Vg_search is "<<Vg_search;
    mat Ve_search(2,2);
    Ve_search.ones();
    Ve_search(0,1) = corr_pheno;
    Ve_search(1,0) = corr_pheno;
    
    cube LDSC_output = join_slices(Ve_search,Vg_search);
    return(LDSC_output);
    
}

double Grid_Search(int ain, int bin,int elem_indexin,int n1in,int n2in,mat LDscorein,mat Zscorein,bool Weightin,mat Vgin,mat Vein){
    int n1 = n1in;
    
    int n2 = n2in;
    
    const mat num_weight_mat ={{n1,sqrt(n1)*sqrt(n2)},{sqrt(n1)*sqrt(n2),n2}};
    
    mat LDscore = LDscorein;
    
    mat Zscore = Zscorein;
    
    mat Vg = Vgin;
    mat Ve = Vein;
    
    const int num_snp = LDscore.n_rows;
    
    
    bool weight = Weightin;
    
    double ain_c = ain;
    double bin_c = bin;
    int elem_index = elem_indexin;
    
    double fa, fb;
    
    double tolerance = 0.01;
    int    num_subintervals = 10;
    Max_Search_Exhaustive_Grid( Log_L, &ain_c, &fa, &bin_c, &fb, num_subintervals,tolerance,elem_index,Vg,Ve,num_weight_mat,LDscore,Zscore,weight);
    // Rcout<<ain_c<<" "<<bin_c;
    return((ain_c+bin_c)/2);
    
    
}


//     double a = -1.0;                                                       //


static int Stopping_Rule(double x0, double x1, double tolerance)
{
    double xm = 0.5 * fabs( x1 + x0 );
    
    if ( xm <= 1.0 ) return ( fabs( x1 - x0 ) < tolerance ) ? 1 : 0;
    return ( fabs( x1 - x0 ) < tolerance * xm ) ? 1 : 0;
}

double Log_L(double grid_value, int elem_index,mat Vg,mat Ve,mat num_weight_mat,mat ldscore,mat zscore,bool weight) {
    
    if(elem_index == 1){
        Vg(0,0) = grid_value;
    }else if(elem_index == 2){
        Vg(1,1) =  grid_value;
    }else if(elem_index == 3){
        Vg(0,1) = grid_value;
        Vg(1,0) = grid_value;
    }else if(elem_index == 4){
        Ve(0,1) = grid_value;
        Ve(1,0) = grid_value;
    }else{
        
    }
    Vg = Vg % num_weight_mat;
    vec ld_weight(ldscore.n_rows);
    int num_snp = ldscore.n_rows;
    if(weight == TRUE){
        ld_weight = 1/ldscore.col(0);
    }else{
        ld_weight.ones();
    }
    
    double LogL = 0;
    for(arma::uword i=0;i<ldscore.n_rows;i++){
        arma::rowvec zj = zscore.row(i);
        arma::mat sigma = 1/ld_weight(i)/num_snp * Vg + Ve;
        arma::mat sigma_inv = inv(sigma);
        LogL-= 0.5 * ld_weight(i)*(trace(zj*sigma_inv*zj.t()) + log(det(sigma)));
    }
    return(LogL);
}


void Max_Search_Exhaustive_Grid( double (*Log_L)(double,int,mat,mat,mat,mat,mat,bool),
                                double* a, double* fa,double* b, double* fb, int num_subintervals, double tolerance,int elem_index,mat Vg,mat Ve,mat num_weight_mat,mat ldscore,mat zscore,bool weight)
{
    double h;
    double x;
    double fx;
    double old_f;
    double xmax;
    double fmax;
    int i;
    int save_next;
    
    // Verify that the tolerance is an acceptable number
    
    // if (tolerance <= 0.0) tolerance = sqrt(DBL_EPSILON) * (b - a);
    
    // Verify that the number of subintervals is at least 3.
    
    // if (num_subintervals < 3) num_subintervals = 3;
    
    // Iterate until the length of the interval is smaller than the tolerance.
    
    while ( ! Stopping_Rule( *a, *b, tolerance) ) {
        h = (*b - *a) / (double) num_subintervals;
        x = *a + h;
        xmax = x;
        fmax = Log_L(x,elem_index,Vg,Ve,num_weight_mat,ldscore,zscore,weight);
        
        save_next = 1;
        for (i = 2; i < num_subintervals; i++) {
            old_f = fx;
            x += h;
            fx = Log_L(x,elem_index,Vg,Ve,num_weight_mat,ldscore,zscore,weight);
            
            //    cout<<i<<" iteration "<<x<<" x value "<<fx;
            if (fx > fmax) { xmax = x; fmax = fx; *fa = old_f; save_next = 1;}
            else if (save_next) {*fb = fx; save_next = 0; }
        }
        *a = xmax - h;
        *b = xmax + h;
    }
    return;
    
}

arma::cube EM(int n1,int n2,int num_snp,arma::cube Vg,arma::mat Ve,arma::mat ldscore,arma::mat zscore,arma::rowvec snp_num_ann,bool weight,bool FixVe) {
    
    
    
    arma::rowvec annotation = ldscore.row(1);
    int num_ann = annotation.size();
    
    arma::vec ld_weight(num_snp);
    ld_weight.zeros();
    if(weight == TRUE){
        for(int i=0;i<num_ann;i++){
            ld_weight+=ldscore.col(i);
        }
        ld_weight = 1/ld_weight;
    }else{
        ld_weight+=1;
    }
    
    arma::mat num_weight_mat(2,2);
    num_weight_mat(0,0) = n1;
    num_weight_mat(1,0) = sqrt(n1)*sqrt(n2);
    num_weight_mat(0,1) = sqrt(n1)*sqrt(n2);
    num_weight_mat(1,1) = n2;
    for(int i=0;i<num_ann;i++){
        Vg.slice(i) = Vg.slice(i) % num_weight_mat;
    }
    
    
    arma::cube par_cube(2,2,num_ann+2);
    par_cube.zeros();
    
    
    if(FixVe == TRUE){
        Ve(0,1) = 0;
        Ve(1,0) = 0;
        Ve(0,0) = 1;
        Ve(1,1) = 1;
        
    }else{
        Ve(0,0) = 1;
        Ve(1,1) = 1;
    }
    
    for(int i=0;i<num_snp;i++){
        
        double snp_weight = ld_weight(i);
        arma::rowvec lj = ldscore.row(i);
        
        arma::rowvec rj = lj/snp_num_ann;
        
        arma::rowvec zj = zscore.row(i);
        arma::mat sigma(2,2);
        sigma.zeros();
        for(int i=0;i<num_ann;i++){
            sigma += rj(i) * Vg.slice(i);
        }
        sigma = sigma + Ve;
        
        arma::mat sigma_inv = inv(sigma);
        
        /************************
         
         Estimate of E_ej_ej
         
         
         ************************/
        
        arma::mat ve_sigma_inv = Ve*sigma_inv;
        arma::colvec E_ej_zj = ve_sigma_inv*zj.t();
        arma::mat V_ej_zj = Ve-ve_sigma_inv*Ve;
        arma::mat E_ej_ej = E_ej_zj*E_ej_zj.t()+V_ej_zj;
        par_cube.slice(0) += E_ej_ej*snp_weight;
        /************************
         
         Estimate of E_uj_uj
         
         
         ************************/
        arma::mat LogL(1,1);
        LogL.zeros();
        
        for(int i=0;i<num_ann;i++){
            arma::mat V_kg_inv = rj(i)*Vg.slice(i)*sigma_inv;
            arma::colvec E_gj_zj = V_kg_inv*zj.t();
            arma::mat V_gj_zj = rj(i)*Vg.slice(i)-rj(i)*V_kg_inv*Vg.slice(i);
            arma::mat E_gj_gj = (E_gj_zj*E_gj_zj.t()+V_gj_zj)*snp_weight;
            
            
            if(rj(i)==0){
                par_cube.slice(i+1) += 0;
                LogL +=0;
            }else{
                par_cube.slice(i+1) += 1/rj(i)*E_gj_gj;
                //   LogL += (-0.5 * trace(E_gj_gj * inv(rj(i) * Vg.slice(i))) - 0.5 * log(det(rj(i) * Vg.slice(i))))*snp_weight;
            }
            
            
        }
        
        //   LogL = LogL - snp_weight*0.5*trace(E_ej_ej*inv(Ve)) - snp_weight*0.5*log(det(Ve));
        LogL = LogL - snp_weight*0.5*trace(zj*sigma_inv*zj.t()) - snp_weight*0.5*log(det(sigma));
        par_cube(0,0,num_ann+1) += LogL(0,0);
        //Rcout<<par_cube;
    }
    
    arma::cube out(2,2,num_ann+2);
    for(int i=0;i<num_ann+1;i++){
        out.slice(i) = par_cube.slice(i)/sum(ld_weight);
    }
    for(int i=1;i<num_ann+1;i++){
        out.slice(i) =  out.slice(i)/num_weight_mat;
    }
    out.slice(num_ann+1) = par_cube.slice(num_ann+1);
    
    return(out);
}
arma::cube Newton_Raphson(int n1,int n2,int num_snp,arma::cube Vg,arma::mat Ve,arma::mat ldscore,arma::mat zscore,arma::rowvec snp_num_ann,bool weight,bool FixVe) {
    
    
    if(FixVe == TRUE){
        Ve(0,1) = 0;
        Ve(1,0) = 0;
        Ve(0,0) = 1;
        Ve(1,1) = 1;
        
    }else{
        Ve(0,0) = 1;
        Ve(1,1) = 1;
    }
    
    arma::rowvec annotation = ldscore.row(1);
    int num_ann = annotation.size();
    
    arma::vec ld_weight(num_snp);
    ld_weight.zeros();
    if(weight == TRUE){
        for(int i=0;i<num_ann;i++){
            ld_weight+=ldscore.col(i);
        }
        ld_weight = 1/ld_weight;
    }else{
        ld_weight+=1;
    }
    
    arma::mat num_weight_mat(2,2);
    num_weight_mat(0,0) = n1;
    num_weight_mat(1,0) = sqrt(n1)*sqrt(n2);
    num_weight_mat(0,1) = sqrt(n1)*sqrt(n2);
    num_weight_mat(1,1) = n2;
    for(int i=0;i<num_ann;i++){
        Vg.slice(i) = Vg.slice(i) % num_weight_mat;
    }
    
    arma::cube par_space(2,2,num_ann+1);
    par_space.slice(0) = Ve;
    for (int i=0;i<num_ann;i++){
        par_space.slice(i+1) = Vg.slice(i);
    }
    arma::rowvec initial_vec(3*num_ann+3);
    
    for (int k=0;k<num_ann+1;k++){
        arma::mat part_cov = par_space.slice(k);
        for (int i=0;i<2;i++){
            for (int j=0;j<=i;j++){
                initial_vec(3*k+i+j) = part_cov(i,j);
            }
        }
    }
    arma::rowvec first_der_vec(3*num_ann+3);
    first_der_vec.zeros();
    arma::mat sec_der_mat(3*num_ann+3,3*num_ann+3);
    sec_der_mat.zeros();
    double LogL = 0;
    
    for(int i=0;i<num_snp;i++){
        double snp_weight = ld_weight(i);
        arma::rowvec lj = ldscore.row(i);
        arma::rowvec rj = lj/snp_num_ann;
        
        arma::rowvec rj_ve(num_ann+1);
        rj_ve(0) = 1;
        for(int i=0;i<num_ann;i++){
            
            rj_ve(i+1) = rj(i);
            
        }
        //Rcout<<'rj is'<<rj_ve<<arma::endl;
        
        arma::rowvec zj = zscore.row(i);
        arma::mat sigma(2,2);
        sigma.zeros();
        for(int i=0;i<lj.size();i++){
            sigma +=rj(i)*Vg.slice(i);
        }
        sigma = sigma + Ve;
        //Rcout<<sigma(0,0)<<','<<sigma(0,1)<<','<<sigma(1,1);
        
        arma::mat sigma_inv = inv(sigma);
        LogL = LogL - snp_weight*0.5*trace(zj*sigma_inv*zj.t()) - snp_weight*0.5*log(det(sigma));
        /************************
         
         First Order Derivative
         
         
         ************************/
        arma::rowvec first_der_num_cons(3);
        
        arma::mat omega_mat_11(2,2);//Omega11
        omega_mat_11.zeros();
        omega_mat_11(0,0) = 1;
        arma::mat sigma_inv_omega_11 = sigma_inv * omega_mat_11;
        arma::mat sigma_inv_omega_11_inv = zj * sigma_inv_omega_11 * sigma_inv * zj.t();
        first_der_num_cons(0) = trace(-0.5*sigma_inv_omega_11) + 0.5 * sigma_inv_omega_11_inv(0,0);
        
        
        arma::mat omega_mat_12(2,2); //Omega12
        omega_mat_12.zeros();
        omega_mat_12(0,1) = 1;
        omega_mat_12(1,0) = 1;
        arma::mat sigma_inv_omega_12 = sigma_inv * omega_mat_12;
        arma::mat sigma_inv_omega_12_inv = zj * sigma_inv_omega_12 * sigma_inv * zj.t();
        first_der_num_cons(1) = trace(-0.5*sigma_inv_omega_12) + 0.5 * sigma_inv_omega_12_inv(0,0);
        
        arma::mat omega_mat_22(2,2); //Omega22
        omega_mat_22.zeros();
        omega_mat_22(1,1) = 1;
        arma::mat sigma_inv_omega_22 = sigma_inv * omega_mat_22;
        arma::mat sigma_inv_omega_22_inv = zj * sigma_inv_omega_22 * sigma_inv * zj.t();
        first_der_num_cons(2) = trace(-0.5*sigma_inv_omega_22) + 0.5 * sigma_inv_omega_22_inv(0,0);
        
        
        for (int i=0;i<num_ann+1;i++){
            first_der_vec(3*i) += first_der_num_cons(0) * rj_ve(i) * snp_weight;
            first_der_vec(3*i+1) += first_der_num_cons(1) * rj_ve(i)* snp_weight;
            first_der_vec(3*i+2) += first_der_num_cons(2) * rj_ve(i)* snp_weight;
        }
        //Rcout<<first_der_vec;
        /************************
         
         Second Order Derivative
         
         
         ************************/
        
        arma::mat rj_ve_rj_ve = rj_ve.t() * rj_ve;//7*7 matrix
        //Rcout<<rj_ve_rj_ve;
        arma::mat num_cons_mat(3,3); // 3*3 matrix of scalar which will times all 7*7 elements above
        num_cons_mat.zeros();
        arma::cube sigma_inv_omega(2,2,3);
        sigma_inv_omega.slice(0) = sigma_inv_omega_11;
        sigma_inv_omega.slice(1) = sigma_inv_omega_12;
        sigma_inv_omega.slice(2) = sigma_inv_omega_22;
        
        for (int i=0;i<3;i++){
            for (int j=0;j<3;j++){
                arma::mat sigma_inv_num_inv_num =  sigma_inv_omega.slice(i) *  sigma_inv_omega.slice(j);
                arma::mat zj_sigma_inv_num_inv_num_zj =  zj * sigma_inv_num_inv_num *  sigma_inv * zj.t();
                double num_cons = 0.5*trace(sigma_inv_num_inv_num) - zj_sigma_inv_num_inv_num_zj(0,0);
                num_cons_mat(i,j) = num_cons;
                
            }
        }
        //  Rcout<<num_cons_mat<<arma::endl;
        sec_der_mat += kron(rj_ve_rj_ve,num_cons_mat)*snp_weight;
        //Rcout<<sec_der_mat(3,3);
        
    }
    arma::rowvec first_der_vec_new(first_der_vec.size()-2);
    first_der_vec_new(0) = first_der_vec(1);
    for (int i=1;i<first_der_vec_new.size();i++){
        first_der_vec_new(i) = first_der_vec(i+2);
    }
    sec_der_mat.shed_row(0);
    sec_der_mat.shed_row(1);
    sec_der_mat.shed_col(0);
    sec_der_mat.shed_col(1);
    
    
    
    //Rcout<<first_der_vec <<sum(first_der_vec) ;
    arma::rowvec update = first_der_vec_new * inv(sec_der_mat);
    update.insert_cols(1,1);
    update.insert_cols(0,1);
    
    for(int i=0;i<num_ann;i++){
        update(3*i+3) =  update(3*i+3)/n1;
        update(3*i+1+3) =  update(3*i+1+3)/sqrt(n1)/sqrt(n2);
        update(3*i+2+3) =  update(3*i+2+3)/n2;
    }
    
    arma::rowvec out = initial_vec - update;
    arma::cube out_cube(2,2,num_ann+2);
    out_cube(0,0,0) = out(0);
    out_cube(1,0,0) = out(1);
    out_cube(0,1,0) = out(1);
    out_cube(1,1,0) = out(2);
    for(int i=1;i<num_ann+1;i++){
        out_cube(0,0,i) = out(3*i)/n1;
        out_cube(1,0,i) = out(3*i+1)/sqrt(n1)/sqrt(n2);
        out_cube(0,1,i) = out(3*i+1)/sqrt(n1)/sqrt(n2);
        out_cube(1,1,i) = out(3*i+2)/n2;
    }
    out_cube(0,0,num_ann+1) = LogL;
    return(out_cube);
    
}
/*** R
 #setwd("C:/Users/Boran Gao/Desktop/GECKO/code")
 #LDscore<-read.table("wtccc_ctrl_no_annot.csv.gz",header = T)
 #LDscore<-as.matrix(LDscore$L2)
 #Z_1<-read.table("summary_-0.4_1_z1.txt",header = T)
 #Z_2<-read.table("summary_-0.4_1_z2.txt",header = T)
 #ain = 0
 #bin = 1
 #elem_indexin = 1
 #LDscorein = LDscore
 #Zscorein = cbind(Z_1$Z,Z_2$Z)
 #Weightin = T
 #Vgin = matrix(c(0.5,0,0,0.5),ncol = 2,nrow = 2)
 #Vein = matrix(c(0.5,0,0,0.5),ncol = 2,nrow = 2)
 #n1in = 2938
 #n2in = 2938
 #ns1n = 2938
 #Fix_Vein = F
 #Annotation_labelin = matrix(rep(1,length(LDscorein)),ncol=1)
 #Test = F
 #result = GECKO(n1in,n2in, ns1n,LDscorein,Zscorein,Annotation_labelin,Weightin,Fix_Vein,Test )
 */

