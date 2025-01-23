#define ARMA_USE_SUPERLU 1

// [[Rcpp::depends(BH)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>



#include "SAIGE_test.hpp"
#include "SPA.hpp"
#include "ER_binary_func.hpp"
#include "UTIL.hpp"
#include "approxfun.hpp"
#include "getMem.hpp"
#include "getMem.hpp"
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/date_time.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <R.h>
#include <Rmath.h>



namespace SAIGE {

SAIGEClass::SAIGEClass(
        arma::mat & t_XVX,
        arma::mat & t_XXVX_inv,
        arma::mat & t_XV,
        arma::mat & t_XVX_inv_XV,
        arma::mat & t_Sigma_iXXSigma_iX,
        arma::mat & t_X,
        arma::mat &  t_S_a,
        arma::mat & t_res,
        arma::mat & t_mu2,
        arma::mat & t_mu,
        arma::mat & t_varRatio_sparse,
        arma::mat & t_varRatio_null,
        arma::mat & t_varRatio_null_sample,
        arma::mat & t_varRatio_null_noXadj,
        arma::mat & t_varRatio_null_eg,
        arma::mat & t_varRatio_sparse_eg,
        arma::vec & t_cateVarRatioMinMACVecExclude,
        arma::vec & t_cateVarRatioMaxMACVecInclude,
        double t_SPA_Cutoff,
        arma::mat & t_tauvec,
        arma::mat & t_varWeightsvec,
        std::vector<std::string> & t_traitType,
        arma::mat & t_y,
        std::string t_impute_method,
        bool t_flagSparseGRM,
        bool t_isnoadjCov,
        double t_pval_cutoff_for_fastTest,
        bool t_isCondition,
        std::vector<uint32_t> & t_condition_genoIndex,
        bool t_is_Firth_beta,
        double t_pCutoffforFirth,
        arma::mat & t_offset,
        arma::mat & t_resout,
        arma::sp_mat & t_SigmaMat_sp,
        float t_tauVal_sp,
         arma::sp_mat & t_Ilongmat,
        arma::vec & t_I_longl_vec,
        arma::sp_mat & t_Tlongmat,
        arma::vec & t_T_longl_vec,
        bool t_is_EmpSPA,
        arma::mat & t_cumul,
	bool t_is_gxe,
	arma::mat &  t_XV_gxe,
   arma::mat &  t_X_gxe,
        arma::mat &  t_XVX_inv_XV_gxe,
        arma::mat &  t_XVX_gxe,
	        arma::mat & t_S_a_gxe,
        arma::mat & t_XXVX_inv_gxe,
        arma::mat & t_y_gxe,
        arma::mat & t_res_gxe,
        arma::mat & t_mu2_gxe,
        arma::mat & t_mu_gxe,
        arma::mat & t_varWeights_gxe){


m_is_gxe = t_is_gxe;
    m_y_mt = t_y;
    m_XV_mt = t_XV;
    m_XXVX_inv_mt = t_XXVX_inv;
    m_XVX_mt = t_XVX;

    m_XVX_inv_XV_mt = t_XVX_inv_XV;
    m_Sigma_iXXSigma_iX_mt = t_Sigma_iXXSigma_iX;
    m_isVarPsadj = false;
    m_X_mt = t_X;
    m_S_a_mt = t_S_a;
    m_res_mt = t_res;
    m_resout_mt = t_resout;
    m_mu2_mt = t_mu2 % t_varWeightsvec;
    m_mu_mt = t_mu;
    m_varRatio_sparse_mt = t_varRatio_sparse;
    m_varRatio_null_mt = t_varRatio_null;
    m_varRatio_null_sample_mt = t_varRatio_null_sample;
    m_varRatio_null_noXadj_mt = t_varRatio_null_noXadj;
    m_varRatio_null_eg_mt = t_varRatio_null_eg;
    m_varRatio_sparse_eg_mt = t_varRatio_sparse_eg;
    m_tauvec_mt = t_tauvec;
    m_varWeightsvec_mt = t_varWeightsvec;
    
    
    m_number_of_traits = t_Ilongmat.n_rows;
    m_n = m_res_mt.n_rows;
    m_tauvec_each_length = t_tauvec.n_rows;

/*
    if(m_Sigma_iXXSigma_iX.n_cols == 1 && m_Sigma_iXXSigma_iX.n_rows == 1){
	m_isVarPsadj = false;
    }else{
	m_isVarPsadj = true;
    }
*/


    m_cateVarRatioMinMACVecExclude = t_cateVarRatioMinMACVecExclude;
    m_cateVarRatioMaxMACVecInclude = t_cateVarRatioMaxMACVecInclude;
    m_traitType_vec = t_traitType;

    //m_case_indices = arma::find(m_y == 1);
    //m_ctrl_indices = arma::find(m_y == 0);
    m_p = t_X.n_cols;
    m_SPA_Cutoff = t_SPA_Cutoff;
    m_impute_method =  t_impute_method;
    m_isCondition = t_isCondition;
    m_condition_genoIndex = t_condition_genoIndex;
    if(m_isCondition){	
	        m_numMarker_cond = t_condition_genoIndex.size();      
    }else{
		m_numMarker_cond = 0;
    }	    
    //m_p = t_X.nrow();

    if(m_traitType == "binary"){
      //if(m_isOutputAFinCaseCtrl){
      if(m_y_mt.n_cols == 1){
        m_case_indices = arma::find(m_y_mt.col(0) == 1);
        m_ctrl_indices = arma::find(m_y_mt.col(0) == 0);
	m_n_case = m_case_indices.n_elem;
	m_n_ctrl = m_ctrl_indices.n_elem;
      } 	
	m_is_Firth_beta = t_is_Firth_beta;
	m_pCutoffforFirth = t_pCutoffforFirth;
    }
    m_offset_mt = t_offset;
    //m_dimNum = t_dimNum;
    m_flagSparseGRM = t_flagSparseGRM;
    //m_isFastTest = t_isFastTest;
    m_isnoadjCov = t_isnoadjCov;
    m_pval_cutoff_for_fastTest = t_pval_cutoff_for_fastTest;
    //if(m_dimNum != 0){
    //	m_locationMat = t_locationMat;
    //	m_valueVec = t_valueVec;
    //}
    //
    //if(t_SigmaMat_sp.n_rows > 2){
	//std::cout << "here m_SigmaMat_sp" << std::endl;
	m_SigmaMat_sp = t_SigmaMat_sp;
    //}	   
   m_tauVal_sp = t_tauVal_sp;
   g_I_longl_mat = t_Ilongmat;
   arma::uvec t_I_longl_vec_new = arma::conv_to< arma::uvec >::from(t_I_longl_vec);
   g_I_longl_vec = t_I_longl_vec_new;
   g_T_longl_mat = t_Tlongmat;
   g_T_longl_vec = t_T_longl_vec;

   //m_res_sample = (g_I_longl_mat.t()) * m_res;
   //m_mu2_sample = (g_I_longl_mat.t()) * m_mu2;
   //m_mu_sample = (g_I_longl_mat.t()) * m_mu;


   m_is_EmpSPA = false;
   if(t_is_EmpSPA){
      m_K_0_emp.setApproxFun(t_cumul.col(0), t_cumul.col(1));
      m_K_1_emp.setApproxFun(t_cumul.col(0), t_cumul.col(2));
      m_K_2_emp.setApproxFun(t_cumul.col(0), t_cumul.col(3));
      m_is_EmpSPA = t_is_EmpSPA;
      m_cumul = t_cumul;
      m_varResid = arma::var(m_res);
   }
      m_X_gxe_mt = t_X_gxe;
      m_XVX_inv_XV_gxe_mt = t_XVX_inv_XV_gxe;
      m_XVX_gxe_mt = t_XVX_gxe;

      m_S_a_gxe_mt = t_S_a_gxe;
      m_XV_gxe_mt = t_XV_gxe;
      m_XXVX_inv_gxe_mt = t_XXVX_inv_gxe;
      m_y_gxe_mt = t_y_gxe;
      m_res_gxe_mt = t_res_gxe;
      m_mu2_gxe_mt = t_mu2_gxe;
      m_mu_gxe_mt = t_mu_gxe;
      m_varWeights_gxe_mt = t_varWeights_gxe;
}    


// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
void SAIGEClass::set_seed(unsigned int seed){
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

void SAIGEClass::scoreTest(arma::vec & t_GVec,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2,
                     arma::vec & t_gtilde,
		     arma::vec & t_P2Vec,
		     double& t_gy, 
		     bool t_is_region,
		     arma::uvec & t_indexForNonZero, 
		     double & t_pval){
    arma::vec Sm, var2m;
    double S, var2;
    //t_GVec.print("t_GVec");
    //t_gtilde.print("t_gtilde before");
    getadjGFast(t_GVec, t_gtilde, t_indexForNonZero);
    //t_gtilde.print("t_gtilde after");
    //getadjG(t_GVec, t_gtilde);


    if(t_is_region && m_traitType == "binary"){
      t_gy = dot(t_gtilde, m_y_mt.col(m_itrait));
    }
    
    S = dot(t_gtilde, m_res_mt.col(m_itrait) % m_varWeightsvec_mt.col(m_itrait));
     
    //std::cout << "S " << S << std::endl;

    S = S/m_tauvec_mt(0,m_itrait); 

    double varRatioVal_var2 = m_varRatioVal;

   if(!m_flagSparseGRM_cur){
      t_P2Vec = t_gtilde % (m_mu2_mt.col(m_itrait)) *(m_tauvec_mt(0,m_itrait)); 
      //t_P2Vec = t_gtilde % m_mu2; 
      var2m = dot(t_P2Vec, t_gtilde);
    }else{
      if(m_SigmaMat_sp.n_rows > 2){	
      //t_P2Vec = arma::spsolve(m_SigmaMat_sp, t_gtilde);
      //std::cout << "m_startin " << m_startin << " m_endin " << m_endin << std::endl;
      //arma::sp_mat m_SigmaMat_sp_firstcol = m_SigmaMat_sp.cols(m_startin, m_endin);
      //m_SigmaMat_sp_firstcol.print("m_SigmaMat_sp_firstcol");
      //t_gtilde.print("t_gtilde");
      t_P2Vec = (m_SigmaMat_sp.cols(m_startin, m_endin)) * t_gtilde;

      var2m = dot(t_P2Vec, t_gtilde);
      //std::cout << "m_isVarPsadj " << m_isVarPsadj << std::endl;

      if(m_isVarPsadj){
	varRatioVal_var2 = 1;     
	var2m = var2m - t_gtilde.t() * m_Sigma_iXXSigma_iX * m_X.t() * t_P2Vec;	
      }
     }else{
	//t_P2Vec = m_sigmainvG_noV;
      //std::cout << t_gtilde.n_elem <<  std::endl;
      //#Sigma_iGE = getSigma_G_multiV(W, tauVecNew, GE_tilde, maxiterPCG, tolPCG, LOCO=FALSE)
      t_P2Vec = getSigma_G_V(t_gtilde, 500, 1e-5);
      var2m = dot(t_P2Vec, t_gtilde);
      if(m_isVarPsadj){
	varRatioVal_var2 = 1;      
	var2m = var2m - t_gtilde.t() * m_Sigma_iXXSigma_iX * m_X.t() * t_P2Vec;	
      }
     }
    }

    var2 = var2m(0,0);
    //std::cout << "var2 " << var2 << std::endl;
    //std::cout << "m_varRatioVal " << m_varRatioVal << std::endl;
    //double var1 = var2 * m_varRatioVal;
    double var1 = var2 * varRatioVal_var2;
    //std::cout << "var1 " << var1 << std::endl;
    double stat = S*S/var1;
    //double t_pval;
    //std::cout << "S " << S << std::endl;    
    //std::cout << "var1 " << var1 << std::endl;    

    //if (var1 <= std::pow(std::numeric_limits<double>::min(), 2)){
    if (var1 <= std::numeric_limits<double>::min()){
        t_pval = 1;
    }else{
        boost::math::chi_squared chisq_dist(1);
        t_pval = boost::math::cdf(complement(chisq_dist, stat));
    }

    char pValueBuf[100];
    if (t_pval != 0)
        sprintf(pValueBuf, "%.6E", t_pval);
    else {
        double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
}


void SAIGEClass::scoreTestFast(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2, 
		     double & t_pval){

    arma::vec g1 = t_GVec.elem(t_indexForNonZero);

    arma::vec m_varWeightsvec_new = m_varWeightsvec.elem(t_indexForNonZero);

    arma::mat X1 = m_X.rows(t_indexForNonZero);
    arma::mat A1 = m_XVX_inv_XV.rows(t_indexForNonZero);
    arma::vec mu21;
    arma::vec res1 = m_res.elem(t_indexForNonZero);
    res1 = res1 % m_varWeightsvec_new;
    arma::vec Z = A1.t() * g1;
    arma::vec B = X1 * Z;
    arma::vec g1_tilde = g1 - B;
    double var1, var2, S, S1, S2, g1tildemu2;
    arma::vec S_a2;
    double Bmu2;
    arma::mat  ZtXVXZ = Z.t() * m_XVX * Z;
    if(m_traitType == "binary" || m_traitType == "count"){
      mu21  = m_mu2.elem(t_indexForNonZero);
      g1tildemu2 = dot(square(g1_tilde), mu21);
      Bmu2 = arma::dot(square(B),  mu21);
      var2 = ZtXVXZ(0,0) - Bmu2 + g1tildemu2;
    }else if(m_traitType == "quantitative" || m_traitType == "count_nb"){
      Bmu2 = dot(g1, B % m_varWeightsvec_new);
      //Bmu2 = dot(g1, B);
      var2 = ZtXVXZ(0,0)*m_tauvec[0] +  dot(g1,g1 % m_varWeightsvec_new) - 2*Bmu2;
      //var2 = ZtXVXZ(0,0)*m_tauvec[0] +  dot(g1,g1) - 2*Bmu2;
    }
    
    var1 = var2 * m_varRatioVal;
    S1 = dot(res1, g1_tilde);
    arma::mat res1X1_temp = (res1.t()) * X1;
    arma::vec res1X1 = res1X1_temp.t();
    S_a2 = m_S_a - res1X1;
    S2 = - arma::dot(S_a2,  Z);
    S = S1 + S2;
    S = S/m_tauvec[0];

    double stat = S*S/var1;
    //double t_pval;
/*
    std::cout << "S FastTest " << S << std::endl;
    std::cout << "var1 " << var1 << std::endl;
    std::cout << "m_varRatioVal " << m_varRatioVal << std::endl;
*/

    //if (var1 <= std::pow(std::numeric_limits<double>::min(), 2)){
    if (var1 <= std::numeric_limits<double>::min()){
        t_pval = 1;
    } else{
      boost::math::chi_squared chisq_dist(1);
      t_pval = boost::math::cdf(complement(chisq_dist, stat));
    }
    

    char pValueBuf[100];
    if (t_pval != 0)
        sprintf(pValueBuf, "%.6E", t_pval);
    else {
        double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
    //std::cout << "t_pval_str scoreTestFast " << t_pval_str << std::endl;
    //std::cout << "end of scoreTestFast" << std::endl;
}



void SAIGEClass::scoreTestFast_noadjCov(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2){

    arma::vec t_GVec_center = t_GVec-arma::mean(t_GVec);
    //arma::vec mu2sub = m_mu2.subvec(t_start,t_end);
    //arma::vec ressub = m_res.subvec(t_start,t_end);
    //arma::vec tauvecsub = m_tauvec.subvec(t_start*m_tauvec_each_length, (t_start+1)*m_tauvec_each_length-1);
    double var2 = dot(m_mu2, pow(t_GVec_center,2));

    double var1 = var2 * m_varRatioVal;
    double S = dot(m_res, t_GVec_center);
    S = S/m_tauvec[0];
    double stat = S*S/var1;
    double t_pval;

    //if (var1 <= std::pow(std::numeric_limits<double>::min(), 2)){
    if (var1 <= std::numeric_limits<double>::min()){
        t_pval = 1;
    } else{
      boost::math::chi_squared chisq_dist(1);
      t_pval = boost::math::cdf(complement(chisq_dist, stat));
       //                                                    t_pval = Rf_pchisq(S*S/var1, 1.0, 0, 0);
    }


    char pValueBuf[100];
    if (t_pval != 0)
        sprintf(pValueBuf, "%.6E", t_pval);
    else {
        double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
    //std::cout << "t_pval_str scoreTestFast " << t_pval_str << std::endl;
    //std::cout << "end of scoreTestFast" << std::endl;
}




void SAIGEClass::scoreTestFast_noadjCov_multiTrait(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2){

    arma::vec t_GVec_center = t_GVec-arma::mean(t_GVec);
    //arma::vec mu2sub = m_mu2.subvec(t_start,t_end);
    //arma::vec ressub = m_res.subvec(t_start,t_end);
    //arma::vec tauvecsub = m_tauvec.subvec(t_start*m_tauvec_each_length, (t_start+1)*m_tauvec_each_length-1);
    double var2 = dot(m_mu2_mt.col(m_itrait), pow(t_GVec_center,2));
    //double var2 = dot(t_GVec_center,t_GVec_center);

    double var1 = var2 * m_varRatioVal;
    //std::cout << "var2 " << var2 << std::endl;
    //std::cout << "var1 " << var1 << std::endl;


    double S = dot(m_res_mt.col(m_itrait), t_GVec_center);
    S = S/m_tauvec_mt(0,m_itrait);
    double stat = S*S/var1;
    double t_pval;
    //std::cout << "S " << S << std::endl; 
    //if (var1 <= std::pow(std::numeric_limits<double>::min(), 2)){
    if (var1 <= std::numeric_limits<double>::min()){
        t_pval = 1;
    } else{
      boost::math::chi_squared chisq_dist(1);
      t_pval = boost::math::cdf(complement(chisq_dist, stat));
       //                                                    t_pval = Rf_pchisq(S*S/var1, 1.0, 0, 0);
    }
    //std::cout << "t_pval " << t_pval << std::endl;

    char pValueBuf[100];
    if (t_pval != 0)
        sprintf(pValueBuf, "%.6E", t_pval);
    else {
        double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
    //std::cout << "t_pval_str scoreTestFast " << t_pval_str << std::endl;
    //std::cout << "end of scoreTestFast" << std::endl;
}




void SAIGEClass::getadjG(arma::vec & t_GVec, arma::vec & g){
   //g = m_XV * t_GVec;
   g = m_XV_mt.rows(m_startip, m_endip) * t_GVec;
   //t_GVec.t().print("t_GVec");
   //std::cout << "sum(t_GVec) " << arma::accu(t_GVec) << std::endl;
   //      m_XV.print("m_XV");
   //g = t_GVec - m_XXVX_inv * g;
   g = t_GVec - m_XXVX_inv_mt.rows(m_startin, m_endin) * g;
    //m_XXVX_inv.print("m_XXVX_inv");
    //g.t().print("g");
}


void SAIGEClass::getadjGFast(arma::vec & t_GVec, arma::vec & g, arma::uvec & iIndex)
{

   arma::vec t_GVec0 = t_GVec;	
   arma::uvec indexNonZeroVec0_arma = iIndex;	
   //if(g_I_longl_mat.n_cols != g_I_longl_mat.n_rows && t_GVec.n_elem != m_XV.n_cols){
      //t_GVec0 = g_I_longl_mat * t_GVec;
      //indexNonZeroVec0_arma = arma::find(t_GVec0 > 0.0);
   //}
   /*else{
      t_GVec0 = t_GVec;
      indexNonZeroVec0_arma = iIndex;
   }*/	   


  // To increase computational efficiency when lots of GVec elements are 0
  arma::vec m_XVG(m_p, arma::fill::zeros);
  for(int i = 0; i < indexNonZeroVec0_arma.n_elem; i++){
      m_XVG += m_XV_mt.rows(m_startip, m_endip).col(indexNonZeroVec0_arma(i)) * t_GVec0(indexNonZeroVec0_arma(i));
  }

   //std::cout << "m_XXVX_inv " << m_XXVX_inv.n_cols << " " << m_XXVX_inv.n_rows << std::endl;
   //std::cout << "m_XVG " << m_XVG.n_cols << " " << m_XVG.n_rows << std::endl;
  //t_GVec0.print("t_GVec0");
  //m_XVG.print("m_XVG");
  //arma::mat m_XXVX_inv_mtsub = m_XXVX_inv_mt.rows(m_startin, m_endin);
  //m_XXVX_inv_mtsub.print("m_XXVX_inv_mtsub");
  
  //g.print("g before");
  g = t_GVec0 - m_XXVX_inv_mt.rows(m_startin, m_endin) * m_XVG;
  //g.print("g after");

}


void SAIGEClass::get_mu(arma::vec & t_mu){
    t_mu = m_mu;
}

void SAIGEClass::getindices(arma::uvec & t_case_indices,
      arma::uvec & t_ctrl_indices){
     t_case_indices = m_case_indices;
     t_ctrl_indices = m_ctrl_indices;
  }


void SAIGEClass::setupSparseMat(int r, arma::umat & locationMatinR, arma::vec & valueVecinR) {
    m_locationMat = locationMatinR;
    m_valueVec = valueVecinR;
    m_dimNum = r;
}



arma::sp_mat SAIGEClass::gen_sp_SigmaMat() {
    //std::cout << "gen_sp_SigmaMat " << std::endl;
    arma::sp_mat resultMat(m_locationMat, m_valueVec, m_dimNum, m_dimNum);
    //std::cout << "gen_sp_SigmaMat2 " << std::endl;
    return resultMat;
}

// revised
// need to add sparse Sigma version 
// This function only uses variance ratio and does not use sparse GRM
void SAIGEClass::getMarkerPval(arma::vec & t_GVec,
			       arma::uvec & iIndex,
			       arma::uvec & iIndexComVec,
                               double& t_Beta, 
                               double& t_seBeta, 
                               double& t_pval,
			       double& t_pval_noSPA, 
                               double t_altFreq,
                               double& t_Tstat,
			       double& t_gy,
			       double& t_var1, 
			       bool & t_isSPAConverge, 
			       arma::vec & t_gtilde, 
			       bool & is_gtilde,
			       bool  is_region,
                               arma::vec & t_P2Vec,
			       bool t_isCondition,
			       double& t_Beta_c,
                           	double& t_seBeta_c,
                           	double& t_pval_c,
                          	double& t_pval_noSPA_c,
                           	double& t_Tstat_c,
                           	double& t_varT_c,
			   	arma::rowvec & t_G1tilde_P_G2tilde, 
				bool & t_isFirth,
				bool & t_isFirthConverge, 
				bool t_isER, 
				bool t_isnoadjCov,
                                bool t_isSparseGRM)
{


  //std::cout << "herehrere" << std::endl;
  t_isFirth = false;
  std::string t_pval_str;
  double t_var2, t_SPApval;
  double altFreq0;
  arma::vec t_GVec0 = t_GVec;
  arma::uvec indexNonZeroVec0_arma = iIndex;
  arma::uvec indexZeroVec0_arma = iIndexComVec;

 if(t_isSparseGRM){
 	t_isnoadjCov = false;
 }
  //std::cout << "t_isnoadjCov " << t_isnoadjCov << std::endl;
  //std::cout << "t_isSparseGRM " << t_isSparseGRM << std::endl;

  if(!t_isnoadjCov){  //t_isSparseGRM
	unsigned int nonzero = indexNonZeroVec0_arma.n_elem;
	unsigned int nGvec = t_GVec0.n_elem;
        is_gtilde = true;
	scoreTest(t_GVec0, t_Beta, t_seBeta, t_pval_str, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde, t_P2Vec, t_gy, is_region, indexNonZeroVec0_arma, t_pval);	
  }else{
  	is_gtilde = false;
	unsigned int nonzero = iIndex.n_elem;
	unsigned int nGvec = t_GVec.n_elem;
	//std::cout << "nGvec " << nGvec << std::endl;
	scoreTestFast_noadjCov_multiTrait(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_str, altFreq0,t_Tstat, t_var1, t_var2);
  }


  double StdStat = std::abs(t_Tstat) / sqrt(t_var1);
  t_isSPAConverge = false;

  double pval_noadj;
  try {
        pval_noadj = std::stod(t_pval_str);
  } catch (const std::invalid_argument&) {
        pval_noadj = 0;
        std::cerr << "Argument is invalid\n";
        //throw;
  } catch (const std::out_of_range&) {
        std::cerr << "Argument is out of range for a double\n";
        //throw;
        pval_noadj = 0;
  }


  t_pval_noSPA = pval_noadj;
  t_pval = pval_noadj;
 double tol1, m1;
if(abs(StdStat) > m_SPA_Cutoff && m_traitType != "quantitative"){
   if(!is_gtilde){
          getadjGFast(t_GVec, t_gtilde, indexNonZeroVec0_arma);
   }
   t_var2 = dot(m_mu2_mt.col(m_itrait), pow(t_gtilde,2));
   arma::vec mu_vec = m_mu_mt.col(m_itrait);
   double q, qinv;
   //= t_Tstat/sqrt(t_var1/var2);
   //double qinv = -q;
   bool logp=false;
   m1 = dot(mu_vec, t_gtilde);
   double tol0 = std::numeric_limits<double>::epsilon();
   tol1 = std::pow(tol0, 0.25);

        if(m_traitType == "binary"){
                q = t_Tstat/sqrt(t_var1/t_var2) + m1;

                if((q-m1) > 0){
                        qinv = -1 * std::abs(q-m1) + m1;
                }else if ((q-m1) == 0){
                        qinv =  m1;
                }else{
                        qinv = std::abs(q-m1) + m1;
                }
        }else if(m_traitType == "survival" || m_traitType == "count"){
                q = t_Tstat/sqrt(t_var1/t_var2);
                qinv = -q;
        }

   SPA(mu_vec, t_gtilde, q, qinv, pval_noadj, tol1, logp, m_traitType, t_SPApval, t_isSPAConverge);
   //if(t_SPApval > pval_noadj){
   t_pval = t_SPApval;
     //t_isSPAConverge = true;
   //}
}


//condition
if(t_isCondition){
	if(!is_gtilde){
        	getadjGFast(t_GVec, t_gtilde, iIndex);
        	is_gtilde = true;
        }

	assignSingleVarianceRatio(m_flagSparseGRM_cur, false, true); //bool isnoXadj, bool issample

        t_G1tilde_P_G2tilde = sqrt(m_varRatioVal) * t_gtilde.t() * (m_P2Mat_cond.cols(m_startic, m_endic));

        arma::vec t_Tstat_ctemp =  t_G1tilde_P_G2tilde * (m_VarInvMat_cond.cols(m_startic, m_endic)) * (m_Tstat_cond.subvec(m_startic, m_endic));
	//std::cout << "here is_gtilde 2" << std::endl;
	arma::mat tempgP2 = t_gtilde.t() * (m_P2Mat_cond.cols(m_startic, m_endic));
	//std::cout << "here is_gtilde 3" << std::endl;

    	t_Tstat_c = t_Tstat - t_Tstat_ctemp(0);
	//std::cout << "here is_gtilde 4" << std::endl;
    	arma::vec t_varT_ctemp = t_G1tilde_P_G2tilde * (m_VarInvMat_cond.cols(m_startic, m_endic)) * (t_G1tilde_P_G2tilde.t());
	//std::cout << "here is_gtilde 5" << std::endl;
    	t_varT_c = t_var1 - t_varT_ctemp(0);
	//std::cout << "here is_gtilde 6" << std::endl;

    double S_c = t_Tstat_c;

    double stat_c = S_c*S_c/t_varT_c;
     if (t_varT_c <= std::numeric_limits<double>::min()){
        t_pval_noSPA_c = 1;
	stat_c = 0;
     }else{
        boost::math::chi_squared chisq_dist(1);
        t_pval_noSPA_c = boost::math::cdf(complement(chisq_dist, stat_c));
     }

    char pValueBuf_c[100];
    if (t_pval_noSPA_c != 0)
        sprintf(pValueBuf_c, "%.6E", t_pval_noSPA_c);
    else {
        double log10p_c = log10(2.0) - M_LOG10E*stat_c/2 - 0.5*log10(stat_c*2*M_PI);
        int exponent_c = floor(log10p_c);
        double fraction_c = pow(10.0, log10p_c - exponent_c);
        if (fraction_c >= 9.95) {
          fraction_c = 1;
           exponent_c++;
         }
        sprintf(pValueBuf_c, "%.1fE%d", fraction_c, exponent_c);
    }
    std::string buffAsStdStr_c = pValueBuf_c;
    std::string& t_pval_noSPA_str_c = buffAsStdStr_c;

    t_Beta_c = S_c/t_varT_c;
    t_seBeta_c = fabs(t_Beta_c) / sqrt(stat_c);
    t_Tstat_c = S_c;
/*
*/

  double pval_noSPA_c;  
  try {
        pval_noSPA_c = std::stod(t_pval_noSPA_str_c);
  } catch (const std::invalid_argument&) {
        pval_noSPA_c = 0;
        std::cerr << "Argument is invalid\n";
        //throw;
  } catch (const std::out_of_range&) {
        std::cerr << "Argument is out of range for a double\n";
        //throw;
        pval_noSPA_c = 0;
  }
  t_pval_noSPA_c = pval_noSPA_c; 

  double StdStat_c = std::abs(t_Tstat_c) / sqrt(t_varT_c);

  if(abs(StdStat_c) > m_SPA_Cutoff && m_traitType != "quantitative" && t_varT_c > std::numeric_limits<double>::min()){
     bool t_isSPAConverge_c;
     double q_c, qinv_c, pval_noadj_c, SPApval_c;
     arma::vec m_mu_vec = m_mu_mt.col(m_itrait);
       t_var2 = dot(m_mu2_mt.col(m_itrait), pow(t_gtilde,2));
      if(m_traitType == "binary"){
                q_c = t_Tstat_c/sqrt(t_varT_c/t_var2) + m1;

                if((q_c-m1) > 0){
                        qinv_c = -1 * std::abs(q_c-m1) + m1;
                }else if ((q_c-m1) == 0){
                        qinv_c =  m1;
                }else{
                        qinv_c = std::abs(q_c-m1) + m1;
                }
        }else if(m_traitType == "count"){
                q_c = t_Tstat_c/sqrt(t_varT_c/t_var2);
                qinv_c = -q_c;
        }	
     bool logp=false;
     SPA(m_mu_vec, t_gtilde, q_c, qinv_c, pval_noadj_c, tol1, logp, m_traitType, SPApval_c, t_isSPAConverge_c);

        boost::math::normal ns;
        t_pval_c = SPApval_c;
        double t_qval_c;
        try {
           t_qval_c = boost::math::quantile(ns, t_pval_c/2);
           t_qval_c = fabs(t_qval_c);
           t_seBeta_c = fabs(t_Beta_c)/t_qval_c;
        }catch (const std::overflow_error&) {
          t_qval_c = std::numeric_limits<double>::infinity();
          t_seBeta_c = 0;
        }
  }else{
     t_pval_c = t_pval_noSPA_c;
  }
} //if(t_isCondition){





   //std::cout << "is_region " << is_region << std::endl;
  if(is_region && !is_gtilde){
	    //arma::vec t_GVec0 = g_I_longl_mat * t_GVec;
	    //arma::uvec iIndex0 = arma::find(t_GVec0 > 0.0);
      getadjGFast(t_GVec, t_gtilde, iIndex);
	//getadjGFast(t_GVec0, t_gtilde, iIndex0);
      is_gtilde = true; 
  }

    //if(is_region && isScoreFast){
    //if(is_region && (t_isnoadjCov || (!t_isnoadjCov && !t_isSparseGRM))){
    if(is_region){


      //if(!m_flagSparseGRM_cur){
      if(!t_isSparseGRM){
        t_P2Vec = t_gtilde % (m_mu2_mt.col(m_itrait)) *(m_tauvec_mt(0,m_itrait));
      }
      //else{
        //if(t_isnoadjCov){
	//  t_P2Vec = t_gtilde % (m_mu2_mt.col(m_itrait)) *(m_tauvec_mt(0,m_itrait));		
        //}
       //if(m_SigmaMat_sp.n_rows > 2){
       // 	t_P2Vec = m_SigmaMat_sp * t_gtilde;
		
       // }else{
       //  	t_P2Vec = getSigma_G_V(t_gtilde, 500, 1e-5);
	//	//m_sigmainvG_noV;
       //}
        //arma::sp_mat m_SigmaMat_sp = gen_sp_SigmaMat();
        //t_P2Vec = arma::spsolve(m_SigmaMat_sp, t_gtilde);
      //}
    }
    //std::cout << "end of getPval" << std::endl;
}


bool SAIGEClass::assignVarianceRatio(double MAC, bool issparseforVR, bool isnoXadj){
    bool hasVarRatio = false;
    arma::vec m_varRatio;
    if(issparseforVR){
	//m_varRatio = m_varRatio_sparse;
	m_varRatio = m_varRatio_sparse_mt.col(m_itrait);
    }else{
	if(!isnoXadj){    
	    //m_varRatio = m_varRatio_null;
	    m_varRatio = m_varRatio_null_mt.col(m_itrait);
	}else{
	//m_varRatio_null_noXadj.print("m_varRatio_null_noXadj");
	    m_varRatio = m_varRatio_null_noXadj_mt.col(m_itrait);
	    //m_varRatio.print("m_varRatio");
	}	
    }
   

    //m_cateVarRatioMinMACVecExclude.print("m_cateVarRatioMinMACVecExclude");
    //m_cateVarRatioMaxMACVecInclude.print("m_cateVarRatioMaxMACVecInclude");
    //m_varRatio.print("m_varRatio");
    //m_varRatio_null_noXadj_mt.print("m_varRatio_null_noXadj_mt");
    for(unsigned int i = 0; i < m_cateVarRatioMaxMACVecInclude.n_elem; i++)
    {
        if(MAC <= m_cateVarRatioMaxMACVecInclude(i) && MAC > m_cateVarRatioMinMACVecExclude(i)){    	    
		m_varRatioVal = m_varRatio(i);
		hasVarRatio = true;
	}	
    }

    if(!hasVarRatio){	
	if(MAC < m_cateVarRatioMinMACVecExclude(0)){
		m_varRatioVal = m_varRatio(0);
		hasVarRatio = true;
	}	
    }

    if(!hasVarRatio){
        if(MAC > m_cateVarRatioMaxMACVecInclude.back()){
		//m_varRatioVal = m_varRatio(a-1);
		m_varRatioVal = m_varRatio.back();
                hasVarRatio = true;
        }
    }

   //m_varRatioVal = m_varRatio(0);
   hasVarRatio = true;
   //std::cout << "hasVarRatio " << hasVarRatio << std::endl;
   return(hasVarRatio);    
}


bool SAIGEClass::assignVarianceRatio(double MAC, bool issparseforVR, bool isnoXadj, bool issample){
    bool hasVarRatio = false;
    arma::vec m_varRatio;
    if(issparseforVR){
        //m_varRatio = m_varRatio_sparse;
        m_varRatio = m_varRatio_sparse_mt.col(m_itrait);
    }else{
        if(!isnoXadj){
            //m_varRatio = m_varRatio_null;
	   if(!issample){
            m_varRatio = m_varRatio_null_mt.col(m_itrait);
           }else{
            m_varRatio = m_varRatio_null_sample_mt.col(m_itrait);
	   }
        }else{
        //m_varRatio_null_noXadj.print("m_varRatio_null_noXadj");
            m_varRatio = m_varRatio_null_noXadj_mt.col(m_itrait);
            //m_varRatio.print("m_varRatio");
        }
    }


    //m_cateVarRatioMinMACVecExclude.print("m_cateVarRatioMinMACVecExclude");
    //m_cateVarRatioMaxMACVecInclude.print("m_cateVarRatioMaxMACVecInclude");
    //m_varRatio.print("m_varRatio");
    //m_varRatio_null_noXadj_mt.print("m_varRatio_null_noXadj_mt");
    for(unsigned int i = 0; i < m_cateVarRatioMaxMACVecInclude.n_elem; i++)
    {
        if(MAC <= m_cateVarRatioMaxMACVecInclude(i) && MAC > m_cateVarRatioMinMACVecExclude(i)){
                m_varRatioVal = m_varRatio(i);
                hasVarRatio = true;
        }
    }

    if(!hasVarRatio){
        if(MAC < m_cateVarRatioMinMACVecExclude(0)){
                m_varRatioVal = m_varRatio(0);
                hasVarRatio = true;
        }
    }

    if(!hasVarRatio){
        if(MAC > m_cateVarRatioMaxMACVecInclude.back()){
                //m_varRatioVal = m_varRatio(a-1);
                m_varRatioVal = m_varRatio.back();
                hasVarRatio = true;
        }
    }

   //m_varRatioVal = m_varRatio(0);
   hasVarRatio = true;
   //std::cout << "hasVarRatio " << hasVarRatio << std::endl;
   return(hasVarRatio);
}




arma::vec SAIGEClass::assignVarianceRatioi_multi(arma::vec & MACvec, bool issparseforVR, bool isnoXadj){
    bool hasVarRatio = false;
    arma::vec m_varRatio;
    if(issparseforVR){
        m_varRatio = m_varRatio_sparse;
    }else{
        if(!isnoXadj){
            m_varRatio = m_varRatio_null;
        }else{
            m_varRatio = m_varRatio_null_noXadj;
        }
    }

    unsigned int q = MACvec.n_elem;
    arma::vec varRatio_vec_multi(q);

 double MAC;
 for(unsigned int j = 0; j < q; j++){ 
    MAC = MACvec(j);
    for(unsigned int i = 0; i < m_cateVarRatioMaxMACVecInclude.n_elem; i++)
    {
        if(MAC <= m_cateVarRatioMaxMACVecInclude(i) && MAC > m_cateVarRatioMinMACVecExclude(i)){
                m_varRatioVal = m_varRatio(i);
                hasVarRatio = true;
        }
    }

    if(!hasVarRatio){
        if(MAC < m_cateVarRatioMinMACVecExclude(0)){
                m_varRatioVal = m_varRatio(0);
                hasVarRatio = true;
        }
    }

    if(!hasVarRatio){
        if(MAC > m_cateVarRatioMaxMACVecInclude.back()){
                //m_varRatioVal = m_varRatio(a-1);
                m_varRatioVal = m_varRatio.back();
                hasVarRatio = true;
        }
    }
    varRatio_vec_multi(j) = m_varRatioVal;
 } 
    return(varRatio_vec_multi);
}








void SAIGEClass::assignSingleVarianceRatio(bool issparseforVR, bool isnoXadj){ 
    arma::rowvec m_varRatio;
    if(issparseforVR){
        m_varRatio = m_varRatio_sparse_mt.row(0);
    }else{
	if(isnoXadj){    
            m_varRatio = m_varRatio_null_noXadj_mt.row(0);
	}else{
	    m_varRatio = m_varRatio_null_mt.row(0);	
	}	
    }	
    m_varRatioVal = m_varRatio(m_itrait);
}


void SAIGEClass::assignSingleVarianceRatio(bool issparseforVR, bool isnoXadj, bool issample){
    arma::rowvec m_varRatio;
    if(issparseforVR){
        m_varRatio = m_varRatio_sparse_mt.row(0);
    }else{
        if(isnoXadj){
            m_varRatio = m_varRatio_null_noXadj_mt.row(0);
        }else{
	  if(!issample){
            m_varRatio = m_varRatio_null_mt.row(0);
	  }else{
	    m_varRatio = m_varRatio_null_sample_mt.row(0);
	  }
        }
    }
    m_varRatioVal = m_varRatio(m_itrait);
}

void SAIGEClass::assignSingleVarianceRatio_withinput(double t_varRatioVal){
        m_varRatioVal = t_varRatioVal;
}


void SAIGEClass::assignConditionFactors(
      arma::mat & t_P2Mat_cond,
      arma::mat & t_VarInvMat_cond,
      arma::mat & t_VarMat_cond,
      arma::vec & t_Tstat_cond,
      arma::mat & t_G2_Weight_cond,
      arma::vec & t_MAF_cond,
      arma::mat & t_qsum_cond,
      arma::mat & t_gsum_cond,
      arma::vec & t_p_cond
      ){
	m_P2Mat_cond = t_P2Mat_cond;
	m_VarInvMat_cond = t_VarInvMat_cond;
	m_VarMat_cond = t_VarMat_cond;
	m_Tstat_cond = t_Tstat_cond;
	m_MAF_cond = t_MAF_cond;
	m_qsum_cond_Mat = t_qsum_cond;
	m_gsum_cond_Mat = t_gsum_cond;
	m_G2_Weight_cond = t_G2_Weight_cond;


	m_p_cond = t_p_cond;
	m_numMarker_cond = t_Tstat_cond.n_elem;
}

void SAIGEClass::assignConditionFactors_scalefactor( 
	arma::vec & t_scalefactor_G2_cond	
		){
	m_scalefactor_G2_cond = t_scalefactor_G2_cond;
	arma::mat scalefactor_G2_cond_Mat = arma::diagmat(arma::sqrt(m_scalefactor_G2_cond));
	arma::mat weightMat_G2_G2 = m_G2_Weight_cond * m_G2_Weight_cond.t(); 
	arma::mat VarMat_cond_scaled = scalefactor_G2_cond_Mat * m_VarMat_cond * scalefactor_G2_cond_Mat;
	arma::mat VarMat_cond_scaled_weighted = VarMat_cond_scaled % weightMat_G2_G2;
	m_VarInvMat_cond_scaled_weighted = VarMat_cond_scaled_weighted.i();
	//m_VarInvMat_cond_region_binary = (1/scalefactor_G2_cond_Mat) * m_VarInvMat_cond	* (1/scalefactor_G2_cond_Mat);
	
}


void SAIGEClass::assignConditionFactors_scalefactor_multiTrait(
        arma::mat & t_scalefactor_G2_cond,
	unsigned int oml
                ){
	unsigned int startoml = oml * (m_VarMat_cond.n_rows);
	unsigned int endoml = (oml+1) * (m_VarMat_cond.n_rows) - 1;
	m_scalefactor_G2_cond_Mat = t_scalefactor_G2_cond;
		
        m_scalefactor_G2_cond = t_scalefactor_G2_cond.col(oml);
        arma::mat scalefactor_G2_cond_Mat = arma::diagmat(arma::sqrt(m_scalefactor_G2_cond));
        arma::mat weightMat_G2_G2 = m_G2_Weight_cond * m_G2_Weight_cond.t();
        arma::mat VarMat_cond_scaled = scalefactor_G2_cond_Mat * (m_VarMat_cond.cols(startoml, endoml)) * scalefactor_G2_cond_Mat;
        arma::mat VarMat_cond_scaled_weighted = VarMat_cond_scaled % weightMat_G2_G2;
        m_VarInvMat_cond_scaled_weighted.cols(startoml, endoml) = VarMat_cond_scaled_weighted.i();
        //m_VarInvMat_cond_region_binary = (1/scalefactor_G2_cond_Mat) * m_VarInvMat_cond       * (1/scalefactor_G2_cond_Mat);
}


void SAIGEClass::extract_XV_XXVX_inv(arma::mat & t_XV, arma::mat & t_XXVX_inv){
	//std::cout << "m_XV_mt.n_rows " << m_XV_mt.n_rows << " " << m_XV_mt.n_cols << std::endl;
	//std::cout << "m_XXVX_inv_mt.n_rows " <<m_XXVX_inv_mt.n_rows << " " << m_XXVX_inv_mt.n_cols  << std::endl;
	t_XV = m_XV_mt.rows(m_startip, m_endip);
	t_XXVX_inv = m_XXVX_inv_mt.rows(m_startin, m_endin);

	//std::cout << " t_XV.n_rows " << t_XV.n_rows << " " << t_XV.n_cols << std::endl;
	//std::cout << " t_XXVX_inv.n_rows " << t_XXVX_inv.n_rows << " " << t_XXVX_inv.n_cols << std::endl;


} 



void SAIGEClass::fast_logistf_fit_simple(arma::mat & x,
                arma::vec & y,
                arma::vec & offset,
                bool firth,
        arma::vec init,
        int maxit,
        int maxstep,
        int maxhs,
        double lconv,
        double gconv,
        double xconv,
        double & beta_G,
        double & sebeta_G, 
	bool & isfirthconverge){
  isfirthconverge = false;	
  int n = x.n_rows;
  int k = x.n_cols;
  arma::vec beta = init;
  int iter = 0;
  arma::vec pi_0 = -x * beta - offset;
  pi_0 = arma::exp(pi_0) + 1;
  arma::vec pi = 1/pi_0;
  int evals = 1;
  arma::vec beta_old;
  arma::mat oneVec(k, 1 , arma::fill::ones);
  arma::mat XX_covs(k, k, arma::fill::zeros);
  while(iter <= maxit){
        beta_old = beta;
        arma::vec wpi = pi % (1 - pi);
        arma::vec W2 = arma::sqrt(wpi);
        //arma::vec wpi_sqrt = arma::sqrt(wpi);
        //arma::vec W2 = weight % wpi_sqrt;
        arma::mat XW2(n, k, arma::fill::zeros);
        for(int j = 0; j < k; j++){
                XW2.col(j) = x.col(j) % W2;
        }

        arma::mat Q;
        arma::mat R;
        arma::qr_econ(Q, R, XW2);
        arma::vec h = Q % Q * oneVec;
        arma::vec U_star(2, arma::fill::zeros);
        arma::vec ypih;
        if(firth){
                ypih = (y - pi) + (h % (0.5 - pi));
        }else{
                ypih = (y - pi);
        }
        //ypih.print();
        arma::vec xcol(n, arma::fill::zeros);
        U_star = x.t() * ypih;

        arma::mat XX_XW2(n, k, arma::fill::zeros);
        for(int j = 0; j < k; j++){
                xcol = x.col(j);
                XX_XW2.col(j) = xcol % W2;
        }
        arma::mat XX_Fisher = XX_XW2.t() * (XX_XW2);
        bool isinv = arma::inv_sympd (XX_covs, XX_Fisher);
        
	if(!isinv){
                break;
        }
        //}
        arma::vec delta = XX_covs * U_star;
        //delta.replace(arma::datum::nan, 0);

        double mx = arma::max(arma::abs(delta))/maxstep;
        if(mx > 1){
                delta = delta/mx;
        }
        evals = evals + 1;
        iter = iter + 1;
        beta = beta + delta;
        pi_0 = -x * beta - offset;
        pi_0 = arma::exp(pi_0) + 1;
        pi = 1/pi_0;
        if((iter == maxit) || ( (arma::max(arma::abs(delta)) <= xconv) & (abs(U_star).is_zero(gconv)))){
		isfirthconverge = true;
                break;
        }
  }
        arma::mat var;
        if(XX_covs.has_nan()){
                var = XX_covs;
                beta_G = arma::datum::nan;
                sebeta_G = arma::datum::nan;
        }else{
                beta_G = beta(1);
                sebeta_G = sqrt(XX_covs(1,1));
        }

	//std::cout << "beta_G " << beta_G << std::endl;
	//std::cout << "sebeta_G " << sebeta_G << std::endl;
        //return beta;
}

void SAIGEClass::set_flagSparseGRM_cur(bool t_flagSparseGRM_cur){
	m_flagSparseGRM_cur = t_flagSparseGRM_cur;
}




arma::vec SAIGEClass::getDiagOfSigma_V(){
	std::cout << "getDiagOfSigma_V() 0" << std::endl;
	double tauVal0 = m_tauvec_mt(0,m_itrait);
	//m_mu2_mt.col(m_itrait)
        int Nnomissing = m_mu2_mt.n_rows;
        arma::vec diagVec(Nnomissing);
        arma::sp_vec diagVecV0;
        arma::vec diagVecG, diagVecV, diagVecG_I, diagVecG_T, diagVecG_IT,diagVecV_I, diagVecV_T, diagVecV_IT;
        unsigned int tauind = 0;
	std::cout << "getDiagOfSigma_V() " << std::endl;
        if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){

             //diagVec = tauVal0/m_mu2;
             diagVec = 1/(m_mu2_mt.col(m_itrait));
             tauind = tauind + 1;
             //diagVecG = spV.diag();
             diagVec = diagVec + m_varRatioVal;
             tauind = tauind + 1;


        }else{ //if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){
                //diagVec = tauVal0/m_mu2;
                diagVec = 1/(m_mu2_mt.col(m_itrait));
                //tauind = tauind + 1;
                //diagVecG = spV.diag();
                //diagVecG_I = diagVecG.elem(g_I_longl_vec);
                //diagVec = diagVec + tauVec(indforV) * diagVecG_I;
                diagVec = diagVec + m_varRatioVal;

                /*
                tauind = tauind + 1 + 3 * (indforV - 1);
                tauind = tauind + 1;

                if(g_T_longl_mat.n_rows > 0){
                  diagVecG_I.ones(Nnomissing);
                  diagVecG_IT = diagVecG_I % g_T_longl_vec;
                  diagVecG_T = diagVecG_IT % g_T_longl_vec;
                  diagVecG_IT = 2 * diagVecG_IT;
                  diagVec = diagVec + tauVec(tauind) * diagVecG_IT;
                  tauind = tauind + 1;
                  diagVec = diagVec + tauVec(tauind) * diagVecG_T;
                  tauind = tauind + 1;
                }
                */
        }

        for(unsigned int i=0; i< Nnomissing; i++){
                if(diagVec(i) < 1e-4){
                        diagVec(i) = 1e-4 ;
                }
        }

        return(diagVec);
}







arma::vec SAIGEClass::getCrossprod_V(arma::vec& bVec){
        //indforV starts from 0
        arma::vec crossProdVec;
        arma::vec crossProd1, GRM_I_bvec, Ibvec, Tbvec, GRM_T_bvec, crossProdGRM_TGIb, crossProdGRM_IGTb, V_I_bvec, V_T_bvec, crossProdV_TGIb, crossProdV_IGTb, crossProdGRM_TIb, crossProdGRM_ITb;
	double tauVal0 = m_tauvec(0); 

        unsigned int tau_ind = 0;
        if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){ //it must have specified GRM
                  crossProd1 = bVec;
                  //crossProdVec = tauVal0*(bVec % (1/m_mu2)) + m_tauVal_sp*crossProd1;
                  crossProdVec = (bVec % (1/m_mu2)) + m_tauVal_sp*crossProd1;
                  tau_ind = tau_ind + 2;
        }else{

                  Ibvec = g_I_longl_mat.t() * bVec;

                  GRM_I_bvec = Ibvec;
                  crossProd1 = g_I_longl_mat * GRM_I_bvec;
                  //crossProdVec = tauVal0*(bVec % (1/m_mu2)) + m_tauVal_sp*crossProd1;
                  crossProdVec = (bVec % (1/m_mu2)) + m_tauVal_sp*crossProd1;
                /*
                  tau_ind = tau_ind + 1 + 3 * (indforV - 1);
                  tau_ind = tau_ind + 1;

                if(g_T_longl_mat.n_rows > 0){
                        Tbvec = g_T_longl_mat.t() * bVec;
                        GRM_T_bvec = spV * Tbvec;
                        crossProdGRM_TGIb = g_T_longl_mat * GRM_I_bvec;
                        crossProdGRM_IGTb = g_I_longl_mat * GRM_T_bvec;
                        crossProdVec = crossProdVec + tauVec(tau_ind) * (crossProdGRM_TGIb + crossProdGRM_IGTb);
                        tau_ind = tau_ind + 1;
                        crossProdVec = crossProdVec + tauVec(tau_ind) * (g_T_longl_mat * GRM_T_bvec);
                        tau_ind = tau_ind + 1;
                }
                */

        }

        return(crossProdVec);
}



arma::vec SAIGEClass::getPCG1ofSigmaAndVector_V(arma::vec& bVec, int maxiterPCG, double tolPCG){
    // Start Timers
    //double wall0 = get_wall_time();
    //double cpu0  = get_cpu_time();
    int Nnomissing = m_mu2.n_elem;
    arma::vec xVec(Nnomissing);
    xVec.zeros();
   std::cout << "getPCG1ofSigmaAndVector_V " << std::endl; 
    //if(g_isStoreSigma){
    //    std::cout << " arma::spsolve(g_spSigma, bVec) 0" << std::endl;
    //    //xVec = arma::spsolve(g_spSigma, bVec);
    //    xVec = arma::spsolve(g_spSigma_V, bVec);
    //    std::cout << " arma::spsolve(g_spSigma, bVec) 1" << std::endl;
    //}else{
        arma::vec rVec = bVec;
        arma::vec r1Vec;
        arma::vec crossProdVec(Nnomissing);
        arma::vec zVec(Nnomissing);
        arma::vec minvVec(Nnomissing);
   std::cout << "getPCG1ofSigmaAndVector_V 2" << std::endl; 

        minvVec = 1/getDiagOfSigma_V();

	std::cout <<  "rVec.n_elem " << rVec.n_elem << std::endl;
	std::cout <<  "minvVec.n_elem " << minvVec.n_elem << std::endl;
        zVec = minvVec % rVec;

        float sumr2 = sum(rVec % rVec);
        arma::vec z1Vec(Nnomissing);
        arma::vec pVec = zVec;

        int iter = 0;

        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                arma::colvec ApVec = getCrossprod_V(pVec);
                arma::vec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                double a = preA(0);
                xVec = xVec + a * pVec;
                r1Vec = rVec - a * ApVec;
                z1Vec = minvVec % r1Vec;

                arma::vec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                double bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;
                sumr2 = sum(rVec % rVec);
        }

        if (iter >= maxiterPCG){
		std::cout << "pcg did not converge. You may increase maxiter number." << std::endl;

        }
	std::cout << "iter from getPCG1ofSigmaAndVector " << iter << std::endl;
//}
        return(xVec);
}


arma::vec SAIGEClass::getSigma_G_V(arma::vec & bVec, int maxiterPCG, double tolPCG){
	arma::vec Sigma_iG;
        Sigma_iG = getPCG1ofSigmaAndVector_V(bVec, maxiterPCG, tolPCG);
        return(Sigma_iG);
}



double SAIGEClass::EmpSPA_K_0(double t, 
             int N0, 
             double adjG0, 
             arma::vec & adjG1)        // adjusted Genotype 
{
    double t_adjG0 = t * adjG0;
    arma::vec t_adjG1 = t * adjG1;
    double out = N0 * m_K_0_emp.getValue(t_adjG0) + arma::sum(m_K_0_emp.getVector(t_adjG1));
    return out;
}

double SAIGEClass::EmpSPA_K_1(double t,
             int N0, 
             double adjG0, 
             arma::vec & adjG1,        // adjusted Genotype
             double q2)
{
    double t_adjG0 = t * adjG0;
    arma::vec t_adjG1 = t * adjG1;
    double out = N0 * adjG0 * m_K_1_emp.getValue(t_adjG0) + arma::sum(adjG1 % m_K_1_emp.getVector(t_adjG1)) - q2;
    return out;
}


double SAIGEClass::EmpSPA_K_2(double t, 
             int N0, 
             double adjG0, 
             arma::vec & adjG1)       // adjusted Genotype
{
    double t_adjG0 = t * adjG0;
    arma::vec t_adjG1 = t * adjG1;
    double out = N0 * pow(adjG0, 2) * m_K_2_emp.getValue(t_adjG0) + arma::sum(pow(adjG1, 2) % m_K_2_emp.getVector(t_adjG1));
    return out;
}

Rcpp::List SAIGEClass::EmpSPA_fastgetroot_K1(double t_initX,
                            int N0, 
                            double adjG0, 
                            arma::vec & adjG1,        // adjusted Genotype
                            double q2)
  {
    double x = t_initX, oldX;
    double K1 = 0, K2 = 0, oldK1;
    double diffX = arma::datum::inf, oldDiffX;
    bool converge = true;
    double tol = 0.001;
    int maxiter = 100;
    int iter = 0;
    
    for(iter = 0; iter < maxiter; iter ++){
      
      oldX = x;
      oldDiffX = diffX;
      oldK1 = K1;
      
      K1 = EmpSPA_K_1(x, N0, adjG0, adjG1, q2);
      K2 = EmpSPA_K_2(x, N0, adjG0, adjG1);
      
      diffX = -1 * K1 / K2;
      
      // Checked on 03/25/2021: Expected!!!!
      // std::cout << "x:\t" << x << std::endl;
      // std::cout << "K1:\t" << K1 << std::endl;
      // std::cout << "arma::sign(K1):\t" << arma::sign(K1) << std::endl;
      // std::cout << "arma::sign(oldK1):\t" << arma::sign(oldK1) << std::endl;
      // std::cout << "K2:\t" << K2 << std::endl;
      // std::cout << "diffX:\t" << diffX << std::endl;
      
      if(!std::isfinite(K1)){
        // checked it on 07/05:
        // if the solution 'x' tends to infinity, 'K2' tends to 0, and 'K1' tends to 0 very slowly.
        // then we can set the one sided p value as 0 (instead of setting converge = F)
        x = arma::datum::inf;
        K2 = 0;
        break;
      }
      
      if(arma::sign(K1) != arma::sign(oldK1)){
        while(std::abs(diffX) > std::abs(oldDiffX) - tol){
          diffX = diffX / 2;
        }
      }
      
      if(std::abs(diffX) < tol) break;
      
      x = oldX + diffX;
    }
    
    if(iter == maxiter) 
      converge = false;

    Rcpp::List yList = Rcpp::List::create(Rcpp::Named("root") = x,
                                          Rcpp::Named("iter") = iter,
                                          Rcpp::Named("converge") = converge,
                                          Rcpp::Named("K2") = K2);
    return yList;
  }

double SAIGEClass::EmpSPA_GetProb_SPA(double adjG0, 
                     arma::vec & adjG1,
                     int N0, 
                     double q2, 
                     bool lowerTail)
  {
    double initX = 0;
    
    // The following initial values are validated on 03/25/2021
    if(q2 > 0) initX = 3;
    if(q2 <= 0) initX = -3;
    
    Rcpp::List rootList = EmpSPA_fastgetroot_K1(initX, N0, adjG0, adjG1, q2);
    double zeta = rootList["root"];
    
    double k1 = EmpSPA_K_0(zeta,  N0, adjG0, adjG1);
    double k2 = EmpSPA_K_2(zeta,  N0, adjG0, adjG1);
    double temp1 = zeta * q2 - k1;
    
    double w = arma::sign(zeta) * sqrt(2 * temp1);
    double v = zeta * sqrt(k2);
    
    double pval = arma::normcdf(arma::sign(lowerTail-0.5) * (w + 1/w * log(v/w)));

    return pval;
  }

double SAIGEClass::EmpSPA_getMarkerPval(arma::vec & t_g,
                       double t_Tstat,
		       double g_altFreq_new, 
		       int g_N)
{

    // estimated variance after adjusting for covariates

    arma::vec adjGVec2 = pow(t_g, 2);
    double VarS = m_varResid * sum(adjGVec2) * m_varRatioVal;
    double VarS_rawG = m_varResid*2*g_altFreq_new*(1-g_altFreq_new)*g_N * m_varRatioVal; 
    //double t_zScore = t_Tstat / sqrt(VarS);
    double t_zScore = t_Tstat;
    //arma::vec adjGVecNorm = t_g / sqrt(VarS);
    arma::vec adjGVecNorm = t_g / sqrt(sum(adjGVec2));
    double pval_noEmpSPA = arma::normcdf(-1*std::abs(t_zScore))*2;
    
    std::cout << "pval_noEmpSPA " << pval_noEmpSPA << std::endl;
    std::cout << "t_Tstat " << t_Tstat << std::endl;
    std::cout << "VarS " << VarS << std::endl;
    std::cout << "VarS_rawG " << VarS_rawG << std::endl;
    std::cout << "m_varResid " << m_varResid << std::endl;
    std::cout << "m_varRatioVal " << m_varRatioVal << std::endl;

    double varg = arma::mean(t_g);
    std::cout << "varg " << varg << std::endl;

    int N0 = 0;
    double adjG0 = 0;   // since N0=0, this value actually does not matter

    double pval1 = EmpSPA_GetProb_SPA(adjG0, adjGVecNorm, N0, std::abs(t_zScore), false);
    double pval2 = EmpSPA_GetProb_SPA(adjG0, adjGVecNorm, N0, -1*std::abs(t_zScore), true);
    double pval = pval1 + pval2;
    std::cout << "pval " << pval << std::endl;
    return pval;
  }


void SAIGEClass::getMarkerPval_multi(arma::mat & t_GMat,
				arma::vec & t_MAF_vec,
				arma::vec & t_MAC_vec,
                               std::vector <double> & t_Beta,
                               std::vector <double> & t_seBeta,
                               std::vector <double> & t_pval,
                               std::vector <double> & t_pval_noSPA,
                               std::vector <double> & t_Tstat,
                               std::vector <double> & t_var1, 
			       std::vector <std::string> & t_pval_str_vec)
{
  arma::vec t_Beta_vec, t_seBeta_vec, t_Tstat_vec, t_var1_vec, t_var2_vec, t_pval_noadj_vec, t_pval_vec;
  arma::uvec t_fastTestind_vec, t_SPAind_vec;
  //std::cout << "t_GMat here2a " << std::endl; 

  arma::vec varRatio_multi_vec;
  if(m_varRatio_null.n_elem == 1){
        assignSingleVarianceRatio(false, true);
  }else{
  	varRatio_multi_vec = assignVarianceRatioi_multi(t_MAC_vec, false, true);
  }


/*
double t_Tstat_compare, t_Tstat_compare_tilde, t_Tstat_compare_mean;
arma::vec t_GMatraw;
arma::vec t_GMatraw_tilde;
arma::vec t_GMatraw_mean;
arma::uvec indexForNonZero_compare;
for(int i = 0; i < 100; i++){
	t_GMatraw = t_GMat.col(i);
	indexForNonZero_compare = arma::find(t_GMatraw != 0); 
    	getadjGFast(t_GMatraw, t_GMatraw_tilde, indexForNonZero_compare);
	t_GMatraw_mean = t_GMatraw - arma::mean(t_GMatraw);
std::cout << "t_GMatraw.n_elem " << t_GMatraw.n_elem << std::endl;
	std::cout << "m_res.n_elem " << m_res.n_elem << std::endl;
	t_Tstat_compare = arma::dot(t_GMatraw,  m_res) / (m_tauvec[0]);
	t_Tstat_compare_tilde = arma::dot(t_GMatraw_tilde,  m_res) / (m_tauvec[0]);
	t_Tstat_compare_mean = arma::dot(t_GMatraw_mean,  m_res) / (m_tauvec[0]);
	std::cout << i << ": " << t_Tstat_compare << ", " << t_Tstat_compare_tilde << ", " << t_Tstat_compare_mean << std::endl;
	indexForNonZero_compare.clear();
}
*/


  //arma::mat t_GMat_center = t_GMat.each_row() - t_MAF_vec.t();
  scoreTestFast_noadjCov_multi(t_GMat, t_MAF_vec, t_Beta_vec, t_seBeta_vec, t_pval_str_vec, t_Tstat_vec, t_var1_vec, t_var2_vec, t_pval_noadj_vec, varRatio_multi_vec);


  //arma::vec StdStat_vec = arma::abs(t_Tstat_vec) / arma::sqrt(t_var1_vec);
  //t_fastTestind_vec = arma::find(t_pval_noadj_vec < m_pval_cutoff_for_fastTest);
  
  //std::cout << "t_GMat here2b " << std::endl; 
  t_pval_vec = t_pval_noadj_vec;

  //t_pval_noadj_vec(t_fastTestind_vec).print("t_pval_noadj_vec");
  //t_fastTestind_vec.print("t_fastTestind_vec");

  //arma::vec timeoutput3_a = getTime();

/*

  if(t_fastTestind_vec.n_elem > 0){

  	if(m_varRatio_null.n_elem == 1){
        	assignSingleVarianceRatio(false, false);
  	}


	arma::vec t_Beta_vec_sub, t_seBeta_vec_sub, t_pval_vec_sub, t_Tstat_vec_sub, t_var1_vec_sub, t_var2_vec_sub, t_gy_vec_sub;
	arma::mat t_P2Mat_sub;
	std::vector <std::string> t_pval_str_vec_sub;
		
	arma::vec t_gtildeVec, t_P2Vec;
	double t_gy;
        t_GMat_sub = t_GMat.cols(t_fastTestind_vec);



        double altFreq0;
        arma::vec t_GVec0;
        arma::uvec indexNonZeroVec0_arma, indexZeroVec0_arma;

	
        if((g_I_longl_mat.n_cols != g_I_longl_mat.n_rows) && (t_GMat_sub.n_rows < m_y.n_elem)){
        	t_GMat0_sub = g_I_longl_mat * t_GMat_sub;
   	}else{
		t_GMat0_sub = t_GMat_sub;
   	}
	getadjGFast_multi(t_GMat_sub, t_GMat_tilde_sub);
}


*/



/*
	arma::vec t_varRatio_vec_sub;
	if(m_varRatio_null.n_elem == 1){
        	assignSingleVarianceRatio(false, true);
  	}else{
		t_varRatio_vec_sub = varRatio_multi_vec(t_fastTestind_vec); //need to correct
  	}
	


	scoreTest_multi(t_GMat_tilde_sub,
			t_Beta_vec_sub,
			t_seBeta_vec_sub,
			t_pval_str_vec_sub,
			t_pval_vec_sub,
			t_Tstat_vec_sub,
			t_var1_vec_sub,
			t_var2_vec_sub,
			t_P2Mat_sub,
			t_gy_vec_sub,
			t_varRatio_vec_sub,
			false);
}


	for(int i = 0; i < t_fastTestind_vec.n_elem; i++){	
		int j = t_fastTestind_vec(i);
		t_GVec0 = t_GMat0_sub.col(i);
		t_altFreq0_val = arma::mean(t_GVec0);
      		iIndex = arma::find(t_GVec0 != 0);
		if(m_varRatio_null.n_elem > 1){	
			hasVarRatio = assignVarianceRatio(t_MAC_vec(j), false, false);	
		}
		scoreTestFast(t_GVec0, iIndex, t_Beta_val, t_seBeta_val, t_pval_str_val, t_altFreq0_val, t_Tstat_val, t_var1_val, t_var2_val, t_pval_val);
		//scoreTest(t_GVec, t_Beta_val, t_seBeta_val, t_pval_str_val, t_altFreq0_val, t_Tstat_val, t_var1_val, t_var2_val, t_gtildeVec, t_P2Vec, t_gy, false, iIndex, t_pval_val);


		t_Beta_vec(j) = t_Beta_val;
		t_seBeta_vec(j) = t_seBeta_val;
		t_Beta_vec(j) = t_Tstat_val;
		t_var1_vec(j) = t_var1_val;
		t_pval_vec(j) = t_pval_val; //convert from string
		t_pval_noadj_vec(j) = t_pval_val;
		std::cout << "t_pval_val " << t_pval_val << std::endl;
	}
   }


*/


  arma::vec StdStat_vec = arma::abs(t_Tstat_vec) / arma::sqrt(t_var1_vec);
  if(m_traitType != "quantitative"){
  double q, qinv, m1, NAmu, NAsigma, tol1, p_iIndexComVecSize;
  arma::uvec iIndexComVec, iIndex;
  bool t_isSPAConverge;
  arma::vec gNB, gNA, muNB, muNA, t_gtilde, t_GVec;
  double gmuNB, t_SPApval, t_Beta_val, t_seBeta_val, t_altFreq0_val, t_Tstat_val, t_var1_val, t_var2_val, t_pval_val;
  std::string t_pval_str_val;
  unsigned int iIndexComVecSize, iIndexSize;
  arma::mat t_GMat_sub, t_GMat_tilde_sub;
  //, t_GMat0_sub, t_GMat_tilde_sub;
  bool hasVarRatio;	
	t_SPAind_vec = arma::find(StdStat_vec > m_SPA_Cutoff);
	if(t_SPAind_vec.n_elem > 0){
		t_GMat_sub = t_GMat.cols(t_SPAind_vec);
		//if((g_I_longl_mat.n_cols != g_I_longl_mat.n_rows) && (t_GMat_sub.n_rows < m_y.n_elem)){
                //	t_GMat0_sub = g_I_longl_mat * t_GMat_sub;
        	//}else{
                //	t_GMat0_sub = t_GMat_sub;
        	//}		
		//getadjGFast_multi(t_GMat0_sub, t_GMat_tilde_sub);
		//t_GMat_sub = t_GMat_center.cols(t_SPAind_vec);
		//t_GMat_sub = t_GMat_center.cols(t_SPAind_vec);
		getadjGFast_multi(t_GMat_sub, t_GMat_tilde_sub);
		//
		//if(m_varRatio_null.n_elem == 1){
		//	assignSingleVarianceRatio(false, false);
		//}
	}
    //std::cout << "hererer" << std::endl; 
    //std::cout << "t_SPAind_vec.n_elem " << t_SPAind_vec.n_elem << std::endl; 
    for(int i = 0; i < t_SPAind_vec.n_elem; i++){
      unsigned int j = t_SPAind_vec(i);

     if(m_varRatio_null.n_elem > 1){
     	hasVarRatio = assignVarianceRatio(t_MAC_vec(j), m_flagSparseGRM_cur, true);
      }
      t_gtilde = t_GMat_tilde_sub.col(i);
      //t_gtilde = t_GMat_sub.col(i);
      t_GVec = t_GMat_sub.col(i);
      //m1 = dot(m_mu_sample, t_GVec);
      iIndex = arma::find(t_GVec != 0);
      iIndexComVec = arma::find(t_GVec == 0);
      iIndexComVecSize = iIndexComVec.n_elem;
      iIndexSize = iIndex.n_elem;
      gNB.set_size(iIndexSize);
      gNA.set_size(iIndexComVecSize);
      muNB.set_size(iIndexSize);
      muNA.set_size(iIndexComVecSize);

      m1 = dot(m_mu, t_gtilde);
      p_iIndexComVecSize = double(iIndexComVecSize)/m_n;

     //if(p_iIndexComVecSize >= 0.5 && !m_flagSparseGRM_cur){
      gNB = t_gtilde(iIndex);
      gNA = t_gtilde(iIndexComVec);
      muNB = m_mu(iIndex);
      muNA = m_mu(iIndexComVec);
      gmuNB = dot(gNB,muNB);
      NAmu= m1-gmuNB;
     //} 



        if(m_traitType == "binary"){
                q = (t_Tstat_vec(j))/sqrt((t_var1_vec(j))/(t_var2_vec(j))) + m1;

                if((q-m1) > 0){
                        qinv = -1 * std::abs(q-m1) + m1;
                }else if ((q-m1) == 0){
                        qinv =  m1;
                }else{
                        qinv = std::abs(q-m1) + m1;
                }
                //if(p_iIndexComVecSize >= 0.5 && !m_flagSparseGRM_cur){
                        NAsigma = t_var2_vec(j) - arma::sum(muNB % (1-muNB) % arma::pow(gNB,2));
                //}
		
        }else if(m_traitType == "survival" || m_traitType == "count"){
		q = (t_Tstat_vec(j))/sqrt((t_var1_vec(j))/(t_var2_vec(j)));
                qinv = -q;
		//if(p_iIndexComVecSize >= 0.5 && !m_flagSparseGRM_cur){
                        NAsigma = t_var2_vec(j) - arma::sum(muNB % arma::pow(gNB,2));
                //}
        }
        bool logp=false;
        double tol0 = std::numeric_limits<double>::epsilon();
        double pval_noadj = t_pval_noadj_vec(j); 
        tol1 = std::pow(tol0, 0.25);



	//if(p_iIndexComVecSize >= 0.5 && !m_flagSparseGRM_cur){
        //        std::cout << "SPA_fast" << std::endl;
        SPA_fast(m_mu, t_gtilde, q, qinv, pval_noadj, logp, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1, m_traitType, t_SPApval, t_isSPAConverge);

  //std::cout << "q " << q << std::endl;
        //m_mu.print("m_mu");
	//
/*
	        std::cout << "iIndex.n_elem " << iIndex.n_elem << std::endl;
	arma::uvec v = { 0, 1, 2, 3, 4 };
        //t_gtilde.print("t_gtilde");
	if(pval_noadj < 0.0001){
	arma::vec muvec = m_mu.elem(v);
	muvec.print("muvec");
	arma::vec gtilevec = t_gtilde.elem(v);
	gtilevec.print("gtilevec");
	std::cout << "NAmu " << NAmu << " NAsigma " << NAsigma << " m1 " << m1 << " gmuNB " << gmuNB << std::endl;
        std::cout << "pval_noadj " << pval_noadj << std::endl;
        std::cout << "t_SPApval " << t_SPApval << std::endl;
	}
*/	
	//}else{
 //               std::cout << "SPA" << std::endl;
        //        SPA(m_mu, t_GVec, q, qinv, pval_noadj, tol1, logp, m_traitType, t_SPApval, t_isSPAConverge);
        //}


        boost::math::normal ns;
        t_pval_val = t_SPApval;
        double t_qval;
        t_Beta_val = t_Beta_vec(j);
        try {
           t_qval = boost::math::quantile(ns, t_pval_val/2);
           t_qval = fabs(t_qval);
           t_seBeta_val = fabs(t_Beta_val)/t_qval;
        }catch (const std::overflow_error&) {
          t_qval = std::numeric_limits<double>::infinity();
          t_seBeta_val = 0;
        }

        t_seBeta_vec(j) = t_seBeta_val;
        t_pval_vec(j) = t_pval_val; 

	}
  }
  
  t_Beta = arma::conv_to < std::vector<double> >::from(t_Beta_vec);
  t_seBeta = arma::conv_to < std::vector<double> >::from(t_seBeta_vec);
  t_Tstat = arma::conv_to < std::vector<double> >::from(t_Tstat_vec);
  t_var1 = arma::conv_to < std::vector<double> >::from(t_var1_vec);
  t_pval = arma::conv_to < std::vector<double> >::from(t_pval_vec);
  t_pval_noSPA = arma::conv_to < std::vector<double> >::from(t_pval_noadj_vec);
}





void SAIGEClass::scoreTestFast_noadjCov_multi(arma::mat & t_GMat,
		     arma::vec & t_MAFvec,
                     arma::vec & t_Beta,
                     arma::vec & t_seBeta,
                     std::vector <std::string> & t_pval_str,
                     arma::vec &t_Tstat,
                     arma::vec &t_var1,
                     arma::vec &t_var2,
                     arma::vec & t_pval_noadj_vec, 
		     arma::vec & t_varRatio_multi_vec)
{


arma::mat GMat_centered;
/*
getadjGFast_multi(t_GMat, GMat_tilde);

std::cout << "xxxxx" << std::endl;
unsigned int x = t_GMat.n_cols;
std::cout << "x " << x  << std::endl;
x = t_GMat.n_rows;
std::cout << "x " << x  << std::endl;
x = GMat_tilde.n_cols;
std::cout << "x " << x  << std::endl;
x = GMat_tilde.n_rows;
std::cout << "x " << x  << std::endl;
*/

//t_Tstat = GMat_tilde.t() * m_res_sample / (m_tauvec[0]);
//t_Tstat = t_GMat.t() * m_res / (m_tauvec[0]);
GMat_centered = t_GMat.each_row() - t_MAFvec.t();
t_Tstat = GMat_centered.t() * m_res / (m_tauvec[0]);
//arma::mat mu2Mat = arma::diag(m_mu2_sample);
arma::mat GMat_centered_weight = (GMat_centered.each_col()) % m_mu2;
GMat_centered_weight = GMat_centered_weight % GMat_centered;
t_var2 = (arma::sum(GMat_centered_weight, 0).t());

if(t_varRatio_multi_vec.n_elem == 0){
	t_var1 = t_var2 * m_varRatioVal;
}else{
	t_var1 = t_var2 % t_varRatio_multi_vec;
}

arma::vec stat_vec = arma::pow(t_Tstat, 2) / t_var1;
double pval, stat;

int numMarkers = t_var1.n_elem;

t_pval_noadj_vec.set_size(numMarkers);

for(int i = 0; i < numMarkers; i++){
        if (t_var1(i) <= std::numeric_limits<double>::min()){
                pval = 1;
                stat = 0;
        }else{
                stat = stat_vec(i);
                pval = Rf_pchisq(stat, 1.0, 0, 0);
        }
        char pValueBuf[100];
        if (pval != 0){
                sprintf(pValueBuf, "%.6E", pval);
        }else {
                double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
                int exponent = floor(log10p);
                double fraction = pow(10.0, log10p - exponent);
                if (fraction >= 9.95) {
                        fraction = 1;
                         exponent++;
                }
                sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
        }
        std::string buffAsStdStr = pValueBuf;
        t_pval_str.push_back(buffAsStdStr);
        t_pval_noadj_vec(i) = pval;
}
t_Beta = t_Tstat/t_var1;
t_seBeta = arma::abs(t_Beta) / sqrt(arma::abs(stat_vec));
}



void SAIGEClass::getadjGFast_multi(arma::mat & t_GMat, arma::mat & t_GMat_tilde)
{
 unsigned int z = m_XV.n_cols;
 std::cout << "z " << z << std::endl;
 z = m_XV.n_rows;
 std::cout << "z " << z << std::endl;


   if(g_I_longl_mat.n_cols != g_I_longl_mat.n_rows && t_GMat.n_rows != m_XV.n_cols){
      t_GMat_tilde = g_I_longl_mat * t_GMat;
   }else{
      t_GMat_tilde = t_GMat;	

   }
  t_GMat_tilde = t_GMat_tilde - m_XXVX_inv * m_XV * t_GMat_tilde;
}



void SAIGEClass::scoreTest_multi(arma::mat & t_GMat_tilde,
                     arma::vec & t_Beta_vec,
                     arma::vec & t_seBeta_vec,
                     std::vector <std::string> & t_pval_str_vec,
                     arma::vec & t_pval_vec,
                     arma::vec & t_Tstat_vec,
                     arma::vec & t_var1_vec,
                     arma::vec & t_var2_vec,
                     arma::mat & t_P2Mat,
                     arma::vec & t_gy_vec,
		     arma::vec & t_varRatio_vec,
                     bool t_is_region){
    arma::vec Sm_vec, var2m_vec, S_vec;
    double var2;

    if(t_is_region && m_traitType == "binary"){
      t_gy_vec = t_GMat_tilde.t() * m_y;
    }

    arma::vec weighted_res_vec = m_res % m_varWeightsvec;

    S_vec = t_GMat_tilde.t() * weighted_res_vec;

    S_vec = S_vec/m_tauvec[0];

    arma::mat t_GMat_tilde_weight;


        double varRatioVal_var2 = m_varRatioVal;



    if(!m_flagSparseGRM_cur){
      t_P2Mat = t_GMat_tilde.each_col() % m_mu2 * m_tauvec[0];
      t_GMat_tilde_weight = t_P2Mat % t_GMat_tilde;
      var2m_vec = (arma::sum(t_GMat_tilde_weight, 0).t());
    }else{
      if(m_SigmaMat_sp.n_rows > 2){
      //t_P2Vec = arma::spsolve(m_SigmaMat_sp, t_gtilde);
      t_P2Mat = t_GMat_tilde.t() * m_SigmaMat_sp;
      t_GMat_tilde_weight = t_P2Mat % t_GMat_tilde;
      var2m_vec = (arma::sum(t_GMat_tilde_weight, 0).t());
      if(m_isVarPsadj){
        varRatioVal_var2 = 1;
        var2m_vec = var2m_vec - t_GMat_tilde.t() * m_Sigma_iXXSigma_iX * m_X.t() * t_P2Mat;
      }
     }else{
        //t_P2Vec = m_sigmainvG_noV;
	for(int i = 0; i < t_GMat_tilde.n_cols; i++){
	arma::vec t_GVec_tilde=t_GMat_tilde.col(i);
          t_P2Mat.col(i) = getSigma_G_V(t_GVec_tilde, 500, 1e-5);
	}
	t_GMat_tilde_weight = t_P2Mat % t_GMat_tilde;
	var2m_vec = (arma::sum(t_GMat_tilde_weight, 0).t());
      if(m_isVarPsadj){
        varRatioVal_var2 = 1;
	var2m_vec = var2m_vec - t_GMat_tilde.t() * m_Sigma_iXXSigma_iX * m_X.t() * t_P2Mat;
      }
     }
    }


    //std::cout << "var2 " << var2 << std::endl;
    //double var1 = var2 * m_varRatioVal;
	if(t_varRatio_vec.n_elem == 0){
        	t_var1_vec = var2m_vec * m_varRatioVal;
	}else{
        	t_var1_vec = var2m_vec % t_varRatio_vec;
	}

    arma::vec stat_vec = arma::pow(S_vec, 2)/t_var1_vec;
    double t_pval, var1, stat;
    t_pval_vec.set_size(t_GMat_tilde.n_cols);


  for(int i = 0; i < t_GMat_tilde.n_cols; i++){
    var1 = t_var1_vec(i);
    stat = stat_vec(i);
    if (var1 <= std::numeric_limits<double>::min()){
        t_pval = 1;
    }else{
        boost::math::chi_squared chisq_dist(1);
        t_pval = boost::math::cdf(complement(chisq_dist, stat));
    }

    char pValueBuf[100];
    if (t_pval != 0)
        sprintf(pValueBuf, "%.6E", t_pval);
    else {
        double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    }

    std::string buffAsStdStr = pValueBuf;
    t_pval_str_vec.push_back(buffAsStdStr);
  }  
    t_Beta_vec = S_vec/t_var1_vec;
    t_seBeta_vec = arma::abs(t_Beta_vec) / arma::sqrt(arma::abs(stat_vec));
    t_Tstat_vec = S_vec;
}


void SAIGEClass::assign_for_trait_i_2(unsigned int itrait){
	//std::cout << "m_n " << m_n << std::endl;
        unsigned int startip = itrait*m_p;
        unsigned int endip = startip + m_p - 1;
        unsigned int startin = itrait*m_n;
        unsigned int endin = startin + m_n - 1;

	//std::cout << startip << " " << endip << " " << startin << " " << endin;
        m_y = m_y_mt.col(itrait);
        //m_offset = m_offset_mt.col(itrait);

        //m_X = m_X_mt.rows(startin, endin);
        //m_S_a = m_S_a_mt.col(itrait);
        m_mu = m_mu_mt.col(itrait);
        //m_resout = m_resout_mt.col(itrait);

        //m_XVX = m_XVX_mt.rows(startip, endip);
        m_XV = m_XV_mt.rows(startip, endip);
        m_XXVX_inv = m_XXVX_inv_mt.rows(startin, endin);
        //m_XVX_inv_XV = m_XVX_inv_XV_mt.rows(startin, endin);
}



void SAIGEClass::assign_for_trait_i(unsigned int itrait){
        //std::cout << "m_n " << m_n << std::endl;
        unsigned int startip = itrait*m_p;
        unsigned int endip = startip + m_p - 1;
        unsigned int startin = itrait*m_n;
        unsigned int endin = startin + m_n - 1;

        //std::cout << startip << " " << endip << " " << startin << " " << endin;


        //m_XVX = m_XVX_mt.rows(startip, endip);
        //m_XV = m_XV_mt.rows(startip, endip);
        //m_XXVX_inv = m_XXVX_inv_mt.rows(startin, endin);
        //m_XVX_inv_XV = m_XVX_inv_XV_mt.rows(startin, endin);
        //m_Sigma_iXXSigma_iX = m_Sigma_iXXSigma_iX_mt.rows(startin, endin); //check
        m_res = m_res_mt.col(itrait);
        m_mu2 = m_mu2_mt.col(itrait);
        m_varRatio_sparse = m_varRatio_sparse_mt.col(itrait);
        m_varRatio_null = m_varRatio_null_mt.col(itrait);
        m_varRatio_null_noXadj = m_varRatio_null_noXadj_mt.col(itrait);
        m_varRatio_null_eg = m_varRatio_null_eg_mt.col(itrait);
        m_varRatio_sparse_eg = m_varRatio_sparse_eg_mt.col(itrait);
        m_tauvec = m_tauvec_mt.col(itrait);
        m_varWeightsvec = m_varWeightsvec_mt.col(itrait);
        m_traitType = m_traitType_vec.at(itrait);
}


void SAIGEClass::assign_for_itrait(unsigned int t_itrait){
	m_itrait = t_itrait;
	m_traitType = m_traitType_vec.at(m_itrait);
	m_startip = m_itrait*m_p;
	m_endip = m_startip + m_p - 1;
	//std::cout << "assign_for_itrait m_is_gxe " << m_is_gxe << std::endl;

	if(m_is_gxe){
	//std::cout << "m_XV_gxe_mt.n_rows " << m_XV_gxe_mt.n_rows << std::endl;
	//std::cout << "m_itrait " << m_itrait << std::endl;
	m_p_gxe = (m_XV_gxe_mt.n_rows) / (m_traitType_vec.size());
	//std::cout << "m_p_gxe " << m_p_gxe << std::endl;
	m_startip_gxe = m_itrait*m_p_gxe;
	m_endip_gxe = m_startip_gxe + m_p_gxe - 1;
	m_n_gxe = (m_XV_gxe_mt.n_cols) / (m_traitType_vec.size());
	m_startin_gxe = m_itrait*m_n_gxe;
	m_endin_gxe = m_startin_gxe + m_n_gxe - 1;
	}
	m_startin = m_itrait*m_n;
	m_endin = m_startin + m_n - 1;
	

	m_startic = m_itrait*m_numMarker_cond;
	m_endic = m_startic + m_numMarker_cond - 1;

}	


void SAIGEClass::assign_for_itrait_binaryindices(unsigned int t_itrait){

        m_case_indices = arma::find(m_y_mt.col(t_itrait) == 1);
        m_ctrl_indices = arma::find(m_y_mt.col(t_itrait) == 0);
        m_n_case = m_case_indices.n_elem;
        m_n_ctrl = m_ctrl_indices.n_elem;
}



// revised
// need to add sparse Sigma version
// This function only uses variance ratio and does not use sparse GRM
void SAIGEClass::getMarkerPval_gxe(arma::vec & t_GVec,
                               arma::uvec & iIndex,
                               arma::uvec & iIndexComVec,
                               double& t_Beta,
                               double& t_seBeta,
                               double& t_pval,
                               double& t_pval_noSPA,
                               double t_altFreq,
                               double& t_Tstat,
                               double& t_gy,
                               double& t_var1,
                               bool & t_isSPAConverge,
                               arma::vec & t_gtilde,
                               bool & is_gtilde,
                               bool  is_region,
                               arma::vec & t_P2Vec,
                               bool t_isCondition,
                               double& t_Beta_c,
                                double& t_seBeta_c,
                                double& t_pval_c,
                                double& t_pval_noSPA_c,
                                double& t_Tstat_c,
                                double& t_varT_c,
                                arma::rowvec & t_G1tilde_P_G2tilde,
                                bool & t_isFirth,
                                bool & t_isFirthConverge,
                                bool t_isER,
                                bool t_isnoadjCov,
                                bool t_isSparseGRM){

  std::string t_pval_str;
  double t_var2, t_SPApval;
  double altFreq0;
  arma::vec t_GVec0 = t_GVec;
  arma::uvec indexNonZeroVec0_arma = iIndex;
  arma::uvec indexZeroVec0_arma = iIndexComVec;


	scoreTest_gxe(t_GVec0, t_Beta, t_seBeta, t_pval_str, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde, t_P2Vec, t_gy, is_region, indexNonZeroVec0_arma, t_pval);

        //scoreTestFast_gxe(t_GVec0, indexNonZeroVec0_arma, t_Beta, t_seBeta, t_pval_str, altFreq0, t_Tstat, t_var1, t_var2, t_pval);


 is_gtilde = true;

 if(!is_gtilde){
        getadjGFast_gxe(t_GVec, t_gtilde, indexNonZeroVec0_arma);
   }


//std::cout << "here is_gtilde" << std::endl;

//std::cout << "t_pval " << t_pval << std::endl;


        t_G1tilde_P_G2tilde = sqrt(m_varRatioVal) * t_gtilde.t() * (m_P2Mat_cond.cols(m_startic, m_endic));
        arma::vec t_Tstat_ctemp =  t_G1tilde_P_G2tilde * (m_VarInvMat_cond.cols(m_startic, m_endic)) * (m_Tstat_cond.subvec(m_startic, m_endic));
//std::cout << "here is_gtilde 2" << std::endl;
        arma::mat tempgP2 = t_gtilde.t() * (m_P2Mat_cond.cols(m_startic, m_endic));
//std::cout << "here is_gtilde 3" << std::endl;
	//std::cout << "t_Tstat " << t_Tstat << std::endl;
	//std::cout << "t_Tstat_ctemp(0) " << t_Tstat_ctemp(0) << std::endl;
	//m_Tstat_cond.print("m_Tstat_cond");
	//t_G1tilde_P_G2tilde.print("t_G1tilde_P_G2tilde");
	//m_VarInvMat_cond.print("m_VarInvMat_cond");
        t_Tstat_c = t_Tstat - t_Tstat_ctemp(0);
//std::cout << "here is_gtilde 4" << std::endl;
        arma::vec t_varT_ctemp = t_G1tilde_P_G2tilde * (m_VarInvMat_cond.cols(m_startic, m_endic)) * (t_G1tilde_P_G2tilde.t());
//std::cout << "here is_gtilde 5" << std::endl;
        t_varT_c = t_var1 - t_varT_ctemp(0);
        //std::cout << "t_varT_c " << t_varT_c << " t_var1 " << t_var1 << " t_varT_ctemp(0) " << t_varT_ctemp(0) << std::endl;
//	std::cout << "t_Tstat_c " << t_Tstat_c << " t_Tstat " << t_Tstat << " t_Tstat_ctemp(0) " << t_Tstat_ctemp(0) << std::endl;
//std::cout << "here is_gtilde 6" << std::endl;

    double S_c = t_Tstat_c;
    double stat_c = S_c*S_c/t_varT_c;
     //std::cout << "here is_gtilde 7" << std::endl;
     //std::cout << "t_varT_c " << t_varT_c << std::endl; 		
     //std::cout << "t_Tstat_c " << t_Tstat_c << std::endl; 		


     //std::cout << "t_var1 " << t_var1 << std::endl; 		
     //std::cout << "t_varT_ctemp(0) " << t_varT_ctemp(0) << std::endl; 	
     //std::cout << "std::numeric_limits<double>::min() " << std::numeric_limits<double>::min() << std::endl;
	double pval_noSPA_c; 
     //std::cout << "std::pow(std::numeric_limits<double>::min()) " << std::numeric_limits<double>::min() << std::endl;
     //std::cout << "t_varT_c " << t_varT_c << std::endl;
     if (t_varT_c <= std::numeric_limits<double>::min() || t_varT_c <= std::pow(10, -5)){
     //if (t_varT_c <= 10^-10){
        t_pval_noSPA_c = 1;
        stat_c = 0;
	t_Beta_c = 0;
	t_seBeta_c = fabs(t_Beta_c) / sqrt(stat_c);
	t_Tstat_c = S_c;
	std::string buffAsStdStr_c = "1";
	std::string& t_pval_noSPA_str_c =  buffAsStdStr_c;
	pval_noSPA_c = 1;
     }else{
        boost::math::chi_squared chisq_dist(1);
        t_pval_noSPA_c = boost::math::cdf(complement(chisq_dist, stat_c));

    //std::cout << "t_pval_noSPA_c " << t_pval_noSPA_c << std::endl;
    char pValueBuf_c[100];
    if (t_pval_noSPA_c != 0)
        sprintf(pValueBuf_c, "%.6E", t_pval_noSPA_c);
    else {
        double log10p_c = log10(2.0) - M_LOG10E*stat_c/2 - 0.5*log10(stat_c*2*M_PI);
        int exponent_c = floor(log10p_c);
        double fraction_c = pow(10.0, log10p_c - exponent_c);
        if (fraction_c >= 9.95) {
          fraction_c = 1;
           exponent_c++;
         }
        sprintf(pValueBuf_c, "%.1fE%d", fraction_c, exponent_c);
    }


    std::string buffAsStdStr_c = pValueBuf_c;
    std::string& t_pval_noSPA_str_c = buffAsStdStr_c;
//std::cout << "here is_gtilde 7" << std::endl;

    t_Beta_c = S_c/t_varT_c;
    t_seBeta_c = fabs(t_Beta_c) / sqrt(stat_c);
    t_Tstat_c = S_c;
  //}
  //double pval_noSPA_c;
  try {
        pval_noSPA_c = std::stod(t_pval_noSPA_str_c);
  } catch (const std::invalid_argument&) {
        pval_noSPA_c = 0;
        std::cerr << "Argument is invalid\n";
        //throw;
  } catch (const std::out_of_range&) {
        std::cerr << "Argument is out of range for a double\n";
        //throw;
        pval_noSPA_c = 0;
  }
  
  
  t_pval_noSPA_c = pval_noSPA_c;
}
  //std::cout << "here is_gtilde 8" << std::endl;
  //std::cout << "pval_noSPA_c " << pval_noSPA_c << std::endl;

double q, qinv, m1, NAmu, NAsigma, tol1, p_iIndexComVecSize;
arma::vec gNB, gNA, muNB, muNA;
  double gmuNB;;

    if((m_traitType == "binary" || m_traitType == "count")  && pval_noSPA_c != 1 && stat_c > std::pow(m_SPA_Cutoff,2)){
/*
  	arma::vec mu_vec = m_mu_gxe_mt.col(m_itrait);
         m1 = dot(mu_vec, t_gtilde);
        bool t_isSPAConverge_c;
        double q_c, qinv_c, pval_noadj_c, SPApval_c;


        if(p_iIndexComVecSize >= 0.5 && !m_flagSparseGRM_cur ){
                unsigned int j1 = 0;
                unsigned int j2 = 0;
                gNB = t_gtilde(iIndex);
                gNA = t_gtilde(iIndexComVec);
                muNB = mu_vec(iIndex);
                //muNB = m_mu_mt(iIndex, m_itrait);
                muNA = mu_vec(iIndexComVec);
                //muNA = m_mu_mt(iIndexComVec, m_itrait);
                gmuNB = dot(gNB,muNB);
                NAmu= m1-gmuNB;
        }else{
                gNA.clear();
                gNB.clear();
                muNA.clear();
                gNB.clear();

        }

        if(m_traitType == "binary"){
                q_c = t_Tstat_c/sqrt(t_varT_c/t_var2) + m1;

                if((q_c-m1) > 0){
                        qinv_c = -1 * std::abs(q_c-m1) + m1;
                }else if ((q_c-m1) == 0){
                        qinv_c =  m1;
                }else{
                        qinv_c = std::abs(q_c-m1) + m1;
                }
        }else if(m_traitType == "count"){
                q_c = t_Tstat_c/sqrt(t_varT_c/t_var2);
                qinv = -q_c;
		NAsigma = t_var2 - arma::sum(muNB % arma::pow(gNB,2));
        }

        bool logp=false;

std::cout << "here is_gtilde 9" << std::endl;


        if(p_iIndexComVecSize >= 0.5 && !m_flagSparseGRM_cur){
                //std::cout << "SPA_fast " << std::endl;
                SPA_fast(mu_vec, t_gtilde, q_c, qinv_c, pval_noadj_c, false, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1, m_traitType, SPApval_c, t_isSPAConverge_c);
        }else{
                std::cout << "SPA " << std::endl;
                SPA(mu_vec, t_gtilde, q_c, qinv_c, pval_noadj_c, tol1, logp, m_traitType, SPApval_c, t_isSPAConverge_c);
        }
*/
	//std::cout << "here is_gtilde 10" << std::endl;
        boost::math::normal ns;
        t_pval_c =pval_noSPA_c;
        //t_pval_c = SPApval_c;
        double t_qval_c;
        try {
           t_qval_c = boost::math::quantile(ns, t_pval_c/2);
           t_qval_c = fabs(t_qval_c);
           t_seBeta_c = fabs(t_Beta_c)/t_qval_c;
        }catch (const std::overflow_error&) {
          t_qval_c = std::numeric_limits<double>::infinity();
          t_seBeta_c = 0;
        }
    }else{
     t_pval_c = t_pval_noSPA_c;
    
    }
   }




void SAIGEClass::getadjGFast_gxe(arma::vec & t_GVec, arma::vec & g, arma::uvec & iIndex)
{

  //arma::vec t_GVec0 = t_GVec;
  arma::uvec indexNonZeroVec0_arma = iIndex;
  // To increase computational efficiency when lots of GVec elements are 0
  arma::vec m_XVG_gxe(m_p_gxe, arma::fill::zeros);
  //arma::mat m_XV_gxe_mt_sub(m_p_gxe, m_XV_gxe_mt.n_cols);
  for(int i = 0; i < indexNonZeroVec0_arma.n_elem; i++){
  	//m_XV_gxe_mt_sub = m_XV_gxe_mt.rows(m_startip_gxe, m_endip_gxe);
	//m_XVG_gxe += m_XV_gxe_mt_sub.col(indexNonZeroVec0_arma(i)) * t_GVec0(indexNonZeroVec0_arma(i)); 
        m_XVG_gxe += m_XV_gxe_mt.rows(m_startip_gxe, m_endip_gxe).col(indexNonZeroVec0_arma(i)) * t_GVec(indexNonZeroVec0_arma(i));
  }
  g = t_GVec - m_XXVX_inv_gxe_mt.rows(m_startin_gxe, m_endin_gxe) * m_XVG_gxe;
}

void SAIGEClass::scoreTest_gxe(arma::vec & t_GVec,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2,
                     arma::vec & t_gtilde,
                     arma::vec & t_P2Vec,
                     double& t_gy,
                     bool t_is_region,
                     arma::uvec & t_indexForNonZero,
                     double & t_pval){
    arma::vec Sm, var2m;
    double S, var2;
    getadjGFast_gxe(t_GVec, t_gtilde, t_indexForNonZero);
    //getadjG(t_GVec, t_gtilde);


    if(t_is_region && m_traitType == "binary"){
      t_gy = dot(t_gtilde, m_y_gxe_mt.col(m_itrait));
    }

    S = dot(t_gtilde, m_res_gxe_mt.col(m_itrait) % m_varWeights_gxe_mt.col(m_itrait));

    //std::cout << "S scoreTest_gxe " << S << std::endl;

    S = S/m_tauvec_mt(0,m_itrait);

    double varRatioVal_var2 = m_varRatioVal;

   if(!m_flagSparseGRM_cur){
      t_P2Vec = t_gtilde % (m_mu2_gxe_mt.col(m_itrait)) *(m_tauvec_mt(0,m_itrait));
      //t_P2Vec = t_gtilde % m_mu2;
      var2m = dot(t_P2Vec, t_gtilde);
    }else{
      if(m_SigmaMat_sp.n_rows > 2){
      //t_P2Vec = arma::spsolve(m_SigmaMat_sp, t_gtilde);
      t_P2Vec = m_SigmaMat_sp * t_gtilde;
      var2m = dot(t_P2Vec, t_gtilde);
      if(m_isVarPsadj){
        varRatioVal_var2 = 1;
        var2m = var2m - t_gtilde.t() * m_Sigma_iXXSigma_iX * m_X.t() * t_P2Vec;
      }
     }else{
        //t_P2Vec = m_sigmainvG_noV;
        t_P2Vec = getSigma_G_V(t_gtilde, 500, 1e-5);
        var2m = dot(t_P2Vec, t_gtilde);
      if(m_isVarPsadj){
        varRatioVal_var2 = 1;
        var2m = var2m - t_gtilde.t() * m_Sigma_iXXSigma_iX * m_X.t() * t_P2Vec;
      }
     }
    }

    var2 = var2m(0,0);
    //std::cout << "var2 Scoretest_gxe " << var2 << std::endl;
    //std::cout << "m_varRatioVal " << m_varRatioVal << std::endl;
    //double var1 = var2 * m_varRatioVal;
    double var1 = var2 * varRatioVal_var2;
    //std::cout << "var1 Scoretest_gxe " << var1 << std::endl;
    double stat = S*S/var1;
    //double t_pval;
    //std::cout << "S " << S << std::endl;
    //std::cout << "var1 " << var1 << std::endl;

    //if (var1 <= std::pow(std::numeric_limits<double>::min(), 2)){
    if (var1 <= std::numeric_limits<double>::min()){
        t_pval = 1;
    }else{
        boost::math::chi_squared chisq_dist(1);
        t_pval = boost::math::cdf(complement(chisq_dist, stat));
    }
    //std::cout << "t_pval Scoretest_gxe " << t_pval << std::endl;
/*
    char pValueBuf[100];
    if (t_pval != 0)
        sprintf(pValueBuf, "%.6E", t_pval);
    else {
        double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
*/
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
}


void SAIGEClass::scoreTestFast_gxe(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2,
                     double & t_pval){

    arma::vec g1 = t_GVec.elem(t_indexForNonZero);
    arma::vec m_varWeightsvec_new_sub = m_varWeights_gxe_mt.col(m_itrait);
    arma::vec m_varWeightsvec_new = m_varWeightsvec_new_sub.elem(t_indexForNonZero);
//    std::cout << "m_startin_gxe " << m_startin_gxe << " " << m_endin_gxe << std::endl;
//    std::cout << "m_startip_gxe " << m_startip_gxe << " " << m_endip_gxe << std::endl;
//    std::cout << "m_X_gxe_mt.rows " << m_X_gxe_mt.n_rows << " " <<  m_X_gxe_mt.n_cols <<  std::endl;
    arma::mat X_gxe_mt = m_X_gxe_mt.rows(m_startin_gxe, m_endin_gxe);
    arma::mat X1 = X_gxe_mt.rows(t_indexForNonZero);
    //arma::mat X1 = (m_X_gxe_mt.rows(m_startip_gxe, m_endip_gxe)).cols(t_indexForNonZero);
//    std::cout << "m_XVX_inv_XV_gxe_mt.rows " << m_XVX_inv_XV_gxe_mt.n_rows << " " <<  m_XVX_inv_XV_gxe_mt.n_cols <<  std::endl;
    arma::mat XVX_inv_XV_gxe_mt = m_XVX_inv_XV_gxe_mt.rows(m_startin_gxe, m_endin_gxe);
    arma::mat A1 = XVX_inv_XV_gxe_mt.rows(t_indexForNonZero);
    //arma::mat A1 = (m_XVX_inv_XV_gxe_mt.rows(m_startip_gxe, m_endip_gxe)).cols(t_indexForNonZero);
    arma::vec mu21;
    arma::vec res_gxe_vec = m_res_gxe_mt.col(m_itrait);
    arma::vec res1 = res_gxe_vec.elem(t_indexForNonZero);
    res1 = res1 % m_varWeightsvec_new;
    arma::vec Z = A1.t() * g1;
    arma::vec B = X1 * Z;
    arma::vec g1_tilde = g1 - B;
    double var1, var2, S, S1, S2, g1tildemu2;
    arma::vec S_a2;
    double Bmu2;
    arma::mat  ZtXVXZ = Z.t() * m_XVX_gxe_mt * Z;
//    std::cout << "OKKKKK\n";
    if(m_traitType == "binary" || m_traitType == "count"){

	arma::vec mu2_gxe_vec = m_mu2_gxe_mt.col(m_itrait);
      mu21  = mu2_gxe_vec.elem(t_indexForNonZero);
//      std::cout << "OKKKKKa\n";
      g1tildemu2 = dot(square(g1_tilde), mu21);
//      std::cout << "OKKKKKb\n";
      Bmu2 = arma::dot(square(B),  mu21);
//      std::cout << "OKKKKKc\n";
      var2 = ZtXVXZ(0,0) - Bmu2 + g1tildemu2;
//      std::cout << "OKKKKK1\n";
    }else if(m_traitType == "quantitative" || m_traitType == "count_nb"){
      Bmu2 = dot(g1, B % m_varWeightsvec_new);
      //Bmu2 = dot(g1, B);
      var2 = ZtXVXZ(0,0)*m_tauvec[0] +  dot(g1,g1 % m_varWeightsvec_new) - 2*Bmu2;
      //var2 = ZtXVXZ(0,0)*m_tauvec[0] +  dot(g1,g1) - 2*Bmu2;
    }

    var1 = var2 * m_varRatioVal;
  //    std::cout << "OKKKKK1\n";
    S1 = dot(res1, g1_tilde);
  //    std::cout << "OKKKKK1a\n";
    arma::mat res1X1_temp = (res1.t()) * X1;
  //    std::cout << "OKKKKK1b\n";
    arma::vec res1X1 = res1X1_temp.t();
  //    std::cout << "OKKKKK1c\n";
  //    m_S_a_gxe_mt.print("m_S_a_gxe_mt");
  //    res1X1.print("res1X1");
    S_a2 = m_S_a_gxe_mt.col(m_itrait) - res1X1;
    S2 = - arma::dot(S_a2,  Z);
  //    std::cout << "OKKKKK1d\n";
    S = S1 + S2;
    S = S/(m_tauvec_mt(0,m_itrait));
  //    std::cout << "OKKKKK2\n";
    double stat = S*S/var1;
    //double t_pval;
/*
    std::cout << "S FastTest " << S << std::endl;
    std::cout << "var1 " << var1 << std::endl;
    std::cout << "m_varRatioVal " << m_varRatioVal << std::endl;
*/

    //if (var1 <= std::pow(std::numeric_limits<double>::min(), 2)){
    if (var1 <= std::numeric_limits<double>::min()){
        t_pval = 1;
    } else{
      boost::math::chi_squared chisq_dist(1);
      t_pval = boost::math::cdf(complement(chisq_dist, stat));
    }


//      std::cout << "OKKKKK3\n";
    char pValueBuf[100];
    if (t_pval != 0)
        sprintf(pValueBuf, "%.6E", t_pval);
    else {
        double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
    //std::cout << "t_pval_str scoreTestFast " << t_pval_str << std::endl;
    //std::cout << "end of scoreTestFast" << std::endl;
}



}
