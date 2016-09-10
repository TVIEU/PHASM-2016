#include <iostream>
#include<fstream>
#include<stdlib.h>
#include<math.h>
#include <sstream>
#include<cstdlib>
#include <iomanip>

using namespace std;

double tab_exp_int[4][1001];
double rotation_axis[365];

//*********************************************************************************************************
//******************** CONSTANTS, GLOBAL VARIABLES ********************************************************
//*********************************************************************************************************

double h_planck = 6.62607004e-34; //Planck constant
double c_light = 299792458; //speed of light
double k_boltzmann = 1.38064852e-23;

double I0 = 4.35*1e13; //intensity
double B_field; //strength of the magnetic field (G)
double psincl; //inclination of the rotation axis to the dipole axis
double u = 0.73; //limb-darkening coefficient 0.73
double temperature = 4500; //temperature at height 400-500 km (typical for the emission of the Sr I line)

double Ju = 1.; //total angular momentum of the upper level
double gu = 1 ; //Landé factor of the upper level
double Jl = 0.; //total angular momentum of the lower level
double A_ul = 2.01*1e8; //Einstein coefficient for spontaneous emission
double nu = c_light/(4607*1e-10); //frequency corresponding to the transition

double w_JuJl0 = 1; //polarizability
double w_JuJl2 = 1; //polarizability

double N_tot=1; //total density (does not matter since it simplifies when computing the fractional Stokes parameters X/I)

double B_lu = 3*c_light*c_light/(2*h_planck*nu*nu*nu)*A_ul; //Einstein coefficient for spontaneous absorption

double B_nu = 2*h_planck/pow(c_light,2)*pow(nu,3)*1/(exp(h_planck*nu/(k_boltzmann*temperature))-1); //Planck function

double factor = 10.1;

//*********************************************************************************************************
//******************** FUNCTIONS **************************************************************************
//*********************************************************************************************************

//compute arctan(sina/cosa) taking into account signs in order to get a result between 0 and 2*PI, not in [-PI/2,PI/2]
double arctan(double cosa , double sina)
{
    double a = atan(sina/cosa) ;
    if (cosa == 0 and sina > 0) a = M_PI/2;
    if (cosa == 0 and sina < 0) a = 3*M_PI/2;
    if (cosa < 0 and sina < 0) a += M_PI;
    if (cosa < 0 and sina > 0) a += M_PI;
    if (cosa > 0 and sina < 0) a += 2*M_PI;
    return(a);
}
//exponential integral of order n
double exp_int(double x,int n)
{
    double res = 0;
    double b=1000,pas=0.01;
    int i ; int ncompt = (b-1)/pas;
    for(i=0;i<=ncompt;i++)
        res += pas*exp(-x*(1+i*pas))/pow((1+i*pas),n);

    return (res);
}
void calcul_file_exp_int()
{
    int n;
    for(n=0;n<4;n++)
    {
        std::stringstream file;
        file << "exp_int/exp_int_" << n+2 << ".txt";
        fstream fun(file.str().c_str(),ios::out);

        int i;
        for(i=0;i<=1000;i++)
        {
            fun << (float)i/(float)1000 << " " << exp_int((double)i/(double)1000,n+2) << endl;
            cout << i<< endl;
        }
    }
}

void init_tab_exp_int()
{
    int n;
    for(n=0;n<4;n++)
    {
        std::stringstream file;
        file << "exp_int/exp_int_" << n+2 << ".txt";
        fstream fun(file.str().c_str(),ios::in);

        double abs;

        int i;
        for(i=0;i<=1000;i++)
            fun >> abs >> tab_exp_int[n][i];
    }
}

void init_tab_rotation()
{
    fstream data("rotation_axis.txt",ios::in);
    int n;double abs;
    for(n=0;n<366;n++) //fill the array of angles jncl from the data file
        data >> abs >> rotation_axis[n];
    data.close();
}

//provides the integration of a function on the projected disk
double integ_tilde(double (*f)(double, double, double, double, double),double small_gamma,double incl)
{
    int n = 30;//defines the precision
    double tab_values[n] = {-0.99689348, -0.98366812, -0.96002186, -0.92620005, -0.88256054,
        -0.82956576, -0.76777743, -0.69785049, -0.62052618, -0.53662415,
        -0.44703377, -0.35270473, -0.25463693, -0.15386991, -0.05147184,
         0.05147184,  0.15386991,  0.25463693,  0.35270473,  0.44703377,
         0.53662415,  0.62052618,  0.69785049,  0.76777743,  0.82956576,
         0.88256054,  0.92620005,  0.96002186,  0.98366812,  0.99689348};
     double tab_weight[n] = {0.00796819,  0.01846647,  0.02878471,  0.03879919,  0.04840267,
         0.05749316,  0.06597423,  0.07375597,  0.0807559 ,  0.08689979,
         0.09212252,  0.09636874,  0.09959342,  0.10176239,  0.10285265,
         0.10285265,  0.10176239,  0.09959342,  0.09636874,  0.09212252,
         0.08689979,  0.0807559 ,  0.07375597,  0.06597423,  0.05749316,
         0.04840267,  0.03879919,  0.02878471,  0.01846647,  0.00796819};

    double lambda,longi,eta,xi,ro,alpha;
    double s=0 ; int i_eta, i_xi;
    for(i_eta=0; i_eta<n; i_eta++)
    {
        for(i_xi=0; i_xi<n; i_xi++)
        {
            eta = tab_values[i_eta]; //eta,xi follow a Gaussian distribution in [-1,1]
            xi = tab_values[i_xi];
            ro = (xi+1)/2; //ro, alpha are polar coordinates
            alpha = (eta+1)*M_PI;
            lambda = fmod(acos(ro*cos(alpha)*sin(incl)-cos(incl)*sqrt(1-ro*ro)),(M_PI)); //lambda, longi are spherical coordinates
            longi = arctan((sin(incl)*sqrt(1-ro*ro)+cos(incl)*cos(alpha)*ro),ro*sin(alpha));
            s+=ro*M_PI/2*tab_weight[i_eta]*tab_weight[i_xi]*f(small_gamma,lambda,incl,longi,sqrt(1-ro*ro));
        }
    }
    return (s);
}

double hydro(double mu) //gives the hydrogen atoms density as a function of mu
{
    double tab_mu[9] = {0.9,0.55,0.34,0.2,0.1,0.065,0.04,0.02,0.01};
    double tab_hydro[9] = {2.33*1e16,1.54*1e16,1.01*1e16,6.54*1e15,4.19*1e15,2.90*1e15,2.08*1e15,1.48*1e15,9.89*1e14};

    int i;
    int ibest=0;
    for(i=1;i<=8;i++)
        if (fabs(mu-tab_mu[i]) <= fabs(mu-tab_mu[ibest])) ibest=i;

    return(tab_hydro[ibest]);
}
//reduced rotation matrix with j=2.
double reduce_rot_mat(double beta , int i , int j)
{
    switch(i)
    {
        case -2 :
            switch(j)
            {
                case -2 :
                    return(pow(cos(beta/2),4));
                    break;
                case -1 :
                    return(1/2.*sin(beta)*(1+cos(beta)));
                    break;
                case 0 :
                    return(sqrt(3/8.)*sin(beta)*sin(beta));
                    break;
                case 1 :
                    return(-1/2.*sin(beta)*(cos(beta)-1));
                    break;
                case 2 :
                    return(pow(sin(beta/2),4));
                    break;
                default:
                    cout << "Error to compute a coefficient of the reduces rotation matrix" << endl;
                    return 0;
            }
            break;
        case -1 :
            switch(j)
            {
                case -2 :
                    return(-1/2.*sin(beta)*(1+cos(beta)));
                    break;
                case -1 :
                    return(1/2.*(2*cos(beta)-1)*(cos(beta)+1));
                    break;
                case 0 :
                    return(sqrt(3/2.)*sin(beta)*cos(beta));
                    break;
                case 1 :
                    return(1/2.*(2.*cos(beta)+1)*(1-cos(beta)));
                    break;
                case 2 :
                    return(-1/2.*sin(beta)*(cos(beta)-1));
                    break;
                default:
                    cout << "Error to compute a coefficient of the reduces rotation matrix" << endl;
                    return 0;
            }
            break;
        case 0 :
            switch(j)
            {
                case -2 :
                    return(sqrt(3/8.)*sin(beta)*sin(beta));
                    break;
                case -1 :
                    return(-sqrt(3/2.)*sin(beta)*cos(beta));
                    break;
                case 0 :
                    return(1/2.*(3.*cos(beta)*cos(beta)-1));
                    break;
                case 1 :
                    return(sqrt(3/2.)*sin(beta)*cos(beta));
                    break;
                case 2 :
                    return(sqrt(3/8.)*sin(beta)*sin(beta));
                    break;
                default:
                    cout << "Error to compute a coefficient of the reduces rotation matrix" << endl;
                    return 0;
            }
            break;
        case 1 :
            switch(j)
            {
                case -2 :
                    return(1/2.*sin(beta)*(cos(beta)-1));
                    break;
                case -1 :
                    return(1/2.*(2.*cos(beta)+1)*(1-cos(beta)));
                    break;
                case 0 :
                    return(-sqrt(3/2.)*sin(beta)*cos(beta));
                    break;
                case 1 :
                    return(1/2.*(2.*cos(beta)-1)*(cos(beta)+1));
                    break;
                case 2 :
                    return(1/2.*sin(beta)*(1+cos(beta)));
                    break;
                default:
                    cout << "Error to compute a coefficient of the reduces rotation matrix" << endl;
                    return 0;
            }
            break;
        case 2 :
            switch(j)
            {
                case -2 :
                    return(pow(sin(beta/2),4));
                    break;
                case -1 :
                    return(1/2.*sin(beta)*(cos(beta)-1));
                    break;
                case 0 :
                    return(sqrt(3/8.)*sin(beta)*sin(beta));
                    break;
                case 1 :
                    return(-1/2.*sin(beta)*(1+cos(beta)));
                    break;
                case 2 :
                    return(pow(cos(beta/2),4));
                    break;
                default:
                    cout << "Error to compute a coefficient of the reduces rotation matrix" << endl;
                    return 0;
            }
            break;
        default:
            cout << "Error to compute a coefficient of the reduces rotation matrix" << endl;
            return 0;
    }
}

//------------- NON INTEGRATED fp gp hp --------------
double func_fp(double gamma_d , double lambda , double incl , double longi)
{
    double thetaB = acos(2*cos(lambda)/sqrt(1+3*(pow(cos(lambda),2))));
    double gamma = gamma_d/2*sqrt(1+3*pow(cos(lambda),2));

    int compt_Qp,compt_Qpp;
    double sum = 0;
    for(compt_Qp=-2;compt_Qp<=2;compt_Qp++)
    {
        for(compt_Qpp=-2;compt_Qpp<=2;compt_Qpp++)
        {
            sum += 1/(1+gamma*gamma*compt_Qp*compt_Qp)*reduce_rot_mat(thetaB,0,compt_Qp)*reduce_rot_mat(-thetaB-lambda,compt_Qp,compt_Qpp)*reduce_rot_mat(-incl,compt_Qpp,0)*(cos(compt_Qpp*longi)-gamma*compt_Qp*sin(compt_Qpp*longi));
        }
    }
    return(sum);
}

double func_gp(double gamma_d ,double lambda , double incl , double longi)
{
    double thetaB = acos(2*cos(lambda)/sqrt(1+3*(pow(cos(lambda),2))));
    double gamma = gamma_d/2*sqrt(1+3*pow(cos(lambda),2));
    int compt_Qp,compt_Qpp;
    double sum = 0;
    for(compt_Qp=-2;compt_Qp<=2;compt_Qp++)
    {
        for(compt_Qpp=-2;compt_Qpp<=2;compt_Qpp++)
        {
            sum += 1/(1+gamma*gamma*compt_Qp*compt_Qp)*reduce_rot_mat(thetaB,0,compt_Qp)*reduce_rot_mat(-thetaB-lambda,compt_Qp,compt_Qpp)*reduce_rot_mat(-incl,compt_Qpp,2)*(cos(compt_Qpp*longi)-gamma*compt_Qp*sin(compt_Qpp*longi));
        }
    }
    return(sum);
}

double func_hp(double gamma_d ,double lambda , double incl , double longi)
{
    double thetaB = acos(2*cos(lambda)/sqrt(1+3*(pow(cos(lambda),2))));
    double gamma = gamma_d/2*sqrt(1+3*pow(cos(lambda),2));
    int compt_Qp,compt_Qpp;
    double sum = 0;
    for(compt_Qp=-2;compt_Qp<=2;compt_Qp++)
    {
        for(compt_Qpp=-2;compt_Qpp<=2;compt_Qpp++)
        {
            sum += -1/(1+gamma*gamma*compt_Qp*compt_Qp)*reduce_rot_mat(thetaB,0,compt_Qp)*reduce_rot_mat(-thetaB-lambda,compt_Qp,compt_Qpp)*reduce_rot_mat(-incl,compt_Qpp,2)*(gamma*compt_Qp*cos(compt_Qpp*longi)+sin(compt_Qpp*longi));
        }
    }
    return(sum);
}

//new functions, to perform the integration extended on J02 and the collision variables which depend on the hydrogen atoms density
double func_fp_tilde(double small_gamma , double lambda , double incl , double longi , double mu)
{
    double Nhydro = hydro(mu);
    double D_2 = 37.0942*1e-9*pow((temperature/5000),0.5247)*pow(0.7198,temperature/5000)*Nhydro;
    double C_ul = 8.42791*1e-9*pow(temperature/5000,0.38)*Nhydro - D_2;
    double epsilon = C_ul/(A_ul+C_ul);
    double gamma_d = small_gamma*(1-epsilon)/(1+D_2/A_ul*(1-epsilon));
    double beta=u/(1-u);
    double J00 = I0*(1-u)/(2+2*beta*mu-2*tab_exp_int[0][int(mu*1000)]+beta*tab_exp_int[1][int(mu*1000)]) ;
    double J02 = J00/(2*sqrt(2))*(tab_exp_int[0][int(mu*1000)]-beta*tab_exp_int[1][int(mu*1000)]-3*tab_exp_int[2][int(mu*1000)]+3*beta*tab_exp_int[3][int(mu*1000)])/(2+2*beta*mu-2*tab_exp_int[0][int(mu*1000)]+beta*tab_exp_int[1][int(mu*1000)]) ;
    return((1-epsilon)/(1+D_2/A_ul*(1-epsilon))*J02*func_fp(gamma_d,lambda,incl,longi));
}

double func_gp_tilde(double small_gamma , double lambda , double incl , double longi , double mu)
{
    double Nhydro = hydro(mu);
    double D_2 = 37.0942*1e-9*pow((temperature/5000),0.5247)*pow(0.7198,temperature/5000)*Nhydro;
    double C_ul = 8.42791*1e-9*pow(temperature/5000,0.38)*Nhydro - D_2;
    double epsilon = C_ul/(A_ul+C_ul);
    double gamma_d = small_gamma*(1-epsilon)/(1+D_2/A_ul*(1-epsilon));
    double beta=u/(1-u);
    double J00 = I0*(1-u)/(2+2*beta*mu-2*tab_exp_int[0][int(mu*1000)]+beta*tab_exp_int[1][int(mu*1000)]) ;
    double J02 = J00/(2*sqrt(2))*(tab_exp_int[0][int(mu*1000)]-beta*tab_exp_int[1][int(mu*1000)]-3*tab_exp_int[2][int(mu*1000)]+3*beta*tab_exp_int[3][int(mu*1000)])/(2+2*beta*mu-2*tab_exp_int[0][int(mu*1000)]+beta*tab_exp_int[1][int(mu*1000)]) ;
    return((1-epsilon)/(1+D_2/A_ul*(1-epsilon))*J02*func_gp(gamma_d,lambda,incl,longi));
}

double func_hp_tilde(double small_gamma , double lambda , double incl , double longi , double mu)
{
    double Nhydro = hydro(mu);
    double D_2 = 37.0942*1e-9*pow((temperature/5000),0.5247)*pow(0.7198,temperature/5000)*Nhydro;
    double C_ul = 8.42791*1e-9*pow(temperature/5000,0.38)*Nhydro - D_2;
    double epsilon = C_ul/(A_ul+C_ul);
    double gamma_d = small_gamma*(1-epsilon)/(1+D_2/A_ul*(1-epsilon));
    double beta=u/(1-u);
    double J00 = I0*(1-u)/(2+2*beta*mu-2*tab_exp_int[0][int(mu*1000)]+beta*tab_exp_int[1][int(mu*1000)]) ;
    double J02 = J00/(2*sqrt(2))*(tab_exp_int[0][int(mu*1000)]-beta*tab_exp_int[1][int(mu*1000)]-3*tab_exp_int[2][int(mu*1000)]+3*beta*tab_exp_int[3][int(mu*1000)])/(2+2*beta*mu-2*tab_exp_int[0][int(mu*1000)]+beta*tab_exp_int[1][int(mu*1000)]) ;
    return((1-epsilon)/(1+D_2/A_ul*(1-epsilon))*J02*func_hp(gamma_d,lambda,incl,longi));
}

double func_op_tilde(double small_gamma , double lambda , double incl , double longi , double mu)
{
    double Nhydro = hydro(mu);
    double D_0 = 43.4158*1e-9*pow((temperature/5000),0.4514)*pow(0.8079,temperature/5000)*Nhydro;
    double D_2 = 37.0942*1e-9*pow((temperature/5000),0.5247)*pow(0.7198,temperature/5000)*Nhydro;
    double C_ul = 8.42791*1e-9*pow(temperature/5000,0.38)*Nhydro - D_2;
    double epsilon = C_ul/(A_ul+C_ul);
    double beta=u/(1-u);
    double J00 = I0*(1-u)/(2+2*beta*mu-2*tab_exp_int[0][int(mu*1000)]+beta*tab_exp_int[1][int(mu*1000)]);
    return(((1-epsilon)*w_JuJl0*J00+epsilon*B_nu)/((1+D_0/A_ul*(1-epsilon))));
}

double func_fp_barre_tilde(double small_gamma, double incl)
{
    if (incl>M_PI/2) incl = M_PI - incl;
    return(integ_tilde(func_fp_tilde,small_gamma,incl));
}

double func_gp_barre_tilde(double small_gamma , double incl)
{
    if (incl>M_PI/2) incl = M_PI - incl;
    return(integ_tilde(func_gp_tilde,small_gamma,incl));
}

double func_hp_barre_tilde(double small_gamma , double incl)
{
    if (incl<M_PI/2) return(integ_tilde(func_hp_tilde,small_gamma,incl));
    else return(-integ_tilde(func_hp_tilde,small_gamma,M_PI-incl));
}

double func_op_barre_tilde(double small_gamma , double incl)
{
    return(integ_tilde(func_op_tilde,small_gamma,incl));
}

//*********************************************************************************************************
//******************** COMPUTATION *****************************************************************
//*********************************************************************************************************

void computing(double B_field, double psincl)
{
    //*********************************************************************************************************
    //******************** VARIABLES **************************************************************************
    //*********************************************************************************************************
    double xi,incl,jncl;
    psincl = psincl*M_PI/180;//inclination of the rotation axis to the dipole axis
    double angf; //rotational phase angle

    double nu_larmor = 1.3996e6*B_field; //Larmor frequency
    double small_gamma = 2*M_PI*nu_larmor*gu/A_ul;

    double ro00, ro02, re_ro22, im_ro22;

    double fun_phi = 1/(4*M_PI*M_PI)*(A_ul+B_lu)/((A_ul+B_lu)*(A_ul+B_lu)/(16*M_PI*M_PI)+nu*nu); //absorption profile
    double epsi_0 = h_planck*nu/(4*M_PI)*A_ul*N_tot*sqrt(2*Ju+1)*fun_phi;

    double I_ref;
    double epsi_I_line, epsi_Q_line, epsi_U_line;
    double Stokes_I, Stokes_Q, Stokes_U;

    double vartemp1; //useful in huge calculations, to make intermediate steps

    //*********************************************************************************************************
    //******************** COMPUTING **************************************************************************
    //*********************************************************************************************************

    double h=1; int n_t = 365, i_t; double day;

    std::stringstream fichier_i;
    fichier_i << "BDD/i/i_" << B_field << "_" << psincl*180/M_PI << ".txt";
    std::stringstream fichier_p;
    fichier_p << "BDD/p/p_" << B_field << "_" << psincl*180/M_PI << ".txt";
    std::stringstream fichier_a;
    fichier_a << "BDD/a/a_" << B_field << "_" << psincl*180/M_PI << ".txt";
    std::stringstream fichier_q;
    fichier_q << "BDD/q/q_" << B_field << "_" << psincl*180/M_PI << ".txt";
    std::stringstream fichier_u;
    fichier_u << "BDD/u/u_" << B_field << "_" << psincl*180/M_PI << ".txt";

    fstream resultat_i(fichier_i.str().c_str(),ios::out);
    fstream resultat_p(fichier_p.str().c_str(),ios::out);
    fstream resultat_a(fichier_a.str().c_str(),ios::out);
    fstream resultat_q(fichier_q.str().c_str(),ios::out);
    fstream resultat_u(fichier_u.str().c_str(),ios::out);


    for(i_t=1;i_t<=n_t;i_t++)
    {
    day = i_t*h;

    //******************************************************************
    //geometric angles
    jncl = (90-rotation_axis[i_t-1])*M_PI/180;
    angf = 2*M_PI*(float)(int(day)%26)/26.; //we suppose that the LOS, omega and e are coplanar at t = 0 (cf paper).
    incl = acos(cos(psincl)*cos(jncl)-sin(psincl)*sin(jncl)*cos(angf));
    if (incl==0) incl = 0.000001;
    xi = arctan(1/sin(incl)*(sin(psincl)*cos(angf)*cos(jncl)+cos(psincl)*sin(jncl)) , -sin(psincl)*sin(angf)/sin(incl) ) ;
    if (sin(incl)==0) xi = M_PI;

    //******************************************************************
    //Calculation of the spherical components of the density matrix
    vartemp1 = 2*h_planck*pow(nu,3)/(c_light*c_light)*(2*Jl+1)/sqrt(2*Ju+1);
    ro00 = func_op_barre_tilde(small_gamma , incl)/vartemp1;
    ro02 = w_JuJl2*func_fp_barre_tilde(small_gamma , incl)/vartemp1;
    re_ro22 = w_JuJl2*func_gp_barre_tilde(small_gamma , incl)/vartemp1;
    im_ro22 = w_JuJl2*func_hp_barre_tilde(small_gamma , incl)/vartemp1;

    //******************************************************************
    //Calculation of  the emissivity coefficients
    epsi_I_line = epsi_0*(ro00+w_JuJl2/sqrt(2)*ro02);
    epsi_Q_line = -epsi_0*w_JuJl2*sqrt(3)*(cos(2*xi)*re_ro22-sin(2*xi)*im_ro22);
    epsi_U_line = epsi_0*w_JuJl2*sqrt(3)*(sin(2*xi)*re_ro22+cos(2*xi)*im_ro22);

    //*****************************************************************
    //Calculation of the Stokes parameters
    Stokes_I = epsi_I_line;
     if (day == 1) I_ref = Stokes_I;//keep in memory the intensity at day=1
    Stokes_Q = Stokes_I*epsi_Q_line/epsi_I_line/factor;
    Stokes_U = Stokes_I*epsi_U_line/epsi_I_line/factor;

    resultat_i << day << " "  << (Stokes_I-I_ref)/I_ref << endl;
    resultat_p << day << " "  << sqrt(Stokes_Q*Stokes_Q+Stokes_U*Stokes_U)/Stokes_I << endl;
    resultat_a << day << " "  << 0.5*arctan(Stokes_Q,Stokes_U) << endl;
    resultat_q << day << " "  << (Stokes_Q/Stokes_I) << endl;
    resultat_u << day << " "  << (Stokes_U/Stokes_I) << endl;

    }
    resultat_i.close();
    resultat_p.close();
    resultat_a.close();
    resultat_q.close();
    resultat_u.close();
}

int main()
{
    int start;
    cout << "*************************  PHASM (T.VIEU 2016)  *************************" << endl; cout << endl;
    cout << "*************************************************************************" << endl; cout << endl;
    cout << "*************************  Welcome to PHASM-DB  *************************" << endl; cout << endl;
    cout << ">>>>> Start the computation (0 / 1) ? <<<<<" << endl; cin >> start; if (start == 0) return 0; cout << endl;

    init_tab_exp_int();
    init_tab_rotation();

    //it is possible to change the steps here
    double stepB = 20;
    double steppsi = 90;
    double Bmin=0,Bmax=20,psimin=0,psimax=90;

    int iB,ipsi;
    int nB=(Bmax-Bmin)/stepB; int npsi=(psimax-psimin)/steppsi;
    for(iB=0;iB<=nB;iB++)
    {
        for(ipsi=0;ipsi<=npsi;ipsi++)
        {
            cout << iB*90/nB+ipsi*100/(10*npsi) << "%  ";
            computing(Bmin+iB*stepB , psimin+ipsi*steppsi);
        }
    }

    cout << endl; cout << endl; cout << "DONE" << endl;
    return 0;
}
