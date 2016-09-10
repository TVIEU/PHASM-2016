#include <iostream>
#include<fstream>
#include<stdlib.h>
#include<math.h>
#include <sstream>
#include<cstdlib>
#include <iomanip>

using namespace std;

double tab_exp_int[4][1001]; //the exponential integrals are tabulated in the "exp_int" folder
double rotation_axis[366]; //the modulation of the rotation axis during the year 2016 is tabulated in the "rotation_axis.txt" file.

//*********************************************************************************************************
//******************** CONSTANTS, GLOBAL VARIABLES ********************************************************
//*********************************************************************************************************

double h_planck = 6.62607004e-34; //Planck constant
double c_light = 299792458; //speed of light
double k_boltzmann = 1.38064852e-23; //Boltzmann constant

double I0 = 4.35*1e13; //intensity
double B_field; //strength of the magnetic field (G)
double psincl;//inclination of the rotation axis to the dipole axis
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

double factor = 10.1; //correction, depolarization factor

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

//exponential integral of order n (NOT USED BECAUSE TOO LONG : THEY ARE PRE-COMPUTED AND TABULATED WITH THE FUNCTION calcul_file_exp_int()
double exp_int(double x,int n)
{
    double res = 0;
    double b=1000,pas=0.01;
    int i ; int ncompt = (b-1)/pas;
    for(i=0;i<=ncompt;i++)
        res += pas*exp(-x*(1+i*pas))/pow((1+i*pas),n);

    return (res);
}

//compute the file containing the tabulated exponential integral. NOT USED, USE THE FOLDER AND FILES ALREADY COMPUTED
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
        fun.close();
    }
}

void init_tab_exp_int() //put the tables of the exponential integral from the pre-computed files to the array.
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

void files_for_Python(double h)// creates two files containing the values of alpha and ro, which is useful for python plots
{
    int n_ro, n_alpha, i_ro, i_alpha;
    n_ro=6/h; n_alpha=(2*M_PI)/h+1;
    fstream res_alpha("res_alpha.txt",ios::out);
    fstream res_ro("res_ro.txt",ios::out);
    for(i_ro=0; i_ro<=n_ro; i_ro++) res_ro << i_ro*h/6 << " " << 0 << endl;
    for(i_alpha=0; i_alpha<=n_alpha; i_alpha++) res_alpha << i_alpha*h << " " << 0 << endl;
    res_alpha.close();
    res_ro.close();
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
            sum += 1/(1+gamma*gamma*compt_Qp*compt_Qp)*reduce_rot_mat(thetaB,0,compt_Qp)*reduce_rot_mat(-thetaB-lambda,compt_Qp,compt_Qpp)*reduce_rot_mat(-incl,compt_Qpp,0)*(cos(compt_Qpp*longi)-gamma*compt_Qp*sin(compt_Qpp*longi));
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
            sum += 1/(1+gamma*gamma*compt_Qp*compt_Qp)*reduce_rot_mat(thetaB,0,compt_Qp)*reduce_rot_mat(-thetaB-lambda,compt_Qp,compt_Qpp)*reduce_rot_mat(-incl,compt_Qpp,2)*(cos(compt_Qpp*longi)-gamma*compt_Qp*sin(compt_Qpp*longi));
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
            sum += -1/(1+gamma*gamma*compt_Qp*compt_Qp)*reduce_rot_mat(thetaB,0,compt_Qp)*reduce_rot_mat(-thetaB-lambda,compt_Qp,compt_Qpp)*reduce_rot_mat(-incl,compt_Qpp,2)*(gamma*compt_Qp*cos(compt_Qpp*longi)+sin(compt_Qpp*longi));
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

//integration using the previous functions and the quadratic integral
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
//******************** MAIN *******************************************************************************
//*********************************************************************************************************

int main()
{
    cout << "*************************  PHASM (T.VIEU 2016)  *************************" << endl; cout << endl;
    cout << "*************************************************************************" << endl; cout << endl;
    cout << "*************************  Welcome to PHASM-MI  *************************" << endl; cout << endl;
    cout << ">>>>> Enter the geometry of a global dipolar magnetic field <<<<<" << endl; cout << endl;

    cout << "Magnetic dipole strength ? "; cin >> B_field;
    cout << "Inclination of the dipole axis to the rotation axis ? "; cin >> psincl; psincl = psincl*M_PI/180; cout << endl;

    //calcul_file_exp_int(); //use only if you want to compute again the tabulated exponential integrals !
    init_tab_exp_int();
    init_tab_rotation();

    //*********************************************************************************************************
    //******************** VARIABLES **************************************************************************
    //*********************************************************************************************************

    double mu;
    double lambda,coslongi,sinlongi,longi;

    double xi;//defines the polarization reference system
    double jncl;//inclination of the rotation axis to the LOS
    double incl;//inclination of the stellar dipole with respect to the LOS.
    double angf; //rotational phase angle

    double nu_larmor = 1.3996e6*B_field; //Larmor frequency
    double Nhydro; //density of hydrogen atoms
    double D_2,C_ul; //collision rates
    double small_gamma = 2*M_PI*nu_larmor*gu/A_ul;
    double gamma_d,epsilon;

    double ro00, ro02, re_ro22, im_ro22;

    double fun_phi = 1/(4*M_PI*M_PI)*(A_ul+B_lu)/((A_ul+B_lu)*(A_ul+B_lu)/(16*M_PI*M_PI)+nu*nu); //absorption profile
    double epsi_0 = h_planck*nu/(4*M_PI)*A_ul*N_tot*sqrt(2*Ju+1)*fun_phi; //does not matter since it simplifies when computing the fractional Stokes parameters X/I

    double I_ref;
    double epsi_I_line, epsi_Q_line, epsi_U_line; //emissivity
    double Stokes_I, Stokes_Q, Stokes_U; //Stokes parameters

    double vartemp1; //useful in huge calculations, to make intermediate steps

    char mode,request;
    cout << ">>>  Which result do you want to get ?  <<<" << endl; //m for map ; p for plot
    cout << "'m' ............  Distribution maps on the projected disk" << endl;
    cout << "'i' ............  Evolution of the integrated signals during the year 2016" << endl;
    cin >> mode; cout << endl;

//*********************************************************************************************************
//******************** DISTRIBUTION MAPS ******************************************************************
//*********************************************************************************************************

if (mode == 'm') //2D maps
    {
    cout << ">>>  Which distribution ?  <<<" << endl;
    cout << "'m' ............  Distribution of gamma" << endl;
    cout << "'i' ............  Stokes parameter I" << endl;
    cout << "'q' ............  Stokes parameter Q" << endl;
    cout << "'u' ............  Stokes parameter U" << endl;
    cout << "'p' ............  Stokes parameter P" << endl;
    cout << "'a' ............  Stokes parameter ALPHA" << endl;
    cout << "'d' ............  Signals along all directions (spread)" << endl;
    cin >> request; cout << endl;

    int day;
    cout << ">>>  Day ? (between 1 and 366)  <<<" << endl;
    cin >> day;
    cout << endl;

    //geometric angles
    jncl = (90-rotation_axis[day-1])*M_PI/180;
    angf = 2*M_PI*(float)(int(day)%26)/26.;
    incl = acos(cos(psincl)*cos(jncl)-sin(psincl)*sin(jncl)*cos(angf));
    if (incl==0) incl = 0.000001;
    xi = arctan(1/sin(incl)*(sin(psincl)*cos(angf)*cos(jncl)+cos(psincl)*sin(jncl)) , -sin(psincl)*sin(angf)/sin(incl)) ;
    if (sin(incl)==0) xi = M_PI;//not expected to happen

    double alpha,ro;
    double h=0.02;
    int n_ro, n_alpha, i_ro, i_alpha;
    n_ro=6/h; n_alpha=(2*M_PI)/h+1;

    files_for_Python(h);// creates two files containing the values of alpha and ro, which is useful for python plots
    fstream resultat("resultat.txt",ios::out);

    for(i_alpha=1; i_alpha<=n_alpha; i_alpha++)
    {
    if(i_alpha*100/n_alpha != (i_alpha-1)*100/n_alpha) cout << i_alpha*100/n_alpha << " " << "% ";
    for(i_ro=0; i_ro<=n_ro; i_ro++)
    {
    alpha = i_alpha*h; ro = i_ro*h/6; mu = sqrt(1-ro*ro);

    //collisions:
    Nhydro = hydro(sqrt(1-ro*ro));
    D_2 = 37.0942*1e-9*pow((temperature/5000),0.5247)*pow(0.7198,temperature/5000)*Nhydro;
    C_ul = 8.42791*1e-9*pow(temperature/5000,0.38)*Nhydro - D_2;
    epsilon = C_ul/(A_ul+C_ul);
    gamma_d = small_gamma*(1-epsilon)/(1+D_2/A_ul*(1-epsilon));

    //spherical coordinates lambda, longi on the solar normalized sphere from polar coordinates ro, alpha from the solar disk:
    lambda = acos(ro*cos(alpha)*sin(incl)-cos(incl)*sqrt(1-ro*ro));
    coslongi = (sin(incl)*sqrt(1-ro*ro)+cos(incl)*cos(alpha)*ro);
    sinlongi = ro*sin(alpha);
    longi = arctan(coslongi,sinlongi);

    //**************************************************************************************
    //Calculation of the spherical components of the density matrix
    vartemp1 = 2*h_planck*pow(nu,3)/(c_light*c_light)*(2*Jl+1)/sqrt(2*Ju+1);
    ro00 = func_op_tilde(small_gamma , lambda , incl , longi, mu)/vartemp1;
    ro02 = w_JuJl2*func_fp_tilde(small_gamma , lambda , incl , longi, mu)/vartemp1;
    re_ro22 = w_JuJl2*func_gp_tilde(small_gamma , lambda , incl , longi, mu)/vartemp1;
    im_ro22 = w_JuJl2*func_hp_tilde(small_gamma , lambda , incl , longi, mu)/vartemp1;

    //*******************************************************************************
    //Calculation of the emissivity coefficients
    epsi_I_line = epsi_0*(ro00+w_JuJl2/sqrt(2)*ro02);
    epsi_Q_line = -epsi_0*w_JuJl2*sqrt(3)*(cos(2*xi)*re_ro22-sin(2*xi)*im_ro22);
    epsi_U_line = epsi_0*w_JuJl2*sqrt(3)*(sin(2*xi)*re_ro22+cos(2*xi)*im_ro22);

    //******************************************************************************
    //Calculation of the Stokes parameters
    Stokes_I = epsi_I_line;
    Stokes_Q = Stokes_I*epsi_Q_line/epsi_I_line/factor;
    Stokes_U = Stokes_I*epsi_U_line/epsi_I_line/factor;

    if (request == 'm') resultat << alpha << " " << ro << " " << gamma_d/2*sqrt(1+3*pow(cos(lambda),2)) << endl;
    if (request == 'i') resultat << alpha << " " << ro << " " << Stokes_I << endl;
    if (request == 'q') resultat << alpha << " " << ro << " " << Stokes_Q/Stokes_I << endl;
    if (request == 'u') resultat << alpha << " " << ro << " " << Stokes_U/Stokes_I << endl;
    if (request == 'p') resultat << alpha << " " << ro << " " << sqrt(Stokes_Q*Stokes_Q+Stokes_U*Stokes_U) / Stokes_I << endl;
    if (request == 'a') resultat << alpha << " " << ro << " " << 0.5*arctan(Stokes_Q,Stokes_U) <<endl;
    if (request == 'd') resultat << alpha << " " << sqrt(1-ro*ro) << " " << sqrt(Stokes_Q*Stokes_Q+Stokes_U*Stokes_U)/Stokes_I  << endl;

    }
    }
    resultat.close();
    }

//************************************************************************************
//******************** INTEGRATED SIGNATURES *****************************************
//************************************************************************************

if (mode == 'i') //integrated signals
    {
    char request;
    cout << ">>>  What do you want ?  <<<" << endl;
    cout << "'i' ............  Stokes parameter I" << endl;
    cout << "'q' ............  Stokes parameter Q" << endl;
    cout << "'u' ............  Stokes parameter U" << endl;
    cout << "'p' ............  Stokes parameter P" << endl;
    cout << "'a' ............  Stokes parameter ALPHA" << endl;
    cout << "'d' ............  Diagramme Q VS U" << endl;
    cin >> request; cout << endl;

    fstream resultat("resultat.txt",ios::out);
    double h=1; int n_t = 366, i_t; double day; //h is the temporal step (in days)
    for(i_t=1;i_t<=n_t;i_t++)
    {
    day = i_t*h; //day of the year
    if(i_t*100/n_t != (i_t-1)*100/n_t) cout << i_t*100/n_t << "% ";

    //******************************************************************
    //geometric angles (modulated with the rotation)
    jncl = (90-rotation_axis[i_t-1])*M_PI/180;
    angf = 2*M_PI*(float)(int(day)%26)/26.;
    incl = acos(cos(psincl)*cos(jncl)-sin(psincl)*sin(jncl)*cos(angf));
    if (incl==0) incl = 0.000001; //border effect, singularity
    xi = arctan(1/sin(incl)*(sin(psincl)*cos(angf)*cos(jncl)+cos(psincl)*sin(jncl)),-sin(psincl)*sin(angf)/sin(incl)) ;
    if (sin(incl)==0) xi = M_PI; //not expected to happen

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

    if (request == 'i') resultat << day << " "  << (Stokes_I-I_ref)/I_ref << endl;
    if (request == 'p') resultat << day << " "  << sqrt(Stokes_Q*Stokes_Q+Stokes_U*Stokes_U)/Stokes_I << endl;
    if (request == 'a') resultat << day << " "  << 0.5*arctan(Stokes_Q,Stokes_U) << endl;
    if (request == 'q') resultat << day << " "  << (Stokes_Q/Stokes_I) << endl;
    if (request == 'u') resultat << day << " "  << (Stokes_U/Stokes_I) << endl;
    if (request == 'd') resultat << Stokes_Q/Stokes_I << " "  << Stokes_U/Stokes_I << endl;

    }
    resultat.close();
    }

//*********************************************************************************************************
//********************  PLOT  *****************************************************************************
//*********************************************************************************************************

//generation of the Python code for the plot, stored in 'plot_python.txt'
    fstream plotpy("plot_python.txt",ios::out);
    plotpy << "import numpy as np" << endl;
    plotpy << "import matplotlib.pyplot as plt" << endl;
    plotpy << "from numpy import *" << endl;
    plotpy << "from matplotlib.pyplot import *" << endl;
    plotpy << "data=loadtxt('resultat.txt')" << endl;
    if (mode == 'm' && request != 'd')
    {
        plotpy << "res_alpha=loadtxt('res_alpha.txt')" << endl;
        plotpy << "res_ro=loadtxt('res_ro.txt')" << endl;
        plotpy << "alpha=res_alpha[:,0]" << endl;
        plotpy << "ro=res_ro[:,0]" << endl;
        plotpy << "r, theta = np.meshgrid(ro, alpha)" << endl;
        plotpy << "values = np.resize(data[:,2],(len(alpha),len(ro)))" << endl;
        plotpy << "fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))" << endl;
        plotpy << "ax.set_theta_zero_location(\"S\")" << endl;
        plotpy << "ax.contourf(theta, r, values)" << endl;
    }
    if(mode == 'm' && request == 'd')
    {
        plotpy << "plt.semilogy()" << endl;
        plotpy << "plt.plot(data[:,1],data[:,2])" << endl;
    }
    if(mode == 'i')
        plotpy << "plt.plot(data[:,0],data[:,1])" << endl;
    plotpy << "plt.show()" << endl;
    plotpy.close();

    int plot;
    cout << endl; cout << endl;
    cout << ">>>  Do you want to plot the result ? (you need the right directory in the system command for Python call at the very end of this main.c++)  <<<" << endl; //m for map ; p for plot
    cout << "0 ............  No" << endl;
    cout << "1 ............  Yes" << endl;
    cin >> plot; cout << endl;

    if(plot == 1) system("\"C:/WinPython-64bit-3.4.4.2/python-3.4.4.amd64/python.exe\" plot_python.txt");
    //put here the directory to your python application. If you can't, compute the results without plotting and look to the instructions in "PYTHON CODES FOR THE PLOTS"
    if(plot == 0) cout << "Use the Python code in plot_python.txt to plot the result. More information in the README file." << endl;
    cout << endl; cout << "DONE" << endl;

    return 0;
}
