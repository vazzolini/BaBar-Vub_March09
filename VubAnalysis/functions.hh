// Polynomials
double f_p1(double *x, double *par);
double f_p2(double *x, double *par);

// Exponential
double f_expo(double *x, double *par);
// Gauss: this is area/mean/sigma
double f_Gauss(double *x, double *par);
double f_2G(double *x, double *par);
// Gauss: This is with max/mean/sigma
double f_gauss(double *x, double *par);
double f_2g(double *x, double *par);
double f_3g(double *x, double *par);

// Crystall Ball
double f_cb(double *x, double *par);
double f_p1acb(double *x, double *par);
double f_p2acb(double *x, double *par);
// Novosibirsk
double f_fnov(double *x, double *par);

// Argus (only)
double f_argus(double *x, double *par);
// Argus and Gauss
double f_aag(double *x, double *par);
// Argus and Crystall Ball
double f_aacb(double *x, double *par);

// Polynomials and Gayss
double f_p0ag(double *x, double *par);
double f_p1ag(double *x, double *par);
double f_p2ag(double *x, double *par);
double f_p0a2g(double *x, double *par);

double f2_chi2ellipsis(double *x, double *par);
