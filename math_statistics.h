/*! \file math_statistics.h
 * \brief Used for statistics acquisition.
 */
// IceGrain: Extraction and parameterization of grain boundary networks of ice
//
// Copyright (c) 2013 Tobias Binder.
// 
// This software was developed at the University of Heidelberg by
// Tobias Binder, Bjoern Andres, Thorsten Beier and Arthur Kuehlwein.
// Enquiries shall be directed to tobias.binder@iwr.uni-heidelberg.de.
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// - Redistributions in binary form must reproduce the above copyright notice, 
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// - All advertising materials mentioning features or use of this software must 
//   display the following acknowledgement: ``This product includes the IceGrain
//   package developed by Tobias Binder and others''.
// - The name of the author must not be used to endorse or promote products 
//   derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED 
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <algorithm>

template<class T>
inline T square(const T& x)
{
    return x*x;
}

bool cmp (float a, float b)
{
    return (a > b);
}

struct grain
{
    int x;
    int y;
    long size;
};

bool grain_cmp (grain a, grain b)
{
    return (a.size > b.size);
}

struct size_index
{
    int index;
    float size;
};

bool size_index_cmp (size_index a, size_index b)
{
    return (a.size > b.size);
}

struct area_range
{
    int x_low;
    int y_low;
    int x_high;
    int y_high;
};

struct data
{
    size_t n;
    double * y;
};

struct correlation
{
    float rho;
    float m;
    float b;
    int nr;
};

struct element_entry
{
    float value;
    int image;
    int element;
};
/*! \struct point3d
 * \brief A 3D point.
 */
struct point3d
{
    float x; /*!< x-coordinate*/
    float y; /*!< y-coordinate*/
    float z; /*!< z-coordinate*/
};

struct param_entry
{
    float size;
    float flat;
    float round;
    float boxflat;
    float width;
    float height;
    float angle;
    float ratio;
};

bool entry_cmp (element_entry a, element_entry b)
{
    return (a.value > b.value);
}

bool entry_cmp_abs (element_entry a, element_entry b)
{
    return (fabs(a.value) > fabs(b.value));
}

bool x_cmp (point3d a, point3d b)
{
    return (a.x > b.x);
}

bool y_cmp (point3d a, point3d b)
{
    return (a.y > b.y);
}

bool size_cmp (param_entry a, param_entry b)
{
    return (a.size > b.size);
}

void calculate_mean_standard_deviation(std::vector<param_entry> values, float & mean, float & standard_deviation)
{
    mean=0;
    for(int i=0;i<(int)values.size();i++)
    {
        mean=mean+values[i].size;
    }
    mean=mean/values.size();

    standard_deviation=0;
    if(values.size()>1)
    {
        for(int i=0;i<(int)values.size();i++)
        {
            standard_deviation=standard_deviation + ((values[i].size-mean)*(values[i].size-mean));
        }
        standard_deviation=sqrt(standard_deviation*(1.0/((float)values.size()-1.0)));
    }
}

void calculate_mean_standard_deviation(std::vector<float> values, float & mean, float & standard_deviation)
{
    mean=0;
    for(int i=0;i<(int)values.size();i++)
    {
        mean=mean+values[i];
    }
    mean=mean/values.size();

    standard_deviation=0;
    if(values.size()>1)
    {
        for(int i=0;i<(int)values.size();i++)
        {
            standard_deviation=standard_deviation + ((values[i]-mean)*(values[i]-mean));
        }
        standard_deviation=sqrt(standard_deviation*(1.0/((float)values.size()-1.0)));
    }
}

float get_length(std::vector<point> pixels)
{
    if (pixels.size()<2) return (float)pixels.size();

    bool horizontal_step=false;
    bool vertical_step=false;
    float length=0.0f;

    for (int i=1; i<pixels.size(); i++)
    {
        if (fabs(pixels[i].x-pixels[i-1].x)==1 && fabs(pixels[i].y-pixels[i-1].y)==0 && vertical_step)
        {
            length+=(sqrt(2)-1.0f);
            horizontal_step=false;
            vertical_step=false;
        }
        else if (fabs(pixels[i].x-pixels[i-1].x)==0 && fabs(pixels[i].y-pixels[i-1].y)==1 && horizontal_step)
        {
            length+=(sqrt(2)-1.0f);
            horizontal_step=false;
            vertical_step=false;
        }
        else
        {
            length+=sqrt((pixels[i].x-pixels[i-1].x)*(pixels[i].x-pixels[i-1].x)+
                (pixels[i].y-pixels[i-1].y)*(pixels[i].y-pixels[i-1].y));

            if (fabs(pixels[i].x-pixels[i-1].x)==1 && fabs(pixels[i].y-pixels[i-1].y)==0) horizontal_step=true;
            else horizontal_step=false;
            if (fabs(pixels[i].x-pixels[i-1].x)==0 && fabs(pixels[i].y-pixels[i-1].y)==1) vertical_step=true;
            else vertical_step=false;
        }
    }
    return length;
}

bool grain_boundary_inverted(std::vector<point> pixels)
{
    if (pixels.size()<2) return false;
    else if (pixels[0].x>pixels.back().x || (pixels[0].x==pixels.back().x && pixels[0].y>pixels.back().y)) return true;
    else return false;
}

float calculate_angle(std::vector< std::vector<point> > points, point junction, point other)
{
    if(points.size()!=2)
    {
        std::cout<<"Error: Two filled vectors needed!"<<std::endl;
        exit(-1);
    }

    if(other.x==-1)
    {
        std::cout<<"Error: No other arcs found!"<<std::endl;
        exit(-1);
    }

    //here we loop over the points of the boundary
    //we have to compute the weighted "center" of pixels
    double weights=0;
    double sum_x=0;
    double sum_y=0;
    double av_x;
    double av_y;

    for(size_t k=0; k<points[0].size(); k++)
    {
        double xx=points[0][k].x;
        double yy=points[0][k].y;
        double weight_factor=1;
        sum_x=sum_x+(xx*weight_factor);
        sum_y=sum_y+(yy*weight_factor);
        weights=weights+weight_factor;
    }

    // normalize
    av_x=sum_x/weights;
    av_y=sum_y/weights;

    float angle_front=atan2((av_y-junction.y),(av_x-junction.x));

    weights=0;
    sum_x=0;
    sum_y=0;

    for(size_t k=0; k<points[1].size(); k++)
    {
        double xx=points[1][k].x;
        double yy=points[1][k].y;
        double weight_factor=1;
        sum_x=sum_x+(xx*weight_factor);
        sum_y=sum_y+(yy*weight_factor);
        weights=weights+weight_factor;
    }

    // normalize
    av_x=sum_x/weights;
    av_y=sum_y/weights;

    float angle_back=atan2((av_y-junction.y),(av_x-junction.x));

    float angle_other=atan2((other.y-junction.y),(other.x-junction.x));

    //if angle of other arc is in between, add 2pi to smaller angle
    if(angle_front<angle_other && angle_other<angle_back) angle_front+=2.0f*PI;
    if(angle_back<angle_other && angle_other<angle_front) angle_back+=2.0f*PI;

    if (angle_front>angle_back) return(angle_front-angle_back);
    else return (angle_back-angle_front);
}

//in this version the number of pixels to conside it automatically determined
float calculate_angle2(std::vector< std::vector<point> > points, point junction, point other)
{
    if(points.size()!=2)
    {
        std::cout<<"Error: Two filled vectors needed!"<<std::endl;
        exit(-1);
    }

    if(other.x==-1)
    {
        std::cout<<"Error: No other arcs found!"<<std::endl;
        exit(-1);
    }

    std::vector <point> RANSAC1;
    std::vector <point> RANSAC2;

    if (points[0].size()>4) 
    {
        std::vector <float> xfit(2);
        std::vector <float> yfit(2);
        std::vector <float> residues;
        float median=0;
        float tolerance;
        float residue;

        for(int i=4;i<points[0].size();i++)
        {
            xfit[0]=((float)points[0][i].x-(float)points[0][0].x)/i;
            xfit[1]=(float)points[0][i].x-i*xfit[0];
            yfit[0]=((float)points[0][i].y-(float)points[0][0].y)/i;
            yfit[1]=(float)points[0][i].y-i*yfit[0];
            for(int j=0;j<=i;j++)
            {
                residues.push_back(sqrt((xfit[0]*j+xfit[1]-points[0][j].x)*(xfit[0]*j+xfit[1]-points[0][j].x)
                    +(yfit[0]*j+yfit[1]-points[0][j].y)*(yfit[0]*j+yfit[1]-points[0][j].y)));
            }
            std::sort(residues.begin(),residues.end());
            tolerance=exp(-(i+4)/sqrt(points[0].size()));

            if(i==4 || residues[floor((float)residues.size()/2)]<=median+tolerance)
            {
                median=residues[floor((float)residues.size()/2)];
                RANSAC1.clear();
                for(int j=0;j<=i;j++)
                {
                    residue=sqrt((xfit[0]*j+xfit[1]-points[0][j].x)*(xfit[0]*j+xfit[1]-points[0][j].x)
                        +(yfit[0]*j+yfit[1]-points[0][j].y)*(yfit[0]*j+yfit[1]-points[0][j].y));
                    if(residue<3) RANSAC1.push_back(points[0][j]);
                }
            }
            residues.clear();
        }
    }
    else for(int i=0;i<points[0].size();i++) RANSAC1.push_back(points[0][i]);

    if (points[1].size()>4)
    {
        std::vector <float> xfit(2);
        std::vector <float> yfit(2);
        std::vector <float> residues;
        float median=0;
        float tolerance;
        float residue;

        for(int i=4;i<points[1].size();i++)
        {
            xfit[0]=((float)points[1][i].x-(float)points[1][0].x)/i;
            xfit[1]=(float)points[1][i].x-i*xfit[0];
            yfit[0]=((float)points[1][i].y-(float)points[1][0].y)/i;
            yfit[1]=(float)points[1][i].y-i*yfit[0];
            for(int j=0;j<=i;j++)
            {
                residues.push_back(sqrt((xfit[0]*j+xfit[1]-points[1][j].x)*(xfit[0]*j+xfit[1]-points[1][j].x)
                    +(yfit[0]*j+yfit[1]-points[1][j].y)*(yfit[0]*j+yfit[1]-points[1][j].y)));
            }
            std::sort(residues.begin(),residues.end());
            tolerance=exp(-(i+4)/sqrt(points[1].size()));

            if(i==4 || residues[floor((float)residues.size()/2)]<=median+tolerance)
            {
                median=residues[floor((float)residues.size()/2)];
                RANSAC2.clear();
                for(int j=0;j<=i;j++)
                {
                    residue=sqrt((xfit[0]*j+xfit[1]-points[1][j].x)*(xfit[0]*j+xfit[1]-points[1][j].x)
                        +(yfit[0]*j+yfit[1]-points[1][j].y)*(yfit[0]*j+yfit[1]-points[1][j].y));
                    if(residue<3) RANSAC2.push_back(points[1][j]);
                }
            }
            residues.clear();
        }
    }
    else for(int i=0;i<points[1].size();i++) RANSAC2.push_back(points[1][i]);

    //here we loop over the points of the boundary
    //we have to compute the weighted "center" of pixels
    double weights=0;
    double sum_x=0;
    double sum_y=0;
    double av_x;
    double av_y;

    for(size_t k=0; k<RANSAC1.size(); k++)
    {
        double xx=RANSAC1[k].x;
        double yy=RANSAC1[k].y;
        double weight_factor=1;
        sum_x=sum_x+(xx*weight_factor);
        sum_y=sum_y+(yy*weight_factor);
        weights=weights+weight_factor;
    }

    // normalize
    av_x=sum_x/weights;
    av_y=sum_y/weights;

    float angle_front=atan2((av_y-junction.y),(av_x-junction.x));

    weights=0;
    sum_x=0;
    sum_y=0;

    for(size_t k=0; k<RANSAC2.size(); k++)
    {
        double xx=RANSAC2[k].x;
        double yy=RANSAC2[k].y;
        double weight_factor=1;
        sum_x=sum_x+(xx*weight_factor);
        sum_y=sum_y+(yy*weight_factor);
        weights=weights+weight_factor;
    }

    // normalize
    av_x=sum_x/weights;
    av_y=sum_y/weights;

    float angle_back=atan2((av_y-junction.y),(av_x-junction.x));

    float angle_other=atan2((other.y-junction.y),(other.x-junction.x));

    //if angle of other arc is in between, add 2pi to smaller angle
    if(angle_front<angle_other && angle_other<angle_back) angle_front+=2.0f*PI;
    if(angle_back<angle_other && angle_other<angle_front) angle_back+=2.0f*PI;

    if (angle_front>angle_back) return(angle_front-angle_back);
    else return (angle_back-angle_front);
}

int expb_f (const gsl_vector * x, void *data, gsl_vector * f)
{
    size_t n = ((struct data *)data)->n;
    double *y = ((struct data *)data)->y;

    double mu = gsl_vector_get(x,0);
    double sigma = gsl_vector_get(x,1);

    size_t i;

    for (i = 0; i < n; i++)
    {
        /* Model Yi = exp(-0.5f * (log(t)-mu)^2 / sigma^2 ) / (t*sigma*sqrt(2.0f*PI)) */
        double t = i+1;
        double Yi = exp(-0.5f * (log(t)-mu)*(log(t)-mu) / (sigma*sigma) ) / (t*sigma*sqrt(2.0f*3.141592654f));
        gsl_vector_set (f, i, (Yi - y[i]));
    }

    return GSL_SUCCESS;
}

int expb_df (const gsl_vector * x, void *data, gsl_matrix * J)
{
    size_t n = ((struct data *)data)->n;

    double mu = gsl_vector_get(x,0);
    double sigma = gsl_vector_get(x,1);

    size_t i;

    for (i = 0; i < n; i++)
    {
        /* Jacobian matrix J(i,j) = dfi / dxj,                                        */
        /* where fi = (Yi - yi),                                                      */
        /*       Yi = exp(-0.5f * (log(t)-mu)^2 / sigma^2 ) / (t*sigma*sqrt(2.0f*PI)) */
        /* and the xj are the parameters (mu,sigma)                                   */
        double t = i+1;
        double sigma2 = sigma*sigma;
        double e = exp(-0.5f * (log(t)-mu)*(log(t)-mu) / (sigma2) ) / (t*sqrt(2.0f*3.141592654f));
        gsl_matrix_set (J, i, 0, 2.0f * e * (log(t)-mu)/sigma ); 
        gsl_matrix_set (J, i, 1, exp(2*sigma2) * e * (4.0f-(1.0f/sigma2)) );
    }
    return GSL_SUCCESS;
}

int expb_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
    expb_f (x, data, f);
    expb_df (x, data, J);

    return GSL_SUCCESS;
}

void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
    printf ("iter: %3u x = % 3.2f % 3.2f |f(x)| = %g\n", (unsigned int) iter,
        gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), gsl_blas_dnrm2 (s->f));
}

float log_normal_fit(const size_t n, double * y, double * x_init, float & max_x, float & stdabw, float & mu, float & sigma)
{
    std::cout<<"log normal fit to histogram values"<<std::endl;

    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    int status;
    unsigned int i, iter = 0;
    const size_t p = 2;

    gsl_matrix *covar = gsl_matrix_alloc (p, p);

    struct data d = { n, y };
    gsl_multifit_function_fdf f;
    gsl_vector_view x = gsl_vector_view_array (x_init, p);
    const gsl_rng_type * type;
    gsl_rng * r;

    gsl_rng_env_setup();

    type = gsl_rng_default;
    r = gsl_rng_alloc (type);

    f.f = &expb_f;
    f.df = &expb_df;
    f.fdf = &expb_fdf;
    f.n = n;
    f.p = p;
    f.params = &d;

    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc (T, n, p);
    gsl_multifit_fdfsolver_set (s, &f, &x.vector);

    print_state (iter, s);

    do
    {
        iter++;
        status = gsl_multifit_fdfsolver_iterate (s);

        //printf ("status = %s\n", gsl_strerror (status));
        print_state (iter, s);

        if (status) break;

        status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
    }
    while (status == GSL_CONTINUE && iter < 500);

    gsl_multifit_covar (s->J, 0.0, covar);

    #define FIT(i) gsl_vector_get(s->x, i)
    #define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 

    //printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
    //printf ("mu    = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
    //printf ("sigma = %.5f +/- %.5f\n", FIT(1), c*ERR(1));

    //printf ("status = %s\n", gsl_strerror (status));

    mu=FIT(0);
    sigma=FIT(1);
    max_x=exp(FIT(0)-FIT(1)*FIT(1));
    stdabw=sqrt(exp(2*FIT(0)+FIT(1)*FIT(1))*(exp(FIT(1)*FIT(1))-1));

    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);
    gsl_rng_free (r);

    if(status==-2) return pow(chi, 2.0) / dof;
    else return -1.0f;
}

void linfit(std::vector<float> x, std::vector<float> y, float & m, float & b)
{
    float sumx = 0.0f, sumy = 0.0f, sumx2 = 0.0f, sumxy = 0.0f;
    float length = (float)y.size();

    if (y.size() <= 1)
    {
        m = 0.0f;
        b = 0.0f;
    }
    else
    {
        for(size_t i=0; i < y.size(); i++)
        {
            sumx += x[i];
            sumy += y[i];
            sumx2 += (x[i] * x[i]);
            sumxy += (x[i] * y[i]);
        }
        
        float divisor = (sumx2 - ((sumx * sumx) / length));

        if (divisor != 0.0f)
        {
            m = (sumxy - ((sumx * sumy) / length)) / divisor;
            b = (sumy - (m * sumx)) / length;
        }
        else
        {
            m = 0.0f;
            b = 0.0f;
        }
    }
}

void grain_correlation(std::vector< std::vector<int> > grain_areas,
                       std::vector< std::vector<float> > grain_box_flattening,
                       std::vector< std::vector<float> > grain_ellipse_flattening,
                       std::vector< std::vector<float> > parameter1,
                       std::vector< std::vector<float> > parameter2,
                       correlation & corr,
                       bool abs1=false,
                       bool abs2=false)
{
    //Calculate means and sums
    float mean1=0.0f;
    float mean2=0.0f;
    float sumx = 0.0f, sumy = 0.0f, sumx2 = 0.0f, sumxy = 0.0f;
    corr.nr=0;

    for (int i=0; i<grain_areas.size(); i++)
    {
        for (int area=0; area<grain_areas[i].size(); area++)
        {
            //ellipse not degenerated and exes found correctly
            if (grain_box_flattening[i][grain_areas[i][area]]>0.0f && grain_ellipse_flattening[i][grain_areas[i][area]]>0.0f)
            {
                if (abs1) parameter1[i][grain_areas[i][area]]=fabs(parameter1[i][grain_areas[i][area]]);
                if (abs2) parameter2[i][grain_areas[i][area]]=fabs(parameter2[i][grain_areas[i][area]]);

                mean1+=parameter1[i][grain_areas[i][area]];
                mean2+=parameter2[i][grain_areas[i][area]];

                sumx += parameter1[i][grain_areas[i][area]];
                sumy += parameter2[i][grain_areas[i][area]];
                sumx2 += (parameter1[i][grain_areas[i][area]] * parameter1[i][grain_areas[i][area]]);
                sumxy += (parameter1[i][grain_areas[i][area]] * parameter2[i][grain_areas[i][area]]);

                corr.nr++;
            }
        }
    }
    mean1/=corr.nr;
    mean2/=corr.nr;

    //Calculate standard deviations and covariance
    float standard_deviation1=0.0f;
    float standard_deviation2=0.0f;
    float covariance=0.0f;

    for (int i=0; i<grain_areas.size(); i++)
    {
        for (int area=0; area<grain_areas[i].size(); area++)
        {
            //ellipse not degenerated and exes found correctly
            if (grain_box_flattening[i][grain_areas[i][area]]>0.0f && grain_ellipse_flattening[i][grain_areas[i][area]]>0.0f)
            {   
                standard_deviation1+=square(parameter1[i][grain_areas[i][area]]-mean1);
                standard_deviation2+=square(parameter2[i][grain_areas[i][area]]-mean2);

                covariance+=(parameter1[i][grain_areas[i][area]]-mean1)*(parameter2[i][grain_areas[i][area]]-mean2);
            }
        }
    }
    standard_deviation1=sqrt(standard_deviation1*(1.0f/((float)corr.nr-1.0f)));
    standard_deviation2=sqrt(standard_deviation2*(1.0f/((float)corr.nr-1.0f)));

    //Calculate correlation
    corr.rho=covariance/((float)corr.nr*standard_deviation1*standard_deviation2);

    float divisor = (sumx2 - ((sumx * sumx) / (float)corr.nr));

    if (divisor != 0.0f)
    {
        corr.m = (sumxy - ((sumx * sumy) / (float)corr.nr)) / divisor;
        corr.b = (sumy - (corr.m * sumx)) / (float)corr.nr;
    }
    else
    {
        corr.m = 0.0f;
        corr.b = 0.0f;
    }
}

void grain_correlation(std::vector< std::vector<float> > parameter1,
                       std::vector< std::vector<float> > parameter2,
                       correlation & corr)
{
    //Calculate means and sums
    float mean1=0.0f;
    float mean2=0.0f;
    float sumx = 0.0f, sumy = 0.0f, sumx2 = 0.0f, sumxy = 0.0f;
    corr.nr=0;

    for (int i=0; i<parameter1.size(); i++)
    {
        for (int entry=0; entry<parameter1[i].size(); entry++)
        {
            mean1+=parameter1[i][entry];
            mean2+=parameter2[i][entry];

            sumx += parameter1[i][entry];
            sumy += parameter2[i][entry];
            sumx2 += (parameter1[i][entry] * parameter1[i][entry]);
            sumxy += (parameter1[i][entry] * parameter2[i][entry]);

            corr.nr++;
        }
    }
    mean1/=corr.nr;
    mean2/=corr.nr;

    //Calculate standard deviations and covariance
    float standard_deviation1=0.0f;
    float standard_deviation2=0.0f;
    float covariance=0.0f;

    for (int i=0; i<parameter1.size(); i++)
    {
        for (int entry=0; entry<parameter1[i].size(); entry++)
        {
            standard_deviation1+=square(parameter1[i][entry]-mean1);
            standard_deviation2+=square(parameter2[i][entry]-mean2);

            covariance+=(parameter1[i][entry]-mean1)*(parameter2[i][entry]-mean2);
        }
    }
    standard_deviation1=sqrt(standard_deviation1*(1.0f/((float)corr.nr-1.0f)));
    standard_deviation2=sqrt(standard_deviation2*(1.0f/((float)corr.nr-1.0f)));

    //Calculate correlation
    corr.rho=covariance/((float)corr.nr*standard_deviation1*standard_deviation2);

    float divisor = (sumx2 - ((sumx * sumx) / (float)corr.nr));

    if (divisor != 0.0f)
    {
        corr.m = (sumxy - ((sumx * sumy) / (float)corr.nr)) / divisor;
        corr.b = (sumy - (corr.m * sumx)) / (float)corr.nr;
    }
    else
    {
        corr.m = 0.0f;
        corr.b = 0.0f;
    }
}
