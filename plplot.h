/*! \file plplot.h
 * \brief Class to generate plots based om plplot library
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
#include "plstream.h"

#ifndef ROUND
#define ROUND( a )    (PLINT) ( ( a ) < 0. ? ( ( a ) - 0.5 ) : ( ( a ) + 0.5 ) )
#endif

class plplot
{
    public:

    void plfbox(PLFLT, PLFLT, float, float);
    void draw_histogram(std::string, std::string, std::string, std::vector<int>, float, float, float, float, float, std::string);
    void draw_histogram_shift(std::string, std::string, std::string, std::vector<int>, float, float, float, float, float, std::string);
    void draw_histogram_rose(std::string, std::string, std::string, std::vector<int>, float, float, float, float, float, std::string);
    void draw_histogram(std::string, std::string, std::string, std::vector<int>, float, float, float, std::string);
    void draw_histogram_numbers(std::string, std::string, std::string, std::vector<int>, float, float, float, float, float, std::string);
    void draw_histogram(std::string, std::string, std::string, std::vector< std::vector<int> >, float, float, float, std::string);
    void draw_multi_histogram(std::string, std::string, std::string, std::vector< std::vector<int> >, float, float, float, std::string);
    void draw_histogram_lognorm(std::string, std::string, std::string, std::vector<int>, float, float, float, float, float&, float&, float, std::string);
    void draw_histogram_left_lognorm(std::string, std::string, std::string, std::vector<int>, float, float, float, float, float&, float&, float&, float,
                                     float&, float&, std::string);
    void draw_histogram_log(std::string, std::string, std::string, std::vector<int>, float, float, float, float, float, std::string);
    void draw_histogram_log(std::string, std::string, std::string, std::vector< std::vector<int> >, float, float, std::vector<float>,
                            std::vector<float>, float, std::string);
    void draw_multi_histogram_log(std::string, std::string, std::string, std::vector< std::vector<int> >, float, float, float, std::string);
    void draw_histogram_log_log(std::string, std::string, std::string, std::vector<float>, float, float, float, float, float, float, std::string);
    void draw_values_errors(std::string, std::string, std::string, std::vector<float>, std::vector<float>, std::string);
    void draw_values_errors_lin_log(std::string, std::string, std::string, std::vector<float>, std::vector<float>, float, float, float, std::string);
    void draw_values_errors_log_lin(std::string, std::string, std::string, std::vector<float>, std::vector<float>, float, float, float, std::string);
    void draw_values_errors_log_log(std::string, std::string, std::string, std::vector<float>, std::vector<float>, float, float, float, float, std::string);
    void draw_values_errors_log_log_fit(std::string, std::string, std::string, std::vector<float>, std::vector<float>, float, float, float, float, std::string);
    void draw_depth(std::string, std::string, std::string, std::vector<float>, std::vector<float>, float, std::string);
    void draw_depth(std::string, std::string, std::string, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, float, std::string);
    void draw_depth(std::string, std::string, std::string, std::vector< std::vector<float> >, std::vector< std::vector<float> >, std::string);
    void draw_depth_nofit(std::string, std::string, std::string, std::vector< std::vector<float> >, std::vector< std::vector<float> >, std::string);
    void draw_depth_nofit2(std::string, std::string, std::string, std::vector< std::vector<float> >, std::vector< std::vector<float> >, std::string,
                          std::vector<float>, float);
    void draw_depth_errors(std::string, std::string, std::string, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>,
                           float, std::string);
    void draw_depth_errors(std::string, std::string, std::string, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>,
                           std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>,float, std::string);
    void draw_depth_3d(std::string, std::string, std::string, std::string, std::vector<float>, std::vector<float>, std::vector< std::vector<float> >,
                       std::string);
    void draw_depth_color(std::string, std::string, std::string, std::string, std::vector<float>, std::vector<float>, std::vector< std::vector<float> >,
                       std::string);
    void draw_curv_cross(std::string, std::string, std::string, std::vector<int>, std::vector<float>, std::string);
    void draw_correlation(std::string, std::string, std::vector< std::vector<int> >, std::vector< std::vector<float> >, std::vector< std::vector<float> >,
                          std::vector< std::vector<float> >, std::vector< std::vector<float> >, correlation, std::string, bool, bool);
    void draw_correlation(std::string, std::string, std::vector< std::vector<float> >, std::vector< std::vector<float> >, correlation, std::string);
    void draw_correlation(std::string, std::string, std::vector<float>, std::vector<float>, std::string);
    void draw_3d_blocks(std::string, std::string, std::string, std::string, std::vector<point3d>, int, int, std::string);
    void draw_3d_histogram_shift(std::string, std::string, std::string, std::string, std::vector<float>, std::vector< std::vector<int> >, float, float,
                                 std::string);

    private:

    plstream *pls;
    static PLFLT pos[], red[], green[], blue[];
    static PLINT black[];
};

//color space for axis (0-255)
PLINT plplot::black[] = {0};

//color space for values (0.0-1.0)
PLFLT plplot::pos[] = {0.0, 0.25, 0.5, 0.75, 1.0};
PLFLT plplot::red[] = {0.0, 0.25, 0.5, 1.0, 1.0};
PLFLT plplot::green[] = {1.0, 0.5, 0.5, 0.5, 1.0};
PLFLT plplot::blue[] = {1.0, 1.0, 0.5, 0.25, 0.0};
//PLFLT plplot::red[] = {0.0, 0.25, 0.5, 0.75, 1.0};
//PLFLT plplot::green[] = {1.0, 0.75, 0.5, 0.25, 0.0};
//PLFLT plplot::blue[] = {1.0, 0.75, 0.5, 0.25, 0.0};

void plplot::plfbox(PLFLT x0, PLFLT y0, float min, float width)//draw one box
{
    PLFLT *x = new PLFLT[5];
    PLFLT *y = new PLFLT[5];

    x[0] = x0;
    y[0] = min;
    x[1] = x0;
    y[1] = y0;
    x[2] = x0 + width;
    y[2] = y0;
    x[3] = x0 + width;
    y[3] = min;
    x[4] = x0;
    y[4] = min;
    pls->fill(5, x, y);
    pls->col0(0);
    pls->lsty(1);
    pls->line(5, x, y);

    delete x;
    delete y;
}

void plplot::draw_histogram(std::string x_axis, std::string y_axis, std::string title, std::vector<int> hist_values, float bin_width, float scaling,
                            float mean, float standard_deviation, float y_max, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    char mean_string[20];
    sprintf(mean_string, "%.3f", mean);

    char standard_deviation_string[20];
    sprintf(standard_deviation_string, "%.3f", standard_deviation);

    title.append(" (");
    title.append(mean_string);
    title.append("+-");
    title.append(standard_deviation_string);
    title.append(")");

    int x_max=hist_values.size();
    float width=bin_width/scaling;

    int sum_value=0;
    for(int i=0; i<hist_values.size(); i++)
    {
        sum_value+=hist_values[i];
    }

    for(int i=0; i<hist_values.size(); i++)
    {
        if(hist_values[i]>y_max*sum_value) y_max=(float)hist_values[i]/(float)sum_value;
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(-width, (x_max+1)*width, 0.0, y_max);//window size
    pls->box("binstv", 0.0, 0, "bcginstv", 0.0, 0);//box around window
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    for (int i = 0; i < hist_values.size(); i++)
    {
        pls->col1(0.1);//bar filler color
        pls->psty(0);//style of filler
        plfbox((float)i*width, (float)hist_values[i]/(float)sum_value, 0.0, width);//histogram values
    }

    delete pls;
}

void plplot::draw_histogram_shift(std::string x_axis, std::string y_axis, std::string title, std::vector<int> hist_values, float bin_width, float scaling,
                            float mean, float standard_deviation, float y_max, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    char mean_string[20];
    sprintf(mean_string, "%.3f", mean);

    char standard_deviation_string[20];
    sprintf(standard_deviation_string, "%.3f", standard_deviation);

    title.append(" (");
    title.append(mean_string);
    title.append("+-");
    title.append(standard_deviation_string);
    title.append(")");

    int x_mid=18;//hist_values.size()/2;
    float width=bin_width/scaling;

    int sum_value=0;
    for(int i=0; i<hist_values.size(); i++)
    {
        sum_value+=hist_values[i];
    }

    for(int i=0; i<hist_values.size(); i++)
    {
        if(hist_values[i]>y_max*sum_value) y_max=(float)hist_values[i]/(float)sum_value;
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(-(x_mid+1)*width, (hist_values.size()-x_mid+1)*width, 0.0, y_max);//window size
    pls->box("binstv", 30.0, 3, "bcginstv", 0.0, 0);//box around window
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    for (int i = 0; i < hist_values.size(); i++)
    {
        pls->col1(0.1);//bar filler color
        pls->psty(0);//style of filler
        plfbox((float)(i-x_mid)*width, (float)hist_values[i]/(float)sum_value, 0.0, width);//histogram values
    }

    delete pls;
}

void plplot::draw_histogram_rose(std::string x_axis, std::string y_axis, std::string title, std::vector<int> hist_values, float bin_width, float scaling,
                float mean, float standard_deviation, float y_max, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    char mean_string[20];
    sprintf(mean_string, "%.3f", mean);

    char standard_deviation_string[20];
    sprintf(standard_deviation_string, "%.3f", standard_deviation);

    title.append(" (");
    title.append(mean_string);
    title.append("+-");
    title.append(standard_deviation_string);
    title.append(")");

    int x_mid=18;//hist_values.size()/2;
    float width=bin_width/scaling;

    int sum_value=0;
    for(int i=0; i<hist_values.size(); i++)
    {
        sum_value+=hist_values[i];
    }

    for(int i=0; i<hist_values.size(); i++)
    {
        if(hist_values[i]>y_max*sum_value) y_max=(float)hist_values[i]/(float)sum_value;
    }

    PLFLT *x0 = new PLFLT[hist_values.size()+1];
    PLFLT *y0 = new PLFLT[hist_values.size()+1];
    PLFLT *x  = new PLFLT[4*hist_values.size()];
    PLFLT *y  = new PLFLT[4*hist_values.size()];

    for (int i = 0; i < hist_values.size()+1; i++)
    {
        x0[i] = cos(PI*(float)(i-x_mid)*width/180.0f);
        y0[i] = sin(PI*(float)(i-x_mid)*width/180.0f);
    }

    pls->scmap0(black,black,black,1);
    pls->col0(0);

    // Set up viewport and window, but do not draw box.
    pls->env( -1.3, 1.3, -1.3, 1.3, 1, -2 );

    // Draw circles for polar grid
    for (int i = 1; i <= 10; i++)
    {
        pls->arc( 0.0, 0.0, 0.1 * i, 0.1 * i, 0.0, 360.0, 0.0, 0 );
    }

    for (int i = 1; i <= 24; i++)
    {
        PLFLT dx = -cos(PI*i/12.0f);
        PLFLT dy = -sin(PI*i/12.0f);
        char  text[5];
        PLFLT theta, offset;

        // Draw radial spokes for polar grid.
        pls->join( 0.0, 0.0, dx, dy );
        if (i>=6 && i<=18) sprintf( text, "%d", (int) ROUND(15.0*i)-180 );
        else if (i<6) sprintf( text, "%d", (int) ROUND(15.0*i));
        else sprintf( text, "%d", (int) ROUND(15.0*i)-360);

        // Write labels for angle.
        if (theta < 9.99) offset = 0.45;
        else if (theta < 99.9) offset = 0.30;
        else offset = 0.15;

        //Slightly off zero to avoid floating point logic flips at 90 and 270 deg.
        if (dx >= -0.00001) pls->ptex(dx, dy, dx, dy, -offset, text);
        else pls->ptex(dx, dy, -dx, -dy, 1. + offset, text);
    }

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    // Draw the graph.
    for (int i = 0; i < hist_values.size(); i++)
    {
        int j=2*i;

        x[j] = 0.95f/y_max * x0[i] * (float)hist_values[i]/(float)sum_value;
        y[j] = 0.95f/y_max * y0[i] * (float)hist_values[i]/(float)sum_value;
        x[4*hist_values.size()-1-j] = -0.95f/y_max * x0[i] * (float)hist_values[i]/(float)sum_value;
        y[4*hist_values.size()-1-j] = -0.95f/y_max * y0[i] * (float)hist_values[i]/(float)sum_value;

        x[j+1] = 0.95f/y_max * x0[i+1] * (float)hist_values[i]/(float)sum_value;
        y[j+1] = 0.95f/y_max * y0[i+1] * (float)hist_values[i]/(float)sum_value;
        x[4*hist_values.size()-2-j] = -0.95f/y_max * x0[i+1] * (float)hist_values[i]/(float)sum_value;
        y[4*hist_values.size()-2-j] = -0.95f/y_max * y0[i+1] * (float)hist_values[i]/(float)sum_value;
    }

    pls->col1(0.1);
    pls->fill(4*hist_values.size(), x, y);
    pls->line(4*hist_values.size(), x, y);

    pls->col0(0);
    pls->mtex( "t", 1.0, 0.5, 0.5, title.c_str());

    delete x;
    delete y;
    delete x0;
    delete y0;

    delete pls;
}

void plplot::draw_histogram(std::string x_axis, std::string y_axis, std::string title, std::vector<int> hist_values, float bin_width, float scaling,
                            float y_max, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    int x_max=hist_values.size();
    float width=bin_width/scaling;

    int sum_value=0;
    for(int i=0; i<hist_values.size(); i++)
    {
        sum_value+=hist_values[i];
    }

    for(int i=0; i<hist_values.size(); i++)
    {
        if(hist_values[i]>y_max*sum_value) y_max=(float)hist_values[i]/(float)sum_value;
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(-width, (x_max+1)*width, 0.0, y_max);//window size
    pls->box("binstv", 0.0, 0, "bcginstv", 0.0, 0);//box around window
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    for (int i = 0; i < hist_values.size(); i++)
    {
        pls->col1(0.1);//bar filler color
        pls->psty(0);//style of filler
        plfbox((float)i*width, (float)hist_values[i]/(float)sum_value, 0.0, width);//histogram values
    }

    delete pls;
}

void plplot::draw_histogram_numbers(std::string x_axis, std::string y_axis, std::string title, std::vector<int> hist_values, float bin_width, float scaling,
                                    float mean, float standard_deviation, float y_max, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    char mean_string[20];
    sprintf(mean_string, "%.3f", mean);

    char standard_deviation_string[20];
    sprintf(standard_deviation_string, "%.3f", standard_deviation);

    title.append(" (");
    title.append(mean_string);
    title.append("+-");
    title.append(standard_deviation_string);
    title.append(")");

    int x_max=hist_values.size();
    float width=bin_width/scaling;

    int sum_value=0;
    for(int i=0; i<hist_values.size(); i++)
    {
        sum_value+=hist_values[i];
    }

    for(int i=0; i<hist_values.size(); i++)
    {
        if(hist_values[i]>y_max*sum_value) y_max=(float)hist_values[i]/(float)sum_value;
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(-width, (x_max+1)*width, 0.0, y_max);//window size
    pls->box("binstv", 0.0, 0, "bcginstv", 0.0, 0);//box around window
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    for (int i = 0; i < hist_values.size(); i++)
    {
        pls->col1(0.1);//bar filler color
        pls->psty(0);//style of filler
        plfbox((float)(i-0.5f)*width, (float)hist_values[i]/(float)sum_value, 0.0, width);//histogram values
    }

    delete pls;
}

void plplot::draw_histogram(std::string x_axis, std::string y_axis, std::string title, std::vector< std::vector<int> > hist_values, float bin_width,
                            float scaling, float y_max, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    int x_max=hist_values[0].size();
    float width=bin_width/scaling;

    int sum_value=0;
    for(int i=0; i<hist_values[1].size(); i++)
    {
        sum_value+=hist_values[1][i];
    }

    for(int i=0; i<hist_values[1].size(); i++)
    {
        if(hist_values[1][i]>y_max*sum_value) y_max=(float)hist_values[1][i]/(float)sum_value;
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(-width, (x_max+1)*width, 0.0, y_max);//window size
    //pls->box("binstv", 100.0, 10, "bcginstv", 0.0, 0);//box around window
    pls->box("binstv", 10.0, 10, "bcginstv", 0.0, 0);//box around window
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    for (int j=hist_values.size(); j>0; j--)
    {
        for (int i = 0; i < hist_values[0].size(); i++)
        {
            if(j==2) pls->col1(0.8f);//bar filler color
            else pls->col1(0.1f);//bar filler color
            pls->psty(0);//style of filler
            plfbox((float)i*width, (float)hist_values[j-1][i]/(float)sum_value, 0.0, width);//histogram values
        }
    }

    delete pls;
}

void plplot::draw_multi_histogram(std::string x_axis, std::string y_axis, std::string title, std::vector< std::vector<int> > hist_values, float bin_width,
                                  float scaling, float y_max, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    int x_max=0;
    float width=bin_width/scaling;
    std::vector<int> sum_value(hist_values.size(),0);

    for(int j=0; j<hist_values.size(); j++)
    {
        if (hist_values[j].size()>x_max) x_max=hist_values[j].size();

        int temp=0;
        for(int i=0; i<hist_values[j].size(); i++)
        {
            temp+=hist_values[j][i];
        }

        if(temp>sum_value[j]) sum_value[j]=temp;
    }

    for(int j=0; j<hist_values.size(); j++)
    {
        for(int i=0; i<hist_values[j].size(); i++)
        {
            if(hist_values[j][i]>y_max*sum_value[j]) y_max=(float)hist_values[j][i]/(float)sum_value[j];
        }
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(-width, (x_max+1)*width, 0.0, y_max);//window size
    pls->box("binstv", 0.0, 0, "bcginstv", 0.0, 0);//box around window
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    for (int j=hist_values.size(); j>0; j--)
    {
        for (int i=0; i<hist_values[j-1].size(); i++)
        {
            if(j==2) pls->col1(0.8f);//bar filler color
            else pls->col1(0.1f);//bar filler color
            pls->psty(0);//style of filler
            plfbox((float)(i+0.4*(j-1))*width, (float)hist_values[j-1][i]/(float)sum_value[j-1], 0.0, 0.4f*width);//histogram values
        }
    }

    delete pls;
}

void plplot::draw_histogram_lognorm(std::string x_axis, std::string y_axis, std::string title, std::vector<int> hist_values, float bin_width, float scaling,
                                    float start_mu, float start_sigma, float & max_x, float & stdabw, float y_max, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    //normalized histogram values
    double * y = new double[hist_values.size()];

    int x_max=hist_values.size();
    float width=bin_width/scaling;

    int sum_value=0;
    for(int i=0; i<hist_values.size(); i++)
    {
        sum_value+=hist_values[i];
    }

    for(int i=0; i<hist_values.size(); i++)
    {
        y[i]=(double)hist_values[i]/(double)sum_value;
        if(y[i]>y_max) y_max=y[i];
    }

    //do lognormal fit
    double x_init[2] = {start_mu, start_sigma};//start values for fit
    PLFLT * x = new PLFLT[5*hist_values.size()];
    PLFLT * fit = new PLFLT[5*hist_values.size()];
    float mu, sigma;
    float chisqdof = log_normal_fit(hist_values.size(), y, x_init, max_x, stdabw, mu, sigma);

    std::cout<<"chisq/dof = "<<chisqdof<<std::endl;
    std::cout<<"max x  = "<<max_x<<std::endl;
    std::cout<<"stdabw = "<<stdabw<<std::endl;
    std::cout<<"mu     = "<<mu<<std::endl;
    std::cout<<"sigma  = "<<sigma<<std::endl;

    //fit line to show in histogram
    for (int i = 0; i < 5*hist_values.size(); i++)
    {
        float t = (float)i/5.0f;
        x[i]=t;
        fit[i] = exp(-0.5f * (log(t)-mu)*(log(t)-mu) / (sigma*sigma) ) / (t*sigma*sqrt(2.0f*PI));
        //std::cout<<"data: "<<t<<" "<<y[(int)t]<<" "<<fit[i]<<std::endl;
        if(fit[i]>y_max) y_max=fit[i];
    }

    char x_mean[20];
    sprintf(x_mean, "%.3f", max_x);

    char standard_deviation_string[20];
    sprintf(standard_deviation_string, "%.3f", stdabw);

    title.append(" (");
    title.append(x_mean);
    title.append("+-");
    title.append(standard_deviation_string);
    title.append(")");

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(-width, (x_max+1)*width, 0.0, y_max);//window size
    pls->box("binstv", 0.0, 0, "bcginstv", 0.0, 0);//box around window
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    for (int i = 0; i < hist_values.size(); i++)
    {
        pls->col1(0.1);//bar filler color
        pls->psty(0);//style of filler
        plfbox((float)(i-0.5f)*width, y[i], 0.0, width);//histogram values
    }

    pls->col1(0.8);//line color
    pls->line(5*hist_values.size(),x,fit);

    delete y;
    delete x;
    delete fit;

    delete pls;
}

void plplot::draw_histogram_left_lognorm(std::string x_axis, std::string y_axis, std::string title, std::vector<int> hist_values, float bin_width,
                                    float scaling, float start_mu, float start_sigma, float & max_x, float & stdabw_low, float & stdabw_high, float y_max,
                                    float & mu, float & sigma, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    //normalized histogram values
    double * y = new double[hist_values.size()];

    int x_max=hist_values.size();
    float width=bin_width/scaling;

    int sum_value=0;
    for(int i=0; i<hist_values.size(); i++)
    {
        sum_value+=hist_values[i];
    }

    for(int i=0; i<hist_values.size(); i++)
    {
        y[i]=(double)hist_values[hist_values.size()-1-i]/(double)sum_value;
        if(y[i]>y_max) y_max=y[i];
    }

    //do lognormal fit
    int resolution=5;
    double x_init[2] = {start_mu, start_sigma};//start values for fit
    PLFLT * x = new PLFLT[(hist_values.size()+1)*resolution];
    PLFLT * fit = new PLFLT[(hist_values.size()+1)*resolution];
    float stdabw;
    float chisqdof = log_normal_fit(hist_values.size(), y, x_init, max_x, stdabw, mu, sigma);

    float temp=max_x;
    max_x=pow(10,(hist_values.size()-max_x+0.5f)/bin_width-log10(scaling));
    stdabw_low=max_x-pow(10,(hist_values.size()-temp-stdabw)/bin_width-log10(scaling));
    stdabw_high=pow(10,(hist_values.size()-temp+stdabw)/bin_width-log10(scaling))-max_x;

    std::cout<<"chisq/dof = "<<chisqdof<<std::endl;
    std::cout<<"max x = "<<max_x<<std::endl;
    std::cout<<"stdabw low = "<<stdabw_low<<std::endl;
    std::cout<<"stdabw high = "<<stdabw_high<<std::endl;
    std::cout<<"mu    = "<<mu<<std::endl;
    std::cout<<"sigma = "<<sigma<<std::endl;

    //fit line to show in histogram
    for (int i=0; i<hist_values.size()+1; i++)
    {
        for (int j=0; j<resolution; j++)
        {
            float t = (float)i + (float)j/(float)resolution;
            x[resolution*i+j] = (float)(hist_values.size()-t+0.5f)/bin_width-log10(scaling);
            fit[resolution*i+j] = exp(-0.5f * (log(t)-mu)*(log(t)-mu) / (sigma*sigma) ) / (t*sigma*sqrt(2.0f*PI));
            //std::cout<<"data: "<<x[resolution*i+j]<<" "<<y[i]<<" "<<fit[resolution*i+j]<<std::endl;
            if(fit[i]>y_max) y_max=fit[i];
        }
    }

    char x_mean[20];
    sprintf(x_mean, "%.3f", max_x);

    char standard_low[20];
    sprintf(standard_low, "%.3f", stdabw_low);

    char standard_high[20];
    sprintf(standard_high, "%.3f", stdabw_high);

    title.append(" (");
    title.append(x_mean);
    title.append(" -");
    title.append(standard_low);
    title.append("/+");
    title.append(standard_high);
    title.append(")");

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(-log10(scaling), (float)(x_max+1)/bin_width-log10(scaling), 0.0, y_max);//window size
    pls->box("bilnstv", 0.0, 0, "bcginstv", 0.0, 0);//box around window
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    for (int i = 0; i < hist_values.size(); i++)
    {
        pls->col1(0.1);//bar filler color
        pls->psty(0);//style of filler
        plfbox((float)i/bin_width-log10(scaling), (float)hist_values[i]/(float)sum_value, 0.0, 1.0f/bin_width);//histogram values
    }

    pls->col1(0.8);//line color
    pls->line((hist_values.size()+1)*resolution,x,fit);

    delete y;
    delete x;
    delete fit;

    delete pls;
}

void plplot::draw_histogram_log(std::string x_axis, std::string y_axis, std::string title, std::vector<int> hist_values, float bin_width, float scaling,
                                float mean, float standard_deviation, float y_max, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    char mean_string[20];
    sprintf(mean_string, "%.3f", mean);

    char standard_deviation_string[20];
    sprintf(standard_deviation_string, "%.3f", standard_deviation);

    title.append(" (");
    title.append(mean_string);
    title.append("+-");
    title.append(standard_deviation_string);
    title.append(")");

    int x_max=hist_values.size();

    int sum_value=0;
    for(int i=0; i<hist_values.size(); i++)
    {
        sum_value+=hist_values[i];
    }

    for(int i=0; i<hist_values.size(); i++)
    {
        if(hist_values[i]>y_max*sum_value) y_max=(float)hist_values[i]/(float)sum_value;
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(-log10(scaling), (float)(x_max+1)/bin_width-log10(scaling), 0.0, y_max);//window size
    pls->box("bilnstv", 0.0, 0, "bcginstv", 0.0, 0);//box around window
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    for (int i = 0; i < hist_values.size(); i++)
    {
        pls->col1(0.1);//bar filler color
        pls->psty(0);//style of filler
        plfbox((float)i/bin_width-log10(scaling), (float)hist_values[i]/(float)sum_value, 0.0, 1.0f/bin_width);//histogram values
    }

    delete pls;
}

void plplot::draw_histogram_log(std::string x_axis, std::string y_axis, std::string title, std::vector< std::vector<int> > hist_values, float bin_width,
                                float scaling, std::vector<float> mean, std::vector<float> standard_deviation, float y_max, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    char mean_string[20];
    sprintf(mean_string, "%.3f", mean.back());

    char standard_deviation_string[20];
    sprintf(standard_deviation_string, "%.3f", standard_deviation.back());

    title.append(" (");
    title.append(mean_string);
    title.append("+-");
    title.append(standard_deviation_string);
    title.append(")");

    int x_max=hist_values[0].size();

    int sum_value=0;
    for(int j=0; j<hist_values.size(); j++)
    {
        int temp=0;
        for(int i=0; i<hist_values[j].size(); i++)
        {
            temp+=hist_values[j][i];
        }

        if(temp>sum_value) sum_value=temp;
    }

    for(int j=0; j<hist_values.size(); j++)
    {
        for(int i=0; i<hist_values[j].size(); i++)
        {
            if(hist_values[j][i]>y_max*sum_value) y_max=(float)hist_values[j][i]/(float)sum_value;
        }
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(-log10(scaling), (float)(x_max+1)/bin_width-log10(scaling), 0.0, y_max);//window size
    pls->box("bilnstv", 0.0, 0, "bcginstv", 0.0, 0);//box around window
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    for (int j=hist_values.size(); j>0; j--)
    {
        for (int i=0; i<hist_values[j-1].size(); i++)
        {
            pls->col1(0.1*(j*10/(int)hist_values.size()-1));//bar filler color
            pls->psty(0);//style of filler
            plfbox((float)i/bin_width-log10(scaling), (float)hist_values[j-1][i]/(float)sum_value, 0.0, 1.0f/bin_width);//histogram values
        }
    }

    delete pls;
}

void plplot::draw_multi_histogram_log(std::string x_axis, std::string y_axis, std::string title, std::vector< std::vector<int> > hist_values,
                                      float bin_width, float scaling,  float y_max, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    int x_max=0;
    std::vector<int> sum_value(hist_values.size(),0);

    for(int j=0; j<hist_values.size(); j++)
    {
        if (hist_values[j].size()>x_max) x_max=hist_values[j].size();

        int temp=0;
        for(int i=0; i<hist_values[j].size(); i++)
        {
            temp+=hist_values[j][i];
        }

        if(temp>sum_value[j]) sum_value[j]=temp;
    }

    for(int j=0; j<hist_values.size(); j++)
    {
        for(int i=0; i<hist_values[j].size(); i++)
        {
            if(hist_values[j][i]>y_max*sum_value[j]) y_max=(float)hist_values[j][i]/(float)sum_value[j];
        }
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(-log10(scaling), (float)(x_max+1)/bin_width-log10(scaling), 0.0, y_max);//window size
    pls->box("bilnstv", 0.0, 0, "bcginstv", 0.0, 0);//box around window
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    for (int j=hist_values.size(); j>0; j--)
    {
        for (int i=0; i<hist_values[j-1].size(); i++)
        {
            if(j==2) pls->col1(0.8f);//bar filler color
            else pls->col1(0.1f);//bar filler color
            pls->psty(0);//style of filler
            plfbox((float)(i+0.4*(j-1))/bin_width-log10(scaling), (float)hist_values[j-1][i]/(float)sum_value[j-1], 0.0, 0.4f/bin_width);//histogram values
        }
    }

    delete pls;
}

void plplot::draw_histogram_log_log(std::string x_axis, std::string y_axis, std::string title, std::vector<float> hist_values, float bin_width,
                                    float scaling, float x_min, float mean, float standard_deviation, float y_max, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    char mean_string[20];
    sprintf(mean_string, "%.2e", mean);

    char standard_deviation_string[20];
    sprintf(standard_deviation_string, "%.2e", standard_deviation);

    title.append(" (");
    title.append(mean_string);
    title.append("+-");
    title.append(standard_deviation_string);
    title.append(")");

    float x_max=hist_values.size()/bin_width+x_min;

    float sum_value=0;
    for(int i=0; i<hist_values.size(); i++)
    {
        sum_value+=hist_values[i];
    }

    float y_min=0.0f;
    for(int i=0; i<hist_values.size(); i++)
    {
        if(hist_values[i]>0)
        {
            float log_value=log10(hist_values[i]/sum_value);
            if(log_value>y_max) y_max=log_value;
            if(log_value<y_min) y_min=log_value;
        }
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(x_min-(1.0f/bin_width)-log10(scaling), x_max+(1.0f/bin_width)-log10(scaling), y_min-1, y_max);//window size
    pls->box("bilnstv", 0.0, 0, "bcgilnstv", 0.0, 0);//box around window
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    for (int i = 0; i < hist_values.size(); i++)
    {
        pls->col1(0.1);//bar filler color
        //pls->col1(std::min(20.0f,(float)i)/20.0f);//bar filler color
        pls->psty(0);//style of filler
        plfbox(x_min + (float)i/bin_width-log10(scaling), log10(hist_values[i]/sum_value), y_min-1, 1.0f/bin_width);//histogram values
    }

    delete pls;
}

void plplot::draw_values_errors(std::string x_axis, std::string y_axis, std::string title, std::vector<float> values, std::vector<float> errors, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    PLFLT * y_min = new PLFLT[values.size()];
    PLFLT * y_max = new PLFLT[values.size()];
    PLFLT * y = new PLFLT[values.size()];
    PLFLT * x = new PLFLT[values.size()];

    float y_high=0;

    for(int i=0; i<values.size(); i++)
    {
        y_min[i]=values[i]-errors[i];
        y_max[i]=values[i]+errors[i];
        y[i]=values[i];
        x[i]=i+1;
        if(y_max[i]*1.1>y_high) y_high=y_max[i]*1.1;
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind( 0.0, values.size()+1, 0.0, y_high);
    pls->box( "bcinst", 2.0, 2.0, "bcinstv", 0, 0 );
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,red,red,NULL);

    pls->col1(0.0);
    pls->poin( values.size(), x , y , 0.0 );
    pls->erry( values.size(), x, y_min, y_max );

    delete y_min;
    delete y_max;
    delete y;
    delete x;

    delete pls;
}

void plplot::draw_values_errors_lin_log(std::string x_axis, std::string y_axis, std::string title, std::vector<float> values, std::vector<float> errors,
                                        float x_bin_width, float y_bin_width, float y_scaling, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    PLFLT * y_min = new PLFLT[values.size()];
    PLFLT * y_max = new PLFLT[values.size()];
    PLFLT * y = new PLFLT[values.size()];
    PLFLT * x = new PLFLT[values.size()];

    float y_high=0.0f;
    for(int i=0; i<values.size(); i++) if (values[i]+errors[i]>y_high) y_high=values[i]+errors[i];
    float y_low=y_high;
    for(int i=0; i<values.size(); i++) if (values[i]-errors[i]<y_low) y_low=values[i]-errors[i];

    float x_min=values.size();

    for(int i=0; i<values.size(); i++)
    {
        y_min[i]=(values[i]-errors[i])/y_bin_width-log10(y_scaling);
        y_max[i]=(values[i]+errors[i])/y_bin_width-log10(y_scaling);
        y[i]=values[i]/y_bin_width-log10(y_scaling);
        x[i]=(float)(i+0.5)*x_bin_width;
        if(values[i]>0.0f && i<x_min) x_min=(float)(i-1);
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind( x_min*x_bin_width, float(values.size()+1)*x_bin_width, (float)(y_low-1)/y_bin_width-log10(y_scaling), (float)(y_high+1)/y_bin_width-log10(y_scaling));
    pls->box( "bcinst", 0.0, 0.0, "bcilnstv", 0.0, 0.0 );
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());
    pls->scmap1l(true,5,pos,red,red,red,NULL);

    pls->col1(0.0);
    pls->poin( values.size(), x , y , 0.0 );
    pls->erry( values.size(), x, y_min, y_max );

    delete y_min;
    delete y_max;
    delete y;
    delete x;

    delete pls;
}

void plplot::draw_values_errors_log_lin(std::string x_axis, std::string y_axis, std::string title, std::vector<float> values, std::vector<float> errors,
                                        float x_bin_width, float x_scaling, float y_bin_width, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    PLFLT * y_min = new PLFLT[values.size()];
    PLFLT * y_max = new PLFLT[values.size()];
    PLFLT * y = new PLFLT[values.size()];
    PLFLT * x = new PLFLT[values.size()];

    float y_high=0.0f;
    for(int i=0; i<values.size(); i++) if (values[i]+errors[i]>y_high) y_high=values[i]+errors[i];
    float y_low=y_high;
    for(int i=0; i<values.size(); i++) if (values[i]-errors[i]<y_low) y_low=values[i]-errors[i];

    int x_max=values.size();

    for(int i=0; i<values.size(); i++)
    {
        y_min[i]=(values[i]-errors[i])*y_bin_width;
        y_max[i]=(values[i]+errors[i])*y_bin_width;
        y[i]=values[i]*y_bin_width;
        x[i]=(float)(i+0.5)/x_bin_width-log10(x_scaling);
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(-log10(x_scaling), (float)(x_max+1)/x_bin_width-log10(x_scaling), (float)(y_low-1)*y_bin_width, (float)(y_high+1)*y_bin_width);
    pls->box( "bcilnst", 0.0, 0.0, "bcinstv", 0.0, 0.0 );
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());
    pls->scmap1l(true,5,pos,red,red,red,NULL);

    pls->col1(0.0);
    pls->poin( values.size(), x , y , 0.0 );
    pls->erry( values.size(), x, y_min, y_max );

    delete y_min;
    delete y_max;
    delete y;
    delete x;

    delete pls;
}

void plplot::draw_values_errors_log_log(std::string x_axis, std::string y_axis, std::string title, std::vector<float> values, std::vector<float> errors,
                                        float x_bin_width, float x_scaling, float y_bin_width, float y_scaling, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    PLFLT * y_min = new PLFLT[values.size()];
    PLFLT * y_max = new PLFLT[values.size()];
    PLFLT * y = new PLFLT[values.size()];
    PLFLT * x = new PLFLT[values.size()];

    float y_high=0.0f;
    for(int i=0; i<values.size(); i++) if (values[i]+errors[i]>y_high) y_high=values[i]+errors[i];
    float y_low=y_high;
    for(int i=0; i<values.size(); i++) if (values[i]-errors[i]<y_low) y_low=values[i]-errors[i];

    int x_max=values.size();

    for(int i=0; i<values.size(); i++)
    {
        y_min[i]=(values[i]-errors[i])/y_bin_width-log10(y_scaling);
        y_max[i]=(values[i]+errors[i])/y_bin_width-log10(y_scaling);
        y[i]=values[i]/y_bin_width-log10(y_scaling);
        x[i]=(float)(i+0.5)/x_bin_width-log10(x_scaling);
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(-log10(x_scaling), (float)(x_max+1)/x_bin_width-log10(x_scaling), (float)(y_low-1)/y_bin_width-log10(y_scaling),
              (float)(y_high+1)/y_bin_width-log10(y_scaling));
    pls->box( "bcilnst", 0.0, 0.0, "bcilnstv", 0.0, 0.0 );
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());
    pls->scmap1l(true,5,pos,red,red,red,NULL);

    pls->col1(0.0);
    pls->poin( values.size(), x , y , 0.0 );
    pls->erry( values.size(), x, y_min, y_max );

    delete y_min;
    delete y_max;
    delete y;
    delete x;

    delete pls;
}

void plplot::draw_values_errors_log_log_fit(std::string x_axis, std::string y_axis, std::string title, std::vector<float> values, std::vector<float> errors,
                                        float x_bin_width, float x_scaling, float y_bin_width, float y_scaling, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    PLFLT * y_min = new PLFLT[values.size()];
    PLFLT * y_max = new PLFLT[values.size()];
    PLFLT * y = new PLFLT[values.size()];
    PLFLT * x = new PLFLT[values.size()];

    float y_high=0.0f;
    for(int i=0; i<values.size(); i++) if (values[i]+errors[i]>y_high) y_high=values[i]+errors[i];
    float y_low=y_high;
    for(int i=0; i<values.size(); i++) if (values[i]-errors[i]<y_low) y_low=values[i]-errors[i];

    int x_max=values.size();
    std::vector<float> xx, yy;

    for(int i=0; i<values.size(); i++)
    {
        y_min[i]=(values[i]-errors[i])/y_bin_width-log10(y_scaling);
        y_max[i]=(values[i]+errors[i])/y_bin_width-log10(y_scaling);
        y[i]=values[i]/y_bin_width-log10(y_scaling);
        x[i]=(float)(i+0.5)/x_bin_width-log10(x_scaling);

        if(values[i]>0.0f)
        {
            yy.push_back(values[i]/y_bin_width-log10(y_scaling));
            xx.push_back((float)(i+0.5)/x_bin_width-log10(x_scaling));
        }
    }

    float m, b;
    linfit(xx,yy,m,b);
    char fit[20];
    sprintf(fit, " (y = %.2f * x^ %.2f)", pow(10.0f,b), m);
    title.append(fit);

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(-log10(x_scaling), (float)(x_max+1)/x_bin_width-log10(x_scaling), (float)(y_low-1)/y_bin_width-log10(y_scaling),
              (float)(y_high+1)/y_bin_width-log10(y_scaling));
    pls->box( "bcilnst", 0.0, 0.0, "bcilnstv", 0.0, 0.0 );
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());
    pls->scmap1l(true,5,pos,red,red,red,NULL);

    pls->col1(0.0);
    pls->poin( values.size(), x , y , 0.0 );
    pls->erry( values.size(), x, y_min, y_max );

    pls->scmap1l(true,5,pos,red,green,blue,NULL);
    pls->col1(0.8);

    delete y;
    delete x;
    y = new PLFLT[2];
    x = new PLFLT[2];

    y[0]=-m*log10(x_scaling)+b;
    x[0]=-log10(x_scaling);
    y[1]=m*((float)(x_max+1)/x_bin_width-log10(x_scaling))+b;
    x[1]=(float)(x_max+1)/x_bin_width-log10(x_scaling);

    pls->lsty(1);
    pls->line(2, x, y);

    delete y_min;
    delete y_max;
    delete y;
    delete x;

    delete pls;
}

void plplot::draw_depth(std::string x_axis, std::string y_axis, std::string title, std::vector<float> x_values, std::vector<float> y_values,
                        float y_minimal, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    PLFLT * y = new PLFLT[y_values.size()];
    PLFLT * x = new PLFLT[y_values.size()];

    float y_high=0;
    float x_max=0;

    for(int i=0; i<y_values.size(); i++)
    {
        y[i]=y_values[i];
        x[i]=x_values[i];
        if(y[i]*1.1>y_high) y_high=y[i]*1.1;
        if(x[i]*1.05>x_max) x_max=x[i]*1.05;
        if(y[i]<0.0f && y[i]*1.1<y_minimal) y_minimal=y[i]*1.1;
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind( 0.0, x_max, y_minimal, y_high);
    if (x_max>500) pls->box( "bcinst", 250.0, 5.0, "bcinstv", 0, 0 );
    else pls->box( "bcinst", 50.0, 5.0, "bcinstv", 0, 0 );
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,red,red,NULL);

    pls->col1(0.0);
    pls->poin( y_values.size(), x , y , 0.0 );

    delete y;
    delete x;

    delete pls;
}

void plplot::draw_depth(std::string x_axis, std::string y_axis, std::string title, std::vector<float> x_values1, std::vector<float> y_values1,
                        std::vector<float> x_values2, std::vector<float> y_values2, float y_minimal, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    PLFLT * y1 = new PLFLT[y_values1.size()];
    PLFLT * x1 = new PLFLT[y_values1.size()];
    PLFLT * y2 = new PLFLT[y_values2.size()];
    PLFLT * x2 = new PLFLT[y_values2.size()];

    float y_high=0;
    float x_max=0;

    for(int i=0; i<y_values1.size(); i++)
    {
        y1[i]=y_values1[i];
        x1[i]=x_values1[i];
        if(y1[i]*1.1>y_high) y_high=y1[i]*1.1;
        if(x1[i]*1.05>x_max) x_max=x1[i]*1.05;
        if(y1[i]<0.0f && y1[i]*1.1<y_minimal) y_minimal=y1[i]*1.1;
    }

    for(int i=0; i<y_values2.size(); i++)
    {
        y2[i]=y_values2[i];
        x2[i]=x_values2[i];
        if(y2[i]*1.1>y_high) y_high=y2[i]*1.1;
        if(x2[i]*1.05>x_max) x_max=x2[i]*1.05;
        if(y2[i]<0.0f && y2[i]*1.1<y_minimal) y_minimal=y2[i]*1.1;
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind( 0.0, x_max, y_minimal, y_high);
    if (x_max>500) pls->box( "bcinst", 250.0, 5.0, "bcinstv", 0, 0 );
    else pls->box( "bcinst", 50.0, 5.0, "bcinstv", 0, 0 );
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    pls->col1(0.1);
    pls->poin( y_values1.size(), x1 , y1 , 0.0 );

    delete y1;
    delete x1;

    pls->col1(0.8);
    pls->poin( y_values2.size(), x2 , y2 , 0.0 );

    delete y2;
    delete x2;

    delete pls;
}

void plplot::draw_depth(std::string x_axis, std::string y_axis, std::string title, std::vector< std::vector<float> > x_values,
                        std::vector< std::vector<float> > y_values, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    float y_high=0;
    float x_min=1000.0f;
    float x_max=0;

    for (int j=0; j<y_values.size(); j++)
    {
        for(int i=0; i<y_values[j].size(); i++)
        {
            if(y_values[j][i]*1.1>y_high) y_high=y_values[j][i]*1.1;
            if(x_values[j][i]<x_min) x_min=x_values[j][i];
            if(x_values[j][i]*1.05>x_max) x_max=x_values[j][i]*1.05;
        }
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind( x_min-0.09*x_max, x_max, 0.0, y_high);
    if (x_max>500) pls->box( "bcinst", 250.0, 5.0, "bcinstv", 0, 0 );
    else pls->box( "bcinst", 50.0, 5.0, "bcinstv", 0, 0 );
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    for (int j=0; j<y_values.size(); j++)
    {
        PLFLT * y = new PLFLT[y_values[j].size()];
        PLFLT * x = new PLFLT[y_values[j].size()];
        std::vector<float> xx, yy;

        for (int i=0; i<y_values[j].size(); i++)
        {
            y[i]=y_values[j][i];
            x[i]=x_values[j][i];
            yy.push_back(y[i]);
            xx.push_back(x[i]);
        }

        //pls->col1(0.1*((j+1)*10/(int)y_values.size()-1));
        if(j==0) pls->col1(0.7);
        if(j==1) pls->col1(0.4);
        if(j==2) pls->col1(0.1);
        if(j==3) pls->col1(0.8);
        pls->poin( y_values[j].size(), x , y , j+2 );
        float m, b;
        linfit(xx,yy,m,b);
        //std::cout<<"Slope is "<<m<<", intercept is "<<b<<std::endl;

        for (int i=0; i<y_values[j].size(); i++)
        {
            x[i]=x_values[j][i];
            y[i]=m * x[i] + b;
        }

        pls->lsty(1);
        pls->line(y_values[j].size(), x, y);

        delete y;
        delete x;
    }

    delete pls;
}

void plplot::draw_depth_nofit(std::string x_axis, std::string y_axis, std::string title, std::vector< std::vector<float> > x_values,
                        std::vector< std::vector<float> > y_values, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    float y_high=0;
    float x_min=1000.0f;
    float x_max=0;

    for (int j=0; j<y_values.size(); j++)
    {
        for(int i=0; i<y_values[j].size(); i++)
        {
            if(y_values[j][i]*1.1>y_high) y_high=y_values[j][i]*1.1;
            if(x_values[j][i]<x_min) x_min=x_values[j][i];
            if(x_values[j][i]*1.05>x_max) x_max=x_values[j][i]*1.05;
        }
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind( x_min-0.09*x_max, x_max, 0.0, y_high);
    if (x_max>500) pls->box( "bcinst", 250.0, 5.0, "bcinstv", 0, 0 );
    else pls->box( "bcinst", 50.0, 5.0, "bcinstv", 0, 0 );
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    for (int j=0; j<y_values.size(); j++)
    {
        PLFLT * y = new PLFLT[y_values[j].size()];
        PLFLT * x = new PLFLT[y_values[j].size()];

        for (int i=0; i<y_values[j].size(); i++)
        {
            y[i]=y_values[j][i];
            x[i]=x_values[j][i];
        }

        //pls->col1(0.1*((j+1)*10/(int)y_values.size()-1));
        if(j==0) pls->col1(0.7);
        if(j==1) pls->col1(0.4);
        if(j==2) pls->col1(0.1);
        if(j==3) pls->col1(0.8);
        pls->poin( y_values[j].size(), x , y , j+2 );

        delete y;
        delete x;
    }

    delete pls;
}

void plplot::draw_depth_nofit2(std::string x_axis, std::string y_axis, std::string title, std::vector< std::vector<float> > x_values,
                        std::vector< std::vector<float> > y_values, std::string dest_path, std::vector<float> errors, float y_minimal=0.0f)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    float y_high=0;
    float x_min=1000.0f;
    float x_max=0;

    for (int j=0; j<y_values.size(); j++)
    {
        for(int i=0; i<y_values[j].size(); i++)
        {
            if(y_values[j][i]*1.1>y_high) y_high=y_values[j][i]*1.1;
            if(x_values[j][i]<x_min) x_min=x_values[j][i];
            if(x_values[j][i]*1.05>x_max) x_max=x_values[j][i]*1.05;
        }
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
//    pls->wind( x_min-0.09*x_max, x_max, y_minimal, y_high);
    pls->wind(0.0, x_max, y_minimal, y_high);
    if (x_max>4000) pls->box( "bcinst", 0.0, 0, "bcinstv", 0, 0 );
    else if (x_max>500) pls->box( "bcinst", 250.0, 5.0, "bcinstv", 0, 0 );
    else pls->box( "bcinst", 50.0, 5.0, "bcinstv", 0, 0 );
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,red,red,NULL);

    int XPTS = 100;
    int YPTS = 100;
    PLFLT **z = new PLFLT *[XPTS];
    for (int i = 0; i < XPTS; i++ )
    {
        z[i] = new PLFLT[YPTS];
        for(int j = 0; j < YPTS; j++) z[i][j] = 0.0;
    }

    for(int i=0; i<y_values[2].size(); i++)
    {
        int xx=XPTS*x_values[2][i]/x_max;

        PLFLT y_min=y_values[2][i]-errors[i];
        PLFLT y_max=y_values[2][i]+errors[i];

        int yy1=YPTS*(std::max(y_minimal,(float)y_min)-y_minimal)/(y_high-y_minimal);
        int yy2=YPTS*(std::min(y_high,(float)y_max)-y_minimal)/(y_high-y_minimal);

        for(int yy = yy1; yy <= yy2 && yy<YPTS; yy++)
        {
            z[xx][yy] = 1.0;
            if(xx>0) z[xx-1][yy] = 1.0;
            if(xx<XPTS-1) z[xx+1][yy] = 1.0;
        }
    }

    PLFLT shade_min = 0.9, shade_max = 1.1, sh_color = 0.65;
    int   sh_cmap   = 1, sh_width = 0;
    int   min_color = 0, min_width = 0, max_color = 0, max_width = 0;

    pls->psty(0);
    pls->shade(z, XPTS, YPTS, NULL, 0, (PLFLT) x_max, (PLFLT) y_minimal, (PLFLT) y_high,
        shade_min, shade_max, sh_cmap, sh_color, sh_width, min_color, min_width, max_color, max_width,
        plstream::fill, false, NULL, NULL);

    for (int j=0; j<y_values.size(); j++)
    {
        PLFLT * y = new PLFLT[y_values[j].size()];
        PLFLT * x = new PLFLT[y_values[j].size()];

        for (int i=0; i<y_values[j].size(); i++)
        {
            y[i]=y_values[j][i];
            x[i]=x_values[j][i];
        }

        if(j==2)
        {
            pls->scmap1l(true,5,pos,red,red,red,NULL);
            pls->col1(0.0);
            pls->poin( y_values[j].size(), x , y , 0.0 );
        }
        else
        {
            pls->scmap1l(true,5,pos,red,green,blue,NULL);
            if(j==0) pls->col1(0.8);
            if(j==1) pls->col1(0.4);
            if(j==0) pls->poin( y_values[j].size(), x , y , 2 );
            if(j==1) pls->poin( y_values[j].size(), x , y , 5 );
        }

        delete y;
        delete x;
    }

    delete pls;
}

void plplot::draw_depth_errors(std::string x_axis, std::string y_axis, std::string title, std::vector<float> x_values,
                               std::vector<float> y_values, std::vector<float> y_errors_low, std::vector<float> y_errors_high,
                               float y_minimal, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    PLFLT * y_min = new PLFLT[y_values.size()];
    PLFLT * y_max = new PLFLT[y_values.size()];
    PLFLT * y = new PLFLT[y_values.size()];
    PLFLT * x = new PLFLT[y_values.size()];

    float y_high=0;
    float x_max=0;

    for(int i=0; i<y_values.size(); i++)
    {
        y_min[i]=y_values[i]-y_errors_low[i];
        y_max[i]=y_values[i]+y_errors_high[i];
        y[i]=y_values[i];
        x[i]=x_values[i];
        if(y_max[i]>y_high)
        {
            if (y_max[i]<1.5f*y[i]) y_high=std::max(y_high,(float)y_max[i]*1.1f);
            else y_high=std::max(y_high,(float)y[i]*1.5f);
            if (y[i]<0.0f)
            {
                if (y_min[i]>1.5f*y[i]) y_minimal=std::min(y_minimal,(float)y_min[i]*1.1f);
                else y_minimal=std::min(y_minimal,(float)y[i]*1.5f);
            }
        }
        if(x[i]*1.05>x_max) x_max=x[i]*1.05;
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind( 0.0, x_max, y_minimal, y_high);
    if (x_max>500) pls->box( "bcinst", 250.0, 5.0, "bcinstv", 0, 0 );
    else if (x_max>50) pls->box( "bcinst", 50.0, 5.0, "bcinstv", 0, 0 );
    else pls->box( "bcinst", 0, 0, "bcinstv", 0, 0 );
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,red,red,NULL);

    pls->col1(0.0);
    pls->poin( y_values.size(), x , y , 0.0 );
    pls->erry( y_values.size(), x, y_min, y_max );

    delete y_min;
    delete y_max;
    delete y;
    delete x;

    delete pls;
}

void plplot::draw_depth_errors(std::string x_axis, std::string y_axis, std::string title, std::vector<float> x_values1,
                               std::vector<float> y_values1, std::vector<float> y_errors1_low, std::vector<float> y_errors1_high,
                               std::vector<float> x_values2, std::vector<float> y_values2, std::vector<float> y_errors2_low,
                               std::vector<float> y_errors2_high, float y_minimal, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    PLFLT * y_min1 = new PLFLT[y_values1.size()];
    PLFLT * y_max1 = new PLFLT[y_values1.size()];
    PLFLT * y1 = new PLFLT[y_values1.size()];
    PLFLT * x1 = new PLFLT[y_values1.size()];
    PLFLT * y_min2 = new PLFLT[y_values2.size()];
    PLFLT * y_max2 = new PLFLT[y_values2.size()];
    PLFLT * y2 = new PLFLT[y_values2.size()];
    PLFLT * x2 = new PLFLT[y_values2.size()];

    float y_high=0;
    float x_max=0;

    for(int i=0; i<y_values1.size(); i++)
    {
        y_min1[i]=y_values1[i]-y_errors1_low[i];
        y_max1[i]=y_values1[i]+y_errors1_high[i];
        y1[i]=y_values1[i];
        x1[i]=x_values1[i];
        if(y_max1[i]>y_high)
        {
            if (y_max1[i]<1.5f*y1[i]) y_high=std::max(y_high,(float)y_max1[i]*1.1f);
            else y_high=std::max(y_high,(float)y1[i]*1.5f);
            if (y1[i]<0.0f)
            {
                if (y_min1[i]>1.5f*y1[i]) y_minimal=std::min(y_minimal,(float)y_min1[i]*1.1f);
                else y_minimal=std::min(y_minimal,(float)y1[i]*1.5f);
            }
        }

        if(x1[i]*1.05>x_max) x_max=x1[i]*1.05;
    }

    for(int i=0; i<y_values2.size(); i++)
    {
        y_min2[i]=y_values2[i]-y_errors2_low[i];
        y_max2[i]=y_values2[i]+y_errors2_high[i];
        y2[i]=y_values2[i];
        x2[i]=x_values2[i];
        if(y_max2[i]>y_high)
        {
            if (y_max2[i]<1.5f*y2[i]) y_high=std::max(y_high,(float)y_max2[i]*1.1f);
            else y_high=std::max(y_high,(float)y2[i]*1.5f);
            if (y2[i]<0.0f)
            {
                if (y_min2[i]>1.5f*y2[i]) y_minimal=std::min(y_minimal,(float)y_min2[i]*1.1f);
                else y_minimal=std::min(y_minimal,(float)y2[i]*1.5f);
            }

        }
        if(x2[i]*1.05>x_max) x_max=x2[i]*1.05;
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind( 0.0, x_max, y_minimal, y_high);
    if (x_max>500) pls->box( "bcinst", 250.0, 5.0, "bcinstv", 0, 0 );
    else if (x_max>50) pls->box( "bcinst", 50.0, 5.0, "bcinstv", 0, 0 );
    else pls->box( "bcinst", 0, 0, "bcinstv", 0, 0 );
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    pls->col1(0.1);
    pls->poin( y_values1.size(), x1 , y1 , 0.0 );
    pls->erry( y_values1.size(), x1, y_min1, y_max1 );

    delete y_min1;
    delete y_max1;
    delete y1;
    delete x1;

    pls->col1(0.8);
    pls->poin( y_values2.size(), x2 , y2 , 0.0 );
    //pls->erry( y_values2.size(), x2, y_min2, y_max2 );

    delete y_min2;
    delete y_max2;
    delete y2;
    delete x2;

    delete pls;
}

void plplot::draw_depth_3d(std::string x_axis, std::string depth_axis, std::string z_axis, std::string title, std::vector<float> x_values,
                           std::vector<float> depth, std::vector< std::vector<float> > values, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    int max_bin_nr=0;
    for(int depth_index=0; depth_index<depth.size(); depth_index++)
    {
        if (values[depth_index].size()>max_bin_nr) max_bin_nr=values[depth_index].size();
    }

    PLFLT * x = new PLFLT[depth.size()*max_bin_nr];
    PLFLT * z = new PLFLT[depth.size()*max_bin_nr];
    PLFLT * depth_values = new PLFLT[depth.size()*max_bin_nr];

    float z_high=0.0f;
    float depth_min=3000.0f;
    float depth_max=0.0f;
    float x_min=0.0f;
    float x_max=0.0f;

    for(int depth_index=0; depth_index<depth.size(); depth_index++)
    {
        for(int i=0; i<values[depth_index].size(); i++)
        {
            int ii=depth_index*max_bin_nr+i;
            depth_values[ii] = depth[depth_index];
            if(depth_values[ii]*1.1>depth_max) depth_max=depth_values[ii]*1.1;
            if(depth_values[ii]*0.9<depth_min) depth_min=depth_values[ii]*0.9;

            z[ii] = values[depth_index][i];
            x[ii] = x_values[i];
            if(z[ii]*1.1>z_high) z_high=z[ii]*1.1;;

            if (x_values[i]*1.05>x_max) x_max=x_values[i]*1.05;
            if (x_values[i]<x_min) x_min=x_values[i];
        }

        for(int i=values[depth_index].size(); i<max_bin_nr; i++)
        {
            int ii=depth_index*max_bin_nr+i;
            depth_values[ii] = depth[depth_index];

            if(max_bin_nr>20) z[ii]=values[depth_index].back();//step
            else z[ii] = 0.0f;
            x[ii] = x_values[i];

            if (x_values[i]*1.05>x_max) x_max=x_values[i]*1.05;
            if (x_values[i]<x_min) x_min=x_values[i];
        }
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);

    int max_bin=max_bin_nr;

    if(x_min<0)//grain sizes and radii
    {
        pls->wind(-3.48, 2.9, -1.3, 1.15);
        pls->w3d(3.5, 4.5, 0.7, -3.0, x_max, depth_min, depth_max, 0.0, z_high, 20.0, -45.0);
        if (depth_max>1000) pls->box3("blnstu", x_axis.c_str(), 0.0, 0, "bnstu", depth_axis.c_str(), 500.0, 5.0, "bnstu", z_axis.c_str(), 0.0, 0);
        else if (depth_max>500) pls->box3("blnstu", x_axis.c_str(), 0.0, 0, "bnstu", depth_axis.c_str(), 250.0, 5.0, "bnstu", z_axis.c_str(), 0.0, 0);
        else pls->box3("blnstu", x_axis.c_str(), 0.0, 0, "bnstu", depth_axis.c_str(), 50.0, 5.0, "bnstu", z_axis.c_str(), 0.0, 0);
    }
    else if(x_max<1)//quantiles
    {
        pls->wind(-2.0, 1.55, -0.75, 1.1);
        pls->w3d(1.0, 6.0, 0.7, x_min, x_max, depth_min, depth_max, 0.0, z_high, 10.0, -20.0);
        if (depth_max>500) pls->box3("bnstu", x_axis.c_str(), 0.05, 5, "bnstu", depth_axis.c_str(), 250.0, 5.0, "bnstu", z_axis.c_str(), 0.0, 0);
        else pls->box3("bnstu", x_axis.c_str(), 0.05, 5, "bnstu", depth_axis.c_str(), 50.0, 5.0, "bnstu", z_axis.c_str(), 0.0, 0);
    }
    else if(x_max>2)//step and radius norm
    {
        if(max_bin>20)
        {
            x_max=x_values[19]*1.05f;
            max_bin=20;
        }

        pls->wind(-2.0, 1.55, -0.75, 1.1);
        pls->w3d(1.0, 6.0, 0.7, x_min, x_max, depth_min, depth_max, 0.0, z_high, 10.0, -20.0);
        if (depth_max>500) pls->box3("bnstu", x_axis.c_str(), 0.0, 0, "bnstu", depth_axis.c_str(), 250.0, 5.0, "bnstu", z_axis.c_str(), 0.0, 0);
        else pls->box3("bnstu", x_axis.c_str(), 0.0, 0, "bnstu", depth_axis.c_str(), 50.0, 5.0, "bnstu", z_axis.c_str(), 0.0, 0);
    }
    else//percent
    {
        pls->wind(-2.0, 1.55, -0.75, 1.1);
        pls->w3d(1.0, 6.0, 0.7, x_min, x_max, depth_min, depth_max, 0.0, z_high, 10.0, -20.0);
        if (depth_max>500) pls->box3("bnstu", x_axis.c_str(), 0.5, 5, "bnstu", depth_axis.c_str(), 250.0, 5.0, "bnstu", z_axis.c_str(), 0.0, 0);
        else pls->box3("bnstu", x_axis.c_str(), 0.5, 5, "bnstu", depth_axis.c_str(), 50.0, 5.0, "bnstu", z_axis.c_str(), 0.0, 0);
    }

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    pls->col1(0.1);

    PLFLT * x_temp = new PLFLT[max_bin];
    PLFLT * z_temp = new PLFLT[max_bin];
    PLFLT * depth_values_temp = new PLFLT[max_bin];

    for(int depth_index=0; depth_index<depth.size(); depth_index++)
    {
        for(int i=0; i<max_bin; i++)
        {
            int ii=depth_index*max_bin_nr+i;
            if(i<values[depth_index].size()) x_temp[i]=x[ii];
            else x_temp[i]=x[depth_index*max_bin_nr+values[depth_index].size()-1];
            z_temp[i]=z[ii];
            depth_values_temp[i]=depth_values[ii];
        }

        pls->line3(max_bin, x_temp, depth_values_temp, z_temp );
    }

    delete x_temp;
    delete z_temp;
    delete depth_values_temp;

    x_temp = new PLFLT[2];
    z_temp = new PLFLT[2];
    depth_values_temp = new PLFLT[2];

    for(int i=0; i<max_bin; i++)
    {
        for(int depth_index=0; depth_index<depth.size()-1; depth_index++)
        {
            int ii=depth_index*max_bin_nr+i;
            x_temp[0]=x[ii];
            z_temp[0]=z[ii];
            depth_values_temp[0]=depth_values[ii];
            
            x_temp[1]=x[ii+max_bin_nr];
            z_temp[1]=z[ii+max_bin_nr];
            depth_values_temp[1]=depth_values[ii+max_bin_nr];

            if(values[depth_index].size()>i && values[depth_index+1].size()>i)
                pls->line3(2, x_temp, depth_values_temp, z_temp );
        }
       
    }
/*
    pls->col1(0.8);

    delete x_temp;
    delete z_temp;
    delete depth_values_temp;

    x_temp = new PLFLT[depth.size()];
    z_temp = new PLFLT[depth.size()];
    depth_values_temp = new PLFLT[depth.size()];

    for(int depth_index=0; depth_index<depth.size(); depth_index++)
    {
        int i=12;
        int ii=depth_index*max_bin_nr+i;
        x_temp[depth_index]=x[ii];
        z_temp[depth_index]=z[ii];
        depth_values_temp[depth_index]=depth_values[ii];
    }

    pls->line3(depth.size(), x_temp, depth_values_temp, z_temp );
*/
    pls->col0(0);
    pls->mtex( "t", 1.0, 0.5, 0.5, title.c_str());

    delete x_temp;
    delete z_temp;
    delete depth_values_temp;

    delete x;
    delete z;
    delete depth_values;

    delete pls;
}

void plplot::draw_depth_color(std::string y_axis, std::string depth_axis, std::string z_axis, std::string title, std::vector<float> y_values,
                              std::vector<float> depth, std::vector< std::vector<float> > values, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    int y_min_bin=y_values.size()-1;
    int y_max_bin=0;
    for(int depth_index=0; depth_index<depth.size(); depth_index++)
    {
        for(int i=0; i<values[depth_index].size(); i++)
        {
            if (i<y_min_bin && values[depth_index][i]>0.0f) y_min_bin=i;
        }
        if (values[depth_index].size()-1>y_max_bin) y_max_bin=values[depth_index].size()-1;
    }

    float y_min=y_values[y_min_bin];
    float y_max=y_values[y_max_bin];
    int YPTS=1+y_max_bin-y_min_bin;

    float depth_max=depth[0];
    int min_depth_bin=3000.0f;
    for(int d1=0; d1<depth.size()-1; d1++)
    {
        int depth_bin=depth[d1+1]-depth[d1];
        if (depth_bin<min_depth_bin) min_depth_bin=std::max(1,depth_bin);
        if (depth[d1+1]>depth_max) depth_max=depth[d1+1];
    }
    int XPTS=1+depth_max/min_depth_bin;

    PLFLT **z = new PLFLT *[XPTS];
    for (int i = 0; i < XPTS; i++ )
    {
        z[i] = new PLFLT[YPTS];
        for(int j = 0; j < YPTS; j++) z[i][j] = 0.0f;
    }

    float z_high=0.0f;
    for(int depth_index=0; depth_index<depth.size(); depth_index++)
    {
        int depth_bin=depth[depth_index]/min_depth_bin;
        for(int i=y_min_bin; i<=y_max_bin && i<values[depth_index].size(); i++)
        {
            z[depth_bin][i-y_min_bin] = values[depth_index][i];
            if (values[depth_index][i]>z_high) z_high=values[depth_index][i];
        }
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind( 0.0, depth_max, y_min, y_max);
    pls->lab(depth_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->spal0( "cmap0_black_on_white.pal" );
    pls->spal1( "cmap1_blue_yellow.pal", true );
    pls->scmap0n( 3 );

    int ns = 50;
    PLFLT *shedge = new PLFLT[ns + 1];
    for (int i = 0; i < ns + 1; i++ )
        shedge[i] = z_high * (PLFLT) i / (PLFLT) (ns-1);

    int fill_width = 2, cont_color = 0, cont_width = 0;

    pls->psty(0);
    pls->shades(z, XPTS, YPTS, NULL, 0, depth_max, y_min, y_max,
                shedge, ns + 1, fill_width, cont_color, cont_width,
                plstream::fill, false, NULL, NULL);

    pls->scmap0(black,black,black,1);
    pls->col0(0);
    pls->box( "bcinst", 250.0, 5.0, "bcilnstu", 0, 0 );//"bcinstv"

    if (depth_max > 2360)
    {
        pls->scmap1l(true,5,pos,pos,pos,pos,NULL);
        pls->col1(1.0);
        plfbox(2214.3, y_max, y_min, 0.0);
    }

    delete pls;
}

void plplot::draw_curv_cross(std::string x_axis, std::string y_axis, std::string title, std::vector<int> x_values, std::vector<float> y_values,
                             std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    PLFLT * y = new PLFLT[y_values.size()];
    PLFLT * x = new PLFLT[y_values.size()];

    float y_high=0.05f;
    int x_min=10;
    int x_max=0;

    for(int i=0; i<y_values.size(); i++)
    {
        x[i]=x_values[i];
        if(x[i]*1.05f>x_max) x_max=x[i]*1.05f;
        if(x[i]<x_min) x_min=x[i];

        y[i]=y_values[i];
        if(y[i]*1.1f>y_high || y[i]*1.1f<-y_high) y_high=fabs(y[i])*1.1f;
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    if (x_min==0) pls->wind( x_min, x_max, -y_high, y_high);//curv
    else if (x_min>0) pls->wind( x_min, x_max, 0, y_high);//deviation
    else pls->wind( x_min, x_max, 0, 255);//cross
    pls->box( "bcinst", 0, 0, "bcinstv", 0, 0 );
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    pls->col1(0.1);
    pls->line(y_values.size(), x, y);

    delete y;
    delete x;

    delete pls;
}

void plplot::draw_correlation(std::string x_axis, std::string y_axis, std::vector< std::vector<int> > grain_areas,
                              std::vector< std::vector<float> > grain_box_flattening, std::vector< std::vector<float> > grain_ellipse_flattening,
                              std::vector< std::vector<float> > x_values, std::vector< std::vector<float> > y_values, correlation corr,
                              std::string dest_path, bool abs1=false, bool abs2=false)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    //plsdev("svg");//sets output device
    plsdev("psc");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    char coeff[20];
    sprintf(coeff, "%.2f", corr.rho);

    char nr[20];
    sprintf(nr, "%i", corr.nr);

    std::string title="Pearson correlation: ";
    title.append(coeff);
    title.append(", nr grains: ");
    title.append(nr);

    PLFLT * y = new PLFLT[corr.nr];
    PLFLT * x = new PLFLT[corr.nr];

    float y_low=1000;
    float y_high=0;
    float x_min=1000;
    float x_max=0;
    int index=0;

    for (int i=0; i<grain_areas.size(); i++)
    {
        for (int area=0; area<grain_areas[i].size(); area++)
        {
            //ellipse not degenerated and exes found correctly
            if (grain_box_flattening[i][grain_areas[i][area]]>0.0f && grain_ellipse_flattening[i][grain_areas[i][area]]>0.0f)
            {
                if (abs1) x_values[i][grain_areas[i][area]]=fabs(x_values[i][grain_areas[i][area]]);
                if (abs2) y_values[i][grain_areas[i][area]]=fabs(y_values[i][grain_areas[i][area]]);

                y[index]=y_values[i][grain_areas[i][area]];
                x[index]=x_values[i][grain_areas[i][area]];
                if(y[index]<y_low) y_low=y[index];
                if(y[index]>y_high) y_high=y[index];
                if(x[index]<x_min) x_min=x[index];
                if(x[index]*1.05f>x_max) x_max=x[index]*1.05f;
                index++;
            }
        }
    }

    x_min-=(x_max-x_min)/10;
    x_max+=(x_max-x_min)/10;
    y_low-=(y_high-y_low)/10;
    y_high+=(y_high-y_low)/10;

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(x_min, x_max, y_low, y_high);
    pls->box("bcinst", 0, 0, "bcinstv", 0, 0);
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    pls->col1(0.1);
    pls->poin(corr.nr, x, y, 1.0);

    if (corr.m != 0.0f || corr.b != 0.0f)
    {
        pls->col1(0.8);

        PLFLT * yy = new PLFLT[2];
        PLFLT * xx = new PLFLT[2];

        if(corr.m>0)
        {
            if (corr.m*x_min+corr.b<y_low) x_min=(y_low-corr.b)/corr.m;
            if (corr.m*x_max+corr.b>y_high) x_max=(y_high-corr.b)/corr.m;
        }
        else
        {
            if (corr.m*x_min+corr.b>y_high) x_min=(y_high-corr.b)/corr.m;
            if (corr.m*x_max+corr.b<y_low) x_max=(y_low-corr.b)/corr.m;
        }

        xx[0]=x_min;
        yy[0]=corr.m*x_min+corr.b;
        xx[1]=x_max,
        yy[1]=corr.m*x_max+corr.b;

        pls->line(2, xx, yy);

        delete yy;
        delete xx;
    }

    delete y;
    delete x;

    delete pls;
}

void plplot::draw_correlation(std::string x_axis, std::string y_axis, std::vector< std::vector<float> > x_values, std::vector< std::vector<float> > y_values,
                              correlation corr, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    //plsdev("svg");//sets output device
    plsdev("psc");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    char coeff[20];
    sprintf(coeff, "%.2f", corr.rho);

    char nr[20];
    sprintf(nr, "%i", corr.nr);

    //std::string title="Pearson correlation: ";
    //title.append(coeff);
    //title.append(", nr boundaries: ");
    //title.append(nr);
    std::string title="Nr boundaries: ";
    title.append(nr);

    PLFLT * y = new PLFLT[corr.nr];
    PLFLT * x = new PLFLT[corr.nr];

    float y_low=1000;
    float y_high=0;
    float x_min=1000;
    float x_max=0;
    int index=0;
    int zeros=0;

    for (int i=0; i<y_values.size(); i++)
    {
        for (int entry=0; entry<y_values[i].size(); entry++)
        {
            if(y_values[i][entry]!=0.0f)
            {
                y[index]=y_values[i][entry];
                x[index]=x_values[i][entry];
                if(y[index]<y_low) y_low=y[index];
                if(y[index]>y_high) y_high=y[index];
                if(x[index]<x_min) x_min=x[index];
                if(x[index]*1.05f>x_max) x_max=x[index]*1.05f;
                index++;
            }
            else zeros++;
        }
    }

    x_min-=(x_max-x_min)/10;
    x_max+=(x_max-x_min)/10;
    y_low-=(y_high-y_low)/10;
    y_high+=(y_high-y_low)/10;

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(x_min, x_max, y_low, y_high);
    pls->box("bcinst", 0, 0, "bcinstv", 0, 0);
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    pls->col1(0.1);
    pls->poin(corr.nr-zeros, x, y, 1.0);
/*
    if (corr.m != 0.0f || corr.b != 0.0f)
    {
        pls->col1(0.8);

        PLFLT * yy = new PLFLT[2];
        PLFLT * xx = new PLFLT[2];

        if(corr.m>0)
        {
            if (corr.m*x_min+corr.b<y_low) x_min=(y_low-corr.b)/corr.m;
            if (corr.m*x_max+corr.b>y_high) x_max=(y_high-corr.b)/corr.m;
        }
        else
        {
            if (corr.m*x_min+corr.b>y_high) x_min=(y_high-corr.b)/corr.m;
            if (corr.m*x_max+corr.b<y_low) x_max=(y_low-corr.b)/corr.m;
        }

        xx[0]=x_min;
        yy[0]=corr.m*x_min+corr.b;
        xx[1]=x_max,
        yy[1]=corr.m*x_max+corr.b;

        pls->line(2, xx, yy);

        delete yy;
        delete xx;
    }
*/
    delete y;
    delete x;

    delete pls;
}

void plplot::draw_correlation(std::string x_axis, std::string y_axis, std::vector<float> x_values, std::vector<float> y_values, std::string dest_path)
{
    correlation corr;

    //Calculate means and sums
    float mean1=0.0f;
    float mean2=0.0f;
    float sumx = 0.0f, sumy = 0.0f, sumx2 = 0.0f, sumxy = 0.0f;
    corr.nr=0;

    for (int i=0; i<x_values.size(); i++)
    {
        if(y_values[i]!=0.0f)
        {
            mean1+=x_values[i];
            mean2+=y_values[i];

            sumx += x_values[i];
            sumy += y_values[i];
            sumx2 += (x_values[i] * x_values[i]);
            sumxy += (x_values[i] * y_values[i]);

            corr.nr++;
        }
    }
    mean1/=corr.nr;
    mean2/=corr.nr;

    //Calculate standard deviations and covariance
    float standard_deviation1=0.0f;
    float standard_deviation2=0.0f;
    float covariance=0.0f;

    for (int i=0; i<y_values.size(); i++)
    {
        if(y_values[i]!=0.0f)
        {
            standard_deviation1+=square(x_values[i]-mean1);
            standard_deviation2+=square(y_values[i]-mean2);

            covariance+=(x_values[i]-mean1)*(y_values[i]-mean2);
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

    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    //plsdev("psc");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    char coeff[20];
    sprintf(coeff, "%.2f", corr.rho);

    char nr[20];
    sprintf(nr, "%i", corr.nr);

    std::string title="Pearson correlation: ";
    title.append(coeff);
    title.append(", nr bags: ");
    title.append(nr);

    PLFLT * y = new PLFLT[corr.nr];
    PLFLT * x = new PLFLT[corr.nr];

    float y_low=1000;
    float y_high=-100;
    float x_min=1000;
    float x_max=0;
    int index=0;

    for (int i=0; i<y_values.size(); i++)
    {
        if(y_values[i]!=0.0f)
        {
            y[index]=y_values[i];
            x[index]=x_values[i];
            if(y[index]<y_low) y_low=y[index];
            if(y[index]>y_high) y_high=y[index];
            if(x[index]<x_min) x_min=x[index];
            if(x[index]*1.05f>x_max) x_max=x[index]*1.05f;
            index++;
        }
    }

    x_min-=(x_max-x_min)/10;
    x_max+=(x_max-x_min)/10;
    y_low-=(y_high-y_low)/10;
    y_high+=(y_high-y_low)/10;

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);
    pls->wind(x_min, x_max, y_low, y_high);
    pls->box("bcinst", 0, 0, "bcinstv", 0, 0);
    pls->lab(x_axis.c_str(), y_axis.c_str(), title.c_str());

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    pls->col1(0.1);
    pls->poin(corr.nr, x, y, 1.0);

    if (corr.m != 0.0f || corr.b != 0.0f)
    {
        pls->col1(0.8);

        PLFLT * yy = new PLFLT[2];
        PLFLT * xx = new PLFLT[2];

        if(corr.m>0)
        {
            if (corr.m*x_min+corr.b<y_low) x_min=(y_low-corr.b)/corr.m;
            if (corr.m*x_max+corr.b>y_high) x_max=(y_high-corr.b)/corr.m;
        }
        else
        {
            if (corr.m*x_min+corr.b>y_high) x_min=(y_high-corr.b)/corr.m;
            if (corr.m*x_max+corr.b<y_low) x_max=(y_low-corr.b)/corr.m;
        }

        xx[0]=x_min;
        yy[0]=corr.m*x_min+corr.b;
        xx[1]=x_max,
        yy[1]=corr.m*x_max+corr.b;

        pls->line(2, xx, yy);

        delete yy;
        delete xx;
    }

    delete y;
    delete x;

    delete pls;
}

void plplot::draw_3d_blocks(std::string x_axis, std::string y_axis, std::string z_axis, std::string title, std::vector<point3d> values, int x_steps, int y_steps,
    std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    float x_max=0;
    float y_max=0;
    float z_min=0;
    float z_max=0;

    for(int i=0; i<values.size(); i++)
    {
        if(values[i].x*1.1f>x_max) x_max=values[i].x*1.1f;
        if(values[i].y*1.1f>y_max) y_max=values[i].y*1.1f;
        if(values[i].z*1.1f>z_max) z_max=values[i].z*1.1f;
        if(values[i].z*1.1f<z_min) z_min=values[i].z*1.1f;
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);

    pls->wind(-1.75, 1.5, -1, 1);
    pls->w3d(1.3, 2.5, 0.5, 0.0, x_max, 0.0, y_max, z_min, z_max, 20.0, -45.0);
    pls->box3("bnstu", x_axis.c_str(), 0.0, 0, "bnstu", y_axis.c_str(), 0.0, 0, "bnstu", z_axis.c_str(), 0.0, 0);

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    pls->col1(0.1);

    PLFLT * x_temp = new PLFLT[y_steps];
    PLFLT * y_temp = new PLFLT[y_steps];
    PLFLT * z_temp = new PLFLT[y_steps];

    //find x values
    std::vector<float> x_values;
    x_values.push_back(values[0].x);

    for(int index=1; index<values.size() && x_values.size()<x_steps; index++)
    {
        bool found=false;

        for(int x=0; x<x_values.size() && !false; x++)
            if (values[index].x==x_values[x]) found=true;

        if (!found) x_values.push_back(values[index].x);
    }

    for(int x=0; x<x_values.size(); x++)
    {
        int i=0;

        for(int index=0; index<values.size(); index++)
        {
            if (values[index].x==x_values[x])
            {
                x_temp[i]=values[index].x;
                y_temp[i]=values[index].y;
                z_temp[i]=values[index].z;
                i++;
            }
        }

        //draw line
        pls->line3(i, x_temp, y_temp, z_temp );
    }

    delete x_temp;
    delete y_temp;
    delete z_temp;

    pls->scmap1l(true,5,pos,red,red,red,NULL);

    pls->col1(0.0);

    x_temp = new PLFLT[2];
    y_temp = new PLFLT[2];
    z_temp = new PLFLT[2];

    for(int x=0; x<x_values.size(); x++)
    {
        //for (float angle=-30.0f; angle<=30.0f; angle+=30.0f)
        float angle=0.0f;
        {
            x_temp[0]=x_values[x];
            y_temp[0]=0;
            z_temp[0]=angle;

            x_temp[1]=x_values[x];
            y_temp[1]=y_max;
            z_temp[1]=angle;

            pls->line3(2, x_temp, y_temp, z_temp );
        }
    }

    pls->col0(0);
    pls->mtex( "t", 1.0, 0.5, 0.5, title.c_str());

    delete x_temp;
    delete y_temp;
    delete z_temp;

    delete pls;
}

void plplot::draw_3d_histogram_shift(std::string x_axis, std::string y_axis, std::string z_axis, std::string title, std::vector<float> y_values,
    std::vector< std::vector<int> > values, float bin_width, float scaling, std::string dest_path)
{
    pls = new plstream();

    plscolbg(255,255,255);//set background white
    plsdev("svg");//sets output device
    plsfnam (dest_path.c_str());//output filename

    // Initialize plplot.
    pls->init();

    int x_mid=values[0].size()/2;
    float width=bin_width/scaling;
    int bin_nr=values[0].size();

    PLFLT * x = new PLFLT[y_values.size()*values[0].size()];
    PLFLT * y = new PLFLT[y_values.size()*values[0].size()];
    PLFLT * z = new PLFLT[y_values.size()*values[0].size()];

    float z_max=0.0f;
    float y_max=0.0f;

    for(int y_index=0; y_index<y_values.size(); y_index++)
    {
        for(int i=0; i<values[0].size(); i++)
        {
            int ii=y_index*values[0].size()+i;
            x[ii] = (float)(i-x_mid)*width;

            y[ii] = y_values[y_index];
            if(y[ii]*1.1>y_max) y_max=y[ii]*1.1;

            z[ii] = values[y_index][i];
            if(z[ii]*1.1>z_max) z_max=z[ii]*1.1;;
        }
    }

    pls->scmap0(black,black,black,1);

    pls->adv(0);
    pls->vsta();
    pls->col0(0);

    pls->wind(-2.35, 1.68, -0.9, 0.95);
    pls->w3d(1.0, 4.0, 0.7, -(x_mid+1)*width, (x_mid+values[0].size()%2+1)*width, 0.0, y_max, 0.0, z_max, 20.0, -50.0);
    pls->box3("bnstu", x_axis.c_str(), 0.0, 0, "bnstu", y_axis.c_str(), 0.0, 0, "bnstu", z_axis.c_str(), 0.0, 0);

    pls->scmap1l(true,5,pos,red,green,blue,NULL);

    pls->col1(0.1);

    PLFLT * x_temp = new PLFLT[values[0].size()];
    PLFLT * y_temp = new PLFLT[values[0].size()];
    PLFLT * z_temp = new PLFLT[values[0].size()];

    for(int y_index=0; y_index<y_values.size(); y_index++)
    {
        for(int i=0; i<values[0].size(); i++)
        {
            int ii=y_index*values[0].size()+i;
            x_temp[i]=x[ii];
            y_temp[i]=y[ii];
            z_temp[i]=z[ii];
        }

        //pls->line3(values[0].size(), x_temp, y_temp, z_temp );
    }

    delete x_temp;
    delete y_temp;
    delete z_temp;

    x_temp = new PLFLT[y_values.size()];
    y_temp = new PLFLT[y_values.size()];
    z_temp = new PLFLT[y_values.size()];

    for(int i=0; i<values[0].size(); i++)
    {
        for(int y_index=0; y_index<y_values.size(); y_index++)
        {
            int ii=y_index*values[0].size()+i;
            x_temp[y_index]=x[ii];
            y_temp[y_index]=y[ii];
            z_temp[y_index]=z[ii];
        }

        pls->line3(y_values.size(), x_temp, y_temp, z_temp );
    }

    pls->col0(0);
    pls->mtex( "t", 1.0, 0.5, 0.5, title.c_str());

    delete x_temp;
    delete y_temp;
    delete z_temp;

    delete x;
    delete y;
    delete z;

    delete pls;
}
