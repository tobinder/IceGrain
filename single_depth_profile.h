/*! \file single_depth_profile.h
 * \brief Single depth profile.
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
void single_lin_depth_profile(std::string path_results, int minimal_bubble_distance, int minimal_grain_size, int grain_size_min, std::string str,
    plplot plot, std::string filename, std::string overview_name, std::string x, std::string y, std::string title, float y_minimal, float scaling=1.0f)
{
    std::vector<float> parameter_depths;
    std::vector<float> parameter_values;

    //string is read from temp file to check whether file is empty
    std::string teststring;

    std::stringstream s;

    if (minimal_bubble_distance==0) s << "." << minimal_grain_size;
    else s << "." << minimal_bubble_distance << "." << minimal_grain_size;
    if (grain_size_min>10) s << "_" << grain_size_min;
    s << ".txt";

    std::string filepath_parameter=path_results;
    filepath_parameter.append(filename.c_str());
    filepath_parameter.append(s.str());

    std::ifstream parameter_file(filepath_parameter.c_str());
    std::ifstream temp_parameter_file(filepath_parameter.c_str());
    temp_parameter_file>>teststring;

    if(parameter_file && teststring.size()!=0)
    {
        std::vector<depth_parameter> entries;

        while(!parameter_file.eof())
        {
            std::string name;
            float mean=-1, stdabw;
            parameter_file>>name>>mean>>stdabw;
            depth_parameter entry;
            entry.name=name;
            entry.mean=mean;
            if(mean!=-1) entries.push_back(entry);
        }

        for (int k=0; k<entries.size(); k++)
        {
            float depth;

            if(entries[k].name=="NGRIP211.bmp") depth=211.0f*0.55f;
            else if(entries[k].name=="NGRIP301.bmp") depth=301.0f*0.55f;
            else if(entries[k].name=="NGRIP401.bmp") depth=401.0f*0.55f;
            else if(entries[k].name=="NGRIP502.bmp") depth=502.0f*0.55f;
            else if(entries[k].name=="NGRIP601.bmp") depth=601.0f*0.55f;
            else if(entries[k].name=="NGRIP701.bmp") depth=701.0f*0.55f;
            else if(entries[k].name=="NGRIP801.bmp") depth=801.0f*0.55f;
            else if(entries[k].name=="NGRIP902.bmp") depth=902.0f*0.55f;
            else if(entries[k].name=="NGRIP1001.bmp") depth=1001.0f*0.55f;
            else if(entries[k].name=="NGRIP1101.bmp") depth=1101.0f*0.55f;
            else if(entries[k].name=="NGRIP1200.bmp") depth=1200.0f*0.55f;
            else if(entries[k].name=="NGRIP1298.bmp") depth=1298.0f*0.55f;
            else if(entries[k].name=="NGRIP1401.bmp") depth=1401.0f*0.55f;
            else if(entries[k].name=="NGRIP1502.bmp") depth=1502.0f*0.55f;
            else if(entries[k].name=="NGRIP1602.bmp") depth=1602.0f*0.55f;
            else if(entries[k].name=="GRIP202.bmp") depth=202.0f*0.55f;
            else if(entries[k].name=="GRIP400.bmp") depth=400.0f*0.55f;
            else if(entries[k].name=="GRIP611.bmp") depth=611.0f*0.55f;
            else if(entries[k].name=="GRIP808.bmp") depth=808.0f*0.55f;
            else if(entries[k].name=="GRIP996.bmp") depth=996.0f*0.55f;
            else if(entries[k].name=="GRIP1208.bmp") depth=1208.0f*0.55f;
            else if(entries[k].name=="GRIP1395.bmp") depth=1395.0f*0.55f;
            else if(entries[k].name=="GRIP1610.bmp") depth=1610.0f*0.55f;
            else if(entries[k].name=="GRIP1819.bmp") depth=1819.0f*0.55f;
            else
            {
                std::cout<<"Error: Depth not defined for "<<entries[k].name<<std::endl;
                exit(-1);
            }

            parameter_depths.push_back(depth);
            parameter_values.push_back(entries[k].mean/scaling);
        }

        parameter_file.close();
        temp_parameter_file.close();
    }

    std::cout<<overview_name.c_str()<<": "<<parameter_values.size()<<std::endl;

    std::string filepath_parameter_out=path_results;
    filepath_parameter_out.append(filename.c_str());
    filepath_parameter_out.append(str.c_str());
    filepath_parameter_out.resize(filepath_parameter_out.size()-3);
    filepath_parameter_out.append("svg");

    if(parameter_values.size()>0) plot.draw_depth(x.c_str(), y.c_str(), title.c_str(), parameter_depths, parameter_values,
        y_minimal, filepath_parameter_out.c_str());
};

void single_depth_profile(std::string path_results, ParameterFile paramFile)
{
   	Parameter<int> low_grain_size;
   	low_grain_size.assign("", "low_grain_size", 0);
    low_grain_size.load(paramFile,"config");

   	Parameter<int> high_grain_size;
   	high_grain_size.assign("", "high_grain_size", 2000);
    high_grain_size.load(paramFile,"config");

   	Parameter<int> grain_size_step;
   	grain_size_step.assign("", "grain_size_step", 5000);
    grain_size_step.load(paramFile,"config");

   	Parameter<int> m_bubble_distance;
   	m_bubble_distance.assign("", "min_bubble_distance", 0);
    m_bubble_distance.load(paramFile,"config");
    int minimal_bubble_distance=m_bubble_distance;

   	Parameter<int> g_step;
   	g_step.assign("", "grain_step", 100);
    g_step.load(paramFile,"config");
    int grain_step=g_step;

    Parameter<int> g_size_min;
    g_size_min.assign("", "grain_size_min", 4);
    g_size_min.load(paramFile,"config");
    int grain_size_min=g_size_min;

    //initialise plplot class
    plplot plot = plplot();

    //scaling from pixels to length
    std::string length_unit="Length in mm (";
    char length_pixels[20];
    sprintf(length_pixels, "%.0f", length_scaling);
    length_unit.append(length_pixels);
    length_unit.append(" pixels)");

    //scaling from pixels to area
    std::string area_unit="Size in mm^2 (";
    char area_pixels[20];
    sprintf(area_pixels, "%.0f", area_scaling);
    area_unit.append(area_pixels);
    area_unit.append(" pixels)");

    float grain_bin_width=5.0f;
    float grain_equiv_radius_bin_width=10.0f;
    float grain_equiv_radius_norm_bin_width=0.2f;

    for(int minimal_grain_size=low_grain_size; minimal_grain_size<high_grain_size+1; minimal_grain_size+=grain_size_step)
    {
        std::cout<<std::endl;
        std::cout<<"Overview for parameters"<<std::endl;
        std::cout<<"Minimal bubble distance: "<<minimal_bubble_distance<<std::endl;
        std::cout<<"Minimal grain size: "<< minimal_grain_size<<std::endl;
        std::cout<<"Nr of different depths found:"<<std::endl;

        std::stringstream s;
        if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
        else s << "." << minimal_grain_size;
        if (grain_size_min>10) s << "_" << grain_size_min;
        s << ".txt";

        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grainshape", "Grain shape", "Depth in m", "Roundness factor (4pi*Area/Perimeter^2)", "Grain roundness profile", 0.28f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "nr_grainarcs", "Number of grain boundaries", "Depth in m", "Nr of grain boundaries",
            "Number of grain boundaries profile", 0.0f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "nr_neighbors", "Number of grain neighbors", "Depth in m", "Nr of grain neighbors", "Number of grain neighbors profile",
             0.0f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "length_grainarcs", "Grain arc length", "Depth in m", length_unit, "Grain boundary length profile", 0.0f,
            length_scaling);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "longest_grainarcs", "Grain longest arc length", "Depth in m", length_unit, "Grain longest boundary length profile",
            0.0f, length_scaling);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "dislocation_densities", "Dislocation density differences", "Depth in m", "Dislocation density difference in m^(-2)",
            "Profile of mean dislocation density difference at grain boundaries", 0.0f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_equivradius", "Equivalent grain radius", "Depth in m", length_unit, "Grain equivalent radius profile", 0.0f,
            length_scaling);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_boxshape", "Grain box shape", "Depth in m", "Box flattening factor", "Box flattening profile", 0.9f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_boxwidth", "Grain box width", "Depth in m", length_unit, "Grain box width profile", 0.0f, length_scaling);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_boxheight", "Grain box height", "Depth in m", length_unit, "Grain box height profile", 0.0f, length_scaling);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_ellipselong", "Grain ellipse long axis", "Depth in m", length_unit, "Ellipse long axis profile", 0.0f,
            length_scaling);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_ellipseshape", "Grain ellipse shape", "Depth in m", "Ellipse flattening factor",
            "Ellipse flattening profile", 0.9f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_ellipseangle", "Grain ellipse orientation", "Depth in m", "Ellipse long axis angle",
            "Ellipse long axis angle profile", 0.0f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grainwidth", "Grain width", "Depth in m", length_unit, "Grain width profile", 0.0f, length_scaling);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grainheight", "Grain height", "Depth in m", length_unit, "Grain height profile", 0.0f, length_scaling);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grainflattening", "Vertical grain flattening", "Depth in m", "Vertical grain flattening factor",
            "Vertical grain flattening profile", 0.9f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "turning_points", "Number of turning points", "Depth in m", "Nr of turning points",
            "Grain boundary turning points profile", 0.0f);

        std::vector<float> grainsize_all_depths;
        std::vector<float> grainsize_all_values;

        std::vector<float> grainsize_fit_depths;
        std::vector<float> grainsize_fit_values;

        std::vector<float> corr_ellipse_box_flattening_depths;
        std::vector<float> corr_ellipse_box_flattening_values;

        std::vector< std::vector<float> > grainsize_step_depths;
        std::vector< std::vector<float> > grainsize_step_values;

        std::vector< std::vector<float> > grainsize_percent_depths(9);
        std::vector< std::vector<float> > grainsize_percent_values(9);

        std::vector< std::vector<float> > grainsize_quantile_depths(10);
        std::vector< std::vector<float> > grainsize_quantile_values(10);

        std::vector< std::vector<float> > grainsize_bin_depths;
        std::vector< std::vector<float> > grainsize_bin_values;

        std::vector< std::vector<float> > grainradius_bin_depths;
        std::vector< std::vector<float> > grainradius_bin_values;
        
        std::vector< std::vector<float> > grainradius_norm_bin_depths;
        std::vector< std::vector<float> > grainradius_norm_bin_values;

        std::vector<float> grainsize_step_profile_depths;
        std::vector<float> grainsize_step_profile_x_values;
        for(int s=0; s<500; s++) grainsize_step_profile_x_values.push_back(s+1);
        std::vector< std::vector<float> > grainsize_step_profile_values;

        std::vector<float> grainsize_percent_profile_depths;
        std::vector<float> grainsize_percent_profile_x_values;
        for(int p=0; p<10; p++) grainsize_percent_profile_x_values.push_back((float)(p+1)/10.0f);
        std::vector< std::vector<float> > grainsize_percent_profile_values;

        std::vector<float> grainsize_quantile_profile_depths;
        std::vector<float> grainsize_quantile_profile_x_values;
        for(int q=0; q<10; q++) grainsize_quantile_profile_x_values.push_back((float)(q+1)/100.0f);
        std::vector< std::vector<float> > grainsize_quantile_profile_values;

        std::vector<float> grainsize_bin_profile_depths;
        std::vector<float> grainsize_bin_profile_x_values;
        for(int bin=0; bin<50; bin++) grainsize_bin_profile_x_values.push_back(((float)bin/grain_bin_width)-log10(area_scaling));
        std::vector< std::vector<float> > grainsize_bin_profile_values;

        std::vector<float> grainradius_bin_profile_depths;
        std::vector<float> grainradius_bin_profile_x_values;
        for(int bin=0; bin<50; bin++)
            grainradius_bin_profile_x_values.push_back(((float)bin/grain_equiv_radius_bin_width)-log10(length_scaling));
        std::vector< std::vector<float> > grainradius_bin_profile_values;

        std::vector<float> grainradius_norm_bin_profile_depths;
        std::vector<float> grainradius_norm_bin_profile_x_values;
        for(int bin=0; bin<100; bin++) grainradius_norm_bin_profile_x_values.push_back((float)bin*grain_equiv_radius_norm_bin_width);
        std::vector< std::vector<float> > grainradius_norm_bin_profile_values;

        std::vector< std::vector<float> > grainsize_combined_depths;
        grainsize_combined_depths.resize(4);
        std::vector< std::vector<float> > grainsize_combined_values;
        grainsize_combined_values.resize(4);

        std::vector<float> dihedral_angles_depths;
        std::vector<float> dihedral_angles_values;
        std::vector<float> dihedral_angles_errors;

        //string is read from temp file to check whether file is empty
        std::string teststring;

        std::string filepath_grainsize_all=path_results;
        filepath_grainsize_all.append("grainsize_all");
        filepath_grainsize_all.append(s.str());

        std::string filepath_grainsize_fit=path_results;
        filepath_grainsize_fit.append("grainsize_fit");
        filepath_grainsize_fit.append(s.str());

        std::string filepath_corr_ellipse_box_flattening=path_results;
        filepath_corr_ellipse_box_flattening.append("corr_ellipse_box_flattening");
        filepath_corr_ellipse_box_flattening.append(s.str());

        std::string filepath_grainsize_step=path_results;
        filepath_grainsize_step.append("grainsize_step");
        filepath_grainsize_step.append(s.str());

        std::string filepath_grainsize_percent=path_results;
        filepath_grainsize_percent.append("grainsize_percent");
        filepath_grainsize_percent.append(s.str());

        std::string filepath_grainsize_quantile=path_results;
        filepath_grainsize_quantile.append("grainsize_quantile");
        filepath_grainsize_quantile.append(s.str());

        std::string filepath_grainsize_bin=path_results;
        filepath_grainsize_bin.append("grainsize_bin");
        filepath_grainsize_bin.append(s.str());

        std::string filepath_grainradius_bin=path_results;
        filepath_grainradius_bin.append("grainradius_bin");
        filepath_grainradius_bin.append(s.str());

        std::string filepath_grainradius_norm_bin=path_results;
        filepath_grainradius_norm_bin.append("grainradius_norm_bin");
        filepath_grainradius_norm_bin.append(s.str());

        //grain size all
        std::ifstream grainsize_all_file(filepath_grainsize_all.c_str());
        std::ifstream temp_grainsize_all_file(filepath_grainsize_all.c_str());
        temp_grainsize_all_file>>teststring;

        if(grainsize_all_file && teststring.size()!=0)
        {
            std::vector<depth_parameter> entries;

	        while(!grainsize_all_file.eof())
	        {
                std::string name;
                float mean=-1, stdabw;
		        grainsize_all_file>>name>>mean>>stdabw;
                depth_parameter entry;
                entry.name=name;
                entry.mean=mean;
                if(mean!=-1) entries.push_back(entry);
	        }

            for (int k=0; k<entries.size(); k++)
            {
                float depth;

                if(entries[k].name=="NGRIP211.bmp") depth=211.0f*0.55f;
                else if(entries[k].name=="NGRIP301.bmp") depth=301.0f*0.55f;
                else if(entries[k].name=="NGRIP401.bmp") depth=401.0f*0.55f;
                else if(entries[k].name=="NGRIP502.bmp") depth=502.0f*0.55f;
                else if(entries[k].name=="NGRIP601.bmp") depth=601.0f*0.55f;
                else if(entries[k].name=="NGRIP701.bmp") depth=701.0f*0.55f;
                else if(entries[k].name=="NGRIP801.bmp") depth=801.0f*0.55f;
                else if(entries[k].name=="NGRIP902.bmp") depth=902.0f*0.55f;
                else if(entries[k].name=="NGRIP1001.bmp") depth=1001.0f*0.55f;
                else if(entries[k].name=="NGRIP1101.bmp") depth=1101.0f*0.55f;
                else if(entries[k].name=="NGRIP1200.bmp") depth=1200.0f*0.55f;
                else if(entries[k].name=="NGRIP1298.bmp") depth=1298.0f*0.55f;
                else if(entries[k].name=="NGRIP1401.bmp") depth=1401.0f*0.55f;
                else if(entries[k].name=="NGRIP1502.bmp") depth=1502.0f*0.55f;
                else if(entries[k].name=="NGRIP1602.bmp") depth=1602.0f*0.55f;
                else if(entries[k].name=="GRIP202.bmp") depth=202.0f*0.55f;
                else if(entries[k].name=="GRIP400.bmp") depth=400.0f*0.55f;
                else if(entries[k].name=="GRIP611.bmp") depth=611.0f*0.55f;
                else if(entries[k].name=="GRIP808.bmp") depth=808.0f*0.55f;
                else if(entries[k].name=="GRIP996.bmp") depth=996.0f*0.55f;
                else if(entries[k].name=="GRIP1208.bmp") depth=1208.0f*0.55f;
                else if(entries[k].name=="GRIP1395.bmp") depth=1395.0f*0.55f;
                else if(entries[k].name=="GRIP1610.bmp") depth=1610.0f*0.55f;
                else if(entries[k].name=="GRIP1819.bmp") depth=1819.0f*0.55f;

                grainsize_all_depths.push_back(depth);
                grainsize_all_values.push_back(entries[k].mean/area_scaling);

                grainsize_combined_depths[0].push_back(depth);
                grainsize_combined_values[0].push_back(entries[k].mean/area_scaling);
            }

            grainsize_all_file.close();
            temp_grainsize_all_file.close();
        }

        //grainsize fit
        std::ifstream grainsize_fit_file(filepath_grainsize_fit.c_str());
        std::ifstream temp_grainsize_fit_file(filepath_grainsize_fit.c_str());
        temp_grainsize_fit_file>>teststring;

        if(grainsize_fit_file && teststring.size()!=0)
        {
            std::vector<depth_parameter> entries;

            while(!grainsize_fit_file.eof())
            {
                std::string name;
                float mean=-1, stdabw_low, stdabw_high;
                grainsize_fit_file>>name>>mean>>stdabw_low>>stdabw_high;
                depth_parameter entry;
                entry.name=name;
                entry.mean=mean;
                if(mean!=-1) entries.push_back(entry);
            }

            for (int k=0; k<entries.size(); k++)
            {
                float depth;

                if(entries[k].name=="NGRIP211.bmp") depth=211.0f*0.55f;
                else if(entries[k].name=="NGRIP301.bmp") depth=301.0f*0.55f;
                else if(entries[k].name=="NGRIP401.bmp") depth=401.0f*0.55f;
                else if(entries[k].name=="NGRIP502.bmp") depth=502.0f*0.55f;
                else if(entries[k].name=="NGRIP601.bmp") depth=601.0f*0.55f;
                else if(entries[k].name=="NGRIP701.bmp") depth=701.0f*0.55f;
                else if(entries[k].name=="NGRIP801.bmp") depth=801.0f*0.55f;
                else if(entries[k].name=="NGRIP902.bmp") depth=902.0f*0.55f;
                else if(entries[k].name=="NGRIP1001.bmp") depth=1001.0f*0.55f;
                else if(entries[k].name=="NGRIP1101.bmp") depth=1101.0f*0.55f;
                else if(entries[k].name=="NGRIP1200.bmp") depth=1200.0f*0.55f;
                else if(entries[k].name=="NGRIP1298.bmp") depth=1298.0f*0.55f;
                else if(entries[k].name=="NGRIP1401.bmp") depth=1401.0f*0.55f;
                else if(entries[k].name=="NGRIP1502.bmp") depth=1502.0f*0.55f;
                else if(entries[k].name=="NGRIP1602.bmp") depth=1602.0f*0.55f;
                else if(entries[k].name=="GRIP202.bmp") depth=202.0f*0.55f;
                else if(entries[k].name=="GRIP400.bmp") depth=400.0f*0.55f;
                else if(entries[k].name=="GRIP611.bmp") depth=611.0f*0.55f;
                else if(entries[k].name=="GRIP808.bmp") depth=808.0f*0.55f;
                else if(entries[k].name=="GRIP996.bmp") depth=996.0f*0.55f;
                else if(entries[k].name=="GRIP1208.bmp") depth=1208.0f*0.55f;
                else if(entries[k].name=="GRIP1395.bmp") depth=1395.0f*0.55f;
                else if(entries[k].name=="GRIP1610.bmp") depth=1610.0f*0.55f;
                else if(entries[k].name=="GRIP1819.bmp") depth=1819.0f*0.55f;

                grainsize_fit_depths.push_back(depth);
                grainsize_fit_values.push_back(entries[k].mean/area_scaling);
            }

            grainsize_fit_file.close();
            temp_grainsize_fit_file.close();
        }

        //corr ellipse box flattening
        std::ifstream corr_ellipse_box_flattening_file(filepath_corr_ellipse_box_flattening.c_str());
        std::ifstream temp_corr_ellipse_box_flattening_file(filepath_corr_ellipse_box_flattening.c_str());
        temp_corr_ellipse_box_flattening_file>>teststring;

        if(corr_ellipse_box_flattening_file && teststring.size()!=0)
        {
            std::vector<depth_parameter> entries;

            while(!corr_ellipse_box_flattening_file.eof())
            {
                std::string name;
                float mean=-1;
                corr_ellipse_box_flattening_file>>name>>mean;
                depth_parameter entry;
                entry.name=name;
                entry.mean=mean;
                if(mean!=-1) entries.push_back(entry);
            }

            for (int k=0; k<entries.size(); k++)
            {
                float depth;

                if(entries[k].name=="NGRIP211.bmp") depth=211.0f*0.55f;
                else if(entries[k].name=="NGRIP301.bmp") depth=301.0f*0.55f;
                else if(entries[k].name=="NGRIP401.bmp") depth=401.0f*0.55f;
                else if(entries[k].name=="NGRIP502.bmp") depth=502.0f*0.55f;
                else if(entries[k].name=="NGRIP601.bmp") depth=601.0f*0.55f;
                else if(entries[k].name=="NGRIP701.bmp") depth=701.0f*0.55f;
                else if(entries[k].name=="NGRIP801.bmp") depth=801.0f*0.55f;
                else if(entries[k].name=="NGRIP902.bmp") depth=902.0f*0.55f;
                else if(entries[k].name=="NGRIP1001.bmp") depth=1001.0f*0.55f;
                else if(entries[k].name=="NGRIP1101.bmp") depth=1101.0f*0.55f;
                else if(entries[k].name=="NGRIP1200.bmp") depth=1200.0f*0.55f;
                else if(entries[k].name=="NGRIP1298.bmp") depth=1298.0f*0.55f;
                else if(entries[k].name=="NGRIP1401.bmp") depth=1401.0f*0.55f;
                else if(entries[k].name=="NGRIP1502.bmp") depth=1502.0f*0.55f;
                else if(entries[k].name=="NGRIP1602.bmp") depth=1602.0f*0.55f;
                else if(entries[k].name=="GRIP202.bmp") depth=202.0f*0.55f;
                else if(entries[k].name=="GRIP400.bmp") depth=400.0f*0.55f;
                else if(entries[k].name=="GRIP611.bmp") depth=611.0f*0.55f;
                else if(entries[k].name=="GRIP808.bmp") depth=808.0f*0.55f;
                else if(entries[k].name=="GRIP996.bmp") depth=996.0f*0.55f;
                else if(entries[k].name=="GRIP1208.bmp") depth=1208.0f*0.55f;
                else if(entries[k].name=="GRIP1395.bmp") depth=1395.0f*0.55f;
                else if(entries[k].name=="GRIP1610.bmp") depth=1610.0f*0.55f;
                else if(entries[k].name=="GRIP1819.bmp") depth=1819.0f*0.55f;

                corr_ellipse_box_flattening_depths.push_back(depth);
                corr_ellipse_box_flattening_values.push_back(entries[k].mean);
            }

            corr_ellipse_box_flattening_file.close();
            temp_corr_ellipse_box_flattening_file.close();
        }

        //grain size step
        std::ifstream grainsize_step_file(filepath_grainsize_step.c_str());
        std::ifstream temp_grainsize_step_file(filepath_grainsize_step.c_str());
        temp_grainsize_step_file>>teststring;

        if(grainsize_step_file && teststring.size()!=0)
        {
            std::vector<depth_parameter3> entries;

	        while(!grainsize_step_file.eof())
	        {
                std::string name;
		        grainsize_step_file>>name;                    
                depth_parameter3 entry;
                entry.name=name;

                char line[2048];
                grainsize_step_file.getline(line,2048);

                int old_pos=1;
                for(int c=1; c<grainsize_step_file.gcount(); c++)
                {
                    if(line[c]==32)
                    {
                        char number[20]="";
                        memmove(number+0,line+old_pos,c-old_pos);
                        old_pos=c+1;
                        if (entry.mean.size()==entry.stdabw.size()) entry.mean.push_back(atof(number));
                        else entry.stdabw.push_back(atof(number));
                    }
                }

                float depth=-1.0f;

                if(entry.name=="NGRIP211.bmp") depth=211.0f*0.55f;
                else if(entry.name=="NGRIP301.bmp") depth=301.0f*0.55f;
                else if(entry.name=="NGRIP401.bmp") depth=401.0f*0.55f;
                else if(entry.name=="NGRIP502.bmp") depth=502.0f*0.55f;
                else if(entry.name=="NGRIP601.bmp") depth=601.0f*0.55f;
                else if(entry.name=="NGRIP701.bmp") depth=701.0f*0.55f;
                else if(entry.name=="NGRIP801.bmp") depth=801.0f*0.55f;
                else if(entry.name=="NGRIP902.bmp") depth=902.0f*0.55f;
                else if(entry.name=="NGRIP1001.bmp") depth=1001.0f*0.55f;
                else if(entry.name=="NGRIP1101.bmp") depth=1101.0f*0.55f;
                else if(entry.name=="NGRIP1200.bmp") depth=1200.0f*0.55f;
                else if(entry.name=="NGRIP1298.bmp") depth=1298.0f*0.55f;
                else if(entry.name=="NGRIP1401.bmp") depth=1401.0f*0.55f;
                else if(entry.name=="NGRIP1502.bmp") depth=1502.0f*0.55f;
                else if(entry.name=="NGRIP1602.bmp") depth=1602.0f*0.55f;
                else if(entry.name=="GRIP202.bmp") depth=202.0f*0.55f;
                else if(entry.name=="GRIP400.bmp") depth=400.0f*0.55f;
                else if(entry.name=="GRIP611.bmp") depth=611.0f*0.55f;
                else if(entry.name=="GRIP808.bmp") depth=808.0f*0.55f;
                else if(entry.name=="GRIP996.bmp") depth=996.0f*0.55f;
                else if(entry.name=="GRIP1208.bmp") depth=1208.0f*0.55f;
                else if(entry.name=="GRIP1395.bmp") depth=1395.0f*0.55f;
                else if(entry.name=="GRIP1610.bmp") depth=1610.0f*0.55f;
                else if(entry.name=="GRIP1819.bmp") depth=1819.0f*0.55f;

                if (depth!=-1.0f)
                {
                    grainsize_step_profile_depths.push_back(depth);

                    std::vector<float> new_entry;
                    grainsize_step_profile_values.push_back(new_entry);

                    for(int s=0; s<entry.mean.size(); s++)
                    {
                        if(s+1>grainsize_step_depths.size()) grainsize_step_depths.resize(s+1);
                        if(s+1>grainsize_step_values.size()) grainsize_step_values.resize(s+1);
                        grainsize_step_depths[s].push_back(depth);
                        grainsize_step_values[s].push_back(entry.mean[s]/area_scaling);

                        grainsize_step_profile_values.back().push_back(entry.mean[s]/area_scaling);
                    }  

                    grainsize_combined_depths[1].push_back(depth);
                    grainsize_combined_values[1].push_back(entry.mean[1]/area_scaling);
                }  
	        }

            grainsize_step_file.close();
            temp_grainsize_step_file.close();
        }

        //grain size percent
        std::ifstream grainsize_percent_file(filepath_grainsize_percent.c_str());
        std::ifstream temp_grainsize_percent_file(filepath_grainsize_percent.c_str());
        temp_grainsize_percent_file>>teststring;

        if(grainsize_percent_file && teststring.size()!=0)
        {
            std::vector<depth_parameter3> entries;

	        while(!grainsize_percent_file.eof())
	        {
                std::string name;
		        grainsize_percent_file>>name;
                depth_parameter3 entry;
                entry.name=name;

                float mean[10];
                float stdabw[10];
                bool end=false;

                for (int p=0; p<10 && !end; p++)
                {
                    mean[p]=-1;
    				    grainsize_percent_file>>mean[p]>>stdabw[p];
                    if(mean[p]!=-1)
                    {
                        entry.mean.push_back(mean[p]);
                    }
                    else end=true;
                }

                if(!end) entries.push_back(entry);
	        }

            for (int k=0; k<entries.size(); k++)
            {
                float depth;

                if(entries[k].name=="NGRIP211.bmp") depth=211.0f*0.55f;
                else if(entries[k].name=="NGRIP301.bmp") depth=301.0f*0.55f;
                else if(entries[k].name=="NGRIP401.bmp") depth=401.0f*0.55f;
                else if(entries[k].name=="NGRIP502.bmp") depth=502.0f*0.55f;
                else if(entries[k].name=="NGRIP601.bmp") depth=601.0f*0.55f;
                else if(entries[k].name=="NGRIP701.bmp") depth=701.0f*0.55f;
                else if(entries[k].name=="NGRIP801.bmp") depth=801.0f*0.55f;
                else if(entries[k].name=="NGRIP902.bmp") depth=902.0f*0.55f;
                else if(entries[k].name=="NGRIP1001.bmp") depth=1001.0f*0.55f;
                else if(entries[k].name=="NGRIP1101.bmp") depth=1101.0f*0.55f;
                else if(entries[k].name=="NGRIP1200.bmp") depth=1200.0f*0.55f;
                else if(entries[k].name=="NGRIP1298.bmp") depth=1298.0f*0.55f;
                else if(entries[k].name=="NGRIP1401.bmp") depth=1401.0f*0.55f;
                else if(entries[k].name=="NGRIP1502.bmp") depth=1502.0f*0.55f;
                else if(entries[k].name=="NGRIP1602.bmp") depth=1602.0f*0.55f;
                else if(entries[k].name=="GRIP202.bmp") depth=202.0f*0.55f;
                else if(entries[k].name=="GRIP400.bmp") depth=400.0f*0.55f;
                else if(entries[k].name=="GRIP611.bmp") depth=611.0f*0.55f;
                else if(entries[k].name=="GRIP808.bmp") depth=808.0f*0.55f;
                else if(entries[k].name=="GRIP996.bmp") depth=996.0f*0.55f;
                else if(entries[k].name=="GRIP1208.bmp") depth=1208.0f*0.55f;
                else if(entries[k].name=="GRIP1395.bmp") depth=1395.0f*0.55f;
                else if(entries[k].name=="GRIP1610.bmp") depth=1610.0f*0.55f;
                else if(entries[k].name=="GRIP1819.bmp") depth=1819.0f*0.55f;

                grainsize_percent_profile_depths.push_back(depth);

                std::vector<float> new_entry;
                grainsize_percent_profile_values.push_back(new_entry);

                for(int p=0; p<9; p++)
                {
                    grainsize_percent_depths[p].push_back(depth);
                    grainsize_percent_values[p].push_back(entries[k].mean[p]/area_scaling);

                    grainsize_percent_profile_values.back().push_back(entries[k].mean[p]/area_scaling);
                } 

                grainsize_percent_profile_values.back().push_back(entries[k].mean[9]/area_scaling);
                grainsize_combined_depths[2].push_back(depth);
                grainsize_combined_values[2].push_back(entries[k].mean[2]/area_scaling);
            }

            grainsize_percent_file.close();
            temp_grainsize_percent_file.close();
        }

        //grain size quantile
        std::ifstream grainsize_quantile_file(filepath_grainsize_quantile.c_str());
        std::ifstream temp_grainsize_quantile_file(filepath_grainsize_quantile.c_str());
        temp_grainsize_quantile_file>>teststring;

        if(grainsize_quantile_file && teststring.size()!=0)
        {
            std::vector<depth_parameter3> entries;

	        while(!grainsize_quantile_file.eof())
	        {
                std::string name;
		        grainsize_quantile_file>>name;
                depth_parameter3 entry;
                entry.name=name;

                float mean[10];
                float stdabw[10];
                bool end=false;

                for (int q=0; q<10 && !end; q++)
                {
                    mean[q]=-1;
    				    grainsize_quantile_file>>mean[q]>>stdabw[q];
                    if(mean[q]!=-1)
                    {
                        entry.mean.push_back(mean[q]);
                    }
                    else end=true;
                }

                if(!end) entries.push_back(entry);
	        }

            for (int k=0; k<entries.size(); k++)
            {
                float depth;

                if(entries[k].name=="NGRIP211.bmp") depth=211.0f*0.55f;
                else if(entries[k].name=="NGRIP301.bmp") depth=301.0f*0.55f;
                else if(entries[k].name=="NGRIP401.bmp") depth=401.0f*0.55f;
                else if(entries[k].name=="NGRIP502.bmp") depth=502.0f*0.55f;
                else if(entries[k].name=="NGRIP601.bmp") depth=601.0f*0.55f;
                else if(entries[k].name=="NGRIP701.bmp") depth=701.0f*0.55f;
                else if(entries[k].name=="NGRIP801.bmp") depth=801.0f*0.55f;
                else if(entries[k].name=="NGRIP902.bmp") depth=902.0f*0.55f;
                else if(entries[k].name=="NGRIP1001.bmp") depth=1001.0f*0.55f;
                else if(entries[k].name=="NGRIP1101.bmp") depth=1101.0f*0.55f;
                else if(entries[k].name=="NGRIP1200.bmp") depth=1200.0f*0.55f;
                else if(entries[k].name=="NGRIP1298.bmp") depth=1298.0f*0.55f;
                else if(entries[k].name=="NGRIP1401.bmp") depth=1401.0f*0.55f;
                else if(entries[k].name=="NGRIP1502.bmp") depth=1502.0f*0.55f;
                else if(entries[k].name=="NGRIP1602.bmp") depth=1602.0f*0.55f;
                else if(entries[k].name=="GRIP202.bmp") depth=202.0f*0.55f;
                else if(entries[k].name=="GRIP400.bmp") depth=400.0f*0.55f;
                else if(entries[k].name=="GRIP611.bmp") depth=611.0f*0.55f;
                else if(entries[k].name=="GRIP808.bmp") depth=808.0f*0.55f;
                else if(entries[k].name=="GRIP996.bmp") depth=996.0f*0.55f;
                else if(entries[k].name=="GRIP1208.bmp") depth=1208.0f*0.55f;
                else if(entries[k].name=="GRIP1395.bmp") depth=1395.0f*0.55f;
                else if(entries[k].name=="GRIP1610.bmp") depth=1610.0f*0.55f;
                else if(entries[k].name=="GRIP1819.bmp") depth=1819.0f*0.55f;

                grainsize_quantile_profile_depths.push_back(depth);

                std::vector<float> new_entry;
                grainsize_quantile_profile_values.push_back(new_entry);

                for(int q=0; q<10; q++)
                {
                    grainsize_quantile_depths[q].push_back(depth);
                    grainsize_quantile_values[q].push_back(entries[k].mean[q]/area_scaling);

                    grainsize_quantile_profile_values.back().push_back(entries[k].mean[q]/area_scaling);
                } 

                grainsize_combined_depths[3].push_back(depth);
                grainsize_combined_values[3].push_back(entries[k].mean[9]/area_scaling);
            }

            grainsize_quantile_file.close();
            temp_grainsize_quantile_file.close();
        }

        //grain size bin
        std::ifstream grainsize_bin_file(filepath_grainsize_bin.c_str());
        std::ifstream temp_grainsize_bin_file(filepath_grainsize_bin.c_str());
        temp_grainsize_bin_file>>teststring;

        if(grainsize_bin_file && teststring.size()!=0)
        {
            std::vector<depth_parameter3> entries;

	        while(!grainsize_bin_file.eof())
	        {
                std::string name;
		        grainsize_bin_file>>name;                    
                depth_parameter3 entry;
                entry.name=name;

                char line[256];
                grainsize_bin_file.getline(line,256);

                int old_pos=1;
                for(int c=1; c<grainsize_bin_file.gcount(); c++)
                {
                    if(line[c]==32)
                    {
                        char number[10]="";
                        memmove(number+0,line+old_pos,c-old_pos);
                        old_pos=c+1;
                        entry.mean.push_back(atof(number));
                    }
                }

                float depth=-1.0f;

                if(entry.name=="NGRIP211.bmp") depth=211.0f*0.55f;
                else if(entry.name=="NGRIP301.bmp") depth=301.0f*0.55f;
                else if(entry.name=="NGRIP401.bmp") depth=401.0f*0.55f;
                else if(entry.name=="NGRIP502.bmp") depth=502.0f*0.55f;
                else if(entry.name=="NGRIP601.bmp") depth=601.0f*0.55f;
                else if(entry.name=="NGRIP701.bmp") depth=701.0f*0.55f;
                else if(entry.name=="NGRIP801.bmp") depth=801.0f*0.55f;
                else if(entry.name=="NGRIP902.bmp") depth=902.0f*0.55f;
                else if(entry.name=="NGRIP1001.bmp") depth=1001.0f*0.55f;
                else if(entry.name=="NGRIP1101.bmp") depth=1101.0f*0.55f;
                else if(entry.name=="NGRIP1200.bmp") depth=1200.0f*0.55f;
                else if(entry.name=="NGRIP1298.bmp") depth=1298.0f*0.55f;
                else if(entry.name=="NGRIP1401.bmp") depth=1401.0f*0.55f;
                else if(entry.name=="NGRIP1502.bmp") depth=1502.0f*0.55f;
                else if(entry.name=="NGRIP1602.bmp") depth=1602.0f*0.55f;
                else if(entry.name=="GRIP202.bmp") depth=202.0f*0.55f;
                else if(entry.name=="GRIP400.bmp") depth=400.0f*0.55f;
                else if(entry.name=="GRIP611.bmp") depth=611.0f*0.55f;
                else if(entry.name=="GRIP808.bmp") depth=808.0f*0.55f;
                else if(entry.name=="GRIP996.bmp") depth=996.0f*0.55f;
                else if(entry.name=="GRIP1208.bmp") depth=1208.0f*0.55f;
                else if(entry.name=="GRIP1395.bmp") depth=1395.0f*0.55f;
                else if(entry.name=="GRIP1610.bmp") depth=1610.0f*0.55f;
                else if(entry.name=="GRIP1819.bmp") depth=1819.0f*0.55f;

                if (depth!=-1.0f)
                {
                    int nr_grains_sum=0;

                    grainsize_bin_profile_depths.push_back(depth);

                    std::vector<float> new_entry;
                    grainsize_bin_profile_values.push_back(new_entry);

                    for(int bin=0; bin<entry.mean.size(); bin++)
                    {
                        nr_grains_sum+=(entry.mean[bin]);
                    } 

                    for(int bin=0; bin<entry.mean.size(); bin++)
                    {
                        if(bin+1>grainsize_bin_depths.size()) grainsize_bin_depths.resize(bin+1);
                        if(bin+1>grainsize_bin_values.size()) grainsize_bin_values.resize(bin+1);
                        grainsize_bin_depths[bin].push_back(depth);
                        grainsize_bin_values[bin].push_back(entry.mean[bin]/nr_grains_sum);
      
                        grainsize_bin_profile_values.back().push_back(entry.mean[bin]/nr_grains_sum);
                    }
                } 
	        }

            grainsize_bin_file.close();
            temp_grainsize_bin_file.close();
        }

        //grain radius bin
        std::ifstream grainradius_bin_file(filepath_grainradius_bin.c_str());
        std::ifstream temp_grainradius_bin_file(filepath_grainradius_bin.c_str());
        temp_grainradius_bin_file>>teststring;

        if(grainradius_bin_file && teststring.size()!=0)
        {
            std::vector<depth_parameter3> entries;

	        while(!grainradius_bin_file.eof())
	        {
                std::string name;
		        grainradius_bin_file>>name;                    
                depth_parameter3 entry;
                entry.name=name;

                char line[256];
                grainradius_bin_file.getline(line,256);

                int old_pos=1;
                for(int c=1; c<grainradius_bin_file.gcount(); c++)
                {
                    if(line[c]==32)
                    {
                        char number[10]="";
                        memmove(number+0,line+old_pos,c-old_pos);
                        old_pos=c+1;
                        entry.mean.push_back(atof(number));
                    }
                }

                float depth=-1.0f;

                if(entry.name=="NGRIP211.bmp") depth=211.0f*0.55f;
                else if(entry.name=="NGRIP301.bmp") depth=301.0f*0.55f;
                else if(entry.name=="NGRIP401.bmp") depth=401.0f*0.55f;
                else if(entry.name=="NGRIP502.bmp") depth=502.0f*0.55f;
                else if(entry.name=="NGRIP601.bmp") depth=601.0f*0.55f;
                else if(entry.name=="NGRIP701.bmp") depth=701.0f*0.55f;
                else if(entry.name=="NGRIP801.bmp") depth=801.0f*0.55f;
                else if(entry.name=="NGRIP902.bmp") depth=902.0f*0.55f;
                else if(entry.name=="NGRIP1001.bmp") depth=1001.0f*0.55f;
                else if(entry.name=="NGRIP1101.bmp") depth=1101.0f*0.55f;
                else if(entry.name=="NGRIP1200.bmp") depth=1200.0f*0.55f;
                else if(entry.name=="NGRIP1298.bmp") depth=1298.0f*0.55f;
                else if(entry.name=="NGRIP1401.bmp") depth=1401.0f*0.55f;
                else if(entry.name=="NGRIP1502.bmp") depth=1502.0f*0.55f;
                else if(entry.name=="NGRIP1602.bmp") depth=1602.0f*0.55f;
                else if(entry.name=="GRIP202.bmp") depth=202.0f*0.55f;
                else if(entry.name=="GRIP400.bmp") depth=400.0f*0.55f;
                else if(entry.name=="GRIP611.bmp") depth=611.0f*0.55f;
                else if(entry.name=="GRIP808.bmp") depth=808.0f*0.55f;
                else if(entry.name=="GRIP996.bmp") depth=996.0f*0.55f;
                else if(entry.name=="GRIP1208.bmp") depth=1208.0f*0.55f;
                else if(entry.name=="GRIP1395.bmp") depth=1395.0f*0.55f;
                else if(entry.name=="GRIP1610.bmp") depth=1610.0f*0.55f;
                else if(entry.name=="GRIP1819.bmp") depth=1819.0f*0.55f;

                if (depth!=-1.0f)
                {
                    int nr_grains_sum=0;

                    grainradius_bin_profile_depths.push_back(depth);

                    std::vector<float> new_entry;
                    grainradius_bin_profile_values.push_back(new_entry);

                    for(int bin=0; bin<entry.mean.size(); bin++)
                    {
                        nr_grains_sum+=(entry.mean[bin]);
                    } 

                    for(int bin=0; bin<entry.mean.size(); bin++)
                    {
                        if(bin+1>grainradius_bin_depths.size()) grainradius_bin_depths.resize(bin+1);
                        if(bin+1>grainradius_bin_values.size()) grainradius_bin_values.resize(bin+1);
                        grainradius_bin_depths[bin].push_back(depth);
                        grainradius_bin_values[bin].push_back(entry.mean[bin]/nr_grains_sum);

                        grainradius_bin_profile_values.back().push_back(entry.mean[bin]/nr_grains_sum);
                    }
                }
	        }

            grainradius_bin_file.close();
            temp_grainradius_bin_file.close();
        }

        //grain radius norm bin
        std::ifstream grainradius_norm_bin_file(filepath_grainradius_norm_bin.c_str());
        std::ifstream temp_grainradius_norm_bin_file(filepath_grainradius_norm_bin.c_str());
        temp_grainradius_norm_bin_file>>teststring;

        if(grainradius_norm_bin_file && teststring.size()!=0)
        {
            std::vector<depth_parameter3> entries;

	        while(!grainradius_norm_bin_file.eof())
	        {
                std::string name;
		        grainradius_norm_bin_file>>name;                    
                depth_parameter3 entry;
                entry.name=name;

                char line[1024];
                grainradius_norm_bin_file.getline(line,1024);

                int old_pos=1;
                for(int c=1; c<grainradius_norm_bin_file.gcount(); c++)
                {
                    if(line[c]==32)
                    {
                        char number[10]="";
                        memmove(number+0,line+old_pos,c-old_pos);
                        old_pos=c+1;
                        entry.mean.push_back(atof(number));
                    }
                }

                float depth=-1.0f;

                if(entry.name=="NGRIP211.bmp") depth=211.0f*0.55f;
                else if(entry.name=="NGRIP301.bmp") depth=301.0f*0.55f;
                else if(entry.name=="NGRIP401.bmp") depth=401.0f*0.55f;
                else if(entry.name=="NGRIP502.bmp") depth=502.0f*0.55f;
                else if(entry.name=="NGRIP601.bmp") depth=601.0f*0.55f;
                else if(entry.name=="NGRIP701.bmp") depth=701.0f*0.55f;
                else if(entry.name=="NGRIP801.bmp") depth=801.0f*0.55f;
                else if(entry.name=="NGRIP902.bmp") depth=902.0f*0.55f;
                else if(entry.name=="NGRIP1001.bmp") depth=1001.0f*0.55f;
                else if(entry.name=="NGRIP1101.bmp") depth=1101.0f*0.55f;
                else if(entry.name=="NGRIP1200.bmp") depth=1200.0f*0.55f;
                else if(entry.name=="NGRIP1298.bmp") depth=1298.0f*0.55f;
                else if(entry.name=="NGRIP1401.bmp") depth=1401.0f*0.55f;
                else if(entry.name=="NGRIP1502.bmp") depth=1502.0f*0.55f;
                else if(entry.name=="NGRIP1602.bmp") depth=1602.0f*0.55f;
                else if(entry.name=="GRIP202.bmp") depth=202.0f*0.55f;
                else if(entry.name=="GRIP400.bmp") depth=400.0f*0.55f;
                else if(entry.name=="GRIP611.bmp") depth=611.0f*0.55f;
                else if(entry.name=="GRIP808.bmp") depth=808.0f*0.55f;
                else if(entry.name=="GRIP996.bmp") depth=996.0f*0.55f;
                else if(entry.name=="GRIP1208.bmp") depth=1208.0f*0.55f;
                else if(entry.name=="GRIP1395.bmp") depth=1395.0f*0.55f;
                else if(entry.name=="GRIP1610.bmp") depth=1610.0f*0.55f;
                else if(entry.name=="GRIP1819.bmp") depth=1819.0f*0.55f;

                if (depth!=-1.0f)
                {
                    int nr_grains_sum=0;

                    grainradius_norm_bin_profile_depths.push_back(depth);

                    std::vector<float> new_entry;
                    grainradius_norm_bin_profile_values.push_back(new_entry);
                
                    for(int bin=0; bin<entry.mean.size(); bin++)
                    {
                        nr_grains_sum+=(entry.mean[bin]);
                    } 

                    for(int bin=0; bin<entry.mean.size(); bin++)
                    {
                        if(bin+1>grainradius_norm_bin_depths.size()) grainradius_norm_bin_depths.resize(bin+1);
                        if(bin+1>grainradius_norm_bin_values.size()) grainradius_norm_bin_values.resize(bin+1);
                        grainradius_norm_bin_depths[bin].push_back(depth);
                        grainradius_norm_bin_values[bin].push_back(entry.mean[bin]/nr_grains_sum);

                        grainradius_norm_bin_profile_values.back().push_back(entry.mean[bin]/nr_grains_sum);
                    }
                }
	        }

            grainradius_norm_bin_file.close();
            temp_grainradius_norm_bin_file.close();
        }

        std::string filepath_dihedral_angles=path_results;
        filepath_dihedral_angles.append("dihedral_angles");
        filepath_dihedral_angles.append(s.str());

        //dihedral angles
        std::ifstream dihedral_angles_file(filepath_dihedral_angles.c_str());
        std::ifstream temp_dihedral_angles_file(filepath_dihedral_angles.c_str());
        temp_dihedral_angles_file>>teststring;

        if(dihedral_angles_file && teststring.size()!=0)
        {
            std::vector<depth_parameter> entries;

            while(!dihedral_angles_file.eof())
            {
                std::string name;
                float mean=-1, stdabw;
	            dihedral_angles_file>>name>>mean>>stdabw;
                depth_parameter entry;
                entry.name=name;
                entry.mean=mean;
                entry.stdabw=stdabw;
                if(mean!=-1) entries.push_back(entry);
            }

            for (int k=0; k<entries.size(); k++)
            {
                float depth;

                if(entries[k].name=="NGRIP211.bmp") depth=211.0f*0.55f;
                else if(entries[k].name=="NGRIP301.bmp") depth=301.0f*0.55f;
                else if(entries[k].name=="NGRIP401.bmp") depth=401.0f*0.55f;
                else if(entries[k].name=="NGRIP502.bmp") depth=502.0f*0.55f;
                else if(entries[k].name=="NGRIP601.bmp") depth=601.0f*0.55f;
                else if(entries[k].name=="NGRIP701.bmp") depth=701.0f*0.55f;
                else if(entries[k].name=="NGRIP801.bmp") depth=801.0f*0.55f;
                else if(entries[k].name=="NGRIP902.bmp") depth=902.0f*0.55f;
                else if(entries[k].name=="NGRIP1001.bmp") depth=1001.0f*0.55f;
                else if(entries[k].name=="NGRIP1101.bmp") depth=1101.0f*0.55f;
                else if(entries[k].name=="NGRIP1200.bmp") depth=1200.0f*0.55f;
                else if(entries[k].name=="NGRIP1298.bmp") depth=1298.0f*0.55f;
                else if(entries[k].name=="NGRIP1401.bmp") depth=1401.0f*0.55f;
                else if(entries[k].name=="NGRIP1502.bmp") depth=1502.0f*0.55f;
                else if(entries[k].name=="NGRIP1602.bmp") depth=1602.0f*0.55f;
                else if(entries[k].name=="GRIP202.bmp") depth=202.0f*0.55f;
                else if(entries[k].name=="GRIP400.bmp") depth=400.0f*0.55f;
                else if(entries[k].name=="GRIP611.bmp") depth=611.0f*0.55f;
                else if(entries[k].name=="GRIP808.bmp") depth=808.0f*0.55f;
                else if(entries[k].name=="GRIP996.bmp") depth=996.0f*0.55f;
                else if(entries[k].name=="GRIP1208.bmp") depth=1208.0f*0.55f;
                else if(entries[k].name=="GRIP1395.bmp") depth=1395.0f*0.55f;
                else if(entries[k].name=="GRIP1610.bmp") depth=1610.0f*0.55f;
                else if(entries[k].name=="GRIP1819.bmp") depth=1819.0f*0.55f;

                dihedral_angles_depths.push_back(depth);
                dihedral_angles_values.push_back(entries[k].mean);
                dihedral_angles_errors.push_back(entries[k].stdabw);
            }

            dihedral_angles_file.close();
            temp_dihedral_angles_file.close();
        }

        std::cout<<"Mean grain size (all grains): "<<grainsize_all_values.size()<<std::endl;

        for (int s=0; s<grainsize_step_values.size(); s++)
        {
            std::cout<<"Mean grain size ("<<(s+1)*grain_step<<" largest grains): "<<grainsize_step_values[s].size()<<std::endl;
        }

        for (int p=0; p<grainsize_percent_values.size(); p++)
        {
            std::cout<<"Mean grain size ("<<(p+1)*10<<" percent largest grains): "<<grainsize_percent_values[p].size()<<std::endl;
        }

        for (int q=0; q<grainsize_quantile_values.size(); q++)
        {
            std::cout<<"Mean grain size ("<<99-q<<" percent quantile): "<<grainsize_quantile_values[q].size()<<std::endl;
        }

        for (int bin=0; bin<grainsize_bin_values.size(); bin++)
        {
            std::cout<<"Grain size bin ratio (bin "<<(bin+1)<<"): "<<grainsize_bin_values[bin].size()<<std::endl;
        }

        for (int bin=0; bin<grainradius_bin_values.size(); bin++)
        {
            std::cout<<"Grain radius bin ratio (bin "<<(bin+1)<<"): "<<grainradius_bin_values[bin].size()<<std::endl;
        }

        for (int bin=0; bin<grainradius_norm_bin_values.size(); bin++)
        {
            std::cout<<"Grain radius norm bin ratio (bin "<<(bin+1)<<"): "<<grainradius_norm_bin_values[bin].size()<<std::endl;
        }

        std::cout<<"Mean grain size (step depth profile): "<<grainsize_step_profile_values.size()<<std::endl;
        std::cout<<"Mean grain size (percent depth profile): "<<grainsize_percent_profile_values.size()<<std::endl;
        std::cout<<"Mean grain size (quantile depth profile): "<<grainsize_quantile_profile_values.size()<<std::endl;
        std::cout<<"Mean grain size (bin depth profile): "<<grainsize_bin_profile_values.size()<<std::endl;
        std::cout<<"Mean grain radius (bin depth profile): "<<grainradius_bin_profile_values.size()<<std::endl;
        std::cout<<"Mean grain radius norm (bin depth profile): "<<grainradius_norm_bin_profile_values.size()<<std::endl;

        for (int i=0; i<grainsize_combined_values.size(); i++)
        {
            std::cout<<"Mean grain size (parameter "<<i+1<<"/"<<grainsize_combined_values.size()<<" in combined profile): "<<
                grainsize_combined_values[i].size()<<std::endl;
        }

        std::cout<<"Grain size fit: "<<grainsize_fit_values.size()<<std::endl;
        std::cout<<"Correlation ellipse box flattening: "<<corr_ellipse_box_flattening_values.size()<<std::endl;
        std::cout<<"Dihedral angle: "<<dihedral_angles_values.size()<<std::endl;

        std::string filepath_grainsize_all_out=path_results;
        filepath_grainsize_all_out.append("grainsize_all");
        filepath_grainsize_all_out.append(s.str());
        filepath_grainsize_all_out.resize(filepath_grainsize_all_out.size()-3);
        filepath_grainsize_all_out.append("svg");

        std::string filepath_grainsize_all_combined_out=path_results;
        filepath_grainsize_all_combined_out.append("grainsize_all_combined");
        filepath_grainsize_all_combined_out.append(s.str());
        filepath_grainsize_all_combined_out.resize(filepath_grainsize_all_combined_out.size()-3);
        filepath_grainsize_all_combined_out.append("svg");

        std::string filepath_grainsize_fit_out=path_results;
        filepath_grainsize_fit_out.append("grainsize_fit");
        filepath_grainsize_fit_out.append(s.str());
        filepath_grainsize_fit_out.resize(filepath_grainsize_fit_out.size()-3);
        filepath_grainsize_fit_out.append("svg");

        std::string filepath_grainsize_fit_combined_out=path_results;
        filepath_grainsize_fit_combined_out.append("grainsize_fit_combined");
        filepath_grainsize_fit_combined_out.append(s.str());
        filepath_grainsize_fit_combined_out.resize(filepath_grainsize_fit_combined_out.size()-3);
        filepath_grainsize_fit_combined_out.append("svg");

        std::string filepath_corr_ellipse_box_flattening_out=path_results;
        filepath_corr_ellipse_box_flattening_out.append("corr_ellipse_box_flattening");
        filepath_corr_ellipse_box_flattening_out.append(s.str());
        filepath_corr_ellipse_box_flattening_out.resize(filepath_corr_ellipse_box_flattening_out.size()-3);
        filepath_corr_ellipse_box_flattening_out.append("svg");

        std::string filepath_corr_ellipse_box_flattening_combined_out=path_results;
        filepath_corr_ellipse_box_flattening_combined_out.append("corr_ellipse_box_flattening_combined");
        filepath_corr_ellipse_box_flattening_combined_out.append(s.str());
        filepath_corr_ellipse_box_flattening_combined_out.resize(filepath_corr_ellipse_box_flattening_combined_out.size()-3);
        filepath_corr_ellipse_box_flattening_combined_out.append("svg");

        std::string filepath_grainsize_step_out=path_results;
        filepath_grainsize_step_out.append("step/grainsize_");

        std::string filepath_grainsize_percent_out=path_results;
        filepath_grainsize_percent_out.append("percent/grainsize_");

        std::string filepath_grainsize_quantile_out=path_results;
        filepath_grainsize_quantile_out.append("quantile/grainsize_");

        std::string filepath_grainsize_bin_out=path_results;
        filepath_grainsize_bin_out.append("bin/grainsize_");

        std::string filepath_grainradius_bin_out=path_results;
        filepath_grainradius_bin_out.append("bin/grainradius_");

        std::string filepath_grainradius_norm_bin_out=path_results;
        filepath_grainradius_norm_bin_out.append("bin/grainradius_norm");

        std::string filepath_grainsize_step_profile_out=path_results;
        filepath_grainsize_step_profile_out.append("grainsize_step_profile");
        filepath_grainsize_step_profile_out.append(s.str());
        filepath_grainsize_step_profile_out.resize(filepath_grainsize_step_profile_out.size()-3);
        filepath_grainsize_step_profile_out.append("svg");

        std::string filepath_grainsize_percent_profile_out=path_results;
        filepath_grainsize_percent_profile_out.append("grainsize_percent_profile");
        filepath_grainsize_percent_profile_out.append(s.str());
        filepath_grainsize_percent_profile_out.resize(filepath_grainsize_percent_profile_out.size()-3);
        filepath_grainsize_percent_profile_out.append("svg");

        std::string filepath_grainsize_quantile_profile_out=path_results;
        filepath_grainsize_quantile_profile_out.append("grainsize_quantile_profile");
        filepath_grainsize_quantile_profile_out.append(s.str());
        filepath_grainsize_quantile_profile_out.resize(filepath_grainsize_quantile_profile_out.size()-3);
        filepath_grainsize_quantile_profile_out.append("svg");

        std::string filepath_grainsize_bin_profile_out=path_results;
        filepath_grainsize_bin_profile_out.append("grainsize_bin_profile");
        filepath_grainsize_bin_profile_out.append(s.str());
        filepath_grainsize_bin_profile_out.resize(filepath_grainsize_bin_profile_out.size()-3);
        filepath_grainsize_bin_profile_out.append("svg");

        std::string filepath_grainradius_bin_profile_out=path_results;
        filepath_grainradius_bin_profile_out.append("grainradius_bin_profile");
        filepath_grainradius_bin_profile_out.append(s.str());
        filepath_grainradius_bin_profile_out.resize(filepath_grainradius_bin_profile_out.size()-3);
        filepath_grainradius_bin_profile_out.append("svg");

        std::string filepath_grainradius_norm_bin_profile_out=path_results;
        filepath_grainradius_norm_bin_profile_out.append("grainradius_norm_bin_profile");
        filepath_grainradius_norm_bin_profile_out.append(s.str());
        filepath_grainradius_norm_bin_profile_out.resize(filepath_grainradius_norm_bin_profile_out.size()-3);
        filepath_grainradius_norm_bin_profile_out.append("svg");

        std::string filepath_grainsize_combined_out=path_results;
        filepath_grainsize_combined_out.append("grainsize_combined");
        filepath_grainsize_combined_out.append(s.str());
        filepath_grainsize_combined_out.resize(filepath_grainsize_combined_out.size()-3);
        filepath_grainsize_combined_out.append("svg");

        std::string filepath_dihedral_angles_out=path_results;
        filepath_dihedral_angles_out.append("dihedral_angles");
        filepath_dihedral_angles_out.append(s.str());
        filepath_dihedral_angles_out.resize(filepath_dihedral_angles_out.size()-3);
        filepath_dihedral_angles_out.append("svg");

        if(grainsize_all_values.size()>0) plot.draw_depth("Depth in m", area_unit,
                                                                 "Mean grain size profile (all grains)", grainsize_all_depths, grainsize_all_values,
                                                                 0.0f, filepath_grainsize_all_out.c_str());

        if(grainsize_fit_values.size()>0) plot.draw_depth("Depth in m", area_unit, 
                                                                 "Grain size profile log-normal fit", grainsize_fit_depths, grainsize_fit_values,
                                                                 0.0f,
                                                                 filepath_grainsize_fit_out.c_str());

        if(corr_ellipse_box_flattening_values.size()>0) plot.draw_depth("Depth in m", "Pearson correlation",
                                                                 "Correlation ellipse box flattening profile", corr_ellipse_box_flattening_depths,
                                                                 corr_ellipse_box_flattening_values, 0.0f,
                                                                 filepath_corr_ellipse_box_flattening_out.c_str());

        for (int st=0; st<grainsize_step_values.size(); st++)
        {
            char step[5];
            sprintf(step, "%i", (st+1)*grain_step);

            std::string titel="Mean grain size profile (";
            titel.append(step);
            titel.append(" largest grains)");

            std::string filepath_out=filepath_grainsize_step_out;
            filepath_out.append(step);
            filepath_out.append("step");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(grainsize_step_values[st].size()>0) plot.draw_depth("Depth in m", area_unit, titel, grainsize_step_depths[st],
                                                                 grainsize_step_values[st], 0.0f,
                                                                 filepath_out.c_str());
        }

        for (int p=0; p<grainsize_percent_values.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Largest grains profile (";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainsize_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(grainsize_percent_values[p].size()>0) plot.draw_depth("Depth in m", area_unit, titel,
                                                                 grainsize_percent_depths[p], grainsize_percent_values[p], 0.0f, filepath_out.c_str());
        }

        for (int q=0; q<grainsize_quantile_values.size(); q++)
        {
            char quantile[3];
            sprintf(quantile, "%i", 99-q);

            std::string titel="Grain size quantiles (";
            titel.append(quantile);
            titel.append(" percent)");

            std::string filepath_out=filepath_grainsize_quantile_out;
            filepath_out.append(quantile);
            filepath_out.append("quantile");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(grainsize_quantile_values[q].size()>0) plot.draw_depth("Depth in m", area_unit, titel,
                                                                 grainsize_quantile_depths[q], grainsize_quantile_values[q], 0.0f, filepath_out.c_str());
        }

        for (int b=0; b<grainsize_bin_values.size(); b++)
        {
            bool value_found=false;

            for (int depth_index=0; depth_index<grainsize_bin_values[b].size(); depth_index++)
            {
                if(grainsize_bin_values[b][depth_index]>0.0f) value_found=true;
            }

            if (!value_found) continue;

            char bin[3];
            char size_low[10];
            char size_high[10];
            sprintf(bin, "%i", b+1);
            sprintf(size_low, "%2.1e", pow(10.0f,(float)b/grain_bin_width)/area_scaling);
            sprintf(size_high, "%2.1e", pow(10.0f,(float)(b+1)/grain_bin_width)/area_scaling);

            std::string titel="Grain size bin profile (";
            titel.append(size_low);
            titel.append("mm^2 - ");
            titel.append(size_high);
            titel.append("mm^2)");

            std::string filepath_out=filepath_grainsize_bin_out;
            filepath_out.append(bin);
            filepath_out.append("bin");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(grainsize_bin_values[b].size()>0) plot.draw_depth("Depth in m", "Relative occurrence", titel,
                                                                 grainsize_bin_depths[b], grainsize_bin_values[b], 0.0f, filepath_out.c_str());
        }

        for (int b=0; b<grainradius_bin_values.size(); b++)
        {
            bool value_found=false;

            for (int depth_index=0; depth_index<grainradius_bin_values[b].size(); depth_index++)
            {
                if(grainradius_bin_values[b][depth_index]>0.0f) value_found=true;
            }

            if (!value_found) continue;

            char bin[3];
            char radius_low[10];
            char radius_high[10];
            sprintf(bin, "%i", b+1);
            sprintf(radius_low, "%2.1e", pow(10.0f,(float)b/grain_equiv_radius_bin_width)/length_scaling);
            sprintf(radius_high, "%2.1e", pow(10.0f,(float)(b+1)/grain_equiv_radius_bin_width)/length_scaling);

            std::string titel="Grain radius bin profile (";
            titel.append(radius_low);
            titel.append("mm - ");
            titel.append(radius_high);
            titel.append("mm)");

            std::string filepath_out=filepath_grainradius_bin_out;
            filepath_out.append(bin);
            filepath_out.append("bin");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(grainradius_bin_values[b].size()>0) plot.draw_depth("Depth in m", "Relative occurrence", titel,
                                                                 grainradius_bin_depths[b], grainradius_bin_values[b], 0.0f, filepath_out.c_str());
        }

        for (int b=0; b<grainradius_norm_bin_values.size(); b++)
        {
            bool value_found=false;

            for (int depth_index=0; depth_index<grainradius_norm_bin_values[b].size(); depth_index++)
            {
                if(grainradius_norm_bin_values[b][depth_index]>0.0f) value_found=true;
            }

            if (!value_found) continue;

            char bin[3];
            char radius_low[10];
            char radius_high[10];
            sprintf(bin, "%i", b+1);
            sprintf(radius_low, "%.2f", (float)b*grain_equiv_radius_norm_bin_width);
            sprintf(radius_high, "%.2f", (float)(b+1)*grain_equiv_radius_norm_bin_width);

            std::string titel="Normalized grain radius bin profile (";
            titel.append(radius_low);
            titel.append(" - ");
            titel.append(radius_high);
            titel.append(")");

            std::string filepath_out=filepath_grainradius_norm_bin_out;
            filepath_out.append(bin);
            filepath_out.append("bin");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(grainradius_norm_bin_values[b].size()>0) plot.draw_depth("Depth in m", "Relative occurrence", titel,
                                                                 grainradius_norm_bin_depths[b], grainradius_norm_bin_values[b], 0.0f, filepath_out.c_str());
        }

        char step[20];
        sprintf(step, "%d", grain_step);

        std::string x_axis="Largest grains * ";
        x_axis.append(step);

        if(grainsize_step_profile_values.size()>0) plot.draw_depth_3d(x_axis, "Depth in m", area_unit, "Largest grains profile (absolute)",
                                                                 grainsize_step_profile_x_values, grainsize_step_profile_depths,
                                                                 grainsize_step_profile_values, filepath_grainsize_step_profile_out.c_str());
        if(grainsize_percent_profile_values.size()>0) plot.draw_depth_3d("Largest grains", "Depth in m", area_unit,
                                                                 "Largest grains profile (relative)",
                                                                 grainsize_percent_profile_x_values, grainsize_percent_profile_depths,
                                                                 grainsize_percent_profile_values, filepath_grainsize_percent_profile_out.c_str());
        if(grainsize_quantile_profile_values.size()>0) plot.draw_depth_3d("Quantile 1-", "Depth in m", area_unit, "Grain size quantiles",
                                                                 grainsize_quantile_profile_x_values, grainsize_quantile_profile_depths,
                                                                 grainsize_quantile_profile_values, filepath_grainsize_quantile_profile_out.c_str());
        if(grainsize_bin_profile_values.size()>0) plot.draw_depth_3d(area_unit, "Depth in m", "Relative occurrence", "Grain size profile",
                                                                 grainsize_bin_profile_x_values, grainsize_bin_profile_depths,
                                                                 grainsize_bin_profile_values, filepath_grainsize_bin_profile_out.c_str());
        if(grainradius_bin_profile_values.size()>0) plot.draw_depth_3d(length_unit, "Depth in m", "Relative occurrence", "Grain radius profile",
                                                                 grainradius_bin_profile_x_values, grainradius_bin_profile_depths,
                                                                 grainradius_bin_profile_values, filepath_grainradius_bin_profile_out.c_str());
        if(grainradius_norm_bin_profile_values.size()>0) plot.draw_depth_3d("Normalized grain radius", "Depth in m", "Relative occurrence",
                                                                 "Normalized grain radius profile", grainradius_norm_bin_profile_x_values,
                                                                 grainradius_norm_bin_profile_depths, grainradius_norm_bin_profile_values,
                                                                 filepath_grainradius_norm_bin_profile_out.c_str());

        if(grainsize_combined_values.size()>0) plot.draw_depth("Depth in m", area_unit, "Different grain size parameters", grainsize_combined_depths,
                                                                 grainsize_combined_values, filepath_grainsize_combined_out.c_str());

        if(dihedral_angles_values.size()>0) plot.draw_depth("Depth in m", "Dihedral angle standard deviation [degree]", "Dihedral angle standard deviation profile",
                                                                 dihedral_angles_depths, dihedral_angles_errors, 10.0f, filepath_dihedral_angles_out.c_str());
    }
};
