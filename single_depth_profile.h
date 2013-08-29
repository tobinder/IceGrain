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
float get_depth(std::string name)
{
    float depth=0.0f;

    std::string core=name;
    core.erase(4,core.size()-3);

    if(core!="NGRI" && core!="GRIP" && core!="EDML" && core!="KCI-")
    {
        std::cout<<"Error: Depth not defined for "<<name<<std::endl;
    }
    else
    {
        if(core=="NGRI")
        {
            name.erase(0,5);
            depth=atof(name.c_str())*0.55f;
        }
        else if(core=="GRIP")
        {
            name.erase(0,4);
            depth=atof(name.c_str())*0.55f;
        }
        else if(core=="EDML")
        {
            name.erase(0,4);
            name.erase(4,name.size()-3);
            depth=atof(name.c_str());
        }
        else if(core=="KCI-")
        {
            if(name=="KCI-35_01-120719.bmp") depth=29.79f;
            else if(name=="KCI-44_02-120718.bmp") depth=37.30f;
            else if(name=="KCI-46_x1-120725.bmp") depth=39.82f;
            else if(name=="KCI-49_02-120719.bmp") depth=41.74f;
            else if(name=="KCI-52_x1-120725.bmp") depth=45.07f;
            else if(name=="KCI-54_02-120719.bmp") depth=45.91f;
            else if(name=="KCI-56_x1-120726.bmp") depth=48.42f;
            else if(name=="KCI-58_05-120718.bmp") depth=49.35f;
            else if(name=="KCI-61_x1-120726.bmp") depth=52.32f;
            else if(name=="KCI-64_06-120719.bmp") depth=53.87f;
            else if(name=="KCI-64_12-120720.bmp") depth=53.97f;
            else std::cout<<"Error: Depth not defined for "<<name<<std::endl;
        }
    }

    return depth;
}

void single_lin_depth_profile(std::string path_results, int minimal_bubble_distance, int minimal_grain_size, int grain_size_min,
    std::string str, plplot plot, std::string filename, std::string overview_name, std::string x, std::string y, std::string title,
    float y_minimal, float scaling=1.0f)
{
    std::vector<float> parameter_depths;
    std::vector<float> parameter_values;
    std::vector<float> parameter_errors;

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
            entry.stdabw=stdabw;
            if(mean!=-1) entries.push_back(entry);
        }

        //check for multiple entries for same image
        for (int k=0; k<entries.size(); k++)
        {
            bool other_entry=false;

            for(int j=k+1; j<entries.size() && !other_entry; j++)
            {
                if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
            }

            float depth=get_depth(entries[k].name);

            if (!other_entry && depth>0.0f)
            {
                parameter_depths.push_back(depth);
                parameter_values.push_back(entries[k].mean/scaling);
                parameter_errors.push_back(entries[k].stdabw/scaling);
            }
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

    if(parameter_values.size()>0) plot.draw_depth_errors(x.c_str(), y.c_str(), title.c_str(), parameter_depths, parameter_values,
        parameter_errors, parameter_errors, y_minimal, filepath_parameter_out.c_str());
};

void single_new_lin_depth_profile(std::string path_results, int minimal_bubble_distance, int minimal_grain_size, int grain_size_min,
    std::string str, plplot plot, std::string filename, std::string overview_name, std::string x, std::string y, std::string title,
    float y_minimal, int nr_values_selection, std::vector<depth_numbers> nr_values, int depth_max, float depth_bin_width,
    float scaling=1.0f)
{
    std::vector<float> parameter_depths((size_t)(1+depth_max/depth_bin_width));
    std::vector<float> parameter_values((size_t)(1+depth_max/depth_bin_width),0.0f);
    std::vector<float> parameter_sum((size_t)(1+depth_max/depth_bin_width),0.0f);

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

        //check for multiple entries for same image
        for (int k=0; k<entries.size(); k++)
        {
            bool other_entry=false;

            for(int j=k+1; j<entries.size() && !other_entry; j++)
            {
                if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
            }

            int depth_bin=get_depth(entries[k].name)/depth_bin_width;

            //nr of values found
            if(!other_entry && depth_bin>0 && entries[k].name.compare(nr_values[k].name)==0 && nr_values_selection>0)
            {
                parameter_values[depth_bin]+=entries[k].mean*(float)nr_values[k].nr[nr_values_selection-1];
                parameter_sum[depth_bin]+=nr_values[k].nr[nr_values_selection-1];
            }
        }

        parameter_file.close();
        temp_parameter_file.close();
    }

    for (int bin=0; bin<parameter_values.size(); bin++)
    {
        parameter_depths[bin]=(bin+0.5f)*depth_bin_width;
    }

    for (int bin=0; bin<parameter_values.size(); bin++)
    {
        if(parameter_sum[bin]>0.0f) parameter_values[bin]/=(parameter_sum[bin]*scaling);
        else
        {
            parameter_depths.erase (parameter_depths.begin()+bin);
            parameter_values.erase (parameter_values.begin()+bin);
            parameter_sum.erase (parameter_sum.begin()+bin);
            bin--;
        }
    }

    std::string filepath_parameter_out=path_results;
    filepath_parameter_out.append("new_");
    filepath_parameter_out.append(filename.c_str());
    filepath_parameter_out.append(str.c_str());
    filepath_parameter_out.resize(filepath_parameter_out.size()-3);
    filepath_parameter_out.append("svg");

    if(parameter_values.size()>0) plot.draw_depth(x.c_str(), y.c_str(), title.c_str(), parameter_depths, parameter_values, y_minimal,
        filepath_parameter_out.c_str());
};

void single_depth_profile(std::string path_results, ParameterFile paramFile, float depth_bin_width=0.0f)
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

    Parameter<int> b_step;
    b_step.assign("", "boundary_step", 100);
    b_step.load(paramFile,"config");
    int boundary_step=b_step;

    Parameter<int> g_size_min;
    g_size_min.assign("", "grain_size_min", 4);
    g_size_min.load(paramFile,"config");
    int grain_size_min=g_size_min;

    Parameter<float> length_scaling;
    length_scaling.assign("", "length_scaling", 193.5f);
    length_scaling.load(paramFile,"config");

    Parameter<float> area_scaling;
    area_scaling.assign("", "area_scaling", 37444.0f);
    area_scaling.load(paramFile,"config");

    //initialise plplot class
    plplot plot = plplot();

    //scaling from pixels to length
    std::string length_unit="Length in mm (";
    char length_pixels[20];
    sprintf(length_pixels, "%.0f", length_scaling());
    length_unit.append(length_pixels);
    length_unit.append(" pixels)");

    //scaling from pixels to area
    std::string area_unit="Size in mm^2 (";
    char area_pixels[20];
    sprintf(area_pixels, "%.0f", area_scaling());
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
            length_scaling());
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "longest_grainarcs", "Grain longest arc length", "Depth in m", length_unit, "Grain longest boundary length profile",
            0.0f, length_scaling());
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "dislocation_densities", "Dislocation density differences", "Depth in m", "Dislocation density difference in m^(-2)",
            "Profile of mean dislocation density difference at grain boundaries", 0.0f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_equivradius", "Equivalent grain radius", "Depth in m", length_unit, "Grain equivalent radius profile", 0.0f,
            length_scaling());
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_boxshape", "Grain box shape", "Depth in m", "Box flattening factor", "Box flattening profile", 0.9f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_boxwidth", "Grain box width", "Depth in m", length_unit, "Grain box width profile", 0.0f, length_scaling());
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_boxheight", "Grain box height", "Depth in m", length_unit, "Grain box height profile", 0.0f, length_scaling());
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_ellipselong", "Grain ellipse long axis", "Depth in m", length_unit, "Ellipse long axis profile", 0.0f,
            length_scaling());
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_ellipseshape", "Grain ellipse shape", "Depth in m", "Ellipse flattening factor",
            "Ellipse flattening profile", 0.9f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_ellipseangle", "Grain ellipse orientation", "Depth in m", "Ellipse long axis angle",
            "Ellipse long axis angle profile", 0.0f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "boundary_orientation", "Boundary orientation", "Depth in m", "Boundary orientation",
            "Boundary orientation profile", 0.0f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grainwidth", "Grain width", "Depth in m", length_unit, "Grain width profile", 0.0f, length_scaling());
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grainheight", "Grain height", "Depth in m", length_unit, "Grain height profile", 0.0f, length_scaling());
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grainflattening", "Vertical grain flattening", "Depth in m", "Vertical grain flattening factor",
            "Vertical grain flattening profile", 0.9f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "turning_points", "Number of turning points", "Depth in m", "Nr of turning points",
            "Grain boundary turning points profile", 0.0f);
        single_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
            "grain_perimeter_ratio", "Grain perimeter ratio", "Depth in m", "Perimeter ratio", "Perimeter ratio profile", 0.7f);

        std::vector<float> grainsize_all_depths;
        std::vector<float> grainsize_all_values;
        std::vector<float> grainsize_all_errors;

        std::vector<float> grainsize_fit_depths;
        std::vector<float> grainsize_fit_values;
        std::vector<float> grainsize_fit_errors_low;
        std::vector<float> grainsize_fit_errors_high;

        std::vector<float> corr_ellipse_box_flattening_depths;
        std::vector<float> corr_ellipse_box_flattening_values;

        std::vector< std::vector<float> > grainsize_step_depths;
        std::vector< std::vector<float> > grainsize_step_values;
        std::vector< std::vector<float> > grainsize_step_errors;

        std::vector< std::vector<float> > grainsize_percent_depths(9);
        std::vector< std::vector<float> > grainsize_percent_values(9);
        std::vector< std::vector<float> > grainsize_percent_errors(9);

        std::vector< std::vector<float> > grainsize_quantile_depths(10);
        std::vector< std::vector<float> > grainsize_quantile_values(10);
        std::vector< std::vector<float> > grainsize_quantile_errors(10);

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
        for(int bin=0; bin<50; bin++) grainsize_bin_profile_x_values.push_back(((float)bin/grain_bin_width)-log10(area_scaling()));
        std::vector< std::vector<float> > grainsize_bin_profile_values;

        std::vector<float> grainradius_bin_profile_depths;
        std::vector<float> grainradius_bin_profile_x_values;
        for(int bin=0; bin<50; bin++)
            grainradius_bin_profile_x_values.push_back(((float)bin/grain_equiv_radius_bin_width)-log10(length_scaling()));
        std::vector< std::vector<float> > grainradius_bin_profile_values;

        std::vector<float> grainradius_norm_bin_profile_depths;
        std::vector<float> grainradius_norm_bin_profile_x_values;
        for(int bin=0; bin<100; bin++) grainradius_norm_bin_profile_x_values.push_back((float)bin*grain_equiv_radius_norm_bin_width);
        std::vector< std::vector<float> > grainradius_norm_bin_profile_values;

        std::vector<float> grain_neighbors_bin_profile_depths;
        std::vector<float> grain_neighbors_bin_profile_x_values;
        for(int bin=0; bin<100; bin++) grain_neighbors_bin_profile_x_values.push_back((float)bin);
        std::vector< std::vector<float> > grain_neighbors_bin_profile_values;

        std::vector< std::vector<float> > grainsize_combined_depths;
        grainsize_combined_depths.resize(4);
        std::vector< std::vector<float> > grainsize_combined_values;
        grainsize_combined_values.resize(4);

        std::vector<float> dihedral_angles_depths;
        std::vector<float> dihedral_angles_values;
        std::vector<float> dihedral_angles_errors;

        std::vector<float> dihedral_angles2_depths;
        std::vector<float> dihedral_angles2_values;
        std::vector<float> dihedral_angles2_errors;

        std::vector< std::vector<float> > disl_dens_longest_boundaries_depths;
        std::vector< std::vector<float> > disl_dens_longest_boundaries_values;
        std::vector< std::vector<float> > disl_dens_longest_boundaries_errors;

        std::vector< std::vector<float> > disl_dens_percent_boundaries_depths(10);
        std::vector< std::vector<float> > disl_dens_percent_boundaries_values(10);
        std::vector< std::vector<float> > disl_dens_percent_boundaries_errors(10);

        std::vector<float> disl_dens_longest_boundaries_profile_depths;
        std::vector<float> disl_dens_longest_boundaries_profile_x_values;
        for(int s=0; s<500; s++) disl_dens_longest_boundaries_profile_x_values.push_back(s+1);
        std::vector< std::vector<float> > disl_dens_longest_boundaries_profile_values;

        std::vector<float> disl_dens_percent_boundaries_profile_depths;
        std::vector<float> disl_dens_percent_boundaries_profile_x_values;
        for(int p=0; p<10; p++) disl_dens_percent_boundaries_profile_x_values.push_back((float)(p+1)/10.0f);
        std::vector< std::vector<float> > disl_dens_percent_boundaries_profile_values;

        std::vector< std::vector<float> > percentage_filled_grains_depths;
        percentage_filled_grains_depths.resize(4);
        std::vector< std::vector<float> > percentage_filled_grains_values;
        percentage_filled_grains_values.resize(4);

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

        std::string filepath_grain_neighbors_bin=path_results;
        filepath_grain_neighbors_bin.append("grain_neighbors_bin");
        filepath_grain_neighbors_bin.append(s.str());

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
                entry.stdabw=stdabw;
                if(mean!=-1) entries.push_back(entry);
	        }

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if (!other_entry && depth>0.0f)
                {
                    grainsize_all_depths.push_back(depth);
                    grainsize_all_values.push_back(entries[k].mean/area_scaling());
                    grainsize_all_errors.push_back(entries[k].stdabw/area_scaling());

                    grainsize_combined_depths[0].push_back(depth);
                    grainsize_combined_values[0].push_back(entries[k].mean/area_scaling());
                }
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
                entry.stdabw=stdabw_low;
                entry.stdabw2=stdabw_high;
                if(mean!=-1) entries.push_back(entry);
            }

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if (!other_entry && depth>0.0f)
                {
                    grainsize_fit_depths.push_back(depth);
                    grainsize_fit_values.push_back(entries[k].mean/area_scaling());
                    grainsize_fit_errors_low.push_back(entries[k].stdabw/area_scaling());
                    grainsize_fit_errors_high.push_back(std::min(1000.0f,entries[k].stdabw2/area_scaling()));
                }
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

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if (!other_entry && depth>0.0f)
                {
                    corr_ellipse_box_flattening_depths.push_back(depth);
                    corr_ellipse_box_flattening_values.push_back(entries[k].mean);
                }
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

                char line[8192];
                grainsize_step_file.getline(line,8192);

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

                if(entry.mean.size()>0)
                {
                    entries.push_back(entry);
                }
            }

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if (!other_entry && depth>0.0f)
                {
                    grainsize_step_profile_depths.push_back(depth);

                    std::vector<float> new_entry;
                    grainsize_step_profile_values.push_back(new_entry);

                    for(int s=0; s<entries[k].mean.size(); s++)
                    {
                        if(s+1>grainsize_step_depths.size()) grainsize_step_depths.resize(s+1);
                        if(s+1>grainsize_step_values.size()) grainsize_step_values.resize(s+1);
                        if(s+1>grainsize_step_errors.size()) grainsize_step_errors.resize(s+1);
                        grainsize_step_depths[s].push_back(depth);
                        grainsize_step_values[s].push_back(entries[k].mean[s]/area_scaling());
                        grainsize_step_errors[s].push_back(entries[k].stdabw[s]/area_scaling());

                        grainsize_step_profile_values.back().push_back(entries[k].mean[s]/area_scaling());
                    }  

                    grainsize_combined_depths[1].push_back(depth);
                    grainsize_combined_values[1].push_back(entries[k].mean[1]/area_scaling());
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
                        entry.stdabw.push_back(stdabw[p]);
                    }
                    else end=true;
                }

                if(!end) entries.push_back(entry);
	        }

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if (!other_entry && depth>0.0f)
                {
                    grainsize_percent_profile_depths.push_back(depth);

                    std::vector<float> new_entry;
                    grainsize_percent_profile_values.push_back(new_entry);

                    for(int p=0; p<9; p++)
                    {
                        grainsize_percent_depths[p].push_back(depth);
                        grainsize_percent_values[p].push_back(entries[k].mean[p]/area_scaling());
                        grainsize_percent_errors[p].push_back(entries[k].stdabw[p]/area_scaling());

                        grainsize_percent_profile_values.back().push_back(entries[k].mean[p]/area_scaling());
                    } 

                    grainsize_percent_profile_values.back().push_back(entries[k].mean[9]/area_scaling());
                    grainsize_combined_depths[2].push_back(depth);
                    grainsize_combined_values[2].push_back(entries[k].mean[2]/area_scaling());
                }
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
                        entry.stdabw.push_back(stdabw[q]);
                    }
                    else end=true;
                }

                if(!end) entries.push_back(entry);
	        }

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if (!other_entry && depth>0.0f)
                {
                    grainsize_quantile_profile_depths.push_back(depth);

                    std::vector<float> new_entry;
                    grainsize_quantile_profile_values.push_back(new_entry);

                    for(int q=0; q<10; q++)
                    {
                        grainsize_quantile_depths[q].push_back(depth);
                        grainsize_quantile_values[q].push_back(entries[k].mean[q]/area_scaling());
                        grainsize_quantile_errors[q].push_back(entries[k].stdabw[q]/area_scaling());

                        grainsize_quantile_profile_values.back().push_back(entries[k].mean[q]/area_scaling());
                    } 

                    grainsize_combined_depths[3].push_back(depth);
                    grainsize_combined_values[3].push_back(entries[k].mean[9]/area_scaling());
                }
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

                if(entry.mean.size()>0)
                {
                    entries.push_back(entry);
                }
            }

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if (!other_entry && depth>0.0f)
                {
                    int nr_grains_sum=0;

                    grainsize_bin_profile_depths.push_back(depth);

                    std::vector<float> new_entry;
                    grainsize_bin_profile_values.push_back(new_entry);

                    for(int bin=0; bin<entries[k].mean.size(); bin++)
                    {
                        nr_grains_sum+=(entries[k].mean[bin]);
                    } 

                    for(int bin=0; bin<entries[k].mean.size(); bin++)
                    {
                        if(bin+1>grainsize_bin_depths.size()) grainsize_bin_depths.resize(bin+1);
                        if(bin+1>grainsize_bin_values.size()) grainsize_bin_values.resize(bin+1);
                        grainsize_bin_depths[bin].push_back(depth);
                        grainsize_bin_values[bin].push_back(entries[k].mean[bin]/nr_grains_sum);
      
                        grainsize_bin_profile_values.back().push_back(entries[k].mean[bin]/nr_grains_sum);
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

                if(entry.mean.size()>0)
                {
                    entries.push_back(entry);
                }
            }

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if (!other_entry && depth>0.0f)
                {
                    int nr_grains_sum=0;

                    grainradius_bin_profile_depths.push_back(depth);

                    std::vector<float> new_entry;
                    grainradius_bin_profile_values.push_back(new_entry);

                    for(int bin=0; bin<entries[k].mean.size(); bin++)
                    {
                        nr_grains_sum+=(entries[k].mean[bin]);
                    } 

                    for(int bin=0; bin<entries[k].mean.size(); bin++)
                    {
                        if(bin+1>grainradius_bin_depths.size()) grainradius_bin_depths.resize(bin+1);
                        if(bin+1>grainradius_bin_values.size()) grainradius_bin_values.resize(bin+1);
                        grainradius_bin_depths[bin].push_back(depth);
                        grainradius_bin_values[bin].push_back(entries[k].mean[bin]/nr_grains_sum);

                        grainradius_bin_profile_values.back().push_back(entries[k].mean[bin]/nr_grains_sum);
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

                if(entry.mean.size()>0)
                {
                    entries.push_back(entry);
                }
            }

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if (!other_entry && depth>0.0f)
                {
                    int nr_grains_sum=0;

                    grainradius_norm_bin_profile_depths.push_back(depth);

                    std::vector<float> new_entry;
                    grainradius_norm_bin_profile_values.push_back(new_entry);
                
                    for(int bin=0; bin<entries[k].mean.size(); bin++)
                    {
                        nr_grains_sum+=(entries[k].mean[bin]);
                    } 

                    for(int bin=0; bin<entries[k].mean.size(); bin++)
                    {
                        if(bin+1>grainradius_norm_bin_depths.size()) grainradius_norm_bin_depths.resize(bin+1);
                        if(bin+1>grainradius_norm_bin_values.size()) grainradius_norm_bin_values.resize(bin+1);
                        grainradius_norm_bin_depths[bin].push_back(depth);
                        grainradius_norm_bin_values[bin].push_back(entries[k].mean[bin]/nr_grains_sum);

                        grainradius_norm_bin_profile_values.back().push_back(entries[k].mean[bin]/nr_grains_sum);
                    }
                }
	        }

            grainradius_norm_bin_file.close();
            temp_grainradius_norm_bin_file.close();
        }

        //grain neighbors bin
        std::ifstream grain_neighbors_bin_file(filepath_grain_neighbors_bin.c_str());
        std::ifstream temp_grain_neighbors_bin_file(filepath_grain_neighbors_bin.c_str());
        temp_grain_neighbors_bin_file>>teststring;

        if(grain_neighbors_bin_file && teststring.size()!=0)
        {
            std::vector<depth_parameter3> entries;

	        while(!grain_neighbors_bin_file.eof())
	        {
                std::string name;
		        grain_neighbors_bin_file>>name;                    
                depth_parameter3 entry;
                entry.name=name;

                char line[1024];
                grain_neighbors_bin_file.getline(line,1024);

                int old_pos=1;
                for(int c=1; c<grain_neighbors_bin_file.gcount(); c++)
                {
                    if(line[c]==32)
                    {
                        char number[10]="";
                        memmove(number+0,line+old_pos,c-old_pos);
                        old_pos=c+1;
                        entry.mean.push_back(atof(number));
                    }
                }

                if(entry.mean.size()>0)
                {
                    entries.push_back(entry);
                }
	        }

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if(!other_entry && depth>0.0f)
                {
                    int nr_grains_sum=0;

                    grain_neighbors_bin_profile_depths.push_back(depth);

                    std::vector<float> new_entry;
                    grain_neighbors_bin_profile_values.push_back(new_entry);

                    for(int bin=0; bin<entries[k].mean.size(); bin++)
                    {
                        nr_grains_sum+=entries[k].mean[bin];
                    } 

                    for(int bin=0; bin<entries[k].mean.size(); bin++)
                    {
                        grain_neighbors_bin_profile_values.back().push_back(entries[k].mean[bin]/nr_grains_sum);
                    }
                }
            }

            grain_neighbors_bin_file.close();
            temp_grain_neighbors_bin_file.close();
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

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if (!other_entry && depth>0.0f)
                {
                    dihedral_angles_depths.push_back(depth);
                    dihedral_angles_values.push_back(entries[k].mean);
                    dihedral_angles_errors.push_back(entries[k].stdabw);
                }
            }

            dihedral_angles_file.close();
            temp_dihedral_angles_file.close();
        }

        std::string filepath_dihedral_angles2=path_results;
        filepath_dihedral_angles2.append("dihedral_angles2");
        filepath_dihedral_angles2.append(s.str());

        //dihedral angles 2nd version
        std::ifstream dihedral_angles2_file(filepath_dihedral_angles2.c_str());
        std::ifstream temp_dihedral_angles2_file(filepath_dihedral_angles2.c_str());
        temp_dihedral_angles2_file>>teststring;

        if(dihedral_angles2_file && teststring.size()!=0)
        {
            std::vector<depth_parameter> entries;

            while(!dihedral_angles2_file.eof())
            {
                std::string name;
                float mean=-1, stdabw;
	            dihedral_angles2_file>>name>>mean>>stdabw;
                depth_parameter entry;
                entry.name=name;
                entry.mean=mean;
                entry.stdabw=stdabw;
                if(mean!=-1) entries.push_back(entry);
            }

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if (!other_entry && depth>0.0f)
                {
                    dihedral_angles2_depths.push_back(depth);
                    dihedral_angles2_values.push_back(entries[k].mean);
                    dihedral_angles2_errors.push_back(entries[k].stdabw);
                }
            }

            dihedral_angles2_file.close();
            temp_dihedral_angles2_file.close();
        }

        std::string filepath_disl_dens_longest_boundaries=path_results;
        filepath_disl_dens_longest_boundaries.append("disl_dens_step");
        filepath_disl_dens_longest_boundaries.append(s.str());

        //disl dens longest boundaries
        std::ifstream disl_dens_longest_boundaries_file(filepath_disl_dens_longest_boundaries.c_str());
        std::ifstream temp_disl_dens_longest_boundaries_file(filepath_disl_dens_longest_boundaries.c_str());
        temp_disl_dens_longest_boundaries_file>>teststring;

        if(disl_dens_longest_boundaries_file && teststring.size()!=0)
        {
            std::vector<depth_parameter3> entries;

	        while(!disl_dens_longest_boundaries_file.eof())
	        {
                std::string name;
		        disl_dens_longest_boundaries_file>>name;                    
                depth_parameter3 entry;
                entry.name=name;

                char line[16384];
                disl_dens_longest_boundaries_file.getline(line,16384);

                int old_pos=1;
                for(int c=1; c<disl_dens_longest_boundaries_file.gcount(); c++)
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

                if(entry.mean.size()>0)
                {
                    entries.push_back(entry);
                }
	        }

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if(!other_entry && depth>0.0f)
                {
                    disl_dens_longest_boundaries_profile_depths.push_back(depth);

                    std::vector<float> new_entry;
                    disl_dens_longest_boundaries_profile_values.push_back(new_entry);

                    for(int s=0; s<entries[k].mean.size(); s++)
                    {
                        if(s+1>disl_dens_longest_boundaries_depths.size()) disl_dens_longest_boundaries_depths.resize(s+1);
                        if(s+1>disl_dens_longest_boundaries_values.size()) disl_dens_longest_boundaries_values.resize(s+1);
                        if(s+1>disl_dens_longest_boundaries_errors.size()) disl_dens_longest_boundaries_errors.resize(s+1);
                        disl_dens_longest_boundaries_depths[s].push_back(depth);
                        disl_dens_longest_boundaries_values[s].push_back(entries[k].mean[s]);
                        disl_dens_longest_boundaries_errors[s].push_back(entries[k].stdabw[s]);

                        disl_dens_longest_boundaries_profile_values.back().push_back(entries[k].mean[s]);
                    } 
                }
            }

            disl_dens_longest_boundaries_file.close();
            temp_disl_dens_longest_boundaries_file.close();
        }

        std::string filepath_disl_dens_percent_boundaries=path_results;
        filepath_disl_dens_percent_boundaries.append("disl_dens_percent");
        filepath_disl_dens_percent_boundaries.append(s.str());

        //disl dens percent boundaries
        std::ifstream disl_dens_percent_boundaries_file(filepath_disl_dens_percent_boundaries.c_str());
        std::ifstream temp_disl_dens_percent_boundaries_file(filepath_disl_dens_percent_boundaries.c_str());
        temp_disl_dens_percent_boundaries_file>>teststring;

        if(disl_dens_percent_boundaries_file && teststring.size()!=0)
        {
            std::vector<depth_parameter3> entries;

	        while(!disl_dens_percent_boundaries_file.eof())
	        {
                std::string name;
		        disl_dens_percent_boundaries_file>>name;
                depth_parameter3 entry;
                entry.name=name;

                float mean[10];
                float stdabw[10];
                bool end=false;

                for (int p=0; p<10 && !end; p++)
                {
                    mean[p]=-1;
    				    disl_dens_percent_boundaries_file>>mean[p]>>stdabw[p];
                    if(mean[p]!=-1)
                    {
                        entry.mean.push_back(mean[p]);
                        entry.stdabw.push_back(stdabw[p]);
                    }
                    else end=true;
                }

                if(!end) entries.push_back(entry);
	        }

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if(!other_entry && depth>0.0f)
                {
                    disl_dens_percent_boundaries_profile_depths.push_back(depth); 

                    std::vector<float> new_entry;
                    disl_dens_percent_boundaries_profile_values.push_back(new_entry);

                    for(int p=0; p<10; p++)
                    {
                        disl_dens_percent_boundaries_depths[p].push_back(depth);
                        disl_dens_percent_boundaries_values[p].push_back(entries[k].mean[p]);
                        disl_dens_percent_boundaries_errors[p].push_back(entries[k].stdabw[p]);

                        disl_dens_percent_boundaries_profile_values.back().push_back(entries[k].mean[p]);
                    } 
                }
            }

            disl_dens_percent_boundaries_file.close();
            temp_disl_dens_percent_boundaries_file.close();
        }

        std::string filepath_percentage_filled_grains=path_results;
        filepath_percentage_filled_grains.append("percentage_filled_grains");
        filepath_percentage_filled_grains.append(s.str());

        //percentage filled grains
        std::ifstream percentage_filled_grains_file(filepath_percentage_filled_grains.c_str());
        std::ifstream temp_percentage_filled_grains_file(filepath_percentage_filled_grains.c_str());
        temp_percentage_filled_grains_file>>teststring;

        if(percentage_filled_grains_file && teststring.size()!=0)
        {
            std::vector<depth_parameter3> entries;

	        while(!percentage_filled_grains_file.eof())
	        {
                std::string name;

		        percentage_filled_grains_file>>name;
                depth_parameter3 entry;
                entry.name=name;

                float mean[13];
                float stdabw[13];
                bool end=false;

                for (int p=0; p<13 && !end; p++)
                {
                    mean[p]=-1;
    				    percentage_filled_grains_file>>mean[p]>>stdabw[p];
                    if(mean[p]!=-1)
                    {
                        entry.mean.push_back(mean[p]);
                    }
                    else end=true;
                }

                if(!end) entries.push_back(entry);
	        }

            //check for multiple entries for same image
            for (int k=0; k<entries.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<entries.size() && !other_entry; j++)
                {
                    if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                }

                float depth=get_depth(entries[k].name);

                if (!other_entry && depth>0.0f)
                {
                    percentage_filled_grains_depths[0].push_back(depth);
                    percentage_filled_grains_values[0].push_back(entries[k].mean[1]);

                    percentage_filled_grains_depths[1].push_back(depth);
                    percentage_filled_grains_values[1].push_back(entries[k].mean[2]);

                    percentage_filled_grains_depths[2].push_back(depth);
                    percentage_filled_grains_values[2].push_back(entries[k].mean[7]);

                    percentage_filled_grains_depths[3].push_back(depth);
                    percentage_filled_grains_values[3].push_back(entries[k].mean[12]);
                }
            }

            percentage_filled_grains_file.close();
            temp_percentage_filled_grains_file.close();
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
        std::cout<<"Number of grain neighbors (bin depth profile): "<<grain_neighbors_bin_profile_values.size()<<std::endl;

        for (int i=0; i<grainsize_combined_values.size(); i++)
        {
            std::cout<<"Mean grain size (parameter "<<i+1<<"/"<<grainsize_combined_values.size()<<" in combined profile): "<<
                grainsize_combined_values[i].size()<<std::endl;
        }

        std::cout<<"Grain size fit: "<<grainsize_fit_values.size()<<std::endl;
        std::cout<<"Correlation ellipse box flattening: "<<corr_ellipse_box_flattening_values.size()<<std::endl;
        std::cout<<"Dihedral angle: "<<dihedral_angles_values.size()<<std::endl;
        std::cout<<"Dihedral angle 2nd version: "<<dihedral_angles2_values.size()<<std::endl;

        for (int s=0; s<disl_dens_longest_boundaries_values.size(); s++)
        {
            std::cout<<"Disl dens diff ("<<(s+1)*boundary_step<<" longest boundaries): "<<
                disl_dens_longest_boundaries_values[s].size()<<std::endl;
        }

        for (int p=0; p<disl_dens_percent_boundaries_values.size(); p++)
        {
            std::cout<<"Disl dens diff ("<<(p+1)*10<<" percent longest boundaries): "<<
                disl_dens_percent_boundaries_values[p].size()<<std::endl;
        }

        std::cout<<"Disl dens diff (step depth profile): "<<disl_dens_longest_boundaries_profile_values.size()<<std::endl;
        std::cout<<"Disl dens diff (percent depth profile): "<<disl_dens_percent_boundaries_profile_values.size()<<std::endl;

        for (int i=0; i<percentage_filled_grains_values.size(); i++)
        {
            std::cout<<"Mean grain size percentage filled (parameter "<<i+1<<"/"<<percentage_filled_grains_values.size()<<
                " in combined profile): "<< percentage_filled_grains_values[i].size()<<std::endl;
        }

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

        std::string filepath_grain_neighbors_bin_profile_out=path_results;
        filepath_grain_neighbors_bin_profile_out.append("grain_neighbors_bin_profile");
        filepath_grain_neighbors_bin_profile_out.append(s.str());
        filepath_grain_neighbors_bin_profile_out.resize(filepath_grain_neighbors_bin_profile_out.size()-3);
        filepath_grain_neighbors_bin_profile_out.append("svg");

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

        std::string filepath_dihedral_angles2_out=path_results;
        filepath_dihedral_angles2_out.append("dihedral_angles2");
        filepath_dihedral_angles2_out.append(s.str());
        filepath_dihedral_angles2_out.resize(filepath_dihedral_angles2_out.size()-3);
        filepath_dihedral_angles2_out.append("svg");

        std::string filepath_disl_dens_longest_boundaries_out=path_results;
        filepath_disl_dens_longest_boundaries_out.append("step/disl_dens_");

        std::string filepath_disl_dens_percent_boundaries_out=path_results;
        filepath_disl_dens_percent_boundaries_out.append("percent/disl_dens_");

        std::string filepath_disl_dens_longest_boundaries_profile_out=path_results;
        filepath_disl_dens_longest_boundaries_profile_out.append("disl_dens_step_profile");
        filepath_disl_dens_longest_boundaries_profile_out.append(s.str());
        filepath_disl_dens_longest_boundaries_profile_out.resize(filepath_disl_dens_longest_boundaries_profile_out.size()-3);
        filepath_disl_dens_longest_boundaries_profile_out.append("svg");

        std::string filepath_disl_dens_percent_boundaries_profile_out=path_results;
        filepath_disl_dens_percent_boundaries_profile_out.append("disl_dens_percent_profile");
        filepath_disl_dens_percent_boundaries_profile_out.append(s.str());
        filepath_disl_dens_percent_boundaries_profile_out.resize(filepath_disl_dens_percent_boundaries_profile_out.size()-3);
        filepath_disl_dens_percent_boundaries_profile_out.append("svg");

        std::string filepath_grainsize_percentage_filled_out=path_results;
        filepath_grainsize_percentage_filled_out.append("grainsize_percentage_filled");
        filepath_grainsize_percentage_filled_out.append(s.str());
        filepath_grainsize_percentage_filled_out.resize(filepath_grainsize_percentage_filled_out.size()-3);
        filepath_grainsize_percentage_filled_out.append("svg");

        if(grainsize_all_values.size()>0) plot.draw_depth_errors("Depth in m", area_unit, "Mean grain size profile (all grains)",
                                                                 grainsize_all_depths, grainsize_all_values, grainsize_all_errors,
                                                                 grainsize_all_errors, 0.0f, filepath_grainsize_all_out.c_str());

        if(grainsize_fit_values.size()>0) plot.draw_depth_errors("Depth in m", area_unit, "Grain size profile log-normal fit",
                                                                 grainsize_fit_depths, grainsize_fit_values,
                                                                 grainsize_fit_errors_low, grainsize_fit_errors_high, 0.0f,
                                                                 filepath_grainsize_fit_out.c_str());

        if(corr_ellipse_box_flattening_values.size()>0) plot.draw_depth("Depth in m", "Pearson correlation",
                                                                 "Correlation ellipse box flattening profile",
                                                                 corr_ellipse_box_flattening_depths,
                                                                 corr_ellipse_box_flattening_values, 0.0f,
                                                                 filepath_corr_ellipse_box_flattening_out.c_str());

        for (int st=0; st<grainsize_step_values.size(); st++)
        {
            char step[6];
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

            if(grainsize_step_values[st].size()>0) plot.draw_depth_errors("Depth in m", area_unit, titel, grainsize_step_depths[st],
                                                                 grainsize_step_values[st], grainsize_step_errors[st],
                                                                 grainsize_step_errors[st], 0.0f, filepath_out.c_str());
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

            if(grainsize_percent_values[p].size()>0) plot.draw_depth_errors("Depth in m", area_unit, titel,
                                                                 grainsize_percent_depths[p], grainsize_percent_values[p],
                                                                 grainsize_percent_errors[p], grainsize_percent_errors[p], 0.0f,
                                                                 filepath_out.c_str());
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

            if(grainsize_quantile_values[q].size()>0) plot.draw_depth_errors("Depth in m", area_unit, titel,
                                                                 grainsize_quantile_depths[q], grainsize_quantile_values[q],
                                                                 grainsize_quantile_errors[q], grainsize_quantile_errors[q], 0.0f,
                                                                 filepath_out.c_str());
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
            sprintf(size_low, "%2.1e", pow(10.0f,(float)b/grain_bin_width)/area_scaling());
            sprintf(size_high, "%2.1e", pow(10.0f,(float)(b+1)/grain_bin_width)/area_scaling());

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

            if(grainsize_bin_values[b].size()>0) plot.draw_depth("Depth in m", "Relative occurrence", titel, grainsize_bin_depths[b],
                                                                 grainsize_bin_values[b], 0.0f, filepath_out.c_str());
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
            sprintf(radius_low, "%2.1e", pow(10.0f,(float)b/grain_equiv_radius_bin_width)/length_scaling());
            sprintf(radius_high, "%2.1e", pow(10.0f,(float)(b+1)/grain_equiv_radius_bin_width)/length_scaling());

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
                                                                 grainradius_bin_depths[b], grainradius_bin_values[b], 0.0f,
                                                                 filepath_out.c_str());
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
                                                                 grainradius_norm_bin_depths[b], grainradius_norm_bin_values[b], 0.0f,
                                                                 filepath_out.c_str());
        }

        char step[20];
        sprintf(step, "%d", grain_step);

        std::string x_axis="Largest grains * ";
        x_axis.append(step);

        if(grainsize_step_profile_values.size()>0) plot.draw_depth_3d(x_axis, "Depth in m", area_unit,
                                                                 "Largest grains profile (absolute)",
                                                                 grainsize_step_profile_x_values, grainsize_step_profile_depths,
                                                                 grainsize_step_profile_values,
                                                                 filepath_grainsize_step_profile_out.c_str());
        if(grainsize_percent_profile_values.size()>0) plot.draw_depth_3d("Largest grains", "Depth in m", area_unit,
                                                                 "Largest grains profile (relative)",
                                                                 grainsize_percent_profile_x_values, grainsize_percent_profile_depths,
                                                                 grainsize_percent_profile_values,
                                                                 filepath_grainsize_percent_profile_out.c_str());
        if(grainsize_quantile_profile_values.size()>0) plot.draw_depth_3d("Quantile 1-", "Depth in m", area_unit,
                                                                 "Grain size quantiles", grainsize_quantile_profile_x_values,
                                                                 grainsize_quantile_profile_depths,
                                                                 grainsize_quantile_profile_values,
                                                                 filepath_grainsize_quantile_profile_out.c_str());
        if(grainsize_bin_profile_values.size()>0) plot.draw_depth_3d(area_unit, "Depth in m", "Relative occurrence",
                                                                 "Grain size profile", grainsize_bin_profile_x_values,
                                                                 grainsize_bin_profile_depths, grainsize_bin_profile_values,
                                                                 filepath_grainsize_bin_profile_out.c_str());
        if(grainradius_bin_profile_values.size()>0) plot.draw_depth_3d(length_unit, "Depth in m", "Relative occurrence",
                                                                 "Grain radius profile", grainradius_bin_profile_x_values,
                                                                 grainradius_bin_profile_depths, grainradius_bin_profile_values,
                                                                 filepath_grainradius_bin_profile_out.c_str());
        if(grainradius_norm_bin_profile_values.size()>0) plot.draw_depth_3d("Normalized grain radius", "Depth in m",
                                                                 "Relative occurrence", "Normalized grain radius profile",
                                                                 grainradius_norm_bin_profile_x_values,
                                                                 grainradius_norm_bin_profile_depths,
                                                                 grainradius_norm_bin_profile_values,
                                                                 filepath_grainradius_norm_bin_profile_out.c_str());
        if(grain_neighbors_bin_profile_values.size()>0) plot.draw_depth_3d("Number of grain neighbors", "Depth in m",
                                                                 "Relative occurrence", "Number of grain neighbors profile",
                                                                 grain_neighbors_bin_profile_x_values,
                                                                 grain_neighbors_bin_profile_depths,
                                                                 grain_neighbors_bin_profile_values,
                                                                 filepath_grain_neighbors_bin_profile_out.c_str());

        if(grainsize_combined_values.size()>0) plot.draw_depth_nofit("Depth in m", area_unit, "Different grain size parameters",
                                                                 grainsize_combined_depths, grainsize_combined_values,
                                                                 filepath_grainsize_combined_out.c_str());

        if(dihedral_angles_values.size()>0) plot.draw_depth("Depth in m", "Dihedral angle standard deviation [degree]",
                                                                 "Dihedral angle standard deviation profile",
                                                                 dihedral_angles_depths, dihedral_angles_errors, 10.0f,
                                                                 filepath_dihedral_angles_out.c_str());
        if(dihedral_angles2_values.size()>0) plot.draw_depth("Depth in m", "Dihedral angle standard deviation [degree]",
                                                                "Dihedral angle (2nd version) standard deviation profile",
                                                                 dihedral_angles2_depths, dihedral_angles2_errors, 10.0f,
                                                                 filepath_dihedral_angles2_out.c_str());

        for (int st=0; st<disl_dens_longest_boundaries_values.size(); st++)
        {
            char step[6];
            sprintf(step, "%i", (st+1)*boundary_step);

            std::string titel="Disl. dens. diff. at longest boundaries profile (";
            titel.append(step);
            titel.append(" longest boundaries)");

            std::string filepath_out=filepath_disl_dens_longest_boundaries_out;
            filepath_out.append(step);
            filepath_out.append("step");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(disl_dens_longest_boundaries_values[st].size()>0) plot.draw_depth_errors("Depth in m",
                                                                 "Dislocation density difference in m^(-2)", titel,
                                                                 disl_dens_longest_boundaries_depths[st],
                                                                 disl_dens_longest_boundaries_values[st],
                                                                 disl_dens_longest_boundaries_errors[st],
                                                                 disl_dens_longest_boundaries_errors[st], 0.0f,
                                                                 filepath_out.c_str());
        }

        for (int p=0; p<disl_dens_percent_boundaries_values.size(); p++)
        {
            char percent[4];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Disl. dens. diff. at longest boundaries profile (";
            titel.append(percent);
            titel.append(" percent longest boundaries)");

            std::string filepath_out=filepath_disl_dens_percent_boundaries_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(disl_dens_percent_boundaries_values[p].size()>0) plot.draw_depth_errors("Depth in m",
                                                                 "Dislocation density difference in m^(-2)", titel,
                                                                disl_dens_percent_boundaries_depths[p],
                                                                disl_dens_percent_boundaries_values[p],
                                                                disl_dens_percent_boundaries_errors[p],
                                                                disl_dens_percent_boundaries_errors[p], 0.0f,
                                                                filepath_out.c_str());
        }

        sprintf(step, "%d", boundary_step);

        x_axis="Longest boundaries * ";
        x_axis.append(step);

        if(disl_dens_longest_boundaries_profile_values.size()>0) plot.draw_depth_3d(x_axis, "Depth in m",
                                                                 "Dislocation density difference in m^(-2)",
                                                                 "Disl. dens. diff. at longest boundaries profile (absolute)",
                                                                 disl_dens_longest_boundaries_profile_x_values,
                                                                 disl_dens_longest_boundaries_profile_depths,
                                                                 disl_dens_longest_boundaries_profile_values,
                                                                 filepath_disl_dens_longest_boundaries_profile_out.c_str());
        if(disl_dens_percent_boundaries_profile_values.size()>0) plot.draw_depth_3d("Longest boundaries", "Depth in m",
                                                                 "Dislocation density difference in m^(-2)",
                                                                 "Disl. dens. diff. at longest boundaries profile (relative)",
                                                                 disl_dens_percent_boundaries_profile_x_values,
                                                                 disl_dens_percent_boundaries_profile_depths,
                                                                 disl_dens_percent_boundaries_profile_values,
                                                                 filepath_disl_dens_percent_boundaries_profile_out.c_str());

        if(percentage_filled_grains_values.size()>0) plot.draw_depth_nofit("Depth in m", area_unit,
                                                                 "Mean grain size for different percentage of area filled by grains",
                                                                 percentage_filled_grains_depths, percentage_filled_grains_values,
                                                                 filepath_grainsize_percentage_filled_out.c_str());

        if(depth_bin_width>0.0f)
        {
            std::cout<<"New plots generated"<<std::endl;

            //load nr values 
            std::vector<depth_numbers> nr_values;

            std::string filepath_nr=path_results;
            filepath_nr.append("nr_values");
            filepath_nr.append(s.str());

            std::ifstream nr_file(filepath_nr.c_str());
            std::ifstream temp_nr_file(filepath_nr.c_str());
            temp_nr_file>>teststring;

            if(nr_file && teststring.size()!=0)
            {
                while(!nr_file.eof())
                {
                    std::string name;
                    int nr_grains=-1;
                    int nr_not_deg_ellipse;
                    int nr_correct_flatt;
                    int nr_correct_long;
                    int nr_correct_length;
                    int nr_correct_longest;
                    int nr_angles;
                    int nr_curved_pixels;
                    int nr_boundaries;
                    int nr_junctions;

                    nr_file>>name>>nr_grains>>nr_not_deg_ellipse>>nr_correct_flatt>>nr_correct_long>>nr_correct_length>>nr_correct_longest>>
                        nr_angles>>nr_curved_pixels>>nr_boundaries>>nr_junctions;
                    depth_numbers entry;
                    entry.name=name;
                    entry.nr[0]=nr_grains;
                    entry.nr[1]=nr_not_deg_ellipse;
                    entry.nr[2]=nr_correct_flatt;
                    entry.nr[3]=nr_correct_long;
                    entry.nr[4]=nr_correct_length;
                    entry.nr[5]=nr_correct_longest;
                    entry.nr[6]=nr_angles;
                    entry.nr[7]=nr_curved_pixels;
                    entry.nr[8]=nr_boundaries;
                    entry.nr[9]=nr_junctions;
                    if(nr_grains!=-1) nr_values.push_back(entry);
                }

                nr_file.close();
                temp_nr_file.close();
            }

            //number histograms
            std::vector<int> nr_grains_histogram(1,0);
            std::vector<int> nr_boundaries_histogram(1,0);
            std::vector<int> nr_junctions_histogram(1,0);
            std::vector<int> nr_images_histogram(1,0);
            float depth_max=0.0f;

            //check for multiple entries for same image
            for (int k=0; k<nr_values.size(); k++)
            {
                bool other_entry=false;

                for(int j=k+1; j<nr_values.size() && !other_entry; j++)
                {
                    if(nr_values[k].name.compare(nr_values[j].name)==0) other_entry=true;
                }

                float depth=get_depth(nr_values[k].name);

                if (!other_entry && depth>0.0f)
                {
                    if (depth>depth_max)
                    {
                        depth_max=depth;

                        nr_grains_histogram.resize(1+depth_max/depth_bin_width,0);
                        nr_boundaries_histogram.resize(1+depth_max/depth_bin_width,0);
                        nr_junctions_histogram.resize(1+depth_max/depth_bin_width,0);
                        nr_images_histogram.resize(1+depth_max/depth_bin_width,0);
                    }
                    
                    nr_grains_histogram[depth/depth_bin_width]+=nr_values[k].nr[0];
                    nr_boundaries_histogram[depth/depth_bin_width]+=nr_values[k].nr[8];
                    nr_junctions_histogram[depth/depth_bin_width]+=nr_values[k].nr[9];
                    nr_images_histogram[depth/depth_bin_width]++;
                }
            }

            std::string filepath_nr_grains=path_results;
            std::string filepath_nr_boundaries=path_results;
            std::string filepath_nr_junctions=path_results;
            std::string filepath_nr_images=path_results;

            filepath_nr_grains.append("nr_grains");
            filepath_nr_boundaries.append("nr_boundaries");
            filepath_nr_junctions.append("nr_junctions");
            filepath_nr_images.append("nr_images.svg");

            filepath_nr_grains.append(s.str());
            filepath_nr_boundaries.append(s.str());
            filepath_nr_junctions.append(s.str());

            filepath_nr_grains.resize(filepath_nr_grains.size()-3);
            filepath_nr_boundaries.resize(filepath_nr_boundaries.size()-3);
            filepath_nr_junctions.resize(filepath_nr_junctions.size()-3);

            filepath_nr_grains.append("svg");
            filepath_nr_boundaries.append("svg");
            filepath_nr_junctions.append("svg");

            plot.draw_histogram("Depth in m", "Relative occurence", "Grain depth distribution", nr_grains_histogram, depth_bin_width,
                                1.0f, 0.0f, filepath_nr_grains);
            plot.draw_histogram("Depth in m", "Relative occurence", "Boundary depth distribution", nr_boundaries_histogram,
                                depth_bin_width, 1.0f, 0.0f, filepath_nr_boundaries);
            plot.draw_histogram("Depth in m", "Relative occurence", "Junctions depth distribution", nr_junctions_histogram,
                                depth_bin_width, 1.0f, 0.0f, filepath_nr_junctions);
            plot.draw_histogram("Depth in m", "Relative occurence", "Images depth distribution", nr_images_histogram,
                                depth_bin_width, 1.0f, 0.0f,filepath_nr_images);

            //************
            //New profiles
            //************

            //load nr values2 
            std::vector<depth_numbers> nr_values2;

            std::string filepath_nr2=path_results;
            filepath_nr2.append("nr_percentage_filled_grains");
            filepath_nr2.append(s.str());

            std::ifstream nr2_file(filepath_nr2.c_str());
            std::ifstream temp_nr2_file(filepath_nr2.c_str());
            temp_nr2_file>>teststring;

            if(nr2_file && teststring.size()!=0)
            {
                while(!nr2_file.eof())
                {
                    std::string name;
                    int nr_70=-1;
                    int nr_80;
                    int nr_90;
                    int nr_91;
                    int nr_92;
                    int nr_93;
                    int nr_94;
                    int nr_95;
                    int nr_96;
                    int nr_97;
                    int nr_98;
                    int nr_99;
                    int nr_100;

                    nr2_file>>name>>nr_70>>nr_80>>nr_90>>nr_91>>nr_92>>nr_93>>nr_94>>nr_95>>nr_96>>nr_97>>nr_98>>nr_99>>nr_100;
                    depth_numbers entry;
                    entry.name=name;
                    entry.nr2[0]=nr_70;
                    entry.nr2[1]=nr_80;
                    entry.nr2[2]=nr_90;
                    entry.nr2[3]=nr_91;
                    entry.nr2[4]=nr_92;
                    entry.nr2[5]=nr_93;
                    entry.nr2[6]=nr_94;
                    entry.nr2[7]=nr_95;
                    entry.nr2[8]=nr_96;
                    entry.nr2[9]=nr_97;
                    entry.nr2[10]=nr_98;
                    entry.nr2[11]=nr_99;
                    entry.nr2[12]=nr_100;
                    if(nr_70!=-1) nr_values2.push_back(entry);
                }

                nr2_file.close();
                temp_nr2_file.close();
            }

            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "grainshape", "Grain shape", "Depth in m", "Roundness factor (4pi*Area/Perimeter^2)", "Grain roundness profile",
                0.28f, 1, nr_values, depth_max, depth_bin_width);
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "nr_grainarcs", "Number of grain boundaries", "Depth in m", "Nr of grain boundaries",
                "Number of grain boundaries profile", 0.0f, 0, nr_values, depth_max, depth_bin_width);
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "nr_neighbors", "Number of grain neighbors", "Depth in m", "Nr of grain neighbors",
                "Number of grain neighbors profile", 0.0f, 0, nr_values, depth_max, depth_bin_width);
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "length_grainarcs", "Grain arc length", "Depth in m", length_unit, "Grain boundary length profile", 0.0f, 5,
                nr_values, depth_max, depth_bin_width, length_scaling());
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "longest_grainarcs", "Grain longest arc length", "Depth in m", length_unit, "Grain longest boundary length profile",
                0.0f, 6, nr_values, depth_max, depth_bin_width, length_scaling());
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "dislocation_densities", "Dislocation density differences", "Depth in m", "Dislocation density difference in m^(-2)",
                "Profile of mean dislocation density difference at grain boundaries", 0.0f, 8, nr_values, depth_max, depth_bin_width);
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "grain_equivradius", "Equivalent grain radius", "Depth in m", length_unit, "Grain equivalent radius profile", 0.0f,
                1, nr_values, depth_max, depth_bin_width, length_scaling());
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "grain_boxshape", "Grain box shape", "Depth in m", "Box flattening factor", "Box flattening profile", 0.9f, 2,
                nr_values, depth_max, depth_bin_width);
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "grain_boxwidth", "Grain box width", "Depth in m", length_unit, "Grain box width profile", 0.0f, 2, nr_values,
                depth_max, depth_bin_width, length_scaling());
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "grain_boxheight", "Grain box height", "Depth in m", length_unit, "Grain box height profile", 0.0f, 2, nr_values,
                depth_max, depth_bin_width, length_scaling());
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "grain_ellipselong", "Grain ellipse long axis", "Depth in m", length_unit, "Ellipse long axis profile", 0.0f, 4,
                nr_values, depth_max, depth_bin_width, length_scaling());
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "grain_ellipseshape", "Grain ellipse shape", "Depth in m", "Ellipse flattening factor",
                "Ellipse flattening profile", 0.9f, 3, nr_values, depth_max, depth_bin_width);
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "grain_ellipseangle", "Grain ellipse orientation", "Depth in m", "Ellipse long axis angle",
                "Ellipse long axis angle profile", 0.0f, 2, nr_values, depth_max, depth_bin_width);
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "boundary_orientation", "Boundary orientation", "Depth in m", "Boundary orientation",
                "Boundary orientation profile", 0.0f, 7, nr_values, depth_max, depth_bin_width);
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "grainwidth", "Grain width", "Depth in m", length_unit, "Grain width profile", 0.0f, 1, nr_values, depth_max,
                depth_bin_width, length_scaling());
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "grainheight", "Grain height", "Depth in m", length_unit, "Grain height profile", 0.0f, 1, nr_values, depth_max,
                depth_bin_width, length_scaling());
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "grainflattening", "Vertical grain flattening", "Depth in m", "Vertical grain flattening factor",
                "Vertical grain flattening profile", 0.9f, 1, nr_values, depth_max, depth_bin_width);
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "turning_points", "Number of turning points", "Depth in m", "Nr of turning points",
                "Grain boundary turning points profile", 0.0f, 9, nr_values, depth_max, depth_bin_width);
            single_new_lin_depth_profile(path_results, minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(), plot,
                "grain_perimeter_ratio", "Grain perimeter ratio", "Depth in m", "Perimeter ratio", "Perimeter ratio profile", 0.7f, 1,
                nr_values, depth_max, depth_bin_width);

            std::vector<float> grainsize_all_depths((size_t)(1+depth_max/depth_bin_width));
            std::vector<float> grainsize_all_values((size_t)(1+depth_max/depth_bin_width),0.0f);
            std::vector<float> grainsize_all_sum((size_t)(1+depth_max/depth_bin_width),0.0f);

            std::vector< std::vector<float> > grainsize_step_depths;
            std::vector< std::vector<float> > grainsize_step_values;
            std::vector< std::vector<float> > grainsize_step_sum(500);
            for(int st=0; st<500; st++) grainsize_step_sum[st].resize((size_t)(1+depth_max/depth_bin_width),0.0f);

            std::vector< std::vector<float> > grainsize_bin_depths;
            std::vector< std::vector<float> > grainsize_bin_values;
            std::vector<int> grainsize_bin_sum((size_t)(1+depth_max/depth_bin_width),0);

            std::vector< std::vector<float> > grainradius_bin_depths;
            std::vector< std::vector<float> > grainradius_bin_values;
            std::vector<int> grainradius_bin_sum((size_t)(1+depth_max/depth_bin_width),0);
            
            std::vector< std::vector<float> > grainradius_norm_bin_depths;
            std::vector< std::vector<float> > grainradius_norm_bin_values;
            std::vector<int> grainradius_norm_bin_sum((size_t)(1+depth_max/depth_bin_width),0);

            std::vector<float> grainsize_step_profile_depths((size_t)(1+depth_max/depth_bin_width));
            std::vector<float> grainsize_step_profile_x_values;
            for(int st=0; st<500; st++) grainsize_step_profile_x_values.push_back(st+1);
            std::vector< std::vector<float> > grainsize_step_profile_values((size_t)(1+depth_max/depth_bin_width));

            std::vector<float> grainsize_bin_profile_depths((size_t)(1+depth_max/depth_bin_width));
            std::vector<float> grainsize_bin_profile_x_values;
            for(int bin=0; bin<50; bin++) grainsize_bin_profile_x_values.push_back
                (((float)bin/grain_bin_width)-log10(area_scaling()));
            std::vector< std::vector<float> > grainsize_bin_profile_values((size_t)(1+depth_max/depth_bin_width));

            std::vector<float> grainradius_bin_profile_depths((size_t)(1+depth_max/depth_bin_width));
            std::vector<float> grainradius_bin_profile_x_values;
            for(int bin=0; bin<50; bin++)
                grainradius_bin_profile_x_values.push_back(((float)bin/grain_equiv_radius_bin_width)-log10(length_scaling()));
            std::vector< std::vector<float> > grainradius_bin_profile_values((size_t)(1+depth_max/depth_bin_width));

            std::vector<float> grainradius_norm_bin_profile_depths((size_t)(1+depth_max/depth_bin_width));
            std::vector<float> grainradius_norm_bin_profile_x_values;
            for(int bin=0; bin<100; bin++) grainradius_norm_bin_profile_x_values.push_back
                ((float)bin*grain_equiv_radius_norm_bin_width);
            std::vector< std::vector<float> > grainradius_norm_bin_profile_values((size_t)(1+depth_max/depth_bin_width));

            std::vector<int> grain_neighbors_bin_sum((size_t)(1+depth_max/depth_bin_width),0);
            std::vector<float> grain_neighbors_bin_profile_depths((size_t)(1+depth_max/depth_bin_width));
            std::vector<float> grain_neighbors_bin_profile_x_values;
            for(int bin=0; bin<100; bin++) grain_neighbors_bin_profile_x_values.push_back((float)bin);
            std::vector< std::vector<float> > grain_neighbors_bin_profile_values((size_t)(1+depth_max/depth_bin_width));

            std::vector<float> dihedral_angles_depths((size_t)(1+depth_max/depth_bin_width));
            std::vector<float> dihedral_angles_errors((size_t)(1+depth_max/depth_bin_width),0.0f);
            std::vector<float> dihedral_angles_sum((size_t)(1+depth_max/depth_bin_width),0.0f);

            std::vector<float> dihedral_angles2_depths((size_t)(1+depth_max/depth_bin_width));
            std::vector<float> dihedral_angles2_errors((size_t)(1+depth_max/depth_bin_width),0.0f);
            std::vector<float> dihedral_angles2_sum((size_t)(1+depth_max/depth_bin_width),0.0f);

            std::vector< std::vector<float> > disl_dens_longest_boundaries_depths;
            std::vector< std::vector<float> > disl_dens_longest_boundaries_values;
            std::vector< std::vector<float> > disl_dens_longest_boundaries_sum(500);
            for(int st=0; st<500; st++) disl_dens_longest_boundaries_sum[st].resize((size_t)(1+depth_max/depth_bin_width),0.0f);

            std::vector<float> disl_dens_longest_boundaries_profile_depths((size_t)(1+depth_max/depth_bin_width));
            std::vector<float> disl_dens_longest_boundaries_profile_x_values;
            for(int st=0; st<500; st++) disl_dens_longest_boundaries_profile_x_values.push_back(st+1);
            std::vector< std::vector<float> > disl_dens_longest_boundaries_profile_values((size_t)(1+depth_max/depth_bin_width));

            std::vector< std::vector<float> > percentage_filled_grains_depths;
            percentage_filled_grains_depths.resize(4);
            for(int p=0; p<4; p++) percentage_filled_grains_depths[p].resize((size_t)(1+depth_max/depth_bin_width));
            std::vector< std::vector<float> > percentage_filled_grains_values;
            percentage_filled_grains_values.resize(4);
            for(int p=0; p<4; p++) percentage_filled_grains_values[p].resize((size_t)(1+depth_max/depth_bin_width),0.0f);
            std::vector< std::vector<float> > percentage_filled_grains_sum;
            percentage_filled_grains_sum.resize(4);
            for(int p=0; p<4; p++) percentage_filled_grains_sum[p].resize((size_t)(1+depth_max/depth_bin_width),0.0f);

            std::string filepath_grainsize_all=path_results;
            filepath_grainsize_all.append("grainsize_all");
            filepath_grainsize_all.append(s.str());

            std::string filepath_grainsize_step=path_results;
            filepath_grainsize_step.append("grainsize_step");
            filepath_grainsize_step.append(s.str());

            std::string filepath_grainsize_bin=path_results;
            filepath_grainsize_bin.append("grainsize_bin");
            filepath_grainsize_bin.append(s.str());

            std::string filepath_grainradius_bin=path_results;
            filepath_grainradius_bin.append("grainradius_bin");
            filepath_grainradius_bin.append(s.str());

            std::string filepath_grainradius_norm_bin=path_results;
            filepath_grainradius_norm_bin.append("grainradius_norm_bin");
            filepath_grainradius_norm_bin.append(s.str());

            std::string filepath_grain_neighbors_bin=path_results;
            filepath_grain_neighbors_bin.append("grain_neighbors_bin");
            filepath_grain_neighbors_bin.append(s.str());

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

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    int depth_bin=get_depth(entries[k].name)/depth_bin_width;

                    if(!other_entry && depth_bin>0 && entries[k].name.compare(nr_values[k].name)==0)//nr of values found
                    {
                        grainsize_all_values[depth_bin]+=entries[k].mean*(float)nr_values[k].nr[0]/area_scaling();
                        grainsize_all_sum[depth_bin]+=nr_values[k].nr[0];
                    }
                }

                grainsize_all_file.close();
                temp_grainsize_all_file.close();
            }

            //grain size step
            std::ifstream grainsize_step_file(filepath_grainsize_step.c_str());
            std::ifstream temp_grainsize_step_file(filepath_grainsize_step.c_str());
            temp_grainsize_step_file>>teststring;

            if(grainsize_step_file && teststring.size()!=0)
            {
                std::vector<depth_parameter3> entries;
                int min_step_nr=1000;

	            while(!grainsize_step_file.eof())
	            {
                    std::string name;
		            grainsize_step_file>>name;                    
                    depth_parameter3 entry;
                    entry.name=name;

                    char line[8192];
                    grainsize_step_file.getline(line,8192);

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

                    if(entry.mean.size()>0)
                    {
                        entries.push_back(entry);
                        if(entry.mean.size()<min_step_nr) min_step_nr=entry.mean.size();
                    }
	            }

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    int depth_bin=get_depth(entries[k].name)/depth_bin_width;

                    if(!other_entry && depth_bin>0)
                    {
                        if (min_step_nr>grainsize_step_values.size()) grainsize_step_values.resize(min_step_nr);
                        for(int s=0; s<min_step_nr; s++)
                        {
                            grainsize_step_values[s].resize((size_t)(1+depth_max/depth_bin_width));
                            grainsize_step_values[s][depth_bin]+=entries[k].mean[s]/area_scaling();

                            grainsize_step_sum[s][depth_bin]++;
                        }

                        if (min_step_nr>grainsize_step_profile_values[depth_bin].size())
                            grainsize_step_profile_values[depth_bin].resize(min_step_nr);

                        for(int s=0; s<min_step_nr; s++)
                        {
                            grainsize_step_profile_values[depth_bin][s]+=entries[k].mean[s]/area_scaling();
                        }
                    }
                }

                grainsize_step_file.close();
                temp_grainsize_step_file.close();
            }

            //grain size bin
            std::ifstream grainsize_bin_file(filepath_grainsize_bin.c_str());
            std::ifstream temp_grainsize_bin_file(filepath_grainsize_bin.c_str());
            temp_grainsize_bin_file>>teststring;

            if(grainsize_bin_file && teststring.size()!=0)
            {
                std::vector<depth_parameter3> entries;
                int min_bin_nr=1000;

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

                    if(entry.mean.size()>0)
                    {
                        entries.push_back(entry);
                        if(entry.mean.size()<min_bin_nr) min_bin_nr=entry.mean.size();
                    }
	            }

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    int depth_bin=get_depth(entries[k].name)/depth_bin_width;

                    if(!other_entry && depth_bin>0)
                    {
                        if (min_bin_nr>grainsize_bin_values.size()) grainsize_bin_values.resize(min_bin_nr);
                        for(int bin=0; bin<min_bin_nr; bin++)
                        {
                            grainsize_bin_values[bin].resize((size_t)(1+depth_max/depth_bin_width));
                            grainsize_bin_values[bin][depth_bin]+=entries[k].mean[bin];

                            grainsize_bin_sum[depth_bin]+=entries[k].mean[bin];
                        }

                        if (min_bin_nr>grainsize_bin_profile_values[depth_bin].size())
                            grainsize_bin_profile_values[depth_bin].resize(min_bin_nr);

                        for(int bin=0; bin<min_bin_nr; bin++)
                        {
                            grainsize_bin_profile_values[depth_bin][bin]+=entries[k].mean[bin];
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
                int min_bin_nr=1000;

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

                    if(entry.mean.size()>0)
                    {
                        entries.push_back(entry);
                        if(entry.mean.size()<min_bin_nr) min_bin_nr=entry.mean.size();
                    }
	            }

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    int depth_bin=get_depth(entries[k].name)/depth_bin_width;

                    if(!other_entry && depth_bin>0)
                    {
                        if (min_bin_nr>grainradius_bin_values.size()) grainradius_bin_values.resize(min_bin_nr);
                        for(int bin=0; bin<min_bin_nr; bin++)
                        {
                            grainradius_bin_values[bin].resize((size_t)(1+depth_max/depth_bin_width));
                            grainradius_bin_values[bin][depth_bin]+=entries[k].mean[bin];

                            grainradius_bin_sum[depth_bin]+=entries[k].mean[bin];
                        }

                        if (min_bin_nr>grainradius_bin_profile_values[depth_bin].size())
                            grainradius_bin_profile_values[depth_bin].resize(min_bin_nr);

                        for(int bin=0; bin<min_bin_nr; bin++)
                        {
                            grainradius_bin_profile_values[depth_bin][bin]+=entries[k].mean[bin];
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
                int min_bin_nr=1000;

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

                    if(entry.mean.size()>0)
                    {
                        entries.push_back(entry);
                        if(entry.mean.size()<min_bin_nr) min_bin_nr=entry.mean.size();
                    }
	            }

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    int depth_bin=get_depth(entries[k].name)/depth_bin_width;

                    if(!other_entry && depth_bin>0)
                    {
                        if (min_bin_nr>grainradius_norm_bin_values.size()) grainradius_norm_bin_values.resize(min_bin_nr);
                        for(int bin=0; bin<min_bin_nr; bin++)
                        {
                            grainradius_norm_bin_values[bin].resize((size_t)(1+depth_max/depth_bin_width));
                            grainradius_norm_bin_values[bin][depth_bin]+=entries[k].mean[bin];

                            grainradius_norm_bin_sum[depth_bin]+=entries[k].mean[bin];
                        }

                        if (min_bin_nr>grainradius_norm_bin_profile_values[depth_bin].size())
                            grainradius_norm_bin_profile_values[depth_bin].resize(min_bin_nr);

                        for(int bin=0; bin<min_bin_nr; bin++)
                        {
                            grainradius_norm_bin_profile_values[depth_bin][bin]+=entries[k].mean[bin];
                        }
                    }
                }

                grainradius_norm_bin_file.close();
                temp_grainradius_norm_bin_file.close();
            }

            //grain neighbors bin
            std::ifstream grain_neighbors_bin_file(filepath_grain_neighbors_bin.c_str());
            std::ifstream temp_grain_neighbors_bin_file(filepath_grain_neighbors_bin.c_str());
            temp_grain_neighbors_bin_file>>teststring;

            if(grain_neighbors_bin_file && teststring.size()!=0)
            {
                std::vector<depth_parameter3> entries;
                int min_bin_nr=1000;

	            while(!grain_neighbors_bin_file.eof())
	            {
                    std::string name;
		            grain_neighbors_bin_file>>name;                    
                    depth_parameter3 entry;
                    entry.name=name;

                    char line[1024];
                    grain_neighbors_bin_file.getline(line,1024);

                    int old_pos=1;
                    for(int c=1; c<grain_neighbors_bin_file.gcount(); c++)
                    {
                        if(line[c]==32)
                        {
                            char number[10]="";
                            memmove(number+0,line+old_pos,c-old_pos);
                            old_pos=c+1;
                            entry.mean.push_back(atof(number));
                        }
                    }

                    if(entry.mean.size()>0)
                    {
                        entries.push_back(entry);
                        if(entry.mean.size()<min_bin_nr) min_bin_nr=entry.mean.size();
                    }
	            }

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    int depth_bin=get_depth(entries[k].name)/depth_bin_width;

                    if(!other_entry && depth_bin>0)
                    {
                        if (min_bin_nr>grain_neighbors_bin_profile_values[depth_bin].size())
                            grain_neighbors_bin_profile_values[depth_bin].resize(min_bin_nr);

                        for(int bin=0; bin<min_bin_nr; bin++)
                        {
                            grain_neighbors_bin_sum[depth_bin]+=entries[k].mean[bin];
                            grain_neighbors_bin_profile_values[depth_bin][bin]+=entries[k].mean[bin];
                        }
                    }
                }

                grain_neighbors_bin_file.close();
                temp_grain_neighbors_bin_file.close();
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
                    entry.stdabw=stdabw;
                    if(mean!=-1) entries.push_back(entry);
                }

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    int depth_bin=get_depth(entries[k].name)/depth_bin_width;

                    if(!other_entry && depth_bin>0 && entries[k].name.compare(nr_values[k].name)==0)//nr of values found
                    {
                        dihedral_angles_errors[depth_bin]+=entries[k].stdabw*nr_values[k].nr[6];
                        dihedral_angles_sum[depth_bin]+=nr_values[k].nr[6];
                    }
                }

                dihedral_angles_file.close();
                temp_dihedral_angles_file.close();
            }

            std::string filepath_dihedral_angles2=path_results;
            filepath_dihedral_angles2.append("dihedral_angles2");
            filepath_dihedral_angles2.append(s.str());

            //dihedral angles 2nd version
            std::ifstream dihedral_angles2_file(filepath_dihedral_angles2.c_str());
            std::ifstream temp_dihedral_angles2_file(filepath_dihedral_angles2.c_str());
            temp_dihedral_angles2_file>>teststring;

            if(dihedral_angles2_file && teststring.size()!=0)
            {
                std::vector<depth_parameter> entries;

                while(!dihedral_angles2_file.eof())
                {
                    std::string name;
                    float mean=-1, stdabw;
	                dihedral_angles2_file>>name>>mean>>stdabw;
                    depth_parameter entry;
                    entry.name=name;
                    entry.stdabw=stdabw;
                    if(mean!=-1) entries.push_back(entry);
                }

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    int depth_bin=get_depth(entries[k].name)/depth_bin_width;

                    if(!other_entry && depth_bin>0 && entries[k].name.compare(nr_values[k].name)==0)//nr of values found
                    {
                        dihedral_angles2_errors[depth_bin]+=entries[k].stdabw*nr_values[k].nr[6];
                        dihedral_angles2_sum[depth_bin]+=nr_values[k].nr[6];
                    }
                }

                dihedral_angles2_file.close();
                temp_dihedral_angles2_file.close();
            }

            std::string filepath_disl_dens_longest_boundaries=path_results;
            filepath_disl_dens_longest_boundaries.append("disl_dens_step");
            filepath_disl_dens_longest_boundaries.append(s.str());

            //disl dens longest boundaries
            std::ifstream disl_dens_longest_boundaries_file(filepath_disl_dens_longest_boundaries.c_str());
            std::ifstream temp_disl_dens_longest_boundaries_file(filepath_disl_dens_longest_boundaries.c_str());
            temp_disl_dens_longest_boundaries_file>>teststring;

            if(disl_dens_longest_boundaries_file && teststring.size()!=0)
            {
                std::vector<depth_parameter3> entries;
                int min_step_nr=1000;

	            while(!disl_dens_longest_boundaries_file.eof())
	            {
                    std::string name;
		            disl_dens_longest_boundaries_file>>name;                    
                    depth_parameter3 entry;
                    entry.name=name;

                    char line[16384];
                    disl_dens_longest_boundaries_file.getline(line,16384);

                    int old_pos=1;
                    for(int c=1; c<disl_dens_longest_boundaries_file.gcount(); c++)
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

                    if(entry.mean.size()>0)
                    {
                        entries.push_back(entry);
                        if(entry.mean.size()<min_step_nr) min_step_nr=entry.mean.size();
                    }
	            }

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    int depth_bin=get_depth(entries[k].name)/depth_bin_width;

                    if(!other_entry && depth_bin>0)
                    {
                        if (min_step_nr>disl_dens_longest_boundaries_values.size())
                            disl_dens_longest_boundaries_values.resize(min_step_nr);
                        for(int s=0; s<min_step_nr; s++)
                        {
                            disl_dens_longest_boundaries_values[s].resize((size_t)(1+depth_max/depth_bin_width));
                            disl_dens_longest_boundaries_values[s][depth_bin]+=entries[k].mean[s];

                            disl_dens_longest_boundaries_sum[s][depth_bin]++;
                        }

                        if (min_step_nr>disl_dens_longest_boundaries_profile_values[depth_bin].size())
                            disl_dens_longest_boundaries_profile_values[depth_bin].resize(min_step_nr);

                        for(int s=0; s<min_step_nr; s++)
                        {
                            disl_dens_longest_boundaries_profile_values[depth_bin][s]+=entries[k].mean[s];
                        }
                    }
                }

                disl_dens_longest_boundaries_file.close();
                temp_disl_dens_longest_boundaries_file.close();
            }

            std::string filepath_percentage_filled_grains=path_results;
            filepath_percentage_filled_grains.append("percentage_filled_grains");
            filepath_percentage_filled_grains.append(s.str());

            //percentage filled grains
            std::ifstream percentage_filled_grains_file(filepath_percentage_filled_grains.c_str());
            std::ifstream temp_percentage_filled_grains_file(filepath_percentage_filled_grains.c_str());
            temp_percentage_filled_grains_file>>teststring;

            if(percentage_filled_grains_file && teststring.size()!=0)
            {
                std::vector<depth_parameter3> entries;

	            while(!percentage_filled_grains_file.eof())
	            {
                    std::string name;

		            percentage_filled_grains_file>>name;
                    depth_parameter3 entry;
                    entry.name=name;

                    float mean[13];
                    float stdabw[13];
                    bool end=false;

                    for (int p=0; p<13 && !end; p++)
                    {
                        mean[p]=-1;
        				    percentage_filled_grains_file>>mean[p]>>stdabw[p];
                        if(mean[p]!=-1)
                        {
                            entry.mean.push_back(mean[p]);
                        }
                        else end=true;
                    }

                    if(!end) entries.push_back(entry);
	            }

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    int depth_bin=get_depth(entries[k].name)/depth_bin_width;

                    if(!other_entry && depth_bin>0 && entries[k].name.compare(nr_values2[k].name)==0)//nr of values found
                    {
                        percentage_filled_grains_values[0][depth_bin]+=entries[k].mean[1]*nr_values2[k].nr2[1];
                        percentage_filled_grains_sum[0][depth_bin]+=nr_values2[k].nr2[1];
                        percentage_filled_grains_values[1][depth_bin]+=entries[k].mean[2]*nr_values2[k].nr2[2];
                        percentage_filled_grains_sum[1][depth_bin]+=nr_values2[k].nr2[2];
                        percentage_filled_grains_values[2][depth_bin]+=entries[k].mean[7]*nr_values2[k].nr2[7];
                        percentage_filled_grains_sum[2][depth_bin]+=nr_values2[k].nr2[7];
                        percentage_filled_grains_values[3][depth_bin]+=entries[k].mean[12]*nr_values2[k].nr2[12];
                        percentage_filled_grains_sum[3][depth_bin]+=nr_values2[k].nr2[12];
                    }
                }

                percentage_filled_grains_file.close();
                temp_percentage_filled_grains_file.close();
            }

            for (int bin=0; bin<grainsize_all_values.size(); bin++)
            {
                grainsize_all_depths[bin]=(bin+0.5f)*depth_bin_width;
            }

            for (int bin=0; bin<grainsize_all_values.size(); bin++)
            {
                if(grainsize_all_sum[bin]>0.0f) grainsize_all_values[bin]/=grainsize_all_sum[bin];
                else
                {
                    grainsize_all_depths.erase(grainsize_all_depths.begin()+bin);
                    grainsize_all_values.erase(grainsize_all_values.begin()+bin);
                    grainsize_all_sum.erase(grainsize_all_sum.begin()+bin);
                    bin--;
                }
            }

            grainsize_step_depths.resize(grainsize_step_values.size());

            for(int grainsize_step=0; grainsize_step<grainsize_step_values.size(); grainsize_step++)
            {
                grainsize_step_depths[grainsize_step].resize(grainsize_step_values[grainsize_step].size());

                for(int bin=0; bin<grainsize_step_values[grainsize_step].size(); bin++)
                {
                    grainsize_step_depths[grainsize_step][bin]=(bin+0.5f)*depth_bin_width;
                }
            }

            for(int grainsize_step=0; grainsize_step<grainsize_step_values.size(); grainsize_step++)
            {
                int count=0;
                for (int bin=0; bin<grainsize_step_values[grainsize_step].size(); bin++)
                {
                    if(grainsize_step_sum[grainsize_step][count]>0.0f)
                        grainsize_step_values[grainsize_step][bin]/=grainsize_step_sum[grainsize_step][count];
                    else
                    {
                        grainsize_step_depths[grainsize_step].erase(grainsize_step_depths[grainsize_step].begin()+bin);
                        grainsize_step_values[grainsize_step].erase(grainsize_step_values[grainsize_step].begin()+bin);
                        bin--;
                    }
                    count++;
                }
            }

            grainsize_bin_depths.resize(grainsize_bin_values.size());

            for(int grainsize_bin=0; grainsize_bin<grainsize_bin_values.size(); grainsize_bin++)
            {
                grainsize_bin_depths[grainsize_bin].resize(grainsize_bin_values[grainsize_bin].size());

                for(int bin=0; bin<grainsize_bin_values[grainsize_bin].size(); bin++)
                {
                    if(grainsize_bin_sum[bin]>0.0f) grainsize_bin_values[grainsize_bin][bin]/=grainsize_bin_sum[bin];
                    grainsize_bin_depths[grainsize_bin][bin]=(bin+0.5f)*depth_bin_width;
                }
            }

            for (int bin=0; bin<grainsize_step_profile_values.size(); bin++)
            {
                grainsize_step_profile_depths[bin]=(bin+0.5f)*depth_bin_width;
            }

            int count=0;
            for (int bin=0; bin<grainsize_step_profile_values.size(); bin++)
            {
                if(grainsize_step_profile_values[bin].size()>0)
                {
                    for(int grainsize_step=0; grainsize_step<grainsize_step_profile_values[bin].size(); grainsize_step++)
                    {
                        if(grainsize_step_sum[grainsize_step][count]>0.0f)
                        {
                            grainsize_step_profile_values[bin][grainsize_step]/=(grainsize_step_sum[grainsize_step][count]);
                        }
                    }
                }
                else
                {
                    grainsize_step_profile_depths.erase(grainsize_step_profile_depths.begin()+bin);
                    grainsize_step_profile_values.erase(grainsize_step_profile_values.begin()+bin);
                    bin--;
                }
                count++;
            }

            for (int bin=0; bin<grainsize_bin_profile_values.size(); bin++)
            {
                grainsize_bin_profile_depths[bin]=(bin+0.5f)*depth_bin_width;
            }

            for (int bin=0; bin<grainsize_bin_profile_values.size(); bin++)
            {
                if(grainsize_bin_sum[bin]>0.0f)
                {
                    for(int grainsize_bin=0; grainsize_bin<grainsize_bin_profile_values[bin].size(); grainsize_bin++)
                        grainsize_bin_profile_values[bin][grainsize_bin]/=grainsize_bin_sum[bin];
                }
                else
                {
                    grainsize_bin_profile_depths.erase(grainsize_bin_profile_depths.begin()+bin);
                    grainsize_bin_profile_values.erase(grainsize_bin_profile_values.begin()+bin);
                    grainsize_bin_sum.erase(grainsize_bin_sum.begin()+bin);
                    bin--;
                }
            }

            grainradius_bin_depths.resize(grainradius_bin_values.size());

            for(int grainradius_bin=0; grainradius_bin<grainradius_bin_values.size(); grainradius_bin++)
            {
                grainradius_bin_depths[grainradius_bin].resize(grainradius_bin_values[grainradius_bin].size());

                for(int bin=0; bin<grainradius_bin_values[grainradius_bin].size(); bin++)
                {
                    if(grainradius_bin_sum[bin]>0.0f) grainradius_bin_values[grainradius_bin][bin]/=grainradius_bin_sum[bin];
                    grainradius_bin_depths[grainradius_bin][bin]=(bin+0.5f)*depth_bin_width;
                }
            }

            for (int bin=0; bin<grainradius_bin_profile_values.size(); bin++)
            {
                grainradius_bin_profile_depths[bin]=(bin+0.5f)*depth_bin_width;
            }

            for (int bin=0; bin<grainradius_bin_profile_values.size(); bin++)
            {
                if(grainradius_bin_sum[bin]>0.0f)
                {
                    for(int grainradius_bin=0; grainradius_bin<grainradius_bin_profile_values[bin].size(); grainradius_bin++)
                        grainradius_bin_profile_values[bin][grainradius_bin]/=grainradius_bin_sum[bin];
                }
                else
                {
                    grainradius_bin_profile_depths.erase(grainradius_bin_profile_depths.begin()+bin);
                    grainradius_bin_profile_values.erase(grainradius_bin_profile_values.begin()+bin);
                    grainradius_bin_sum.erase(grainradius_bin_sum.begin()+bin);
                    bin--;
                }
            }

            grainradius_norm_bin_depths.resize(grainradius_norm_bin_values.size());

            for(int grainradius_norm_bin=0; grainradius_norm_bin<grainradius_norm_bin_values.size(); grainradius_norm_bin++)
            {
                grainradius_norm_bin_depths[grainradius_norm_bin].resize(grainradius_norm_bin_values[grainradius_norm_bin].size());

                for(int bin=0; bin<grainradius_norm_bin_values[grainradius_norm_bin].size(); bin++)
                {
                    if(grainradius_norm_bin_sum[bin]>0.0f)
                        grainradius_norm_bin_values[grainradius_norm_bin][bin]/=grainradius_norm_bin_sum[bin];
                    grainradius_norm_bin_depths[grainradius_norm_bin][bin]=(bin+0.5f)*depth_bin_width;
                }
            }

            for (int bin=0; bin<grainradius_norm_bin_profile_values.size(); bin++)
            {
                grainradius_norm_bin_profile_depths[bin]=(bin+0.5f)*depth_bin_width;
            }

            for (int bin=0; bin<grainradius_norm_bin_profile_values.size(); bin++)
            {
                if(grainradius_norm_bin_sum[bin]>0.0f)
                {
                    for(int grainradius_norm_bin=0; grainradius_norm_bin<grainradius_norm_bin_profile_values[bin].size();
                        grainradius_norm_bin++)
                        grainradius_norm_bin_profile_values[bin][grainradius_norm_bin]/=grainradius_norm_bin_sum[bin];
                }
                else
                {
                    grainradius_norm_bin_profile_depths.erase(grainradius_norm_bin_profile_depths.begin()+bin);
                    grainradius_norm_bin_profile_values.erase(grainradius_norm_bin_profile_values.begin()+bin);
                    grainradius_norm_bin_sum.erase(grainradius_norm_bin_sum.begin()+bin);
                    bin--;
                }
            }

            for (int bin=0; bin<grain_neighbors_bin_profile_values.size(); bin++)
            {
                grain_neighbors_bin_profile_depths[bin]=(bin+0.5f)*depth_bin_width;
            }

            for (int bin=0; bin<grain_neighbors_bin_profile_values.size(); bin++)
            {
                if(grain_neighbors_bin_sum[bin]>0.0f)
                {
                    for(int grain_neighbors_bin=0; grain_neighbors_bin<grain_neighbors_bin_profile_values[bin].size();
                        grain_neighbors_bin++)
                        grain_neighbors_bin_profile_values[bin][grain_neighbors_bin]/=grain_neighbors_bin_sum[bin];
                }
                else
                {
                    grain_neighbors_bin_profile_depths.erase(grain_neighbors_bin_profile_depths.begin()+bin);
                    grain_neighbors_bin_profile_values.erase(grain_neighbors_bin_profile_values.begin()+bin);
                    grain_neighbors_bin_sum.erase(grain_neighbors_bin_sum.begin()+bin);
                    bin--;
                }
            }

            for (int bin=0; bin<dihedral_angles_errors.size(); bin++)
            {
                dihedral_angles_depths[bin]=(bin+0.5f)*depth_bin_width;
            }

            for (int bin=0; bin<dihedral_angles_errors.size(); bin++)
            {
                if(dihedral_angles_sum[bin]>0.0f) dihedral_angles_errors[bin]/=dihedral_angles_sum[bin];
                else
                {
                    dihedral_angles_depths.erase(dihedral_angles_depths.begin()+bin);
                    dihedral_angles_errors.erase(dihedral_angles_errors.begin()+bin);
                    dihedral_angles_sum.erase(dihedral_angles_sum.begin()+bin);
                    bin--;
                }
            }

            for (int bin=0; bin<dihedral_angles2_errors.size(); bin++)
            {
                dihedral_angles2_depths[bin]=(bin+0.5f)*depth_bin_width;
            }

            for (int bin=0; bin<dihedral_angles2_errors.size(); bin++)
            {
                if(dihedral_angles2_sum[bin]>0.0f) dihedral_angles2_errors[bin]/=dihedral_angles2_sum[bin];
                else
                {
                    dihedral_angles2_depths.erase(dihedral_angles2_depths.begin()+bin);
                    dihedral_angles2_errors.erase(dihedral_angles2_errors.begin()+bin);
                    dihedral_angles2_sum.erase(dihedral_angles2_sum.begin()+bin);
                    bin--;
                }
            }

            disl_dens_longest_boundaries_depths.resize(disl_dens_longest_boundaries_values.size());

            for(int step=0; step<disl_dens_longest_boundaries_values.size(); step++)
            {
                disl_dens_longest_boundaries_depths[step].resize(disl_dens_longest_boundaries_values[step].size());

                for(int bin=0; bin<disl_dens_longest_boundaries_values[step].size(); bin++)
                {
                    disl_dens_longest_boundaries_depths[step][bin]=(bin+0.5f)*depth_bin_width;
                }
            }

            for(int step=0; step<disl_dens_longest_boundaries_values.size(); step++)
            {
                int count=0;
                for (int bin=0; bin<disl_dens_longest_boundaries_values[step].size(); bin++)
                {
                    if(disl_dens_longest_boundaries_sum[step][count]>0.0f)
                        disl_dens_longest_boundaries_values[step][bin]/=disl_dens_longest_boundaries_sum[step][count];
                    else
                    {
                        disl_dens_longest_boundaries_depths[step].erase(disl_dens_longest_boundaries_depths[step].begin()+bin);
                        disl_dens_longest_boundaries_values[step].erase(disl_dens_longest_boundaries_values[step].begin()+bin);
                        bin--;
                    }
                    count++;
                }
            }

            for (int bin=0; bin<disl_dens_longest_boundaries_profile_values.size(); bin++)
            {
                disl_dens_longest_boundaries_profile_depths[bin]=(bin+0.5f)*depth_bin_width;
            }

            count=0;
            for (int bin=0; bin<disl_dens_longest_boundaries_profile_values.size(); bin++)
            {
                if(disl_dens_longest_boundaries_profile_values[bin].size()>0)
                {
                    for(int step=0; step<disl_dens_longest_boundaries_profile_values[bin].size(); step++)
                    {
                        if(disl_dens_longest_boundaries_sum[step][count]>0.0f)
                        {
                            disl_dens_longest_boundaries_profile_values[bin][step]/=(disl_dens_longest_boundaries_sum[step][count]);
                        }
                    }
                }
                else
                {
                    disl_dens_longest_boundaries_profile_depths.erase(disl_dens_longest_boundaries_profile_depths.begin()+bin);
                    disl_dens_longest_boundaries_profile_values.erase(disl_dens_longest_boundaries_profile_values.begin()+bin);
                    bin--;
                }
                count++;
            }

            for(int p=0; p<4; p++)
            {
                for (int bin=0; bin<percentage_filled_grains_values[p].size(); bin++)
                {
                    percentage_filled_grains_depths[p][bin]=(bin+0.5f)*depth_bin_width;
                }

                for (int bin=0; bin<percentage_filled_grains_values[p].size(); bin++)
                {
                    if(percentage_filled_grains_sum[p][bin]>0.0f)
                        percentage_filled_grains_values[p][bin]/=percentage_filled_grains_sum[p][bin];
                    else
                    {
                        percentage_filled_grains_depths[p].erase(percentage_filled_grains_depths[p].begin()+bin);
                        percentage_filled_grains_values[p].erase(percentage_filled_grains_values[p].begin()+bin);
                        percentage_filled_grains_sum[p].erase(percentage_filled_grains_sum[p].begin()+bin);
                        bin--;
                    }
                }
            }

            std::string filepath_grainsize_all_out=path_results;
            filepath_grainsize_all_out.append("new_grainsize_all");
            filepath_grainsize_all_out.append(s.str());
            filepath_grainsize_all_out.resize(filepath_grainsize_all_out.size()-3);
            filepath_grainsize_all_out.append("svg");

            std::string filepath_grainsize_all_combined_out=path_results;
            filepath_grainsize_all_combined_out.append("new_grainsize_all_combined");
            filepath_grainsize_all_combined_out.append(s.str());
            filepath_grainsize_all_combined_out.resize(filepath_grainsize_all_combined_out.size()-3);
            filepath_grainsize_all_combined_out.append("svg");

            std::string filepath_grainsize_step_out=path_results;
            filepath_grainsize_step_out.append("step/new_grainsize_");

            std::string filepath_grainsize_bin_out=path_results;
            filepath_grainsize_bin_out.append("bin/new_grainsize_");

            std::string filepath_grainradius_bin_out=path_results;
            filepath_grainradius_bin_out.append("bin/new_grainradius_");

            std::string filepath_grainradius_norm_bin_out=path_results;
            filepath_grainradius_norm_bin_out.append("bin/new_grainradius_norm");

            std::string filepath_grainsize_step_profile_out=path_results;
            filepath_grainsize_step_profile_out.append("new_grainsize_step_profile");
            filepath_grainsize_step_profile_out.append(s.str());
            filepath_grainsize_step_profile_out.resize(filepath_grainsize_step_profile_out.size()-3);
            filepath_grainsize_step_profile_out.append("svg");

            std::string filepath_grainsize_bin_profile_out=path_results;
            filepath_grainsize_bin_profile_out.append("new_grainsize_bin_profile");
            filepath_grainsize_bin_profile_out.append(s.str());
            filepath_grainsize_bin_profile_out.resize(filepath_grainsize_bin_profile_out.size()-3);
            filepath_grainsize_bin_profile_out.append("svg");

            std::string filepath_grainradius_bin_profile_out=path_results;
            filepath_grainradius_bin_profile_out.append("new_grainradius_bin_profile");
            filepath_grainradius_bin_profile_out.append(s.str());
            filepath_grainradius_bin_profile_out.resize(filepath_grainradius_bin_profile_out.size()-3);
            filepath_grainradius_bin_profile_out.append("svg");

            std::string filepath_grainradius_norm_bin_profile_out=path_results;
            filepath_grainradius_norm_bin_profile_out.append("new_grainradius_norm_bin_profile");
            filepath_grainradius_norm_bin_profile_out.append(s.str());
            filepath_grainradius_norm_bin_profile_out.resize(filepath_grainradius_norm_bin_profile_out.size()-3);
            filepath_grainradius_norm_bin_profile_out.append("svg");

            std::string filepath_grain_neighbors_bin_profile_out=path_results;
            filepath_grain_neighbors_bin_profile_out.append("new_grain_neighbors_bin_profile");
            filepath_grain_neighbors_bin_profile_out.append(s.str());
            filepath_grain_neighbors_bin_profile_out.resize(filepath_grain_neighbors_bin_profile_out.size()-3);
            filepath_grain_neighbors_bin_profile_out.append("svg");

            std::string filepath_dihedral_angles_out=path_results;
            filepath_dihedral_angles_out.append("new_dihedral_angles");
            filepath_dihedral_angles_out.append(s.str());
            filepath_dihedral_angles_out.resize(filepath_dihedral_angles_out.size()-3);
            filepath_dihedral_angles_out.append("svg");

            std::string filepath_dihedral_angles2_out=path_results;
            filepath_dihedral_angles2_out.append("new_dihedral_angles2");
            filepath_dihedral_angles2_out.append(s.str());
            filepath_dihedral_angles2_out.resize(filepath_dihedral_angles2_out.size()-3);
            filepath_dihedral_angles2_out.append("svg");

            std::string filepath_disl_dens_longest_boundaries_out=path_results;
            filepath_disl_dens_longest_boundaries_out.append("step/new_disl_dens_");

            std::string filepath_disl_dens_longest_boundaries_profile_out=path_results;
            filepath_disl_dens_longest_boundaries_profile_out.append("new_disl_dens_step_profile");
            filepath_disl_dens_longest_boundaries_profile_out.append(s.str());
            filepath_disl_dens_longest_boundaries_profile_out.resize(filepath_disl_dens_longest_boundaries_profile_out.size()-3);
            filepath_disl_dens_longest_boundaries_profile_out.append("svg");

            std::string filepath_grainsize_percentage_filled_out=path_results;
            filepath_grainsize_percentage_filled_out.append("new_grainsize_percentage_filled");
            filepath_grainsize_percentage_filled_out.append(s.str());
            filepath_grainsize_percentage_filled_out.resize(filepath_grainsize_percentage_filled_out.size()-3);
            filepath_grainsize_percentage_filled_out.append("svg");

            if(grainsize_all_values.size()>0) plot.draw_depth("Depth in m", area_unit, "Mean grain size profile (all grains)",
                grainsize_all_depths, grainsize_all_values, 0.0f, filepath_grainsize_all_out.c_str());

            for (int st=0; st<grainsize_step_values.size(); st++)
            {
                char step[6];
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
                    grainsize_step_values[st], 0.0f, filepath_out.c_str());
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
                sprintf(size_low, "%2.1e", pow(10.0f,(float)b/grain_bin_width)/area_scaling());
                sprintf(size_high, "%2.1e", pow(10.0f,(float)(b+1)/grain_bin_width)/area_scaling());

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
                sprintf(radius_low, "%2.1e", pow(10.0f,(float)b/grain_equiv_radius_bin_width)/length_scaling());
                sprintf(radius_high, "%2.1e", pow(10.0f,(float)(b+1)/grain_equiv_radius_bin_width)/length_scaling());

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

            if(grainsize_step_profile_values.size()>0) plot.draw_depth_3d(x_axis, "Depth in m", area_unit,
                "Largest grains profile (absolute)", grainsize_step_profile_x_values, grainsize_step_profile_depths,
                grainsize_step_profile_values, filepath_grainsize_step_profile_out.c_str());
            if(grainsize_bin_profile_values.size()>0) plot.draw_depth_3d(area_unit, "Depth in m", "Relative occurrence",
                "Grain size profile", grainsize_bin_profile_x_values, grainsize_bin_profile_depths, grainsize_bin_profile_values,
                filepath_grainsize_bin_profile_out.c_str());
            if(grainradius_bin_profile_values.size()>0) plot.draw_depth_3d(length_unit, "Depth in m", "Relative occurrence",
                "Grain radius profile", grainradius_bin_profile_x_values, grainradius_bin_profile_depths,
                grainradius_bin_profile_values, filepath_grainradius_bin_profile_out.c_str());
            if(grainradius_norm_bin_profile_values.size()>0) plot.draw_depth_3d("Normalized grain radius", "Depth in m",
                "Relative occurrence", "Normalized grain radius profile", grainradius_norm_bin_profile_x_values,
                grainradius_norm_bin_profile_depths, grainradius_norm_bin_profile_values,
                filepath_grainradius_norm_bin_profile_out.c_str());
            if(grain_neighbors_bin_profile_values.size()>0) plot.draw_depth_3d("Number of grain neighbors", "Depth in m",
                "Relative occurrence", "Number of grain neighbors profile", grain_neighbors_bin_profile_x_values,
                grain_neighbors_bin_profile_depths, grain_neighbors_bin_profile_values,
                filepath_grain_neighbors_bin_profile_out.c_str());
            if(dihedral_angles_errors.size()>0) plot.draw_depth("Depth in m", "Dihedral angle standard deviation [degree]",
                "Dihedral angle standard deviation profile", dihedral_angles_depths, dihedral_angles_errors,
                 10.0f, filepath_dihedral_angles_out.c_str());
            if(dihedral_angles2_errors.size()>0) plot.draw_depth("Depth in m", "Dihedral angle standard deviation [degree]",
                "Dihedral angle (2nd version) standard deviation profile", dihedral_angles2_depths, dihedral_angles2_errors, 10.0f,
                filepath_dihedral_angles2_out.c_str());

            for (int st=0; st<disl_dens_longest_boundaries_values.size(); st++)
            {
                char step[6];
                sprintf(step, "%i", (st+1)*boundary_step);

                std::string titel="Disl. dens. diff. at longest boundaries profile (";
                titel.append(step);
                titel.append(" longest boundaries)");

                std::string filepath_out=filepath_disl_dens_longest_boundaries_out;
                filepath_out.append(step);
                filepath_out.append("step");
                filepath_out.append(s.str());
                filepath_out.resize(filepath_out.size()-3);
                filepath_out.append("svg");

                if(disl_dens_longest_boundaries_values[st].size()>0) plot.draw_depth("Depth in m",
                    "Dislocation density difference in m^(-2)", titel, disl_dens_longest_boundaries_depths[st],
                    disl_dens_longest_boundaries_values[st], 0.0f, filepath_out.c_str());
            }

            sprintf(step, "%d", boundary_step);

            x_axis="Longest boundaries * ";
            x_axis.append(step);

            if(disl_dens_longest_boundaries_profile_values.size()>0) plot.draw_depth_3d(x_axis, "Depth in m",
                "Dislocation density difference in m^(-2)", "Disl. dens. diff. at longest boundaries profile (absolute)",
                disl_dens_longest_boundaries_profile_x_values, disl_dens_longest_boundaries_profile_depths,
                disl_dens_longest_boundaries_profile_values, filepath_disl_dens_longest_boundaries_profile_out.c_str());
            if(disl_dens_percent_boundaries_profile_values.size()>0) plot.draw_depth_3d("Longest boundaries", "Depth in m",
                "Dislocation density difference in m^(-2)", "Disl. dens. diff. at longest boundaries profile (relative)",
                disl_dens_percent_boundaries_profile_x_values, disl_dens_percent_boundaries_profile_depths,
                disl_dens_percent_boundaries_profile_values, filepath_disl_dens_percent_boundaries_profile_out.c_str());
            if(percentage_filled_grains_values.size()>0) plot.draw_depth_nofit("Depth in m", area_unit,
                "Mean grain size for different percentage of area filled by grains", percentage_filled_grains_depths,
                percentage_filled_grains_values, filepath_grainsize_percentage_filled_out.c_str());
        }
    }
};
