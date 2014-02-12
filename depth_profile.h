/*! \file depth_profile.h
 *  \brief Depth profile calculation.
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
struct depth_parameter
{
    std::string name;
    float mean;
    float stdabw;
    float stdabw2;
};

struct depth_parameter3
{
    std::string name;
    std::vector<float> mean;
    std::vector<float> stdabw;
};

struct depth_numbers
{
    std::string name;
    int nr[10];
    int nr2[13];
};

void lin_depth_profile(std::string path_results, int* bagnr, int nr_depths, int minimal_bubble_distance, int minimal_grain_size,
    int grain_size_min, std::string str, plplot plot, std::string filename, std::string overview_name, std::string x, std::string y,
    std::string title, float y_minimal, int nr_values_selection, std::vector<std::vector<depth_numbers> > nr_values,
    std::vector<std::vector<depth_numbers> > nr_values_relax, float scaling=1.0f)
{
    int * depth = new int[(size_t)nr_depths];

    std::vector<float> parameter_depths;
    std::vector<float> parameter_values;
    std::vector<float> parameter_errors;

    std::vector<float> parameter_depths_relax;
    std::vector<float> parameter_values_relax;
    std::vector<float> parameter_errors_relax;

    for(int i=0; i<nr_depths; i++)
    {
        int iterations=1;
        if (minimal_bubble_distance>0 && i<67) iterations=2;

        depth[i]=0.55*bagnr[i];

        std::stringstream label;
        label << bagnr[i] << "/";

        //string is read from temp file to check whether file is empty
        std::string teststring;

        for(int iter=0; iter<iterations; iter++)
        {
            std::stringstream s;

            if (iter==0) s << "." << minimal_grain_size;
            else s << "." << minimal_bubble_distance << "." << minimal_grain_size;
            if (grain_size_min>10) s << "_" << grain_size_min;
            s << ".txt";

            std::string filepath_parameter=path_results;
            filepath_parameter.append(label.str());
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

                float mean_sum=0;
                float stdabw_sum=0;
                int nr_sum=0;

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    if(!other_entry)
                    {
                        if (iter==0)
                        {
                            if(entries[k].name.compare(nr_values[i][k].name)==0 && nr_values_selection>0)//nr of values found
                            {
                                mean_sum+=entries[k].mean*nr_values[i][k].nr[nr_values_selection-1];
                                stdabw_sum+=entries[k].stdabw*nr_values[i][k].nr[nr_values_selection-1];
                                nr_sum+=nr_values[i][k].nr[nr_values_selection-1];
                            }
                            else
                            {
                                mean_sum+=entries[k].mean;
                                stdabw_sum+=entries[k].stdabw;
                                nr_sum++;
                            }
                        }
                        else
                        {
                            if(entries[k].name.compare(nr_values_relax[i][k].name)==0 && nr_values_selection>0)//nr of values found
                            {
                                mean_sum+=entries[k].mean*nr_values_relax[i][k].nr[nr_values_selection-1];
                                stdabw_sum+=entries[k].stdabw*nr_values_relax[i][k].nr[nr_values_selection-1];
                                nr_sum+=nr_values_relax[i][k].nr[nr_values_selection-1];
                            }
                            else
                            {
                                mean_sum+=entries[k].mean;
                                stdabw_sum+=entries[k].stdabw;
                                nr_sum++;
                            }
                        }
                    }
                }

                if (iter==0)
                {
                    parameter_depths.push_back(depth[i]);
                    parameter_values.push_back(mean_sum/(nr_sum*scaling));
                    parameter_errors.push_back(stdabw_sum/(nr_sum*scaling));
                }
                else
                {
                    parameter_depths_relax.push_back(depth[i]);
                    parameter_values_relax.push_back(mean_sum/(nr_sum*scaling));
                    parameter_errors_relax.push_back(stdabw_sum/(nr_sum*scaling));
                }

                parameter_file.close();
                temp_parameter_file.close();
            }
        }

    }

    std::cout<<overview_name.c_str()<<": "<<parameter_values.size()<<std::endl;

    std::string filepath_parameter_out=path_results;
    filepath_parameter_out.append(filename.c_str());
    filepath_parameter_out.append(str.c_str());
    filepath_parameter_out.resize(filepath_parameter_out.size()-3);
    filepath_parameter_out.append("svg");

    std::string filepath_parameter_combined_out=path_results;
    filepath_parameter_combined_out.append(filename.c_str());
    filepath_parameter_combined_out.append("_combined");
    filepath_parameter_combined_out.append(str.c_str());
    filepath_parameter_combined_out.resize(filepath_parameter_combined_out.size()-3);
    filepath_parameter_combined_out.append("svg");

    if(minimal_bubble_distance==0 && parameter_values.size()>0) plot.draw_depth_errors(x.c_str(), y.c_str(), title.c_str(),
                                                             parameter_depths, parameter_values, parameter_errors, parameter_errors,
                                                             y_minimal, filepath_parameter_out.c_str());
    else if(parameter_values_relax.size()>0) plot.draw_depth_errors(x.c_str(), y.c_str(), title.c_str(), parameter_depths_relax,
                                                             parameter_values_relax, parameter_errors_relax, parameter_errors_relax,
                                                             y_minimal, filepath_parameter_out.c_str());
    if(parameter_values.size()>0 && parameter_values_relax.size()>0) plot.draw_depth_errors(x.c_str(), y.c_str(), title.c_str(),
                                                             parameter_depths, parameter_values, parameter_errors, parameter_errors,
                                                             parameter_depths_relax, parameter_values_relax, parameter_errors_relax,
                                                             parameter_errors_relax, y_minimal,
                                                             filepath_parameter_combined_out.c_str());

    delete depth;
};

void new_lin_depth_profile(std::string path_results, int* bagnr, int nr_depths, int minimal_bubble_distance, int minimal_grain_size,
    int grain_size_min, std::string str, plplot plot, std::string filename, std::string overview_name, std::string x, std::string y,
    std::string title, float y_minimal, int nr_values_selection, std::vector<std::vector<depth_numbers> > nr_values,
    std::vector<std::vector<depth_numbers> > nr_values_relax, int depth_max, float depth_bin_width, float scaling=1.0f)
{
    int * depth = new int[(size_t)nr_depths];

    std::vector<float> parameter_depths((size_t)(1+depth_max/depth_bin_width));
    std::vector<float> parameter_values((size_t)(1+depth_max/depth_bin_width),0.0f);
    std::vector<float> parameter_sum((size_t)(1+depth_max/depth_bin_width),0.0f);

    std::vector<float> parameter_depths_relax((size_t)(1+depth_max/depth_bin_width));
    std::vector<float> parameter_values_relax((size_t)(1+depth_max/depth_bin_width),0.0f);
    std::vector<float> parameter_sum_relax((size_t)(1+depth_max/depth_bin_width),0.0f);

    for(int i=0; i<nr_depths; i++)
    {
        int iterations=1;
        if (minimal_bubble_distance>0 && i<67) iterations=2;

        depth[i]=0.55*bagnr[i];

        std::stringstream label;
        label << bagnr[i] << "/";

        //string is read from temp file to check whether file is empty
        std::string teststring;

        for(int iter=0; iter<iterations; iter++)
        {
            std::stringstream s;

            if (iter==0) s << "." << minimal_grain_size;
            else s << "." << minimal_bubble_distance << "." << minimal_grain_size;
            if (grain_size_min>10) s << "_" << grain_size_min;
            s << ".txt";

            std::string filepath_parameter=path_results;
            filepath_parameter.append(label.str());
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

                    if (iter==0)
                    {
                        //nr of values found
                        if(!other_entry && entries[k].name.compare(nr_values[i][k].name)==0 && nr_values_selection>0)
                        {
                            parameter_values[depth[i]/depth_bin_width]+=entries[k].mean*nr_values[i][k].nr[nr_values_selection-1];
                            parameter_sum[depth[i]/depth_bin_width]+=nr_values[i][k].nr[nr_values_selection-1];
                        }
                    }
                    else
                    {
                        //nr of values found
                        if(!other_entry && entries[k].name.compare(nr_values_relax[i][k].name)==0 && nr_values_selection>0)
                        {
                            parameter_values_relax[depth[i]/depth_bin_width]+=
                                entries[k].mean*nr_values_relax[i][k].nr[nr_values_selection-1];
                            parameter_sum_relax[depth[i]/depth_bin_width]+=nr_values_relax[i][k].nr[nr_values_selection-1];
                        }
                    }
                }

                parameter_file.close();
                temp_parameter_file.close();
            }
        }
    }

    for (int bin=0; bin<parameter_values.size(); bin++)
    {
        parameter_depths[bin]=(bin+0.5f)*depth_bin_width;
    }

    for (int bin=0; bin<parameter_values_relax.size(); bin++)
    {
        parameter_depths_relax[bin]=(bin+0.5f)*depth_bin_width;
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

    for (int bin=0; bin<parameter_values_relax.size(); bin++)
    {
        if(parameter_sum_relax[bin]>0.0f) parameter_values_relax[bin]/=(parameter_sum_relax[bin]*scaling);
        else
        {
            parameter_depths_relax.erase (parameter_depths_relax.begin()+bin);
            parameter_values_relax.erase (parameter_values_relax.begin()+bin);
            parameter_sum_relax.erase (parameter_sum_relax.begin()+bin);
            bin--;
        }
    }

    std::string filepath_parameter_out=path_results;
    filepath_parameter_out.append("new_");
    filepath_parameter_out.append(filename.c_str());
    filepath_parameter_out.append(str.c_str());
    filepath_parameter_out.resize(filepath_parameter_out.size()-3);
    filepath_parameter_out.append("svg");

    std::string filepath_parameter_combined_out=path_results;
    filepath_parameter_combined_out.append("new_");
    filepath_parameter_combined_out.append(filename.c_str());
    filepath_parameter_combined_out.append("_combined");
    filepath_parameter_combined_out.append(str.c_str());
    filepath_parameter_combined_out.resize(filepath_parameter_combined_out.size()-3);
    filepath_parameter_combined_out.append("svg");

    if(minimal_bubble_distance==0 && parameter_values.size()>0) plot.draw_depth(x.c_str(), y.c_str(), title.c_str(), parameter_depths,
                                                             parameter_values, y_minimal, filepath_parameter_out.c_str());
    else if(parameter_values_relax.size()>0) plot.draw_depth(x.c_str(), y.c_str(), title.c_str(), parameter_depths_relax,
                                                             parameter_values_relax, y_minimal, filepath_parameter_out.c_str());
    if(parameter_values.size()>0 && parameter_values_relax.size()>0) plot.draw_depth(x.c_str(), y.c_str(), title.c_str(),
                                                             parameter_depths, parameter_values,parameter_depths_relax,
                                                             parameter_values_relax, y_minimal,
                                                             filepath_parameter_combined_out.c_str());

    delete depth;
};

void depth_profile(std::string path_results, ParameterFile paramFile, float depth_bin_width=0.0f)
{
    int bagnr[153] = {186, 226, 266, 304, 306, 346, 386, 426, 466, 506, 508, 521, 546, 586, 606, 656, 706, 756, 806, 846, 906, 956,
                      1006, 1056, 1105, 1147, 1187, 1227, 1306, 1346, 1375, 1426, 1496, 1546, 1616, 1676, 1766, 1815, 1866, 1885,
                      1896, 1906, 1946, 1997, 2007, 2028, 2046, 2067, 2068, 2069, 2085, 2106, 2146, 2156, 2204, 2216, 2227, 2235,
                      2246, 2275, 2283, 2295, 2327, 2336, 2386, 2466, 2506, 2546, 2579, 2580, 2586, 2606, 2646, 2686, 2726, 2766,
                      2806, 2846, 2926, 3196, 3206, 3214, 3216, 3226, 3237, 3245, 3246, 3249, 3256, 3266, 3276, 3278, 3286, 3293,
                      3296, 3306, 3307, 3313, 3316, 3326, 3336, 3346, 3356, 3396, 3436, 3451, 3476, 3516, 3556, 3596, 3636, 3676,
                      3726, 3756, 3796, 3817, 3836, 3849, 3856, 3876, 3906, 3917, 3956, 4006, 4016, 4026, 4077, 4106, 4116, 4126,
                      4156, 4196, 4236, 4276, 4297, 4299, 4300, 4306, 4316, 4356, 4397, 4441, 4442, 4476, 4514, 4515, 4516, 4517,
                      4518, 4519, 4520, 4556, 4606};

   	Parameter<int> nr_depths;
   	nr_depths.assign("", "nr_depths", 125);
    nr_depths.load(paramFile,"config");

    int * depth = new int[(size_t)nr_depths];

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

        //load nr values 
        std::vector<std::vector<depth_numbers> > nr_values;
        std::vector<std::vector<depth_numbers> > nr_values_relax;
        std::vector<std::vector<depth_numbers> > nr_values_new;
        std::vector<std::vector<depth_numbers> > nr_values_new_relax;

        for(int i=0; i<nr_depths; i++)
        {
            int iterations=1;
            if (minimal_bubble_distance>0 && i<67) iterations=2;

            std::stringstream label;
            label << bagnr[i] << "/";

            //string is read from temp file to check whether file is empty
            std::string teststring;

            for(int iter=0; iter<iterations; iter++)
            {
                std::stringstream s;

                if (iter==0) s << "." << minimal_grain_size;
                else s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                if (grain_size_min>10) s << "_" << grain_size_min;
                s << ".txt";

                std::string filepath_parameter=path_results;
                filepath_parameter.append(label.str());
                filepath_parameter.append("nr_values");
                filepath_parameter.append(s.str());

                std::ifstream parameter_file(filepath_parameter.c_str());
                std::ifstream temp_parameter_file(filepath_parameter.c_str());
                temp_parameter_file>>teststring;

                std::vector<depth_numbers> entries;

                if(parameter_file && teststring.size()!=0)
                {
                    while(!parameter_file.eof())
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

                        parameter_file>>name>>nr_grains>>nr_not_deg_ellipse>>nr_correct_flatt>>nr_correct_long>>nr_correct_length>>nr_correct_longest>>
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
                        if(nr_grains!=-1) entries.push_back(entry);
                    }

                    parameter_file.close();
                    temp_parameter_file.close();
                }

                if (iter==0) nr_values.push_back(entries);
                else nr_values_relax.push_back(entries);

                std::string filepath_parameter_new=path_results;
                filepath_parameter_new.append(label.str());
                filepath_parameter_new.append("nr_values2");
                filepath_parameter_new.append(s.str());

                std::ifstream parameter_new_file(filepath_parameter_new.c_str());
                std::ifstream temp_parameter_new_file(filepath_parameter_new.c_str());
                temp_parameter_new_file>>teststring;

                entries.clear();

                if(parameter_new_file && teststring.size()!=0)
                {
                    while(!parameter_new_file.eof())
                    {
                        std::string name;
                        int nr_correct_length=-1;
                        int nr_curved_pixels;

                        parameter_new_file>>name>>nr_correct_length>>nr_curved_pixels;
                        depth_numbers entry;
                        entry.name=name;
                        entry.nr[4]=nr_correct_length;
                        entry.nr[7]=nr_curved_pixels;
                        if(nr_correct_length!=-1) entries.push_back(entry);
                    }

                    parameter_new_file.close();
                    temp_parameter_new_file.close();
                }

                if (iter==0) nr_values_new.push_back(entries);
                else nr_values_new_relax.push_back(entries);
            }
        }

        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "grainshape", "Grain shape", "Depth in m", "Roundness factor (4pi*Area/Perimeter^2)", "Grain roundness profile",
            0.35f, 1, nr_values, nr_values_relax);
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "nr_grainarcs", "Number of grain boundaries", "Depth in m", "Nr of grain boundaries",
            "Number of grain boundaries profile", 0.0f, 0, nr_values, nr_values_relax);
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "nr_neighbors", "Number of grain neighbors", "Depth in m", "Nr of grain neighbors",
            "Number of grain neighbors profile", 0.0f, 0, nr_values, nr_values_relax);
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "length_grainarcs", "Grain arc length", "Depth in m", length_unit, "Grain boundary length profile", 0.0f, 5,
            nr_values, nr_values_relax, length_scaling());
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "length_grainarcs2", "Grain arc length", "Depth in m", length_unit, "Grain boundary length profile", 0.0f, 5,
            nr_values_new, nr_values_new_relax, length_scaling());
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "longest_grainarcs", "Grain longest arc length", "Depth in m", length_unit, "Grain longest boundary length profile",
            0.0f, 6, nr_values, nr_values_relax, length_scaling());
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "dislocation_densities", "Dislocation density differences", "Depth in m",
            "Dislocation density difference in m^(-2)", "Profile of mean dislocation density difference at grain boundaries", 0.0f,
            8, nr_values, nr_values_relax);
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "dislocation_densities2", "Dislocation density differences", "Depth in m",
            "Dislocation density difference in m^(-2)", "Profile of mean dislocation density difference at grain boundaries", 0.0f,
            8, nr_values_new, nr_values_new_relax);
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "grain_equivradius", "Equivalent grain radius", "Depth in m", length_unit, "Grain equivalent radius profile", 0.0f,
            1, nr_values, nr_values_relax, length_scaling());
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "grain_boxshape", "Grain box shape", "Depth in m", "Box flattening factor", "Box flattening profile", 0.9f, 2,
            nr_values, nr_values_relax);
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "grain_boxwidth", "Grain box width", "Depth in m", length_unit, "Grain box width profile", 0.0f, 2, nr_values,
            nr_values_relax, length_scaling());
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "grain_boxheight", "Grain box height", "Depth in m", length_unit, "Grain box height profile", 0.0f, 2, nr_values,
            nr_values_relax, length_scaling());
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "grain_ellipselong", "Grain ellipse long axis", "Depth in m", length_unit, "Ellipse long axis profile", 0.0f, 4,
            nr_values, nr_values_relax, length_scaling());
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "grain_ellipseshape", "Grain ellipse shape", "Depth in m", "Ellipse flattening factor",
            "Ellipse flattening profile", 0.9f, 3, nr_values, nr_values_relax);
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "grain_ellipseangle", "Grain ellipse orientation", "Depth in m", "Ellipse long axis angle",
            "Ellipse long axis angle profile", 0.0f, 2, nr_values, nr_values_relax);
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "boundary_orientation", "Boundary orientation", "Depth in m", "Boundary orientation",
            "Boundary orientation profile", 0.0f, 7, nr_values, nr_values_relax);
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "grainwidth", "Grain width", "Depth in m", length_unit, "Grain width profile", 0.0f, 1, nr_values, nr_values_relax,
            length_scaling());
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "grainheight", "Grain height", "Depth in m", length_unit, "Grain height profile", 0.0f, 1, nr_values,
            nr_values_relax, length_scaling());
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "grainflattening", "Vertical grain flattening", "Depth in m", "Vertical grain flattening factor",
            "Vertical grain flattening profile", 0.9f, 1, nr_values, nr_values_relax);
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "turning_points", "Number of turning points", "Depth in m", "Nr of turning points",
            "Grain boundary turning points profile", 0.0f, 9, nr_values, nr_values_relax);
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "grain_perimeter_ratio", "Grain perimeter ratio", "Depth in m", "Perimeter ratio", "Perimeter ratio profile", 0.7f,
            1, nr_values, nr_values_relax);
        lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min, s.str(),
            plot, "subGB_densities", "Subgrain boundary density", "Depth in m", "Subgrain boundary density (1/mm)",
            "Subgrain boundary density profile", 0.0f, 0, nr_values, nr_values_relax, 1.0f/length_scaling());

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

        std::vector<float> grainsize_all_depths_relax;
        std::vector<float> grainsize_all_values_relax;
        std::vector<float> grainsize_all_errors_relax;

        std::vector<float> grainsize_fit_depths_relax;
        std::vector<float> grainsize_fit_values_relax;
        std::vector<float> grainsize_fit_errors_low_relax;
        std::vector<float> grainsize_fit_errors_high_relax;

        std::vector<float> corr_ellipse_box_flattening_depths_relax;
        std::vector<float> corr_ellipse_box_flattening_values_relax;

        std::vector< std::vector<float> > grainsize_step_depths_relax;
        std::vector< std::vector<float> > grainsize_step_values_relax;
        std::vector< std::vector<float> > grainsize_step_errors_relax;

        std::vector< std::vector<float> > grainsize_percent_depths_relax(9);
        std::vector< std::vector<float> > grainsize_percent_values_relax(9);
        std::vector< std::vector<float> > grainsize_percent_errors_relax(9);

        std::vector< std::vector<float> > grainsize_quantile_depths_relax(10);
        std::vector< std::vector<float> > grainsize_quantile_values_relax(10);
        std::vector< std::vector<float> > grainsize_quantile_errors_relax(10);

        std::vector< std::vector<float> > grainsize_bin_depths_relax;
        std::vector< std::vector<float> > grainsize_bin_values_relax;

        std::vector< std::vector<float> > grainradius_bin_depths_relax;
        std::vector< std::vector<float> > grainradius_bin_values_relax;

        std::vector< std::vector<float> > grainradius_norm_bin_depths_relax;
        std::vector< std::vector<float> > grainradius_norm_bin_values_relax;

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

        std::vector<float> grain_perimeter_ratio_depths;
        std::vector<float> grain_perimeter_ratio_values;
        std::vector<float> grain_perimeter_ratio_errors;

        std::vector<float> grain_perimeter_ratio_depths_bubble;
        std::vector<float> grain_perimeter_ratio_values_bubble;
        std::vector<float> grain_perimeter_ratio_errors_bubble;

        for(int i=0; i<nr_depths; i++)
        {
            int iterations=1;
            if (minimal_bubble_distance>0 && i<67) iterations=2;

            depth[i]=0.55*bagnr[i];

            std::stringstream label;
            label << bagnr[i] << "/";

            //string is read from temp file to check whether file is empty
            std::string teststring;

            for(int iter=0; iter<iterations; iter++)
            {
                std::stringstream s;

                if (iter==0) s << "." << minimal_grain_size;
                else s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                if (grain_size_min>10) s << "_" << grain_size_min;
                s << ".txt";

                std::string filepath_grainsize_all=path_results;
                filepath_grainsize_all.append(label.str());
                filepath_grainsize_all.append("grainsize_all");
                filepath_grainsize_all.append(s.str());

                std::string filepath_grainsize_fit=path_results;
                filepath_grainsize_fit.append(label.str());
                filepath_grainsize_fit.append("grainsize_fit");
                filepath_grainsize_fit.append(s.str());

                std::string filepath_corr_ellipse_box_flattening=path_results;
                filepath_corr_ellipse_box_flattening.append(label.str());
                filepath_corr_ellipse_box_flattening.append("corr_ellipse_box_flattening");
                filepath_corr_ellipse_box_flattening.append(s.str());

                std::string filepath_grainsize_step=path_results;
                filepath_grainsize_step.append(label.str());
                filepath_grainsize_step.append("grainsize_step");
                filepath_grainsize_step.append(s.str());

                std::string filepath_grainsize_percent=path_results;
                filepath_grainsize_percent.append(label.str());
                filepath_grainsize_percent.append("grainsize_percent");
                filepath_grainsize_percent.append(s.str());

                std::string filepath_grainsize_quantile=path_results;
                filepath_grainsize_quantile.append(label.str());
                filepath_grainsize_quantile.append("grainsize_quantile");
                filepath_grainsize_quantile.append(s.str());

                std::string filepath_grainsize_bin=path_results;
                filepath_grainsize_bin.append(label.str());
                filepath_grainsize_bin.append("grainsize_bin");
                filepath_grainsize_bin.append(s.str());

                std::string filepath_grainradius_bin=path_results;
                filepath_grainradius_bin.append(label.str());
                filepath_grainradius_bin.append("grainradius_bin");
                filepath_grainradius_bin.append(s.str());

                std::string filepath_grainradius_norm_bin=path_results;
                filepath_grainradius_norm_bin.append(label.str());
                filepath_grainradius_norm_bin.append("grainradius_norm_bin");
                filepath_grainradius_norm_bin.append(s.str());

                std::string filepath_grain_neighbors_bin=path_results;
                filepath_grain_neighbors_bin.append(label.str());
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

                    float mean_sum=0;
                    float stdabw_sum=0;
                    int nr_sum=0;

                    //check for multiple entries for same image
                    for (int k=0; k<entries.size(); k++)
                    {
                        bool other_entry=false;

                        for(int j=k+1; j<entries.size() && !other_entry; j++)
                        {
                            if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                        }

                        if(!other_entry)
                        {
                            if (iter==0)
                            {
                                if(entries[k].name.compare(nr_values[i][k].name)==0)//nr of values found
                                {
                                    mean_sum+=entries[k].mean*nr_values[i][k].nr[0];
                                    stdabw_sum+=entries[k].stdabw*nr_values[i][k].nr[0];
                                    nr_sum+=nr_values[i][k].nr[0];
                                }
                                else
                                {
                                    mean_sum+=entries[k].mean;
                                    stdabw_sum+=entries[k].stdabw;
                                    nr_sum++;
                                }
                            }
                            else
                            {
                                if(entries[k].name.compare(nr_values_relax[i][k].name)==0)//nr of values found
                                {
                                    mean_sum+=entries[k].mean*nr_values_relax[i][k].nr[0];
                                    stdabw_sum+=entries[k].stdabw*nr_values_relax[i][k].nr[0];
                                    nr_sum+=nr_values_relax[i][k].nr[0];
                                }
                                else
                                {
                                    mean_sum+=entries[k].mean;
                                    stdabw_sum+=entries[k].stdabw;
                                    nr_sum++;
                                }
                            }
                        }
                    }

                    if (iter==0)
                    {
                        grainsize_all_depths.push_back(depth[i]);
                        grainsize_all_values.push_back(mean_sum/(nr_sum*area_scaling()));
                        grainsize_all_errors.push_back(stdabw_sum/(nr_sum*area_scaling()));
                    }
                    else
                    {
                        grainsize_all_depths_relax.push_back(depth[i]);
                        grainsize_all_values_relax.push_back(mean_sum/(nr_sum*area_scaling()));
                        grainsize_all_errors_relax.push_back(stdabw_sum/(nr_sum*area_scaling()));
                    }

                    grainsize_all_file.close();
                    temp_grainsize_all_file.close();

                    if (iter==iterations-1)
                    {
                        grainsize_combined_depths[0].push_back(depth[i]);
                        grainsize_combined_values[0].push_back(mean_sum/(nr_sum*area_scaling()));
                    }
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

                    float mean_sum=0;
                    float stdabw_low_sum=0;
                    float stdabw_high_sum=0;
                    int nr_sum=0;

                    //check for multiple entries for same image
                    for (int k=0; k<entries.size(); k++)
                    {
                        bool other_entry=false;

                        for(int j=k+1; j<entries.size() && !other_entry; j++)
                        {
                            if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                        }

                        if(!other_entry)
                        {
                            mean_sum+=entries[k].mean;
                            stdabw_low_sum+=entries[k].stdabw;
                            stdabw_high_sum+=entries[k].stdabw2;
                            nr_sum++;
                        }
                    }

                    if (iter==0)
                    {
                        grainsize_fit_depths.push_back(depth[i]);
                        grainsize_fit_values.push_back(mean_sum/(nr_sum*area_scaling()));
                        grainsize_fit_errors_low.push_back(stdabw_low_sum/(nr_sum*area_scaling()));
                        grainsize_fit_errors_high.push_back(std::min(1000.0f,stdabw_high_sum/(nr_sum*area_scaling())));
                    }
                    else
                    {
                        grainsize_fit_depths_relax.push_back(depth[i]);
                        grainsize_fit_values_relax.push_back(mean_sum/(nr_sum*area_scaling()));
                        grainsize_fit_errors_low_relax.push_back(stdabw_low_sum/(nr_sum*area_scaling()));
                        grainsize_fit_errors_high_relax.push_back(std::min(1000.0f,stdabw_high_sum/(nr_sum*area_scaling())));
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

                    float mean_sum=0;
                    int nr_sum=0;

                    //check for multiple entries for same image
                    for (int k=0; k<entries.size(); k++)
                    {
                        bool other_entry=false;

                        for(int j=k+1; j<entries.size() && !other_entry; j++)
                        {
                            if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                        }

                        if(!other_entry)
                        {
                            mean_sum+=entries[k].mean;
                            nr_sum++;
                        }
                    }

                    if (iter==0)
                    {
                        corr_ellipse_box_flattening_depths.push_back(depth[i]);
                        corr_ellipse_box_flattening_values.push_back(mean_sum/nr_sum);
                    }
                    else
                    {
                        corr_ellipse_box_flattening_depths_relax.push_back(depth[i]);
                        corr_ellipse_box_flattening_values_relax.push_back(mean_sum/nr_sum);
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

                    std::vector<float> mean_sum(min_step_nr,0.0f);
                    std::vector<float> stdabw_sum(min_step_nr,0.0f);
                    int nr_sum=0;

                    //check for multiple entries for same image
                    for (int k=0; k<entries.size(); k++)
                    {
                        bool other_entry=false;

                        for(int j=k+1; j<entries.size() && !other_entry; j++)
                        {
                            if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                        }

                        if(!other_entry)
                        {
                            for(int s=0; s<min_step_nr; s++)
                            {
                                mean_sum[s]+=entries[k].mean[s];
                                stdabw_sum[s]+=entries[k].stdabw[s];
                            } 
                            nr_sum++;
                        }
                    }

                    if (nr_sum>0)
                    {
                        if (iter==iterations-1)
                        {
                            grainsize_step_profile_depths.push_back(depth[i]);

                            std::vector<float> new_entry;
                            grainsize_step_profile_values.push_back(new_entry);
                        }

                        for(int s=0; s<min_step_nr; s++)
                        {
                            if (iter==0)
                            {
                                if(s+1>grainsize_step_depths.size()) grainsize_step_depths.resize(s+1);
                                if(s+1>grainsize_step_values.size()) grainsize_step_values.resize(s+1);
                                if(s+1>grainsize_step_errors.size()) grainsize_step_errors.resize(s+1);
                                grainsize_step_depths[s].push_back(depth[i]);
                                grainsize_step_values[s].push_back(mean_sum[s]/(nr_sum*area_scaling()));
                                grainsize_step_errors[s].push_back(stdabw_sum[s]/(nr_sum*area_scaling()));
                            }
                            else
                            {
                                if(s+1>grainsize_step_depths_relax.size()) grainsize_step_depths_relax.resize(s+1);
                                if(s+1>grainsize_step_values_relax.size()) grainsize_step_values_relax.resize(s+1);
                                if(s+1>grainsize_step_errors_relax.size()) grainsize_step_errors_relax.resize(s+1);
                                grainsize_step_depths_relax[s].push_back(depth[i]);
                                grainsize_step_values_relax[s].push_back(mean_sum[s]/(nr_sum*area_scaling()));
                                grainsize_step_errors_relax[s].push_back(stdabw_sum[s]/(nr_sum*area_scaling()));
                            }
                            if (iter==iterations-1)
                            {
                                grainsize_step_profile_values.back().push_back(mean_sum[s]/(nr_sum*area_scaling()));
                            }
                        }

                        if (iter==iterations-1)
                        {
                            grainsize_combined_depths[1].push_back(depth[i]);
                            grainsize_combined_values[1].push_back(mean_sum[1]/(nr_sum*area_scaling()));
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
                        if(entry.mean.size()==0) break;
			        }

                    std::vector<float> mean_sum(10,0.0f);
                    std::vector<float> stdabw_sum(10,0.0f);
                    int nr_sum=0;

                    //check for multiple entries for same image
                    for (int k=0; k<entries.size(); k++)
                    {
                        bool other_entry=false;

                        for(int j=k+1; j<entries.size() && !other_entry; j++)
                        {
                            if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                        }

                        if(!other_entry)
                        {
                            for(int p=0; p<10; p++)
                            {
                                mean_sum[p]+=entries[k].mean[p];
                                stdabw_sum[p]+=entries[k].stdabw[p];
                            } 
                            nr_sum++;
                        }
                    }

                    if (iter==iterations-1)
                    {
                        grainsize_percent_profile_depths.push_back(depth[i]);

                        std::vector<float> new_entry;
                        grainsize_percent_profile_values.push_back(new_entry);
                    }

                    for(int p=0; p<9; p++)
                    {
                        if (iter==0)
                        {
                            grainsize_percent_depths[p].push_back(depth[i]);
                            grainsize_percent_values[p].push_back(mean_sum[p]/(nr_sum*area_scaling()));
                            grainsize_percent_errors[p].push_back(stdabw_sum[p]/(nr_sum*area_scaling()));
                        }
                        else
                        {
                            grainsize_percent_depths_relax[p].push_back(depth[i]);
                            grainsize_percent_values_relax[p].push_back(mean_sum[p]/(nr_sum*area_scaling()));
                            grainsize_percent_errors_relax[p].push_back(stdabw_sum[p]/(nr_sum*area_scaling()));
                        }
                        if (iter==iterations-1)
                        {
                            grainsize_percent_profile_values.back().push_back(mean_sum[p]/(nr_sum*area_scaling()));
                        }
                    }

                    grainsize_percent_file.close();
                    temp_grainsize_percent_file.close();

                    if (iter==iterations-1)
                    {
                        grainsize_percent_profile_values.back().push_back(mean_sum[9]/(nr_sum*area_scaling()));
                        grainsize_combined_depths[2].push_back(depth[i]);
                        grainsize_combined_values[2].push_back(mean_sum[2]/(nr_sum*area_scaling()));
                    }
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
                        if(entry.mean.size()==0) break;
			        }

                    std::vector<float> mean_sum(10,0.0f);
                    std::vector<float> stdabw_sum(10,0.0f);
                    int nr_sum=0;

                    //check for multiple entries for same image
                    for (int k=0; k<entries.size(); k++)
                    {
                        bool other_entry=false;

                        for(int j=k+1; j<entries.size() && !other_entry; j++)
                        {
                            if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                        }

                        if(!other_entry)
                        {
                            for(int q=0; q<10; q++)
                            {
                                mean_sum[q]+=entries[k].mean[q];
                                stdabw_sum[q]+=entries[k].stdabw[q];
                            } 
                            nr_sum++;
                        }
                    }

                    if (iter==iterations-1)
                    {
                        grainsize_quantile_profile_depths.push_back(depth[i]);

                        std::vector<float> new_entry;
                        grainsize_quantile_profile_values.push_back(new_entry);
                    }

                    for(int q=0; q<10; q++)
                    {
                        if (iter==0)
                        {
                            grainsize_quantile_depths[q].push_back(depth[i]);
                            grainsize_quantile_values[q].push_back(mean_sum[q]/(nr_sum*area_scaling()));
                            grainsize_quantile_errors[q].push_back(stdabw_sum[q]/(nr_sum*area_scaling()));
                        }
                        else
                        {
                            grainsize_quantile_depths_relax[q].push_back(depth[i]);
                            grainsize_quantile_values_relax[q].push_back(mean_sum[q]/(nr_sum*area_scaling()));
                            grainsize_quantile_errors_relax[q].push_back(stdabw_sum[q]/(nr_sum*area_scaling()));
                        }
                        if (iter==iterations-1)
                        {
                            grainsize_quantile_profile_values.back().push_back(mean_sum[q]/(nr_sum*area_scaling()));
                        }
                    }

                    grainsize_quantile_file.close();
                    temp_grainsize_quantile_file.close();

                    if (iter==iterations-1)
                    {
                        grainsize_combined_depths[3].push_back(depth[i]);
                        grainsize_combined_values[3].push_back(mean_sum[9]/(nr_sum*area_scaling()));
                    }
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

                    std::vector<float> nr_grains_bin(min_bin_nr,0.0f);//nr of grains for each grain size bin
                    int nr_grains_sum=0;

                    //check for multiple entries for same image
                    for (int k=0; k<entries.size(); k++)
                    {
                        bool other_entry=false;

                        for(int j=k+1; j<entries.size() && !other_entry; j++)
                        {
                            if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                        }

                        if(!other_entry)
                        {
                            for(int bin=0; bin<min_bin_nr; bin++)
                            {
                                nr_grains_bin[bin]+=entries[k].mean[bin];
                                nr_grains_sum+=entries[k].mean[bin];
                            } 
                        }
                    }

                    if (iter==iterations-1)
                    {
                        grainsize_bin_profile_depths.push_back(depth[i]);

                        std::vector<float> new_entry;
                        grainsize_bin_profile_values.push_back(new_entry);
                    }

                    for(int bin=0; bin<min_bin_nr; bin++)
                    {
                        if (iter==0)
                        {
                            if(bin+1>grainsize_bin_depths.size()) grainsize_bin_depths.resize(bin+1);
                            if(bin+1>grainsize_bin_values.size()) grainsize_bin_values.resize(bin+1);
                            grainsize_bin_depths[bin].push_back(depth[i]);
                            grainsize_bin_values[bin].push_back(nr_grains_bin[bin]/nr_grains_sum);
                        }
                        else
                        {
                            if(bin+1>grainsize_bin_depths_relax.size()) grainsize_bin_depths_relax.resize(bin+1);
                            if(bin+1>grainsize_bin_values_relax.size()) grainsize_bin_values_relax.resize(bin+1);
                            grainsize_bin_depths_relax[bin].push_back(depth[i]);
                            grainsize_bin_values_relax[bin].push_back(nr_grains_bin[bin]/nr_grains_sum);
                        }
                        if (iter==iterations-1)
                        {
                            grainsize_bin_profile_values.back().push_back(nr_grains_bin[bin]/nr_grains_sum);
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

                    std::vector<float> nr_grains_bin(min_bin_nr,0.0f);//nr of grains for each grain radius bin
                    int nr_grains_sum=0;

                    //check for multiple entries for same image
                    for (int k=0; k<entries.size(); k++)
                    {
                        bool other_entry=false;

                        for(int j=k+1; j<entries.size() && !other_entry; j++)
                        {
                            if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                        }

                        if(!other_entry)
                        {
                            for(int bin=0; bin<min_bin_nr; bin++)
                            {
                                nr_grains_bin[bin]+=entries[k].mean[bin];
                                nr_grains_sum+=entries[k].mean[bin];
                            } 
                        }
                    }

                    if (iter==iterations-1)
                    {
                        grainradius_bin_profile_depths.push_back(depth[i]);

                        std::vector<float> new_entry;
                        grainradius_bin_profile_values.push_back(new_entry);
                    }

                    for(int bin=0; bin<min_bin_nr; bin++)
                    {
                        if (iter==0)
                        {
                            if(bin+1>grainradius_bin_depths.size()) grainradius_bin_depths.resize(bin+1);
                            if(bin+1>grainradius_bin_values.size()) grainradius_bin_values.resize(bin+1);
                            grainradius_bin_depths[bin].push_back(depth[i]);
                            grainradius_bin_values[bin].push_back(nr_grains_bin[bin]/nr_grains_sum);
                        }
                        else
                        {
                            if(bin+1>grainradius_bin_depths_relax.size()) grainradius_bin_depths_relax.resize(bin+1);
                            if(bin+1>grainradius_bin_values_relax.size()) grainradius_bin_values_relax.resize(bin+1);
                            grainradius_bin_depths_relax[bin].push_back(depth[i]);
                            grainradius_bin_values_relax[bin].push_back(nr_grains_bin[bin]/nr_grains_sum);
                        }
                        if (iter==iterations-1)
                        {
                            grainradius_bin_profile_values.back().push_back(nr_grains_bin[bin]/nr_grains_sum);
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

                    std::vector<float> nr_grains_bin(min_bin_nr,0.0f);//nr of grains for each grain radius bin
                    int nr_grains_sum=0;

                    //check for multiple entries for same image
                    for (int k=0; k<entries.size(); k++)
                    {
                        bool other_entry=false;

                        for(int j=k+1; j<entries.size() && !other_entry; j++)
                        {
                            if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                        }

                        if(!other_entry)
                        {
                            for(int bin=0; bin<min_bin_nr; bin++)
                            {
                                nr_grains_bin[bin]+=entries[k].mean[bin];
                                nr_grains_sum+=entries[k].mean[bin];
                            } 
                        }
                    }

                    if (iter==iterations-1)
                    {
                        grainradius_norm_bin_profile_depths.push_back(depth[i]);

                        std::vector<float> new_entry;
                        grainradius_norm_bin_profile_values.push_back(new_entry);
                    }

                    for(int bin=0; bin<min_bin_nr; bin++)
                    {
                        if (iter==0)
                        {
                            if(bin+1>grainradius_norm_bin_depths.size()) grainradius_norm_bin_depths.resize(bin+1);
                            if(bin+1>grainradius_norm_bin_values.size()) grainradius_norm_bin_values.resize(bin+1);
                            grainradius_norm_bin_depths[bin].push_back(depth[i]);
                            grainradius_norm_bin_values[bin].push_back(nr_grains_bin[bin]/nr_grains_sum);
                        }
                        else
                        {
                            if(bin+1>grainradius_norm_bin_depths_relax.size()) grainradius_norm_bin_depths_relax.resize(bin+1);
                            if(bin+1>grainradius_norm_bin_values_relax.size()) grainradius_norm_bin_values_relax.resize(bin+1);
                            grainradius_norm_bin_depths_relax[bin].push_back(depth[i]);
                            grainradius_norm_bin_values_relax[bin].push_back(nr_grains_bin[bin]/nr_grains_sum);
                        }
                        if (iter==iterations-1)
                        {
                            grainradius_norm_bin_profile_values.back().push_back(nr_grains_bin[bin]/nr_grains_sum);
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

                    std::vector<float> nr_grains_bin(min_bin_nr,0.0f);//nr of grains for each grain radius bin
                    int nr_grains_sum=0;

                    //check for multiple entries for same image
                    for (int k=0; k<entries.size(); k++)
                    {
                        bool other_entry=false;

                        for(int j=k+1; j<entries.size() && !other_entry; j++)
                        {
                            if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                        }

                        if(!other_entry)
                        {
                            for(int bin=0; bin<min_bin_nr; bin++)
                            {
                                nr_grains_bin[bin]+=entries[k].mean[bin];
                                nr_grains_sum+=entries[k].mean[bin];
                            } 
                        }
                    }

                    if (iter==iterations-1)
                    {
                        grain_neighbors_bin_profile_depths.push_back(depth[i]);

                        std::vector<float> new_entry;
                        grain_neighbors_bin_profile_values.push_back(new_entry);
                    }

                    for(int bin=0; bin<min_bin_nr; bin++)
                    {
                        if (iter==iterations-1)
                        {
                            grain_neighbors_bin_profile_values.back().push_back(nr_grains_bin[bin]/nr_grains_sum);
                        }
                    }

                    grain_neighbors_bin_file.close();
                    temp_grain_neighbors_bin_file.close();
                }
            }

            std::string filepath_dihedral_angles=path_results;
            filepath_dihedral_angles.append(label.str());
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

                float mean_sum=0;
                float stdabw_sum=0;
                int nr_sum=0;

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    if(!other_entry)
                    {
                        if(iterations==1)
                        {
                            if(entries[k].name.compare(nr_values[i][k].name)==0)//nr of values found
                            {
                                mean_sum+=entries[k].mean*nr_values[i][k].nr[6];
                                stdabw_sum+=entries[k].stdabw*nr_values[i][k].nr[6];
                                nr_sum+=nr_values[i][k].nr[6];
                            }
                            else
                            {
                                mean_sum+=entries[k].mean;
                                stdabw_sum+=entries[k].stdabw;
                                nr_sum++;
                            }
                        }
                        else
                        {
                            if(entries[k].name.compare(nr_values_relax[i][k].name)==0)//nr of values found
                            {
                                mean_sum+=entries[k].mean*nr_values_relax[i][k].nr[6];
                                stdabw_sum+=entries[k].stdabw*nr_values_relax[i][k].nr[6];
                                nr_sum+=nr_values_relax[i][k].nr[6];
                            }
                            else
                            {
                                mean_sum+=entries[k].mean;
                                stdabw_sum+=entries[k].stdabw;
                                nr_sum++;
                            }
                        }
                    }
                }

                dihedral_angles_depths.push_back(depth[i]);
                dihedral_angles_values.push_back(mean_sum/nr_sum);
                dihedral_angles_errors.push_back(stdabw_sum/nr_sum);
                dihedral_angles_file.close();
                temp_dihedral_angles_file.close();
            }

            std::string filepath_dihedral_angles2=path_results;
            filepath_dihedral_angles2.append(label.str());
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

                float mean_sum=0;
                float stdabw_sum=0;
                int nr_sum=0;

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    if(!other_entry)
                    {
                        if(iterations==1)
                        {
                            if(entries[k].name.compare(nr_values[i][k].name)==0)//nr of values found
                            {
                                mean_sum+=entries[k].mean*nr_values[i][k].nr[6];
                                stdabw_sum+=entries[k].stdabw*nr_values[i][k].nr[6];
                                nr_sum+=nr_values[i][k].nr[6];
                            }
                            else
                            {
                                mean_sum+=entries[k].mean;
                                stdabw_sum+=entries[k].stdabw;
                                nr_sum++;
                            }
                        }
                        else
                        {
                            if(entries[k].name.compare(nr_values_relax[i][k].name)==0)//nr of values found
                            {
                                mean_sum+=entries[k].mean*nr_values_relax[i][k].nr[6];
                                stdabw_sum+=entries[k].stdabw*nr_values_relax[i][k].nr[6];
                                nr_sum+=nr_values_relax[i][k].nr[6];
                            }
                            else
                            {
                                mean_sum+=entries[k].mean;
                                stdabw_sum+=entries[k].stdabw;
                                nr_sum++;
                            }
                        }
                    }
                }

                dihedral_angles2_depths.push_back(depth[i]);
                dihedral_angles2_values.push_back(mean_sum/nr_sum);
                dihedral_angles2_errors.push_back(stdabw_sum/nr_sum);
                dihedral_angles2_file.close();
                temp_dihedral_angles2_file.close();
            }

            std::string filepath_disl_dens_longest_boundaries=path_results;
            filepath_disl_dens_longest_boundaries.append(label.str());
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

                std::vector<float> mean_sum(min_step_nr,0.0f);
                std::vector<float> stdabw_sum(min_step_nr,0.0f);
                int nr_sum=0;

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    if(!other_entry)
                    {
                        for(int s=0; s<min_step_nr; s++)
                        {
                            mean_sum[s]+=entries[k].mean[s];
                            stdabw_sum[s]+=entries[k].stdabw[s];
                        } 
                        nr_sum++;
                    }
                }

                disl_dens_longest_boundaries_profile_depths.push_back(depth[i]);

                std::vector<float> new_entry;
                disl_dens_longest_boundaries_profile_values.push_back(new_entry);

                for(int s=0; s<min_step_nr; s++)
                {
                    if(s+1>disl_dens_longest_boundaries_depths.size()) disl_dens_longest_boundaries_depths.resize(s+1);
                    if(s+1>disl_dens_longest_boundaries_values.size()) disl_dens_longest_boundaries_values.resize(s+1);
                    if(s+1>disl_dens_longest_boundaries_errors.size()) disl_dens_longest_boundaries_errors.resize(s+1);
                    disl_dens_longest_boundaries_depths[s].push_back(depth[i]);
                    disl_dens_longest_boundaries_values[s].push_back(mean_sum[s]/nr_sum);
                    disl_dens_longest_boundaries_errors[s].push_back(stdabw_sum[s]/nr_sum);

                    disl_dens_longest_boundaries_profile_values.back().push_back(mean_sum[s]/nr_sum);
                }

                disl_dens_longest_boundaries_file.close();
                temp_disl_dens_longest_boundaries_file.close();
            }

            std::string filepath_disl_dens_percent_boundaries=path_results;
            filepath_disl_dens_percent_boundaries.append(label.str());
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
                    if(entry.mean.size()==0) break;
		        }

                std::vector<float> mean_sum(10,0.0f);
                std::vector<float> stdabw_sum(10,0.0f);
                int nr_sum=0;

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    if(!other_entry)
                    {
                        for(int p=0; p<10; p++)
                        {
                            mean_sum[p]+=entries[k].mean[p];
                            stdabw_sum[p]+=entries[k].stdabw[p];
                        } 
                        nr_sum++;
                    }
                }

                disl_dens_percent_boundaries_profile_depths.push_back(depth[i]);

                std::vector<float> new_entry;
                disl_dens_percent_boundaries_profile_values.push_back(new_entry);

                for(int p=0; p<10; p++)
                {
                    disl_dens_percent_boundaries_depths[p].push_back(depth[i]);
                    disl_dens_percent_boundaries_values[p].push_back(mean_sum[p]/nr_sum);
                    disl_dens_percent_boundaries_errors[p].push_back(stdabw_sum[p]/nr_sum);

                    disl_dens_percent_boundaries_profile_values.back().push_back(mean_sum[p]/nr_sum);
                }

                disl_dens_percent_boundaries_file.close();
                temp_disl_dens_percent_boundaries_file.close();
            }

            std::string filepath_percentage_filled_grains=path_results;
            filepath_percentage_filled_grains.append(label.str());
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
                    if(entry.mean.size()==0) break;
		        }

                std::vector<float> mean_sum(13,0.0f);
                int nr_sum=0;

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    if(!other_entry)
                    {
                        for(int p=0; p<13; p++)
                        {
                            mean_sum[p]+=entries[k].mean[p];
                        } 
                        nr_sum++;
                    }
                }

                percentage_filled_grains_depths[0].push_back(depth[i]);
                percentage_filled_grains_values[0].push_back(mean_sum[1]/nr_sum);

                percentage_filled_grains_depths[1].push_back(depth[i]);
                percentage_filled_grains_values[1].push_back(mean_sum[2]/nr_sum);

                percentage_filled_grains_depths[2].push_back(depth[i]);
                percentage_filled_grains_values[2].push_back(mean_sum[7]/nr_sum);

                percentage_filled_grains_depths[3].push_back(depth[i]);
                percentage_filled_grains_values[3].push_back(mean_sum[12]/nr_sum);

                percentage_filled_grains_file.close();
                temp_percentage_filled_grains_file.close();
            }

            std::string filepath_grain_perimeter_ratio=path_results;
            filepath_grain_perimeter_ratio.append(label.str());
            filepath_grain_perimeter_ratio.append("grain_perimeter_ratio");
            filepath_grain_perimeter_ratio.append(s.str());

            //grain perimeter ratio
            std::ifstream grain_perimeter_ratio_file(filepath_grain_perimeter_ratio.c_str());
            std::ifstream temp_grain_perimeter_ratio_file(filepath_grain_perimeter_ratio.c_str());
            temp_grain_perimeter_ratio_file>>teststring;

            if(grain_perimeter_ratio_file && teststring.size()!=0)
            {
                std::vector<depth_parameter> entries;

	            while(!grain_perimeter_ratio_file.eof())
	            {
                    std::string name;
                    float mean=-1, stdabw;
		            grain_perimeter_ratio_file>>name>>mean>>stdabw;
                    depth_parameter entry;
                    entry.name=name;
                    entry.mean=mean;
                    entry.stdabw=stdabw;
                    if(mean!=-1) entries.push_back(entry);
	            }

                float mean_sum=0;
                float stdabw_sum=0;
                int nr_sum=0;

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    if(!other_entry)
                    {
                        if(iterations==1)
                        {
                            if(entries[k].name.compare(nr_values[i][k].name)==0)//nr of values found
                            {
                                mean_sum+=entries[k].mean*nr_values[i][k].nr[0];
                                stdabw_sum+=entries[k].stdabw*nr_values[i][k].nr[0];
                                nr_sum+=nr_values[i][k].nr[0];
                            }
                            else
                            {
                                mean_sum+=entries[k].mean;
                                stdabw_sum+=entries[k].stdabw;
                                nr_sum++;
                            }
                        }
                        else
                        {
                            if(entries[k].name.compare(nr_values_relax[i][k].name)==0)//nr of values found
                            {
                                mean_sum+=entries[k].mean*nr_values_relax[i][k].nr[0];
                                stdabw_sum+=entries[k].stdabw*nr_values_relax[i][k].nr[0];
                                nr_sum+=nr_values_relax[i][k].nr[0];
                            }
                            else
                            {
                                mean_sum+=entries[k].mean;
                                stdabw_sum+=entries[k].stdabw;
                                nr_sum++;
                            }
                        }
                    }
                }

                grain_perimeter_ratio_depths.push_back(depth[i]);
                grain_perimeter_ratio_values.push_back(mean_sum/nr_sum);
                grain_perimeter_ratio_errors.push_back(stdabw_sum/nr_sum);
                grain_perimeter_ratio_file.close();
                temp_grain_perimeter_ratio_file.close();
            }

            std::string filepath_grain_perimeter_ratio2=path_results;
            filepath_grain_perimeter_ratio2.append(label.str());
            filepath_grain_perimeter_ratio2.append("grain_perimeter_ratio2");
            filepath_grain_perimeter_ratio2.append(s.str());

            //grain perimeter ratio
            std::ifstream grain_perimeter_ratio_file2(filepath_grain_perimeter_ratio2.c_str());
            std::ifstream temp_grain_perimeter_ratio_file2(filepath_grain_perimeter_ratio2.c_str());
            temp_grain_perimeter_ratio_file2>>teststring;

            if(grain_perimeter_ratio_file2 && teststring.size()!=0)
            {
                std::vector<depth_parameter> entries;

	            while(!grain_perimeter_ratio_file2.eof())
	            {
                    std::string name;
                    float mean=-1, stdabw;
		            grain_perimeter_ratio_file2>>name>>mean>>stdabw;
                    depth_parameter entry;
                    entry.name=name;
                    entry.mean=mean;
                    entry.stdabw=stdabw;
                    if(mean!=-1) entries.push_back(entry);
	            }

                float mean_sum=0;
                float stdabw_sum=0;
                int nr_sum=0;

                //check for multiple entries for same image
                for (int k=0; k<entries.size(); k++)
                {
                    bool other_entry=false;

                    for(int j=k+1; j<entries.size() && !other_entry; j++)
                    {
                        if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                    }

                    if(!other_entry)
                    {
                        if(iterations==1)
                        {
                            if(entries[k].name.compare(nr_values[i][k].name)==0)//nr of values found
                            {
                                mean_sum+=entries[k].mean*nr_values[i][k].nr[0];
                                stdabw_sum+=entries[k].stdabw*nr_values[i][k].nr[0];
                                nr_sum+=nr_values[i][k].nr[0];
                            }
                            else
                            {
                                mean_sum+=entries[k].mean;
                                stdabw_sum+=entries[k].stdabw;
                                nr_sum++;
                            }
                        }
                        else
                        {
                            if(entries[k].name.compare(nr_values_relax[i][k].name)==0)//nr of values found
                            {
                                mean_sum+=entries[k].mean*nr_values_relax[i][k].nr[0];
                                stdabw_sum+=entries[k].stdabw*nr_values_relax[i][k].nr[0];
                                nr_sum+=nr_values_relax[i][k].nr[0];
                            }
                            else
                            {
                                mean_sum+=entries[k].mean;
                                stdabw_sum+=entries[k].stdabw;
                                nr_sum++;
                            }
                        }
                    }
                }

                grain_perimeter_ratio_depths_bubble.push_back(depth[i]);
                grain_perimeter_ratio_values_bubble.push_back(mean_sum/nr_sum);
                grain_perimeter_ratio_errors_bubble.push_back(stdabw_sum/nr_sum);
                grain_perimeter_ratio_file2.close();
                temp_grain_perimeter_ratio_file2.close();
            }
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

        std::string filepath_grain_perimeter_ratio_combined_out=path_results;
        filepath_grain_perimeter_ratio_combined_out.append("grain_perimeter_ratio_combined");
        filepath_grain_perimeter_ratio_combined_out.append(s.str());
        filepath_grain_perimeter_ratio_combined_out.resize(filepath_grain_perimeter_ratio_combined_out.size()-3);
        filepath_grain_perimeter_ratio_combined_out.append("svg");

        if(minimal_bubble_distance==0 && grainsize_all_values.size()>0) plot.draw_depth_errors("Depth in m", area_unit,
                                                                 "Mean grain size profile (all grains)", grainsize_all_depths,
                                                                 grainsize_all_values, grainsize_all_errors, grainsize_all_errors,
                                                                 0.0f, filepath_grainsize_all_out.c_str());
        else if(grainsize_all_values_relax.size()>0) plot.draw_depth_errors("Depth in m", area_unit,
                                                                 "Mean grain size profile (all grains)", grainsize_all_depths_relax,
                                                                 grainsize_all_values_relax, grainsize_all_errors_relax,
                                                                 grainsize_all_errors_relax, 0.0f,
                                                                 filepath_grainsize_all_out.c_str());
        if(grainsize_all_values.size()>0 && grainsize_all_values_relax.size()>0) plot.draw_depth_errors("Depth in m", area_unit, 
                                                                 "Mean grain size profile (all grains)", grainsize_all_depths,
                                                                 grainsize_all_values, grainsize_all_errors, grainsize_all_errors,
                                                                 grainsize_all_depths_relax, grainsize_all_values_relax,
                                                                 grainsize_all_errors_relax, grainsize_all_errors_relax, 0.0f,
                                                                 filepath_grainsize_all_combined_out.c_str());

        if(minimal_bubble_distance==0 && grainsize_fit_values.size()>0) plot.draw_depth_errors("Depth in m", area_unit, 
                                                                 "Grain size profile log-normal fit", grainsize_fit_depths,
                                                                 grainsize_fit_values, grainsize_fit_errors_low,
                                                                 grainsize_fit_errors_high, 0.0f, filepath_grainsize_fit_out.c_str());
        else if(grainsize_fit_values_relax.size()>0) plot.draw_depth_errors("Depth in m", area_unit,
                                                                 "Grain size profile log-normal fit", grainsize_fit_depths_relax,
                                                                 grainsize_fit_values_relax, grainsize_fit_errors_low_relax,
                                                                 grainsize_fit_errors_high_relax, 0.0f,
                                                                 filepath_grainsize_fit_out.c_str());
        if(grainsize_fit_values.size()>0 && grainsize_fit_values_relax.size()>0) plot.draw_depth_errors("Depth in m", area_unit,
                                                                 "Grain size profile log-normal fit", grainsize_fit_depths,
                                                                 grainsize_fit_values, grainsize_fit_errors_low,
                                                                 grainsize_fit_errors_high, grainsize_fit_depths_relax,
                                                                 grainsize_fit_values_relax, grainsize_fit_errors_low_relax,
                                                                 grainsize_fit_errors_high_relax, 0.0f,
                                                                 filepath_grainsize_fit_combined_out.c_str());

        if(minimal_bubble_distance==0 && corr_ellipse_box_flattening_values.size()>0) plot.draw_depth("Depth in m",
                                                                 "Pearson correlation", "Correlation ellipse box flattening profile",
                                                                 corr_ellipse_box_flattening_depths,
                                                                 corr_ellipse_box_flattening_values, 0.0f,
                                                                 filepath_corr_ellipse_box_flattening_out.c_str());
        else if(corr_ellipse_box_flattening_values_relax.size()>0) plot.draw_depth("Depth in m", "Pearson correlation",
                                                                 "Corelation ellipse box flattening profile",
                                                                 corr_ellipse_box_flattening_depths_relax,
                                                                 corr_ellipse_box_flattening_values_relax, 0.0f,
                                                                 filepath_corr_ellipse_box_flattening_out.c_str());
        if(corr_ellipse_box_flattening_values.size()>0 && corr_ellipse_box_flattening_values_relax.size()>0) plot.draw_depth(
                                                                "Depth in m", "Pearson correlation",
                                                                "Correlation ellipse box profile", corr_ellipse_box_flattening_depths,
                                                                 corr_ellipse_box_flattening_values,
                                                                 corr_ellipse_box_flattening_depths_relax,
                                                                 corr_ellipse_box_flattening_values_relax, 0.0f,
                                                                 filepath_corr_ellipse_box_flattening_combined_out.c_str());

        if(minimal_bubble_distance==0)
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
        else
        for (int st=0; st<grainsize_step_values_relax.size(); st++)
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

            std::string filepath_grainsize_step_combined_out=filepath_grainsize_step_out;
            filepath_grainsize_step_combined_out.append(step);
            filepath_grainsize_step_combined_out.append("step_combined");
            filepath_grainsize_step_combined_out.append(s.str());
            filepath_grainsize_step_combined_out.resize(filepath_grainsize_step_combined_out.size()-3);
            filepath_grainsize_step_combined_out.append("svg");

            if(grainsize_step_values_relax[st].size()>0) plot.draw_depth_errors("Depth in m", area_unit, titel,
                                                                 grainsize_step_depths_relax[st], grainsize_step_values_relax[st],
                                                                 grainsize_step_errors_relax[st], grainsize_step_errors_relax[st],
                                                                 0.0f, filepath_out.c_str());

            if(grainsize_step_values.size()>0)
                if(grainsize_step_values[st].size()>0 && grainsize_step_values_relax[st].size()>0) plot.draw_depth_errors(
                                                                "Depth in m", area_unit, titel, grainsize_step_depths[st],
                                                                 grainsize_step_values[st], grainsize_step_errors[st],
                                                                 grainsize_step_errors[st], grainsize_step_depths_relax[st],
                                                                 grainsize_step_values_relax[st], grainsize_step_errors_relax[st],
                                                                 grainsize_step_errors_relax[st], 0.0f,
                                                                 filepath_grainsize_step_combined_out.c_str());
        }

        if(minimal_bubble_distance==0)
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

            if(minimal_bubble_distance==0 && grainsize_percent_values[p].size()>0) plot.draw_depth_errors("Depth in m", area_unit,
                                                                 titel, grainsize_percent_depths[p], grainsize_percent_values[p],
                                                                 grainsize_percent_errors[p], grainsize_percent_errors[p], 0.0f,
                                                                 filepath_out.c_str());
        }
        else
        for (int p=0; p<grainsize_percent_values_relax.size(); p++)
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

            std::string filepath_grainsize_percent_combined_out=filepath_grainsize_percent_out;
            filepath_grainsize_percent_combined_out.append(percent);
            filepath_grainsize_percent_combined_out.append("percent_combined");
            filepath_grainsize_percent_combined_out.append(s.str());
            filepath_grainsize_percent_combined_out.resize(filepath_grainsize_percent_combined_out.size()-3);
            filepath_grainsize_percent_combined_out.append("svg");

            if(grainsize_percent_values_relax[p].size()>0) plot.draw_depth_errors("Depth in m", area_unit, titel,
                                                                 grainsize_percent_depths_relax[p],
                                                                 grainsize_percent_values_relax[p], grainsize_percent_errors_relax[p],
                                                                 grainsize_percent_errors_relax[p], 0.0f, filepath_out.c_str());

            if(grainsize_percent_values.size()>0)
                if(grainsize_percent_values[p].size()>0 && grainsize_percent_values_relax[p].size()>0) plot.draw_depth_errors(
                                                                 "Depth in m", area_unit, titel, grainsize_percent_depths[p],
                                                                 grainsize_percent_values[p], grainsize_percent_errors[p],
                                                                 grainsize_percent_errors[p], grainsize_percent_depths_relax[p],
                                                                 grainsize_percent_values_relax[p],
                                                                 grainsize_percent_errors_relax[p],
                                                                 grainsize_percent_errors_relax[p], 0.0f,
                                                                 filepath_grainsize_percent_combined_out.c_str());
        }

        if(minimal_bubble_distance==0)
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

            if(minimal_bubble_distance==0 && grainsize_quantile_values[q].size()>0) plot.draw_depth_errors("Depth in m", area_unit,
                                                                 titel, grainsize_quantile_depths[q], grainsize_quantile_values[q],
                                                                 grainsize_quantile_errors[q], grainsize_quantile_errors[q], 0.0f,
                                                                 filepath_out.c_str());
        }
        else
        for (int q=0; q<grainsize_quantile_values_relax.size(); q++)
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

            std::string filepath_grainsize_quantile_combined_out=filepath_grainsize_quantile_out;
            filepath_grainsize_quantile_combined_out.append(quantile);
            filepath_grainsize_quantile_combined_out.append("quantile_combined");
            filepath_grainsize_quantile_combined_out.append(s.str());
            filepath_grainsize_quantile_combined_out.resize(filepath_grainsize_quantile_combined_out.size()-3);
            filepath_grainsize_quantile_combined_out.append("svg");

            if(grainsize_quantile_values_relax[q].size()>0) plot.draw_depth_errors("Depth in m", area_unit, titel,
                                                                 grainsize_quantile_depths_relax[q],
                                                                 grainsize_quantile_values_relax[q],
                                                                 grainsize_quantile_errors_relax[q],
                                                                 grainsize_quantile_errors_relax[q], 0.0f, filepath_out.c_str());

            if(grainsize_quantile_values.size()>0)
                if(grainsize_quantile_values[q].size()>0 && grainsize_quantile_values_relax[q].size()>0) plot.draw_depth_errors(
                                                                 "Depth in m", area_unit, titel, grainsize_quantile_depths[q],
                                                                 grainsize_quantile_values[q], grainsize_quantile_errors[q],
                                                                 grainsize_quantile_errors[q], grainsize_quantile_depths_relax[q],
                                                                 grainsize_quantile_values_relax[q],
                                                                 grainsize_quantile_errors_relax[q],
                                                                 grainsize_quantile_errors_relax[q], 0.0f,
                                                                 filepath_grainsize_quantile_combined_out.c_str());
        }

        for (int b=0; b<std::max(grainsize_bin_values.size(),grainsize_bin_values_relax.size()); b++)
        {
            bool value_found=false;

            if(minimal_bubble_distance==0)
            {
                if (b>=grainsize_bin_values.size()) break;
                for (int depth_index=0; depth_index<grainsize_bin_values[b].size(); depth_index++)
                {
                    if(grainsize_bin_values[b][depth_index]>0.0f) value_found=true;
                }
            }
            else
            {
                if (b>=grainsize_bin_values_relax.size()) break;
                for (int depth_index=0; depth_index<grainsize_bin_values_relax[b].size(); depth_index++)
                {
                    if(grainsize_bin_values_relax[b][depth_index]>0.0f) value_found=true;
                }
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

            if(minimal_bubble_distance==0 && grainsize_bin_values[b].size()>0) plot.draw_depth("Depth in m", "Relative occurrence",
                                                                 titel, grainsize_bin_depths[b], grainsize_bin_values[b], 0.0f,
                                                                 filepath_out.c_str());
            else
            {
                if(grainsize_bin_values_relax[b].size()>0) plot.draw_depth("Depth in m", "Relative occurrence", titel,
                                                                 grainsize_bin_depths_relax[b], grainsize_bin_values_relax[b], 0.0f,
                                                                 filepath_out.c_str());

                std::string filepath_grainsize_bin_combined_out=filepath_grainsize_bin_out;
                filepath_grainsize_bin_combined_out.append(bin);
                filepath_grainsize_bin_combined_out.append("bin_combined");
                filepath_grainsize_bin_combined_out.append(s.str());
                filepath_grainsize_bin_combined_out.resize(filepath_grainsize_bin_combined_out.size()-3);
                filepath_grainsize_bin_combined_out.append("svg");

                if(grainsize_bin_values.size()>0)
                    if(grainsize_bin_values[b].size()>0 && grainsize_bin_values_relax[b].size()>0) plot.draw_depth("Depth in m",
                                                                 "Relative occurrence", titel, grainsize_bin_depths[b],
                                                                 grainsize_bin_values[b], grainsize_bin_depths_relax[b],
                                                                 grainsize_bin_values_relax[b], 0.0f,
                                                                 filepath_grainsize_bin_combined_out.c_str());
            }
        }

        for (int b=0; b<std::max(grainradius_bin_values.size(),grainradius_bin_values_relax.size()); b++)
        {
            bool value_found=false;

            if(minimal_bubble_distance==0)
            {
                if (b>=grainradius_bin_values.size()) break;
                for (int depth_index=0; depth_index<grainradius_bin_values[b].size(); depth_index++)
                {
                    if(grainradius_bin_values[b][depth_index]>0.0f) value_found=true;
                }
            }
            else
            {
                if (b>=grainradius_bin_values_relax.size()) break;
                for (int depth_index=0; depth_index<grainradius_bin_values_relax[b].size(); depth_index++)
                {
                    if(grainradius_bin_values_relax[b][depth_index]>0.0f) value_found=true;
                }
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

            if(minimal_bubble_distance==0 && grainradius_bin_values[b].size()>0) plot.draw_depth("Depth in m", "Relative occurrence",
                                                                 titel, grainradius_bin_depths[b], grainradius_bin_values[b], 0.0f,
                                                                 filepath_out.c_str());
            else
            {
                if(grainradius_bin_values_relax[b].size()>0) plot.draw_depth("Depth in m", "Relative occurrence", titel,
                                                                 grainradius_bin_depths_relax[b], grainradius_bin_values_relax[b],
                                                                 0.0f, filepath_out.c_str());

                std::string filepath_grainradius_bin_combined_out=filepath_grainradius_bin_out;
                filepath_grainradius_bin_combined_out.append(bin);
                filepath_grainradius_bin_combined_out.append("bin_combined");
                filepath_grainradius_bin_combined_out.append(s.str());
                filepath_grainradius_bin_combined_out.resize(filepath_grainradius_bin_combined_out.size()-3);
                filepath_grainradius_bin_combined_out.append("svg");

                if(grainradius_bin_values.size()>0)
                    if(grainradius_bin_values[b].size()>0 && grainradius_bin_values_relax[b].size()>0) plot.draw_depth("Depth in m",
                                                                 "Relative occurrence", titel, grainradius_bin_depths[b],
                                                                 grainradius_bin_values[b], grainradius_bin_depths_relax[b],
                                                                 grainradius_bin_values_relax[b], 0.0f,
                                                                 filepath_grainradius_bin_combined_out.c_str());
            }
        }

        for (int b=0; b<std::max(grainradius_norm_bin_values.size(),grainradius_norm_bin_values_relax.size()); b++)
        {
            bool value_found=false;

            if(minimal_bubble_distance==0)
            {
                if (b>=grainradius_norm_bin_values.size()) break;
                for (int depth_index=0; depth_index<grainradius_norm_bin_values[b].size(); depth_index++)
                {
                    if(grainradius_norm_bin_values[b][depth_index]>0.0f) value_found=true;
                }
            }
            else
            {
                if (b>=grainradius_norm_bin_values_relax.size()) break;
                for (int depth_index=0; depth_index<grainradius_norm_bin_values_relax[b].size(); depth_index++)
                {
                    if(grainradius_norm_bin_values_relax[b][depth_index]>0.0f) value_found=true;
                }
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

            if(minimal_bubble_distance==0 && grainradius_norm_bin_values[b].size()>0) plot.draw_depth("Depth in m",
                                                                 "Relative occurrence", titel, grainradius_norm_bin_depths[b],
                                                                 grainradius_norm_bin_values[b], 0.0f, filepath_out.c_str());
            else
            {
                if(grainradius_norm_bin_values_relax[b].size()>0) plot.draw_depth("Depth in m", "Relative occurrence", titel,
                                                                 grainradius_norm_bin_depths_relax[b],
                                                                 grainradius_norm_bin_values_relax[b], 0.0f, filepath_out.c_str());

                std::string filepath_grainradius_norm_bin_combined_out=filepath_grainradius_norm_bin_out;
                filepath_grainradius_norm_bin_combined_out.append(bin);
                filepath_grainradius_norm_bin_combined_out.append("bin_combined");
                filepath_grainradius_norm_bin_combined_out.append(s.str());
                filepath_grainradius_norm_bin_combined_out.resize(filepath_grainradius_norm_bin_combined_out.size()-3);
                filepath_grainradius_norm_bin_combined_out.append("svg");

                if(grainradius_norm_bin_values.size()>0)
                    if(grainradius_norm_bin_values[b].size()>0 && grainradius_norm_bin_values_relax[b].size()>0) plot.draw_depth(
                                                                 "Depth in m", "Relative occurrence", titel,
                                                                 grainradius_norm_bin_depths[b], grainradius_norm_bin_values[b],
                                                                 grainradius_norm_bin_depths_relax[b],
                                                                 grainradius_norm_bin_values_relax[b], 0.0f,
                                                                 filepath_grainradius_norm_bin_combined_out.c_str());
            }
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

        if(grain_perimeter_ratio_values.size()>0 && grain_perimeter_ratio_values_bubble.size()>0) plot.draw_depth_errors(
                                                                 "Depth in m", "Perimeter ratio", "Perimeter ratio profile",
                                                                 grain_perimeter_ratio_depths, grain_perimeter_ratio_values,
                                                                 grain_perimeter_ratio_errors, grain_perimeter_ratio_errors,
                                                                 grain_perimeter_ratio_depths_bubble,
                                                                 grain_perimeter_ratio_values_bubble,
                                                                 grain_perimeter_ratio_errors_bubble,
                                                                 grain_perimeter_ratio_errors_bubble, 0.7f,
                                                                 filepath_grain_perimeter_ratio_combined_out.c_str());

        if(depth_bin_width>0.0f)
        {
            std::cout<<"New plots generated"<<std::endl;

            //number histograms
            std::vector<int> nr_grains_histogram(1,0);
            std::vector<int> nr_boundaries_histogram(1,0);
            std::vector<int> nr_junctions_histogram(1,0);
            std::vector<int> nr_images_histogram(1,0);
            int depth_max=0;
            int depth_max_relax=depth[66];

            for(int i=0; i<nr_depths; i++)
            {
                if (depth[i]>depth_max)
                {
                    depth_max=depth[i];

                    nr_grains_histogram.resize(1+depth_max/depth_bin_width,0);
                    nr_boundaries_histogram.resize(1+depth_max/depth_bin_width,0);
                    nr_junctions_histogram.resize(1+depth_max/depth_bin_width,0);
                    nr_images_histogram.resize(1+depth_max/depth_bin_width,0);
                }

                if(minimal_bubble_distance==0)
                    for (int k=0; k<nr_values[i].size(); k++)
                    {
                        bool other_entry=false;

                        //check for multiple entries for same image
                        for(int j=k+1; j<nr_values[i].size() && !other_entry; j++)
                        {
                            if(nr_values[i][k].name.compare(nr_values[i][j].name)==0) other_entry=true;
                        }

                        if (!other_entry)
                        {
                            nr_grains_histogram[depth[i]/depth_bin_width]+=nr_values[i][k].nr[0];
                            nr_boundaries_histogram[depth[i]/depth_bin_width]+=nr_values[i][k].nr[8];
                            nr_junctions_histogram[depth[i]/depth_bin_width]+=nr_values[i][k].nr[9];
                            nr_images_histogram[depth[i]/depth_bin_width]++;
                        }
                    }
                else if(i<67)
                    for (int k=0; k<nr_values_relax[i].size(); k++)
                    {
                        bool other_entry=false;

                        //check for multiple entries for same image
                        for(int j=k+1; j<nr_values[i].size() && !other_entry; j++)
                        {
                            if(nr_values[i][k].name.compare(nr_values[i][j].name)==0) other_entry=true;
                        }

                        if (!other_entry)
                        {
                            nr_grains_histogram[depth[i]/depth_bin_width]+=nr_values_relax[i][k].nr[0];
                            nr_boundaries_histogram[depth[i]/depth_bin_width]+=nr_values_relax[i][k].nr[8];
                            nr_junctions_histogram[depth[i]/depth_bin_width]+=nr_values_relax[i][k].nr[9];
                            nr_images_histogram[depth[i]/depth_bin_width]++;
                        }
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

            std::stringstream s;
            if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
            else s << "." << minimal_grain_size;
            if (grain_size_min>10) s << "_" << grain_size_min;
            s << ".txt";

            //load nr values2
            std::vector<std::vector<depth_numbers> > nr_values2;
            std::vector<std::vector<depth_numbers> > nr_values2_relax;

            for(int i=0; i<nr_depths; i++)
            {
                int iterations=1;
                if (minimal_bubble_distance>0 && i<67) iterations=2;

                std::stringstream label;
                label << bagnr[i] << "/";

                //string is read from temp file to check whether file is empty
                std::string teststring;

                for(int iter=0; iter<iterations; iter++)
                {
                    std::stringstream s;

                    if (iter==0) s << "." << minimal_grain_size;
                    else s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                    if (grain_size_min>10) s << "_" << grain_size_min;
                    s << ".txt";

                    std::string filepath_parameter=path_results;
                    filepath_parameter.append(label.str());
                    filepath_parameter.append("nr_percentage_filled_grains");
                    filepath_parameter.append(s.str());

                    std::ifstream parameter_file(filepath_parameter.c_str());
                    std::ifstream temp_parameter_file(filepath_parameter.c_str());
                    temp_parameter_file>>teststring;

                    std::vector<depth_numbers> entries;

                    if(parameter_file && teststring.size()!=0)
                    {
                        while(!parameter_file.eof())
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

                            parameter_file>>name>>nr_70>>nr_80>>nr_90>>nr_91>>nr_92>>nr_93>>nr_94>>nr_95>>nr_96>>nr_97>>nr_98>>nr_99>>nr_100;
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
                            if(nr_70!=-1) entries.push_back(entry);
                        }

                        parameter_file.close();
                        temp_parameter_file.close();
                    }

                    if (iter==0) nr_values2.push_back(entries);
                    else nr_values2_relax.push_back(entries);
                }
            }

            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "grainshape", "grain shape", "Depth in m", "Roundness factor (4pi*Area/Perimeter^2)",
                "Grain roundness profile", 0.35f, 1, nr_values, nr_values_relax, depth_max, depth_bin_width);
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "length_grainarcs", "grain arc length", "Depth in m", length_unit, "Grain boundary length profile",
                0.0f, 5, nr_values, nr_values_relax, depth_max, depth_bin_width, length_scaling());
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "length_grainarcs2", "grain arc length", "Depth in m", length_unit, "Grain boundary length profile",
                0.0f, 5, nr_values_new, nr_values_new_relax, depth_max, depth_bin_width, length_scaling());
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "longest_grainarcs", "grain longest arc length", "Depth in m", length_unit,
                "Grain longest boundary length profile", 0.0f, 6, nr_values, nr_values_relax, depth_max, depth_bin_width,
                length_scaling());
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "dislocation_densities", "dislocation density differences", "Depth in m",
                "Dislocation density difference in m^(-2)", "Profile of mean dislocation density difference at grain boundaries",
                0.0f, 8, nr_values, nr_values_relax, depth_max, depth_bin_width);
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "dislocation_densities2", "dislocation density differences", "Depth in m",
                "Dislocation density difference in m^(-2)", "Profile of mean dislocation density difference at grain boundaries",
                0.0f, 8, nr_values_new, nr_values_new_relax, depth_max, depth_bin_width);
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "grain_equivradius", "equivalent grain radius", "Depth in m", length_unit,
                "Grain equivalent radius profile", 0.0f, 1, nr_values, nr_values_relax, depth_max, depth_bin_width, length_scaling());
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "grain_boxshape", "grain box shape", "Depth in m", "Box flattening factor", "Box flattening profile",
                0.9f, 2, nr_values, nr_values_relax, depth_max, depth_bin_width);
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "grain_boxwidth", "grain box width", "Depth in m", length_unit, "Grain box width profile", 0.0f, 2,
                nr_values, nr_values_relax, depth_max, depth_bin_width, length_scaling());
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "grain_boxheight", "grain box height", "Depth in m", length_unit, "Grain box height profile", 0.0f, 2,
                nr_values, nr_values_relax, depth_max, depth_bin_width, length_scaling());
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "grain_ellipselong", "grain ellipse long axis", "Depth in m", length_unit,
                "Ellipse long axis profile", 0.0f, 4, nr_values, nr_values_relax, depth_max, depth_bin_width, length_scaling());
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "grain_ellipseshape", "grain ellipse shape", "Depth in m", "Ellipse flattening factor",
                "Ellipse flattening profile", 0.9f, 3, nr_values, nr_values_relax, depth_max, depth_bin_width);
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "grain_ellipseangle", "grain ellipse orientation", "Depth in m", "Ellipse long axis angle",
                "Ellipse long axis angle profile", 0.0f, 2, nr_values, nr_values_relax, depth_max, depth_bin_width);
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "boundary_orientation", "boundary orientation", "Depth in m", "Boundary orientation",
                "Boundary orientation profile", 0.0f, 7, nr_values, nr_values_relax, depth_max, depth_bin_width);
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "grainwidth", "grain width", "Depth in m", length_unit, "Grain width profile", 0.0f, 1, nr_values,
                nr_values_relax, depth_max, depth_bin_width, length_scaling());
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "grainheight", "grain height", "Depth in m", length_unit, "Grain height profile", 0.0f, 1, nr_values,
                nr_values_relax, depth_max, depth_bin_width, length_scaling());
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "grainflattening", "vertical grain flattening", "Depth in m", "Vertical grain flattening factor",
                "Vertical grain flattening profile", 0.9f, 1, nr_values, nr_values_relax, depth_max, depth_bin_width);
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "turning_points", "number of turning points", "Depth in m", "Nr of turning points",
                "Grain boundary turning points profile", 0.0f, 9, nr_values, nr_values_relax, depth_max, depth_bin_width);
            new_lin_depth_profile(path_results, bagnr, nr_depths(), minimal_bubble_distance, minimal_grain_size, grain_size_min,
                s.str(), plot, "grain_perimeter_ratio", "perimeter ratio", "Depth in m", "Perimeter ratio",
                "Perimeter ratio profile", 0.7f, 1, nr_values, nr_values_relax, depth_max, depth_bin_width);

            std::vector<float> grainsize_all_depths((size_t)(1+depth_max/depth_bin_width));
            std::vector<float> grainsize_all_values((size_t)(1+depth_max/depth_bin_width),0.0f);
            std::vector<float> grainsize_all_sum((size_t)(1+depth_max/depth_bin_width),0.0f);

            std::vector< std::vector<float> > grainsize_step_depths;
            std::vector< std::vector<float> > grainsize_step_values;
            std::vector< std::vector<float> > grainsize_step_sum(500);
            for(int s=0; s<500; s++) grainsize_step_sum[s].resize((size_t)(1+depth_max/depth_bin_width),0.0f);

            std::vector< std::vector<float> > grainsize_bin_depths;
            std::vector< std::vector<float> > grainsize_bin_values;
            std::vector<int> grainsize_bin_sum((size_t)(1+depth_max/depth_bin_width),0);

            std::vector< std::vector<float> > grainradius_bin_depths;
            std::vector< std::vector<float> > grainradius_bin_values;
            std::vector<int> grainradius_bin_sum((size_t)(1+depth_max/depth_bin_width),0);
            
            std::vector< std::vector<float> > grainradius_norm_bin_depths;
            std::vector< std::vector<float> > grainradius_norm_bin_values;
            std::vector<int> grainradius_norm_bin_sum((size_t)(1+depth_max/depth_bin_width),0);

            std::vector<float> grainsize_all_depths_relax((size_t)(1+depth_max/depth_bin_width));
            std::vector<float> grainsize_all_values_relax((size_t)(1+depth_max/depth_bin_width),0.0f);
            std::vector<float> grainsize_all_sum_relax((size_t)(1+depth_max/depth_bin_width),0.0f);

            std::vector< std::vector<float> > grainsize_step_depths_relax;
            std::vector< std::vector<float> > grainsize_step_values_relax;
            std::vector< std::vector<float> > grainsize_step_sum_relax(500);
            for(int s=0; s<500; s++) grainsize_step_sum_relax[s].resize((size_t)(1+depth_max_relax/depth_bin_width),0.0f);

            std::vector< std::vector<float> > grainsize_bin_depths_relax;
            std::vector< std::vector<float> > grainsize_bin_values_relax;
            std::vector<int> grainsize_bin_sum_relax((size_t)(1+depth_max_relax/depth_bin_width),0);

            std::vector< std::vector<float> > grainradius_bin_depths_relax;
            std::vector< std::vector<float> > grainradius_bin_values_relax;
            std::vector<int> grainradius_bin_sum_relax((size_t)(1+depth_max_relax/depth_bin_width),0);

            std::vector< std::vector<float> > grainradius_norm_bin_depths_relax;
            std::vector< std::vector<float> > grainradius_norm_bin_values_relax;
            std::vector<int> grainradius_norm_bin_sum_relax((size_t)(1+depth_max_relax/depth_bin_width),0);

            std::vector<float> grainsize_step_profile_depths((size_t)(1+depth_max/depth_bin_width));
            std::vector<float> grainsize_step_profile_x_values;
            for(int s=0; s<500; s++) grainsize_step_profile_x_values.push_back(s+1);
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
            for(int s=0; s<500; s++) disl_dens_longest_boundaries_sum[s].resize((size_t)(1+depth_max/depth_bin_width),0.0f);

            std::vector<float> disl_dens_longest_boundaries_profile_depths((size_t)(1+depth_max/depth_bin_width));
            std::vector<float> disl_dens_longest_boundaries_profile_x_values;
            for(int s=0; s<500; s++) disl_dens_longest_boundaries_profile_x_values.push_back(s+1);
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

            std::vector< std::vector<float> > grain_perimeter_ratio_depths;
            grain_perimeter_ratio_depths.resize(2);
            for(int p=0; p<2; p++) grain_perimeter_ratio_depths[p].resize((size_t)(1+depth_max/depth_bin_width));
            std::vector< std::vector<float> > grain_perimeter_ratio_values;
            grain_perimeter_ratio_values.resize(2);
            for(int p=0; p<2; p++) grain_perimeter_ratio_values[p].resize((size_t)(1+depth_max/depth_bin_width),0.0f);
            std::vector< std::vector<float> > grain_perimeter_ratio_sum;
            grain_perimeter_ratio_sum.resize(2);
            for(int p=0; p<2; p++) grain_perimeter_ratio_sum[p].resize((size_t)(1+depth_max/depth_bin_width),0.0f);

            for(int i=0; i<nr_depths; i++)
            {
                int iterations=1;
                if (minimal_bubble_distance>0 && i<67) iterations=2;

                depth[i]=0.55*bagnr[i];

                std::stringstream label;
                label << bagnr[i] << "/";

                //string is read from temp file to check whether file is empty
                std::string teststring;

                for(int iter=0; iter<iterations; iter++)
                {
                    std::stringstream s;

                    if (iter==0) s << "." << minimal_grain_size;
                    else s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                    if (grain_size_min>10) s << "_" << grain_size_min;
                    s << ".txt";

                    std::string filepath_grainsize_all=path_results;
                    filepath_grainsize_all.append(label.str());
                    filepath_grainsize_all.append("grainsize_all");
                    filepath_grainsize_all.append(s.str());

                    std::string filepath_grainsize_step=path_results;
                    filepath_grainsize_step.append(label.str());
                    filepath_grainsize_step.append("grainsize_step");
                    filepath_grainsize_step.append(s.str());

                    std::string filepath_grainsize_bin=path_results;
                    filepath_grainsize_bin.append(label.str());
                    filepath_grainsize_bin.append("grainsize_bin");
                    filepath_grainsize_bin.append(s.str());

                    std::string filepath_grainradius_bin=path_results;
                    filepath_grainradius_bin.append(label.str());
                    filepath_grainradius_bin.append("grainradius_bin");
                    filepath_grainradius_bin.append(s.str());

                    std::string filepath_grainradius_norm_bin=path_results;
                    filepath_grainradius_norm_bin.append(label.str());
                    filepath_grainradius_norm_bin.append("grainradius_norm_bin");
                    filepath_grainradius_norm_bin.append(s.str());

                    std::string filepath_grain_neighbors_bin=path_results;
                    filepath_grain_neighbors_bin.append(label.str());
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

                            if (iter==0)
                            {
                                if(!other_entry && entries[k].name.compare(nr_values[i][k].name)==0)//nr of values found
                                {
                                    grainsize_all_values[depth[i]/depth_bin_width]+=
                                        entries[k].mean*(float)nr_values[i][k].nr[0]/area_scaling();
                                    grainsize_all_sum[depth[i]/depth_bin_width]+=nr_values[i][k].nr[0];
                                }
                            }
                            else
                            {
                                if(!other_entry && entries[k].name.compare(nr_values_relax[i][k].name)==0)//nr of values found
                                {
                                    grainsize_all_values_relax[depth[i]/depth_bin_width]+=
                                        entries[k].mean*(float)nr_values_relax[i][k].nr[0]/area_scaling();
                                    grainsize_all_sum_relax[depth[i]/depth_bin_width]+=nr_values_relax[i][k].nr[0];
                                }
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

                            if(!other_entry)
                            {
                                if (iter==0)
                                {
                                    if (min_step_nr>grainsize_step_values.size()) grainsize_step_values.resize(min_step_nr);
                                    for(int s=0; s<min_step_nr; s++)
                                    {
                                        grainsize_step_values[s].resize((size_t)(1+depth_max/depth_bin_width));
                                        grainsize_step_values[s][depth[i]/depth_bin_width]+=entries[k].mean[s]/area_scaling();

                                        grainsize_step_sum[s][depth[i]/depth_bin_width]++;
                                    }
                                }
                                else
                                {
                                    if (min_step_nr>grainsize_step_values_relax.size())
                                        grainsize_step_values_relax.resize(min_step_nr);
                                    for(int s=0; s<min_step_nr; s++)
                                    {
                                        grainsize_step_values_relax[s].resize((size_t)(1+depth_max/depth_bin_width));
                                        grainsize_step_values_relax[s][depth[i]/depth_bin_width]+=entries[k].mean[s]/area_scaling();

                                        grainsize_step_sum_relax[s][depth[i]/depth_bin_width]++;
                                    }
                                }

                                if (iter==iterations-1)
                                {
                                    if (min_step_nr>grainsize_step_profile_values[depth[i]/depth_bin_width].size())
                                        grainsize_step_profile_values[depth[i]/depth_bin_width].resize(min_step_nr);

                                    for(int s=0; s<min_step_nr; s++)
                                    {
                                        grainsize_step_profile_values[depth[i]/depth_bin_width][s]+=entries[k].mean[s]/area_scaling();
                                    }
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

                            if(!other_entry)
                            {
                                if (iter==0)
                                {
                                    if (min_bin_nr>grainsize_bin_values.size()) grainsize_bin_values.resize(min_bin_nr);
                                    for(int bin=0; bin<min_bin_nr; bin++)
                                    {
                                        grainsize_bin_values[bin].resize((size_t)(1+depth_max/depth_bin_width));
                                        grainsize_bin_values[bin][depth[i]/depth_bin_width]+=entries[k].mean[bin];

                                        grainsize_bin_sum[depth[i]/depth_bin_width]+=entries[k].mean[bin];
                                    }
                                }
                                else
                                {
                                    if (min_bin_nr>grainsize_bin_values_relax.size()) grainsize_bin_values_relax.resize(min_bin_nr);
                                    for(int bin=0; bin<min_bin_nr; bin++)
                                    {
                                        grainsize_bin_values_relax[bin].resize((size_t)(1+depth_max/depth_bin_width));
                                        grainsize_bin_values_relax[bin][depth[i]/depth_bin_width]+=entries[k].mean[bin];

                                        grainsize_bin_sum_relax[depth[i]/depth_bin_width]+=entries[k].mean[bin];
                                    }
                                }

                                if (iter==iterations-1)
                                {
                                    if (min_bin_nr>grainsize_bin_profile_values[depth[i]/depth_bin_width].size())
                                        grainsize_bin_profile_values[depth[i]/depth_bin_width].resize(min_bin_nr);

                                    for(int bin=0; bin<min_bin_nr; bin++)
                                    {
                                        grainsize_bin_profile_values[depth[i]/depth_bin_width][bin]+=entries[k].mean[bin];
                                    }
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

                            if(!other_entry)
                            {
                                if (iter==0)
                                {
                                    if (min_bin_nr>grainradius_bin_values.size()) grainradius_bin_values.resize(min_bin_nr);
                                    for(int bin=0; bin<min_bin_nr; bin++)
                                    {
                                        grainradius_bin_values[bin].resize((size_t)(1+depth_max/depth_bin_width));
                                        grainradius_bin_values[bin][depth[i]/depth_bin_width]+=entries[k].mean[bin];

                                        grainradius_bin_sum[depth[i]/depth_bin_width]+=entries[k].mean[bin];
                                    }
                                }
                                else
                                {
                                    if (min_bin_nr>grainradius_bin_values_relax.size())
                                        grainradius_bin_values_relax.resize(min_bin_nr);
                                    for(int bin=0; bin<min_bin_nr; bin++)
                                    {
                                        grainradius_bin_values_relax[bin].resize((size_t)(1+depth_max/depth_bin_width));
                                        grainradius_bin_values_relax[bin][depth[i]/depth_bin_width]+=entries[k].mean[bin];

                                        grainradius_bin_sum_relax[depth[i]/depth_bin_width]+=entries[k].mean[bin];
                                    }
                                }

                                if (iter==iterations-1)
                                {
                                    if (min_bin_nr>grainradius_bin_profile_values[depth[i]/depth_bin_width].size())
                                        grainradius_bin_profile_values[depth[i]/depth_bin_width].resize(min_bin_nr);

                                    for(int bin=0; bin<min_bin_nr; bin++)
                                    {
                                        grainradius_bin_profile_values[depth[i]/depth_bin_width][bin]+=entries[k].mean[bin];
                                    }
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

                            if(!other_entry)
                            {
                                if (iter==0)
                                {
                                    if (min_bin_nr>grainradius_norm_bin_values.size()) grainradius_norm_bin_values.resize(min_bin_nr);
                                    for(int bin=0; bin<min_bin_nr; bin++)
                                    {
                                        grainradius_norm_bin_values[bin].resize((size_t)(1+depth_max/depth_bin_width));
                                        grainradius_norm_bin_values[bin][depth[i]/depth_bin_width]+=entries[k].mean[bin];

                                        grainradius_norm_bin_sum[depth[i]/depth_bin_width]+=entries[k].mean[bin];
                                    }
                                }
                                else
                                {
                                    if (min_bin_nr>grainradius_norm_bin_values_relax.size())
                                        grainradius_norm_bin_values_relax.resize(min_bin_nr);
                                    for(int bin=0; bin<min_bin_nr; bin++)
                                    {
                                        grainradius_norm_bin_values_relax[bin].resize((size_t)(1+depth_max/depth_bin_width));
                                        grainradius_norm_bin_values_relax[bin][depth[i]/depth_bin_width]+=entries[k].mean[bin];

                                        grainradius_norm_bin_sum_relax[depth[i]/depth_bin_width]+=entries[k].mean[bin];
                                    }
                                }

                                if (iter==iterations-1)
                                {
                                    if (min_bin_nr>grainradius_norm_bin_profile_values[depth[i]/depth_bin_width].size())
                                        grainradius_norm_bin_profile_values[depth[i]/depth_bin_width].resize(min_bin_nr);

                                    for(int bin=0; bin<min_bin_nr; bin++)
                                    {
                                        grainradius_norm_bin_profile_values[depth[i]/depth_bin_width][bin]+=entries[k].mean[bin];
                                    }
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

                            if(!other_entry)
                            {
                                if (iter==iterations-1)
                                {
                                    if (min_bin_nr>grain_neighbors_bin_profile_values[depth[i]/depth_bin_width].size())
                                        grain_neighbors_bin_profile_values[depth[i]/depth_bin_width].resize(min_bin_nr);

                                    for(int bin=0; bin<min_bin_nr; bin++)
                                    {
                                        grain_neighbors_bin_sum[depth[i]/depth_bin_width]+=entries[k].mean[bin];
                                        grain_neighbors_bin_profile_values[depth[i]/depth_bin_width][bin]+=entries[k].mean[bin];
                                    }
                                }
                            }
                        }

                        grain_neighbors_bin_file.close();
                        temp_grain_neighbors_bin_file.close();
                    }

                }

                std::string filepath_dihedral_angles=path_results;
                filepath_dihedral_angles.append(label.str());
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

                        if(iterations==1)
                        {
                            if(!other_entry && entries[k].name.compare(nr_values[i][k].name)==0)//nr of values found
                            {
                                dihedral_angles_errors[depth[i]/depth_bin_width]+=entries[k].stdabw*nr_values[i][k].nr[6];
                                dihedral_angles_sum[depth[i]/depth_bin_width]+=nr_values[i][k].nr[6];
                            }
                        }
                        else if(!other_entry && entries[k].name.compare(nr_values_relax[i][k].name)==0)//nr of values found
                        {
                            dihedral_angles_errors[depth[i]/depth_bin_width]+=entries[k].stdabw*nr_values_relax[i][k].nr[6];
                            dihedral_angles_sum[depth[i]/depth_bin_width]+=nr_values_relax[i][k].nr[6];
                        }
                    }

                    dihedral_angles_file.close();
                    temp_dihedral_angles_file.close();
                }

                std::string filepath_dihedral_angles2=path_results;
                filepath_dihedral_angles2.append(label.str());
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

                        if(iterations==1)
                        {
                            if(!other_entry && entries[k].name.compare(nr_values[i][k].name)==0)//nr of values found
                            {
                                dihedral_angles2_errors[depth[i]/depth_bin_width]+=entries[k].stdabw*nr_values[i][k].nr[6];
                                dihedral_angles2_sum[depth[i]/depth_bin_width]+=nr_values[i][k].nr[6];
                            }
                        }
                        else if(!other_entry && entries[k].name.compare(nr_values_relax[i][k].name)==0)//nr of values found
                        {
                            dihedral_angles2_errors[depth[i]/depth_bin_width]+=entries[k].stdabw*nr_values_relax[i][k].nr[6];
                            dihedral_angles2_sum[depth[i]/depth_bin_width]+=nr_values_relax[i][k].nr[6];
                        }
                    }

                    dihedral_angles2_file.close();
                    temp_dihedral_angles2_file.close();
                }

                std::string filepath_disl_dens_longest_boundaries=path_results;
                filepath_disl_dens_longest_boundaries.append(label.str());
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

                        if(!other_entry)
                        {
                            if (min_step_nr>disl_dens_longest_boundaries_values.size())
                                disl_dens_longest_boundaries_values.resize(min_step_nr);
                            for(int s=0; s<min_step_nr; s++)
                            {
                                disl_dens_longest_boundaries_values[s].resize((size_t)(1+depth_max/depth_bin_width));
                                disl_dens_longest_boundaries_values[s][depth[i]/depth_bin_width]+=entries[k].mean[s];

                                disl_dens_longest_boundaries_sum[s][depth[i]/depth_bin_width]++;
                            }

                            if (min_step_nr>disl_dens_longest_boundaries_profile_values[depth[i]/depth_bin_width].size())
                                disl_dens_longest_boundaries_profile_values[depth[i]/depth_bin_width].resize(min_step_nr);

                            for(int s=0; s<min_step_nr; s++)
                            {
                                disl_dens_longest_boundaries_profile_values[depth[i]/depth_bin_width][s]+=entries[k].mean[s];
                            }
                        }
                    }

                    disl_dens_longest_boundaries_file.close();
                    temp_disl_dens_longest_boundaries_file.close();
                }

                std::string filepath_percentage_filled_grains=path_results;
                filepath_percentage_filled_grains.append(label.str());
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
                        if(entry.mean.size()==0) break;
		            }

                    //check for multiple entries for same image
                    for (int k=0; k<entries.size(); k++)
                    {
                        bool other_entry=false;

                        for(int j=k+1; j<entries.size() && !other_entry; j++)
                        {
                            if(entries[k].name.compare(entries[j].name)==0) other_entry=true;
                        }

                        if(iterations==1)
                        {
                            if(!other_entry && entries[k].name.compare(nr_values2[i][k].name)==0)//nr of values found
                            {
                                percentage_filled_grains_values[0][depth[i]/depth_bin_width]+=
                                    entries[k].mean[1]*nr_values2[i][k].nr2[1];
                                percentage_filled_grains_sum[0][depth[i]/depth_bin_width]+=nr_values2[i][k].nr2[1];
                                percentage_filled_grains_values[1][depth[i]/depth_bin_width]+=
                                    entries[k].mean[2]*nr_values2[i][k].nr2[2];
                                percentage_filled_grains_sum[1][depth[i]/depth_bin_width]+=nr_values2[i][k].nr2[2];
                                percentage_filled_grains_values[2][depth[i]/depth_bin_width]+=
                                    entries[k].mean[7]*nr_values2[i][k].nr2[7];
                                percentage_filled_grains_sum[2][depth[i]/depth_bin_width]+=nr_values2[i][k].nr2[7];
                                percentage_filled_grains_values[3][depth[i]/depth_bin_width]+=
                                    entries[k].mean[12]*nr_values2[i][k].nr2[12];
                                percentage_filled_grains_sum[3][depth[i]/depth_bin_width]+=nr_values2[i][k].nr2[12];
                            }
                        }
                        else if(!other_entry && entries[k].name.compare(nr_values2_relax[i][k].name)==0)//nr of values found
                        {
                            percentage_filled_grains_values[0][depth[i]/depth_bin_width]+=
                                entries[k].mean[1]*nr_values2_relax[i][k].nr2[1];
                            percentage_filled_grains_sum[0][depth[i]/depth_bin_width]+=nr_values2_relax[i][k].nr2[1];
                            percentage_filled_grains_values[1][depth[i]/depth_bin_width]+=
                                entries[k].mean[2]*nr_values2_relax[i][k].nr2[2];
                            percentage_filled_grains_sum[1][depth[i]/depth_bin_width]+=nr_values2_relax[i][k].nr2[2];
                            percentage_filled_grains_values[2][depth[i]/depth_bin_width]+=
                                entries[k].mean[7]*nr_values2_relax[i][k].nr2[7];
                            percentage_filled_grains_sum[2][depth[i]/depth_bin_width]+=nr_values2_relax[i][k].nr2[7];
                            percentage_filled_grains_values[3][depth[i]/depth_bin_width]+=
                                entries[k].mean[12]*nr_values2_relax[i][k].nr2[12];
                            percentage_filled_grains_sum[3][depth[i]/depth_bin_width]+=nr_values2_relax[i][k].nr2[12];
                        }
                    }

                    percentage_filled_grains_file.close();
                    temp_percentage_filled_grains_file.close();
                }

                std::string filepath_grain_perimeter_ratio=path_results;
                filepath_grain_perimeter_ratio.append(label.str());
                filepath_grain_perimeter_ratio.append("grain_perimeter_ratio");
                filepath_grain_perimeter_ratio.append(s.str());

                //grain perimeter ratio
                std::ifstream grain_perimeter_ratio_file(filepath_grain_perimeter_ratio.c_str());
                std::ifstream temp_grain_perimeter_ratio_file(filepath_grain_perimeter_ratio.c_str());
                temp_grain_perimeter_ratio_file>>teststring;

                if(grain_perimeter_ratio_file && teststring.size()!=0)
                {
                    std::vector<depth_parameter> entries;

	                while(!grain_perimeter_ratio_file.eof())
	                {
                        std::string name;
                        float mean=-1, stdabw;
		                grain_perimeter_ratio_file>>name>>mean>>stdabw;
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

                        if(iterations==1)
                        {
                            if(!other_entry && entries[k].name.compare(nr_values[i][k].name)==0)//nr of values found
                            {
                                grain_perimeter_ratio_values[0][depth[i]/depth_bin_width]+=entries[k].mean*nr_values[i][k].nr[0];
                                grain_perimeter_ratio_sum[0][depth[i]/depth_bin_width]+=nr_values[i][k].nr[0];
                            }
                        }
                        else if(!other_entry && entries[k].name.compare(nr_values_relax[i][k].name)==0)//nr of values found
                        {
                            grain_perimeter_ratio_values[0][depth[i]/depth_bin_width]+=entries[k].mean*nr_values_relax[i][k].nr[0];
                            grain_perimeter_ratio_sum[0][depth[i]/depth_bin_width]+=nr_values_relax[i][k].nr[0];
                        }
                    }

                    grain_perimeter_ratio_file.close();
                    temp_grain_perimeter_ratio_file.close();
                }

                std::string filepath_grain_perimeter_ratio2=path_results;
                filepath_grain_perimeter_ratio2.append(label.str());
                filepath_grain_perimeter_ratio2.append("grain_perimeter_ratio2");
                filepath_grain_perimeter_ratio2.append(s.str());

                //grain perimeter ratio
                std::ifstream grain_perimeter_ratio_file2(filepath_grain_perimeter_ratio2.c_str());
                std::ifstream temp_grain_perimeter_ratio_file2(filepath_grain_perimeter_ratio2.c_str());
                temp_grain_perimeter_ratio_file2>>teststring;

                if(grain_perimeter_ratio_file2 && teststring.size()!=0)
                {
                    std::vector<depth_parameter> entries;

	                while(!grain_perimeter_ratio_file2.eof())
	                {
                        std::string name;
                        float mean=-1, stdabw;
		                grain_perimeter_ratio_file2>>name>>mean>>stdabw;
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

                        if(iterations==1)
                        {
                            if(!other_entry && entries[k].name.compare(nr_values[i][k].name)==0)//nr of values found
                            {
                                grain_perimeter_ratio_values[1][depth[i]/depth_bin_width]+=entries[k].mean*nr_values[i][k].nr[0];
                                grain_perimeter_ratio_sum[1][depth[i]/depth_bin_width]+=nr_values[i][k].nr[0];
                            }
                        }
                        else if(!other_entry && entries[k].name.compare(nr_values_relax[i][k].name)==0)//nr of values found
                        {
                            grain_perimeter_ratio_values[1][depth[i]/depth_bin_width]+=entries[k].mean*nr_values_relax[i][k].nr[0];
                            grain_perimeter_ratio_sum[1][depth[i]/depth_bin_width]+=nr_values_relax[i][k].nr[0];
                        }
                    }

                    grain_perimeter_ratio_file2.close();
                    temp_grain_perimeter_ratio_file2.close();
                }
            }

            for (int bin=0; bin<grainsize_all_values.size(); bin++)
            {
                grainsize_all_depths[bin]=(bin+0.5f)*depth_bin_width;
            }

            for (int bin=0; bin<grainsize_all_values_relax.size(); bin++)
            {
                grainsize_all_depths_relax[bin]=(bin+0.5f)*depth_bin_width;
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

            for (int bin=0; bin<grainsize_all_values_relax.size(); bin++)
            {
                if(grainsize_all_sum_relax[bin]>0.0f) grainsize_all_values_relax[bin]/=grainsize_all_sum_relax[bin];
                else
                {
                    grainsize_all_depths_relax.erase(grainsize_all_depths_relax.begin()+bin);
                    grainsize_all_values_relax.erase(grainsize_all_values_relax.begin()+bin);
                    grainsize_all_sum_relax.erase(grainsize_all_sum_relax.begin()+bin);
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

            grainsize_step_depths_relax.resize(grainsize_step_values_relax.size());

            for(int grainsize_step=0; grainsize_step<grainsize_step_values_relax.size(); grainsize_step++)
            {
                grainsize_step_depths_relax[grainsize_step].resize(grainsize_step_values_relax[grainsize_step].size());

                for(int bin=0; bin<grainsize_step_values_relax[grainsize_step].size(); bin++)
                {
                    grainsize_step_depths_relax[grainsize_step][bin]=(bin+0.5f)*depth_bin_width;
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

            for(int grainsize_step=0; grainsize_step<grainsize_step_values_relax.size(); grainsize_step++)
            {
                int count=0;
                for (int bin=0; bin<grainsize_step_values_relax[grainsize_step].size(); bin++)
                {
                    if(grainsize_step_sum_relax[grainsize_step][count]>0.0f)
                        grainsize_step_values_relax[grainsize_step][bin]/=grainsize_step_sum_relax[grainsize_step][count];
                    else
                    {
                        grainsize_step_depths_relax[grainsize_step].erase(grainsize_step_depths_relax[grainsize_step].begin()+bin);
                        grainsize_step_values_relax[grainsize_step].erase(grainsize_step_values_relax[grainsize_step].begin()+bin);
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

            grainsize_bin_depths_relax.resize(grainsize_bin_values_relax.size());

            for(int grainsize_bin=0; grainsize_bin<grainsize_bin_values_relax.size(); grainsize_bin++)
            {
                grainsize_bin_depths_relax[grainsize_bin].resize(grainsize_bin_values_relax[grainsize_bin].size());

                for(int bin=0; bin<grainsize_bin_values_relax[grainsize_bin].size(); bin++)
                {
                    if(grainsize_bin_sum_relax[bin]>0.0f)
                        grainsize_bin_values_relax[grainsize_bin][bin]/=grainsize_bin_sum_relax[bin];
                    grainsize_bin_depths_relax[grainsize_bin][bin]=(bin+0.5f)*depth_bin_width;
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

            if (minimal_bubble_distance>0)
            {
                int end_bin=0;
                for (int bin=0; bin<grainsize_bin_sum_relax.size(); bin++)
                {
                    if(grainsize_bin_sum_relax[bin]>0.0f)
                    {
                        for(int grainsize_bin=0; grainsize_bin<grainsize_bin_profile_values[bin].size(); grainsize_bin++)
                            grainsize_bin_profile_values[bin][grainsize_bin]/=grainsize_bin_sum_relax[bin];
                    }
                    else
                    {
                        grainsize_bin_profile_depths.erase(grainsize_bin_profile_depths.begin()+bin);
                        grainsize_bin_profile_values.erase(grainsize_bin_profile_values.begin()+bin);
                        grainsize_bin_sum.erase(grainsize_bin_sum.begin()+bin);
                        grainsize_bin_sum_relax.erase(grainsize_bin_sum_relax.begin()+bin);
                        bin--;
                    }
                    end_bin=bin;
                }

                for (int bin=end_bin; bin<grainsize_bin_profile_values.size(); bin++)
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
            }
            else for (int bin=0; bin<grainsize_bin_profile_values.size(); bin++)
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

            grainradius_bin_depths_relax.resize(grainradius_bin_values_relax.size());

            for(int grainradius_bin=0; grainradius_bin<grainradius_bin_values_relax.size(); grainradius_bin++)
            {
                grainradius_bin_depths_relax[grainradius_bin].resize(grainradius_bin_values_relax[grainradius_bin].size());

                for(int bin=0; bin<grainradius_bin_values_relax[grainradius_bin].size(); bin++)
                {
                    if(grainradius_bin_sum_relax[bin]>0.0f)
                        grainradius_bin_values_relax[grainradius_bin][bin]/=grainradius_bin_sum_relax[bin];
                    grainradius_bin_depths_relax[grainradius_bin][bin]=(bin+0.5f)*depth_bin_width;
                }
            }

            for (int bin=0; bin<grainradius_bin_profile_values.size(); bin++)
            {
                grainradius_bin_profile_depths[bin]=(bin+0.5f)*depth_bin_width;
            }

            if (minimal_bubble_distance>0)
            {
                int end_bin=0;
                for (int bin=0; bin<grainradius_bin_sum_relax.size(); bin++)
                {
                    if(grainradius_bin_sum_relax[bin]>0.0f)
                    {
                        for(int grainradius_bin=0; grainradius_bin<grainradius_bin_profile_values[bin].size(); grainradius_bin++)
                            grainradius_bin_profile_values[bin][grainradius_bin]/=grainradius_bin_sum_relax[bin];
                    }
                    else
                    {
                        grainradius_bin_profile_depths.erase(grainradius_bin_profile_depths.begin()+bin);
                        grainradius_bin_profile_values.erase(grainradius_bin_profile_values.begin()+bin);
                        grainradius_bin_sum.erase(grainradius_bin_sum.begin()+bin);
                        grainradius_bin_sum_relax.erase(grainradius_bin_sum_relax.begin()+bin);
                        bin--;
                    }
                    end_bin=bin;
                }

                for (int bin=end_bin; bin<grainradius_bin_profile_values.size(); bin++)
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
            }
            else for (int bin=0; bin<grainradius_bin_profile_values.size(); bin++)
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

            grainradius_norm_bin_depths_relax.resize(grainradius_norm_bin_values_relax.size());

            for(int grainradius_norm_bin=0; grainradius_norm_bin<grainradius_norm_bin_values_relax.size(); grainradius_norm_bin++)
            {
                grainradius_norm_bin_depths_relax[grainradius_norm_bin].resize
                    (grainradius_norm_bin_values_relax[grainradius_norm_bin].size());

                for(int bin=0; bin<grainradius_norm_bin_values_relax[grainradius_norm_bin].size(); bin++)
                {
                    if(grainradius_norm_bin_sum_relax[bin]>0.0f)
                        grainradius_norm_bin_values_relax[grainradius_norm_bin][bin]/=grainradius_norm_bin_sum_relax[bin];
                    grainradius_norm_bin_depths_relax[grainradius_norm_bin][bin]=(bin+0.5f)*depth_bin_width;
                }
            }

            for (int bin=0; bin<grainradius_norm_bin_profile_values.size(); bin++)
            {
                grainradius_norm_bin_profile_depths[bin]=(bin+0.5f)*depth_bin_width;
            }

            if (minimal_bubble_distance>0)
            {
                int end_bin=0;
                for (int bin=0; bin<grainradius_norm_bin_sum_relax.size(); bin++)
                {
                    if(grainradius_norm_bin_sum_relax[bin]>0.0f)
                    {
                        for(int grainradius_norm_bin=0; grainradius_norm_bin<grainradius_norm_bin_profile_values[bin].size();
                            grainradius_norm_bin++)
                            grainradius_norm_bin_profile_values[bin][grainradius_norm_bin]/=grainradius_norm_bin_sum_relax[bin];
                    }
                    else
                    {
                        grainradius_norm_bin_profile_depths.erase(grainradius_norm_bin_profile_depths.begin()+bin);
                        grainradius_norm_bin_profile_values.erase(grainradius_norm_bin_profile_values.begin()+bin);
                        grainradius_norm_bin_sum.erase(grainradius_norm_bin_sum.begin()+bin);
                        grainradius_norm_bin_sum_relax.erase(grainradius_norm_bin_sum_relax.begin()+bin);
                        bin--;
                    }
                    end_bin=bin;
                }

                for (int bin=end_bin; bin<grainradius_norm_bin_profile_values.size(); bin++)
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
            }
            else for (int bin=0; bin<grainradius_norm_bin_profile_values.size(); bin++)
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

            for(int p=0; p<2; p++)
            {
                for (int bin=0; bin<grain_perimeter_ratio_values[p].size(); bin++)
                {
                    grain_perimeter_ratio_depths[p][bin]=(bin+0.5f)*depth_bin_width;
                }

                for (int bin=0; bin<grain_perimeter_ratio_values[p].size(); bin++)
                {
                    if(grain_perimeter_ratio_sum[p][bin]>0.0f)
                        grain_perimeter_ratio_values[p][bin]/=grain_perimeter_ratio_sum[p][bin];
                    else
                    {
                        grain_perimeter_ratio_depths[p].erase(grain_perimeter_ratio_depths[p].begin()+bin);
                        grain_perimeter_ratio_values[p].erase(grain_perimeter_ratio_values[p].begin()+bin);
                        grain_perimeter_ratio_sum[p].erase(grain_perimeter_ratio_sum[p].begin()+bin);
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

            std::string filepath_grain_perimeter_ratio_combined_out=path_results;
            filepath_grain_perimeter_ratio_combined_out.append("new_grain_perimeter_ratio_combined");
            filepath_grain_perimeter_ratio_combined_out.append(s.str());
            filepath_grain_perimeter_ratio_combined_out.resize(filepath_grain_perimeter_ratio_combined_out.size()-3);
            filepath_grain_perimeter_ratio_combined_out.append("svg");

            if(minimal_bubble_distance==0 && grainsize_all_values.size()>0) plot.draw_depth("Depth in m", area_unit,
                "Mean grain size profile (all grains)", grainsize_all_depths, grainsize_all_values, 0.0f,
                filepath_grainsize_all_out.c_str());
            else if(grainsize_all_values_relax.size()>0) plot.draw_depth("Depth in m", area_unit,
                "Mean grain size profile (all grains)", grainsize_all_depths_relax, grainsize_all_values_relax, 0.0f,
                filepath_grainsize_all_out.c_str());
            if(grainsize_all_values.size()>0 && grainsize_all_values_relax.size()>0) plot.draw_depth("Depth in m", area_unit,
                "Mean grain size profile (all grains)", grainsize_all_depths, grainsize_all_values, grainsize_all_depths_relax,
                grainsize_all_values_relax, 0.0f, filepath_grainsize_all_combined_out.c_str());

            if(minimal_bubble_distance==0)
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
            else
            for (int st=0; st<grainsize_step_values_relax.size(); st++)
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

                std::string filepath_grainsize_step_combined_out=filepath_grainsize_step_out;
                filepath_grainsize_step_combined_out.append(step);
                filepath_grainsize_step_combined_out.append("step_combined");
                filepath_grainsize_step_combined_out.append(s.str());
                filepath_grainsize_step_combined_out.resize(filepath_grainsize_step_combined_out.size()-3);
                filepath_grainsize_step_combined_out.append("svg");

                if(grainsize_step_values_relax[st].size()>0) plot.draw_depth("Depth in m", area_unit, titel, 
                    grainsize_step_depths_relax[st], grainsize_step_values_relax[st], 0.0f, filepath_out.c_str());

                if(grainsize_step_values.size()>0)
                    if(grainsize_step_values[st].size()>0 && grainsize_step_values_relax[st].size()>0) plot.draw_depth("Depth in m",
                        area_unit, titel, grainsize_step_depths[st], grainsize_step_values[st], grainsize_step_depths_relax[st],
                        grainsize_step_values_relax[st], 0.0f, filepath_grainsize_step_combined_out.c_str());
            }

            for (int b=0; b<std::max(grainsize_bin_values.size(),grainsize_bin_values_relax.size()); b++)
            {
                bool value_found=false;

                if(minimal_bubble_distance==0)
                {
                    if (b>=grainsize_bin_values.size()) break;
                    for (int depth_index=0; depth_index<grainsize_bin_values[b].size(); depth_index++)
                    {
                        if(grainsize_bin_values[b][depth_index]>0.0f) value_found=true;
                    }
                }
                else
                {
                    if (b>=grainsize_bin_values_relax.size()) break;
                    for (int depth_index=0; depth_index<grainsize_bin_values_relax[b].size(); depth_index++)
                    {
                        if(grainsize_bin_values_relax[b][depth_index]>0.0f) value_found=true;
                    }
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

                if(minimal_bubble_distance==0 && grainsize_bin_values[b].size()>0) plot.draw_depth("Depth in m",
                    "Relative occurrence", titel, grainsize_bin_depths[b], grainsize_bin_values[b], 0.0f, filepath_out.c_str());
                else
                {
                    if(grainsize_bin_values_relax[b].size()>0) plot.draw_depth("Depth in m", "Relative occurrence", titel,
                        grainsize_bin_depths_relax[b], grainsize_bin_values_relax[b], 0.0f, filepath_out.c_str());

                    std::string filepath_grainsize_bin_combined_out=filepath_grainsize_bin_out;
                    filepath_grainsize_bin_combined_out.append(bin);
                    filepath_grainsize_bin_combined_out.append("bin_combined");
                    filepath_grainsize_bin_combined_out.append(s.str());
                    filepath_grainsize_bin_combined_out.resize(filepath_grainsize_bin_combined_out.size()-3);
                    filepath_grainsize_bin_combined_out.append("svg");

                    if(grainsize_bin_values.size()>0)
                        if(grainsize_bin_values[b].size()>0 && grainsize_bin_values_relax[b].size()>0) plot.draw_depth("Depth in m",
                            "Relative occurrence", titel, grainsize_bin_depths[b], grainsize_bin_values[b],
                            grainsize_bin_depths_relax[b], grainsize_bin_values_relax[b], 0.0f,
                            filepath_grainsize_bin_combined_out.c_str());
                }
            }

            for (int b=0; b<std::max(grainradius_bin_values.size(),grainradius_bin_values_relax.size()); b++)
            {
                bool value_found=false;

                if(minimal_bubble_distance==0)
                {
                    if (b>=grainradius_bin_values.size()) break;
                    for (int depth_index=0; depth_index<grainradius_bin_values[b].size(); depth_index++)
                    {
                        if(grainradius_bin_values[b][depth_index]>0.0f) value_found=true;
                    }
                }
                else
                {
                    if (b>=grainradius_bin_values_relax.size()) break;
                    for (int depth_index=0; depth_index<grainradius_bin_values_relax[b].size(); depth_index++)
                    {
                        if(grainradius_bin_values_relax[b][depth_index]>0.0f) value_found=true;
                    }
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

                if(minimal_bubble_distance==0 && grainradius_bin_values[b].size()>0) plot.draw_depth("Depth in m",
                    "Relative occurrence", titel, grainradius_bin_depths[b], grainradius_bin_values[b], 0.0f, filepath_out.c_str());
                else
                {
                    if(grainradius_bin_values_relax[b].size()>0) plot.draw_depth("Depth in m", "Relative occurrence", titel,
                        grainradius_bin_depths_relax[b], grainradius_bin_values_relax[b], 0.0f, filepath_out.c_str());

                    std::string filepath_grainradius_bin_combined_out=filepath_grainradius_bin_out;
                    filepath_grainradius_bin_combined_out.append(bin);
                    filepath_grainradius_bin_combined_out.append("bin_combined");
                    filepath_grainradius_bin_combined_out.append(s.str());
                    filepath_grainradius_bin_combined_out.resize(filepath_grainradius_bin_combined_out.size()-3);
                    filepath_grainradius_bin_combined_out.append("svg");

                    if(grainradius_bin_values.size()>0)
                        if(grainradius_bin_values[b].size()>0 && grainradius_bin_values_relax[b].size()>0) plot.draw_depth(
                            "Depth in m", "Relative occurrence", titel, grainradius_bin_depths[b], grainradius_bin_values[b],
                            grainradius_bin_depths_relax[b], grainradius_bin_values_relax[b], 0.0f,
                            filepath_grainradius_bin_combined_out.c_str());
                }
            }

            for (int b=0; b<std::max(grainradius_norm_bin_values.size(),grainradius_norm_bin_values_relax.size()); b++)
            {
                bool value_found=false;

                if(minimal_bubble_distance==0)
                {
                    if (b>=grainradius_norm_bin_values.size()) break;
                    for (int depth_index=0; depth_index<grainradius_norm_bin_values[b].size(); depth_index++)
                    {
                        if(grainradius_norm_bin_values[b][depth_index]>0.0f) value_found=true;
                    }
                }
                else
                {
                    if (b>=grainradius_norm_bin_values_relax.size()) break;
                    for (int depth_index=0; depth_index<grainradius_norm_bin_values_relax[b].size(); depth_index++)
                    {
                        if(grainradius_norm_bin_values_relax[b][depth_index]>0.0f) value_found=true;
                    }
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

                if(minimal_bubble_distance==0 && grainradius_norm_bin_values[b].size()>0) plot.draw_depth("Depth in m",
                    "Relative occurrence", titel, grainradius_norm_bin_depths[b], grainradius_norm_bin_values[b], 0.0f,
                    filepath_out.c_str());
                else
                {
                    if(grainradius_norm_bin_values_relax[b].size()>0) plot.draw_depth("Depth in m", "Relative occurrence", titel,
                        grainradius_norm_bin_depths_relax[b], grainradius_norm_bin_values_relax[b], 0.0f, filepath_out.c_str());

                    std::string filepath_grainradius_norm_bin_combined_out=filepath_grainradius_norm_bin_out;
                    filepath_grainradius_norm_bin_combined_out.append(bin);
                    filepath_grainradius_norm_bin_combined_out.append("bin_combined");
                    filepath_grainradius_norm_bin_combined_out.append(s.str());
                    filepath_grainradius_norm_bin_combined_out.resize(filepath_grainradius_norm_bin_combined_out.size()-3);
                    filepath_grainradius_norm_bin_combined_out.append("svg");

                    if(grainradius_norm_bin_values.size()>0)
                        if(grainradius_norm_bin_values[b].size()>0 && grainradius_norm_bin_values_relax[b].size()>0)
                            plot.draw_depth("Depth in m", "Relative occurrence", titel, grainradius_norm_bin_depths[b],
                                grainradius_norm_bin_values[b], grainradius_norm_bin_depths_relax[b],
                                grainradius_norm_bin_values_relax[b], 0.0f, filepath_grainradius_norm_bin_combined_out.c_str());
                }
            }

            char step[20];
            sprintf(step, "%d", grain_step);

            std::string x_axis="Largest grains * ";
            x_axis.append(step);

            if(grainsize_step_profile_values.size()>0) plot.draw_depth_3d(x_axis, "Depth in m", area_unit,
                "Largest grains profile (absolute)", grainsize_step_profile_x_values, grainsize_step_profile_depths,
                grainsize_step_profile_values, filepath_grainsize_step_profile_out.c_str());
            if(grainsize_bin_profile_values.size()>0) plot.draw_depth_color(area_unit, "Depth in m", "Relative occurrence",
                "Grain size profile", grainsize_bin_profile_x_values, grainsize_bin_profile_depths, grainsize_bin_profile_values,
                filepath_grainsize_bin_profile_out.c_str());
            if(grainradius_bin_profile_values.size()>0) plot.draw_depth_color(length_unit, "Depth in m", "Relative occurrence",
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
            if(grain_perimeter_ratio_values.size()>0) plot.draw_depth_nofit("Depth in m", "Perimeter ratio",
                "Perimeter ratio profile", grain_perimeter_ratio_depths, grain_perimeter_ratio_values,
                filepath_grain_perimeter_ratio_combined_out.c_str(), 0.7f);
        }
    }

    delete depth;
}
