/*! \file analyze.h
 *  \brief Used for image analyzation.
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
#pragma once
#include <time.h>

void analyze(std::string filepath_list, std::string path_rf_predictions, std::string path_parameters, std::string param_file,
             ParameterFile paramFile, int minimal_grain_size, int minimal_bubble_distance, std::string path_plots, int mode,
             int combo1, int spin, int combo2, std::string path_to_ws_image, bool correct_suffix)
{
    Parameter<float> length_scaling;
    length_scaling.assign("", "length_scaling", 193.5f);
    length_scaling.load(paramFile,"config");

    Parameter<float> area_scaling;
    area_scaling.assign("", "area_scaling", 37444.0f);
    area_scaling.load(paramFile,"config");

    //initialise plplot class
    plplot plot = plplot();

    std::ifstream list_file(filepath_list.c_str());

    if(!list_file)
    {
        std::cout<<"Error: List file could not be found!"<<std::endl;
        return;
    }

    std::ifstream temp_list_file(filepath_list.c_str());
    std::string teststring;
    temp_list_file>>teststring;

    if(teststring.size()==0)
    {
        std::cout<<"Error: List file is empty!"<<std::endl;
        list_file.close();
        temp_list_file.close();
        return;
    }

    //These data structures are extracted and saved for all images in list file
    std::vector< std::vector<int> > grain_junctions;
    std::vector< std::vector<point> > grain_area_center_mass;
    std::vector< std::vector<float> > grain_size_values;
    std::vector<std::string> filenames;
    std::vector< std::vector<int> > grain_areas;
    std::vector< std::vector<float> > grain_roundness;
    std::vector< std::vector<float> > grain_box_flattening;
    std::vector< std::vector<float> > grain_box_width;
    std::vector< std::vector<float> > grain_box_height;
    std::vector< std::vector<float> > grain_ellipse_flattening;
    std::vector< std::vector<float> > grain_ellipse_long_axis;
    std::vector< std::vector<float> > grain_ellipse_long_axis_angle;
    std::vector< std::vector<float> > grain_area_flattening;
    std::vector< std::vector<float> > grain_area_width;
    std::vector< std::vector<float> > grain_area_height;
    std::vector< std::vector<float> > grain_arc_number;
    std::vector< std::vector<float> > grain_neighbors;
    std::vector< std::vector<float> > grain_longest_arc_length;
    std::vector< std::vector<float> > grain_boundary_length;
    std::vector< std::vector<float> > turning_point;
    std::vector< std::vector<float> > stddev_dihedral_angles;
    std::vector< std::vector<float> > grain_perimeter_ratio;

    std::string old="--";
    int line=1;
    
    //Store all file names in the filenames vector
    while(!list_file.eof())
    {
        std::string filename;
        list_file >> filename;

        if(filename == old || filename == "")
        {
            break;
        }
        line++;
        old = filename;
        filenames.push_back(filename);
    }

    //Plot correlations
    if(mode == 0) 
    {  
        //Calculate correlations
        correlation corr_roundness_box_flattening;
        correlation corr_roundness_vertical_flattening;
        correlation corr_ellipse_box_flattening;
        correlation corr_ellipse_long_box_width;
        correlation corr_flattening_abs_ellipse_angle;
        correlation corr_roundness_perimeter_ratio;
        correlation corr_vertical_box_flattening;

        //Loop over every image specified in the list file
        for(int i = 0; i < filenames.size(); i++)
        {
            std::vector<float> float_vector;
            std::string filename = filenames[i];

            //extract suffix
            std::string suffix1=get_path(filename);

            if (suffix1.size()>0)
            {
                suffix1.resize(suffix1.size()-1);
                suffix1=get_filename(suffix1.c_str());
                suffix1.append("/");
            }
            
            std::string param_file_name=get_filename(param_file);
            if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
            else param_file_name = "";

            std::string filepath_new_classification=path_rf_predictions;
            if (correct_suffix) filepath_new_classification.append(suffix1);
            filepath_new_classification.append(param_file_name.c_str());

            std::stringstream s;
            if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
            else s << "." << minimal_grain_size;
            std::string suffix2=s.str();

            filepath_new_classification.append(get_filename(filename));
            if (minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
            filepath_new_classification.append(".h5");

            //IMPORT RESULTS FROM HDF5 file
            gbn grainBoundNet;
            grainBoundNet.load_final_structure(filepath_new_classification);                

            std::string filepath_parameters=path_parameters;
            if (correct_suffix) filepath_parameters.append(suffix1);
            filepath_parameters.append(get_filename(filename));
            filepath_parameters.append(suffix2.c_str());
            filepath_parameters.append(".h5");
            
            param para;

            //gbn attribute is required as param attribute
            para.grain_junctions=grainBoundNet.grain_junctions; 
            para.load_extracted_parameters(filepath_parameters);            

            grain_area_flattening.push_back(float_vector);
            grain_areas.push_back(para.grain_areas);
            grain_box_flattening.push_back(para.grain_box_flattening);    
            grain_ellipse_flattening.push_back(para.ellipse_flattening); 
            grain_roundness.push_back(para.grain_roundness); 
            grain_area_width.push_back(para.grain_area_width); 
            grain_area_height.push_back(para.grain_area_height);
            grain_box_width.push_back(para.grain_box_width);
            grain_ellipse_long_axis.push_back(para.ellipse_long_axis);
            grain_ellipse_long_axis_angle.push_back(para.ellipse_long_axis_angle);
            grain_perimeter_ratio.push_back(para.grain_perimeter_ratio);  

            grain_area_flattening.back().resize(grain_areas.back().back()+1);
            for(int area=0; area<grain_areas.back().size(); area++)
            {
                grain_area_flattening.back()[grain_areas.back()[area]]=
                    grain_area_width.back()[grain_areas.back()[area]]/grain_area_height.back()[grain_areas.back()[area]];
            }
        }

        grain_correlation(grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_box_flattening, grain_roundness,
            corr_roundness_box_flattening);
        grain_correlation(grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_area_flattening, grain_roundness,
            corr_roundness_vertical_flattening);
        grain_correlation(grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_box_flattening, grain_ellipse_flattening, 
            corr_ellipse_box_flattening);
        grain_correlation(grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_box_width, grain_ellipse_long_axis,
            corr_ellipse_long_box_width);
        grain_correlation(grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_ellipse_long_axis_angle, grain_area_flattening, 
            corr_flattening_abs_ellipse_angle, 1, 0);
        grain_correlation(grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_perimeter_ratio, grain_roundness,
            corr_roundness_perimeter_ratio);
        grain_correlation(grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_box_flattening, grain_area_flattening,
            corr_vertical_box_flattening);

        std::stringstream s;
        //if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size <<".svg";
        //else s << "." << minimal_grain_size <<".svg";
        if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size <<".ps";
        else s << "." << minimal_grain_size <<".ps";

        std::string filepath_roundness_box_flattening = path_plots;
        std::string filepath_roundness_vertical_flattening = path_plots;
        std::string filepath_ellipse_box_flattening = path_plots;
        std::string filepath_ellipse_long_box_width = path_plots;
        std::string filepath_ellipse_abs_angle_flattening = path_plots;
        std::string filepath_roundness_perimeter_ratio = path_plots;
        std::string filepath_vertical_box_flattening = path_plots;

        filepath_roundness_box_flattening.append("corr_roundness_box_flattening");
        filepath_roundness_vertical_flattening.append("corr_roundness_vertical_flattening");
        filepath_ellipse_box_flattening.append("corr_ellipse_box_flattening");
        filepath_ellipse_long_box_width.append("corr_ellipse_long_box_width");
        filepath_ellipse_abs_angle_flattening.append("corr_flattening_abs_ellipse_angle");
        filepath_roundness_perimeter_ratio.append("corr_roundness_perimeter_ratio");
        filepath_vertical_box_flattening.append("corr_vertical_box_flattening");

        filepath_roundness_box_flattening.append(s.str());
        filepath_roundness_vertical_flattening.append(s.str());
        filepath_ellipse_box_flattening.append(s.str());
        filepath_ellipse_long_box_width.append(s.str());
        filepath_ellipse_abs_angle_flattening.append(s.str());
        filepath_roundness_perimeter_ratio.append(s.str());
        filepath_vertical_box_flattening.append(s.str());

        plot.draw_correlation("Grain box flattening", "Grain roundness", grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_box_flattening,
            grain_roundness, corr_roundness_box_flattening, filepath_roundness_box_flattening);
        plot.draw_correlation("Vertical grain flattening", "Grain roundness", grain_areas, grain_box_flattening, grain_ellipse_flattening,
            grain_area_flattening, grain_roundness, corr_roundness_vertical_flattening, filepath_roundness_vertical_flattening);
        plot.draw_correlation("Grain box flattening", "Grain ellipse flattening", grain_areas, grain_box_flattening, grain_ellipse_flattening,
            grain_box_flattening, grain_ellipse_flattening, corr_ellipse_box_flattening, filepath_ellipse_box_flattening);
        plot.draw_correlation("Grain box width [mm]", "Grain ellipse long axis [mm]", grain_areas, grain_box_flattening, grain_ellipse_flattening,
            grain_box_width, grain_ellipse_long_axis, corr_ellipse_long_box_width, filepath_ellipse_long_box_width);
        plot.draw_correlation("Absolute grain orientation angle [degree]", "Vertical grain flattening", grain_areas, grain_box_flattening,
            grain_ellipse_flattening, grain_ellipse_long_axis_angle, grain_area_flattening, corr_flattening_abs_ellipse_angle,
            filepath_ellipse_abs_angle_flattening, 1, 0);
        plot.draw_correlation("Grain perimeter ratio", "Grain roundness", grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_perimeter_ratio,
            grain_roundness, corr_roundness_perimeter_ratio, filepath_roundness_perimeter_ratio);
        plot.draw_correlation("Grain box flattening", "Vertical grain flattening", grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_box_flattening,
            grain_area_flattening, corr_vertical_box_flattening, filepath_vertical_box_flattening);
    }
    
    //Analyze parameters
    else if(mode == 1 || mode == 3)
    {
        std::vector< std::vector<float> > paramtr;
        std::string paramtr_name;
        std::vector<element_entry> entry;
        int plot_element=0;
        int plot_type=0;

        //Parameter: Grain Size
        if(combo1 == 0)
        {        
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::vector<float> float_vector;
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);        

                size_t & nr_areas = grainBoundNet.nr_new_areas;    
                long * & grain_area_size = grainBoundNet.grain_area_size;

                grain_size_values.push_back(float_vector);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);

                for(int area=0; area<nr_areas; area++) 
                {
                    grain_size_values.back().push_back((float)grain_area_size[area]/area_scaling());
                }
            }                                   

            paramtr = grain_size_values;
            paramtr_name = "grain_size";
            std::cout << "Parameter: grain_size" << std::endl;
        }
        
        //Parameter: Grain Roundness
        else if(combo1 == 1)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);
                    
                grain_roundness.push_back(para.grain_roundness);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);
            }      
                           
            paramtr = grain_roundness;
            paramtr_name = "grain_roundness";
            std::cout << "Parameter: grain roundness" << std::endl;
        }

        //Parameter: Grain Box Flattening
        else if(combo1 == 2)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);
                                       
                grain_box_flattening.push_back(para.grain_box_flattening);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);
            }                                     

            paramtr = grain_box_flattening;
            paramtr_name = "grain_box_flattening";
            std::cout << "Parameter: grain box flattening" << std::endl;
            plot_type = 1;
        }
        
        //Parameter: Grain Box Width
        else if(combo1 == 3)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);
                                       
                grain_box_width.push_back(para.grain_box_width);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);
            }  
          
            paramtr = grain_box_width;
            paramtr_name = "grain_box_width";
            std::cout << "Parameter: grain box width" << std:: endl;
            plot_type = 1;
        }

        //Parameter: Grain Box Height
        else if(combo1 == 4)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);
                                       
                grain_box_height.push_back(para.grain_box_height);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);
            }  

            paramtr = grain_box_height;
            paramtr_name = "grain_box_height"; 
            std::cout << "Parameter: grain box width" << std::endl;
            plot_type = 1;            
        }

        //Parameter: Grain Ellipse Flattening
        else if(combo1 == 5)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);
                                       
                grain_ellipse_flattening.push_back(para.ellipse_flattening);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);
            }  
    
            paramtr = grain_ellipse_flattening;
            paramtr_name = "grain_ellipse_flattening";
            std::cout << "Parameter: grain ellipse flattening" << std::endl;
            plot_type = 2;
        }         

        //Parameter: Grain Ellipse Long Axis
        else if(combo1 == 6)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);
                                       
                grain_ellipse_long_axis.push_back(para.ellipse_long_axis);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);
            }  

            paramtr = grain_ellipse_long_axis;
            paramtr_name = "grain_ellipse_long_axis_angle";
            std::cout << "Parameter: grain ellipse angle" << std::endl;
            plot_type = 2;
        }

        //Parameter: Grain Ellipse Long Axis Angle
        else if(combo1 == 7)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);
                                       
                grain_ellipse_long_axis_angle.push_back(para.ellipse_long_axis_angle);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);
            }

            paramtr = grain_ellipse_long_axis_angle;
            paramtr_name = "grain_ellipse_long_axis_angle";
            std::cout << "Parameter: grain ellipse angle" << std::endl;
            plot_type = 2;
        }

        //Parameter: Grain Area Flattening
        else if(combo1 == 8)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::vector<float> float_vector;
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);
                
                grain_area_flattening.push_back(float_vector);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);

                //Calculate grain area flattening
                grain_area_flattening.back().resize(para.grain_areas.back()+1);
                for(int area=0; area<para.grain_areas.size(); area++)
                {
                    grain_area_flattening.back()[para.grain_areas[area]]=
                        para.grain_area_width[para.grain_areas[area]]/para.grain_area_height[para.grain_areas[area]];
                }                                          
            }
            
            paramtr = grain_area_flattening;
            paramtr_name = "grain_area_flattening";
            std::cout << "Parameter: grain area flattening" << std::endl;
            plot_type = 4;
        }
        
        //Parameter: Grain Area Width
        else if(combo1 == 9)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);
                                       
                grain_area_width.push_back(para.grain_area_width);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);
            }
            
            paramtr = grain_area_width;
            paramtr_name = "grain_area_width";
            std::cout << "Parameter: grain area width" << std::endl;
            plot_type = 4;
        }

        //Parameter: Grain Area Height
        else if(combo1 == 10)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);
                                       
                grain_area_height.push_back(para.grain_area_height);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);
            }

            paramtr = grain_area_height;
            paramtr_name = "grain_area_height";
            std::cout << "Parameter: grain area height" << std::endl;
            plot_type = 4;
        }

        //Parameter: Grain Arc Number
        else if(combo1 == 11)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);
                                       
                grain_arc_number.push_back(para.grain_arc_number);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);
            }        

            paramtr = grain_arc_number;
            paramtr_name = "grain_boundaries";
            std::cout << "Parameter: number grain boundaries" << std::endl;
        }
        
        //Parameter: Grain Neighbors
        else if(combo1 == 12)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);
                                       
                grain_neighbors.push_back(para.grain_neighbors);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);
            }    

            paramtr = grain_neighbors;
            paramtr_name = "grain_neighbors";
            std::cout << "Parameter: number grain neighbors" << std::endl;
        }   

        //Parameter: Grain Longest Grain Boundary
        else if(combo1 == 13)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);
                                       
                grain_longest_arc_length.push_back(para.grain_longest_arc_length);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);
            }

            paramtr = grain_longest_arc_length;
            paramtr_name = "grain_longest_grain_boundary";
            std::cout << "Parameter: grain longest grain boundary" << std::endl;
        }

        //Parameter: Grain Boundary Length
        else if(combo1 == 14)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::vector<float> float_vector;
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                std::string filepath_to_ws_image=path_to_ws_image;
                if (correct_suffix) filepath_to_ws_image.append(suffix1);
                filepath_to_ws_image.append(get_filename(filename));
                filepath_to_ws_image.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                seg segment(false);
                segment.load_cgp_data_structure(filepath_to_ws_image);

                marray::Marray<unsigned int> & two_boundings = segment.two_boundings;
                std::vector< std::vector<point> > & arcs = segment.arcs;

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                long * & grain_area_size = grainBoundNet.grain_area_size;
                long * & bubble_area_size = grainBoundNet.bubble_area_size;
                std::vector<size_t> & region_labels = grainBoundNet.region_labels;
                std::vector<bool> & grain_arc = grainBoundNet.grain_arc;
                grain_arc.resize(arcs.size(),false);

                //adjust twoCellNeighbors
                marray::Marray<unsigned int> twoCellNeighbors;
                size_t size[] = {arcs.size(),2};
                twoCellNeighbors.resize(size,size+2);
               
                for(size_t arc=0;arc<arcs.size();arc++)
                {                    
                    //FIND THE INDEX OF THE NEIGBOUR REGIONS (2)
                    size_t n0=two_boundings(arc,0);
                    size_t n1=two_boundings(arc,1);                   
                    twoCellNeighbors(arc,0)=region_labels[n0-1];
                    twoCellNeighbors(arc,1)=region_labels[n1-1];
                }

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");
                
                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                //seg attribute is required to resize param attribute
                para.grain_junction.resize(segment.junctions.size(),false);
                para.load_extracted_parameters(filepath_parameters);               

                std::vector< std::vector<int> > & grain_boundary_index = para.grain_boundary_index;
                grain_boundary_length.push_back(float_vector);

                //Calculate grain boundary length
                grain_boundary_length.back().resize(grain_boundary_index.size(), 0.0f);
                for (int boundary=0; boundary<grain_boundary_index.size(); boundary++)
                {
                    bool inner_boundary = true;

                    int arc_index=fabs(grain_boundary_index[boundary][0])-1;

                    if ((grain_area_size[twoCellNeighbors(arc_index,1)-1]==0 && bubble_area_size[twoCellNeighbors(arc_index,1)-1]==0) ||
                        (grain_area_size[twoCellNeighbors(arc_index,0)-1]==0 && bubble_area_size[twoCellNeighbors(arc_index,0)-1]==0))
                        inner_boundary=false;

                    if (!inner_boundary) grain_boundary_length.back()[boundary]=0.0f;
                    else for(int arc=0; arc<grain_boundary_index[boundary].size(); arc++)
                        grain_boundary_length.back()[boundary]+=arcs[fabs(grain_boundary_index[boundary][arc])-1].size()/length_scaling();
                }                                
            }

            paramtr = grain_boundary_length;
            paramtr_name = "grain_boundary_length";
            std::cout << "Parameter: grain boundary length" << std::endl;
            plot_element = 1;
        }

        //Parameter: Number of Turning Points
        else if(combo1 == 15)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::vector<float> float_vector;
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                std::string filepath_to_ws_image=path_to_ws_image;
                if (correct_suffix) filepath_to_ws_image.append(suffix1);
                filepath_to_ws_image.append(get_filename(filename));
                filepath_to_ws_image.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                seg segment(false);
                segment.load_cgp_data_structure(filepath_to_ws_image);

                marray::Marray<unsigned int> & two_boundings = segment.two_boundings;
                std::vector< std::vector<point> > & arcs = segment.arcs;

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);

                long * & grain_area_size = grainBoundNet.grain_area_size;
                long * & bubble_area_size = grainBoundNet.bubble_area_size;
                std::vector<size_t> & region_labels = grainBoundNet.region_labels;
                std::vector<bool> & grain_arc = grainBoundNet.grain_arc;
                grain_arc.resize(arcs.size(),false);

                //adjust twoCellNeighbors
                marray::Marray<unsigned int> twoCellNeighbors;
                size_t size[] = {arcs.size(),2};
                twoCellNeighbors.resize(size,size+2);
               
                for(size_t arc=0;arc<arcs.size();arc++)
                {                    
                    //FIND THE INDEX OF THE NEIGBOUR REGIONS (2)
                    size_t n0=two_boundings(arc,0);
                    size_t n1=two_boundings(arc,1);                   
                    twoCellNeighbors(arc,0)=region_labels[n0-1];
                    twoCellNeighbors(arc,1)=region_labels[n1-1];
                }

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");
                
                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);   

                std::vector< std::vector<int> > & grain_boundary_index = para.grain_boundary_index;
                grain_boundary_length.push_back(float_vector);
                turning_point.push_back(para.turning_point);

                //Calculate grain boundary length
                grain_boundary_length.back().resize(grain_boundary_index.size(), 0.0f);
                for (int boundary=0; boundary<grain_boundary_index.size(); boundary++)
                {
                    bool inner_boundary = true;

                    int arc_index=fabs(grain_boundary_index[boundary][0])-1;

                    if ((grain_area_size[twoCellNeighbors(arc_index,1)-1]==0 && bubble_area_size[twoCellNeighbors(arc_index,1)-1]==0) ||
                        (grain_area_size[twoCellNeighbors(arc_index,0)-1]==0 && bubble_area_size[twoCellNeighbors(arc_index,0)-1]==0))
                        inner_boundary=false;

                    if (!inner_boundary) grain_boundary_length.back()[boundary]=0.0f;
                    else for(int arc=0; arc<grain_boundary_index[boundary].size(); arc++)
                        grain_boundary_length.back()[boundary]+=arcs[fabs(grain_boundary_index[boundary][arc])-1].size()/length_scaling();
                }  
            }

            paramtr = turning_point;
            paramtr_name = "turning_point";
            std::cout << "Parameter: number of turning points" << std::endl;
            plot_element = 1;
        }

        //Parameter: Standard Deviation Dihedral Angles
        else if(combo1 == 16)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::vector<float> float_vector;
                std::string filename = filenames[i];
              
                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if (minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                std::string filepath_to_ws_image=path_to_ws_image;
                if (correct_suffix) filepath_to_ws_image.append(suffix1);
                filepath_to_ws_image.append(get_filename(filename));
                filepath_to_ws_image.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                seg segment(false);
                segment.load_cgp_data_structure(filepath_to_ws_image);

                std::vector< std::vector<point> > & arcs = segment.arcs;
                std::vector<point> & junctions = segment.junctions;

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);

                grain_junctions.push_back(grainBoundNet.grain_junctions);
                std::vector<bool> & grain_arc = grainBoundNet.grain_arc;
                grain_arc.resize(arcs.size(),false);

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions=grainBoundNet.grain_junctions;
                //seg attribute is required to resize param attribute
                para.grain_junction.resize(segment.junctions.size(),false);
                para.load_extracted_parameters(filepath_parameters);

                std::vector< std::vector<int> > & grain_boundary_index = para.grain_boundary_index;
                std::vector< std::vector<int> > & grain_junction_index = para.grain_junction_index;
                stddev_dihedral_angles.push_back(float_vector);

                //fill grain boundary pixels
                std::vector< std::vector<point> > grain_boundary_pixels(grain_boundary_index.size());

                for (int boundary=0; boundary<grain_boundary_index.size(); boundary++)
                {
                    for(int arc=0; arc<grain_boundary_index[boundary].size(); arc++)
                    {
                        if (grain_boundary_index[boundary][arc]>0)
                        {
                            for (int p=0; p<arcs[fabs(grain_boundary_index[boundary][arc])-1].size(); p++)
                            {
                                int x=arcs[fabs(grain_boundary_index[boundary][arc])-1][p].x;
                                int y=arcs[fabs(grain_boundary_index[boundary][arc])-1][p].y;

                                point po;
                                po.x=x;
                                po.y=y;

                                grain_boundary_pixels[boundary].push_back(po);
                            }
                        }
                        else
                        {
                            for (int p=arcs[fabs(grain_boundary_index[boundary][arc])-1].size()-1; p>=0; p--)
                            {
                                int x=arcs[fabs(grain_boundary_index[boundary][arc])-1][p].x;
                                int y=arcs[fabs(grain_boundary_index[boundary][arc])-1][p].y;

                                point po;
                                po.x=x;
                                po.y=y;

                                grain_boundary_pixels[boundary].push_back(po);
                            }
                        }
                    }
                }

                for(int j=0; j<grain_junctions.back().size(); j++)
                {
                    for (int i=0; i<grain_junction_index[grain_junctions.back()[j]].size(); i++)
                    {
                        //check which side of grain boundary fits to junction
                        point front=grain_boundary_pixels[grain_junction_index[grain_junctions.back()[j]][i]][0];
                        point back=grain_boundary_pixels[grain_junction_index[grain_junctions.back()[j]][i]].back();

                        if((fabs(front.x-junctions[grain_junctions.back()[j]].x)+fabs(front.y-junctions[grain_junctions.back()[j]].y))<2)
                        {
                            grain_junction_index[grain_junctions.back()[j]][i]=grain_junction_index[grain_junctions.back()[j]][i]+1;
                        }
                        else if((fabs(back.x-junctions[grain_junctions.back()[j]].x)+fabs(back.y-junctions[grain_junctions.back()[j]].y))<2)
                        {
                            grain_junction_index[grain_junctions.back()[j]][i]=-grain_junction_index[grain_junctions.back()[j]][i]-1;
                        }
                        else
                        {
                            std::cout<<"Error: Junction do not fit to boundaries!"<<std::endl;
                            exit(-1);
                        }
                    }
                }

                //calculate standard deviation dihedral angles
                for (int junction=0; junction<grain_junctions.back().size(); junction++)
                {
                    int pixels_average=15;
                    int * boundary_pos = new int[grain_junction_index[grain_junctions.back()[junction]].size()];
                    point * average_pos = new point [grain_junction_index[grain_junctions.back()[junction]].size()];

                    for (int i=0; i<grain_junction_index[grain_junctions.back()[junction]].size(); i++)
                    {
                        boundary_pos[i]=std::min(pixels_average,(int)grain_boundary_pixels[
                            fabs(grain_junction_index[grain_junctions.back()[junction]][i])-1].size());
                    }

                    //compute the "center" of pixels
                    for (int i=0; i<grain_junction_index[grain_junctions.back()[junction]].size(); i++)
                    {
                        int number=0;
                        int sum_x=0;
                        int sum_y=0;

                        for(int p=0; p<boundary_pos[i]; p++)
                        {
                            int index;
                            if(grain_junction_index[grain_junctions.back()[junction]][i]>0) index=p;
                            else index=grain_boundary_pixels[fabs(grain_junction_index[grain_junctions.back()[junction]][i])-1].size()-1-p;

                            sum_x+=grain_boundary_pixels[fabs(grain_junction_index[grain_junctions.back()[junction]][i])-1][index].x;
                            sum_y+=grain_boundary_pixels[fabs(grain_junction_index[grain_junctions.back()[junction]][i])-1][index].y;
                            number++;
                        }

                        average_pos[i].x=sum_x/number;
                        average_pos[i].y=sum_y/number;
                    }

                    std::vector<float> angles;
                    float mean, standard_deviation;

                    for (int i=0; i<3; i++)
                    {
                        float angle_front=atan2((average_pos[i].y-junctions[grain_junctions.back()[junction]].y),
                                                (average_pos[i].x-junctions[grain_junctions.back()[junction]].x));
                        float angle_back=atan2((average_pos[(i+1)%3].y-junctions[grain_junctions.back()[junction]].y),
                                               (average_pos[(i+1)%3].x-junctions[grain_junctions.back()[junction]].x));
                        float angle_other=atan2((average_pos[(i+2)%3].y-junctions[grain_junctions.back()[junction]].y),
                                                (average_pos[(i+2)%3].x-junctions[grain_junctions.back()[junction]].x));

                        //if angle of other arc is in between, add 2pi to smaller angle
                        if(angle_front<angle_other && angle_other<angle_back) angle_front+=2.0f*PI;
                        if(angle_back<angle_other && angle_other<angle_front) angle_back+=2.0f*PI;

                        if (angle_front>angle_back) angles.push_back(180.0f*(angle_front-angle_back)/PI);
                        else angles.push_back(180.0f*(angle_back-angle_front)/PI);
                    }

                    calculate_mean_standard_deviation(angles, mean, standard_deviation);
                    stddev_dihedral_angles.back().push_back(standard_deviation);

                    delete average_pos;
                    delete boundary_pos;
                }
            }
            
            paramtr = stddev_dihedral_angles;
            paramtr_name = "stddev_dihedral_angles";
            std::cout << "Parameter: standard deviation dihedral angles" << std::endl;
            plot_element = 2;
        }

        //Parameter: Grain Perimeter Ratio
        else if(combo1 == 17)
        {
            //Loop over every image in the filename vector
            for(int i = 0; i < filenames.size(); i++)
            {
                std::string filename = filenames[i];

                //extract suffix
                std::string suffix1=get_path(filename);

                if (suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }

                std::string param_file_name=get_filename(param_file);
                if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if (correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);  

                std::string filepath_parameters=path_parameters;
                if (correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters);
                    
                grain_perimeter_ratio.push_back(para.grain_perimeter_ratio);
                grain_area_center_mass.push_back(grainBoundNet.grain_area_center_mass);
            }      
                           
            paramtr = grain_perimeter_ratio;
            paramtr_name = "grain_perimeter_ratio";
            std::cout << "Parameter: grain perimeter ratio" << std::endl;
            plot_type = 5;
        }
        else
        {
            std::cout << "Error: Selected parameter is not defined!" << std::endl;
            exit(-1);
        }

        //Grain areas
        if(plot_element == 0)
        {
            for(int j = 0; j < filenames.size(); j++)
            {
                std::string filename = filenames[j];  
                
                //Extract suffix
                std::string suffix1=get_path(filename);

                if(suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    suffix1.append("/");
                }
                
                std::string param_file_name=get_filename(param_file);
                if(param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                else param_file_name = "";

                std::string filepath_new_classification=path_rf_predictions;
                if(correct_suffix) filepath_new_classification.append(suffix1);
                filepath_new_classification.append(param_file_name.c_str());

                std::stringstream s;
                if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
                else s << "." << minimal_grain_size;
                std::string suffix2=s.str();

                filepath_new_classification.append(get_filename(filename));
                if(minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
                filepath_new_classification.append(".h5");

                //IMPORT RESULTS FROM HDF5 file
                gbn grainBoundNet;
                grainBoundNet.load_final_structure(filepath_new_classification);

                std::string filepath_parameters=path_parameters;
                if(correct_suffix) filepath_parameters.append(suffix1);
                filepath_parameters.append(get_filename(filename));
                filepath_parameters.append(suffix2.c_str());
                filepath_parameters.append(".h5");

                param para;

                //gbn attribute is required as param attribute
                para.grain_junctions = grainBoundNet.grain_junctions;
                para.load_extracted_parameters(filepath_parameters); 

                grain_areas.push_back(para.grain_areas);

                //Prevent double loading (might cause errors otherwise)
                if(combo1 != 2)
                {
                    grain_box_flattening.push_back(para.grain_box_flattening);
                }
                if(combo1 != 5)
                {   
                    grain_ellipse_flattening.push_back(para.ellipse_flattening);
                } 
                if(combo1 != 11)
                {
                    grain_arc_number.push_back(para.grain_arc_number); 
                }
            }
            
            for(int i=0; i<grain_areas.size(); i++)
            {                
                for(int area=0; area<grain_areas[i].size(); area++)
                {                    
                    if(combo1>1 && combo1<8)
                    {                        
                        //Ellipse not degenerated and axes found correctly
                        if(grain_box_flattening[i][grain_areas[i][area]]>0.0f && grain_ellipse_flattening[i][grain_areas[i][area]]>0.0f)
                        {
                            element_entry new_grain;
                            
                            new_grain.value=paramtr[i][grain_areas[i][area]];
                            new_grain.image=i;
                            new_grain.element=area;

                            entry.push_back(new_grain);
                        }
                    }                   
                    else if(combo1>10 && combo1<14)
                    {                        
                        //Only inner entry
                        if(grain_arc_number[i][grain_areas[i][area]]>0.0f)
                        {
                            element_entry new_grain;
                            
                            new_grain.value=paramtr[i][grain_areas[i][area]];
                            new_grain.image=i;
                            new_grain.element=area;

                            entry.push_back(new_grain);
                        }
                    }
                    else 
                    {
                        element_entry new_grain;
                        
                        new_grain.value=paramtr[i][grain_areas[i][area]];
                        new_grain.image=i;
                        new_grain.element=area;

                        entry.push_back(new_grain);
                    }
                }
            }
        } 

        //Grain boundaries    
        else if(plot_element == 1)
        {           
            for(int i=0; i<grain_boundary_length.size(); i++)
            {
                for(int boundary=0; boundary<grain_boundary_length[i].size(); boundary++)
                {
                    //Only inner entry
                    if (grain_boundary_length[i][boundary]>0.0f)
                    {
                        element_entry new_grain;
                        
                        new_grain.value=paramtr[i][boundary];
                        new_grain.image=i;
                        new_grain.element=boundary;

                        entry.push_back(new_grain);
                    }
                }
            }
        }

        //Grain triple junctions
        else if(plot_element == 2)
        {            
            for(int i=0; i<grain_junctions.size(); i++)
            {
                for(int junction=0; junction<grain_junctions[i].size(); junction++)
                {
                    element_entry new_grain;
                    
                    new_grain.value=paramtr[i][junction];
                    new_grain.image=i;
                    new_grain.element=junction;

                    entry.push_back(new_grain);
                }
            }
        }
        
        //Sort elements
        std::sort(entry.begin(),entry.end(),entry_cmp);
        
        //Lowest values
        if(combo2 == 0)
        {
            std::cout<<"lowest values: "<<std::endl;

            for(int j=entry.size()-1; j>=0 && j>entry.size()-1-spin; j--)
            {
                std::stringstream s;

                if(plot_element==0) s<<get_filename(filenames[entry[j].image].c_str())<<", grain "<<
                    grain_areas[entry[j].image][entry[j].element]<<": "<<paramtr[entry[j].image][grain_areas[entry[j].image][entry[j].element]];
                else if(plot_element==1) s<<get_filename(filenames[entry[j].image].c_str())<<", grain boundary "<<entry[j].element<<": "<<
                    paramtr[entry[j].image][entry[j].element];
                else if(plot_element==2) s<<get_filename(filenames[entry[j].image].c_str())<<", grain triple junction "<<
                    grain_junctions[entry[j].image][entry[j].element]<<": "<<paramtr[entry[j].image][entry[j].element];

                std::string output_string = s.str();
                std::cout<<output_string<<std::endl;

                std::string filepath_to_feature_file=filenames[entry[j].image];
                filepath_to_feature_file.append(".bin");

                //extract suffix
                std::string path_ws_image=path_to_ws_image;
                std::string path_predictions=path_rf_predictions;
                std::string path_extract_params=path_parameters;
                std::string suffix=get_path(filenames[entry[j].image].c_str());

                if(suffix.size()>0)
                {
                    suffix.resize(suffix.size()-1);
                    suffix=get_filename(suffix.c_str());
                    suffix.append("/");
                }

                if(correct_suffix)
                {
                    path_ws_image.append(suffix);
                    path_predictions.append(suffix);
                    path_extract_params.append(suffix);
                }

                std::string analyze_plot_path=path_plots;
                analyze_plot_path.append(paramtr_name);
                analyze_plot_path.append("_lowest/");

                if(!directory_exists(analyze_plot_path.c_str())) mkdir(analyze_plot_path.c_str(),0777);
                
                //plot lowest values
                if (mode==1)
                {
                    if (plot_element==0) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions, analyze_plot_path,
                        path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                        grain_areas[entry[j].image][entry[j].element], true, plot_type);
                    else if (plot_element==1) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions,
                        analyze_plot_path, path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                        entry[j].element, true, plot_type);
                    else if (plot_element==2) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions,
                        analyze_plot_path, path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                        grain_junctions[entry[j].image][entry[j].element], true, plot_type);
                }
                //lowest values to file
                else
                {
                    std::string filepath_list=analyze_plot_path;
                    filepath_list.append("list.txt");
                    std::ofstream file(filepath_list.c_str(), std::ios_base::out | std::ios_base::app);
                    if (plot_element == 0) file << output_string << " at position (" <<
                        grain_area_center_mass[entry[j].image][grain_areas[entry[j].image][entry[j].element]].x <<"," << 
                        grain_area_center_mass[entry[j].image][grain_areas[entry[j].image][entry[j].element]].y << ")" << "\n";
                    else file << output_string << "\n";
                    file.close();
                }
            }
        }
    
        //Median values
        else if(combo2 == 1)
        {
            std::cout<<"median values: "<<std::endl;
            int start=std::max(0,((int)entry.size()-1+spin)/2);
            for (int j=start; j>=0 && j>start-spin; j--)
            {
                std::stringstream s;

                if (plot_element==0) s<<get_filename(filenames[entry[j].image].c_str())<<", grain "<<
                    grain_areas[entry[j].image][entry[j].element]<<": "<<paramtr[entry[j].image][grain_areas[entry[j].image][entry[j].element]];
                else if (plot_element==1) s<<get_filename(filenames[entry[j].image].c_str())<<", grain boundary "<<entry[j].element<<": "<<
                    paramtr[entry[j].image][entry[j].element];
                else if (plot_element==2) s<<get_filename(filenames[entry[j].image].c_str())<<", grain triple junction "<<
                    grain_junctions[entry[j].image][entry[j].element]<<": "<<paramtr[entry[j].image][entry[j].element];

                std::string output_string = s.str();
                std::cout<<output_string<<std::endl;

                std::string filepath_to_feature_file=filenames[entry[j].image];
                filepath_to_feature_file.append(".bin");

                //extract suffix
                std::string path_ws_image=path_to_ws_image;
                std::string path_predictions=path_rf_predictions;
                std::string path_extract_params=path_parameters;
                std::string suffix=get_path(filenames[entry[j].image].c_str());

                if (suffix.size()>0)
                {
                    suffix.resize(suffix.size()-1);
                    suffix=get_filename(suffix.c_str());
                    suffix.append("/");
                }

                if (correct_suffix)
                {
                    path_ws_image.append(suffix);
                    path_predictions.append(suffix);
                    path_extract_params.append(suffix);
                }

                std::string analyze_plot_path=path_plots;
                analyze_plot_path.append(paramtr_name);
                analyze_plot_path.append("_median/");

                if (!directory_exists(analyze_plot_path.c_str())) mkdir(analyze_plot_path.c_str(),0777);

                //plot median values
                if (mode==1)
                {
                    if (plot_element==0) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions, analyze_plot_path,
                        path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                        grain_areas[entry[j].image][entry[j].element], true, plot_type);
                    else if (plot_element==1) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions,
                        analyze_plot_path, path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                        entry[j].element, true, plot_type);
                    else if (plot_element==2) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions,
                        analyze_plot_path, path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                        grain_junctions[entry[j].image][entry[j].element], true, plot_type);
                }
                //median values to file
                else
                {
                    std::string filepath_list=analyze_plot_path;
                    filepath_list.append("list.txt");
                    std::ofstream file(filepath_list.c_str(), std::ios_base::out | std::ios_base::app);
                    if (plot_element == 0) file << output_string << " at position (" <<
                        grain_area_center_mass[entry[j].image][grain_areas[entry[j].image][entry[j].element]].x <<"," << 
                        grain_area_center_mass[entry[j].image][grain_areas[entry[j].image][entry[j].element]].y << ")" << "\n";
                    else file << output_string << "\n";
                    file.close();
                }
            }
        }

        //Highest values
        else if(combo2 == 2)
        {
            std::cout<<"highest values: "<<std::endl;
            for (int j=0; j<entry.size() && j<spin; j++)
            {
                std::stringstream s;

                if (plot_element==0) s<<get_filename(filenames[entry[j].image].c_str())<<", grain "<<
                    grain_areas[entry[j].image][entry[j].element]<<": "<<paramtr[entry[j].image][grain_areas[entry[j].image][entry[j].element]];
                else if (plot_element==1) s<<get_filename(filenames[entry[j].image].c_str())<<", grain boundary "<<entry[j].element<<": "<<
                    paramtr[entry[j].image][entry[j].element];
                else if (plot_element==2) s<<get_filename(filenames[entry[j].image].c_str())<<", grain triple junction "<<
                    grain_junctions[entry[j].image][entry[j].element]<<": "<<paramtr[entry[j].image][entry[j].element];

                std::string output_string = s.str();
                std::cout<<output_string<<std::endl;

                std::string filepath_to_feature_file=filenames[entry[j].image];
                filepath_to_feature_file.append(".bin");

                //extract suffix
                std::string path_ws_image=path_to_ws_image;
                std::string path_predictions=path_rf_predictions;
                std::string path_extract_params=path_parameters;
                std::string suffix=get_path(filenames[entry[j].image].c_str());

                if (suffix.size()>0)
                {
                    suffix.resize(suffix.size()-1);
                    suffix=get_filename(suffix.c_str());
                    suffix.append("/");
                }

                if (correct_suffix)
                {
                    path_ws_image.append(suffix);
                    path_predictions.append(suffix);
                    path_extract_params.append(suffix);
                }

                std::string analyze_plot_path=path_plots;
                analyze_plot_path.append(paramtr_name);
                analyze_plot_path.append("_highest/");

                if (!directory_exists(analyze_plot_path.c_str())) mkdir(analyze_plot_path.c_str(),0777);
                
                //plot highest values
                if (mode==1)
                {
                    if (plot_element==0) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions, analyze_plot_path,
                        path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                        grain_areas[entry[j].image][entry[j].element], true, plot_type);
                    else if (plot_element==1) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions,
                        analyze_plot_path, path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                        entry[j].element, true, plot_type);
                    else if (plot_element==2) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions,
                        analyze_plot_path, path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                        grain_junctions[entry[j].image][entry[j].element], true, plot_type);
                }
                //highest values to file
                else
                {
                    std::string filepath_list=analyze_plot_path;
                    filepath_list.append("list.txt");
                    std::ofstream file(filepath_list.c_str(), std::ios_base::out | std::ios_base::app);
                    if (plot_element == 0) file << output_string << " at position (" <<
                        grain_area_center_mass[entry[j].image][grain_areas[entry[j].image][entry[j].element]].x <<"," << 
                        grain_area_center_mass[entry[j].image][grain_areas[entry[j].image][entry[j].element]].y << ")" << "\n";
                    else file << output_string << "\n";
                    file.close();
                }
            }
        }
    }

    //Analyze correlations
    else if(mode == 2) 
    {
        //Calculate correlations
        correlation corr_roundness_box_flattening;
        correlation corr_roundness_vertical_flattening;
        correlation corr_ellipse_box_flattening;
        correlation corr_ellipse_long_box_width;
        correlation corr_flattening_abs_ellipse_angle;
        correlation corr_roundness_perimeter_ratio;
        correlation corr_vertical_box_flattening;

        //Loop over every image specified in the list file
        for(int i = 0; i < filenames.size(); i++)
        {
            std::vector<float> float_vector;
            std::string filename = filenames[i];

            //extract suffix
            std::string suffix1=get_path(filename);

            if (suffix1.size()>0)
            {
                suffix1.resize(suffix1.size()-1);
                suffix1=get_filename(suffix1.c_str());
                suffix1.append("/");
            }
            
            std::string param_file_name=get_filename(param_file);
            if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
            else param_file_name = "";

            std::string filepath_new_classification=path_rf_predictions;
            if (correct_suffix) filepath_new_classification.append(suffix1);
            filepath_new_classification.append(param_file_name.c_str());

            std::stringstream s;
            if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
            else s << "." << minimal_grain_size;
            std::string suffix2=s.str();

            filepath_new_classification.append(get_filename(filename));
            if (minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix2.c_str());
            filepath_new_classification.append(".h5");

            //IMPORT RESULTS FROM HDF5 file
            gbn grainBoundNet;
            grainBoundNet.load_final_structure(filepath_new_classification);                

            std::string filepath_parameters=path_parameters;
            if (correct_suffix) filepath_parameters.append(suffix1);
            filepath_parameters.append(get_filename(filename));
            filepath_parameters.append(suffix2.c_str());
            filepath_parameters.append(".h5");
            
            param para;

            //gbn attribute is required as param attribute
            para.grain_junctions=grainBoundNet.grain_junctions;
            para.load_extracted_parameters(filepath_parameters);            

            grain_area_flattening.push_back(float_vector);
            grain_areas.push_back(para.grain_areas);
            grain_box_flattening.push_back(para.grain_box_flattening);
            grain_ellipse_flattening.push_back(para.ellipse_flattening);
            grain_roundness.push_back(para.grain_roundness);
            grain_area_width.push_back(para.grain_area_width);
            grain_area_height.push_back(para.grain_area_height);
            grain_box_width.push_back(para.grain_box_width);
            grain_ellipse_long_axis.push_back(para.ellipse_long_axis);
            grain_ellipse_long_axis_angle.push_back(para.ellipse_long_axis_angle);
            grain_perimeter_ratio.push_back(para.grain_perimeter_ratio);

            grain_area_flattening.back().resize(grain_areas.back().back()+1);
            for(int area=0; area<grain_areas.back().size(); area++)
            {
                grain_area_flattening.back()[grain_areas.back()[area]]=
                    grain_area_width.back()[grain_areas.back()[area]]/grain_area_height.back()[grain_areas.back()[area]];
            }
        }

        grain_correlation(grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_box_flattening, grain_roundness,
            corr_roundness_box_flattening);
        grain_correlation(grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_area_flattening, grain_roundness,
            corr_roundness_vertical_flattening);
        grain_correlation(grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_box_flattening, grain_ellipse_flattening, 
            corr_ellipse_box_flattening);
        grain_correlation(grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_box_width, grain_ellipse_long_axis,
            corr_ellipse_long_box_width);
        grain_correlation(grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_ellipse_long_axis_angle, grain_area_flattening, 
            corr_flattening_abs_ellipse_angle, 1, 0);
        grain_correlation(grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_perimeter_ratio, grain_roundness,
            corr_roundness_perimeter_ratio);
        grain_correlation(grain_areas, grain_box_flattening, grain_ellipse_flattening, grain_box_flattening, grain_area_flattening,
            corr_vertical_box_flattening);

        std::vector< std::vector<float> > param1;
        std::vector< std::vector<float> > param2;
        std::string corr_name;
        correlation corr;
        std::vector<element_entry> entry;
        int plot_element = 0;
        int plot_type = 0;

        //Correlation: Ellipse/Box Flattening
        if(combo1 == 0)
        {
            
            param1 = grain_box_flattening;
            param2 = grain_ellipse_flattening;
            corr_name = "corr_ellipse_box_flattening";
            corr = corr_ellipse_box_flattening;
            std::cout<<"Correlation ellipse/box flattening"<<std::endl;
            plot_type=3;
        }
    
        //Correlation: Ellipse Long Axis/Box Width
        else if(combo1==1)
        {
            param1=grain_box_width;
            param2=grain_ellipse_long_axis;
            corr_name="corr_ellipse_long_box_width";
            corr=corr_ellipse_long_box_width;
            std::cout<<"Correlation ellipse long axis/box width"<<std::endl;
            plot_type=3;
        }
        
        //Correlation: Vertical Flattening/Absoulute Grain Orientation Angle
        else if(combo1==2)
        {
            param1=grain_ellipse_long_axis_angle;
            param2=grain_area_flattening;
            corr_name="corr_flattening_abs_ellipse_angle";
            corr=corr_flattening_abs_ellipse_angle;
            std::cout<<"Correlation vertical flattening/absolute grain orientation angle"<<std::endl;
            plot_type=2;
        }

        //Correlation: Roundess/Box Flattening
        else if(combo1==3)
        {
            param1=grain_box_flattening;
            param2=grain_roundness;
            corr_name="corr_roundness_box_flattening";
            corr=corr_roundness_box_flattening;
            std::cout<<"Correlation roundness/box flattening"<<std::endl;
            plot_type=1;
        }

        //Correlation: Roundness/Vertical Flattening
        else if(combo1==4)
        {
            param1=grain_area_flattening;
            param2=grain_roundness;
            corr_name="corr_roundness_vertical_flattening";
            corr=corr_roundness_vertical_flattening;
            std::cout<<"Correlation roundness/vertical flattening"<<std::endl;
            plot_type=4;
        }

        //Correlation: Roundness/Perimeter Ratio
        else if(combo1==5)
        {
            param1=grain_perimeter_ratio;
            param2=grain_roundness;
            corr_name="corr_roundness_perimeter_ratio";
            corr=corr_roundness_perimeter_ratio;
            std::cout<<"Correlation roundness/perimeter ratio"<<std::endl;
            plot_type=5;
        }

        //Correlation: Vertical/Box Flattening
        else if(combo1==6)
        {
            param1=grain_box_flattening;
            param2=grain_area_flattening;
            corr_name="corr_vertical_box_flattening";
            corr=corr_vertical_box_flattening;
            std::cout<<"Correlation vertical/box flattening"<<std::endl;
            plot_type=1;
        }
        else
        {
            std::cout<<"Error: Selected correlation is not defined!"<<std::endl;
            exit(-1);
        }

        for (int i=0; i<grain_areas.size(); i++)
        {
            for (int area=0; area<grain_areas[i].size(); area++)
            {
                //ellipse not degenerated and axes found correctly
                if (grain_box_flattening[i][grain_areas[i][area]]>0.0f && grain_ellipse_flattening[i][grain_areas[i][area]]>0.0f)
                {
                    if (combo1==2) param1[i][grain_areas[i][area]]=fabs(param1[i][grain_areas[i][area]]);                 

                    element_entry new_grain;
                    
                    new_grain.value=param2[i][grain_areas[i][area]]-(corr.m*param1[i][grain_areas[i][area]]+corr.b);
                    new_grain.image=i;
                    new_grain.element=area;

                    entry.push_back(new_grain);
                }
            }
        }

        std::sort(entry.begin(),entry.end(),entry_cmp_abs);

        std::cout<<"Pearson correlation: "<<corr.rho<<", nr grains: "<<corr.nr<<std::endl;

        if (combo2==0)
        {
            std::cout<<"normal case values: "<<std::endl;
            for (int j=entry.size()-1; j>=0 && j>entry.size()-1-spin; j--)
            {
                if (plot_element==0) std::cout<<get_filename(filenames[entry[j].image].c_str())<<", grain "<<grain_areas[entry[j].image][entry[j].element]<<
                    ": ("<<param1[entry[j].image][grain_areas[entry[j].image][entry[j].element]]<<","<<
                    param2[entry[j].image][grain_areas[entry[j].image][entry[j].element]]<<"), residuum: "<<entry[j].value<<std::endl;
                else if (plot_element==1) std::cout<<get_filename(filenames[entry[j].image].c_str())<<", grain boundary "<<entry[j].element<<": ("<<
                    param1[entry[j].image][entry[j].element]<<","<<param2[entry[j].image][entry[j].element]<<"), residuum: "<<entry[j].value<<std::endl;
                else if (plot_element==2) std::cout<<get_filename(filenames[entry[j].image].c_str())<<", grain triple junction "<<
                    grain_junctions[entry[j].image][entry[j].element]<<": ("<<param1[entry[j].image][entry[j].element]<<","<<
                    param2[entry[j].image][entry[j].element]<<"), residuum: "<<entry[j].value<<std::endl;

                std::string filepath_to_feature_file=filenames[entry[j].image];
                filepath_to_feature_file.append(".bin");

                //extract suffix
                std::string path_ws_image=path_to_ws_image;
                std::string path_predictions=path_rf_predictions;
                std::string path_extract_params=path_parameters;
                std::string suffix=get_path(filenames[entry[j].image].c_str());

                if (suffix.size()>0)
                {
                    suffix.resize(suffix.size()-1);
                    suffix=get_filename(suffix.c_str());
                    suffix.append("/");
                }

                if (correct_suffix)
                {
                    path_ws_image.append(suffix);
                    path_predictions.append(suffix);
                    path_extract_params.append(suffix);
                }

                std::string analyze_plot_path=path_plots;
                analyze_plot_path.append(corr_name);
                analyze_plot_path.append("_normal_case/");

                if (!directory_exists(analyze_plot_path.c_str())) mkdir(analyze_plot_path.c_str(),0777);

                if (plot_element==0) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions, analyze_plot_path,
                    path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                    grain_areas[entry[j].image][entry[j].element], true, plot_type);
                else if (plot_element==1) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions,
                    analyze_plot_path, path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                    entry[j].element, true, plot_type);
                else if (plot_element==2) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions,
                    analyze_plot_path, path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                    grain_junctions[entry[j].image][entry[j].element], true, plot_type);
            }
        }
        else if (combo2==1)
        {
            std::cout<<"most outlying values: "<<std::endl;
            for (int j=0; j<entry.size() && j<spin; j++)
            {
                if (plot_element==0) std::cout<<get_filename(filenames[entry[j].image].c_str())<<", grain "<<grain_areas[entry[j].image][entry[j].element]<<
                    ": ("<<param1[entry[j].image][grain_areas[entry[j].image][entry[j].element]]<<","<<
                    param2[entry[j].image][grain_areas[entry[j].image][entry[j].element]]<<"), residuum: "<<entry[j].value<<std::endl;
                else if (plot_element==1) std::cout<<get_filename(filenames[entry[j].image].c_str())<<", grain boundary "<<entry[j].element<<": ("<<
                    param1[entry[j].image][entry[j].element]<<","<<param2[entry[j].image][entry[j].element]<<"), residuum: "<<entry[j].value<<std::endl;
                else if (plot_element==2) std::cout<<get_filename(filenames[entry[j].image].c_str())<<", grain triple junction "<<
                    grain_junctions[entry[j].image][entry[j].element]<<": ("<<param1[entry[j].image][entry[j].element]<<","<<
                    param2[entry[j].image][entry[j].element]<<"), residuum: "<<entry[j].value<<std::endl;

                std::string filepath_to_feature_file=filenames[entry[j].image];
                filepath_to_feature_file.append(".bin");

                //extract suffix
                std::string path_ws_image=path_to_ws_image;
                std::string path_predictions=path_rf_predictions;
                std::string path_extract_params=path_parameters;
                std::string suffix=get_path(filenames[entry[j].image].c_str());

                if (suffix.size()>0)
                {
                    suffix.resize(suffix.size()-1);
                    suffix=get_filename(suffix.c_str());
                    suffix.append("/");
                }

                if (correct_suffix)
                {
                    path_ws_image.append(suffix);
                    path_predictions.append(suffix);
                    path_extract_params.append(suffix);
                }

                std::string analyze_plot_path=path_plots;
                analyze_plot_path.append(corr_name);
                analyze_plot_path.append("_outlying/");

                if (!directory_exists(analyze_plot_path.c_str())) mkdir(analyze_plot_path.c_str(),0777);

                if (plot_element==0) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions, analyze_plot_path,
                    path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                    grain_areas[entry[j].image][entry[j].element], true, plot_type);
                else if (plot_element==1) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions,
                    analyze_plot_path, path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                    entry[j].element, true, plot_type);
                else if (plot_element==2) view(filepath_to_feature_file, get_path(filenames[entry[j].image]), path_ws_image, path_predictions,
                    analyze_plot_path, path_extract_params, param_file, paramFile, minimal_grain_size, minimal_bubble_distance, plot_element,
                    grain_junctions[entry[j].image][entry[j].element], true, plot_type);
            }
        }
        
    }

    list_file.close();
    temp_list_file.close();
}

void new_depth_profiles(std::string filepath_list, std::string path_rf_predictions, std::string param_file, ParameterFile paramFile,
    int minimal_bubble_distance, std::string path_results, bool correct_suffix, int low_grain_size, int high_grain_size,
    int grain_size_step, float depth_bin_width)
{
    if(depth_bin_width<=0.0f)
    {
        std::cout<<"Error: Positive depth bin width required!"<<std::endl;
        return;
    }

    //initialise plplot class
    plplot plot = plplot();

    std::ifstream list_file(filepath_list.c_str());

    if(!list_file)
    {
        std::cout<<"Error: List file could not be found!"<<std::endl;
        return;
    }

    std::ifstream temp_list_file(filepath_list.c_str());
    std::string teststring;
    temp_list_file>>teststring;
    if(teststring.size()==0)
    {
        std::cout<<"Error: List file is empty!"<<std::endl;
        list_file.close();
        temp_list_file.close();
        return;
    }

    Parameter<float> length_scaling;
    length_scaling.assign("", "length_scaling", 193.5f);
    length_scaling.load(paramFile,"config");

    Parameter<float> area_scaling;
    area_scaling.assign("", "area_scaling", 37444.0f);
    area_scaling.load(paramFile,"config");

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

    int grain_size_min=1;
    float grain_bin_width=5.0f;
    float grain_size_start_mu=2.0f;
    float grain_size_start_sigma=0.4f;
    float grain_size_y_max=0.0f;

    float dislocation_density_scaling=length_scaling/16024666666867000.0f;

    for(int minimal_grain_size=low_grain_size; minimal_grain_size<high_grain_size+1; minimal_grain_size+=grain_size_step)
    {
        std::vector<float> grainsize_fit_depths;
        std::vector<float> grainsize_fit_values;
        std::vector<float> grainsize_fit_errors_low;
        std::vector<float> grainsize_fit_errors_high;

        std::vector<float> grainsize_fit_depths_relax;
        std::vector<float> grainsize_fit_values_relax;
        std::vector<float> grainsize_fit_errors_low_relax;
        std::vector<float> grainsize_fit_errors_high_relax;

        std::vector<std::vector<float> > percent_depths(9);
        std::vector<std::vector<float> > grainsize_percent_values(9);
        std::vector<std::vector<float> > grainflattening_percent_values(9);
        std::vector<std::vector<float> > grainround_percent_values(9);
        std::vector<std::vector<float> > grainboxflat_percent_values(9);
        std::vector<std::vector<float> > grainwidth_percent_values(9);
        std::vector<std::vector<float> > grainheight_percent_values(9);
        std::vector<std::vector<float> > grainangle_percent_values(9);
        std::vector<std::vector<float> > grainabsangle_percent_values(10);
        std::vector<std::vector<float> > grainratio_percent_values(9);

        std::vector<std::vector<float> > percent_depths_relax(9);
        std::vector<std::vector<float> > grainsize_percent_values_relax(9);
        std::vector<std::vector<float> > grainflattening_percent_values_relax(9);
        std::vector<std::vector<float> > grainround_percent_values_relax(9);
        std::vector<std::vector<float> > grainboxflat_percent_values_relax(9);
        std::vector<std::vector<float> > grainwidth_percent_values_relax(9);
        std::vector<std::vector<float> > grainheight_percent_values_relax(9);
        std::vector<std::vector<float> > grainangle_percent_values_relax(9);
        std::vector<std::vector<float> > grainabsangle_percent_values_relax(10);
        std::vector<std::vector<float> > grainratio_percent_values_relax(9);

        std::vector<std::vector<float> > grainsize_quantile_depths(10);
        std::vector<std::vector<float> > grainsize_quantile_values(10);

        std::vector<std::vector<float> > grainsize_quantile_depths_relax(10);
        std::vector<std::vector<float> > grainsize_quantile_values_relax(10);

        std::vector<float> grainsize_percent_profile_depths;
        std::vector<float> grainsize_percent_profile_x_values;
        for(int p=0; p<10; p++) grainsize_percent_profile_x_values.push_back((float)(p+1)/10.0f);
        std::vector< std::vector<float> > grainsize_percent_profile_values;

        std::vector<float> grainsize_quantile_profile_depths;
        std::vector<float> grainsize_quantile_profile_x_values;
        for(int q=0; q<10; q++) grainsize_quantile_profile_x_values.push_back((float)(q+1)/100.0f);
        std::vector< std::vector<float> > grainsize_quantile_profile_values;

        std::vector< std::vector<float> > grainsize_combined_depths;
        grainsize_combined_depths.resize(4);
        std::vector< std::vector<float> > grainsize_combined_values;
        grainsize_combined_values.resize(4);

        std::vector<float> grainsize_median_depths;
        std::vector<float> grainsize_median_values;

        std::vector<float> grainsize_median_depths_relax;
        std::vector<float> grainsize_median_values_relax;

        std::vector<float> grainsize_5_largest_depths;
        std::vector<float> grainsize_5_largest_values;
        std::vector<float> grainsize_5_largest_errors;

        std::vector<std::vector<float> > disl_dens_percent_boundaries_depths(10);
        std::vector<std::vector<float> > disl_dens_percent_boundaries_values(10);

        std::vector<float> disl_dens_percent_boundaries_profile_depths;
        std::vector<float> disl_dens_percent_boundaries_profile_x_values;
        for(int p=0; p<10; p++) disl_dens_percent_boundaries_profile_x_values.push_back((float)(p+1)/10.0f);
        std::vector< std::vector<float> > disl_dens_percent_boundaries_profile_values;

        int iterations=1;
        if (minimal_bubble_distance>0) iterations=2;

        for(int iter=0; iter<iterations; iter++)
        {
            //reopen file for second iteration
            if(!list_file.is_open()) list_file.open(filepath_list.c_str());

            //These data structures are saved in depths bins for all images in list file
            std::vector<std::vector<param_entry> > grain_sizes;
            std::vector<std::vector<int> > grain_size_histogram;
            std::vector<long> grain_size_max;
            std::vector<float> grain_size_mean;
            std::vector<float> grain_size_max_x;

            std::vector<std::vector<param_entry> > boundary_sizes;

            std::string old="--";
            int line=1;
            float depth_max=0.0f;

            //loop over images in list file
            while(!list_file.eof())
            {
                std::string filename;
                list_file>>filename;

                if(filename==old || filename=="") break;//avoids reading last line twice

                float depth=0.0f;

                //extract suffix
                std::string suffix1=get_path(filename);

                if (correct_suffix && suffix1.size()>0)
                {
                    suffix1.resize(suffix1.size()-1);
                    suffix1=get_filename(suffix1.c_str());
                    int bagnr=atoi(suffix1.c_str());
                    suffix1.append("/");

                    if(iter==1 && bagnr>2506) break;
                    depth=(float)bagnr*0.55f;
                }
                else
                {
                    depth=get_depth(get_filename(filename));
                }

                if (depth>0.0f)
                {
                    if (depth>depth_max)
                    {
                        depth_max=depth;
                        grain_sizes.resize(1+depth_max/depth_bin_width);
                        grain_size_histogram.resize(1+depth_max/depth_bin_width);
                        grain_size_max.resize(1+depth_max/depth_bin_width,0);
                        grain_size_mean.resize(1+depth_max/depth_bin_width);
                        grain_size_max_x.resize(1+depth_max/depth_bin_width);
                        boundary_sizes.resize(1+depth_max/depth_bin_width);
                    }
                    
                    gbn grainBoundNet;
                    size_t & nr_areas = grainBoundNet.nr_new_areas;
                    long * & grain_area_size = grainBoundNet.grain_area_size;

                    std::string param_file_name=get_filename(param_file);
                    if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
                    else param_file_name = "";

                    std::string filepath_new_classification=path_rf_predictions;
                    if (correct_suffix) filepath_new_classification.append(suffix1);
                    filepath_new_classification.append(param_file_name.c_str());

                    std::stringstream s;
                    if (iter==0 && minimal_grain_size>0) s << "." << minimal_grain_size;
                    else if (iter>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;

                    filepath_new_classification.append(get_filename(filename));
                    filepath_new_classification.append(s.str());
                    filepath_new_classification.append(".h5");

                    //IMPORT RESULTS FROM HDF5 file
                    grainBoundNet.load_grain_sizes(filepath_new_classification);

                    param para;
                    std::vector<int> & grain_areas = para.grain_areas;
                    std::vector<float> & grain_roundness = para.grain_roundness;
                    std::vector<float> & grain_box_flattening = para.grain_box_flattening;
                    std::vector<float> & ellipse_long_axis_angle = para.ellipse_long_axis_angle;
                    std::vector<float> & grain_area_width = para.grain_area_width;
                    std::vector<float> & grain_area_height = para.grain_area_height;
                    std::vector< std::vector<float> > & grain_boundary_curvs = para.grain_boundary_curvs;
                    std::vector<float> & grain_perimeter_ratio = para.grain_perimeter_ratio;

                    std::stringstream s2;
                    if (iter==0) s2 << "." << minimal_grain_size << ".h5";
                    else s2 << "." << minimal_bubble_distance << "." << minimal_grain_size <<".h5";

                    std::string filepath_parameters=path_results;
                    if (correct_suffix) filepath_parameters.append(suffix1);
                    filepath_parameters.append(get_filename(filename));
                    filepath_parameters.append(s2.str());

                    //load grain width and height
                    para.load_selected_extracted_parameters(filepath_parameters);

                    for (int area_index=0; area_index<grain_areas.size(); area_index++)
                    {
                        int area=grain_areas[area_index];

                        if (grain_area_size[area]>grain_size_max[depth/depth_bin_width])
                        {
                            grain_size_histogram[depth/depth_bin_width].resize((int)grain_bin_width*log10(grain_area_size[area])+1,0);
                            grain_size_max[depth/depth_bin_width]=grain_area_size[area];
                        }

                        if (grain_area_size[area]>grain_size_min)
                        {
                            grain_size_histogram[depth/depth_bin_width][(int)grain_bin_width*log10(grain_area_size[area])]++;
                            param_entry entry;
                            entry.size=grain_area_size[area]/area_scaling();
                            entry.flat=grain_area_width[area]/grain_area_height[area];
                            entry.round=grain_roundness[area];
                            entry.boxflat=grain_box_flattening[area];
                            entry.width=grain_area_width[area];
                            entry.height=grain_area_height[area];
                            entry.angle=ellipse_long_axis_angle[area];
                            entry.ratio=grain_perimeter_ratio[area];
                            grain_sizes[depth/depth_bin_width].push_back(entry);
                        }
                    }

                    delete grain_area_size;
                    
                    for (int boundary=0; boundary<grain_boundary_curvs.size(); boundary++)
                    {
                        std::vector<float> curv_values;
                        float mean=0.0f;
                        float standard_deviation;

                        for (int p=0; p<grain_boundary_curvs[boundary].size(); p++)
                        {
                            float curv=fabs(grain_boundary_curvs[boundary][p])/0.75f;//gauging, has to be validated!!

                            if (curv<0.15)//CUT OFF BY CURV 0.15!!
                            {
                                curv_values.push_back(curv);
                            }
                        }

                        if(curv_values.size()>0) calculate_mean_standard_deviation(curv_values, mean, standard_deviation);

                        param_entry entry;
                        entry.size=curv_values.size();
                        entry.flat=mean; //mean curv
                        boundary_sizes[depth/depth_bin_width].push_back(entry);
                    }
                }

                line++;
                old=filename;
            }

            for(int p=0; p<9; p++)
            {
                if (iter==0)
                {
                    percent_depths[p].resize(1+depth_max/depth_bin_width);
                    grainsize_percent_values[p].resize(1+depth_max/depth_bin_width);
                    grainflattening_percent_values[p].resize(1+depth_max/depth_bin_width);
                    grainround_percent_values[p].resize(1+depth_max/depth_bin_width);
                    grainboxflat_percent_values[p].resize(1+depth_max/depth_bin_width);
                    grainwidth_percent_values[p].resize(1+depth_max/depth_bin_width);
                    grainheight_percent_values[p].resize(1+depth_max/depth_bin_width);
                    grainangle_percent_values[p].resize(1+depth_max/depth_bin_width);
                    grainratio_percent_values[p].resize(1+depth_max/depth_bin_width);
                }
                else
                {
                    percent_depths_relax[p].resize(1+depth_max/depth_bin_width);
                    grainsize_percent_values_relax[p].resize(1+depth_max/depth_bin_width);
                    grainflattening_percent_values_relax[p].resize(1+depth_max/depth_bin_width);
                    grainround_percent_values_relax[p].resize(1+depth_max/depth_bin_width);
                    grainboxflat_percent_values_relax[p].resize(1+depth_max/depth_bin_width);
                    grainwidth_percent_values_relax[p].resize(1+depth_max/depth_bin_width);
                    grainheight_percent_values_relax[p].resize(1+depth_max/depth_bin_width);
                    grainangle_percent_values_relax[p].resize(1+depth_max/depth_bin_width);
                    grainratio_percent_values_relax[p].resize(1+depth_max/depth_bin_width);
                }

            }

            for(int q=0; q<10; q++)
            {
                if (iter==0)
                {
                    grainsize_quantile_depths[q].resize(1+depth_max/depth_bin_width);
                    grainsize_quantile_values[q].resize(1+depth_max/depth_bin_width);
                    grainabsangle_percent_values[q].resize(1+depth_max/depth_bin_width);
                    if (q==0)
                    {
                        grainsize_median_depths.resize(1+depth_max/depth_bin_width);
                        grainsize_median_values.resize(1+depth_max/depth_bin_width);
                        grainsize_5_largest_depths.resize(1+depth_max/depth_bin_width);
                        grainsize_5_largest_values.resize(1+depth_max/depth_bin_width);
                        grainsize_5_largest_errors.resize(1+depth_max/depth_bin_width);
                    }
                }
                else
                {
                    grainsize_quantile_depths_relax[q].resize(1+depth_max/depth_bin_width);
                    grainsize_quantile_values_relax[q].resize(1+depth_max/depth_bin_width);
                    grainabsangle_percent_values_relax[q].resize(1+depth_max/depth_bin_width);
                    if (q==0) grainsize_median_depths_relax.resize(1+depth_max/depth_bin_width);
                    if (q==0) grainsize_median_values_relax.resize(1+depth_max/depth_bin_width);
                }
            }

            for (int bin=0; bin<grain_sizes.size(); bin++)
            {
                if(grain_sizes[bin].size()==0) continue;

                if (iter==iterations-1)
                {
                    //grain size histogram
                    float standard_deviation;
                    calculate_mean_standard_deviation(grain_sizes[bin], grain_size_mean[bin], standard_deviation);
                    grainsize_combined_depths[0].push_back((bin+0.5f)*depth_bin_width);
                    grainsize_combined_values[0].push_back(grain_size_mean[bin]);

                    grainsize_percent_profile_depths.push_back((bin+0.5f)*depth_bin_width);
                    grainsize_quantile_profile_depths.push_back((bin+0.5f)*depth_bin_width);

                    std::vector<float> new_entry;
                    grainsize_percent_profile_values.push_back(new_entry);
                    grainsize_quantile_profile_values.push_back(new_entry);
                }

                std::sort(grain_sizes[bin].begin(),grain_sizes[bin].end(),size_cmp);

                std::vector<float> grain_size_region_values;
                std::vector<float> grain_flattening_region_values;
                std::vector<float> grain_round_region_values;
                std::vector<float> grain_boxflat_region_values;
                std::vector<float> grain_width_region_values;
                std::vector<float> grain_height_region_values;
                std::vector<float> grain_angle_region_values;
                std::vector<float> grain_absangle_region_values;
                std::vector<float> grain_ratio_region_values;
                int p=0;
                int q=0;

                int quantile_step_rounded = grain_sizes[bin].size()/100;
                int percent_step_rounded = grain_sizes[bin].size()/10;

                //grain size quantiles and percents
                for (int area_number=0; area_number<grain_sizes[bin].size(); ++area_number)
                {
                    grain_size_region_values.push_back(grain_sizes[bin][area_number].size);
                    grain_flattening_region_values.push_back(grain_sizes[bin][area_number].flat);
                    grain_round_region_values.push_back(grain_sizes[bin][area_number].round);
                    grain_boxflat_region_values.push_back(grain_sizes[bin][area_number].boxflat);
                    grain_width_region_values.push_back(grain_sizes[bin][area_number].width);
                    grain_height_region_values.push_back(grain_sizes[bin][area_number].height);
                    grain_angle_region_values.push_back(grain_sizes[bin][area_number].angle);
                    grain_absangle_region_values.push_back(fabs(grain_sizes[bin][area_number].angle));
                    grain_ratio_region_values.push_back(grain_sizes[bin][area_number].ratio);
                    if (quantile_step_rounded>0)
                        if((area_number+1)%quantile_step_rounded==0 && q<10)
                        {
                            if (iter==0)
                            {
                                grainsize_quantile_values[q][bin]=grain_sizes[bin][area_number].size;
                                grainsize_quantile_depths[q][bin]=(bin+0.5f)*depth_bin_width;
                            }
                            else
                            {
                                grainsize_quantile_values_relax[q][bin]=grain_sizes[bin][area_number].size;
                                grainsize_quantile_depths_relax[q][bin]=(bin+0.5f)*depth_bin_width;
                            }
                            q++;

                            if (iter==iterations-1)
                            {
                                grainsize_quantile_profile_values.back().push_back(grain_sizes[bin][area_number].size);

                                if (q==9)
                                {
                                    grainsize_combined_depths[3].push_back((bin+0.5f)*depth_bin_width);
                                    grainsize_combined_values[3].push_back(grain_sizes[bin][area_number].size);
                                }
                            }
                        }
                    if (percent_step_rounded>0)
                        if ((area_number+1)%percent_step_rounded==0)
                        {
                            float grain_size_region_mean;
                            float grain_size_region_standard_deviation;
                            calculate_mean_standard_deviation(grain_size_region_values, grain_size_region_mean, grain_size_region_standard_deviation);

                            float grain_flattening_region_mean;
                            float grain_flattening_region_standard_deviation;
                            calculate_mean_standard_deviation(grain_flattening_region_values, grain_flattening_region_mean,
                                grain_flattening_region_standard_deviation);

                            float grain_round_region_mean;
                            float grain_round_region_standard_deviation;
                            calculate_mean_standard_deviation(grain_round_region_values, grain_round_region_mean,
                                grain_round_region_standard_deviation);

                            float grain_boxflat_region_mean;
                            float grain_boxflat_region_standard_deviation;
                            calculate_mean_standard_deviation(grain_boxflat_region_values, grain_boxflat_region_mean,
                                grain_boxflat_region_standard_deviation);

                            float grain_width_region_mean;
                            float grain_width_region_standard_deviation;
                            calculate_mean_standard_deviation(grain_width_region_values, grain_width_region_mean,
                                grain_width_region_standard_deviation);

                            float grain_height_region_mean;
                            float grain_height_region_standard_deviation;
                            calculate_mean_standard_deviation(grain_height_region_values, grain_height_region_mean,
                                grain_height_region_standard_deviation);

                            float grain_angle_region_mean;
                            float grain_angle_region_standard_deviation;
                            calculate_mean_standard_deviation(grain_angle_region_values, grain_angle_region_mean,
                                grain_angle_region_standard_deviation);

                            float grain_absangle_region_mean;
                            float grain_absangle_region_standard_deviation;
                            calculate_mean_standard_deviation(grain_absangle_region_values, grain_absangle_region_mean,
                                grain_absangle_region_standard_deviation);

                            float grain_ratio_region_mean;
                            float grain_ratio_region_standard_deviation;
                            calculate_mean_standard_deviation(grain_ratio_region_values, grain_ratio_region_mean,
                                grain_ratio_region_standard_deviation);

                            if(p<9)
                            {
                                if (iter==0)
                                {
                                    percent_depths[p][bin]=(bin+0.5f)*depth_bin_width;
                                    grainsize_percent_values[p][bin]=grain_size_region_mean;
                                    grainflattening_percent_values[p][bin]=grain_flattening_region_mean;
                                    grainround_percent_values[p][bin]=grain_round_region_mean;
                                    grainboxflat_percent_values[p][bin]=grain_boxflat_region_mean;
                                    grainwidth_percent_values[p][bin]=grain_width_region_mean;
                                    grainheight_percent_values[p][bin]=grain_height_region_mean;
                                    grainangle_percent_values[p][bin]=grain_angle_region_mean;
                                    grainabsangle_percent_values[p][bin]=grain_absangle_region_mean;
                                    grainratio_percent_values[p][bin]=grain_ratio_region_mean;
                                }
                                else
                                {
                                    percent_depths_relax[p][bin]=(bin+0.5f)*depth_bin_width;
                                    grainsize_percent_values_relax[p][bin]=grain_size_region_mean;
                                    grainflattening_percent_values_relax[p][bin]=grain_flattening_region_mean;
                                    grainround_percent_values_relax[p][bin]=grain_round_region_mean;
                                    grainboxflat_percent_values_relax[p][bin]=grain_boxflat_region_mean;
                                    grainwidth_percent_values_relax[p][bin]=grain_width_region_mean;
                                    grainheight_percent_values_relax[p][bin]=grain_height_region_mean;
                                    grainangle_percent_values_relax[p][bin]=grain_angle_region_mean;
                                    grainabsangle_percent_values_relax[p][bin]=grain_absangle_region_mean;
                                    grainratio_percent_values_relax[p][bin]=grain_ratio_region_mean;
                                }
                                p++;
                            }

                            if (iter==iterations-1)
                            {
                                grainsize_percent_profile_values.back().push_back(grain_size_region_mean);

                                if (p==2)
                                {
                                    grainsize_combined_depths[2].push_back((bin+0.5f)*depth_bin_width);
                                    grainsize_combined_values[2].push_back(grain_size_region_mean);
                                }
                            }
                        }

                    if (area_number==grain_sizes[bin].size()-1)
                    {
                        float grain_absangle_region_mean;
                        float grain_absangle_region_standard_deviation;
                        calculate_mean_standard_deviation(grain_absangle_region_values, grain_absangle_region_mean,
                            grain_absangle_region_standard_deviation);
                        grainabsangle_percent_values[9][bin]=grain_absangle_region_mean;
                    }

                    if (area_number+1==grain_sizes[bin].size()/2)//median
                    {
                        if (iter==0)
                        {
                            grainsize_median_values[bin]=grain_sizes[bin][area_number].size;
                            grainsize_median_depths[bin]=(bin+0.5f)*depth_bin_width;
                        }
                        else
                        {
                            grainsize_median_values_relax[bin]=grain_sizes[bin][area_number].size;
                            grainsize_median_depths_relax[bin]=(bin+0.5f)*depth_bin_width;
                        }
                    }

                    if (area_number+1==5)//5 largest
                    {
                        float mean;
                        float standard_deviation;
                        calculate_mean_standard_deviation(grain_size_region_values, mean, standard_deviation);

                        grainsize_5_largest_values[bin]=mean;
                        grainsize_5_largest_depths[bin]=(bin+0.5f)*depth_bin_width;
                        grainsize_5_largest_errors[bin]=standard_deviation;
                    }
                }

                //grain size fit
                std::string filepath_grain_size_fit="fit.svg";

                float fit_mu, fit_sigma, stdabw_low, stdabw_high;

                plot.draw_histogram_left_lognorm(area_unit, "Relative occurrence","Grain size distribution", grain_size_histogram[bin], grain_bin_width,
                    area_scaling(), grain_size_start_mu, grain_size_start_sigma, grain_size_max_x[bin], stdabw_low, stdabw_high, grain_size_y_max, fit_mu,
                    fit_sigma, filepath_grain_size_fit.c_str());

                remove(filepath_grain_size_fit.c_str());

                if (iter==0)
                {
                    grainsize_fit_depths.push_back((bin+0.5f)*depth_bin_width);
                    grainsize_fit_values.push_back(grain_size_max_x[bin]);
                    grainsize_fit_errors_low.push_back(stdabw_low);
                    grainsize_fit_errors_high.push_back(std::min(1000.0f,stdabw_high));
                }
                else
                {
                    grainsize_fit_depths_relax.push_back((bin+0.5f)*depth_bin_width);
                    grainsize_fit_values_relax.push_back(grain_size_max_x[bin]);
                    grainsize_fit_errors_low_relax.push_back(stdabw_low);
                    grainsize_fit_errors_high_relax.push_back(std::min(1000.0f,stdabw_high));

                }
            }

            //consider longest grain boundaries
            if (iter==iterations-1)
            {
                for(int p=0; p<10; p++)
                {
                    disl_dens_percent_boundaries_depths[p].resize(1+depth_max/depth_bin_width);
                    disl_dens_percent_boundaries_values[p].resize(1+depth_max/depth_bin_width);
                }

                for (int bin=0; bin<boundary_sizes.size(); bin++)
                {
                    if(boundary_sizes[bin].size()==0) continue;

                    disl_dens_percent_boundaries_profile_depths.push_back((bin+0.5f)*depth_bin_width);

                    std::vector<float> new_entry;
                    disl_dens_percent_boundaries_profile_values.push_back(new_entry);

                    std::sort(boundary_sizes[bin].begin(),boundary_sizes[bin].end(),size_cmp);

                    std::vector<float> boundary_length_region_values;
                    std::vector<float> dislocation_density_values;
                    int p=0;

                    int percent_step_rounded = boundary_sizes[bin].size()/10;

                    //dislocation density percents
                    for (int boundary_number=0; boundary_number<boundary_sizes[bin].size();++boundary_number)
                    {
                        boundary_length_region_values.push_back(boundary_sizes[bin][boundary_number].size);
                        for(int rep=0; rep<boundary_sizes[bin][boundary_number].size; rep++)
                            dislocation_density_values.push_back(boundary_sizes[bin][boundary_number].flat);

                        if (percent_step_rounded>0)
                            if ((boundary_number+1)%percent_step_rounded==0)
                            {
                                float mean;
                                float standard_deviation;
                                calculate_mean_standard_deviation(dislocation_density_values, mean, standard_deviation);

                                disl_dens_percent_boundaries_depths[p][bin]=(bin+0.5f)*depth_bin_width;
                                disl_dens_percent_boundaries_values[p][bin]=mean/dislocation_density_scaling;
                                disl_dens_percent_boundaries_profile_values.back().push_back(mean/dislocation_density_scaling);
                                p++;
                            }
                    }
                }
            }

            list_file.close();
            if(temp_list_file.is_open()) temp_list_file.close();
        }

        std::stringstream s;
        if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size <<".svg";
        else s << "." << minimal_grain_size <<".svg";

        std::string filepath_grainsize_fit_out=path_results;
        filepath_grainsize_fit_out.append("new_grainsize_fit");
        filepath_grainsize_fit_out.append(s.str());

        std::string filepath_grainsize_fit_combined_out=path_results;
        filepath_grainsize_fit_combined_out.append("new_grainsize_fit_combined");
        filepath_grainsize_fit_combined_out.append(s.str());

        std::string filepath_grainsize_percent_out=path_results;
        filepath_grainsize_percent_out.append("percent/new_grainsize_");
        std::string filepath_grainflattening_percent_out=path_results;
        filepath_grainflattening_percent_out.append("percent/new_grainflattening_");
        std::string filepath_grainround_percent_out=path_results;
        filepath_grainround_percent_out.append("percent/new_grainshape_");
        std::string filepath_grainboxflat_percent_out=path_results;
        filepath_grainboxflat_percent_out.append("percent/new_grain_boxshape_");
        std::string filepath_grainwidth_percent_out=path_results;
        filepath_grainwidth_percent_out.append("percent/new_grainwidth_");
        std::string filepath_grainheight_percent_out=path_results;
        filepath_grainheight_percent_out.append("percent/new_grainheight_");
        std::string filepath_grainangle_percent_out=path_results;
        filepath_grainangle_percent_out.append("percent/new_grain_ellipseangle_");
        std::string filepath_grainabsangle_percent_out=path_results;
        filepath_grainabsangle_percent_out.append("percent/new_grain_abs_ellipseangle_");
        std::string filepath_grainratio_percent_out=path_results;
        filepath_grainratio_percent_out.append("percent/new_grain_perimeter_ratio_");

        std::string filepath_grainsize_quantile_out=path_results;
        filepath_grainsize_quantile_out.append("quantile/new_grainsize_");

        std::string filepath_grainsize_percent_profile_out=path_results;
        filepath_grainsize_percent_profile_out.append("new_grainsize_percent_profile");
        filepath_grainsize_percent_profile_out.append(s.str());

        std::string filepath_grainsize_quantile_profile_out=path_results;
        filepath_grainsize_quantile_profile_out.append("new_grainsize_quantile_profile");
        filepath_grainsize_quantile_profile_out.append(s.str());

        std::string filepath_grainsize_combined_out=path_results;
        filepath_grainsize_combined_out.append("new_grainsize_combined");
        filepath_grainsize_combined_out.append(s.str());

        std::string filepath_grainsize_median_out=path_results;
        filepath_grainsize_median_out.append("new_grainsize_median");
        filepath_grainsize_median_out.append(s.str());

        std::string filepath_grainsize_median_combined_out=path_results;
        filepath_grainsize_median_combined_out.append("new_grainsize_median_combined");
        filepath_grainsize_median_combined_out.append(s.str());

        std::string filepath_grainsize_5_largest_out=path_results;
        filepath_grainsize_5_largest_out.append("new_grainsize_5_largest");
        filepath_grainsize_5_largest_out.append(s.str());

        std::string filepath_disl_dens_percent_boundaries_out=path_results;
        filepath_disl_dens_percent_boundaries_out.append("percent/new_disl_dens_");

        std::string filepath_disl_dens_percent_boundaries_profile_out=path_results;
        filepath_disl_dens_percent_boundaries_profile_out.append("new_disl_dens_percent_profile");
        filepath_disl_dens_percent_boundaries_profile_out.append(s.str());

        if(minimal_bubble_distance==0 && grainsize_fit_values.size()>0) plot.draw_depth_errors("Depth in m", area_unit, 
                                                                 "Grain size profile log-normal fit", grainsize_fit_depths, grainsize_fit_values,
                                                                 grainsize_fit_errors_low, grainsize_fit_errors_high, 0.0f,
                                                                 filepath_grainsize_fit_out.c_str());
        else if(grainsize_fit_values_relax.size()>0) plot.draw_depth_errors("Depth in m", area_unit, "Grain size profile log-normal fit",
                                                                 grainsize_fit_depths_relax, grainsize_fit_values_relax, grainsize_fit_errors_low_relax,
                                                                 grainsize_fit_errors_high_relax, 0.0f, filepath_grainsize_fit_out.c_str());
        if(grainsize_fit_values.size()>0 && grainsize_fit_values_relax.size()>0) plot.draw_depth_errors("Depth in m", area_unit,
                                                                 "Grain size profile log-normal fit", grainsize_fit_depths, grainsize_fit_values,
                                                                 grainsize_fit_errors_low, grainsize_fit_errors_high, grainsize_fit_depths_relax,
                                                                 grainsize_fit_values_relax, grainsize_fit_errors_low_relax,
                                                                 grainsize_fit_errors_high_relax, 0.0f, filepath_grainsize_fit_combined_out.c_str());

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

            if(minimal_bubble_distance==0 && grainsize_percent_values[p].size()>0) plot.draw_depth("Depth in m", area_unit, titel,
                                                                 percent_depths[p], grainsize_percent_values[p], 0.0f, filepath_out.c_str());
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

            if(grainsize_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m", area_unit, titel, percent_depths_relax[p],
                                                                 grainsize_percent_values_relax[p], 0.0f, filepath_out.c_str());

            if(grainsize_percent_values.size()>0)
                if(grainsize_percent_values[p].size()>0 && grainsize_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m", area_unit, 
                                                                 titel, percent_depths[p], grainsize_percent_values[p],
                                                                 percent_depths_relax[p], grainsize_percent_values_relax[p], 0.0f,
                                                                 filepath_grainsize_percent_combined_out.c_str());
        }

        if(minimal_bubble_distance==0)
        for (int p=0; p<grainflattening_percent_values.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Vertical grain flattening profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainflattening_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(minimal_bubble_distance==0 && grainflattening_percent_values[p].size()>0) plot.draw_depth("Depth in m", "Vertical grain flattening factor",
                                                                 titel, percent_depths[p], grainflattening_percent_values[p], 0.9f, filepath_out.c_str());
        }
        else
        for (int p=0; p<grainflattening_percent_values_relax.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Vertical grain flattening profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainflattening_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            std::string filepath_grainflattening_percent_combined_out=filepath_grainflattening_percent_out;
            filepath_grainflattening_percent_combined_out.append(percent);
            filepath_grainflattening_percent_combined_out.append("percent_combined");
            filepath_grainflattening_percent_combined_out.append(s.str());
            filepath_grainflattening_percent_combined_out.resize(filepath_grainflattening_percent_combined_out.size()-3);
            filepath_grainflattening_percent_combined_out.append("svg");

            if(grainflattening_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m", "Vertical grain flattening factor", titel, percent_depths_relax[p],
                                                                 grainflattening_percent_values_relax[p], 0.9f, filepath_out.c_str());

            if(grainflattening_percent_values.size()>0)
                if(grainflattening_percent_values[p].size()>0 && grainflattening_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m",
                                                                 "Vertical grain flattening factor", titel, percent_depths[p], grainflattening_percent_values[p],
                                                                 percent_depths_relax[p], grainflattening_percent_values_relax[p], 0.9f,
                                                                 filepath_grainflattening_percent_combined_out.c_str());
        }

        if(minimal_bubble_distance==0)
        for (int p=0; p<grainround_percent_values.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Grain roundness profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainround_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(minimal_bubble_distance==0 && grainround_percent_values[p].size()>0) plot.draw_depth("Depth in m", "Roundness factor (4pi*Area/Perimeter^2)",
                                                                 titel, percent_depths[p], grainround_percent_values[p], 0.35f, filepath_out.c_str());
        }
        else
        for (int p=0; p<grainround_percent_values_relax.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Grain roundness profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainround_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            std::string filepath_grainround_percent_combined_out=filepath_grainround_percent_out;
            filepath_grainround_percent_combined_out.append(percent);
            filepath_grainround_percent_combined_out.append("percent_combined");
            filepath_grainround_percent_combined_out.append(s.str());
            filepath_grainround_percent_combined_out.resize(filepath_grainround_percent_combined_out.size()-3);
            filepath_grainround_percent_combined_out.append("svg");

            if(grainround_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m", "Roundness factor (4pi*Area/Perimeter^2)", titel, percent_depths_relax[p],
                                                                 grainround_percent_values_relax[p], 0.35f, filepath_out.c_str());

            if(grainround_percent_values.size()>0)
                if(grainround_percent_values[p].size()>0 && grainround_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m",
                                                                 "Roundness factor (4pi*Area/Perimeter^2)", titel, percent_depths[p], grainround_percent_values[p],
                                                                 percent_depths_relax[p], grainround_percent_values_relax[p], 0.35f,
                                                                 filepath_grainround_percent_combined_out.c_str());
        }

        if(minimal_bubble_distance==0)
        for (int p=0; p<grainboxflat_percent_values.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Box flattening profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainboxflat_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(minimal_bubble_distance==0 && grainboxflat_percent_values[p].size()>0) plot.draw_depth("Depth in m", "Box flattening factor",
                                                                 titel, percent_depths[p], grainboxflat_percent_values[p], 0.9f, filepath_out.c_str());
        }
        else
        for (int p=0; p<grainboxflat_percent_values_relax.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Box flattening profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainboxflat_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            std::string filepath_boxflat_percent_combined_out=filepath_grainboxflat_percent_out;
            filepath_boxflat_percent_combined_out.append(percent);
            filepath_boxflat_percent_combined_out.append("percent_combined");
            filepath_boxflat_percent_combined_out.append(s.str());
            filepath_boxflat_percent_combined_out.resize(filepath_boxflat_percent_combined_out.size()-3);
            filepath_boxflat_percent_combined_out.append("svg");

            if(grainboxflat_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m", "Box flattening factor", titel, percent_depths_relax[p],
                                                                 grainboxflat_percent_values_relax[p], 0.9f, filepath_out.c_str());

            if(grainboxflat_percent_values.size()>0)
                if(grainboxflat_percent_values[p].size()>0 && grainboxflat_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m",
                                                                 "Box flattening factor", titel, percent_depths[p], grainboxflat_percent_values[p],
                                                                 percent_depths_relax[p], grainboxflat_percent_values_relax[p], 0.9f,
                                                                 filepath_boxflat_percent_combined_out.c_str());
        }

        if(minimal_bubble_distance==0)
        for (int p=0; p<grainwidth_percent_values.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Grain width profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainwidth_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(minimal_bubble_distance==0 && grainwidth_percent_values[p].size()>0) plot.draw_depth("Depth in m", length_unit,
                                                                 titel, percent_depths[p], grainwidth_percent_values[p], 0.0f, filepath_out.c_str());
        }
        else
        for (int p=0; p<grainwidth_percent_values_relax.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Grain width profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainwidth_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            std::string filepath_grainwidth_percent_combined_out=filepath_grainwidth_percent_out;
            filepath_grainwidth_percent_combined_out.append(percent);
            filepath_grainwidth_percent_combined_out.append("percent_combined");
            filepath_grainwidth_percent_combined_out.append(s.str());
            filepath_grainwidth_percent_combined_out.resize(filepath_grainwidth_percent_combined_out.size()-3);
            filepath_grainwidth_percent_combined_out.append("svg");

            if(grainwidth_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m", length_unit, titel, percent_depths_relax[p],
                                                                 grainwidth_percent_values_relax[p], 0.0f, filepath_out.c_str());

            if(grainwidth_percent_values.size()>0)
                if(grainwidth_percent_values[p].size()>0 && grainwidth_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m",
                                                                 length_unit, titel, percent_depths[p], grainwidth_percent_values[p],
                                                                 percent_depths_relax[p], grainwidth_percent_values_relax[p], 0.0f,
                                                                 filepath_grainwidth_percent_combined_out.c_str());
        }

        if(minimal_bubble_distance==0)
        for (int p=0; p<grainheight_percent_values.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Grain height profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainheight_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(minimal_bubble_distance==0 && grainheight_percent_values[p].size()>0) plot.draw_depth("Depth in m", length_unit,
                                                                 titel, percent_depths[p], grainheight_percent_values[p], 0.0f, filepath_out.c_str());
        }
        else
        for (int p=0; p<grainheight_percent_values_relax.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Grain height profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainheight_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            std::string filepath_grainheight_percent_combined_out=filepath_grainheight_percent_out;
            filepath_grainheight_percent_combined_out.append(percent);
            filepath_grainheight_percent_combined_out.append("percent_combined");
            filepath_grainheight_percent_combined_out.append(s.str());
            filepath_grainheight_percent_combined_out.resize(filepath_grainheight_percent_combined_out.size()-3);
            filepath_grainheight_percent_combined_out.append("svg");

            if(grainheight_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m", length_unit, titel, percent_depths_relax[p],
                                                                 grainheight_percent_values_relax[p], 0.0f, filepath_out.c_str());

            if(grainheight_percent_values.size()>0)
                if(grainheight_percent_values[p].size()>0 && grainheight_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m",
                                                                 length_unit, titel, percent_depths[p], grainheight_percent_values[p],
                                                                 percent_depths_relax[p], grainheight_percent_values_relax[p], 0.0f,
                                                                 filepath_grainheight_percent_combined_out.c_str());
        }

        if(minimal_bubble_distance==0)
        for (int p=0; p<grainangle_percent_values.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Grain angle profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainangle_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(minimal_bubble_distance==0 && grainangle_percent_values[p].size()>0) plot.draw_depth("Depth in m", "Ellipse long axis angle",
                                                                 titel, percent_depths[p], grainangle_percent_values[p], 0.0f, filepath_out.c_str());
        }
        else
        for (int p=0; p<grainangle_percent_values_relax.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Grain angle profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainangle_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            std::string filepath_grainangle_percent_combined_out=filepath_grainangle_percent_out;
            filepath_grainangle_percent_combined_out.append(percent);
            filepath_grainangle_percent_combined_out.append("percent_combined");
            filepath_grainangle_percent_combined_out.append(s.str());
            filepath_grainangle_percent_combined_out.resize(filepath_grainangle_percent_combined_out.size()-3);
            filepath_grainangle_percent_combined_out.append("svg");

            if(grainangle_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m", "Ellipse long axis angle", titel, percent_depths_relax[p],
                                                                 grainangle_percent_values_relax[p], 0.0f, filepath_out.c_str());

            if(grainangle_percent_values.size()>0)
                if(grainangle_percent_values[p].size()>0 && grainangle_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m",
                                                                 "Ellipse long axis angle", titel, percent_depths[p], grainangle_percent_values[p],
                                                                 percent_depths_relax[p], grainangle_percent_values_relax[p], 0.0f,
                                                                 filepath_grainangle_percent_combined_out.c_str());
        }

        if(minimal_bubble_distance==0)
        for (int p=0; p<grainabsangle_percent_values.size(); p++)
        {
            char percent[4];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Absolute grain angle profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainabsangle_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(minimal_bubble_distance==0 && grainabsangle_percent_values[p].size()>0) plot.draw_depth("Depth in m", "Absolut ellipse long axis angle",
                                                                 titel, grainsize_quantile_depths[p], grainabsangle_percent_values[p], 0.0f, filepath_out.c_str());
        }
        else
        for (int p=0; p<grainabsangle_percent_values_relax.size(); p++)
        {
            char percent[4];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Absolute grain angle profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainabsangle_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            std::string filepath_grainabsangle_percent_combined_out=filepath_grainabsangle_percent_out;
            filepath_grainabsangle_percent_combined_out.append(percent);
            filepath_grainabsangle_percent_combined_out.append("percent_combined");
            filepath_grainabsangle_percent_combined_out.append(s.str());
            filepath_grainabsangle_percent_combined_out.resize(filepath_grainabsangle_percent_combined_out.size()-3);
            filepath_grainabsangle_percent_combined_out.append("svg");

            if(grainabsangle_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m", "Absolut ellipse long axis angle", titel,
                                                                 grainsize_quantile_depths_relax[p], grainabsangle_percent_values_relax[p], 0.0f,
                                                                 filepath_out.c_str());

            if(grainabsangle_percent_values.size()>0)
                if(grainabsangle_percent_values[p].size()>0 && grainabsangle_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m",
                                                                 "Absolut ellipse long axis angle", titel, grainsize_quantile_depths[p],
                                                                 grainabsangle_percent_values[p], grainsize_quantile_depths_relax[p],
                                                                 grainabsangle_percent_values_relax[p], 0.0f, filepath_grainabsangle_percent_combined_out.c_str());
        }

        if(minimal_bubble_distance==0)
        for (int p=0; p<grainratio_percent_values.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Perimeter ratio profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainratio_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            if(minimal_bubble_distance==0 && grainratio_percent_values[p].size()>0) plot.draw_depth("Depth in m", "Perimeter ratio",
                                                                 titel, percent_depths[p], grainratio_percent_values[p], 0.7f, filepath_out.c_str());
        }
        else
        for (int p=0; p<grainratio_percent_values_relax.size(); p++)
        {
            char percent[3];
            sprintf(percent, "%i", (p+1)*10);

            std::string titel="Perimeter ratio profile (considering ";
            titel.append(percent);
            titel.append(" percent largest grains)");

            std::string filepath_out=filepath_grainratio_percent_out;
            filepath_out.append(percent);
            filepath_out.append("percent");
            filepath_out.append(s.str());
            filepath_out.resize(filepath_out.size()-3);
            filepath_out.append("svg");

            std::string filepath_grainratio_percent_combined_out=filepath_grainratio_percent_out;
            filepath_grainratio_percent_combined_out.append(percent);
            filepath_grainratio_percent_combined_out.append("percent_combined");
            filepath_grainratio_percent_combined_out.append(s.str());
            filepath_grainratio_percent_combined_out.resize(filepath_grainratio_percent_combined_out.size()-3);
            filepath_grainratio_percent_combined_out.append("svg");

            if(grainratio_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m", "Perimeter ratio", titel, percent_depths_relax[p],
                                                                 grainratio_percent_values_relax[p], 0.7f, filepath_out.c_str());

            if(grainratio_percent_values.size()>0)
                if(grainratio_percent_values[p].size()>0 && grainratio_percent_values_relax[p].size()>0) plot.draw_depth("Depth in m",
                                                                 "Perimeter ratio", titel, percent_depths[p], grainratio_percent_values[p],
                                                                 percent_depths_relax[p], grainratio_percent_values_relax[p], 0.7f,
                                                                 filepath_grainratio_percent_combined_out.c_str());
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

            if(minimal_bubble_distance==0 && grainsize_quantile_values[q].size()>0) plot.draw_depth("Depth in m", area_unit, titel,
                                                                 grainsize_quantile_depths[q], grainsize_quantile_values[q], 0.0f, filepath_out.c_str());
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

            if(grainsize_quantile_values_relax[q].size()>0) plot.draw_depth("Depth in m", area_unit, titel, grainsize_quantile_depths_relax[q],
                                                                 grainsize_quantile_values_relax[q], 0.0f, filepath_out.c_str());

            if(grainsize_quantile_values.size()>0)
                if(grainsize_quantile_values[q].size()>0 && grainsize_quantile_values_relax[q].size()>0) plot.draw_depth("Depth in m", area_unit, 
                                                                 titel, grainsize_quantile_depths[q], grainsize_quantile_values[q],
                                                                 grainsize_quantile_depths_relax[q], grainsize_quantile_values_relax[q], 0.0f,
                                                                 filepath_grainsize_quantile_combined_out.c_str());
        }

        if(grainsize_percent_profile_values.size()>0) plot.draw_depth_3d("Largest grains", "Depth in m", area_unit, "Largest grains profile (relative)",
                                                                 grainsize_percent_profile_x_values, grainsize_percent_profile_depths,
                                                                 grainsize_percent_profile_values, filepath_grainsize_percent_profile_out.c_str());
        if(grainsize_quantile_profile_values.size()>0) plot.draw_depth_3d("Quantile 1-", "Depth in m", area_unit, "Grain size quantiles",
                                                                 grainsize_quantile_profile_x_values, grainsize_quantile_profile_depths,
                                                                 grainsize_quantile_profile_values, filepath_grainsize_quantile_profile_out.c_str());

        if(grainsize_combined_values.size()>0) plot.draw_depth("Depth in m", area_unit, "Different grain size parameters", grainsize_combined_depths,
                                                                 grainsize_combined_values, filepath_grainsize_combined_out.c_str());

        if(minimal_bubble_distance==0 && grainsize_median_values.size()>0) plot.draw_depth("Depth in m", area_unit, "Grain size medians", grainsize_median_depths,
                                                                 grainsize_median_values, 0.0f, filepath_grainsize_median_out.c_str());
        else
        {
            if(grainsize_median_values_relax.size()>0) plot.draw_depth("Depth in m", area_unit, "Grain size medians", grainsize_median_depths_relax,
                                                                 grainsize_median_values_relax, 0.0f, filepath_grainsize_median_out.c_str());

            if(grainsize_median_values.size()>0)
                if(grainsize_median_values.size()>0 && grainsize_median_values_relax.size()>0) plot.draw_depth("Depth in m", area_unit, 
                                                                 "Grain size medians", grainsize_median_depths, grainsize_median_values,
                                                                 grainsize_median_depths_relax, grainsize_median_values_relax, 0.0f,
                                                                 filepath_grainsize_median_combined_out.c_str());
        }

        if(grainsize_5_largest_values.size()>0) plot.draw_depth_errors("Depth in m", area_unit, "Mean grain size profile (5 largest grains)",
                                                                 grainsize_5_largest_depths, grainsize_5_largest_values, grainsize_5_largest_errors,
                                                                 grainsize_5_largest_errors,  0.0f, filepath_grainsize_5_largest_out.c_str());

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

            if(disl_dens_percent_boundaries_values[p].size()>0) plot.draw_depth("Depth in m", "Dislocation density difference in m^(-2)", titel,
                                                                disl_dens_percent_boundaries_depths[p], disl_dens_percent_boundaries_values[p], 0.0f,
                                                                filepath_out.c_str());
        }

        if(disl_dens_percent_boundaries_profile_values.size()>0) plot.draw_depth_3d("Longest boundaries", "Depth in m", "Dislocation density difference in m^(-2)",
                                                                 "Disl. dens. diff. at longest boundaries profile (relative)",
                                                                 disl_dens_percent_boundaries_profile_x_values, disl_dens_percent_boundaries_profile_depths,
                                                                 disl_dens_percent_boundaries_profile_values,
                                                                 filepath_disl_dens_percent_boundaries_profile_out.c_str());
    }
}
