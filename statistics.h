/*! \file statistics.h
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
#pragma once
#include <list>
#include <vigra/seededregiongrowing.hxx>

#include "math_statistics.h"
#include "structures_statistics.h"
#include "plplot.h"
#include "FitEllipse.h"
#include "gbn.h"
#include "param.h"

#include "ConvexPerimeter.h"

void do_statistics(std::string filepath_to_feature_file, std::string path_to_ws_image, std::string path_rf_predictions, std::string path_plots,
                   std::string path_results, std::string param_file, ParameterFile paramFile)
{
    std::string param_file_name=get_filename(param_file);
    if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
    else param_file_name = "";

    param para;

    Parameter<int> low_grain_size;
    low_grain_size.assign("", "low_grain_size", 0);
    low_grain_size.load(paramFile,"config");

    Parameter<int> high_grain_size;
    high_grain_size.assign("", "high_grain_size", 20000);
    high_grain_size.load(paramFile,"config");

    Parameter<int> grain_size_step;
    grain_size_step.assign("", "grain_size_step", 5000);
    grain_size_step.load(paramFile,"config");

    Parameter<int> minimal_bubble_distance;
    minimal_bubble_distance.assign("", "min_bubble_distance", 0);
    minimal_bubble_distance.load(paramFile,"config");

    Parameter<int> close_bubble_grain_size;
    close_bubble_grain_size.assign("", "close_bubble_grain_size", 1000);
    close_bubble_grain_size.load(paramFile,"config");

    Parameter<int> grain_step;
    grain_step.assign("", "grain_step", 100);
    grain_step.load(paramFile,"config");

    Parameter<int> boundary_step;
    boundary_step.assign("", "boundary_step", 100);
    boundary_step.load(paramFile,"config");

    Parameter<int> grain_size_min;
    grain_size_min.assign("", "grain_size_min", 4);
    grain_size_min.load(paramFile,"config");

    Parameter<float> percentage_grains;
    percentage_grains.assign("", "percentage_grains", 0.95f);
    percentage_grains.load(paramFile,"config");

    Parameter<float> length_scaling;
    length_scaling.assign("", "length_scaling", 193.5f);
    length_scaling.load(paramFile,"config");

    Parameter<float> area_scaling;
    area_scaling.assign("", "area_scaling", 37444.0f);
    area_scaling.load(paramFile,"config");

    std::string filepath_new_classification=path_rf_predictions;
    filepath_new_classification.append(param_file_name.c_str());

    std::string filename=get_filename(filepath_to_feature_file);
    //remove the ".bin"
    filename.resize(filename.size()-4);

    filepath_new_classification.append(filename);

    std::string filepath_plots=path_plots;
    filepath_plots.append(filename);

    std::string filepath_statistics_results=filepath_new_classification;
    filepath_statistics_results.append(".txt");
    filepath_new_classification.append(".h5");

    //IMPORT RESULTS FROM HDF5 file
    gbn GrainBoundNet;
    GrainBoundNet.load_final_structure(filepath_new_classification);

    size_t & nr_areas = GrainBoundNet.nr_new_areas;
    long * & bubble_area_size = GrainBoundNet.bubble_area_size;
    long * & grain_area_size = GrainBoundNet.grain_area_size;
    std::vector< std::vector< std::vector<int> > > grain_arc_index;
    std::vector< std::vector<int> > & old_grain_arc_index = GrainBoundNet.grain_arc_index;
    std::vector< std::vector<int> > & bubble_arc_index = GrainBoundNet.bubble_arc_index;
    std::vector<size_t> & region_labels = GrainBoundNet.region_labels;
    std::vector<int> & grain_junctions = GrainBoundNet.grain_junctions;
    std::vector<int> & grain_bubble_junctions = GrainBoundNet.grain_bubble_junctions;
    std::vector<point> & grain_area_center_mass = GrainBoundNet.grain_area_center_mass;
    std::vector<point> & bubble_area_center_mass = GrainBoundNet.bubble_area_center_mass;
    std::vector<bool> & grain_arc = GrainBoundNet.grain_arc;
    std::vector<int> & found_bubble_areas = GrainBoundNet.found_bubble_areas;

    //create new grain arc index with segments
    grain_arc_index.resize(nr_areas);

    for (int area=0; area<nr_areas; area++)
    {
        std::vector<int> segment_index;

        if(old_grain_arc_index[area].size()>0)
        {
            grain_arc_index[area].push_back(segment_index);//all arcs are added to first segment
            grain_arc_index[area][0]=old_grain_arc_index[area];
        }
    }

    old_grain_arc_index.clear();

    std::string filepath_to_ws_image=path_to_ws_image;
    filepath_to_ws_image.append(filename);
    filepath_to_ws_image.append(".h5");

    //IMPORT RESULTS FROM HDF5 file
    seg segment(false);
    segment.load_cgp_data_structure(filepath_to_ws_image);

    vigra::BasicImage<unsigned int> & ws_region_image = segment.ws_region_image;
    marray::Marray<unsigned int> & one_boundings = segment.one_boundings;
    marray::Marray<unsigned int> & two_boundings = segment.two_boundings;
    std::vector< std::vector<point> > & arcs = segment.arcs;
    std::vector<point> & junctions = segment.junctions;
    int & dim_x = segment.dim_x;
    int & dim_y = segment.dim_y;

    grain_arc.resize(arcs.size(),false);

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

    //grain size if certain percentage of area is filled
    std::vector<float> percentage_filled_grains_mean;
    std::vector<float> percentage_filled_grains_standard_deviation;
    std::vector<float> percentage_values;
    std::vector<float> nr_percentage_filled_grains;

    //grain size histogram
    int grain_size_max=0;
    float grain_bin_width=5.0f;
    std::vector<int> grain_size_histogram;
    std::vector<float> grain_size_values;
    float grain_size_mean;
    float grain_size_standard_deviation;
    float grain_size_max_x;
    float grain_size_stdabw_low;
    float grain_size_stdabw_high;
    float grain_size_start_mu=2.0f;
    float grain_size_start_sigma=0.4f;
    float grain_size_y_max=0.0f;

    //grain equiv radius histogram
    int grain_equiv_radius_min=0;
    int grain_equiv_radius_max=0;
    float grain_equiv_radius_bin_width=10.0f;
    std::vector<int> grain_equiv_radius_histogram;
    std::vector<float> grain_equiv_radius_values;
    float grain_equiv_radius_mean;
    float grain_equiv_radius_standard_deviation;
    float grain_equiv_radius_y_max=0.0f;

    //grain equiv radius norm histogram
    float grain_equiv_radius_norm_max=0.0f;
    float grain_equiv_radius_norm_bin_width=0.2f;
    std::vector<int> grain_equiv_radius_norm_histogram;
    float grain_equiv_radius_norm_y_max=0.0f;

    //grain roundness histogram
    float grain_roundness_max=1.0f;
    float grain_roundness_bin_width=0.025f;
    std::vector<int> grain_roundness_histogram;
    std::vector<float> & grain_roundness_values = para.grain_roundness;
    float grain_roundness_mean;
    float grain_roundness_standard_deviation;
    float grain_roundness_y_max=0.0f;

    //grain box flattening histogram
    float grain_box_flattening_max=0.0f;
    float grain_box_flattening_bin_width=0.25f;
    std::vector<int> grain_box_flattening_histogram;
    std::vector<float> & grain_box_flattening_values = para.grain_box_flattening;
    std::vector<float> grain_box_flattening_values2;
    float grain_box_flattening_mean;
    float grain_box_flattening_standard_deviation;
    float grain_box_flattening_y_max=0.0f;

    //grain box width histogram
    int grain_box_width_max=0;
    float grain_box_width_bin_width=7.5f;
    std::vector<int> grain_box_width_histogram;
    std::vector<float> & grain_box_width_values = para.grain_box_width;
    std::vector<float> grain_box_width_values2;
    float grain_box_width_mean;
    float grain_box_width_standard_deviation;
    float grain_box_width_y_max=0.0f;

    //grain box height histogram
    int grain_box_height_max=0;
    float grain_box_height_bin_width=7.5f;
    std::vector<int> grain_box_height_histogram;
    std::vector<float> & grain_box_height_values = para.grain_box_height;
    std::vector<float> grain_box_height_values2;
    float grain_box_height_mean;
    float grain_box_height_standard_deviation;
    float grain_box_height_y_max=0.0f;

    //parameters ellipse fit
    std::vector< std::vector<float> > & grain_ellipse_params = para.ellipse_params;
    FitEllipse fitEllipse;//initialize fit class

    //Initialize convex perimeter class
    ConvexPerimeter convPerimeter; 

    //grain ellipse long axis histogram
    int grain_ellipse_long_axis_min=0;
    int grain_ellipse_long_axis_max=0;
    float grain_ellipse_long_axis_bin_width=7.5f;
    std::vector<int> grain_ellipse_long_axis_histogram;
    std::vector<float> & grain_ellipse_long_axis_values = para.ellipse_long_axis;
    std::vector<float> grain_ellipse_long_axis_values2;
    float grain_ellipse_long_axis_mean;
    float grain_ellipse_long_axis_standard_deviation;
    float grain_ellipse_long_axis_y_max=0.0f;

    //grain ellipse flattening histogram
    float grain_ellipse_flattening_min=0.0f;
    float grain_ellipse_flattening_max=0.0f;
    float grain_ellipse_flattening_bin_width=0.25f;
    std::vector<int> grain_ellipse_flattening_histogram;
    std::vector<float> & grain_ellipse_flattening_values = para.ellipse_flattening;
    std::vector<float> grain_ellipse_flattening_values2;
    float grain_ellipse_flattening_mean;
    float grain_ellipse_flattening_standard_deviation;
    float grain_ellipse_flattening_y_max=0.0f;

    //grain ellipse long axis angle histogram
    float angle_scaling=2.0f*PI/360;
    float grain_ellipse_long_axis_angle_max=0.0f;
    float grain_ellipse_long_axis_angle_bin_width=PI/36.0f;
    std::vector<int> grain_ellipse_long_axis_angle_histogram;
    std::vector<float> & grain_ellipse_long_axis_angle_values = para.ellipse_long_axis_angle;
    std::vector<float> grain_ellipse_long_axis_angle_values2;
    float grain_ellipse_long_axis_angle_mean;
    float grain_ellipse_long_axis_angle_standard_deviation;
    float grain_ellipse_long_axis_angle_y_max=0.0f;

    //grain perimeter ratio histogram
    float grain_perimeter_ratio_min=0.0f;
    float grain_perimeter_ratio_max=1.0f;
    float grain_perimeter_ratio_bin_width=0.025f;
    std::vector<int> grain_perimeter_ratio_histogram;
    std::vector<float> & grain_perimeter_ratio_values = para.grain_perimeter_ratio;
    float grain_perimeter_ratio_mean;
    float grain_perimeter_ratio_standard_deviation;
    float grain_perimeter_ratio_y_max=0.0f;

    //grain area width histogram
    int grain_area_width_min=0;
    int grain_area_width_max=0;
    float grain_area_width_bin_width=7.5f;
    std::vector<int> grain_area_width_histogram;
    std::vector<float> & grain_area_width_values = para.grain_area_width;
    float grain_area_width_mean;
    float grain_area_width_standard_deviation;
    float grain_area_width_y_max=0.0f;

    //grain area height histogram
    int grain_area_height_min=0;
    int grain_area_height_max=0;
    float grain_area_height_bin_width=7.5f;
    std::vector<int> grain_area_height_histogram;
    std::vector<float> & grain_area_height_values = para.grain_area_height;
    float grain_area_height_mean;
    float grain_area_height_standard_deviation;
    float grain_area_height_y_max=0.0f;

    //grain area flattening histogram
    float grain_area_flattening_min=0.0f;
    float grain_area_flattening_max=0.0f;
    float grain_area_flattening_bin_width=0.25f;
    std::vector<int> grain_area_flattening_histogram;
    std::vector<float> grain_area_flattening_values;
    float grain_area_flattening_mean;
    float grain_area_flattening_standard_deviation;
    float grain_area_flattening_y_max=0.0f;

    //min_bubble_distance histogram
    int min_bubble_distance_min=1;
    int min_bubble_distance_max=0;
    float min_bubble_distance_bin_width=10.0f;
    int nr_size_regions=5;
    std::vector< std::vector<int> > min_bubble_distance_histogram(nr_size_regions);
    std::vector< std::vector<float> > min_bubble_distance_values(nr_size_regions);
    std::vector<float> min_bubble_distance_mean(nr_size_regions);
    std::vector<float> min_bubble_distance_standard_deviation(nr_size_regions);
    float min_bubble_distance_y_max=0.0f;

    //grain arc number histogram
    int grain_arc_number_min=0;
    int grain_arc_number_max=0;
    float grain_arc_number_bin_width=1.0f;
    std::vector<int> grain_arc_number_histogram;
    std::vector<float> & grain_arc_number_values = para.grain_arc_number;
    std::vector<float> grain_arc_number_values2;
    float grain_arc_number_max_x;
    float grain_arc_number_standard_deviation;
    float grain_arc_number_start_mu=1.8f;
    float grain_arc_number_start_sigma=0.58f;
    float grain_arc_number_y_max=0.0f;

    //grain neighbors histogram
    int grain_neighbors_min=0;
    int grain_neighbors_max=0;
    float grain_neighbors_bin_width=1.0f;
    std::vector<int> grain_neighbors_histogram;
    std::vector<float> & grain_neighbors_values = para.grain_neighbors;
    std::vector<float> grain_neighbors_values2;
    float grain_neighbors_max_x;
    float grain_neighbors_standard_deviation;
    float grain_neighbors_start_mu=1.8f;
    float grain_neighbors_start_sigma=0.58f;
    float grain_neighbors_y_max=0.0f;

    //grain arc length histogram
    int grain_arc_length_min=0;
    int grain_arc_length_max=0;
    float grain_arc_length_bin_width=7.5f;
    std::vector<int> grain_arc_length_histogram;
    std::vector<float> grain_arc_length_values;
    float grain_arc_length_mean;
    float grain_arc_length_standard_deviation;
    float grain_arc_length_y_max=0.0f;

    //grain longest arc length histogram
    int grain_longest_arc_length_min=0;
    int grain_longest_arc_length_max=0;
    float grain_longest_arc_length_bin_width=9.0f;
    std::vector<int> grain_longest_arc_length_histogram;
    std::vector<float> & grain_longest_arc_length_values = para.grain_longest_arc_length;
    std::vector<float> grain_longest_arc_length_values2;
    float grain_longest_arc_length_mean;
    float grain_longest_arc_length_standard_deviation;
    float grain_longest_arc_length_y_max=0.0f;

    //dihedral angle histogram
    float dihedral_angle_max=0.0f;
    float dihedral_angle_bin_width=PI/24.0f;
    std::vector<int> dihedral_angle_histogram;
    std::vector<float> dihedral_angle_values;
    float dihedral_angle_mean;
    float dihedral_angle_standard_deviation;
    float dihedral_angle_y_max=0.0f;
    int pixels_average=15;

    //dihedral angle 2 histogram
    float dihedral_angle2_max=0.0f;
    float dihedral_angle2_bin_width=PI/24.0f;
    std::vector<int> dihedral_angle2_histogram;
    std::vector<float> dihedral_angle2_values;
    float dihedral_angle2_mean;
    float dihedral_angle2_standard_deviation;
    float dihedral_angle2_y_max=0.0f;

    //dislocation density histogram
    float dislocation_density_scaling=1.0f/82814814815850.0f;
    float curvature_min=0.005f;
    float curvature_max=0.0f;
    float dislocation_density_bin_width=20.0f;
    std::vector<float> dislocation_density_histogram;
    std::vector<float> dislocation_density_values;
    float dislocation_density_mean;
    float dislocation_density_standard_deviation;
    float dislocation_density_y_max=0.0f;

    //boundary orientation histogram
    float phi_max=0.0f;
    float boundary_orientation_bin_width=PI/36.0f;
    std::vector<int> boundary_orientation_histogram;
    std::vector<float> boundary_orientation_values;
    float boundary_orientation_mean;
    float boundary_orientation_standard_deviation;
    float boundary_orientation_y_max=0.0f;

    //bubble size histogram
    int bubble_size_min=1;
    int bubble_size_max=0;
    float bubble_bin_width=5.0f;
    std::vector<int> bubble_size_histogram;
    std::vector<float> bubble_size_values;
    float bubble_size_mean;
    float bubble_size_standard_deviation;
    float bubble_size_y_max=0.0f;

    //bubble grain neighbors histogram
    int bubble_grain_neighbors_min=0;
    int bubble_grain_neighbors_max=0;
    float bubble_grain_neighbors_bin_width=1.0f;
    std::vector<int> bubble_grain_neighbors_histogram;
    std::vector<float> bubble_grain_neighbors_values;
    float bubble_grain_neighbors_mean;
    float bubble_grain_neighbors_standard_deviation;
    float bubble_grain_neighbors_y_max=0.0f;

    //grain center of mass histogram
    int grain_center_of_mass_min=0;
    int grain_center_of_mass_max=0;
    float grain_center_of_mass_bin_width=700.0f;
    std::vector<int> grain_center_of_mass_histogram;
    std::vector<float> grain_center_of_mass_values;
    float grain_center_of_mass_mean;
    float grain_center_of_mass_standard_deviation;
    float grain_center_of_mass_y_max=0.0f;

    //grain size region/percent/quantile plots
    std::vector<float> grain_size_region_values;
    std::vector<float> grain_region_values;
    std::vector<float> grain_region_errors;
    std::vector<float> grain_percent_values;
    std::vector<float> grain_percent_errors;
    std::vector<float> grain_quantile_values;
    std::vector<float> grain_quantile_errors;

    //relaxation plots
    std::vector< std::vector<float> > relaxation_values;
    std::vector<float> relaxation_mean;
    std::vector<float> relaxation_standard_deviation;
    std::vector< std::vector<float> > relaxation2_values;
    std::vector<float> relaxation2_mean;
    std::vector<float> relaxation2_standard_deviation;

    //turning points histogram
    float turning_point_max=0.0f;
    std::vector<int> turning_point_histogram;
    std::vector<float> & turning_point_values = para.turning_point;
    float turning_point_mean;
    float turning_point_standard_deviation;

    //boundary size dislocation density region/percent plots
    std::vector<float> disl_dens_boundary_region_values;
    std::vector<float> disl_dens_boundary_region_errors;
    std::vector<float> disl_dens_boundary_percent_values;
    std::vector<float> disl_dens_boundary_percent_errors;

    //correlations
    float correlation_ellipse_box_flattening;

    //variables for merging
    std::vector< std::vector<point> > areas(nr_areas);
    std::vector<area_range> area_ranges(nr_areas);
    std::list<int> found_border_areas;
    std::vector<int> arc_state(arcs.size());
    std::vector< std::vector<int> > region_arc_index;
    region_arc_index.resize(nr_areas);
    marray::Marray<unsigned int> twoCellNeighbors;
    vigra::MultiArray<2, double> not_used;
    std::vector<int> close_bubble_areas;

    //global structures created for statistics
    std::vector< std::vector<int> > grain_area_junctions(nr_areas);
    std::vector<bool> & grain_junction = para.grain_junction;
    grain_junction.resize(one_boundings.shape(0),false);
    std::vector<std::vector<int> > arc_junctions(two_boundings.shape(0));
    std::vector<int> grain_perimeter(nr_areas,0);
    std::vector<int> min_bubble_distance(nr_areas);
    std::vector<int> grain_longest_arc_length(nr_areas);
    std::vector< std::vector<float> > grain_junction_angles(nr_areas);
    std::vector< std::vector<float> > grain_junction_angles2(nr_areas);
    std::vector< std::vector<point> > grain_boundary_pixels;
    std::vector< std::vector<int> > & grain_boundary_index = para.grain_boundary_index;
    std::vector< std::vector<float> > & grain_boundary_phis = para.grain_boundary_phis;
    std::vector< std::vector<float> > & grain_boundary_curvs = para.grain_boundary_curvs;
    std::vector< std::vector<int> > & grain_junction_index = para.grain_junction_index;
    grain_junction_index.resize(junctions.size());

    //region image contains current regions and will be updated during merging
    vigra::BasicImage<unsigned int> region_image(dim_x,dim_y);

    //this vector is used to know which arcs have been used for segment combination
    std::vector<bool> arc_to_segment(arcs.size(),false);

    initialize_structures(segment,
                          GrainBoundNet,
                          para,
                          grain_arc_index,
                          minimal_bubble_distance(),
                          close_bubble_grain_size,
                          pixels_average,
                          areas,
                          area_ranges,
                          found_border_areas,
                          arc_state,
                          region_arc_index,
                          twoCellNeighbors,
                          close_bubble_areas,
                          grain_area_junctions,
                          arc_junctions,
                          grain_perimeter,
                          min_bubble_distance,
                          grain_longest_arc_length,
                          grain_junction_angles,
                          grain_junction_angles2,
                          grain_boundary_pixels,
                          region_image,
                          arc_to_segment);

    //*******************************************************************
    //LOOP OVER MINIMAL GRAIN SIZE, COMBINED WITH MINIMAL BUBBLE DISTANCE
    //*******************************************************************

    for(int minimal_grain_size=low_grain_size; minimal_grain_size<high_grain_size+1; minimal_grain_size+=grain_size_step)
    {
        if(minimal_grain_size>0 || minimal_bubble_distance>0)
        {
            if (minimal_grain_size>0) std::cout<<"merge areas with minimal grain size "<<minimal_grain_size<<std::endl;
            if (minimal_bubble_distance>0)
            {
                std::cout<<"merge areas with minimal distance to bubble "<<minimal_bubble_distance<<std::endl;
                std::cout<<"stop merging for grains larger than "<<close_bubble_grain_size<<std::endl;
            }

            //we merge regions until nr of regions stays the same
            size_t old_nr_areas=0;
            size_t diff=1;
            while (diff>0)
            {
                merge_areas(region_image, arcs, dim_x, dim_y, arc_state, twoCellNeighbors, nr_areas, areas, region_arc_index, minimal_grain_size,
                            found_bubble_areas, found_border_areas, close_bubble_areas, not_used, 2);
                //std::cout<<"nr of regions after merging: "<<nr_areas<<std::endl;
                diff=nr_areas-old_nr_areas;
                old_nr_areas=nr_areas;

                //avoid to large areas to be merged
                int i=0;
                while (i<close_bubble_areas.size())
                {
                    if(areas[close_bubble_areas[i]].size()>close_bubble_grain_size)
                    {
                        close_bubble_areas.erase(close_bubble_areas.begin()+i);
                    }
                    else i++;
                }
            }
        }

        if(minimal_grain_size>0 || minimal_bubble_distance>0)
        {
            std::cout<<"nr of grains after merging: "<<nr_areas-found_bubble_areas.size()-found_border_areas.size()<<std::endl;
            update_structures(segment,
                              GrainBoundNet,
                              para,
                              grain_arc_index,
                              minimal_bubble_distance(),
                              close_bubble_grain_size,
                              pixels_average,
                              areas,
                              area_ranges,
                              found_border_areas,
                              arc_state,
                              region_arc_index,
                              twoCellNeighbors,
                              close_bubble_areas,
                              grain_area_junctions,
                              arc_junctions,
                              grain_perimeter,
                              min_bubble_distance,
                              grain_longest_arc_length,
                              grain_junction_angles,
                              grain_junction_angles2,
                              grain_boundary_pixels,
                              region_image,
                              arc_to_segment);

            //update bubble_arc_index
            std::vector< std::vector<int> > old_bubble_arc_index;

            for (int bubble=0; bubble<bubble_arc_index.size(); bubble++)
            {
                if (bubble_arc_index[bubble].size()>0) old_bubble_arc_index.push_back(bubble_arc_index[bubble]);
            }

            bubble_arc_index.clear();
            bubble_arc_index.resize(nr_areas);

            for (int bubble=0; bubble<found_bubble_areas.size(); bubble++)
            {
                for (int arc=0; arc<old_bubble_arc_index[bubble].size(); arc++)
                {
                    bubble_arc_index[found_bubble_areas[bubble]-1].push_back(old_bubble_arc_index[bubble][arc]);
                }
            }

            old_bubble_arc_index.clear();
        }

        //create old grain arc index without segments
        old_grain_arc_index.clear();
        old_grain_arc_index.resize(nr_areas);

        for (int area=0; area<nr_areas; area++)
        {
            for(int segment=0; segment<grain_arc_index[area].size(); segment++)
            {
                for(int a=0; a<grain_arc_index[area][segment].size(); a++)
                {
                    old_grain_arc_index[area].push_back(grain_arc_index[area][segment][a]);
                }
            }

            std::sort(old_grain_arc_index[area].begin(), old_grain_arc_index[area].end());
        }

        //***************************************
        //FILL HISTOGRAMS FOR THIS MIN GRAIN SIZE
        //***************************************

        if (grain_size_min>10) std::cout<<"Exclude grains smaller than "<<grain_size_min<<" pixels"<<std::endl;

        std::vector<size_index> area_sizes;

        for (int area=0; area<nr_areas; area++)
        {
            if (grain_area_size[area]>grain_size_min) //criterion to exclude small grains instead of merging
            {
                size_index grain;
                grain.size=grain_area_size[area]/area_scaling();
                grain.index=area;
                area_sizes.push_back(grain);
            }
        }

        std::sort(area_sizes.begin(),area_sizes.end(),size_index_cmp);

        grain_size_values.clear();
        grain_size_region_values.clear();
        std::vector<int> & grain_areas = para.grain_areas;
        grain_areas.clear();
        percentage_filled_grains_mean.clear();
        percentage_filled_grains_standard_deviation.clear();
        percentage_values.clear();
        nr_percentage_filled_grains.clear();

        for (int area=0; area<area_sizes.size(); area++)
        {
            grain_size_values.push_back(area_sizes[area].size);
            grain_areas.push_back(area_sizes[area].index);//index of grain areas
        }
        area_sizes.clear();

        double accumulated_grain_sizes=0;
        bool seventy_reached=false;
        bool eighty_reached=false;
        bool ninety_reached=false;
        bool nineties_reached[9]={false,false,false,false,false,false,false,false,false};
        bool percentage_reached=false;
        int percentage_nr_grains=0;

        double sum_grain_size_values=0;
        for (int area_number=0; area_number<grain_size_values.size();++area_number)
            sum_grain_size_values+=grain_size_values[area_number];

        for (int area_number=0; area_number<grain_size_values.size();++area_number)
        {
            accumulated_grain_sizes+=grain_size_values[area_number];

            grain_size_region_values.push_back(grain_size_values[area_number]);

            if (!seventy_reached && accumulated_grain_sizes/sum_grain_size_values>0.7)
            {
                float mean, standard_deviation;
                calculate_mean_standard_deviation(grain_size_region_values, mean, standard_deviation);
                percentage_filled_grains_mean.push_back(mean);
                percentage_filled_grains_standard_deviation.push_back(standard_deviation);
                percentage_values.push_back(70.0f);
                nr_percentage_filled_grains.push_back(area_number+1);
                seventy_reached=true;
            }

            if (!eighty_reached && accumulated_grain_sizes/sum_grain_size_values>0.8)
            {
                float mean, standard_deviation;
                calculate_mean_standard_deviation(grain_size_region_values, mean, standard_deviation);
                percentage_filled_grains_mean.push_back(mean);
                percentage_filled_grains_standard_deviation.push_back(standard_deviation);
                percentage_values.push_back(80.0f);
                nr_percentage_filled_grains.push_back(area_number+1);
                eighty_reached=true;
            }

            if (!ninety_reached && accumulated_grain_sizes/sum_grain_size_values>0.9)
            {
                float mean, standard_deviation;
                calculate_mean_standard_deviation(grain_size_region_values, mean, standard_deviation);
                percentage_filled_grains_mean.push_back(mean);
                percentage_filled_grains_standard_deviation.push_back(standard_deviation);
                percentage_values.push_back(90.0f);
                nr_percentage_filled_grains.push_back(area_number+1);
                ninety_reached=true;
            }

            for (int n=0; n<9; n++)
                if (!nineties_reached[n] && accumulated_grain_sizes/sum_grain_size_values>(0.91+n*0.01))
                {
                    float mean, standard_deviation;
                    calculate_mean_standard_deviation(grain_size_region_values, mean, standard_deviation);
                    percentage_filled_grains_mean.push_back(mean);
                    percentage_filled_grains_standard_deviation.push_back(standard_deviation);
                    percentage_values.push_back(91.0f+n*1.0f);
                    nr_percentage_filled_grains.push_back(area_number+1);
                    nineties_reached[n]=true;
                }

            if (!percentage_reached && accumulated_grain_sizes/sum_grain_size_values>percentage_grains)
            {
                percentage_nr_grains=area_number+1;
                percentage_reached=true;
            }
        }

        float mean, standard_deviation;
        calculate_mean_standard_deviation(grain_size_region_values, mean, standard_deviation);
        percentage_filled_grains_mean.push_back(mean);
        percentage_filled_grains_standard_deviation.push_back(standard_deviation);
        percentage_values.push_back(100.0f);
        nr_percentage_filled_grains.push_back(grain_size_region_values.size());

        if (percentage_grains<1.0f)
        {
            std::cout<<"Parameter percentage_grains: "<<percentage_grains()<<std::endl;
            std::cout<<"Only "<<percentage_nr_grains<<" of "<<grain_size_values.size()<<" grains are considered."<<std::endl;

            //reduce grain_size_values and grain areas
            grain_size_values.resize(percentage_nr_grains);
            grain_areas.resize(percentage_nr_grains);
        }

        std::sort(grain_areas.begin(),grain_areas.end());

        //fill grain size, equivalent radius, roundness histogram, width/height and flattening histograms
        grain_size_max=0;
        grain_size_histogram.clear();

        grain_roundness_max=0.0f;
        grain_roundness_histogram.clear();
        //allocation is necessary if no min_value is set, i.e. 0 as value is allowed
        grain_roundness_histogram.resize((int)(grain_roundness_max/grain_roundness_bin_width)+1);
        grain_roundness_values.clear();

        grain_equiv_radius_max=0.0f;
        grain_equiv_radius_histogram.clear();
        grain_equiv_radius_values.clear();

        grain_equiv_radius_norm_max=0.0f;
        grain_equiv_radius_norm_histogram.clear();
        grain_equiv_radius_norm_histogram.resize((int)(grain_equiv_radius_norm_max/grain_equiv_radius_norm_bin_width)+1);
        grain_equiv_radius_values.clear();

        grain_area_width_max=0;
        grain_area_width_histogram.clear();
        grain_area_width_values.clear();

        grain_area_height_max=0;
        grain_area_height_histogram.clear();
        grain_area_height_values.clear();

        grain_area_flattening_max=0;
        grain_area_flattening_histogram.clear();
        grain_area_flattening_values.clear();
        
        for (int area_index=0; area_index<grain_areas.size(); area_index++)//only selected grains
        {
            int area=grain_areas[area_index];

            if (grain_area_size[area]>grain_size_max)
            {
                grain_size_histogram.resize((int)grain_bin_width*log10(grain_area_size[area])+1);
                grain_size_max=grain_area_size[area];
            }
            //no min criterion required as alredy excluded
            grain_size_histogram[(int)grain_bin_width*log10(grain_area_size[area])]++;

            float grain_area_flattening_value=(float)std::max(1, area_ranges[area].x_high-area_ranges[area].x_low)/
                (float)std::max(1, area_ranges[area].y_high-area_ranges[area].y_low);

            //grain roundness
            float roundness_factor=(float)(4.0f*PI*grain_area_size[area])/(float)(grain_perimeter[area]*grain_perimeter[area]);

            if (roundness_factor>grain_roundness_max)
            {
                grain_roundness_histogram.resize((int)(roundness_factor/grain_roundness_bin_width)+1);
                grain_roundness_max=roundness_factor;
            }
            grain_roundness_histogram[(int)(roundness_factor/grain_roundness_bin_width)]++;
            grain_roundness_values.push_back(roundness_factor);

            //grain equivalent radius
            float equiv_radius=sqrt((float)grain_area_size[area]/PI);

            if (equiv_radius>grain_equiv_radius_max)
            {
                grain_equiv_radius_histogram.resize((int)grain_equiv_radius_bin_width*log10(equiv_radius)+1);
                grain_equiv_radius_max=equiv_radius;
            }
            if (equiv_radius>grain_equiv_radius_min) 
            {
                grain_equiv_radius_histogram[(int)grain_equiv_radius_bin_width*log10(equiv_radius)]++;
            }
            grain_equiv_radius_values.push_back(equiv_radius/length_scaling());

            //grain area width
            float grain_area_width=std::max(1, area_ranges[area].x_high-area_ranges[area].x_low);

            if (grain_area_width>grain_area_width_max)
            {
                grain_area_width_histogram.resize((int)grain_area_width_bin_width*log10(grain_area_width)+1);
                grain_area_width_max=grain_area_width;
            }
            if (grain_area_width>grain_area_width_min) 
            {
                grain_area_width_histogram[(int)grain_area_width_bin_width*log10(grain_area_width)]++;
            }
            grain_area_width_values.push_back(grain_area_width/length_scaling());

            //grain area height
            float grain_area_height=std::max(1, area_ranges[area].y_high-area_ranges[area].y_low);

            if (grain_area_height>grain_area_height_max)
            {
                grain_area_height_histogram.resize((int)grain_area_height_bin_width*log10(grain_area_height)+1);
                grain_area_height_max=grain_area_height;
            }
            if (grain_area_height>grain_area_height_min) 
            {
                grain_area_height_histogram[(int)grain_area_height_bin_width*log10(grain_area_height)]++;
            }
            grain_area_height_values.push_back(grain_area_height/length_scaling());

            //grain area flattening
            float grain_area_flattening=grain_area_width/grain_area_height;

            if (grain_area_flattening>grain_area_flattening_max)
            {
                grain_area_flattening_histogram.resize((int)(grain_area_flattening/grain_area_flattening_bin_width)+1);
                grain_area_flattening_max=grain_area_flattening;
            }
            if (grain_area_flattening>grain_area_flattening_min)
            {
                grain_area_flattening_histogram[(int)(grain_area_flattening/grain_area_flattening_bin_width)]++;
            }
            grain_area_flattening_values.push_back(grain_area_flattening);
        }
        calculate_mean_standard_deviation(grain_size_values, grain_size_mean, grain_size_standard_deviation);
        calculate_mean_standard_deviation(grain_roundness_values, grain_roundness_mean, grain_roundness_standard_deviation);
        calculate_mean_standard_deviation(grain_equiv_radius_values, grain_equiv_radius_mean, grain_equiv_radius_standard_deviation);
        calculate_mean_standard_deviation(grain_area_width_values, grain_area_width_mean, grain_area_width_standard_deviation);
        calculate_mean_standard_deviation(grain_area_height_values, grain_area_height_mean, grain_area_height_standard_deviation);
        calculate_mean_standard_deviation(grain_area_flattening_values, grain_area_flattening_mean, grain_area_flattening_standard_deviation);

        //normalize grain equivalent radius
        for (int area_index=0; area_index<grain_equiv_radius_values.size(); area_index++)
        {
            float normed_values=grain_equiv_radius_values[area_index]/grain_equiv_radius_mean;

            if (normed_values>grain_equiv_radius_norm_max)
            {
                grain_equiv_radius_norm_histogram.resize((int)(normed_values/grain_equiv_radius_norm_bin_width)+1);
                grain_equiv_radius_norm_max=normed_values;
            }
            grain_equiv_radius_norm_histogram[(int)(normed_values/grain_equiv_radius_norm_bin_width)]++;
        }

        //fit ellipse, fill grain ellipse and box flattening histogram
        grain_ellipse_params.clear();
        grain_ellipse_params.resize(grain_size_values.size());

        grain_ellipse_long_axis_max=0;
        grain_ellipse_long_axis_histogram.clear();
        grain_ellipse_long_axis_values.clear();
        grain_ellipse_long_axis_values.resize(grain_size_values.size(), 0.0f);
        grain_ellipse_long_axis_values2.clear();

        grain_ellipse_flattening_max=0.0f;
        grain_ellipse_flattening_histogram.clear();
        grain_ellipse_flattening_values.clear();
        grain_ellipse_flattening_values.resize(grain_size_values.size(), 0.0f);
        grain_ellipse_flattening_values2.clear();

        grain_ellipse_long_axis_angle_max=0.0f;
        grain_ellipse_long_axis_angle_histogram.clear();
        grain_ellipse_long_axis_angle_values.clear();
        grain_ellipse_long_axis_angle_values.resize(grain_size_values.size(), 0.0f);
        grain_ellipse_long_axis_angle_values2.clear();

        grain_box_width_max=0;
        grain_box_width_histogram.clear();
        grain_box_width_values.clear();
        grain_box_width_values.resize(grain_size_values.size(), 0.0f);
        grain_box_width_values2.clear();

        grain_box_height_max=0;
        grain_box_height_histogram.clear();
        grain_box_height_values.clear();
        grain_box_height_values.resize(grain_size_values.size(), 0.0f);
        grain_box_height_values2.clear();

        grain_box_flattening_max=0.0f;
        grain_box_flattening_histogram.clear();
        grain_box_flattening_values.clear();
        grain_box_flattening_values.resize(grain_size_values.size(), 0.0f);
        grain_box_flattening_values2.clear();

        grain_perimeter_ratio_max=0.0f;
        grain_perimeter_ratio_histogram.clear();
        grain_perimeter_ratio_values.clear();
        grain_perimeter_ratio_values.resize(grain_size_values.size(), 0.0f);

        for (int area_index=0; area_index<grain_areas.size(); area_index++)//only selected grains
        {
            //grain boundary pixels
            std::vector<point> all_boundary_pixels;
            all_boundary_pixels.clear();

            for(int arc=0; arc<old_grain_arc_index[grain_areas[area_index]].size(); arc++)
            {
                for(int p=0; p<arcs[old_grain_arc_index[grain_areas[area_index]][arc]].size(); p++)
                {
                    all_boundary_pixels.push_back(arcs[old_grain_arc_index[grain_areas[area_index]][arc]][p]);
                }
            }

            //fit ellipse to grain boundary pixels
            std::vector<float> pvec;
            fitEllipse.fit(all_boundary_pixels,pvec);
            grain_ellipse_params[area_index]=pvec;

            //Fit convex perimeter to grain boundary pixels
            float perimeterLength = convPerimeter.fit(all_boundary_pixels);

            if(fabs(1.0f-grain_ellipse_params[area_index][5])<0.1f)//ellipse fitting gave result
            {
                float axis_1=2.0f*sqrt((2.0f*(grain_ellipse_params[area_index][0]*(grain_ellipse_params[area_index][4]/2.0f)*
                    (grain_ellipse_params[area_index][4]/2.0f)+grain_ellipse_params[area_index][2]*
                    (grain_ellipse_params[area_index][3]/2.0f)*(grain_ellipse_params[area_index][3]/2.0f)+
                    grain_ellipse_params[area_index][5]*(grain_ellipse_params[area_index][1]/2.0f)*
                    (grain_ellipse_params[area_index][1]/2.0f)-2.0f*(grain_ellipse_params[area_index][1]/2.0f)*
                    (grain_ellipse_params[area_index][3]/2.0f)*(grain_ellipse_params[area_index][4]/2.0f)-
                    grain_ellipse_params[area_index][0]*grain_ellipse_params[area_index][2]*
                    grain_ellipse_params[area_index][5]))/(((grain_ellipse_params[area_index][1]/2.0f)*
                    (grain_ellipse_params[area_index][1]/2.0f)-grain_ellipse_params[area_index][0]*
                    grain_ellipse_params[area_index][2])*(-sqrt((grain_ellipse_params[area_index][0]-
                    grain_ellipse_params[area_index][2])*(grain_ellipse_params[area_index][0]-
                    grain_ellipse_params[area_index][2])+4.0f*(grain_ellipse_params[area_index][1]/2.0f)*
                    (grain_ellipse_params[area_index][1]/2.0f))-(grain_ellipse_params[area_index][0]+
                    grain_ellipse_params[area_index][2]))));

                float axis_2=2.0f*sqrt((2.0f*(grain_ellipse_params[area_index][0]*(grain_ellipse_params[area_index][4]/2.0f)*
                    (grain_ellipse_params[area_index][4]/2.0f)+grain_ellipse_params[area_index][2]*
                    (grain_ellipse_params[area_index][3]/2.0f)*(grain_ellipse_params[area_index][3]/2.0f)+
                    grain_ellipse_params[area_index][5]*(grain_ellipse_params[area_index][1]/2.0f)*
                    (grain_ellipse_params[area_index][1]/2.0f)-2.0f*(grain_ellipse_params[area_index][1]/2.0f)*
                    (grain_ellipse_params[area_index][3]/2.0f)*(grain_ellipse_params[area_index][4]/2.0f)-
                    grain_ellipse_params[area_index][0]*grain_ellipse_params[area_index][2]*
                    grain_ellipse_params[area_index][5]))/(((grain_ellipse_params[area_index][1]/2.0f)*
                    (grain_ellipse_params[area_index][1]/2.0f)-grain_ellipse_params[area_index][0]*
                    grain_ellipse_params[area_index][2])*(sqrt((grain_ellipse_params[area_index][0]-
                    grain_ellipse_params[area_index][2])*(grain_ellipse_params[area_index][0]-
                    grain_ellipse_params[area_index][2])+4.0f*(grain_ellipse_params[area_index][1]/2.0f)*
                    (grain_ellipse_params[area_index][1]/2.0f))-(grain_ellipse_params[area_index][0]+
                    grain_ellipse_params[area_index][2]))));

                float long_axis=std::max(axis_1,axis_2);
                float short_axis=std::min(axis_1,axis_2);

                //grain ellipse long axis
                if (long_axis>grain_ellipse_long_axis_max)
                {
                    grain_ellipse_long_axis_histogram.resize((int)grain_ellipse_long_axis_bin_width*log10(long_axis)+1);
                    grain_ellipse_long_axis_max=long_axis;
                }
                if (long_axis>grain_ellipse_long_axis_min) 
                {
                    grain_ellipse_long_axis_histogram[(int)grain_ellipse_long_axis_bin_width*log10(long_axis)]++;
                    grain_ellipse_long_axis_values2.push_back(long_axis/length_scaling());
                }
                grain_ellipse_long_axis_values[area_index]=long_axis/length_scaling();

                //grain ellipse flattening
                if (long_axis>0.0f && short_axis>0.0f) grain_ellipse_flattening_values[area_index]=long_axis/short_axis;//axes found correctly

                if (grain_ellipse_flattening_values[area_index]>grain_ellipse_flattening_max)
                {
                    grain_ellipse_flattening_histogram.resize((int)(grain_ellipse_flattening_values[area_index]/grain_ellipse_flattening_bin_width)+1);
                    grain_ellipse_flattening_max=grain_ellipse_flattening_values[area_index];
                }
                if(grain_ellipse_flattening_values[area_index]>grain_ellipse_flattening_min)//axes found correctly
                {
                    grain_ellipse_flattening_histogram[(int)(grain_ellipse_flattening_values[area_index]/grain_ellipse_flattening_bin_width)]++;
                    grain_ellipse_flattening_values2.push_back(grain_ellipse_flattening_values[area_index]);
                }

                if(grain_ellipse_params[area_index][0]!=grain_ellipse_params[area_index][2])//ellipse is not degenerated
                {   
                    float angle=0.0f;

                    if (grain_ellipse_params[area_index][1]==0.0f &&
                        grain_ellipse_params[area_index][0]<grain_ellipse_params[area_index][2]) angle=0.0f;
                    else if (grain_ellipse_params[area_index][1]==0.0f &&
                        grain_ellipse_params[area_index][0]>grain_ellipse_params[area_index][2]) angle=PI/2.0f;
                    else if (grain_ellipse_params[area_index][1]!=0.0f &&
                        grain_ellipse_params[area_index][0]<grain_ellipse_params[area_index][2])
                        angle=0.5f*atan(grain_ellipse_params[area_index][1]/(grain_ellipse_params[area_index][0]-
                            grain_ellipse_params[area_index][2]));
                    else if (grain_ellipse_params[area_index][1]!=0.0f &&
                        grain_ellipse_params[area_index][0]>grain_ellipse_params[area_index][2])
                        angle=PI/2.0f + 0.5f*atan(grain_ellipse_params[area_index][1]/(grain_ellipse_params[area_index][0]-
                            grain_ellipse_params[area_index][2]));
                    angle=PI-angle;
                    if (angle>PI) angle-=PI;
                    if (angle>0.5f*PI) angle-=PI;

                    //grain ellipse long axis angle
                    if (angle+0.5f*PI>grain_ellipse_long_axis_angle_max)
                    {
                        grain_ellipse_long_axis_angle_histogram.resize((int)((angle+0.5f*PI)/grain_ellipse_long_axis_angle_bin_width)+1);
                        grain_ellipse_long_axis_angle_max=angle+0.5f*PI;
                    }
                    grain_ellipse_long_axis_angle_histogram[(int)((angle+0.5f*PI)/grain_ellipse_long_axis_angle_bin_width)]++;
                    grain_ellipse_long_axis_angle_values[area_index]=angle*180.0f/PI;
                    grain_ellipse_long_axis_angle_values2.push_back(angle*180.0f/PI);

                    area_range box_range;
                    box_range.x_low=dim_x;
                    box_range.y_low=dim_y;
                    box_range.x_high=0;
                    box_range.y_high=0;

                    for(int p=0; p<all_boundary_pixels.size(); p++)
                    {
                        point pixel;
                        pixel.x=all_boundary_pixels[p].x-grain_area_center_mass[grain_areas[area_index]].x;
                        pixel.y=all_boundary_pixels[p].y-grain_area_center_mass[grain_areas[area_index]].y;

                        int x=grain_area_center_mass[grain_areas[area_index]].x+pixel.x*cos(angle)-pixel.y*sin(angle);
                        int y=grain_area_center_mass[grain_areas[area_index]].y+pixel.x*sin(angle)+pixel.y*cos(angle);

                        //find ranges
                        if(x<box_range.x_low) box_range.x_low=x;
                        if(y<box_range.y_low) box_range.y_low=y;
                        if(x>box_range.x_high) box_range.x_high=x;
                        if(y>box_range.y_high) box_range.y_high=y;
                    }

                    if (box_range.x_high-box_range.x_low<box_range.y_high-box_range.y_low)
                    {
                        int temp_low=box_range.x_low;
                        int temp_high=box_range.x_high;
                        box_range.x_low=box_range.y_low;
                        box_range.x_high=box_range.y_high;
                        box_range.y_low=temp_low;
                        box_range.y_high=temp_high;
                    }

                    //grain box width
                    grain_box_width_values[area_index]=box_range.x_high-box_range.x_low;

                    if (grain_box_width_values[area_index]>grain_box_width_max)
                    {
                        grain_box_width_histogram.resize((int)grain_box_width_bin_width*log10(grain_box_width_values[area_index])+1);
                        grain_box_width_max=grain_box_width_values[area_index];
                    }
                    grain_box_width_histogram[(int)grain_box_width_bin_width*log10(grain_box_width_values[area_index])]++;
                    grain_box_width_values[area_index]/=length_scaling();
                    grain_box_width_values2.push_back(grain_box_width_values[area_index]);

                    //grain box height
                    grain_box_height_values[area_index]=box_range.y_high-box_range.y_low;

                    if (grain_box_height_values[area_index]>grain_box_height_max)
                    {
                        grain_box_height_histogram.resize((int)grain_box_height_bin_width*log10(grain_box_height_values[area_index])+1);
                        grain_box_height_max=grain_box_height_values[area_index];
                    }
                    grain_box_height_histogram[(int)grain_box_height_bin_width*log10(grain_box_height_values[area_index])]++;
                    grain_box_height_values[area_index]/=length_scaling();
                    grain_box_height_values2.push_back(grain_box_height_values[area_index]);

                    //grain box flattening
                    grain_box_flattening_values[area_index]=grain_box_width_values[area_index]/grain_box_height_values[area_index];

                    if (grain_box_flattening_values[area_index]>grain_box_flattening_max)
                    {
                        grain_box_flattening_histogram.resize((int)(grain_box_width_values[area_index]/
                            (grain_box_height_values[area_index]*grain_box_flattening_bin_width))+1);
                        grain_box_flattening_max=grain_box_flattening_values[area_index];
                    }
                    grain_box_flattening_histogram[(int)(grain_box_width_values[area_index]/(grain_box_height_values[area_index]*grain_box_flattening_bin_width))]++;
                    grain_box_flattening_values2.push_back(grain_box_flattening_values[area_index]);

                    //Grain perimeter ratio
                    grain_perimeter_ratio_values[area_index] = perimeterLength/(float)grain_perimeter[grain_areas[area_index]];                
                    if(grain_perimeter_ratio_values[area_index] > grain_perimeter_ratio_max)
                    {                        
                        grain_perimeter_ratio_histogram.resize((int)(perimeterLength/((float)grain_perimeter[grain_areas[area_index]]*grain_perimeter_ratio_bin_width))+1);
                        grain_perimeter_ratio_max = grain_perimeter_ratio_values[area_index];
                    }
                    grain_perimeter_ratio_histogram[(int)(perimeterLength/((float)grain_perimeter[grain_areas[area_index]]*grain_perimeter_ratio_bin_width))]++;                
                }
            }
        }

        calculate_mean_standard_deviation(grain_box_flattening_values2, grain_box_flattening_mean, grain_box_flattening_standard_deviation);
        calculate_mean_standard_deviation(grain_box_width_values2, grain_box_width_mean, grain_box_width_standard_deviation);
        calculate_mean_standard_deviation(grain_box_height_values2, grain_box_height_mean, grain_box_height_standard_deviation);
        calculate_mean_standard_deviation(grain_ellipse_long_axis_values2, grain_ellipse_long_axis_mean, grain_ellipse_long_axis_standard_deviation);
        calculate_mean_standard_deviation(grain_ellipse_flattening_values2, grain_ellipse_flattening_mean, grain_ellipse_flattening_standard_deviation);
        calculate_mean_standard_deviation(grain_ellipse_long_axis_angle_values2, grain_ellipse_long_axis_angle_mean,
                                          grain_ellipse_long_axis_angle_standard_deviation);
        calculate_mean_standard_deviation(grain_perimeter_ratio_values, grain_perimeter_ratio_mean, grain_perimeter_ratio_standard_deviation);

        //calculate correlations
        float covariance_ellipse_box_flattening=0.0f;
        int nr_values=0;

        for (int area_index=0; area_index<grain_areas.size(); area_index++)
        {
            if (grain_box_flattening_values[area_index]>0.0f && grain_ellipse_flattening_values[area_index]>0.0f)//ellipse not degenerated and exes found correctly
            {
                covariance_ellipse_box_flattening+=(grain_box_flattening_values[area_index]-grain_box_flattening_mean)*
                    (grain_ellipse_flattening_values[area_index]-grain_ellipse_flattening_mean);
                    nr_values++;
            }
        }
        if(grain_areas.size() > 0)
        {
            correlation_ellipse_box_flattening=covariance_ellipse_box_flattening/((float)nr_values*grain_box_flattening_standard_deviation*
                grain_ellipse_flattening_standard_deviation);
        }

        std::cout<<"Correlation ellipse box flattening: "<<correlation_ellipse_box_flattening<<std::endl;

        //fill grain region/percent/quantile plot
        grain_size_region_values.clear();
        grain_region_values.clear();
        grain_region_errors.clear();
        grain_percent_values.clear();
        grain_percent_errors.clear();
        grain_quantile_values.clear();
        grain_quantile_errors.clear();

        int quantile_step_rounded = grain_size_values.size()/100;
        int percent_step_rounded = grain_size_values.size()/10;

        for (int area_number=0; area_number<grain_size_values.size();++area_number)
        {
            grain_size_region_values.push_back(grain_size_values[area_number]);

            if ((area_number+1)%grain_step()==0)
            {
                float grain_size_region_mean;
                float grain_size_region_standard_deviation;
                calculate_mean_standard_deviation(grain_size_region_values, grain_size_region_mean, grain_size_region_standard_deviation);
                grain_region_values.push_back(grain_size_region_mean);
                grain_region_errors.push_back(grain_size_region_standard_deviation);
            }
            if (quantile_step_rounded>0)
                if((area_number+1)%quantile_step_rounded==0 && grain_quantile_values.size()<10)
                {
                    grain_quantile_values.push_back(grain_size_values[area_number]);
                    //grain_quantile_errors.push_back((float)(grain_size_values[area_number-1]-
                    //    grain_size_values[area_number+1])*(float)grain_size_values[area_number]/2.0f);
                }
            if (percent_step_rounded>0)
                if ((area_number+1)%percent_step_rounded==0)
                {
                    float grain_size_region_mean;
                    float grain_size_region_standard_deviation;
                    calculate_mean_standard_deviation(grain_size_region_values, grain_size_region_mean, grain_size_region_standard_deviation);
                    grain_percent_values.push_back(grain_size_region_mean);
                    grain_percent_errors.push_back(grain_size_region_standard_deviation);
                }
        }

        //fill bubble distance histogram and relaxation plots
        min_bubble_distance_max=0;

        for(int i=0; i<min_bubble_distance_histogram.size(); i++)
        {
            min_bubble_distance_histogram[i].clear();
            min_bubble_distance_values[i].clear();
        }

        relaxation_values.clear();
        relaxation2_values.clear();

        relaxation_values.resize((int)grain_bin_width*log10(grain_size_max)+1);
        relaxation_mean.resize((int)grain_bin_width*log10(grain_size_max)+1);
        relaxation_standard_deviation.resize((int)grain_bin_width*log10(grain_size_max)+1);

        float size_regions[nr_size_regions];
        for(int i=0; i<nr_size_regions; i++) size_regions[i]=0.0f;

        for (int i=0; i<grain_size_values.size();++i)
        {
            int index=nr_size_regions-1-(nr_size_regions*i/(int)grain_size_values.size());
            if(size_regions[index]==0.0f) size_regions[index]=grain_size_values[i];
        }

        for (int area_index=0; area_index<grain_areas.size(); area_index++)//only grains considered for previous histograms
        {
            if (min_bubble_distance[grain_areas[area_index]]>min_bubble_distance_max)
            {
                int new_size=min_bubble_distance_bin_width*log10(min_bubble_distance[grain_areas[area_index]])+1;

                for(int i=0; i<min_bubble_distance_histogram.size(); i++) min_bubble_distance_histogram[i].resize(new_size);
                min_bubble_distance_max=min_bubble_distance[grain_areas[area_index]];

                relaxation2_values.resize(new_size);
                relaxation2_mean.resize(new_size);
                relaxation2_standard_deviation.resize(new_size);
            }
            for(int i=0; i<min_bubble_distance_histogram.size(); i++)
            {
                if (min_bubble_distance[grain_areas[area_index]]>min_bubble_distance_min &&
                    grain_area_size[grain_areas[area_index]]/area_scaling()<size_regions[i]) 
                {
                    min_bubble_distance_histogram[i][(int)min_bubble_distance_bin_width*log10(min_bubble_distance[grain_areas[area_index]])]++;
                    min_bubble_distance_values[i].push_back(min_bubble_distance[grain_areas[area_index]]/length_scaling());

                    relaxation_values[(int)grain_bin_width*log10(grain_area_size[grain_areas[area_index]])].push_back
                        (min_bubble_distance_bin_width*log10(min_bubble_distance[grain_areas[area_index]]));

                    relaxation2_values[(int)min_bubble_distance_bin_width*log10(min_bubble_distance[grain_areas[area_index]])].push_back
                        (grain_bin_width*log10(grain_area_size[grain_areas[area_index]]));
                }
            }
        }

        for(int i=0; i<nr_size_regions; i++)
        {
            calculate_mean_standard_deviation(min_bubble_distance_values[i], min_bubble_distance_mean[i],
                min_bubble_distance_standard_deviation[i]);
        }
        for(int i=0; i<relaxation_values.size(); i++)
        {
            calculate_mean_standard_deviation(relaxation_values[i], relaxation_mean[i], relaxation_standard_deviation[i]);
        }
        for(int i=0; i<relaxation2_values.size(); i++)
        {
            calculate_mean_standard_deviation(relaxation2_values[i], relaxation2_mean[i], relaxation2_standard_deviation[i]);
        }

        //fill grain arc length histogram
        grain_arc_length_max=0;
        grain_arc_length_histogram.clear();
        grain_arc_length_values.clear();
        std::vector<size_index> boundary_sizes;

        for (int segment=0; segment<grain_boundary_pixels.size(); segment++)
        {
            size_index bound;
            bound.size=grain_boundary_pixels[segment].size()/length_scaling();
            bound.index=segment;
            boundary_sizes.push_back(bound);

            bool border_contact=false;
            int arc_index=fabs(grain_boundary_index[segment][0])-1;

            if ((grain_area_size[twoCellNeighbors(arc_index,1)-1]==0 && bubble_area_size[twoCellNeighbors(arc_index,1)-1]==0) ||
                (grain_area_size[twoCellNeighbors(arc_index,0)-1]==0 && bubble_area_size[twoCellNeighbors(arc_index,0)-1]==0))
                border_contact=true;

            int grain_arc_length;

            if (border_contact)//grain boundaries with border on one side can not handled correctly
            {
                grain_arc_length=0;
            }
            else grain_arc_length=get_length(grain_boundary_pixels[segment]);

            if (grain_arc_length>grain_arc_length_max)
            {
                grain_arc_length_histogram.resize((int)grain_arc_length_bin_width*log10(grain_arc_length)+1);
                grain_arc_length_max=grain_arc_length;
            }
            if (grain_arc_length>grain_arc_length_min) 
            {
                grain_arc_length_histogram[(int)grain_arc_length_bin_width*log10(grain_arc_length)]++;
                grain_arc_length_values.push_back(grain_arc_length/length_scaling());
            }
        }

        calculate_mean_standard_deviation(grain_arc_length_values, grain_arc_length_mean, grain_arc_length_standard_deviation);

        //fill grain arc number, nr of neighbors histogram, grain longest arc length and dihedral angle histogram
        grain_arc_number_max=0;
        grain_arc_number_histogram.clear();
        grain_arc_number_values.clear();
        grain_arc_number_values2.clear();

        grain_neighbors_max=0;
        grain_neighbors_histogram.clear();
        grain_neighbors_values.clear();
        grain_neighbors_values2.clear();

        grain_longest_arc_length_max=0;
        grain_longest_arc_length_histogram.clear();
        grain_longest_arc_length_values.clear();
        grain_longest_arc_length_values2.clear();

        dihedral_angle_max=0.0f;
        dihedral_angle_histogram.clear();
        dihedral_angle_histogram.resize((int)(dihedral_angle_max/dihedral_angle_bin_width)+1);
        dihedral_angle_values.clear();

        dihedral_angle2_max=0.0f;
        dihedral_angle2_histogram.clear();
        dihedral_angle2_histogram.resize((int)(dihedral_angle2_max/dihedral_angle2_bin_width)+1);
        dihedral_angle2_values.clear();

        for (int area_index=0; area_index<grain_areas.size(); area_index++)//only selected grains
        {
            int area=grain_areas[area_index];

            //find all neighbors
            std::list<int> neighbors;
            bool border_contact=false;

            for(int a=0; a<grain_arc_index[area].size() && !border_contact; a++)
            {
                int arc_index=grain_arc_index[area][a][0];

                if(twoCellNeighbors(arc_index,0)==area+1 && grain_area_size[twoCellNeighbors(arc_index,1)-1]>0)
                    neighbors.push_back(twoCellNeighbors(arc_index,1));
                else if(twoCellNeighbors(arc_index,1)==area+1 && grain_area_size[twoCellNeighbors(arc_index,0)-1]>0)
                    neighbors.push_back(twoCellNeighbors(arc_index,0));
                if ((twoCellNeighbors(arc_index,0)==area+1 && grain_area_size[twoCellNeighbors(arc_index,1)-1]==0 &&
                     bubble_area_size[twoCellNeighbors(arc_index,1)-1]==0) || (twoCellNeighbors(arc_index,1)==area+1 &&
                     grain_area_size[twoCellNeighbors(arc_index,0)-1]==0 && bubble_area_size[twoCellNeighbors(arc_index,0)-1]==0))
                    border_contact=true;
            }

            if (border_contact)//grains with contact to border can not handled correctly
            {
                neighbors.clear();
                grain_longest_arc_length[area]=0;
            }

            int nr_arcs=(int)neighbors.size();
            neighbors.sort();
            neighbors.unique();
            int nr_neighbors=(int)neighbors.size();

            //nr of arcs
            if (nr_arcs>grain_arc_number_max)
            {
                grain_arc_number_histogram.resize((int)(nr_arcs/grain_arc_number_bin_width)+1);
                grain_arc_number_max=nr_arcs;
            }
            if (nr_arcs>grain_arc_number_min) 
            {
                grain_arc_number_histogram[(int)(nr_arcs/grain_arc_number_bin_width)]++;
                grain_arc_number_values2.push_back(nr_arcs);
            }

            //nr of neighbors
            if (nr_neighbors>grain_neighbors_max)
            {
                grain_neighbors_histogram.resize((int)(nr_neighbors/grain_neighbors_bin_width)+1);
                grain_neighbors_max=nr_neighbors;
            }
            if (nr_neighbors>grain_neighbors_min)
            {
                grain_neighbors_histogram[(int)(nr_neighbors/grain_neighbors_bin_width)]++;
                grain_neighbors_values2.push_back(nr_neighbors);
            }

            //longest arc length
            if (grain_longest_arc_length[area]>grain_longest_arc_length_max)
            {
                grain_longest_arc_length_histogram.resize((int)grain_longest_arc_length_bin_width*log10(grain_longest_arc_length[area])+1);
                grain_longest_arc_length_max=grain_longest_arc_length[area];
            }
            if (grain_longest_arc_length[area]>grain_longest_arc_length_min) 
            {
                grain_longest_arc_length_histogram[(int)grain_longest_arc_length_bin_width*log10(grain_longest_arc_length[area])]++;
                grain_longest_arc_length_values2.push_back(grain_longest_arc_length[area]/length_scaling());
            }

            //dihedral angles
            for (int junction=0; junction<grain_junction_angles[area].size(); junction++)
            {
                if (grain_junction_angles[area][junction]>dihedral_angle_max)
                {
                    dihedral_angle_histogram.resize((int)(grain_junction_angles[area][junction]/dihedral_angle_bin_width)+1);
                    dihedral_angle_max=grain_junction_angles[area][junction];
                }

                dihedral_angle_histogram[(int)(grain_junction_angles[area][junction]/dihedral_angle_bin_width)]++;
                dihedral_angle_values.push_back(grain_junction_angles[area][junction]/angle_scaling);
            }

            //dihedral angles 2
            for (int junction=0; junction<grain_junction_angles2[area].size(); junction++)
            {
                if (grain_junction_angles2[area][junction]>dihedral_angle2_max)
                {
                    dihedral_angle2_histogram.resize((int)(grain_junction_angles2[area][junction]/dihedral_angle2_bin_width)+1);
                    dihedral_angle2_max=grain_junction_angles2[area][junction];
                }

                dihedral_angle2_histogram[(int)(grain_junction_angles2[area][junction]/dihedral_angle2_bin_width)]++;
                dihedral_angle2_values.push_back(grain_junction_angles2[area][junction]/angle_scaling);
            }


            //values for hdf5 parameter file
            if (grain_area_size[area]>grain_size_min)
            {
                grain_arc_number_values.push_back(nr_arcs);
                grain_neighbors_values.push_back(nr_neighbors);
                grain_longest_arc_length_values.push_back(grain_longest_arc_length[area]/length_scaling());
            }

        }
        calculate_mean_standard_deviation(grain_longest_arc_length_values2, grain_longest_arc_length_mean, grain_longest_arc_length_standard_deviation);
        calculate_mean_standard_deviation(dihedral_angle_values, dihedral_angle_mean, dihedral_angle_standard_deviation);
        calculate_mean_standard_deviation(dihedral_angle2_values, dihedral_angle2_mean, dihedral_angle2_standard_deviation);

        //fill dislocation density and turning point histogram
        curvature_max=0.0f;
        dislocation_density_histogram.clear();
        dislocation_density_values.clear();
        float dislocation_in_histogram=0.0f;

        phi_max=0.0f;
        boundary_orientation_histogram.clear();
        boundary_orientation_values.clear();

        turning_point_max=0.0f;
        turning_point_histogram.clear();
        turning_point_histogram.resize((int)turning_point_max+1);
        turning_point_values.clear();

        disl_dens_boundary_region_values.clear();
        disl_dens_boundary_region_errors.clear();
        disl_dens_boundary_percent_values.clear();
        disl_dens_boundary_percent_errors.clear();

        std::sort(boundary_sizes.begin(),boundary_sizes.end(),size_index_cmp);

        percent_step_rounded = boundary_sizes.size()/10;

        for (int boundary_number=0; boundary_number<boundary_sizes.size();++boundary_number)
        {
            int old_pos=-2;
            int not_pos_interval=0;
            bool first_non_zero_found=false;
            bool first_positive=false, last_positive=false;

            for (int p=0; p<grain_boundary_curvs[boundary_sizes[boundary_number].index].size(); p++)
            {
                float curv=fabs(grain_boundary_curvs[boundary_sizes[boundary_number].index][p])/0.75f;//gauging, has to be validated!!

                if (curv>curvature_max && curv>curvature_min && curv<0.15)//CUT OFF BY CURV 0.15!!
                {
                    dislocation_density_histogram.resize((int)(dislocation_density_bin_width*log10(curv/curvature_min))+1);
                    curvature_max=curv;
                }
                if (curv<0.15)
                {
                    //dislocation density histogram is a curvature histogram
                    if (curv>curvature_min)
                    {
                        dislocation_density_histogram[(int)(dislocation_density_bin_width*log10(curv/curvature_min))]++;
                        dislocation_in_histogram++;
                    }
                    dislocation_density_values.push_back(curv/dislocation_density_scaling);
                }

                curv=grain_boundary_curvs[boundary_sizes[boundary_number].index][p];

                if(curv<=0.0f)//find nr of intervals with not positive curvature
                {
                    if(p-old_pos>1) not_pos_interval++;
                    old_pos=p;
                }

                if(!first_non_zero_found && fabs(curv)>0.001f)//check whether first non zero value is positive or negative
                {
                    if(curv<0.0f) first_positive=false;
                    else first_positive=true;
                    first_non_zero_found=true;
                }

                //at the end last_positive determines whether the last non zero value is positive or negative
                if(curv<-0.001f) last_positive=false;
                else if(curv>0.001f) last_positive=true;
            }

            if(grain_boundary_curvs[boundary_sizes[boundary_number].index].size()<5)
            {
                turning_point_values.push_back(0.0f);
            }
            else if(last_positive)
            {
                if(first_positive) turning_point_values.push_back(2*not_pos_interval-4);
                else turning_point_values.push_back(2*not_pos_interval-3);
            }
            else
            {
                if(first_positive) turning_point_values.push_back(2*not_pos_interval-3);
                else turning_point_values.push_back(2*not_pos_interval-2);
            }

            //fill turning point histogram
            if (turning_point_values.back()>turning_point_max)
            {
                turning_point_histogram.resize((int)turning_point_values.back()+1);
                turning_point_max=turning_point_values.back();
            }

            if(turning_point_values.back()<0) std::cout<<"Error: Number of turning points is negative!"<<boundary_sizes[boundary_number].index<<std::endl;
            turning_point_histogram[std::max(0.0f,turning_point_values.back())]++;

            for (int p=0; p<grain_boundary_phis[boundary_sizes[boundary_number].index].size(); p++)
            {
                float phi=grain_boundary_phis[boundary_sizes[boundary_number].index][p];
                while (phi>0.5f*PI) phi-=PI;
                while (phi<-0.5f*PI) phi+=PI;
                phi=-phi;

                if (phi+0.5f*PI>phi_max)
                {
                    boundary_orientation_histogram.resize((int)((phi+0.5f*PI)/boundary_orientation_bin_width)+1);
                    phi_max=phi+0.5f*PI;
                }

                //boundary orientation histogram is a phi histogram
                boundary_orientation_histogram[(int)((phi+0.5f*PI)/boundary_orientation_bin_width)]++;
                boundary_orientation_values.push_back(phi*180.0f/PI);
            }

            //fill dislocation density region/precent plots
            if ((boundary_number+1)%boundary_step()==0)
            {
                float mean;
                float standard_deviation;
                calculate_mean_standard_deviation(dislocation_density_values, mean, standard_deviation);
                disl_dens_boundary_region_values.push_back(mean);
                disl_dens_boundary_region_errors.push_back(standard_deviation);
            }
            if (percent_step_rounded>0)
                if ((boundary_number+1)%percent_step_rounded==0)
                {
                    float mean;
                    float standard_deviation;
                    calculate_mean_standard_deviation(dislocation_density_values, mean, standard_deviation);
                    disl_dens_boundary_percent_values.push_back(mean);
                    disl_dens_boundary_percent_errors.push_back(standard_deviation);
                }
        }
        calculate_mean_standard_deviation(dislocation_density_values, dislocation_density_mean, dislocation_density_standard_deviation);
        dislocation_in_histogram=dislocation_in_histogram/(float)dislocation_density_values.size();
        calculate_mean_standard_deviation(turning_point_values, turning_point_mean, turning_point_standard_deviation);
        calculate_mean_standard_deviation(boundary_orientation_values, boundary_orientation_mean, boundary_orientation_standard_deviation);

        std::stringstream s;
        if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
        else s << "." << minimal_grain_size;
        if (grain_size_min>10) s << "_" << grain_size_min;
        s << ".svg";

        std::string filepath_largest_grains=filepath_plots;
        std::string filepath_percent_grains=filepath_plots;
        std::string filepath_quantile_grains=filepath_plots;
        std::string filepath_grain_size=filepath_plots;
        std::string filepath_grain_size_fit=filepath_plots;
        std::string filepath_grain_roundness=filepath_plots;
        std::string filepath_grain_equiv_radius=filepath_plots;
        std::string filepath_grain_equiv_radius_norm=filepath_plots;
        std::string filepath_grain_box_flattening=filepath_plots;
        std::string filepath_grain_box_width=filepath_plots;
        std::string filepath_grain_box_height=filepath_plots;
        std::string filepath_grain_ellipse_long_axis=filepath_plots;
        std::string filepath_grain_ellipse_flattening=filepath_plots;
        std::string filepath_grain_ellipse_long_axis_angle=filepath_plots;
        std::string filepath_grain_width=filepath_plots;
        std::string filepath_grain_height=filepath_plots;
        std::string filepath_grain_flattening=filepath_plots;
        std::string filepath_grainperimeter_ratio=filepath_plots;
        std::string filepath_min_bubble_distance=filepath_plots;
        std::string filepath_relaxation=filepath_plots;
        std::string filepath_relaxation2=filepath_plots;
        std::string filepath_grain_arc_number=filepath_plots;
        std::string filepath_grain_neighbors=filepath_plots;
        std::string filepath_grain_arc_length=filepath_plots;
        std::string filepath_grain_longest_arc_length=filepath_plots;
        std::string filepath_dislocation_density=filepath_plots;
        std::string filepath_boundary_orientation=filepath_plots;
        std::string filepath_dihedral_angle=filepath_plots;
        std::string filepath_dihedral_angle2=filepath_plots;
        std::string filepath_turning_point=filepath_plots;
        std::string filepath_disl_dens_longest_boundaries=filepath_plots;
        std::string filepath_disl_dens_percent_boundaries=filepath_plots;
        std::string filepath_percentage_filled_grains=filepath_plots;       

        filepath_largest_grains.append(".largest_grains");
        filepath_percent_grains.append(".percent_grains");
        filepath_quantile_grains.append(".quantile_grains");
        filepath_grain_size.append(".grain_size");
        filepath_grain_size_fit.append(".grain_size_fit");
        filepath_grain_roundness.append(".grain_roundness");
        filepath_grain_equiv_radius.append(".grain_equiv_radius");
        filepath_grain_equiv_radius_norm.append(".grain_equiv_radius_norm");
        filepath_grain_box_flattening.append(".grain_box_flattening");
        filepath_grain_box_width.append(".grain_box_width");
        filepath_grain_box_height.append(".grain_box_height");
        filepath_grain_ellipse_long_axis.append(".grain_ellipse_long_axis");
        filepath_grain_ellipse_flattening.append(".grain_ellipse_flattening");
        filepath_grain_ellipse_long_axis_angle.append(".grain_ellipse_long_axis_angle");
        filepath_grain_width.append(".grain_width");
        filepath_grain_height.append(".grain_height");
        filepath_grain_flattening.append(".grain_flattening");
        filepath_grainperimeter_ratio.append(".grain_perimeter_ratio");
        filepath_min_bubble_distance.append(".min_bubble_distance");
        filepath_relaxation.append(".relaxation");
        filepath_relaxation2.append(".relaxation2");
        filepath_grain_arc_number.append(".grain_arc_number");
        filepath_grain_neighbors.append(".grain_neighbors");
        filepath_grain_arc_length.append(".grain_arc_length");
        filepath_grain_longest_arc_length.append(".grain_longest_arc_length");
        filepath_dislocation_density.append(".dislocation_density");
        filepath_boundary_orientation.append(".boundary_orientation");
        filepath_dihedral_angle.append(".dihedral_angle");
        filepath_dihedral_angle2.append(".dihedral_angle2");
        filepath_turning_point.append(".turning_point");
        filepath_disl_dens_longest_boundaries.append(".disl_dens_longest_boundaries");
        filepath_disl_dens_percent_boundaries.append(".disl_dens_percent_boundaries");
        filepath_percentage_filled_grains.append(".percentage_filled_grains");        

        filepath_largest_grains.append(s.str());
        filepath_percent_grains.append(s.str());
        filepath_quantile_grains.append(s.str());
        filepath_grain_size.append(s.str());
        filepath_grain_size_fit.append(s.str());
        filepath_grain_roundness.append(s.str());
        filepath_grain_equiv_radius.append(s.str());
        filepath_grain_equiv_radius_norm.append(s.str());
        filepath_grain_box_flattening.append(s.str());
        filepath_grain_box_width.append(s.str());
        filepath_grain_box_height.append(s.str());
        filepath_grain_ellipse_long_axis.append(s.str());
        filepath_grain_ellipse_flattening.append(s.str());
        filepath_grain_ellipse_long_axis_angle.append(s.str());
        filepath_grain_width.append(s.str());
        filepath_grain_height.append(s.str());
        filepath_grain_flattening.append(s.str());
        filepath_grainperimeter_ratio.append(s.str());
        filepath_min_bubble_distance.append(s.str());
        filepath_relaxation.append(s.str());
        filepath_relaxation2.append(s.str());
        filepath_grain_arc_number.append(s.str());
        filepath_grain_neighbors.append(s.str());
        filepath_grain_arc_length.append(s.str());
        filepath_grain_longest_arc_length.append(s.str());
        filepath_dislocation_density.append(s.str());
        filepath_boundary_orientation.append(s.str());
        filepath_dihedral_angle.append(s.str());
        filepath_dihedral_angle2.append(s.str());
        filepath_turning_point.append(s.str());
        filepath_disl_dens_longest_boundaries.append(s.str());
        filepath_disl_dens_percent_boundaries.append(s.str());
        filepath_percentage_filled_grains.append(s.str());       

        char step[20];
        sprintf(step, "%d", grain_step());

        std::string y_axis="Mean grain size of ";
        y_axis.append(step);
        y_axis.append("*n Largest grains [mm^2]");

        std::string x_axis="Largest grains * ";
        x_axis.append(step);

        char percent_in_histogram[20];
        sprintf(percent_in_histogram, "%1.2f", dislocation_in_histogram);
            std::string dislocation_title="Relative occurrence x";
        dislocation_title.append(percent_in_histogram);

        plot.draw_values_errors(x_axis, y_axis, "Largest grains distribution", grain_region_values, grain_region_errors, filepath_largest_grains.c_str());
        plot.draw_values_errors("* 10 percent largest grains", "Mean grain size of 10*n percent largest grains [mm^2]", "Largest grains distribution",
                                grain_percent_values, grain_percent_errors, filepath_percent_grains.c_str());
        plot.draw_histogram_log(area_unit,"Relative occurrence","Grain size distribution", grain_size_histogram, grain_bin_width, area_scaling(),
                                grain_size_mean, grain_size_standard_deviation, grain_size_y_max, filepath_grain_size.c_str());

        float fit_mu, fit_sigma;

        plot.draw_histogram_left_lognorm(area_unit, "Relative occurrence","Grain size distribution", grain_size_histogram, grain_bin_width, area_scaling(),
                                grain_size_start_mu, grain_size_start_sigma, grain_size_max_x, grain_size_stdabw_low, grain_size_stdabw_high,
                                grain_size_y_max, fit_mu, fit_sigma, filepath_grain_size_fit.c_str());

        for (int j=0; j<grain_quantile_values.size(); j++)
        {
            float i = (float)grain_size_histogram.size() - log10(grain_quantile_values[j]*area_scaling())*grain_bin_width +0.5f;
            float grains = exp(-0.5f * (log(i)-fit_mu)*(log(i)-fit_mu) / (fit_sigma*fit_sigma) ) / (i*fit_sigma*sqrt(2.0f*PI))*
                (nr_areas-found_bubble_areas.size()-found_border_areas.size());
            float size_diff = pow(10.0f,(log10(grain_quantile_values[j]*area_scaling())*grain_bin_width+0.5f)/grain_bin_width-log10(area_scaling()))-
                pow(10.0f,(log10(grain_quantile_values[j]*area_scaling())*grain_bin_width-0.5f)/grain_bin_width-log10(area_scaling()));
            grain_quantile_errors.push_back(size_diff*grain_quantile_values[j]/grains);
        }

        plot.draw_values_errors("Higher percent quantile positions", "Grain size quantiles [mm^2]", "Grain size quantiles distribution",
                                grain_quantile_values, grain_quantile_errors, filepath_quantile_grains.c_str());
        plot.draw_histogram("Roundness factor (4pi*Area/Perimeter^2)","Relative occurrence","Grain roundness distribution", grain_roundness_histogram,
                            grain_roundness_bin_width, 1.0f, grain_roundness_mean, grain_roundness_standard_deviation, grain_roundness_y_max,
                            filepath_grain_roundness.c_str());
        plot.draw_histogram_log(length_unit,"Relative occurrence","Grain equivalent radius distribution", grain_equiv_radius_histogram,
                                grain_equiv_radius_bin_width, length_scaling(), grain_equiv_radius_mean, grain_equiv_radius_standard_deviation,
                                grain_equiv_radius_y_max, filepath_grain_equiv_radius.c_str());
        plot.draw_histogram("Normalized grain equivalent radius","Relative occurrence","Normalized grain equivalent radius distribution",
                            grain_equiv_radius_norm_histogram, grain_equiv_radius_norm_bin_width, 1.0f, grain_equiv_radius_norm_y_max,
                            filepath_grain_equiv_radius_norm.c_str());
        plot.draw_histogram("Box flattening factor","Relative occurrence","Box flattening distribution", grain_box_flattening_histogram,
                            grain_box_flattening_bin_width, 1.0f, grain_box_flattening_mean, grain_box_flattening_standard_deviation,
                            grain_box_flattening_y_max, filepath_grain_box_flattening.c_str());
        plot.draw_histogram_log(length_unit,"Relative occurrence","Grain box width", grain_box_width_histogram, grain_box_width_bin_width, length_scaling(),
                                grain_box_width_mean, grain_box_width_standard_deviation, grain_box_width_y_max, filepath_grain_box_width.c_str());
        plot.draw_histogram_log(length_unit,"Relative occurrence","Grain box height", grain_box_height_histogram, grain_box_height_bin_width, length_scaling(),
                                grain_box_height_mean, grain_box_height_standard_deviation, grain_box_height_y_max, filepath_grain_box_height.c_str());
        plot.draw_histogram_log(length_unit,"Relative occurrence","Ellipse long axis", grain_ellipse_long_axis_histogram,
                                grain_ellipse_long_axis_bin_width, length_scaling(), grain_ellipse_long_axis_mean,
                                grain_ellipse_long_axis_standard_deviation, grain_ellipse_long_axis_y_max, filepath_grain_ellipse_long_axis.c_str());
        plot.draw_histogram("Ellipse flattening factor","Relative occurrence","Ellipse flattening distribution", grain_ellipse_flattening_histogram,
                            grain_ellipse_flattening_bin_width, 1.0f, grain_ellipse_flattening_mean, grain_ellipse_flattening_standard_deviation,
                            grain_ellipse_flattening_y_max, filepath_grain_ellipse_flattening.c_str());
        plot.draw_histogram_rose("Ellipse long axis angle","Relative occurrence","Ellipse long axis angle distribution",
                            grain_ellipse_long_axis_angle_histogram, grain_ellipse_long_axis_angle_bin_width, angle_scaling,
                            grain_ellipse_long_axis_angle_mean, grain_ellipse_long_axis_angle_standard_deviation, grain_ellipse_long_axis_angle_y_max,
                            filepath_grain_ellipse_long_axis_angle.c_str());
        plot.draw_histogram_log(length_unit,"Relative occurrence","Grain width", grain_area_width_histogram, grain_area_width_bin_width, length_scaling(),
                                grain_area_width_mean, grain_area_width_standard_deviation, grain_area_width_y_max, filepath_grain_width.c_str());
        plot.draw_histogram_log(length_unit,"Relative occurrence","Grain height", grain_area_height_histogram, grain_area_height_bin_width, length_scaling(),
                                grain_area_height_mean, grain_area_height_standard_deviation, grain_area_height_y_max, filepath_grain_height.c_str());
        plot.draw_histogram("Grain vertical flattening factor","Relative occurrence","Grain vertical flattening distribution", grain_area_flattening_histogram,
                            grain_area_flattening_bin_width, 1.0f, grain_area_flattening_mean, grain_area_flattening_standard_deviation,
                            grain_area_flattening_y_max, filepath_grain_flattening.c_str());
        plot.draw_histogram("Perimeter ratio","Relative occurrence","Grain perimeter ratio distribution", grain_perimeter_ratio_histogram,
                            grain_perimeter_ratio_bin_width, 1.0f, grain_perimeter_ratio_mean, grain_perimeter_ratio_standard_deviation,
                            grain_perimeter_ratio_y_max, filepath_grainperimeter_ratio.c_str());
        plot.draw_histogram_log(length_unit,"Relative occurrence","Distance to next bubble", min_bubble_distance_histogram,
                                min_bubble_distance_bin_width, length_scaling(), min_bubble_distance_mean, min_bubble_distance_standard_deviation,
                                min_bubble_distance_y_max, filepath_min_bubble_distance.c_str());
        plot.draw_values_errors_log_log(area_unit, length_unit, "Distance to next bubble vs. mean grain size",relaxation_mean,
                                        relaxation_standard_deviation, grain_bin_width, area_scaling(), min_bubble_distance_bin_width, length_scaling(),
                                        filepath_relaxation.c_str());
        plot.draw_values_errors_log_log(length_unit, area_unit, "Mean grain size vs. distance to next bubble",relaxation2_mean,
                                        relaxation2_standard_deviation, min_bubble_distance_bin_width, length_scaling(), grain_bin_width, area_scaling(),
                                        filepath_relaxation2.c_str());
        if(grain_arc_number_histogram.size()>2)
            plot.draw_histogram_lognorm("Number of grain boundaries","Relative occurrence","Number of grain boundaries distribution",
                            grain_arc_number_histogram, grain_arc_number_bin_width, 1.0f, grain_arc_number_start_mu, grain_arc_number_start_sigma,
                            grain_arc_number_max_x, grain_arc_number_standard_deviation, grain_arc_number_y_max,filepath_grain_arc_number.c_str());
        else
        {
            grain_arc_number_max_x=0.0f;
            grain_arc_number_standard_deviation=0.0f;
        }
        if(grain_neighbors_histogram.size()>2)
            plot.draw_histogram_lognorm("Number of neighbors","Relative occurrence","Number of neighbors distribution", grain_neighbors_histogram,
                            grain_neighbors_bin_width, 1.0f, grain_neighbors_start_mu, grain_neighbors_start_sigma, grain_neighbors_max_x,
                            grain_neighbors_standard_deviation, grain_neighbors_y_max, filepath_grain_neighbors.c_str());
        else
        {
            grain_neighbors_max_x=0.0f;
            grain_neighbors_standard_deviation=0.0f;
        }
        plot.draw_histogram_log(length_unit,"Relative occurrence","Grain boundary length distribution", grain_arc_length_histogram,
                                grain_arc_length_bin_width, length_scaling(), grain_arc_length_mean, grain_arc_length_standard_deviation,
                                grain_arc_length_y_max, filepath_grain_arc_length.c_str());
        plot.draw_histogram_log(length_unit,"Relative occurrence","Grain longest boundary length distribution", grain_longest_arc_length_histogram,
                                grain_longest_arc_length_bin_width, length_scaling(), grain_longest_arc_length_mean,
                                grain_longest_arc_length_standard_deviation, grain_longest_arc_length_y_max, filepath_grain_longest_arc_length.c_str());
        plot.draw_histogram_log_log("Dislocation density difference in m^(-2)", dislocation_title,
                                    "Estimation of dislocation density difference at grain boundaries", dislocation_density_histogram,
                                    dislocation_density_bin_width, dislocation_density_scaling, log10(curvature_min), dislocation_density_mean,
                                    dislocation_density_standard_deviation, dislocation_density_y_max, filepath_dislocation_density.c_str());
        plot.draw_histogram_rose("Boundary orientation","Relative occurrence","Boundary orientation distribution",
                            boundary_orientation_histogram, boundary_orientation_bin_width, angle_scaling,
                            boundary_orientation_mean, boundary_orientation_standard_deviation, boundary_orientation_y_max,
                            filepath_boundary_orientation.c_str());
        plot.draw_histogram("Dihedral angle [degree]","Relative occurrence", "Dihedral angle distribution", dihedral_angle_histogram,
                            dihedral_angle_bin_width,angle_scaling,dihedral_angle_mean,dihedral_angle_standard_deviation, dihedral_angle_y_max,
                            filepath_dihedral_angle.c_str());
        plot.draw_histogram("Dihedral angle [degree]","Relative occurrence", "Dihedral angle (2nd version) distribution", dihedral_angle2_histogram,
                            dihedral_angle2_bin_width,angle_scaling,dihedral_angle2_mean,dihedral_angle2_standard_deviation, dihedral_angle2_y_max,
                            filepath_dihedral_angle2.c_str());
        plot.draw_histogram_numbers("Nr of turning points","Relative occurrence", "Grain boundary turning points distribution", turning_point_histogram,
                            1.0f, 1.0f, turning_point_mean, turning_point_standard_deviation, 1.0f, filepath_turning_point.c_str());

        sprintf(step, "%d", boundary_step());

        x_axis="Longest grain boundaries * ";
        x_axis.append(step);

        plot.draw_values_errors(x_axis, "Dislocation density difference in m^(-2)", "Dislocation density difference of longest grain boundaries",
                                disl_dens_boundary_region_values, disl_dens_boundary_region_errors, filepath_disl_dens_longest_boundaries.c_str());
        plot.draw_values_errors("* 10 percent longest grain boundaries", "Dislocation density difference in m^(-2)",
                                "Dislocation density difference of longest grain boundaries", disl_dens_boundary_percent_values,
                                disl_dens_boundary_percent_errors, filepath_disl_dens_percent_boundaries.c_str());

        plot.draw_depth_errors("Percentage of area filled by grains", "Mean grain size [mm^2]",
                               "Mean grain size vs. percentage of area filled by grains distribution", percentage_values, percentage_filled_grains_mean,
                               percentage_filled_grains_standard_deviation, percentage_filled_grains_standard_deviation, 0.0f,
                               filepath_percentage_filled_grains.c_str());

        //*******************************
        //SAVE NEW MERGED FINAL STRUCTURE
        //*******************************

        if (minimal_grain_size>0 || minimal_bubble_distance>0)
        {
            std::string filepath_new_classification2=filepath_new_classification;
            filepath_new_classification2.resize(filepath_new_classification2.size()-3);

            std::stringstream str;
            if(minimal_bubble_distance>0) str << "." << minimal_bubble_distance << "." << minimal_grain_size;
            else str << "." << minimal_grain_size;
            filepath_new_classification2.append(str.str());

            GrainBoundNet.paramFileGBN = &paramFile;
            GrainBoundNet.save_final_structure(filepath_new_classification2);
        }

        //*************
        //VISUALISATION
        //*************
        /*
        //visualisation of merging result
        vigra::FImage grain_image(dim_x,dim_y);
     
        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
            {
                grain_image(x,y)=3;
            }
        }
        
        //bubble areas to output
        for (int a=0; a<found_bubble_areas.size(); a++)
        {
            for(size_t k=0; k<areas[found_bubble_areas[a]-1].size(); ++k)
            {
                int x=areas[found_bubble_areas[a]-1][k].x;
                int y=areas[found_bubble_areas[a]-1][k].y;
                grain_image(x,y)=0;
            }
        }
        /*        
        //DUBUG close bubble areas to output, NOTE these areas will be erase by merging
        for(int b=0; b<close_bubble_areas.size(); b++)
        {
            int area=close_bubble_areas[b];

            for(size_t k=0; k<areas[area].size(); ++k)
            {
                int x=areas[area][k].x;
                int y=areas[area][k].y;
                grain_image(x,y)=1;
            }

            if(old_grain_arc_index[area].size()>0)
            {
                for (int a=0;a<(int)old_grain_arc_index[area].size();a++)//grain boundaries
                {
                    //now we loop over the points in this arc
                    for(int p=0;p<(int)arcs[old_grain_arc_index[area][a]].size();p++)
                    {
                        int x=arcs[old_grain_arc_index[area][a]][p].x;
                        int y=arcs[old_grain_arc_index[area][a]][p].y;
                        grain_image(x,y)=0;
                    }
                }
            }
        }
        * /
        //grain areas and boundaries to output
        for(size_t area=0;area<nr_areas;area++)
        {
            if(old_grain_arc_index[area].size()>0)
            {
                for(size_t k=0; k<areas[area].size(); ++k)
                {
                    int x=areas[area][k].x;
                    int y=areas[area][k].y;
                    grain_image(x,y)=1;
                }

                for (int a=0;a<(int)old_grain_arc_index[area].size();a++)//grain boundaries
                {
                    //now we loop over the points in this arc
                    for(int p=0;p<(int)arcs[old_grain_arc_index[area][a]].size();p++)
                    {
                        int x=arcs[old_grain_arc_index[area][a]][p].x;
                        int y=arcs[old_grain_arc_index[area][a]][p].y;
                        grain_image(x,y)=0;
                    }
                }
            }
        }
        /*
        //grain junctions to output image
        for(int j=0; j<grain_junctions.size(); j++)
        {
            int xx=junctions[grain_junctions[j]].x;
            int yy=junctions[grain_junctions[j]].y;
            grain_image(xx-1,yy-1)=3;
            grain_image(xx,yy-1)=3;
            grain_image(xx+1,yy-1)=3;
            grain_image(xx-1,yy)=3;
            grain_image(xx,yy)=3;
            grain_image(xx+1,yy)=3;
            grain_image(xx-1,yy+1)=3;
            grain_image(xx,yy+1)=3;
            grain_image(xx+1,yy+1)=3;
        }
        */
        /*
        //grain bubble junction to output image
        for(int j=0; j<grain_bubble_junctions.size(); j++)
        {
            int xx=junctions[grain_bubble_junctions[j]].x;
            int yy=junctions[grain_bubble_junctions[j]].y;
            grain_image(xx-1,yy-1)=3;
            grain_image(xx,yy-1)=3;
            grain_image(xx+1,yy-1)=3;
            grain_image(xx-1,yy)=3;
            grain_image(xx,yy)=3;
            grain_image(xx+1,yy)=3;
            grain_image(xx-1,yy+1)=3;
            grain_image(xx,yy+1)=3;
            grain_image(xx+1,yy+1)=3;
        }
        * /

        //grain center of mass positions to output image        
        for (int area=0; area<nr_areas; area++)
        {
            if (grain_area_size[area]>0)
            {
                int xx=grain_area_center_mass[area].x;
                int yy=grain_area_center_mass[area].y;
                grain_image(xx-1,yy-1)=3;
                grain_image(xx,yy-1)=3;
                grain_image(xx+1,yy-1)=3;
                grain_image(xx-1,yy)=3;
                grain_image(xx,yy)=3;
                grain_image(xx+1,yy)=3;
                grain_image(xx-1,yy+1)=3;
                grain_image(xx,yy+1)=3;
                grain_image(xx+1,yy+1)=3;
            }
        }

        /*
        //bubble center of mass positions to output image
        for (int bubble=0; bubble<bubble_area_center_mass.size(); bubble++)
        {
            int xx=bubble_area_center_mass[bubble].x;
            int yy=bubble_area_center_mass[bubble].y;
            grain_image(xx-1,yy-1)=3;
            grain_image(xx,yy-1)=3;
            grain_image(xx+1,yy-1)=3;
            grain_image(xx-1,yy)=3;
            grain_image(xx,yy)=3;
            grain_image(xx+1,yy)=3;
            grain_image(xx-1,yy+1)=3;
            grain_image(xx,yy+1)=3;
            grain_image(xx+1,yy+1)=3;
        }
        * /
        //export visualisation
        std::string filepath_grains=path_rf_predictions;
        filepath_grains.append("grains/");
        filepath_grains.append(param_file_name.c_str());

        filepath_grains.append(filename);

        std::stringstream sts;
        if(minimal_bubble_distance>0) sts << "." << minimal_bubble_distance << "." << minimal_grain_size <<".jpg";
        else sts << "." << minimal_grain_size <<".jpg";

        filepath_grains.append(sts.str());
        exportImage(srcImageRange(grain_image), vigra::ImageExportInfo(filepath_grains.c_str()));
        */
        //********************************
        //OUTPUT VALUES FOR DEPTH PROFILES
        //********************************

        std::string filepath_grainsize_all=path_results;
        filepath_grainsize_all.append("grainsize_all");
        std::string filepath_grainsize_fit=path_results;
        filepath_grainsize_fit.append("grainsize_fit");
        std::string filepath_grainsize_step=path_results;
        filepath_grainsize_step.append("grainsize_step");
        std::string filepath_grainsize_percent=path_results;
        filepath_grainsize_percent.append("grainsize_percent");
        std::string filepath_grainsize_quantile=path_results;
        filepath_grainsize_quantile.append("grainsize_quantile");
        std::string filepath_grainsize_bin=path_results;
        filepath_grainsize_bin.append("grainsize_bin");
        std::string filepath_grainradius_bin=path_results;
        filepath_grainradius_bin.append("grainradius_bin");
        std::string filepath_grainradius_norm_bin=path_results;
        filepath_grainradius_norm_bin.append("grainradius_norm_bin");
        std::string filepath_grain_neighbors_bin=path_results;
        filepath_grain_neighbors_bin.append("grain_neighbors_bin");
        std::string filepath_grainshape=path_results;
        filepath_grainshape.append("grainshape");
        std::string filepath_grain_equivradius=path_results;
        filepath_grain_equivradius.append("grain_equivradius");
        std::string filepath_grain_boxshape=path_results;
        filepath_grain_boxshape.append("grain_boxshape");
        std::string filepath_grain_boxwidth=path_results;
        filepath_grain_boxwidth.append("grain_boxwidth");
        std::string filepath_grain_boxheight=path_results;
        filepath_grain_boxheight.append("grain_boxheight");
        std::string filepath_grain_ellipselong=path_results;
        filepath_grain_ellipselong.append("grain_ellipselong");
        std::string filepath_grain_ellipseshape=path_results;
        filepath_grain_ellipseshape.append("grain_ellipseshape");
        std::string filepath_grain_ellipseangle=path_results;
        filepath_grain_ellipseangle.append("grain_ellipseangle");
        std::string filepath_grainwidth=path_results;
        filepath_grainwidth.append("grainwidth");
        std::string filepath_grainheight=path_results;
        filepath_grainheight.append("grainheight");
        std::string filepath_grainflattening=path_results;
        filepath_grainflattening.append("grainflattening");
        std::string filepath_nr_grainarcs=path_results;
        filepath_nr_grainarcs.append("nr_grainarcs");
        std::string filepath_nr_neighbors=path_results;
        filepath_nr_neighbors.append("nr_neighbors");
        std::string filepath_length_grainarcs=path_results;
        filepath_length_grainarcs.append("length_grainarcs");
        std::string filepath_longest_grainarcs=path_results;
        filepath_longest_grainarcs.append("longest_grainarcs");
        std::string filepath_dihedral_angles=path_results;
        filepath_dihedral_angles.append("dihedral_angles");
        std::string filepath_dihedral_angles2=path_results;
        filepath_dihedral_angles2.append("dihedral_angles2");
        std::string filepath_dislocation_densities=path_results;
        filepath_dislocation_densities.append("dislocation_densities");
        std::string filepath_boundaryorientation=path_results;
        filepath_boundaryorientation.append("boundary_orientation");
        std::string filepath_turning_points=path_results;
        filepath_turning_points.append("turning_points");
        std::string filepath_correlation_ellipse_box_flattening=path_results;
        filepath_correlation_ellipse_box_flattening.append("corr_ellipse_box_flattening");
        std::string filepath_nr_values=path_results;
        filepath_nr_values.append("nr_values");
        std::string filepath_disl_dens_step=path_results;
        filepath_disl_dens_step.append("disl_dens_step");
        std::string filepath_disl_dens_percent=path_results;
        filepath_disl_dens_percent.append("disl_dens_percent");
        std::string filepath_percentage_filled=path_results;
        filepath_percentage_filled.append("percentage_filled_grains");
        std::string filepath_nr_percentage_filled=path_results;
        filepath_nr_percentage_filled.append("nr_percentage_filled_grains");
        std::string filepath_grain_perimeter_ratio = path_results;
        filepath_grain_perimeter_ratio.append("grain_perimeter_ratio");

        std::string suffix=s.str();
        suffix.resize(suffix.size()-3);
        suffix.append("txt");

        filepath_grainsize_all.append(suffix.c_str());
        filepath_grainsize_fit.append(suffix.c_str());
        filepath_grainsize_step.append(suffix.c_str());
        filepath_grainsize_percent.append(suffix.c_str());
        filepath_grainsize_quantile.append(suffix.c_str());
        filepath_grainsize_bin.append(suffix.c_str());
        filepath_grainradius_bin.append(suffix.c_str());
        filepath_grainradius_norm_bin.append(suffix.c_str());
        filepath_grain_neighbors_bin.append(suffix.c_str());
        filepath_grainshape.append(suffix.c_str());
        filepath_grain_equivradius.append(suffix.c_str());
        filepath_grain_boxshape.append(suffix.c_str());
        filepath_grain_boxwidth.append(suffix.c_str());
        filepath_grain_boxheight.append(suffix.c_str());
        filepath_grain_ellipselong.append(suffix.c_str());
        filepath_grain_ellipseshape.append(suffix.c_str());
        filepath_grain_ellipseangle.append(suffix.c_str());
        filepath_grainwidth.append(suffix.c_str());
        filepath_grainheight.append(suffix.c_str());
        filepath_grainflattening.append(suffix.c_str());
        filepath_nr_grainarcs.append(suffix.c_str());
        filepath_nr_neighbors.append(suffix.c_str());
        filepath_length_grainarcs.append(suffix.c_str());
        filepath_longest_grainarcs.append(suffix.c_str());
        filepath_dihedral_angles.append(suffix.c_str());
        filepath_dihedral_angles2.append(suffix.c_str());
        filepath_dislocation_densities.append(suffix.c_str());
        filepath_boundaryorientation.append(suffix.c_str());
        filepath_turning_points.append(suffix.c_str());
        filepath_correlation_ellipse_box_flattening.append(suffix.c_str());
        filepath_nr_values.append(suffix.c_str());
        filepath_disl_dens_step.append(suffix.c_str());
        filepath_disl_dens_percent.append(suffix.c_str());
        filepath_percentage_filled.append(suffix.c_str());
        filepath_nr_percentage_filled.append(suffix.c_str());
        filepath_grain_perimeter_ratio.append(suffix.c_str());

        std::ofstream grainsize_all_file(filepath_grainsize_all.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grainsize_fit_file(filepath_grainsize_fit.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grainsize_step_file(filepath_grainsize_step.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grainsize_percent_file(filepath_grainsize_percent.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grainsize_quantile_file(filepath_grainsize_quantile.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grainsize_bin_file(filepath_grainsize_bin.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grainradius_bin_file(filepath_grainradius_bin.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grainradius_norm_bin_file(filepath_grainradius_norm_bin.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grain_neighbors_bin_file(filepath_grain_neighbors_bin.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grainshape_file(filepath_grainshape.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grain_equivradius_file(filepath_grain_equivradius.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grain_boxshape_file(filepath_grain_boxshape.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grain_boxwidth_file(filepath_grain_boxwidth.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grain_boxheight_file(filepath_grain_boxheight.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grain_ellipselong_file(filepath_grain_ellipselong.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grain_ellipseshape_file(filepath_grain_ellipseshape.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grain_ellipseangle_file(filepath_grain_ellipseangle.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grainwidth_file(filepath_grainwidth.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grainheight_file(filepath_grainheight.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grainflattening_file(filepath_grainflattening.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream nr_grainarcs_file(filepath_nr_grainarcs.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream nr_neighbors_file(filepath_nr_neighbors.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream length_grainarcs_file(filepath_length_grainarcs.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream longest_grainarcs_file(filepath_longest_grainarcs.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream dihedral_angles_file(filepath_dihedral_angles.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream dihedral_angles2_file(filepath_dihedral_angles2.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream dislocation_densities_file(filepath_dislocation_densities.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream boundary_orientation_file(filepath_boundaryorientation.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream turning_points_file(filepath_turning_points.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream correlation_ellipse_box_flattening_file(filepath_correlation_ellipse_box_flattening.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream nr_values_file(filepath_nr_values.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream disl_dens_step_file(filepath_disl_dens_step.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream disl_dens_percent_file(filepath_disl_dens_percent.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream percentage_filled_grains_file(filepath_percentage_filled.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream nr_percentage_filled_grains_file(filepath_nr_percentage_filled.c_str(), std::ios_base::out | std::ios_base::app);
        std::ofstream grain_perimeter_ratio_file(filepath_grain_perimeter_ratio.c_str(), std::ios_base::out | std::ios_base::app);

        grainsize_all_file <<filename<<" "<<grain_size_mean*area_scaling()<<" "<<grain_size_standard_deviation*area_scaling()<< "\n";
        grainsize_fit_file <<filename<<" "<<grain_size_max_x*area_scaling()<<" "<<grain_size_stdabw_low*area_scaling()<<" "<<
                                            grain_size_stdabw_high*area_scaling()<< "\n";

        grainsize_step_file <<filename<<" ";
        for(int step=0; step<grain_region_values.size(); step++) grainsize_step_file
            <<grain_region_values[step]*area_scaling()<<" "<<grain_region_errors[step]*area_scaling()<<" ";
        grainsize_step_file << "\n";

        grainsize_percent_file <<filename<<" ";
        for (int p=0; p<grain_percent_values.size(); p++) grainsize_percent_file
            <<grain_percent_values[p]*area_scaling()<<" "<<grain_percent_errors[p]*area_scaling()<<" ";
        grainsize_percent_file<< "\n";

        grainsize_quantile_file <<filename<<" ";
        for (int q=0; q<grain_quantile_values.size(); q++) grainsize_quantile_file
            <<grain_quantile_values[q]*area_scaling()<<" "<<grain_quantile_errors[q]*area_scaling()<<" ";
        grainsize_quantile_file << "\n";

        grainsize_bin_file <<filename<<" ";
        for(int bin=0; bin<grain_size_histogram.size(); bin++) grainsize_bin_file
            <<grain_size_histogram[bin]<<" ";
        grainsize_bin_file << "\n";

        grainradius_bin_file <<filename<<" ";
        for(int bin=0; bin<grain_equiv_radius_histogram.size(); bin++) grainradius_bin_file
            <<grain_equiv_radius_histogram[bin]<<" ";
        grainradius_bin_file << "\n";

        grainradius_norm_bin_file <<filename<<" ";
        for(int bin=0; bin<grain_equiv_radius_norm_histogram.size(); bin++) grainradius_norm_bin_file
            <<grain_equiv_radius_norm_histogram[bin]<<" ";
        grainradius_norm_bin_file << "\n";

        grain_neighbors_bin_file <<filename<<" ";
        for(int bin=0; bin<grain_neighbors_histogram.size(); bin++) grain_neighbors_bin_file
            <<grain_neighbors_histogram[bin]<<" ";
        grain_neighbors_bin_file << "\n";

        grain_perimeter_ratio_file << filename << " " << grain_perimeter_ratio_mean << " " << grain_perimeter_ratio_standard_deviation << "\n";        

        grainshape_file <<filename<<" "<<grain_roundness_mean<<" "<<grain_roundness_standard_deviation<< "\n";
        grain_equivradius_file <<filename<<" "<<grain_equiv_radius_mean*length_scaling()<<" "<<grain_equiv_radius_standard_deviation*length_scaling()<< "\n";
        grain_boxshape_file <<filename<<" "<<grain_box_flattening_mean<<" "<<grain_box_flattening_standard_deviation<< "\n";
        grain_boxwidth_file <<filename<<" "<<grain_box_width_mean*length_scaling()<<" "<<grain_box_width_standard_deviation*length_scaling()<< "\n";
        grain_boxheight_file <<filename<<" "<<grain_box_height_mean*length_scaling()<<" "<<grain_box_height_standard_deviation*length_scaling()<< "\n";
        grain_ellipselong_file <<filename<<" "<<grain_ellipse_long_axis_mean*length_scaling()<<" "<<grain_ellipse_long_axis_standard_deviation*length_scaling()<<
            "\n";
        grain_ellipseshape_file <<filename<<" "<<grain_ellipse_flattening_mean<<" "<<grain_ellipse_flattening_standard_deviation<< "\n";
        grain_ellipseangle_file <<filename<<" "<<grain_ellipse_long_axis_angle_mean<<" "<<grain_ellipse_long_axis_angle_standard_deviation<< "\n";
        grainwidth_file <<filename<<" "<<grain_area_width_mean*length_scaling()<<" "<<grain_area_width_standard_deviation*length_scaling()<< "\n";
        grainheight_file <<filename<<" "<<grain_area_height_mean*length_scaling()<<" "<<grain_area_height_standard_deviation*length_scaling()<< "\n";
        grainflattening_file <<filename<<" "<<grain_area_flattening_mean<<" "<<grain_area_flattening_standard_deviation<< "\n";
        nr_grainarcs_file <<filename<<" "<<grain_arc_number_max_x<<" "<<grain_arc_number_standard_deviation<< "\n";
        nr_neighbors_file <<filename<<" "<<grain_neighbors_max_x<<" "<<grain_neighbors_standard_deviation<< "\n";
        length_grainarcs_file <<filename<<" "<<grain_arc_length_mean*length_scaling()<<" "<<grain_arc_length_standard_deviation*length_scaling()<< "\n";
        longest_grainarcs_file <<filename<<" "<<grain_longest_arc_length_mean*length_scaling()<<" "<<
            grain_longest_arc_length_standard_deviation*length_scaling()<< "\n";
        dihedral_angles_file <<filename<<" "<<dihedral_angle_mean<<" "<<dihedral_angle_standard_deviation<< "\n";
        dihedral_angles2_file <<filename<<" "<<dihedral_angle2_mean<<" "<<dihedral_angle2_standard_deviation<< "\n";
        dislocation_densities_file <<filename<<" "<<dislocation_density_mean<<" "<<dislocation_density_standard_deviation<< "\n";
        boundary_orientation_file <<filename<<" "<<boundary_orientation_mean<<" "<<boundary_orientation_standard_deviation<< "\n";
        turning_points_file <<filename<<" "<<turning_point_mean<<" "<<turning_point_standard_deviation<< "\n";
        correlation_ellipse_box_flattening_file <<filename<<" "<<correlation_ellipse_box_flattening<< "\n";
        nr_values_file<< filename<<" "<<grain_size_values.size()<<" "<<grain_box_flattening_values2.size()<<" "<<grain_ellipse_flattening_values2.size()<<
            " "<<grain_ellipse_long_axis_values2.size()<<" "<<grain_arc_length_values.size()<<" "<<grain_longest_arc_length_values2.size()<<" "<<
            dihedral_angle_values.size()<<" "<<dislocation_density_values.size()<<" "<<turning_point_values.size()<<" "<<grain_junctions.size()<< "\n";

        disl_dens_step_file <<filename<<" ";
        for(int step=0; step<disl_dens_boundary_region_values.size(); step++) disl_dens_step_file
            <<disl_dens_boundary_region_values[step]<<" "<<disl_dens_boundary_region_errors[step]<<" ";
        disl_dens_step_file << "\n";

        disl_dens_percent_file <<filename<<" ";
        for (int p=0; p<disl_dens_boundary_percent_values.size(); p++) disl_dens_percent_file
            <<disl_dens_boundary_percent_values[p]<<" "<<disl_dens_boundary_percent_errors[p]<<" ";
        disl_dens_percent_file<< "\n";

        percentage_filled_grains_file <<filename<<" ";
        for (int p=0; p<percentage_filled_grains_mean.size(); p++) percentage_filled_grains_file
            <<percentage_filled_grains_mean[p]<<" "<<percentage_filled_grains_standard_deviation[p]<<" ";
        percentage_filled_grains_file<< "\n";

        nr_percentage_filled_grains_file <<filename<<" ";
        for (int p=0; p<nr_percentage_filled_grains.size(); p++) nr_percentage_filled_grains_file
            <<nr_percentage_filled_grains[p]<<" ";
        nr_percentage_filled_grains_file<< "\n";

        grainsize_all_file.close();
        grainsize_fit_file.close();
        grainsize_step_file.close();
        grainsize_percent_file.close();
        grainsize_quantile_file.close();
        grainsize_bin_file.close();
        grainradius_bin_file.close();
        grainradius_norm_bin_file.close();
        grain_neighbors_bin_file.close();
        grainshape_file.close();
        grain_equivradius_file.close();
        grain_boxshape_file.close();
        grain_boxwidth_file.close();
        grain_boxheight_file.close();
        grain_ellipselong_file.close();
        grain_ellipseshape_file.close();
        grain_ellipseangle_file.close();
        grainwidth_file.close();
        grainheight_file.close();
        grainflattening_file.close();
        nr_grainarcs_file.close();
        nr_neighbors_file.close();
        length_grainarcs_file.close();
        longest_grainarcs_file.close();
        dihedral_angles_file.close();
        dihedral_angles2_file.close();
        dislocation_densities_file.close();
        boundary_orientation_file.close();
        turning_points_file.close();
        correlation_ellipse_box_flattening_file.close();
        nr_values_file.close();
        disl_dens_step_file.close();
        disl_dens_percent_file.close();
        percentage_filled_grains_file.close();
        nr_percentage_filled_grains_file.close();
        grain_perimeter_ratio_file.close();

        //***********************************
        //EXTRACT ALL PARAMETERS FOR ANALYSIS
        //***********************************

        std::string filepath_parameters=path_results;
        filepath_parameters.append(filename);

        suffix.resize(suffix.size()-4);
        filepath_parameters.append(suffix.c_str());
        
        para.save_extracted_parameters(filepath_parameters);
    }

    //*********************
    //FILL OTHER HISTOGRAMS
    //*********************

    //fill bubble size histogram
    for (int area=0; area<nr_areas; area++)
    {
        if (bubble_area_size[area]>bubble_size_max)
        {
            bubble_size_histogram.resize((int)bubble_bin_width*log10(bubble_area_size[area])+1);
            bubble_size_max=bubble_area_size[area];
        }
        if (bubble_area_size[area]>bubble_size_min) 
        {
            bubble_size_histogram[(int)bubble_bin_width*log10(bubble_area_size[area])]++;
            bubble_size_values.push_back(bubble_area_size[area]/area_scaling());
        }
    }
    calculate_mean_standard_deviation(bubble_size_values, bubble_size_mean, bubble_size_standard_deviation);

    for (int area=0; area<nr_areas; area++)
    {
        if (bubble_area_size[area]>0) 
        {
            //find all neighbors
            std::list<int> grain_neighbors;
            std::list<int> bubble_neighbors;
            bool border_contact=false;

            for(int a=0; a<bubble_arc_index[area].size() && !border_contact; a++)
            {
                int arc_index=bubble_arc_index[area][a];

                if(twoCellNeighbors(arc_index,0)==area+1 && grain_area_size[twoCellNeighbors(arc_index,1)-1]>10000)
                    grain_neighbors.push_back(twoCellNeighbors(arc_index,1));
                else if(twoCellNeighbors(arc_index,1)==area+1 && grain_area_size[twoCellNeighbors(arc_index,0)-1]>10000)
                    grain_neighbors.push_back(twoCellNeighbors(arc_index,0));

                if(twoCellNeighbors(arc_index,0)==area+1 && bubble_area_size[twoCellNeighbors(arc_index,1)-1]>0)
                    bubble_neighbors.push_back(twoCellNeighbors(arc_index,1));
                else if(twoCellNeighbors(arc_index,1)==area+1 && bubble_area_size[twoCellNeighbors(arc_index,0)-1]>0)
                    bubble_neighbors.push_back(twoCellNeighbors(arc_index,0));

                if ((twoCellNeighbors(arc_index,0)==area+1 && grain_area_size[twoCellNeighbors(arc_index,1)-1]==0 &&
                     bubble_area_size[twoCellNeighbors(arc_index,1)-1]==0) || (twoCellNeighbors(arc_index,1)==area+1 &&
                     grain_area_size[twoCellNeighbors(arc_index,0)-1]==0 && bubble_area_size[twoCellNeighbors(arc_index,0)-1]==0))
                    border_contact=true;
            }

            if (border_contact || bubble_neighbors.size()>0)//grains with contact to border or other bubbles can not handled correctly
            {
                grain_neighbors.clear();
                bubble_neighbors.clear();
            }

            grain_neighbors.sort();
            grain_neighbors.unique();

            int nr_grain_neighbors=(int)grain_neighbors.size();

            //nr of bubble grain neighbors
            if (nr_grain_neighbors>bubble_grain_neighbors_max)
            {
                bubble_grain_neighbors_histogram.resize((int)(nr_grain_neighbors/bubble_grain_neighbors_bin_width)+1);
                bubble_grain_neighbors_max=nr_grain_neighbors;
            }
            if (nr_grain_neighbors>grain_neighbors_min)
            {
                bubble_grain_neighbors_histogram[(int)(nr_grain_neighbors/bubble_grain_neighbors_bin_width)]++;
                bubble_grain_neighbors_values.push_back(nr_grain_neighbors);
            }
        }
    }
    calculate_mean_standard_deviation(bubble_grain_neighbors_values, bubble_grain_neighbors_mean, bubble_grain_neighbors_standard_deviation);

    //fill grain center of mass histogram
    for (int area=0; area<nr_areas; area++)
    {
        if (grain_area_center_mass[area].y>grain_center_of_mass_max)
        {
            grain_center_of_mass_histogram.resize((int)(grain_area_center_mass[area].y/grain_center_of_mass_bin_width)+1);
            grain_center_of_mass_max=grain_area_center_mass[area].y;
        }
        if (grain_area_center_mass[area].y>grain_center_of_mass_min)
        {
            grain_center_of_mass_histogram[(int)(grain_area_center_mass[area].y/grain_center_of_mass_bin_width)]++;
            grain_center_of_mass_values.push_back(grain_area_center_mass[area].y);
        }
    }
    calculate_mean_standard_deviation(grain_center_of_mass_values, grain_center_of_mass_mean, grain_center_of_mass_standard_deviation);

    std::string filepath_bubble_size=filepath_plots;
    std::string filepath_bubble_grain_neighbors=filepath_plots;
    std::string filepath_grain_center_of_mass=filepath_plots;

    filepath_bubble_size.append(".bubble_size.svg");
    filepath_bubble_grain_neighbors.append(".bubble_grain_neighbors.svg");
    filepath_grain_center_of_mass.append(".grain_center_of_mass.svg");

    plot.draw_histogram_log(area_unit,"Relative occurrence","Bubble size distribution", bubble_size_histogram, bubble_bin_width, area_scaling(),
                            bubble_size_mean, bubble_size_standard_deviation, bubble_size_y_max, filepath_bubble_size.c_str());
    plot.draw_histogram("Number of neighboring grains","Relative occurrence","Number of grain neighbors distribution",
                        bubble_grain_neighbors_histogram, bubble_grain_neighbors_bin_width, 1.0f, bubble_grain_neighbors_mean,
                        bubble_grain_neighbors_standard_deviation, bubble_grain_neighbors_y_max, filepath_bubble_grain_neighbors.c_str());
    plot.draw_histogram("Y-coordinate","Relative occurrence","Grain center of mass distribution", grain_center_of_mass_histogram,
                        grain_center_of_mass_bin_width, 1.0f, grain_center_of_mass_mean, grain_center_of_mass_standard_deviation, grain_center_of_mass_y_max,
                        filepath_grain_center_of_mass.c_str());

    delete bubble_area_size;
    delete grain_area_size;
}

void cloudy_bands(std::string file_a1, std::string file_a2, std::string file_b1, std::string file_b2, std::string file_c1, std::string file_c2,
                  std::string file_d1, std::string file_d2, std::string file_e1, std::string file_e2, std::string file_f1, std::string file_f2,
                  std::string path_rf_predictions, std::string path_plots)
{
    int dim_x=8192;

    size_t nr_areas;
    long * bubble_area_size;
    
    std::string filepath_new_classification_a1=path_rf_predictions;
    std::string filepath_new_classification_a2=path_rf_predictions;
    std::string filepath_new_classification_b1=path_rf_predictions;
    std::string filepath_new_classification_b2=path_rf_predictions;
    std::string filepath_new_classification_c1=path_rf_predictions;
    std::string filepath_new_classification_c2=path_rf_predictions;
    std::string filepath_new_classification_d1=path_rf_predictions;
    std::string filepath_new_classification_d2=path_rf_predictions;
    std::string filepath_new_classification_e1=path_rf_predictions;
    std::string filepath_new_classification_e2=path_rf_predictions;
    std::string filepath_new_classification_f1=path_rf_predictions;
    std::string filepath_new_classification_f2=path_rf_predictions;

    filepath_new_classification_a1.append(file_a1);
    filepath_new_classification_a2.append(file_a2);
    filepath_new_classification_b1.append(file_b1);
    filepath_new_classification_b2.append(file_b2);
    filepath_new_classification_c1.append(file_c1);
    filepath_new_classification_c2.append(file_c2);
    filepath_new_classification_d1.append(file_d1);
    filepath_new_classification_d2.append(file_d2);
    filepath_new_classification_e1.append(file_e1);
    filepath_new_classification_e2.append(file_e2);
    filepath_new_classification_f1.append(file_f1);
    filepath_new_classification_f2.append(file_f2);

    filepath_new_classification_a1.append(".h5");
    filepath_new_classification_a2.append(".h5");
    filepath_new_classification_b1.append(".h5");
    filepath_new_classification_b2.append(".h5");
    filepath_new_classification_c1.append(".h5");
    filepath_new_classification_c2.append(".h5");
    filepath_new_classification_d1.append(".h5");
    filepath_new_classification_d2.append(".h5");
    filepath_new_classification_e1.append(".h5");
    filepath_new_classification_e2.append(".h5");
    filepath_new_classification_f1.append(".h5");
    filepath_new_classification_f2.append(".h5");

    //IMPORT RESULTS FROM HDF5 file
    gbn grainBoundNet_a1;
    grainBoundNet_a1.load_final_structure(filepath_new_classification_a1);    
    delete bubble_area_size;

    gbn grainBoundNet_a2;
    grainBoundNet_a2.load_final_structure(filepath_new_classification_a2);
    delete bubble_area_size;

    gbn grainBoundNet_b1;
    grainBoundNet_b1.load_final_structure(filepath_new_classification_b1);
    delete bubble_area_size;

    gbn grainBoundNet_b2;
    grainBoundNet_b2.load_final_structure(filepath_new_classification_b1);
    delete bubble_area_size;

    gbn grainBoundNet_c1;
    grainBoundNet_c1.load_final_structure(filepath_new_classification_b1);
    delete bubble_area_size;

    gbn grainBoundNet_c2;
    grainBoundNet_c2.load_final_structure(filepath_new_classification_b1);
    delete bubble_area_size;

    gbn grainBoundNet_d1;
    grainBoundNet_d1.load_final_structure(filepath_new_classification_b1);
    delete bubble_area_size;

    gbn grainBoundNet_d2;
    grainBoundNet_d2.load_final_structure(filepath_new_classification_b1);
    delete bubble_area_size;

    gbn grainBoundNet_e1;
    grainBoundNet_e1.load_final_structure(filepath_new_classification_b1);
    delete bubble_area_size;

    gbn grainBoundNet_e2;
    grainBoundNet_e2.load_final_structure(filepath_new_classification_b1);
    delete bubble_area_size;

    gbn grainBoundNet_f1;
    grainBoundNet_f1.load_final_structure(filepath_new_classification_b1);
    delete bubble_area_size;

    gbn grainBoundNet_f2;
    grainBoundNet_f2.load_final_structure(filepath_new_classification_b1);
    delete bubble_area_size;

    //Set up references to the attributes
    long * & grain_area_size_a1 = grainBoundNet_a1.grain_area_size;
    long * & grain_area_size_a2 = grainBoundNet_a2.grain_area_size;
    long * & grain_area_size_b1 = grainBoundNet_b1.grain_area_size;
    long * & grain_area_size_b2 = grainBoundNet_b2.grain_area_size;
    long * & grain_area_size_c1 = grainBoundNet_c1.grain_area_size;
    long * & grain_area_size_c2 = grainBoundNet_c2.grain_area_size;
    long * & grain_area_size_d1 = grainBoundNet_d1.grain_area_size;
    long * & grain_area_size_d2 = grainBoundNet_d2.grain_area_size;
    long * & grain_area_size_e1 = grainBoundNet_e1.grain_area_size;
    long * & grain_area_size_e2 = grainBoundNet_e2.grain_area_size;
    long * & grain_area_size_f1 = grainBoundNet_f1.grain_area_size;
    long * & grain_area_size_f2 = grainBoundNet_f2.grain_area_size;    

    std::vector<point> & grain_area_center_mass_a1 = grainBoundNet_a1.grain_area_center_mass;
    std::vector<point> & grain_area_center_mass_a2 = grainBoundNet_a2.grain_area_center_mass;
    std::vector<point> & grain_area_center_mass_b1 = grainBoundNet_b1.grain_area_center_mass;
    std::vector<point> & grain_area_center_mass_b2 = grainBoundNet_b2.grain_area_center_mass;
    std::vector<point> & grain_area_center_mass_c1 = grainBoundNet_c1.grain_area_center_mass;
    std::vector<point> & grain_area_center_mass_c2 = grainBoundNet_c2.grain_area_center_mass;
    std::vector<point> & grain_area_center_mass_d1 = grainBoundNet_d1.grain_area_center_mass;
    std::vector<point> & grain_area_center_mass_d2 = grainBoundNet_d2.grain_area_center_mass;
    std::vector<point> & grain_area_center_mass_e1 = grainBoundNet_e1.grain_area_center_mass;
    std::vector<point> & grain_area_center_mass_e2 = grainBoundNet_e2.grain_area_center_mass;
    std::vector<point> & grain_area_center_mass_f1 = grainBoundNet_f1.grain_area_center_mass;
    std::vector<point> & grain_area_center_mass_f2 = grainBoundNet_f2.grain_area_center_mass;

    //initialise plplot class
    plplot plot = plplot();

    //grain center of mass histogram
    int grain_center_of_mass_min=0;
    int grain_center_of_mass_max=0;
    float grain_center_of_mass_bin_width=700.0f;
    int nr_x_regions=2;
    std::vector< std::vector<int> > grain_center_of_mass_histogram(nr_x_regions);
    float grain_center_of_mass_y_max=0.0f;

    std::vector<grain> right_side;
    std::vector<grain> left_side;
    int right, x_diff, nr_grains;

    for (int area=0; area<grain_area_center_mass_a1.size(); area++)
    {
        grain one_entry;
        one_entry.x=grain_area_center_mass_a1[area].x;
        one_entry.y=grain_area_center_mass_a1[area].y;
        one_entry.size=grain_area_size_a1[area];
        if (one_entry.size>0) right_side.push_back(one_entry);
    }

    for (int area=0; area<grain_area_center_mass_a2.size(); area++)
    {
        grain one_entry;
        one_entry.x=grain_area_center_mass_a2[area].x;
        one_entry.y=grain_area_center_mass_a2[area].y;
        one_entry.size=grain_area_size_a2[area];
        if (one_entry.size>0) left_side.push_back(one_entry);
    }

    std::sort(right_side.begin(),right_side.end(),grain_cmp);
    std::sort(left_side.begin(),left_side.end(),grain_cmp);

    right=0;
    x_diff=0;
    nr_grains=0;

    for(int left=0; left<left_side.size(); left++)
    {
        while (right_side[right].size-left_side[left].size>9 && right<right_side.size()) right++;
        if(fabs(left_side[left].size-right_side[right].size)<10 && fabs(left_side[left].y-right_side[right].y)<100 &&
            left_side[left].x-right_side[right].x>5700 && left_side[left].x-right_side[right].x<6300)
            {
                x_diff+=left_side[left].x-right_side[right].x;
                nr_grains++;
            }
        if (right>0) right--;
    }
    std::cout<<"found x diff: "<<dim_x-x_diff/nr_grains<<std::endl;

    for (int area=0; area<grain_area_center_mass_a1.size(); area++)
    {
        if (grain_area_center_mass_a1[area].y>grain_center_of_mass_max)
        {
            grain_center_of_mass_histogram[0].resize((int)(grain_area_center_mass_a1[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_histogram[1].resize((int)(grain_area_center_mass_a1[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_max=grain_area_center_mass_a1[area].y;
        }
        if (grain_area_center_mass_a1[area].x>grain_center_of_mass_min && grain_area_center_mass_a1[area].x<0.5*(dim_x+x_diff/nr_grains))
        {
            grain_center_of_mass_histogram[0][(int)(grain_area_center_mass_a1[area].y/grain_center_of_mass_bin_width)]++;
            grain_center_of_mass_histogram[1][(int)(grain_area_center_mass_a1[area].y/grain_center_of_mass_bin_width)]++;
        }
    }

    for (int area=0; area<grain_area_center_mass_a2.size(); area++)
    {
        if (grain_area_center_mass_a2[area].y>grain_center_of_mass_max)
        {
            grain_center_of_mass_histogram[0].resize((int)(grain_area_center_mass_a2[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_histogram[1].resize((int)(grain_area_center_mass_a2[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_max=grain_area_center_mass_a2[area].y;
        }
        if (grain_area_center_mass_a2[area].x>0.5*(dim_x-x_diff/nr_grains))
        {
            grain_center_of_mass_histogram[1][(int)(grain_area_center_mass_a2[area].y/grain_center_of_mass_bin_width)]++;
        }
    }

    right_side.clear();
    left_side.clear();

    for (int area=0; area<grain_area_center_mass_b1.size(); area++)
    {
        grain one_entry;
        one_entry.x=grain_area_center_mass_b1[area].x;
        one_entry.y=grain_area_center_mass_b1[area].y;
        one_entry.size=grain_area_size_b1[area];
        if (one_entry.size>0) right_side.push_back(one_entry);
    }

    for (int area=0; area<grain_area_center_mass_b2.size(); area++)
    {
        grain one_entry;
        one_entry.x=grain_area_center_mass_b2[area].x;
        one_entry.y=grain_area_center_mass_b2[area].y;
        one_entry.size=grain_area_size_b2[area];
        if (one_entry.size>0) left_side.push_back(one_entry);
    }

    std::sort(right_side.begin(),right_side.end(),grain_cmp);
    std::sort(left_side.begin(),left_side.end(),grain_cmp);

    right=0;
    x_diff=0;
    nr_grains=0;

    for(int left=0; left<left_side.size(); left++)
    {
        while (right_side[right].size-left_side[left].size>9 && right<right_side.size()) right++;
        if(fabs(left_side[left].size-right_side[right].size)<10 && fabs(left_side[left].y-right_side[right].y)<100 &&
            left_side[left].x-right_side[right].x>5700 && left_side[left].x-right_side[right].x<6300)
            {
                x_diff+=left_side[left].x-right_side[right].x;
                nr_grains++;
            }
        if (right>0) right--;
    }
    std::cout<<"found x diff: "<<dim_x-x_diff/nr_grains<<std::endl;

    for (int area=0; area<grain_area_center_mass_b1.size(); area++) grain_area_center_mass_b1[area].y+=grain_center_of_mass_max
                                                                                                     +grain_center_of_mass_bin_width;
    for (int area=0; area<grain_area_center_mass_b2.size(); area++) grain_area_center_mass_b2[area].y+=grain_center_of_mass_max
                                                                                                     +grain_center_of_mass_bin_width;

    for (int area=0; area<grain_area_center_mass_b1.size(); area++)
    {
        if (grain_area_center_mass_b1[area].y>grain_center_of_mass_max)
        {
            grain_center_of_mass_histogram[0].resize((int)(grain_area_center_mass_b1[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_histogram[1].resize((int)(grain_area_center_mass_b1[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_max=grain_area_center_mass_b1[area].y;
        }
        if (grain_area_center_mass_b1[area].x>grain_center_of_mass_min && grain_area_center_mass_b1[area].x<0.5*(dim_x+x_diff/nr_grains))
        {
            grain_center_of_mass_histogram[0][(int)(grain_area_center_mass_b1[area].y/grain_center_of_mass_bin_width)]++;
            grain_center_of_mass_histogram[1][(int)(grain_area_center_mass_b1[area].y/grain_center_of_mass_bin_width)]++;
        }
    }

    for (int area=0; area<grain_area_center_mass_b2.size(); area++)
    {
        if (grain_area_center_mass_b2[area].y>grain_center_of_mass_max)
        {
            grain_center_of_mass_histogram[0].resize((int)(grain_area_center_mass_b2[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_histogram[1].resize((int)(grain_area_center_mass_b2[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_max=grain_area_center_mass_b2[area].y;
        }
        if (grain_area_center_mass_b2[area].x>0.5*(dim_x-x_diff/nr_grains))
        {
            grain_center_of_mass_histogram[1][(int)(grain_area_center_mass_b2[area].y/grain_center_of_mass_bin_width)]++;
        }
    }

    right_side.clear();
    left_side.clear();

    for (int area=0; area<grain_area_center_mass_c1.size(); area++)
    {
        grain one_entry;
        one_entry.x=grain_area_center_mass_c1[area].x;
        one_entry.y=grain_area_center_mass_c1[area].y;
        one_entry.size=grain_area_size_c1[area];
        if (one_entry.size>0) right_side.push_back(one_entry);
    }

    for (int area=0; area<grain_area_center_mass_c2.size(); area++)
    {
        grain one_entry;
        one_entry.x=grain_area_center_mass_c2[area].x;
        one_entry.y=grain_area_center_mass_c2[area].y;
        one_entry.size=grain_area_size_c2[area];
        if (one_entry.size>0) left_side.push_back(one_entry);
    }

    std::sort(right_side.begin(),right_side.end(),grain_cmp);
    std::sort(left_side.begin(),left_side.end(),grain_cmp);

    right=0;
    x_diff=0;
    nr_grains=0;

    for(int left=0; left<left_side.size(); left++)
    {
        while (right_side[right].size-left_side[left].size>9 && right<right_side.size()) right++;
        if(fabs(left_side[left].size-right_side[right].size)<10 && fabs(left_side[left].y-right_side[right].y)<100 &&
            left_side[left].x-right_side[right].x>5700 && left_side[left].x-right_side[right].x<6300)
            {
                x_diff+=left_side[left].x-right_side[right].x;
                nr_grains++;
            }
        if (right>0) right--;
    }
    std::cout<<"found x diff: "<<dim_x-x_diff/nr_grains<<std::endl;

    for (int area=0; area<grain_area_center_mass_c1.size(); area++) grain_area_center_mass_c1[area].y+=grain_center_of_mass_max
                                                                                                     +grain_center_of_mass_bin_width;
    for (int area=0; area<grain_area_center_mass_c2.size(); area++) grain_area_center_mass_c2[area].y+=grain_center_of_mass_max
                                                                                                     +grain_center_of_mass_bin_width;

    for (int area=0; area<grain_area_center_mass_c1.size(); area++)
    {
        if (grain_area_center_mass_c1[area].y>grain_center_of_mass_max)
        {
            grain_center_of_mass_histogram[0].resize((int)(grain_area_center_mass_c1[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_histogram[1].resize((int)(grain_area_center_mass_c1[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_max=grain_area_center_mass_c1[area].y;
        }
        if (grain_area_center_mass_c1[area].x>grain_center_of_mass_min && grain_area_center_mass_c1[area].x<0.5*(dim_x+x_diff/nr_grains))
        {
            grain_center_of_mass_histogram[0][(int)(grain_area_center_mass_c1[area].y/grain_center_of_mass_bin_width)]++;
            grain_center_of_mass_histogram[1][(int)(grain_area_center_mass_c1[area].y/grain_center_of_mass_bin_width)]++;
        }
    }

    for (int area=0; area<grain_area_center_mass_c2.size(); area++)
    {
        if (grain_area_center_mass_c2[area].y>grain_center_of_mass_max)
        {
            grain_center_of_mass_histogram[0].resize((int)(grain_area_center_mass_c2[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_histogram[1].resize((int)(grain_area_center_mass_c2[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_max=grain_area_center_mass_c2[area].y;
        }
        if (grain_area_center_mass_c2[area].x>0.5*(dim_x-x_diff/nr_grains))
        {
            grain_center_of_mass_histogram[1][(int)(grain_area_center_mass_c2[area].y/grain_center_of_mass_bin_width)]++;
        }
    }

    right_side.clear();
    left_side.clear();

    for (int area=0; area<grain_area_center_mass_d1.size(); area++)
    {
        grain one_entry;
        one_entry.x=grain_area_center_mass_d1[area].x;
        one_entry.y=grain_area_center_mass_d1[area].y;
        one_entry.size=grain_area_size_d1[area];
        if (one_entry.size>0) right_side.push_back(one_entry);
    }

    for (int area=0; area<grain_area_center_mass_d2.size(); area++)
    {
        grain one_entry;
        one_entry.x=grain_area_center_mass_d2[area].x;
        one_entry.y=grain_area_center_mass_d2[area].y;
        one_entry.size=grain_area_size_d2[area];
        if (one_entry.size>0) left_side.push_back(one_entry);
    }

    std::sort(right_side.begin(),right_side.end(),grain_cmp);
    std::sort(left_side.begin(),left_side.end(),grain_cmp);

    right=0;
    x_diff=0;
    nr_grains=0;

    for(int left=0; left<left_side.size(); left++)
    {
        while (right_side[right].size-left_side[left].size>9 && right<right_side.size()) right++;
        if(fabs(left_side[left].size-right_side[right].size)<10 && fabs(left_side[left].y-right_side[right].y)<100 &&
            left_side[left].x-right_side[right].x>5700 && left_side[left].x-right_side[right].x<6300)
            {
                x_diff+=left_side[left].x-right_side[right].x;
                nr_grains++;
            }
        if (right>0) right--;
    }
    std::cout<<"found x diff: "<<dim_x-x_diff/nr_grains<<std::endl;

    for (int area=0; area<grain_area_center_mass_d1.size(); area++) grain_area_center_mass_d1[area].y+=grain_center_of_mass_max
                                                                                                     +grain_center_of_mass_bin_width;
    for (int area=0; area<grain_area_center_mass_d2.size(); area++) grain_area_center_mass_d2[area].y+=grain_center_of_mass_max
                                                                                                     +grain_center_of_mass_bin_width;

    for (int area=0; area<grain_area_center_mass_d1.size(); area++)
    {
        if (grain_area_center_mass_d1[area].y>grain_center_of_mass_max)
        {
            grain_center_of_mass_histogram[0].resize((int)(grain_area_center_mass_d1[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_histogram[1].resize((int)(grain_area_center_mass_d1[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_max=grain_area_center_mass_d1[area].y;
        }
        if (grain_area_center_mass_d1[area].x>grain_center_of_mass_min && grain_area_center_mass_d1[area].x<0.5*(dim_x+x_diff/nr_grains))
        {
            grain_center_of_mass_histogram[0][(int)(grain_area_center_mass_d1[area].y/grain_center_of_mass_bin_width)]++;
            grain_center_of_mass_histogram[1][(int)(grain_area_center_mass_d1[area].y/grain_center_of_mass_bin_width)]++;
        }
    }

    for (int area=0; area<grain_area_center_mass_d2.size(); area++)
    {
        if (grain_area_center_mass_d2[area].y>grain_center_of_mass_max)
        {
            grain_center_of_mass_histogram[0].resize((int)(grain_area_center_mass_d2[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_histogram[1].resize((int)(grain_area_center_mass_d2[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_max=grain_area_center_mass_d2[area].y;
        }
        if (grain_area_center_mass_d2[area].x>0.5*(dim_x-x_diff/nr_grains))
        {
            grain_center_of_mass_histogram[1][(int)(grain_area_center_mass_d2[area].y/grain_center_of_mass_bin_width)]++;
        }
    }

    right_side.clear();
    left_side.clear();

    for (int area=0; area<grain_area_center_mass_e1.size(); area++)
    {
        grain one_entry;
        one_entry.x=grain_area_center_mass_e1[area].x;
        one_entry.y=grain_area_center_mass_e1[area].y;
        one_entry.size=grain_area_size_e1[area];
        if (one_entry.size>0) right_side.push_back(one_entry);
    }

    for (int area=0; area<grain_area_center_mass_e2.size(); area++)
    {
        grain one_entry;
        one_entry.x=grain_area_center_mass_e2[area].x;
        one_entry.y=grain_area_center_mass_e2[area].y;
        one_entry.size=grain_area_size_e2[area];
        if (one_entry.size>0) left_side.push_back(one_entry);
    }

    std::sort(right_side.begin(),right_side.end(),grain_cmp);
    std::sort(left_side.begin(),left_side.end(),grain_cmp);

    right=0;
    x_diff=0;
    nr_grains=0;

    for(int left=0; left<left_side.size(); left++)
    {
        while (right_side[right].size-left_side[left].size>9 && right<right_side.size()) right++;
        if(fabs(left_side[left].size-right_side[right].size)<10 && fabs(left_side[left].y-right_side[right].y)<100 &&
            left_side[left].x-right_side[right].x>5700 && left_side[left].x-right_side[right].x<6300)
            {
                x_diff+=left_side[left].x-right_side[right].x;
                nr_grains++;
            }
        if (right>0) right--;
    }
    std::cout<<"found x diff: "<<dim_x-x_diff/nr_grains<<std::endl;

    for (int area=0; area<grain_area_center_mass_e1.size(); area++) grain_area_center_mass_e1[area].y+=grain_center_of_mass_max
                                                                                                     +grain_center_of_mass_bin_width;
    for (int area=0; area<grain_area_center_mass_e2.size(); area++) grain_area_center_mass_e2[area].y+=grain_center_of_mass_max
                                                                                                     +grain_center_of_mass_bin_width;

    for (int area=0; area<grain_area_center_mass_e1.size(); area++)
    {
        if (grain_area_center_mass_e1[area].y>grain_center_of_mass_max)
        {
            grain_center_of_mass_histogram[0].resize((int)(grain_area_center_mass_e1[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_histogram[1].resize((int)(grain_area_center_mass_e1[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_max=grain_area_center_mass_e1[area].y;
        }
        if (grain_area_center_mass_e1[area].x>grain_center_of_mass_min && grain_area_center_mass_e1[area].x<0.5*(dim_x+x_diff/nr_grains))
        {
            grain_center_of_mass_histogram[0][(int)(grain_area_center_mass_e1[area].y/grain_center_of_mass_bin_width)]++;
            grain_center_of_mass_histogram[1][(int)(grain_area_center_mass_e1[area].y/grain_center_of_mass_bin_width)]++;
        }
    }

    for (int area=0; area<grain_area_center_mass_e2.size(); area++)
    {
        if (grain_area_center_mass_e2[area].y>grain_center_of_mass_max)
        {
            grain_center_of_mass_histogram[0].resize((int)(grain_area_center_mass_e2[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_histogram[1].resize((int)(grain_area_center_mass_e2[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_max=grain_area_center_mass_e2[area].y;
        }
        if (grain_area_center_mass_e2[area].x>0.5*(dim_x-x_diff/nr_grains))
        {
            grain_center_of_mass_histogram[1][(int)(grain_area_center_mass_e2[area].y/grain_center_of_mass_bin_width)]++;
        }
    }

    right_side.clear();
    left_side.clear();

    for (int area=0; area<grain_area_center_mass_f1.size(); area++)
    {
        grain one_entry;
        one_entry.x=grain_area_center_mass_f1[area].x;
        one_entry.y=grain_area_center_mass_f1[area].y;
        one_entry.size=grain_area_size_f1[area];
        if (one_entry.size>0) right_side.push_back(one_entry);
    }

    for (int area=0; area<grain_area_center_mass_f2.size(); area++)
    {
        grain one_entry;
        one_entry.x=grain_area_center_mass_f2[area].x;
        one_entry.y=grain_area_center_mass_f2[area].y;
        one_entry.size=grain_area_size_f2[area];
        if (one_entry.size>0) left_side.push_back(one_entry);
    }

    std::sort(right_side.begin(),right_side.end(),grain_cmp);
    std::sort(left_side.begin(),left_side.end(),grain_cmp);

    right=0;
    x_diff=0;
    nr_grains=0;

    for(int left=0; left<left_side.size(); left++)
    {
        while (right_side[right].size-left_side[left].size>9 && right<right_side.size()) right++;
        if(fabs(left_side[left].size-right_side[right].size)<10 && fabs(left_side[left].y-right_side[right].y)<100 &&
            left_side[left].x-right_side[right].x>5700 && left_side[left].x-right_side[right].x<6300)
            {
                x_diff+=left_side[left].x-right_side[right].x;
                nr_grains++;
            }
        if (right>0) right--;
    }
    std::cout<<"found x diff: "<<dim_x-x_diff/nr_grains<<std::endl;

    for (int area=0; area<grain_area_center_mass_f1.size(); area++) grain_area_center_mass_f1[area].y+=grain_center_of_mass_max
                                                                                                     +grain_center_of_mass_bin_width;
    for (int area=0; area<grain_area_center_mass_f2.size(); area++) grain_area_center_mass_f2[area].y+=grain_center_of_mass_max
                                                                                                     +grain_center_of_mass_bin_width;

    for (int area=0; area<grain_area_center_mass_f1.size(); area++)
    {
        if (grain_area_center_mass_f1[area].y>grain_center_of_mass_max)
        {
            grain_center_of_mass_histogram[0].resize((int)(grain_area_center_mass_f1[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_histogram[1].resize((int)(grain_area_center_mass_f1[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_max=grain_area_center_mass_f1[area].y;
        }
        if (grain_area_center_mass_f1[area].x>grain_center_of_mass_min && grain_area_center_mass_f1[area].x<0.5*(dim_x+x_diff/nr_grains))
        {
            grain_center_of_mass_histogram[0][(int)(grain_area_center_mass_f1[area].y/grain_center_of_mass_bin_width)]++;
            grain_center_of_mass_histogram[1][(int)(grain_area_center_mass_f1[area].y/grain_center_of_mass_bin_width)]++;
        }
    }

    for (int area=0; area<grain_area_center_mass_f2.size(); area++)
    {
        if (grain_area_center_mass_f2[area].y>grain_center_of_mass_max)
        {
            grain_center_of_mass_histogram[0].resize((int)(grain_area_center_mass_f2[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_histogram[1].resize((int)(grain_area_center_mass_f2[area].y/grain_center_of_mass_bin_width)+1,0);
            grain_center_of_mass_max=grain_area_center_mass_f2[area].y;
        }
        if (grain_area_center_mass_f2[area].x>0.5*(dim_x-x_diff/nr_grains))
        {
            grain_center_of_mass_histogram[1][(int)(grain_area_center_mass_f2[area].y/grain_center_of_mass_bin_width)]++;
        }
    }

    //scaling from pixels to length
    float length_scaling=2072.3;//one mm is 193.5 Pixels

    std::string filepath_grain_center_of_mass=path_plots;
    filepath_grain_center_of_mass.append("grain_center_of_mass.svg");

    plot.draw_histogram("Depth in cm (NEEM ice core at 1867 m)","Relative occurrence in %","Grain center of mass points profile",
        grain_center_of_mass_histogram, grain_center_of_mass_bin_width, length_scaling, grain_center_of_mass_y_max, filepath_grain_center_of_mass.c_str());

    delete grain_area_size_a1;
    delete grain_area_size_a2;
    delete grain_area_size_b1;
    delete grain_area_size_b2;
    delete grain_area_size_c1;
    delete grain_area_size_c2;
    delete grain_area_size_d1;
    delete grain_area_size_d2;
    delete grain_area_size_e1;
    delete grain_area_size_e2;
    delete grain_area_size_f1;
    delete grain_area_size_f2;
}

void double_grain_size(std::string path_results, std::string suffix1, std::string suffix2, int minimal_bubble_distance, int minimal_grain_size, std::string path_plots)
{
    std::stringstream s;
    if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size <<".txt";
    else s << "." << minimal_grain_size <<".txt";

    std::string filepath_to_results1=path_results;
    std::string filepath_to_results2=path_results;
    std::string filepath_to_results3=path_results;
    std::string filepath_to_results4=path_results;
    std::string filepath_to_results5=path_results;
    std::string filepath_to_results6=path_results;

    filepath_to_results1.append(suffix1.c_str());
    filepath_to_results2.append(suffix2.c_str());
    filepath_to_results3.append(suffix1.c_str());
    filepath_to_results4.append(suffix2.c_str());
    filepath_to_results5.append(suffix1.c_str());
    filepath_to_results6.append(suffix2.c_str());

    filepath_to_results1.append("grainsize_bin");
    filepath_to_results2.append("grainsize_bin");
    filepath_to_results3.append("grainradius_bin");
    filepath_to_results4.append("grainradius_bin");
    filepath_to_results5.append("grainradius_norm_bin");
    filepath_to_results6.append("grainradius_norm_bin");

    filepath_to_results1.append(s.str());
    filepath_to_results2.append(s.str());
    filepath_to_results3.append(s.str());
    filepath_to_results4.append(s.str());
    filepath_to_results5.append(s.str());
    filepath_to_results6.append(s.str());

    //initialise plplot class
    plplot plot = plplot();

    //scaling from pixels to length
    float length_scaling=193.5;//one mm is 193.5 Pixels
    float area_scaling=37444.0;//one mm is 193.5 Pixels

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

    //grain size histogram
    int grain_size_max=0;
    float grain_bin_width=5.0f;
    std::vector< std::vector<int> > grain_size_histogram;
    grain_size_histogram.resize(2);
    float grain_size_y_max=0.0f;

    //grain equiv radius histogram
    int grain_equiv_radius_max=0;
    float grain_equiv_radius_bin_width=10.0f;
    std::vector< std::vector<int> > grain_equiv_radius_histogram;
    grain_equiv_radius_histogram.resize(2);
    float grain_equiv_radius_y_max=0.0f;

    //grain equiv radius norm histogram
    float grain_equiv_radius_norm_max=0.0f;
    float grain_equiv_radius_norm_bin_width=0.2f;
    std::vector< std::vector<int> > grain_equiv_radius_norm_histogram;
    grain_equiv_radius_norm_histogram.resize(2);
    float grain_equiv_radius_norm_y_max=0.0f;

    std::string teststring;

    //first file
    std::ifstream grainsize_bin_file1(filepath_to_results1.c_str());
    std::ifstream temp_grainsize_bin_file1(filepath_to_results1.c_str());
    temp_grainsize_bin_file1>>teststring;

    if(grainsize_bin_file1 && teststring.size()!=0)
    {
        //read histogram values from first line
        std::string name;
        grainsize_bin_file1>>name;                    

        char line[256];
        grainsize_bin_file1.getline(line,256);

        int old_pos=1;
        for(int c=1; c<grainsize_bin_file1.gcount(); c++)
        {
            if(line[c]==32)
            {
                char number[10]="";
                memmove(number+0,line+old_pos,c-old_pos);
                old_pos=c+1;
                grain_size_histogram[0].push_back(atof(number));
            }
        }

        grainsize_bin_file1.close();
        temp_grainsize_bin_file1.close();
    }

    //second file
    std::ifstream grainsize_bin_file2(filepath_to_results2.c_str());
    std::ifstream temp_grainsize_bin_file2(filepath_to_results2.c_str());
    temp_grainsize_bin_file2>>teststring;

    if(grainsize_bin_file2 && teststring.size()!=0)
    {
        //read histogram values from first line
        std::string name;
        grainsize_bin_file2>>name;                    

        char line[256];
        grainsize_bin_file2.getline(line,256);

        int old_pos=1;
        for(int c=1; c<grainsize_bin_file2.gcount(); c++)
        {
            if(line[c]==32)
            {
                char number[10]="";
                memmove(number+0,line+old_pos,c-old_pos);
                old_pos=c+1;
                grain_size_histogram[1].push_back(atof(number));
            }
        }

        grainsize_bin_file2.close();
        temp_grainsize_bin_file2.close();
    }

    teststring="";

    //first file
    std::ifstream grainradius_bin_file1(filepath_to_results3.c_str());
    std::ifstream temp_grainradius_bin_file1(filepath_to_results3.c_str());
    temp_grainradius_bin_file1>>teststring;

    if(grainradius_bin_file1 && teststring.size()!=0)
    {
        //read histogram values from first line
        std::string name;
        grainradius_bin_file1>>name;                    

        char line[256];
        grainradius_bin_file1.getline(line,256);

        int old_pos=1;
        for(int c=1; c<grainradius_bin_file1.gcount(); c++)
        {
            if(line[c]==32)
            {
                char number[10]="";
                memmove(number+0,line+old_pos,c-old_pos);
                old_pos=c+1;
                grain_equiv_radius_histogram[0].push_back(atof(number));
            }
        }

        grainradius_bin_file1.close();
        temp_grainradius_bin_file1.close();
    }

    //second file
    std::ifstream grainradius_bin_file2(filepath_to_results4.c_str());
    std::ifstream temp_grainradius_bin_file2(filepath_to_results4.c_str());
    temp_grainradius_bin_file2>>teststring;

    if(grainradius_bin_file2 && teststring.size()!=0)
    {
        //read histogram values from first line
        std::string name;
        grainradius_bin_file2>>name;                    

        char line[256];
        grainradius_bin_file2.getline(line,256);

        int old_pos=1;
        for(int c=1; c<grainradius_bin_file2.gcount(); c++)
        {
            if(line[c]==32)
            {
                char number[10]="";
                memmove(number+0,line+old_pos,c-old_pos);
                old_pos=c+1;
                grain_equiv_radius_histogram[1].push_back(atof(number));
            }
        }

        grainradius_bin_file2.close();
        temp_grainradius_bin_file2.close();
    }

    teststring="";

    //first file
    std::ifstream grainradius_norm_bin_file1(filepath_to_results5.c_str());
    std::ifstream temp_grainradius_norm_bin_file1(filepath_to_results5.c_str());
    temp_grainradius_norm_bin_file1>>teststring;

    if(grainradius_norm_bin_file1 && teststring.size()!=0)
    {
        //read histogram values from first line
        std::string name;
        grainradius_norm_bin_file1>>name;                    

        char line[256];
        grainradius_norm_bin_file1.getline(line,256);

        int old_pos=1;
        for(int c=1; c<grainradius_norm_bin_file1.gcount(); c++)
        {
            if(line[c]==32)
            {
                char number[10]="";
                memmove(number+0,line+old_pos,c-old_pos);
                old_pos=c+1;
                grain_equiv_radius_norm_histogram[0].push_back(atof(number));
            }
        }

        grainradius_norm_bin_file1.close();
        temp_grainradius_norm_bin_file1.close();
    }

    //second file
    std::ifstream grainradius_norm_bin_file2(filepath_to_results6.c_str());
    std::ifstream temp_grainradius_norm_bin_file2(filepath_to_results6.c_str());
    temp_grainradius_norm_bin_file2>>teststring;

    if(grainradius_norm_bin_file2 && teststring.size()!=0)
    {
        //read histogram values from first line
        std::string name;
        grainradius_norm_bin_file2>>name;                    

        char line[256];
        grainradius_norm_bin_file2.getline(line,256);

        int old_pos=1;
        for(int c=1; c<grainradius_norm_bin_file2.gcount(); c++)
        {
            if(line[c]==32)
            {
                char number[10]="";
                memmove(number+0,line+old_pos,c-old_pos);
                old_pos=c+1;
                grain_equiv_radius_norm_histogram[1].push_back(atof(number));
            }
        }

        grainradius_norm_bin_file2.close();
        temp_grainradius_norm_bin_file2.close();
    }

    //set filepath
    std::string filepath_grain_size=path_plots;
    std::string filepath_grain_radius=path_plots;
    std::string filepath_grain_radius_norm=path_plots;
    filepath_grain_size.append("combined_grain_size.svg");
    filepath_grain_radius.append("combined_grain_radius.svg");
    filepath_grain_radius_norm.append("combined_grain_radius_norm.svg");

    plot.draw_multi_histogram_log(area_unit,"Relative occurrence","", grain_size_histogram, grain_bin_width, area_scaling,
                                  grain_size_y_max, filepath_grain_size.c_str());
    plot.draw_multi_histogram_log("Equivalent radius in mm","Relative occurrence","", grain_equiv_radius_histogram, grain_equiv_radius_bin_width,
                                  length_scaling, grain_equiv_radius_y_max, filepath_grain_radius.c_str());
    plot.draw_multi_histogram("Normalized grain equivalent radius","Relative occurrence","", grain_equiv_radius_norm_histogram,
                              grain_equiv_radius_norm_bin_width, 1.0f, grain_equiv_radius_norm_y_max, filepath_grain_radius_norm.c_str());
}
