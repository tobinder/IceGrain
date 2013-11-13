/*! \file subgrain.cpp
 * \brief Used for subgrain_extraction.
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
#include "remove_subgrain_arc.cpp"

int get_nr_true(std::vector<bool> bool_vector)
{
    int nr=0;

    for (int i=0; i<bool_vector.size(); i++)
    {
        if(bool_vector[i]==true) nr++;
    }

    return nr;
}

//quantile
float get_quantile(std::vector<float> values,float q)
{

    int size=values.size();
    if(size!=1)
    {
        //sort the "values"-vector
        sort( values.begin(), values.end() );

        //f_index is the index of the lower quartil
        float f_index=((float)size-1)*q;
        //if we make f_index-floor(f_index) we know how to round
        //the quantile
        float round_value=f_index-floor(f_index);

        if(round_value==0)
        {
            return values[(int)f_index];
        }
        else
        {
            float return_value=0;
            int   low_index=floor(f_index);
            return_value=((1-round_value)*values[low_index]+ (round_value)*values[low_index+1]);
            return return_value;
        }
    }
    else
    {
        return values[1];
    }
}

void find_subgrains(std::string filepath_to_feature_file,std::string path_to_image,std::string path_to_ws_image,
                    std::string filepath_rf_predictions, std::string path_subGB_rf_predictions, std::string path_parameters,
                    int minimal_grain_size, int minimal_bubble_distance, std::string param_file, ParameterFile paramFile)
{
    /*vigra::BasicImage<unsigned int> ws_region_image;
    marray::Marray<unsigned int> one_boundings;
    marray::Marray<unsigned int> two_boundings;
    std::vector< std::vector<point> > arcs;
    std::vector< std::vector<point> > double_pixel_arcs;
    std::vector<point> junctions;
    int dim_x, dim_y; */ 

    std::string filepath_to_ws_image=path_to_ws_image;
    filepath_to_ws_image.append(get_filename(filepath_to_feature_file));
    //remove the ".bin"
    filepath_to_ws_image.resize(filepath_to_ws_image.size()-4);
    filepath_to_ws_image.append(".h5");

    //IMPORT RESULTS FROM HDF5 file
    //all variables represent segmentation and don't need changes here
    seg segment(true);
    segment.load_cgp_data_structure(filepath_to_ws_image);

    vigra::BasicImage<unsigned int> & ws_region_image = segment.ws_region_image;
    marray::Marray<unsigned int> & one_boundings = segment.one_boundings;
    marray::Marray<unsigned int> & two_boundings = segment.two_boundings;
    //std::vector< std::vector<point> > & arcs = segment.arcs;
    std::vector< std::vector<point> > & double_pixel_arcs = segment.arcs;
    std::vector<point> & junctions = segment.junctions;
    int & dim_x = segment.dim_x;
    int & dim_y = segment.dim_y;   

    //Seperate instance of seg for arcs
    seg segment2(false);
    std::vector< std::vector<point> > & arcs = segment2.arcs;

    //we have to load arcs again since double_pixel_arcs contains two pixel arcs
    filepath_to_ws_image.append("a");//dirty
    segment2.load_arcs(filepath_to_ws_image, arcs);

	gbn grainBoundNet;
    size_t & nr_areas = grainBoundNet.nr_new_areas;
    long * & bubble_area_size = grainBoundNet.bubble_area_size;
    long * & grain_area_size = grainBoundNet.grain_area_size;
    std::vector< std::vector< std::vector<int> > > grain_arc_index;
    std::vector< std::vector<int> > & old_grain_arc_index = grainBoundNet.grain_arc_index;
    std::vector< std::vector<int> > & bubble_arc_index = grainBoundNet.bubble_arc_index;
    std::vector<size_t> & region_labels = grainBoundNet.region_labels;
    std::vector<int> & grain_junctions = grainBoundNet.grain_junctions;
    std::vector<int> & grain_bubble_junctions = grainBoundNet.grain_bubble_junctions;
    std::vector<point> & grain_area_center_mass = grainBoundNet.grain_area_center_mass;
    std::vector<point> & bubble_area_center_mass = grainBoundNet.bubble_area_center_mass;
    std::vector<point> reduced_bubble_area_center_mass;
    std::vector<bool> & grain_arc = grainBoundNet.grain_arc;
    std::vector<bool> & subgrain = grainBoundNet.subgrain_arcs;

    std::string param_file_name=get_filename(param_file);
    if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
    else param_file_name = "";

    std::stringstream s;
    if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
    else s << "." << minimal_grain_size;
    std::string suffix=s.str();

    filepath_rf_predictions.append(param_file_name.c_str());
    filepath_to_ws_image.resize(filepath_to_ws_image.size()-4);
    filepath_rf_predictions.append(get_filename(filepath_to_ws_image));
    if (minimal_bubble_distance>0 || minimal_grain_size>0) filepath_rf_predictions.append(suffix.c_str());
    filepath_rf_predictions.append(".h5");

    //check if file exists
    FILE *new_classification;
    new_classification=fopen(filepath_rf_predictions.c_str(),"rb");
    if(new_classification==NULL)//file does NOT exist
    {
        std::cout<<"Error: IceGrain prediction "<<filepath_rf_predictions<<" does not exist."<<std::endl;
        exit(-1);
    }
    else fclose(new_classification);

    //IMPORT RESULTS FROM HDF5 file
    //all variables need to be updated!!
	grainBoundNet.load_final_structure(filepath_rf_predictions);

    //datastructures needed for merging/visualization
    std::vector< std::vector<point> > areas(nr_areas);
    std::vector<area_range> area_ranges(nr_areas);
    std::vector<int> found_bubble_areas;
    std::list<int> found_border_areas;
    std::vector<int> arc_state(arcs.size());
    std::vector< std::vector<point> > grain_boundary_pixels;
    
    //if parameters have been extracted these data structures can be loaded from file
    //otherwise the "junction_and_outer_circle"-function is called
    param para;
    std::vector< std::vector<int> > & grain_boundary_index = para.grain_boundary_index;
    std::vector< std::vector<float> > & grain_boundary_phis = para.grain_boundary_phis;
    std::vector< std::vector<float> > & grain_boundary_curvs = para.grain_boundary_curvs;    

    //create mapping from new to old labels
    std::vector< std::vector<int> > region_labels_reverse;

    for(int old_label=0; old_label<region_labels.size(); old_label++)
    {
        if (region_labels[old_label]>region_labels_reverse.size())
        {
            region_labels_reverse.resize(region_labels[old_label]);
        }

        region_labels_reverse[region_labels[old_label]-1].push_back(old_label+1);
    }

    //set arc state
    for(int arc=0; arc<arcs.size(); arc++)
    {
        if (arc<grain_arc.size())
        {
            if(grain_arc[arc]) arc_state[arc]=1;//grain boundary
            else arc_state[arc]=0;//no boundary
        }
        else
        {
            arc_state[arc]=0;//no boundary
            grain_arc.push_back(false);
        }

        if (arc>=subgrain.size())
        {
            subgrain.push_back(false);
        }
    }

    //set found bubble areas
    for (int area=0; area<nr_areas; area++)
    {
        if (bubble_arc_index[area].size()>0)
        {
            found_bubble_areas.push_back(area);
        }

        for (int a=0; a<bubble_arc_index[area].size(); a++)      
        {
            arc_state[bubble_arc_index[area][a]]=2;//bubble boundary
        }

        //check for area outside the selected area
        if (bubble_arc_index[area].size()==0 && old_grain_arc_index[area].size()==0)
        {
            found_border_areas.push_back(area+1);//outside
        }
    }

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

    for(int area=0; area<nr_areas; area++)
    {
        area_ranges[area].x_low=dim_x;
        area_ranges[area].y_low=dim_y;
        area_ranges[area].x_high=0;
        area_ranges[area].y_high=0;
    }

    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
            int area=region_labels[ws_region_image(x,y)-1];

            //fill areas vector
            point p;
            p.x=x;
            p.y=y;
            areas[area-1].push_back(p);

            //fill areas ranges vector
            if(p.x<area_ranges[area-1].x_low) area_ranges[area-1].x_low=p.x;
            if(p.y<area_ranges[area-1].y_low) area_ranges[area-1].y_low=p.y;
            if(p.x>area_ranges[area-1].x_high) area_ranges[area-1].x_high=p.x;
            if(p.y>area_ranges[area-1].y_high) area_ranges[area-1].y_high=p.y;
        }
    }

    //create reduced vector of bubble center of mass positions
    for (int i=0; i<bubble_area_center_mass.size(); i++)
    {
        if(bubble_area_center_mass[i].x!=0 || bubble_area_center_mass[i].y!=0)
        {
            reduced_bubble_area_center_mass.push_back(bubble_area_center_mass[i]);
        }
    }

    std::string filepath_parameters=path_parameters;
    filepath_parameters.append(get_filename(filepath_to_ws_image));
    filepath_parameters.append(suffix.c_str());
    filepath_parameters.append(".h5");

    //check if file exists
    FILE *extracted_params;
    extracted_params=fopen(filepath_parameters.c_str(),"rb");

    if(extracted_params!=NULL)//file exists
    {
        fclose(extracted_params);
        std::vector<int> & grain_areas = para.grain_areas;
        std::vector<float> & grain_roundness = para.grain_roundness;
        std::vector<float> & grain_box_flattening = para.grain_box_flattening;
        std::vector<float> & grain_box_width = para.grain_box_width;
        std::vector<float> & grain_box_height = para.grain_box_height;
        std::vector< std::vector<float> > & ellipse_params = para.ellipse_params;
        std::vector<float> & ellipse_long_axis = para.ellipse_long_axis;
        std::vector<float> & ellipse_flattening = para.ellipse_flattening;
        std::vector<float> & ellipse_long_axis_angle = para.ellipse_long_axis_angle;
        std::vector<float> & grain_area_width = para.grain_area_width;
        std::vector<float> & grain_area_height = para.grain_area_height;
        std::vector<float> & grain_arc_number = para.grain_arc_number;
        std::vector<float> & grain_neighbors = para.grain_neighbors;
        std::vector<float> & grain_longest_arc_length = para.grain_longest_arc_length;
        std::vector< std::vector<int> > grain_junction_index0(junctions.size());
		para.grain_junction_index = grain_junction_index0;
		std::vector< std::vector<int> > & grain_junction_index = para.grain_junction_index;
        std::vector<float> & turning_point = para.turning_point;
        std::vector< std::vector<unsigned int> > & grain_area_boundaries = para.grain_area_boundaries;

        para.grain_junctions=grainBoundNet.grain_junctions;//gbn attribute is required as param attribute
        para.grain_junction.resize(segment.junctions.size(),false);//seg attribute is required to resize param attribute
        para.load_extracted_parameters(filepath_parameters);

        grain_boundary_pixels.resize(grain_boundary_index.size());

        for(int boundary=0; boundary<grain_boundary_index.size(); boundary++)
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
    }
    else
    {
        std::cout<<"Extracted parameters not found!"<<std::endl;
        std::cout<<"Find outer circles, combine arcs to segments, calculate curvature"<<std::endl;

        std::vector<int> close_bubble_areas;
        std::vector< std::vector<int> > grain_area_junctions(nr_areas);
        std::vector<bool> grain_junction(one_boundings.shape(0),false);
        std::vector<std::vector<int> > arc_junctions(two_boundings.shape(0));
        std::vector<float> grain_perimeter(nr_areas,0.0f);
        std::vector<float> grain_perimeter2(nr_areas,0.0f);
        std::vector<int> min_bubble_distance(nr_areas);
        std::vector<float> grain_longest_arc_length(nr_areas);
        std::vector< std::vector<float> > grain_junction_angles(nr_areas);
        std::vector< std::vector<float> > grain_junction_angles2(nr_areas);
        std::vector< std::vector<int> > grain_junction_index(junctions.size());
        std::vector< std::vector<unsigned int> > grain_area_boundaries(nr_areas);

        //this vector is used to know which arcs have been used for segment combination
        std::vector<bool> arc_to_segment(arcs.size(),false);

        //all junctions for all arcs
        for(int y=0;y<(int)one_boundings.shape(0);y++)
        {
            for(int x=0;x<(int)one_boundings.shape(1);x++)
            {
                if (one_boundings(y,x)!=0) arc_junctions[one_boundings(y,x)-1].push_back(y);//start and end junction of every arc
            }
        }

        //find an outer circle for each grain and combine arcs to segments, calculate curvature of segments
        junctions_and_outer_circle(nr_areas,
                                   grain_area_size,
                                   grain_arc_index,
                                   grain_junctions,
                                   grain_bubble_junctions,
                                   grain_area_center_mass,
                                   reduced_bubble_area_center_mass,
                                   grain_arc,
                                   0,
                                   0,
                                   one_boundings,
                                   two_boundings,
                                   arcs,
                                   junctions,
                                   15,
                                   areas,
                                   area_ranges,
                                   found_bubble_areas,
                                   found_border_areas,
                                   close_bubble_areas,
                                   grain_area_junctions,
                                   grain_junction,
                                   arc_junctions,
                                   grain_perimeter,
                                   grain_perimeter2,
                                   min_bubble_distance,
                                   grain_longest_arc_length,
                                   grain_junction_angles,
                                   grain_junction_angles2,
                                   grain_boundary_pixels,
                                   grain_boundary_index,
                                   grain_boundary_phis,
                                   grain_boundary_curvs,
                                   grain_area_boundaries,
                                   grain_junction_index,
                                   arc_to_segment,
                                   ws_region_image,
                                   region_labels);

        std::cout<<"...done"<<std::endl;
    }
 
    color_type * color = new color_type[7];
    color[0][0] = 0;    color[0][1] = 255;  color[0][2] = 0;
    color[1][0] = 255;  color[1][1] = 0;    color[1][2] = 0;
    color[2][0] = 0;    color[2][1] = 100;  color[2][2] = 255;
    color[3][0] = 255;  color[3][1] = 255;  color[3][2] = 255;//white
    color[4][0] = 0;    color[4][1] = 0;    color[4][2] = 0;//black
    color[5][0] = 255;  color[5][1] = 255;  color[5][2] = 0;//subgrain
    color[6][0] = 0;    color[6][1] = 255;  color[6][2] = 255;//uncertain

    //images for GUI, will be updated
    std::string filepath_to_image=path_to_image;
    filepath_to_image.append(get_filename(filepath_to_feature_file));
    //remove the ".bin"
    filepath_to_image.resize(filepath_to_image.size()-4);

    //temp CImg file is used to avoid an error
	cimg_library::CImg<unsigned char> temp_image(filepath_to_image.c_str());
    cimg_library::CImg<unsigned char> image=temp_image;
    cimg_library::CImg<unsigned char> unmarked_image=temp_image;

    for (int arc=0; arc<arcs.size(); arc++)
        if (grain_arc[arc])
            for (int p=0; p<arcs[arc].size(); p++)
            {
                image(arcs[arc][p].x,arcs[arc][p].y,0,0)=255;
                image(arcs[arc][p].x,arcs[arc][p].y,0,1)=255;
                image(arcs[arc][p].x,arcs[arc][p].y,0,2)=255;
            }

    for(int area=0; area<nr_areas; area++)
        for (int arc=0; arc<bubble_arc_index[area].size(); arc++)
            for (int p=0; p<arcs[bubble_arc_index[area][arc]].size(); p++)
            {
                image(arcs[bubble_arc_index[area][arc]][p].x,arcs[bubble_arc_index[area][arc]][p].y,0,0)=0;
                image(arcs[bubble_arc_index[area][arc]][p].x,arcs[bubble_arc_index[area][arc]][p].y,0,1)=0;
                image(arcs[bubble_arc_index[area][arc]][p].x,arcs[bubble_arc_index[area][arc]][p].y,0,2)=255;
            }

    cimg_library::CImg<unsigned char> original_image=image;//all grain boundaries white

	Parameter<int> display_x;
	display_x.assign("", "display_x", 900);
    display_x.load(paramFile,"config");

	Parameter<int> display_y;
	display_y.assign("", "display_y", 600);
    display_y.load(paramFile,"config");

    display_x=std::min(display_x(),original_image.dimx());
    display_y=std::min(display_y(),original_image.dimy());

    int nr_grain_boundaries=grain_boundary_index.size();
    std::vector<area_range> boundary_ranges(nr_grain_boundaries);

    for(int boundary=0; boundary<nr_grain_boundaries; boundary++)
    {
        boundary_ranges[boundary].x_low=dim_x;
        boundary_ranges[boundary].y_low=dim_y;
        boundary_ranges[boundary].x_high=0;
        boundary_ranges[boundary].y_high=0;
    }

    for(int boundary=0; boundary<nr_grain_boundaries; boundary++)
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

                    //fill boundary ranges vector
                    if(po.x<boundary_ranges[boundary].x_low) boundary_ranges[boundary].x_low=po.x;
                    if(po.y<boundary_ranges[boundary].y_low) boundary_ranges[boundary].y_low=po.y;
                    if(po.x>boundary_ranges[boundary].x_high) boundary_ranges[boundary].x_high=po.x;
                    if(po.y>boundary_ranges[boundary].y_high) boundary_ranges[boundary].y_high=po.y;
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

                    //fill boundary ranges vector
                    if(po.x<boundary_ranges[boundary].x_low) boundary_ranges[boundary].x_low=po.x;
                    if(po.y<boundary_ranges[boundary].y_low) boundary_ranges[boundary].y_low=po.y;
                    if(po.x>boundary_ranges[boundary].x_high) boundary_ranges[boundary].x_high=po.x;
                    if(po.y>boundary_ranges[boundary].y_high) boundary_ranges[boundary].y_high=po.y;
                }
            }
        }
    }

    //*****************************
    //grain boundary classification
    //*****************************

    //this is correct even if subgrain boundaries have been loaded from file as they are not part of GBN
    std::vector<bool> subgrain_boundaries(nr_grain_boundaries, false);
    std::vector<int> uncertain_boundaries;

    //calculate grain boundary features
    std::vector<float> mean_abs_curvature(nr_grain_boundaries);
    std::vector<float> quantile_abs_curvature(nr_grain_boundaries);
    std::vector<float> max_abs_curvature(nr_grain_boundaries);
    std::vector<float> mean_grayvalue(nr_grain_boundaries);
    std::vector<float> stdabw_grayvalue(nr_grain_boundaries);
    std::vector<float> quantile_grayvalue(nr_grain_boundaries);
    //std::vector<float> mean_cross_grayvalue_diff(nr_grain_boundaries);

    for (int boundary=0; boundary<nr_grain_boundaries; boundary++)
    {
        float curv_sum=0.0f;
        int grayvalue_sum=0;
        std::vector<float> abs_curv;
        std::vector<float> grayvalues;
        float curv_max=0.0f;
        //float sum_cross_section_diff=0.0f;
        //int cross_section_pixels=0;

        for (int p=0; p<grain_boundary_pixels[boundary].size(); p++)
        {
            int x=grain_boundary_pixels[boundary][p].x;
            int y=grain_boundary_pixels[boundary][p].y;

            curv_sum+=fabs(grain_boundary_curvs[boundary][p]);
            abs_curv.push_back(fabs(grain_boundary_curvs[boundary][p]));
            grayvalue_sum+=unmarked_image(x,y,0,0);
            grayvalues.push_back(unmarked_image(x,y,0,0));
            if(fabs(grain_boundary_curvs[boundary][p])>curv_max && grain_boundary_pixels[boundary].size()>10) curv_max=fabs(grain_boundary_curvs[boundary][p]);
            /*
            for(int k=-2; k<3; k++)
            {
                int xx =  (int)(x + cos(grain_boundary_phis[boundary][p])*k);
                int yy =  (int)(y + sin(grain_boundary_phis[boundary][p])*k);

                if((xx != x || yy != y) && xx >= 0 && xx + 1 < dim_x && yy >= 0 && yy + 1 < dim_y)
                {
                    cross_section_pixels++;
                    sum_cross_section_diff+=(unmarked_image(xx,yy,0,0)-unmarked_image(x,y,0,0));
                }
            }*/
        }

        mean_abs_curvature[boundary]=curv_sum/(float)grain_boundary_pixels[boundary].size();
        quantile_abs_curvature[boundary]=get_quantile(abs_curv,0.9);
        max_abs_curvature[boundary]=curv_max;
        mean_grayvalue[boundary]=(float)grayvalue_sum/(float)grain_boundary_pixels[boundary].size();
        quantile_grayvalue[boundary]=get_quantile(grayvalues,0.75);
        //mean_cross_grayvalue_diff[boundary]=sum_cross_section_diff/cross_section_pixels;

        float grayvalue_stdabw=0.0f;

        for (int p=0; p<grain_boundary_pixels[boundary].size(); p++)
        {
            int x=grain_boundary_pixels[boundary][p].x;
            int y=grain_boundary_pixels[boundary][p].y;

            grayvalue_stdabw+=square(unmarked_image(x,y,0,0)-mean_grayvalue[boundary]);
        }

        stdabw_grayvalue[boundary]=sqrt(grayvalue_stdabw*(1.0f/((float)grain_boundary_pixels[boundary].size()-1.0f)));
    }

    //predict probabilities, necessary or is hard coded solution enough?

    //TODO

    //if there are high subgrain probabilities classify as subgrain boundary
    for (int boundary=0; boundary<nr_grain_boundaries; boundary++)
    {
//        if ((max_abs_curvature[boundary]>0.06f && mean_grayvalue[boundary]>85.0f*0.6666f) ||
//            quantile_grayvalue[boundary]>85.0f || mean_grayvalue[boundary]>85.0f)
        if /*(quantile_grayvalue[boundary]>120.0f || */(mean_grayvalue[boundary]>140.0f)
        {
            subgrain_boundaries[boundary]=true;
            remove_subgrain_arc(fabs(grain_boundary_index[boundary][0])-1,
                                two_boundings,
                                nr_areas,
                                bubble_area_size,
                                grain_area_size,
                                old_grain_arc_index,
                                bubble_arc_index,
                                region_labels,
                                grain_area_center_mass,
                                grain_arc,
                                subgrain,
                                areas,
                                found_bubble_areas,
                                arc_state,
                                region_labels_reverse);
        }
        else if (false) uncertain_boundaries.push_back(boundary);
    }

    //*************************************************************
    //show uncertain boundaries, classification is used as labeling
    //*************************************************************

    if (uncertain_boundaries.size()>0)
    {
        int startnr_subgrains=get_nr_true(subgrain_boundaries);
        std::vector<bool> labeled(nr_grain_boundaries, false);

        // Initialization of menu
        cimg_library::CImg<unsigned char> menu(700,42,1,3,0);
        cimg_library::CImg<unsigned char> men=menu;

        int boundary_index=0;
        int selected_boundary=uncertain_boundaries[boundary_index];

        int diff_x=boundary_ranges[selected_boundary].x_high-boundary_ranges[selected_boundary].x_low;
        int diff_y=boundary_ranges[selected_boundary].y_high-boundary_ranges[selected_boundary].y_low;

        float scaling=1.0f;
        while (diff_x+50>scaling*display_x || diff_y+50>scaling*display_y) scaling+=0.1f;

        // canvas for our gui to draw on
        // CImg<type> name(dimx, dimy, dimz, colors)
        cimg_library::CImg<unsigned char> canvas(display_x, display_y, 1, 3);

        int posx=std::max(0.0f,boundary_ranges[selected_boundary].x_low - (scaling*display_x-diff_x)/2.0f);
        int posy=std::max(0.0f,boundary_ranges[selected_boundary].y_low - (scaling*display_y-diff_y)/2.0f);

        if (posx+scaling*display_x>image.dimx() && posy+scaling*display_y>image.dimy()) canvas.assign(image.dimx()-posx,image.dimy()-posy, 1, 3);
        else if (posy+scaling*display_y>image.dimy()) canvas.assign(display_x,image.dimy()-posy, 1, 3);
        else if (posx+scaling*display_x>image.dimx()) canvas.assign(image.dimx()-posx,display_y, 1, 3);

        for (int p=0; p<grain_boundary_pixels[selected_boundary].size(); p++)
        {
            int x=grain_boundary_pixels[selected_boundary][p].x;
            int y=grain_boundary_pixels[selected_boundary][p].y;
            image(x,y,0,0)=0;
            image(x,y,0,1)=255;
            image(x,y,0,2)=255;
        }

        // draw the current selection on canvas
        cimg_forXY(canvas,x,y)
        {
            if( x + posx >= 0 && x + posx <  image.dimx() && y + posy >= 0 && y + posy < image.dimy() )
            {
                canvas(x,y,0,0) = image(scaling*x + posx,scaling*y + posy,0,0);
                canvas(x,y,0,1) = image(scaling*x + posx,scaling*y + posy,0,1);
                canvas(x,y,0,2) = image(scaling*x + posx,scaling*y + posy,0,2);
            }
        }

        if (scaling>1.0f)
            for (int p=0; p<grain_boundary_pixels[selected_boundary].size(); p++)
            {
                int x=std::max(0,grain_boundary_pixels[selected_boundary][p].x-posx)/scaling;
                int y=std::max(0,grain_boundary_pixels[selected_boundary][p].y-posy)/scaling;
                canvas(x,y,0,0)=0;
                canvas(x,y,0,1)=255;
                canvas(x,y,0,2)=255;
            }

        cimg_library::CImgDisplay main_menu(menu,"Menu");

        // show the current canvas
        cimg_library::CImgDisplay selectionDisplay(canvas,"View single grain boundary");

        main_menu.move(0,0);
        selectionDisplay.move(0,main_menu.dimy()+57);

        bool end = false;
        bool unmarked = false;
        bool trainingset = false;

        while(!end)
        {
       		//this disables to change the menu size
    		    if (main_menu.is_resized)
            {
    			    main_menu.resize(main_menu);
            }

            //MENU STUFF
            men=menu;

            std::ostringstream Str;
            Str << selected_boundary;
            std::string boundary_string(Str.str());

            men.draw_text(2,0,"Classify boundary or change selection!",color[3],0,1,11,1,1);
            if (!labeled[selected_boundary]) men.draw_text(150,10,"undefined",color[6],0,1,11,1,1);
            else if (subgrain_boundaries[selected_boundary]) men.draw_text(150,10,"sub-grain boundary",color[5],0,1,11,1,1);
            else men.draw_text(150,10,"grain boundary",color[1],0,1,11,1,1);
            if (unmarked) men.draw_text(280,10, "Boundaries are not shown",color[3],0,1,11,1,1);
            else men.draw_text(280,10, "Boundaries are shown",color[3],0,1,11,1,1);
            if (trainingset) men.draw_text(280,30, "Trainingset saved",color[0],0,1,11,1,1);
            else men.draw_text(280,30, "Trainingset not saved",color[1],0,1,11,1,1);

            std::string name="Selected boundary: ";
            name.append(boundary_string);
            men.draw_text(2,10,name.c_str(),color[3],0,1,11,1,1);

            {
                std::ostringstream Str;
                Str << scaling;
                std::string temp_string(Str.str());
                std::string name="Scaling: ";
                name.append(temp_string);
                men.draw_text(280,0,name.c_str(),color[0],0,1,11,1,1);
            }

            {   
                std::ostringstream Str;

                for (int arc=0; arc<grain_boundary_index[selected_boundary].size(); arc++)
                {
                    if (arc<10) Str << fabs(grain_boundary_index[selected_boundary][arc])-1 << " ";
                }
                if (grain_boundary_index[selected_boundary].size()>10) Str << "...";

                std::string temp_string(Str.str());
                std::string name="Arcs: ";
                name.append(temp_string);
                men.draw_text(2,20,name.c_str(),color[0],0,1,11,1,1);
            }

            {   
                std::ostringstream Str;
            
                Str << get_nr_true(labeled)+startnr_subgrains-get_nr_true(subgrain_boundaries);
                std::string temp_string(Str.str());
                std::string name="Grain boundaries: ";
                name.append(temp_string);
                men.draw_text(2,30,name.c_str(),color[1],0,1,11,1,1);
            }

            {   
                std::ostringstream Str;
                Str << get_nr_true(subgrain_boundaries)-startnr_subgrains;
                std::string temp_string(Str.str());
                std::string name="SG boundaries: ";
                name.append(temp_string);
                men.draw_text(150,30,name.c_str(),color[5],0,1,11,1,1);
            }

            {   
                std::ostringstream Str;
                Str << std::fixed << std::setprecision(2) << 100.0f * mean_abs_curvature[selected_boundary] <<"/100 [rad/pixel]";
                std::string temp_string(Str.str());
                std::string name="Mean curvature: ";
                name.append(temp_string);
                men.draw_text(440,0,name.c_str(),color[3],0,1,11,1,1);
            }

            {   
                std::ostringstream Str;
                Str << std::fixed << std::setprecision(2) << 100.0f * max_abs_curvature[selected_boundary] <<"/100 [rad/pixel]";
                std::string temp_string(Str.str());
                std::string name="Max curvature: ";
                name.append(temp_string);
                men.draw_text(440,10,name.c_str(),color[3],0,1,11,1,1);
            }

            {   
                std::ostringstream Str;
                Str << std::fixed << std::setprecision(2) << mean_grayvalue[selected_boundary];
                std::string temp_string(Str.str());
                std::string name="Mean grayvalue: ";
                name.append(temp_string);
                men.draw_text(440,20,name.c_str(),color[3],0,1,11,1,1);
            }

            {   
                std::ostringstream Str;
                Str << std::fixed << std::setprecision(2) << quantile_grayvalue[selected_boundary];
                std::string temp_string(Str.str());
                std::string name="Grayvalue quantile: ";
                name.append(temp_string);
                men.draw_text(440,30,name.c_str(),color[3],0,1,11,1,1);
            }

            main_menu.display(men);

            canvas.display(selectionDisplay);

            selectionDisplay.wait();

            // change selected grain boundary
            if (selectionDisplay.is_keyARROWRIGHT)
            {
                boundary_index++;
                if (boundary_index>uncertain_boundaries.size()-1) boundary_index=uncertain_boundaries.size()-1;
                selected_boundary=uncertain_boundaries[boundary_index];
            }
            if (selectionDisplay.is_keyARROWLEFT)
            {
                boundary_index--;
                if (boundary_index<0) boundary_index=0;
                selected_boundary=uncertain_boundaries[boundary_index];
            }

            // press "b" to change selection for marked boundaries
        		if(selectionDisplay.key==98 || selectionDisplay.key==66)
            {
                if (!unmarked) unmarked=true;
                else unmarked=false;
            }

            // press "q" for quit
        		if(selectionDisplay.key==113 || selectionDisplay.key==81 || selectionDisplay.is_closed)
            {
                for (int boundary=0; boundary<labeled.size(); boundary++)
                {
                    if (subgrain_boundaries[boundary])
                        remove_subgrain_arc(fabs(grain_boundary_index[boundary][0])-1,
                                            two_boundings,
                                            nr_areas,
                                            bubble_area_size,
                                            grain_area_size,
                                            old_grain_arc_index,
                                            bubble_arc_index,
                                            region_labels,
                                            grain_area_center_mass,
                                            grain_arc,
                                            subgrain,
                                            areas,
                                            found_bubble_areas,
                                            arc_state,
                                            region_labels_reverse);
                }

                end = true;
            }

            // press "s" for subgrain boundary
        		if(selectionDisplay.key==115 || selectionDisplay.key==83)
            {
                if (labeled[selected_boundary]==false || subgrain_boundaries[selected_boundary]==false) trainingset=false;
                labeled[selected_boundary]=true;
                subgrain_boundaries[selected_boundary]=true;

            }  

            // press "g" for grain boundary
        		if(selectionDisplay.key==103 || selectionDisplay.key==71)
            {
                if (labeled[selected_boundary]==false || subgrain_boundaries[selected_boundary]==true) trainingset=false;
                labeled[selected_boundary]=true;
                subgrain_boundaries[selected_boundary]=false;
            }  

            // press "c" for clear
        		if(selectionDisplay.key==99 || selectionDisplay.key==67)
            {
                if (labeled[selected_boundary]==true) trainingset=false;
                labeled[selected_boundary]=false;
                subgrain_boundaries[selected_boundary]=false;
            }  

            // press "t" for trainingset
        		if(selectionDisplay.key==116 || selectionDisplay.key==84)
            {
                //create trainingset
                if (get_nr_true(labeled)>0)
                {
                    std::string out_filename=get_filename(filepath_to_ws_image);//TODO set path
                    out_filename.append(".dat");

                    std::ofstream training_file(out_filename.c_str());

                    training_file << "boundary  label  mean_abs_curv  max_abs_cur  mean_gray  quantile_gray  cross_gray_diff\n";

                    for (int boundary=0; boundary<nr_grain_boundaries; boundary++)
                    {
                        if (subgrain_boundaries[boundary])
                        {
                            training_file << boundary << " sGB " << mean_abs_curvature[boundary] << " " << max_abs_curvature[boundary] << " " <<
                                mean_grayvalue[boundary] << " " << quantile_grayvalue[boundary] << /*" " << mean_cross_grayvalue_diff[boundary] <<*/ "\n";
                        }
                        else if (labeled[boundary])
                        {
                            training_file << boundary << "  GB " << mean_abs_curvature[boundary] << " " << max_abs_curvature[boundary] << " " <<
                                mean_grayvalue[boundary] << " " << quantile_grayvalue[boundary] << /*" " << mean_cross_grayvalue_diff[boundary] <<*/ "\n";
                        }
                    }

                    training_file.close();
                    trainingset=true;
                }  
            }
            
            //SELECTION DISPLAY
            diff_x=boundary_ranges[selected_boundary].x_high-boundary_ranges[selected_boundary].x_low;
            diff_y=boundary_ranges[selected_boundary].y_high-boundary_ranges[selected_boundary].y_low;

            scaling=1.0f;
            while (diff_x+50>scaling*display_x || diff_y+50>scaling*display_y) scaling+=0.1f;

            if (!unmarked) image=original_image;
            else image=unmarked_image;

            posx=std::max(0.0f,boundary_ranges[selected_boundary].x_low - (scaling*display_x-diff_x)/2.0f);
            posy=std::max(0.0f,boundary_ranges[selected_boundary].y_low - (scaling*display_y-diff_y)/2.0f);

            if (posx+scaling*display_x>image.dimx() && posy+scaling*display_y>image.dimy()) canvas.assign(image.dimx()-posx,image.dimy()-posy,1,3);
            else if (posy+scaling*display_y>image.dimy()) canvas.assign(display_x,image.dimy()-posy,1,3);
            else if (posx+scaling*display_x>image.dimx()) canvas.assign(image.dimx()-posx,display_y,1,3);
            else canvas.assign(display_x,display_y,1,3);

            if (!unmarked)
            {
                for (int p=0; p<grain_boundary_pixels[selected_boundary].size(); p++)
                {
                    int x=grain_boundary_pixels[selected_boundary][p].x;
                    int y=grain_boundary_pixels[selected_boundary][p].y;

                    if (!labeled[selected_boundary])
                    {
                        image(x,y,0,0)=0;
                        image(x,y,0,1)=255;
                        image(x,y,0,2)=255;
                    }
                    else if (subgrain_boundaries[selected_boundary])
                    {
                        image(x,y,0,0)=255;
                        image(x,y,0,1)=255;
                        image(x,y,0,2)=0;
                    }
                    else
                    {
                        image(x,y,0,0)=255;
                        image(x,y,0,1)=0;
                        image(x,y,0,2)=0;
                    }
                }
            }

            cimg_forXY(canvas,x,y)
            {
                if( x + posx >= 0 && x + posx <  image.dimx() && y + posy >= 0 && y + posy < image.dimy() )
                {
                    canvas(x,y,0,0) = image(scaling*x + posx,scaling*y + posy,0,0);
                    canvas(x,y,0,1) = image(scaling*x + posx,scaling*y + posy,0,1);
                    canvas(x,y,0,2) = image(scaling*x + posx,scaling*y + posy,0,2);
                }
            }

            if (!unmarked && scaling>1.0f)
            {
                for (int p=0; p<grain_boundary_pixels[selected_boundary].size(); p++)
                {
                    int x=std::max(0,grain_boundary_pixels[selected_boundary][p].x-posx)/scaling;
                    int y=std::max(0,grain_boundary_pixels[selected_boundary][p].y-posy)/scaling;

                    if (!labeled[selected_boundary])
                    {
                        canvas(x,y,0,0)=0;
                        canvas(x,y,0,1)=255;
                        canvas(x,y,0,2)=255;
                    }
                    else if (subgrain_boundaries[selected_boundary])
                    {
                        canvas(x,y,0,0)=255;
                        canvas(x,y,0,1)=255;
                        canvas(x,y,0,2)=0;
                    }
                    else
                    {
                        canvas(x,y,0,0)=255;
                        canvas(x,y,0,1)=0;
                        canvas(x,y,0,2)=0;
                    }
                }
            }

            selectionDisplay.resize(canvas.dimx(),canvas.dimy());
        }
    }

    //***************************************
    //calculation of new structure and saving
    //***************************************

    std::cout<<"Calculation of new structure and saving..."<<std::endl;

    //create new junction lists
    grain_junctions.clear();
    grain_bubble_junctions.clear();
    
    for(int y=0;y<(int)one_boundings.shape(0);y++)//loop over junctions
    {
        int nr_grain_arcs=0;
        int nr_bubble_arcs=0;

        for(int x=0;x<(int)one_boundings.shape(1);x++)
        {
            bool found=false;

            for(size_t area=0;area<nr_areas && !found;area++)
            {
                for (int a=0;a<(int)bubble_arc_index[area].size() && !found;a++)//bubble boundaries
                {
                    int arc_index=bubble_arc_index[area][a]+1;
                    if (one_boundings(y,x)==arc_index)
                    {
                        nr_bubble_arcs++;
                        found=true;
                    }
                }
            }

            for(size_t area=0;area<nr_areas && !found;area++)
            {
                for (int a=0;a<(int)old_grain_arc_index[area].size() && !found;a++)//grain boundaries
                {
                    int arc_index=old_grain_arc_index[area][a]+1;
                    if (one_boundings(y,x)==arc_index)
                    {
                        nr_grain_arcs++;
                        found=true;
                    }
                }
            }

        }
        //grain junctions to vector and to output image
        if (nr_grain_arcs==3 || nr_grain_arcs==4)
        {
            grain_junctions.push_back(y);
        }
        //grain bubble junction to vector and to output image
        if (nr_bubble_arcs==2 && (nr_grain_arcs==1 || nr_grain_arcs==2))
        {
            grain_bubble_junctions.push_back(y);
        }
    }
    
    //create output images
    vigra::FImage grain_image(dim_x,dim_y);//grain areas gray, bubbles and grain boundaries black

    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
            grain_image(x,y)=3;
        }
    }

    //loop over areas, bubble/grain areas and grain center of mass to output image
    for(size_t area=0;area<nr_areas;area++)
    {
        if (bubble_area_size[area]>0)//bubble area
        {
            for(size_t k=0; k<areas[area].size(); ++k)
            {
                size_t x=areas[area][k].x;
                size_t y=areas[area][k].y;
                grain_image(x,y)=0;
            }
        }
        else if (grain_area_size[area]>0)//grain area
        {
            for(size_t k=0; k<areas[area].size(); ++k)
            {
                size_t x=areas[area][k].x;
                size_t y=areas[area][k].y;
                grain_image(x,y)=1;
            }

            //center of mass positions
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

    //loop over arc lists, boundaries of bubble/grain boundaries to output image
    for (int area=0;area<nr_areas;area++)//loop over areas
    {
        for (int a=0;a<(int)old_grain_arc_index[area].size();a++)//grain boundaries
        {
            //now we loop over the points in this arc
            for(int p=0;p<(int)double_pixel_arcs[old_grain_arc_index[area][a]].size();p++)
            {
                int x=double_pixel_arcs[old_grain_arc_index[area][a]][p].x;             
                int y=double_pixel_arcs[old_grain_arc_index[area][a]][p].y;
                grain_image(x,y)=0;
            }
        }
    }

    //if subgrain boundaries are labeled show them in output grain image
    for(int arc=0; arc<double_pixel_arcs.size(); arc++)
    {
        if (subgrain[arc])
        {
            for (int i=0; i<double_pixel_arcs[arc].size(); i++)
            {
                int x=double_pixel_arcs[arc][i].x;
                int y=double_pixel_arcs[arc][i].y;
                grain_image(x,y)=2;
            }
        }
    }

    std::string filepath_grains=path_subGB_rf_predictions;
    filepath_rf_predictions=path_subGB_rf_predictions;

    filepath_grains.append("grains/");

    filepath_rf_predictions.append(param_file_name.c_str());
    filepath_grains.append(param_file_name.c_str());

    filepath_grains.append(get_filename(filepath_to_ws_image));
    filepath_rf_predictions.append(get_filename(filepath_to_ws_image));
    if (minimal_bubble_distance>0 || minimal_grain_size>0) filepath_rf_predictions.append(suffix.c_str());

    filepath_grains.append(".jpg");
                    
    exportImage(srcImageRange(grain_image), vigra::ImageExportInfo(filepath_grains.c_str()));

    grainBoundNet.paramFileGBN = &paramFile;
	grainBoundNet.save_final_structure(filepath_rf_predictions);

    std::cout<<"...done"<<std::endl;
}
