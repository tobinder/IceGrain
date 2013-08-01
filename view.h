/*! \file view.h
 * \brief Used for network visualization.
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
void grain_menu(cimg_library::CImg<unsigned char> & men, cimg_library::CImg<unsigned char> & output_men, int selected_grain, bool only_inner_grains,
    std::vector<bool> inner_grain_areas, std::string grain_string, int unmarked, color_type * color, int diff_x, int diff_y, float scaling,
    std::vector< std::vector<point> > arcs, int dim_x, int dim_y, long * grain_area_size, std::vector< std::vector<int> >  grain_arc_index,
    std::vector<point> grain_area_center_mass, std::vector<int> grain_areas, std::vector<float> grain_roundness, std::vector<float> grain_box_flattening,
    std::vector<float> grain_box_width, std::vector<float> grain_box_height, std::vector< std::vector<float> > ellipse_params, std::vector<float> ellipse_long_axis,
    std::vector<float> ellipse_flattening, std::vector<float> ellipse_long_axis_angle, std::vector<float> grain_arc_number, std::vector<float> grain_neighbors,
    std::vector<float> grain_longest_arc_length);

void show_ellipse_box(cimg_library::CImg<unsigned char> & image, int selected_grain, int unmarked, float angle, area_range & range, FitEllipse fitEllipse,
    std::vector<point> & ellipse_points, std::vector< std::vector<point> > arcs, int dim_x, int dim_y, std::vector< std::vector<int> > grain_arc_index,
    std::vector<point> grain_area_center_mass, std::vector<bool> grain_arc, std::vector<int> grain_areas, std::vector< std::vector<float> > ellipse_params);

void print_grain(float scaling, int dim_x, int dim_y, cimg_library::CImg<unsigned char> output_men, cimg_library::CImg<unsigned char> image, int posx, int posy,
    std::string filepath_plots, int unmarked, std::string grain_string, std::string suffix);

void boundary_menu(cimg_library::CImg<unsigned char> & men, cimg_library::CImg<unsigned char> & output_men, int selected_boundary, bool only_inner_boundaries,
    std::vector<bool> inner_boundary, std::string boundary_string, int unmarked, color_type * color, float scaling,
    std::vector< std::vector<int> > grain_boundary_index, std::vector< std::vector<point> > grain_boundary_pixels, std::vector<float> turning_point);

void print_boundary(float scaling, int dim_x, int dim_y, cimg_library::CImg<unsigned char> output_men, cimg_library::CImg<unsigned char> image, int posx, int posy,
    std::string filepath_plots, int unmarked, std::string boundary_string, int selected_boundary, std::vector< std::vector<float> > grain_boundary_curvs, plplot plot,
    std::string suffix);

void junction_menu(cimg_library::CImg<unsigned char> & men, cimg_library::CImg<unsigned char> & output_men, int selected_junction, std::string junction_string,
    int unmarked, color_type * color, float scaling, std::vector<int> grain_junctions, std::vector<point> junctions, std::vector< std::vector<int> > grain_junction_index,
    int pixels_average, point * average_pos);

void show_orientation(cimg_library::CImg<unsigned char> & image, std::vector<int> grain_junctions, std::vector<point> junctions, int dim_x, int dim_y,
    std::vector< std::vector<int> > grain_junction_index, std::vector< std::vector<point> > grain_boundary_pixels, int selected_junction, int posx, int posy,
    int pixels_average, point * average_pos, int dimx, int dimy);

void print_junction(float scaling, int dim_x, int dim_y, cimg_library::CImg<unsigned char> output_men, cimg_library::CImg<unsigned char> image, int posx, int posy,
    std::string filepath_plots, int unmarked, std::string junction_string, std::string suffix);

void view(std::string filepath_to_feature_file, std::string path_to_image, std::string path_to_ws_image, std::string path_rf_predictions,
          std::string path_plots, std::string path_parameters, std::string param_file, ParameterFile paramFile, int minimal_grain_size, 
          int minimal_bubble_distance, int mode, int start_value=-1, bool auto_print=false, int mark=0)
{
    gbn grainBoundNet;
    size_t & nr_areas =                                     grainBoundNet.nr_new_areas;
    long * & bubble_area_size =                             grainBoundNet.bubble_area_size;
    long * & grain_area_size =                              grainBoundNet.grain_area_size;
    std::vector< std::vector<int> > & grain_arc_index =     grainBoundNet.grain_arc_index;
    std::vector< std::vector<int> > & bubble_arc_index =    grainBoundNet.bubble_arc_index;
    std::vector<size_t> & region_labels =                   grainBoundNet.region_labels;
    std::vector<int> & grain_junctions =                    grainBoundNet.grain_junctions;
    std::vector<int> & grain_bubble_junctions =             grainBoundNet.grain_bubble_junctions;
    std::vector<point> & grain_area_center_mass =           grainBoundNet.grain_area_center_mass;
    std::vector<point> & bubble_area_center_mass =          grainBoundNet.bubble_area_center_mass;
    std::vector<bool> & grain_arc =                         grainBoundNet.grain_arc;

    std::string param_file_name=get_filename(param_file);
    if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
    else param_file_name = "";

    Parameter<int> display_x;
    display_x.assign("", "display_x", 900);
    display_x.load(paramFile,"config");

    Parameter<int> display_y;
    display_y.assign("", "display_y", 600);
    display_y.load(paramFile,"config");

    std::string filepath_new_classification=path_rf_predictions;
    filepath_new_classification.append(param_file_name.c_str());

    std::string filename=get_filename(filepath_to_feature_file);
    //remove the ".bin"
    filename.resize(filename.size()-4);

    std::stringstream s;
    if(minimal_bubble_distance>0) s << "." << minimal_bubble_distance << "." << minimal_grain_size;
    else s << "." << minimal_grain_size;
    std::string suffix=s.str();

    filepath_new_classification.append(filename);
    if (minimal_grain_size>0 || minimal_bubble_distance>0) filepath_new_classification.append(suffix.c_str());
    filepath_new_classification.append(".h5");

    std::string filepath_to_image=path_to_image;
    filepath_to_image.append(filename);

    std::string filepath_plots=path_plots;
    filepath_plots.append(filename);

    //IMPORT RESULTS FROM HDF5 file
    grainBoundNet.load_final_structure(filepath_new_classification);

    std::string filepath_to_ws_image=path_to_ws_image;
    filepath_to_ws_image.append(filename);
    filepath_to_ws_image.append(".h5");

    //IMPORT RESULTS FROM HDF5 file
    seg segment(false);
    segment.load_cgp_data_structure(filepath_to_ws_image);
    
    vigra::BasicImage<unsigned int> & ws_region_image = segment.ws_region_image;
    marray::Marray<unsigned int> & one_boundings =      segment.one_boundings;
    marray::Marray<unsigned int> & two_boundings =      segment.two_boundings;
    std::vector< std::vector<point> > & arcs =          segment.arcs;
    std::vector<point> & junctions =                    segment.junctions;
    int & dim_x =                                       segment.dim_x;
    int & dim_y =                                       segment.dim_y;    
    
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
    filepath_parameters.append(filename);
    filepath_parameters.append(suffix.c_str());
    filepath_parameters.append(".h5");

    param para;
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
    para.grain_junction_index.resize(junctions.size());
    std::vector< std::vector<int> > & grain_junction_index = para.grain_junction_index;
    std::vector<float> & turning_point = para.turning_point;
    std::vector< std::vector<int> > & grain_boundary_index = para.grain_boundary_index;    
    std::vector< std::vector<float> > & grain_boundary_phis = para.grain_boundary_phis;
    std::vector< std::vector<float> > & grain_boundary_curvs = para.grain_boundary_curvs;

    FitEllipse fitEllipse;//initialize fit class

    para.grain_junctions=grainBoundNet.grain_junctions;//gbn attribute is required as param attribute
    para.grain_junction.resize(segment.junctions.size(),false);//seg attribute is required to resize param attribute
    para.load_extracted_parameters(filepath_parameters);

    if (grain_junctions.size()==0 && (mode==2 || mode==3))
    {
        std::cout<<"Error: No junctions found!";
        return;
    }

    color_type * color = new color_type[5];
    color[0][0] = 0;    color[0][1] = 255;  color[0][2] = 0;
    color[1][0] = 255;  color[1][1] = 0;    color[1][2] = 0;
    color[2][0] = 0;    color[2][1] = 100;  color[2][2] = 255;
    color[3][0] = 255;  color[3][1] = 255;  color[3][2] = 255;//white
    color[4][0] = 0;    color[4][1] = 0;    color[4][2] = 0;//black

    int posx=0;
    int posy=0;

    //temp CImg file is used to avoid an error
    cimg_library::CImg<unsigned char> temp_image(filepath_to_image.c_str());
    cimg_library::CImg<unsigned char> image=temp_image;
    cimg_library::CImg<unsigned char> unmarked_image=temp_image;

    std::string check_image=filename;
    check_image.resize(4);

    if (check_image=="GRIP" || check_image=="NGRI")
        for (int x=0; x<image.dimx(); x++)
            for (int y=0; y<image.dimy(); y++)
            {
                image(x,y,0,0)=0;
                image(x,y,0,1)=0;
                image(x,y,0,2)=0;

            }

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

    //Import the original gray value image
    vigra::ImageImportInfo info(filepath_to_image.c_str());
    vigra::BasicImage <unsigned int> gray_image;
    gray_image.resize(image.dimx(),image.dimy());

    if (check_image=="GRIP" || check_image=="NGRI")
        for (int x=0; x<image.dimx(); x++)
            for (int y=0; y<image.dimy(); y++)
            {
                gray_image(x,y)=0;
            }
    else importImage(info, destImage(gray_image));

    if (mode==0) //show grains
    {
        std::vector<area_range> area_ranges(nr_areas);

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

                point p;
                p.x=x;
                p.y=y;

                //fill areas ranges vector
                if(p.x<area_ranges[area-1].x_low) area_ranges[area-1].x_low=p.x;
                if(p.y<area_ranges[area-1].y_low) area_ranges[area-1].y_low=p.y;
                if(p.x>area_ranges[area-1].x_high) area_ranges[area-1].x_high=p.x;
                if(p.y>area_ranges[area-1].y_high) area_ranges[area-1].y_high=p.y;
            }
        }

        bool only_inner_grains = true;
        std::vector<bool> inner_grain_areas(grain_areas.size(), true);

        for (int grain=0; grain<grain_areas.size(); grain++)
            if(grain_arc_number[grain_areas[grain]]==0) inner_grain_areas[grain]=false;

        int selected_grain=0;
        int selected_inner_grain=0;

        if (start_value==-1)
        {
            while (!inner_grain_areas[selected_inner_grain] && selected_inner_grain<grain_areas.size()-1) selected_inner_grain++;

            if (!inner_grain_areas[selected_inner_grain]) //no inner grain
            {
                selected_inner_grain=-1;
                only_inner_grains = false;
            }
            else selected_grain=selected_inner_grain;
        }
        else
        {
            bool found = false;
            for (int i=0; i<grain_areas.size() && !found; i++)
                if (grain_areas[i]==start_value)
                {
                    selected_grain=i;
                    found=true;
                }

            if (found)
            {
                if (!inner_grain_areas[selected_grain])
                {
                    only_inner_grains = false;

                    while (!inner_grain_areas[selected_inner_grain] && selected_inner_grain<grain_areas.size()-1) selected_inner_grain++;

                    if (!inner_grain_areas[selected_inner_grain]) //no inner grain
                    {
                        selected_inner_grain=-1;
                    }
                }
                else selected_inner_grain=selected_grain;
            }
            else
            {
                std::cout<<"Error! Grain "<<start_value<< " does not exist!"<<std::endl;
                exit(-1);
            }
        }

        display_x=std::min(display_x(),image.dimx());
        display_y=std::min(display_y(),image.dimy());

        int diff_x=area_ranges[grain_areas[selected_grain]].x_high-area_ranges[grain_areas[selected_grain]].x_low;
        int diff_y=area_ranges[grain_areas[selected_grain]].y_high-area_ranges[grain_areas[selected_grain]].y_low;

        float scaling=1.0f;
        while (diff_x+50>scaling*display_x || diff_y+50>scaling*display_y) scaling+=0.1f;

        // canvas for our gui to draw on
        // CImg<type> name(dimx, dimy, dimz, colors)
        cimg_library::CImg<unsigned char> canvas(display_x, display_y, 1, 3);

        posx=std::max(0.0f,area_ranges[grain_areas[selected_grain]].x_low - (scaling*display_x-diff_x)/2.0f);
        posy=std::max(0.0f,area_ranges[grain_areas[selected_grain]].y_low - (scaling*display_y-diff_y)/2.0f);

        if (posx+scaling*display_x>image.dimx() && posy+scaling*display_y>image.dimy()) canvas.assign(image.dimx()-posx,image.dimy()-posy, 1, 3);
        else if (posy+scaling*display_y>image.dimy()) canvas.assign(display_x,image.dimy()-posy, 1, 3);
        else if (posx+scaling*display_x>image.dimx()) canvas.assign(image.dimx()-posx,display_y, 1, 3);

        for (int arc=0; arc<grain_arc_index[grain_areas[selected_grain]].size(); arc++)
            for (int p=0; p<arcs[grain_arc_index[grain_areas[selected_grain]][arc]].size(); p++)
            {
                int x=arcs[grain_arc_index[grain_areas[selected_grain]][arc]][p].x;
                int y=arcs[grain_arc_index[grain_areas[selected_grain]][arc]][p].y;
                if (grain_arc[grain_arc_index[grain_areas[selected_grain]][arc]])
                {
                    image(x,y,0,0)=255;
                    image(x,y,0,1)=0;
                    image(x,y,0,2)=0;
                }
                else
                {
                    image(x,y,0,0)=0;
                    image(x,y,0,1)=0;
                    image(x,y,0,2)=255;
                }
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
        {
            for (int arc=0; arc<grain_arc_index[grain_areas[selected_grain]].size(); arc++)
                for (int p=0; p<arcs[grain_arc_index[grain_areas[selected_grain]][arc]].size(); p++)
                {
                    int x=std::max(0,arcs[grain_arc_index[grain_areas[selected_grain]][arc]][p].x-posx)/scaling;
                    int y=std::max(0,arcs[grain_arc_index[grain_areas[selected_grain]][arc]][p].y-posy)/scaling;
                    if (grain_arc[grain_arc_index[grain_areas[selected_grain]][arc]])
                    {
                        canvas(x,y,0,0)=255;
                        canvas(x,y,0,1)=0;
                        canvas(x,y,0,2)=0;
                    }
                    else
                    {
                        canvas(x,y,0,0)=0;
                        canvas(x,y,0,1)=0;
                        canvas(x,y,0,2)=255;
                    }
                }
        }

        // Initialization of menu
        cimg_library::CImg<unsigned char> menu(700,52,1,3,0);
        cimg_library::CImg<unsigned char> men=menu;
        cimg_library::CImg<unsigned char> output_menu(700,32,1,3,0);
        cimg_library::CImg<unsigned char> output_men=output_menu;

        if (auto_print)
        {
            std::ostringstream Str;
            Str << grain_areas[selected_grain];
            std::string grain_string(Str.str());

            grain_menu(men, output_men, selected_grain, only_inner_grains, inner_grain_areas, grain_string, mark, color, diff_x, diff_y, scaling,
                arcs, dim_x, dim_y, grain_area_size, grain_arc_index, grain_area_center_mass, grain_areas, grain_roundness,
                grain_box_flattening, grain_box_width, grain_box_height, ellipse_params, ellipse_long_axis, ellipse_flattening, ellipse_long_axis_angle,
                grain_arc_number, grain_neighbors, grain_longest_arc_length);

            area_range range;
            range.x_low=dim_x;
            range.y_low=dim_y;
            range.x_high=0;
            range.y_high=0;

            std::vector<point> ellipse_points;
            float angle=0.0f;
            if (mark!=4) angle=ellipse_long_axis_angle[grain_areas[selected_grain]];

            //print unmarked
            print_grain(scaling, canvas.dimx(), canvas.dimy(), output_men, unmarked_image, posx, posy, filepath_plots, 5, grain_string, suffix);

            //print marked
            if (mark<5)
            {
                show_ellipse_box(image, selected_grain, mark, angle, range, fitEllipse, ellipse_points, arcs, dim_x, dim_y, grain_arc_index,
                    grain_area_center_mass, grain_arc, grain_areas, ellipse_params);

                print_grain(scaling, canvas.dimx(), canvas.dimy(), output_men, image, posx, posy, filepath_plots, mark, grain_string, suffix);
            }

            return;
        }

        cimg_library::CImgDisplay main_menu(menu,"Menu");

        // show the current canvas
        cimg_library::CImgDisplay selectionDisplay(canvas,"View single grain");

        main_menu.move(0,0);
        selectionDisplay.move(0,main_menu.dimy()+57);

        bool end = false;
        int unmarked = 0;

        while(!end)
        {
            area_range range;
            range.x_low=dim_x;
            range.y_low=dim_y;
            range.x_high=0;
            range.y_high=0;

            //this disables to change the menu size
            if (main_menu.is_resized)
            {
                main_menu.resize(main_menu);
            }

            //MENU STUFF
            men=menu;
            output_men=output_menu;

            std::ostringstream Str;
            Str << grain_areas[selected_grain];
            std::string grain_string(Str.str());

            grain_menu(men, output_men, selected_grain, only_inner_grains, inner_grain_areas, grain_string, unmarked, color, diff_x, diff_y, scaling,
                arcs, dim_x, dim_y, grain_area_size, grain_arc_index, grain_area_center_mass, grain_areas, grain_roundness,
                grain_box_flattening, grain_box_width, grain_box_height, ellipse_params, ellipse_long_axis, ellipse_flattening, ellipse_long_axis_angle,
                grain_arc_number, grain_neighbors, grain_longest_arc_length);

            main_menu.display(men);
            canvas.display(selectionDisplay);

            selectionDisplay.wait();

            // change selected grain
            if (selectionDisplay.is_keyARROWRIGHT)
            {
                if (only_inner_grains)
                {
                    selected_inner_grain++;
                    while (!inner_grain_areas[selected_inner_grain] && selected_inner_grain<grain_areas.size()-1) selected_inner_grain++;
                    if (!inner_grain_areas[selected_inner_grain]) selected_inner_grain=selected_grain;
                    selected_grain=selected_inner_grain;
                }
                else
                {
                    selected_grain++;
                    if (selected_grain>grain_areas.size()-1) selected_grain=grain_areas.size()-1;
                    if (inner_grain_areas[selected_grain]) selected_inner_grain=selected_grain;
                }
            }
            if (selectionDisplay.is_keyARROWLEFT)
            {
                if (only_inner_grains)
                {
                    selected_inner_grain--;
                    while (!inner_grain_areas[selected_inner_grain] && selected_inner_grain>0) selected_inner_grain--;
                    if (!inner_grain_areas[selected_inner_grain]) selected_inner_grain=selected_grain;
                    selected_grain=selected_inner_grain;
                }
                else
                {
                    selected_grain--;
                    if (selected_grain<0) selected_grain=0;
                    if (inner_grain_areas[selected_grain]) selected_inner_grain=selected_grain;
                }
            }

            // press "b" to change selection for marked boundaries
                if(selectionDisplay.key==98 || selectionDisplay.key==66)
            {
                unmarked=(unmarked+1)%6;
            }

            // press "i" to change selection for inner grains only
            if(selectionDisplay.key==105 || selectionDisplay.key==73)
            {
                if (!only_inner_grains && selected_inner_grain!=-1)
                {
                    only_inner_grains=true;
                    selected_grain=selected_inner_grain;
                }
                else only_inner_grains=false;
            }

            // press "q" for quit
            if(selectionDisplay.key==113 || selectionDisplay.key==81 || selectionDisplay.is_closed)
            {
                end = true;
            }      

            // press "ENTER" to print
            if (selectionDisplay.is_keyENTER) print_grain(scaling, selectionDisplay.dimx(), selectionDisplay.dimy(), output_men, image, posx, posy,
                filepath_plots, unmarked, grain_string, suffix);

            //SELECTION DISPLAY
            diff_x=area_ranges[grain_areas[selected_grain]].x_high-area_ranges[grain_areas[selected_grain]].x_low;
            diff_y=area_ranges[grain_areas[selected_grain]].y_high-area_ranges[grain_areas[selected_grain]].y_low;

            scaling=1.0f;
            while (diff_x+50>scaling*display_x || diff_y+50>scaling*display_y) scaling+=0.1f;

            if (unmarked<5) image=original_image;
            else image=unmarked_image;

            posx=std::max(0.0f,area_ranges[grain_areas[selected_grain]].x_low - (scaling*display_x-diff_x)/2.0f);
            posy=std::max(0.0f,area_ranges[grain_areas[selected_grain]].y_low - (scaling*display_y-diff_y)/2.0f);

            if (posx+scaling*display_x>image.dimx() && posy+scaling*display_y>image.dimy()) canvas.assign(image.dimx()-posx,image.dimy()-posy,1,3);
            else if (posy+scaling*display_y>image.dimy()) canvas.assign(display_x,image.dimy()-posy,1,3);
            else if (posx+scaling*display_x>image.dimx()) canvas.assign(image.dimx()-posx,display_y,1,3);
            else canvas.assign(display_x,display_y,1,3);

            std::vector<point> ellipse_points;
            float angle=0.0f;
            if (unmarked!=4) angle=ellipse_long_axis_angle[grain_areas[selected_grain]];

            if (unmarked<5) show_ellipse_box(image, selected_grain, unmarked, angle, range, fitEllipse, ellipse_points, arcs, dim_x, dim_y, grain_arc_index,
                grain_area_center_mass, grain_arc, grain_areas, ellipse_params);

            cimg_forXY(canvas,x,y)
            {
                if( x + posx >= 0 && x + posx <  image.dimx() && y + posy >= 0 && y + posy < image.dimy() )
                {
                    canvas(x,y,0,0) = image(scaling*x + posx,scaling*y + posy,0,0);
                    canvas(x,y,0,1) = image(scaling*x + posx,scaling*y + posy,0,1);
                    canvas(x,y,0,2) = image(scaling*x + posx,scaling*y + posy,0,2);
                }
            }

            if (unmarked<5 && scaling>1.0f)
            {
                for (int arc=0; arc<grain_arc_index[grain_areas[selected_grain]].size(); arc++)
                    for (int p=0; p<arcs[grain_arc_index[grain_areas[selected_grain]][arc]].size(); p++)
                    {
                        int x=std::max(0,arcs[grain_arc_index[grain_areas[selected_grain]][arc]][p].x-posx)/scaling;
                        int y=std::max(0,arcs[grain_arc_index[grain_areas[selected_grain]][arc]][p].y-posy)/scaling;
                        if (grain_arc[grain_arc_index[grain_areas[selected_grain]][arc]])
                        {
                            canvas(x,y,0,0)=255;
                            canvas(x,y,0,1)=0;
                            canvas(x,y,0,2)=0;
                        }
                        else
                        {
                            canvas(x,y,0,0)=0;
                            canvas(x,y,0,1)=0;
                            canvas(x,y,0,2)=255;
                        }
                    }

                if (unmarked==1 || unmarked==3 || unmarked==4)
                {
                    //show box
                    for(int x=range.x_low; x<=range.x_high; x++)
                    {
                        point pixel;
                        pixel.x=x-grain_area_center_mass[grain_areas[selected_grain]].x;
                        pixel.y=range.y_low-grain_area_center_mass[grain_areas[selected_grain]].y;

                        int xx=(grain_area_center_mass[grain_areas[selected_grain]].x+
                            pixel.x*cos(PI*angle/180.0f)+
                            pixel.y*sin(PI*angle/180.0f)-posx)/scaling;
                        int yy=(grain_area_center_mass[grain_areas[selected_grain]].y-
                            pixel.x*sin(PI*angle/180.0f)+
                            pixel.y*cos(PI*angle/180.0f)-posy)/scaling;

                        if (xx>=0 && xx<canvas.dimx() && yy>=0 && yy<canvas.dimy())
                        {
                            canvas(xx,yy,0,0)=0;
                            canvas(xx,yy,0,1)=255;
                            canvas(xx,yy,0,2)=0;
                        }

                        pixel.x=x-grain_area_center_mass[grain_areas[selected_grain]].x;
                        pixel.y=range.y_high-grain_area_center_mass[grain_areas[selected_grain]].y;

                        xx=(grain_area_center_mass[grain_areas[selected_grain]].x+
                            pixel.x*cos(PI*angle/180.0f)+
                            pixel.y*sin(PI*angle/180.0f)-posx)/scaling;
                        yy=(grain_area_center_mass[grain_areas[selected_grain]].y-
                            pixel.x*sin(PI*angle/180.0f)+
                            pixel.y*cos(PI*angle/180.0f)-posy)/scaling;

                        if (xx>=0 && xx<canvas.dimx() && yy>=0 && yy<canvas.dimy())
                        {
                            canvas(xx,yy,0,0)=0;
                            canvas(xx,yy,0,1)=255;
                            canvas(xx,yy,0,2)=0;
                        }
                    }

                    for(int y=range.y_low; y<=range.y_high; y++)
                    {
                        point pixel;
                        pixel.x=range.x_low-grain_area_center_mass[grain_areas[selected_grain]].x;
                        pixel.y=y-grain_area_center_mass[grain_areas[selected_grain]].y;

                        int xx=(grain_area_center_mass[grain_areas[selected_grain]].x+
                            pixel.x*cos(PI*angle/180.0f)+
                            pixel.y*sin(PI*angle/180.0f)-posx)/scaling;
                        int yy=(grain_area_center_mass[grain_areas[selected_grain]].y-
                            pixel.x*sin(PI*angle/180.0f)+
                            pixel.y*cos(PI*angle/180.0f)-posy)/scaling;

                        if (xx>=0 && xx<canvas.dimx() && yy>=0 && yy<canvas.dimy())
                        {
                            canvas(xx,yy,0,0)=0;
                            canvas(xx,yy,0,1)=255;
                            canvas(xx,yy,0,2)=0;
                        }
                        pixel.x=range.x_high-grain_area_center_mass[grain_areas[selected_grain]].x;
                        pixel.y=y-grain_area_center_mass[grain_areas[selected_grain]].y;

                        xx=(grain_area_center_mass[grain_areas[selected_grain]].x+
                            pixel.x*cos(PI*angle/180.0f)+
                            pixel.y*sin(PI*angle/180.0f)-posx)/scaling;
                        yy=(grain_area_center_mass[grain_areas[selected_grain]].y-
                            pixel.x*sin(PI*angle/180.0f)+
                            pixel.y*cos(PI*angle/180.0f)-posy)/scaling;

                        if (xx>=0 && xx<canvas.dimx() && yy>=0 && yy<canvas.dimy())
                        {
                            canvas(xx,yy,0,0)=0;
                            canvas(xx,yy,0,1)=255;
                            canvas(xx,yy,0,2)=0;
                        }
                    }
                }

                if ((unmarked==2 || unmarked==3) && fabs(1.0f-ellipse_params[grain_areas[selected_grain]][5])<0.1f)
                {
                    for (int i=0; i<ellipse_points.size(); i++)
                    {
                        if (ellipse_points[i].x-posx>=0 && (ellipse_points[i].x-posx)/scaling<canvas.dimx() &&
                            ellipse_points[i].y-posy>=0 && (ellipse_points[i].y-posy)/scaling<canvas.dimy())
                        {
                            canvas((ellipse_points[i].x-posx)/scaling,(ellipse_points[i].y-posy)/scaling,0,0)=0;
                            canvas((ellipse_points[i].x-posx)/scaling,(ellipse_points[i].y-posy)/scaling,0,1)=255;
                            canvas((ellipse_points[i].x-posx)/scaling,(ellipse_points[i].y-posy)/scaling,0,2)=0;
                        }
                    }
                }
            }

            selectionDisplay.resize(canvas.dimx(),canvas.dimy());
        }
    }

    else if (mode==1) //show grain boundaries
    {
        //initialise plplot class
        plplot plot = plplot();

        int nr_grain_boundaries=grain_boundary_index.size();

        std::vector< std::vector<point> > grain_boundary_pixels(nr_grain_boundaries);
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

                        grain_boundary_pixels[boundary].push_back(po);

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

                        grain_boundary_pixels[boundary].push_back(po);

                        //fill boundary ranges vector
                        if(po.x<boundary_ranges[boundary].x_low) boundary_ranges[boundary].x_low=po.x;
                        if(po.y<boundary_ranges[boundary].y_low) boundary_ranges[boundary].y_low=po.y;
                        if(po.x>boundary_ranges[boundary].x_high) boundary_ranges[boundary].x_high=po.x;
                        if(po.y>boundary_ranges[boundary].y_high) boundary_ranges[boundary].y_high=po.y;
                    }
                }
            }
        }

        bool only_inner_boundaries = true;
        std::vector<bool> inner_boundary(grain_boundary_index.size(), true);

        for (int boundary=0; boundary<grain_boundary_index.size(); boundary++)
        {
            int arc_index=fabs(grain_boundary_index[boundary][0])-1;

            if ((grain_area_size[twoCellNeighbors(arc_index,1)-1]==0 && bubble_area_size[twoCellNeighbors(arc_index,1)-1]==0) ||
                (grain_area_size[twoCellNeighbors(arc_index,0)-1]==0 && bubble_area_size[twoCellNeighbors(arc_index,0)-1]==0))
                inner_boundary[boundary]=false;
        }

        int selected_boundary=0;
        int selected_inner_boundary=0;

        if (start_value==-1)
        {
            while (!inner_boundary[selected_inner_boundary] && selected_inner_boundary<grain_boundary_index.size()-1) selected_inner_boundary++;

            if (!inner_boundary[selected_inner_boundary]) //no inner boundary
            {
                selected_inner_boundary=-1;
                only_inner_boundaries = false;
            }
            else selected_boundary=selected_inner_boundary;
        }
        else
        {
            selected_boundary=start_value;

            if (selected_boundary<0 || selected_boundary>grain_boundary_index.size()-1)
            {
                std::cout<<"Error! Grain boundary "<<start_value<< " does not exist!"<<std::endl;
                exit(-1);
            }

            if (!inner_boundary[selected_boundary])
            {
                only_inner_boundaries = false;

                while (!inner_boundary[selected_inner_boundary] && selected_inner_boundary<grain_boundary_index.size()-1) selected_inner_boundary++;

                if (!inner_boundary[selected_inner_boundary]) //no inner grain
                {
                    selected_inner_boundary=-1;
                }
            }
            else selected_inner_boundary=selected_boundary;
        }

        display_x=std::min(display_x(),image.dimx());
        display_y=std::min(display_y(),image.dimy());

        int diff_x=boundary_ranges[selected_boundary].x_high-boundary_ranges[selected_boundary].x_low;
        int diff_y=boundary_ranges[selected_boundary].y_high-boundary_ranges[selected_boundary].y_low;

        float scaling=1.0f;
        while (diff_x+50>scaling*display_x || diff_y+50>scaling*display_y) scaling+=0.1f;

        // canvas for our gui to draw on
        // CImg<type> name(dimx, dimy, dimz, colors)
        cimg_library::CImg<unsigned char> canvas(display_x, display_y, 1, 3);

        posx=std::max(0.0f,boundary_ranges[selected_boundary].x_low - (scaling*display_x-diff_x)/2.0f);
        posy=std::max(0.0f,boundary_ranges[selected_boundary].y_low - (scaling*display_y-diff_y)/2.0f);

        if (posx+scaling*display_x>image.dimx() && posy+scaling*display_y>image.dimy()) canvas.assign(image.dimx()-posx,image.dimy()-posy, 1, 3);
        else if (posy+scaling*display_y>image.dimy()) canvas.assign(display_x,image.dimy()-posy, 1, 3);
        else if (posx+scaling*display_x>image.dimx()) canvas.assign(image.dimx()-posx,display_y, 1, 3);

        for (int p=0; p<grain_boundary_pixels[selected_boundary].size(); p++)
        {
            int x=grain_boundary_pixels[selected_boundary][p].x;
            int y=grain_boundary_pixels[selected_boundary][p].y;
            image(x,y,0,0)=255;
            image(x,y,0,1)=0;
            image(x,y,0,2)=0;
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
                canvas(x,y,0,0)=255;
                canvas(x,y,0,1)=0;
                canvas(x,y,0,2)=0;
            }

        // Initialization of menu
        cimg_library::CImg<unsigned char> menu(500,42,1,3,0);
        cimg_library::CImg<unsigned char> men=menu;
        cimg_library::CImg<unsigned char> output_menu(500,22,1,3,0);
        cimg_library::CImg<unsigned char> output_men=output_menu;

        if (auto_print)
        {
            std::ostringstream Str;
            Str << selected_boundary;
            std::string boundary_string(Str.str());

            boundary_menu(men, output_men, selected_boundary, only_inner_boundaries, inner_boundary, boundary_string, mark, color, scaling,
                grain_boundary_index, grain_boundary_pixels, turning_point);

            //print unmarked
            print_boundary(scaling, canvas.dimx(), canvas.dimy(), output_men, unmarked_image, posx, posy, filepath_plots, true, boundary_string,
                selected_boundary, grain_boundary_curvs, plot, suffix);

            //print marked
            if (!mark) print_boundary(scaling, canvas.dimx(), canvas.dimy(), output_men, image, posx, posy, filepath_plots, false, boundary_string,
                selected_boundary, grain_boundary_curvs, plot, suffix);

            return;
        }

        cimg_library::CImgDisplay main_menu(menu,"Menu");
        cimg_library::CImgDisplay curv_disp(250,250,"Curvature",0);
        cimg_library::CImgDisplay cross_disp(250,250,"Cross-section",0);
        cimg_library::CImgDisplay phi_disp(250,250,"Orientation",0);

        // show the current canvas
        cimg_library::CImgDisplay selectionDisplay(canvas,"View single grain boundary");

        main_menu.move(0,0);
        selectionDisplay.move(0,main_menu.dimy()+57);
        phi_disp.move(selectionDisplay.dimx()+10,0);
        curv_disp.move(selectionDisplay.dimx()+10,phi_disp.dimy()+33);
        cross_disp.move(selectionDisplay.dimx()+10,phi_disp.dimy()+curv_disp.dimy()+66);


        bool end = false;
        bool unmarked = false;
        int boundary_pos = grain_boundary_pixels[selected_boundary].size()/2;

        while(!end)
        {
               //this disables to change the menu size
                if (main_menu.is_resized)
            {
                    main_menu.resize(main_menu);
            }
            if (curv_disp.is_resized)
            {
                curv_disp.resize(curv_disp);
            }
            if (cross_disp.is_resized)
            {
                cross_disp.resize(cross_disp);
            }
            if (phi_disp.is_resized)
            {
                phi_disp.resize(phi_disp);
            }

            //MENU STUFF
            men=menu;
            output_men=output_menu;

            std::ostringstream Str;
            Str << selected_boundary;
            std::string boundary_string(Str.str());

            boundary_menu(men, output_men, selected_boundary, only_inner_boundaries, inner_boundary, boundary_string, unmarked, color, scaling,
                grain_boundary_index, grain_boundary_pixels, turning_point);

            main_menu.display(men);

            //CURV PLOT
            cimg_library::CImg<unsigned char> curv_plot(curv_disp.dimx(),curv_disp.dimy(),1,3,0);
            float curv_scaling = (float)grain_boundary_curvs[selected_boundary].size()/(float)(curv_disp.dimx()-3);
            if (curv_scaling<0.5f) curv_scaling = 0.5f;

            float y_high=0.05f;

            for(int i=0; i<grain_boundary_curvs[selected_boundary].size(); i++)
            {
                if(grain_boundary_curvs[selected_boundary][i]*1.1f>y_high || grain_boundary_curvs[selected_boundary][i]*1.1f<-y_high)
                    y_high=fabs(grain_boundary_curvs[selected_boundary][i])*1.1f;
            }

            for (int x=0; x<curv_disp.dimx() && (x-3)*curv_scaling<grain_boundary_curvs[selected_boundary].size(); x++)
            {
                curv_plot(x,0.5f*(float)curv_disp.dimy(),0,0)=255;
                curv_plot(x,0.5f*(float)curv_disp.dimy(),0,1)=255;
                curv_plot(x,0.5f*(float)curv_disp.dimy(),0,2)=255;

                if (x==2) for (int y=0; y<curv_disp.dimy(); y++)
                {
                    curv_plot(x,y,0,0)=255;
                    curv_plot(x,y,0,1)=255;
                    curv_plot(x,y,0,2)=255;
                }
                
                if (x>2) curv_plot(x,0.5f*(float)curv_disp.dimy()*(1.0f-(float)grain_boundary_curvs[selected_boundary][(x-3)*curv_scaling]/y_high),0,1)=255;
            }

            {
                std::ostringstream Str;
                Str << std::fixed << std::setprecision(2) <<  y_high << " [rad/pixel]";
                std::string high(Str.str());
                std::string low="-";
                low.append(high);

                curv_plot.draw_text(5,0,high.c_str(),color[3],0,1,11,1,1);
                curv_plot.draw_text(5,curv_disp.dimy()-11,low.c_str(),color[3],0,1,11,1,1);
            }

            for (int y=(curv_disp.dimy()/2)-5; y<(curv_disp.dimy()/2)+5; y++)
            {
                curv_plot(3+boundary_pos/curv_scaling,y,0,0)=255;
                curv_plot(3+boundary_pos/curv_scaling,y,0,1)=255;
                curv_plot(3+boundary_pos/curv_scaling,y,0,2)=255;
            }

            //CROSS_PLOT
            cimg_library::CImg<unsigned char> cross_plot(cross_disp.dimx(),cross_disp.dimy(),1,3,0);
            std::vector<float> cross;

            for(int i=-10; i<11; i++)
            {
                int x=grain_boundary_pixels[selected_boundary][boundary_pos].x + cos(grain_boundary_phis[selected_boundary][boundary_pos]+PI/2.0f)*i;
                int y=grain_boundary_pixels[selected_boundary][boundary_pos].y + sin(grain_boundary_phis[selected_boundary][boundary_pos]+PI/2.0f)*i;

                cross_plot(25+10*(i+10),(float)cross_disp.dimy()*(1.0f-gray_image(x,y)/255.0f),0,1)=255;

                if (!unmarked && (x-posx-1)/scaling>=0 && (x-posx+1)/scaling<canvas.dimx() && (y-posy-1)/scaling>=0 && (y-posy+1)/scaling<canvas.dimy())
                {
                    canvas((x-posx)/scaling,(y-posy)/scaling,0,0)=0;
                    canvas((x-1-posx)/scaling,(y-posy)/scaling,0,0)=0;
                    canvas((x+1-posx)/scaling,(y-posy)/scaling,0,0)=0;
                    canvas((x-posx)/scaling,(y-1-posy)/scaling,0,0)=0;
                    canvas((x-posx)/scaling,(y+1-posy)/scaling,0,0)=0;

                    canvas((x-posx)/scaling,(y-posy)/scaling,0,1)=0;
                    canvas((x-1-posx)/scaling,(y-posy)/scaling,0,1)=0;
                    canvas((x+1-posx)/scaling,(y-posy)/scaling,0,1)=0;
                    canvas((x-posx)/scaling,(y-1-posy)/scaling,0,1)=0;
                    canvas((x-posx)/scaling,(y+1-posy)/scaling,0,1)=0;

                    canvas((x-posx)/scaling,(y-posy)/scaling,0,2)=255;
                    canvas((x-1-posx)/scaling,(y-posy)/scaling,0,2)=255;
                    canvas((x+1-posx)/scaling,(y-posy)/scaling,0,2)=255;
                    canvas((x-posx)/scaling,(y-1-posy)/scaling,0,2)=255;
                    canvas((x-posx)/scaling,(y+1-posy)/scaling,0,2)=255;
                }

                if (i>-10)
                {
                    int diff=gray_image(x,y)-cross.back();

                    for (int j=1; j<10; j++)
                    {
                        int value = cross.back()+diff*j/10;
                        cross_plot(15+10*(i+10)+j,(float)cross_disp.dimy()*(1.0f-value/255.0f),0,1)=255;
                    }
                }

                cross.push_back(gray_image(x,y));
            }

            for (int y=0; y<cross_disp.dimy(); y++)
            {
                cross_plot(125,y,0,0)=255;
                cross_plot(125,y,0,1)=255;
                cross_plot(125,y,0,2)=255;
            }

            cross_plot.draw_text(128,0,"255",color[3],0,1,11,1,1);
            cross_plot.draw_text(128,cross_disp.dimy()-11,"0",color[3],0,1,11,1,1);

            //PHI PLOT
            cimg_library::CImg<unsigned char> phi_plot(phi_disp.dimx(),phi_disp.dimy(),1,3,0);
            float phi_scaling = (float)grain_boundary_phis[selected_boundary].size()/(float)(phi_disp.dimx()-3);
            if (phi_scaling<0.5f) phi_scaling = 0.5f;

            y_high=3.15f;

            for(int i=0; i<grain_boundary_phis[selected_boundary].size(); i++)
            {
                if(grain_boundary_phis[selected_boundary][i]*1.1f>y_high || grain_boundary_phis[selected_boundary][i]*1.1f<-y_high)
                    y_high=fabs(grain_boundary_phis[selected_boundary][i])*1.1f;
            }

            for (int x=0; x<phi_disp.dimx() && (x-3)*phi_scaling<grain_boundary_phis[selected_boundary].size(); x++)
            {
                phi_plot(x,0.5f*(float)phi_disp.dimy(),0,0)=255;
                phi_plot(x,0.5f*(float)phi_disp.dimy(),0,1)=255;
                phi_plot(x,0.5f*(float)phi_disp.dimy(),0,2)=255;

                if (x==2) for (int y=0; y<phi_disp.dimy(); y++)
                {
                    phi_plot(x,y,0,0)=255;
                    phi_plot(x,y,0,1)=255;
                    phi_plot(x,y,0,2)=255;
                }
                
                if (x>2) phi_plot(x,0.5f*(float)phi_disp.dimy()*(1.0f-(float)grain_boundary_phis[selected_boundary][(x-3)*phi_scaling]/y_high),0,1)=255;
            }

            {
                std::ostringstream Str;
                Str << std::fixed << std::setprecision(0) <<  y_high*180.0f/PI << " degrees";
                std::string high(Str.str());
                std::string low="-";
                low.append(high);

                phi_plot.draw_text(5,0,high.c_str(),color[3],0,1,11,1,1);
                phi_plot.draw_text(5,phi_disp.dimy()-11,low.c_str(),color[3],0,1,11,1,1);
            }

            for (int y=(phi_disp.dimy()/2)-5; y<(phi_disp.dimy()/2)+5; y++)
            {
                phi_plot(3+boundary_pos/phi_scaling,y,0,0)=255;
                phi_plot(3+boundary_pos/phi_scaling,y,0,1)=255;
                phi_plot(3+boundary_pos/phi_scaling,y,0,2)=255;
            }

            curv_disp.display(curv_plot);
            cross_disp.display(cross_plot);
            phi_disp.display(phi_plot);
            canvas.display(selectionDisplay);

            selectionDisplay.wait();

            // change boundary position for cross section
            if (selectionDisplay.is_keyARROWUP)
            {
                boundary_pos++;
                if (boundary_pos>grain_boundary_pixels[selected_boundary].size()-1) boundary_pos=grain_boundary_pixels[selected_boundary].size()-1;
            }
            if (selectionDisplay.is_keyARROWDOWN)
            {
                boundary_pos--;
                if (boundary_pos<0) boundary_pos=0;
            }

            // change selected grain boundary
            if (selectionDisplay.is_keyARROWRIGHT)
            {
                if (only_inner_boundaries)
                {
                    selected_inner_boundary++;
                    while (!inner_boundary[selected_inner_boundary] && selected_inner_boundary<grain_boundary_index.size()-1) selected_inner_boundary++;
                    if (!inner_boundary[selected_inner_boundary]) selected_inner_boundary=selected_boundary;
                    selected_boundary=selected_inner_boundary;
                }
                else
                {
                    selected_boundary++;
                    if (selected_boundary>grain_boundary_index.size()-1) selected_boundary=grain_boundary_index.size()-1;
                    if (inner_boundary[selected_boundary]) selected_inner_boundary=selected_boundary;
                }

                boundary_pos = grain_boundary_pixels[selected_boundary].size()/2;
            }
            if (selectionDisplay.is_keyARROWLEFT)
            {
                if (only_inner_boundaries)
                {
                    selected_inner_boundary--;
                    while (!inner_boundary[selected_inner_boundary] && selected_inner_boundary>0) selected_inner_boundary--;
                    if (!inner_boundary[selected_inner_boundary]) selected_inner_boundary=selected_boundary;
                    selected_boundary=selected_inner_boundary;
                }
                else
                {
                    selected_boundary--;
                    if (selected_boundary<0) selected_boundary=0;
                    if (inner_boundary[selected_boundary]) selected_inner_boundary=selected_boundary;
                }

                boundary_pos = grain_boundary_pixels[selected_boundary].size()/2;
            }

            // press "b" to change selection for marked boundaries
                if(selectionDisplay.key==98 || selectionDisplay.key==66)
            {
                if (!unmarked) unmarked=true;
                else unmarked=false;
            }

            // press "i" to change selection for inner boundaries only
                if(selectionDisplay.key==105 || selectionDisplay.key==73)
            {
                if (!only_inner_boundaries && selected_inner_boundary!=-1)
                {
                    only_inner_boundaries=true;
                    selected_boundary=selected_inner_boundary;
                }
                else only_inner_boundaries=false;
            }

            // press "q" for quit
                if(selectionDisplay.key==113 || selectionDisplay.key==81 || selectionDisplay.is_closed)
            {
                end = true;
            }      

            // press "ENTER" to print
            if (selectionDisplay.is_keyENTER)
            {
                print_boundary(scaling, selectionDisplay.dimx(), selectionDisplay.dimy(), output_men, image, posx, posy,
                filepath_plots, unmarked, boundary_string, selected_boundary, grain_boundary_curvs, plot, suffix);

                std::string filepath_output=filepath_plots;

                filepath_output.append("_boundary");
                filepath_output.append(boundary_string);
                filepath_output.append("_cross");
                filepath_output.append(suffix.c_str());
                filepath_output.append(".svg");

                std::vector<int> x_values;
                for (int i=-10; i<11; i++) x_values.push_back(i);

                std::cout<<"Export current cross-section to: "<<filepath_output<<std::endl;
                plot.draw_curv_cross("Position [pixel]", "Gray value", "Cross-section", x_values, cross, filepath_output.c_str());
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
                    image(x,y,0,0)=255;
                    image(x,y,0,1)=0;
                    image(x,y,0,2)=0;
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
                    canvas(x,y,0,0)=255;
                    canvas(x,y,0,1)=0;
                    canvas(x,y,0,2)=0;
                }
            }

            selectionDisplay.resize(canvas.dimx(),canvas.dimy());
        }
    }

    else if (mode==2) //show grain triple junctions
    {
        std::vector<bool> grain_junction(one_boundings.shape(0),false);

        //fill grain junction
        for (int j=0; j<grain_junctions.size(); j++)
        {
            grain_junction[grain_junctions[j]]=true;
        }

        int nr_grain_boundaries=grain_boundary_index.size();

        std::vector< std::vector<point> > grain_boundary_pixels(nr_grain_boundaries);
        std::vector<area_range> boundary_ranges(nr_grain_boundaries);
        std::vector<area_range> junction_ranges(grain_junctions.size());

        for(int boundary=0; boundary<nr_grain_boundaries; boundary++)
        {
            boundary_ranges[boundary].x_low=dim_x;
            boundary_ranges[boundary].y_low=dim_y;
            boundary_ranges[boundary].x_high=0;
            boundary_ranges[boundary].y_high=0;
        }

        for(int junction=0; junction<grain_junctions.size(); junction++)
        {
            junction_ranges[junction].x_low=dim_x;
            junction_ranges[junction].y_low=dim_y;
            junction_ranges[junction].x_high=0;
            junction_ranges[junction].y_high=0;
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

                        //fill grain boundary pixels
                        grain_boundary_pixels[boundary].push_back(po);

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

                        grain_boundary_pixels[boundary].push_back(po);

                        //fill boundary ranges vector
                        if(po.x<boundary_ranges[boundary].x_low) boundary_ranges[boundary].x_low=po.x;
                        if(po.y<boundary_ranges[boundary].y_low) boundary_ranges[boundary].y_low=po.y;
                        if(po.x>boundary_ranges[boundary].x_high) boundary_ranges[boundary].x_high=po.x;
                        if(po.y>boundary_ranges[boundary].y_high) boundary_ranges[boundary].y_high=po.y;
                    }
                }
            }
        }

        for(int j=0; j<grain_junctions.size(); j++)
        {
            for (int i=0; i<grain_junction_index[grain_junctions[j]].size(); i++)
            {
                //fill junction ranges vector
                if(boundary_ranges[grain_junction_index[grain_junctions[j]][i]].x_low<junction_ranges[j].x_low)
                    junction_ranges[j].x_low=boundary_ranges[grain_junction_index[grain_junctions[j]][i]].x_low;
                if(boundary_ranges[grain_junction_index[grain_junctions[j]][i]].y_low<junction_ranges[j].y_low)
                    junction_ranges[j].y_low=boundary_ranges[grain_junction_index[grain_junctions[j]][i]].y_low;
                if(boundary_ranges[grain_junction_index[grain_junctions[j]][i]].x_high>junction_ranges[j].x_high)
                    junction_ranges[j].x_high=boundary_ranges[grain_junction_index[grain_junctions[j]][i]].x_high;
                if(boundary_ranges[grain_junction_index[grain_junctions[j]][i]].y_high>junction_ranges[j].y_high)
                    junction_ranges[j].y_high=boundary_ranges[grain_junction_index[grain_junctions[j]][i]].y_high;

                //check which side of grain boundary fits to junction
                point front=grain_boundary_pixels[grain_junction_index[grain_junctions[j]][i]][0];
                point back=grain_boundary_pixels[grain_junction_index[grain_junctions[j]][i]].back();

                if((fabs(front.x-junctions[grain_junctions[j]].x)+fabs(front.y-junctions[grain_junctions[j]].y))<2)
                {
                    grain_junction_index[grain_junctions[j]][i]=grain_junction_index[grain_junctions[j]][i]+1;
                }
                else if((fabs(back.x-junctions[grain_junctions[j]].x)+fabs(back.y-junctions[grain_junctions[j]].y))<2)
                {
                    grain_junction_index[grain_junctions[j]][i]=-grain_junction_index[grain_junctions[j]][i]-1;
                }
                else
                {
                    std::cout<<"Error: Junction do not fit to boundaries!"<<std::endl;
                    exit(-1);
                }
            }
        }

        int selected_junction=0;

        if (start_value!=-1)
        {
            bool found = false;
            for (int i=0; i<grain_junctions.size() && !found; i++)
                if (grain_junctions[i]==start_value)
                {
                    selected_junction=i;
                    found=true;
                }

            if (!found)
            {
                std::cout<<"Error! Grain junction "<<start_value<< " does not exist!"<<std::endl;
                exit(-1);
            }
        }

        display_x=std::min(display_x(),image.dimx());
        display_y=std::min(display_y(),image.dimy());

        int pixels_average=15;
        int * boundary_pos = new int[grain_junction_index[grain_junctions[selected_junction]].size()];
        point * average_pos = new point [grain_junction_index[grain_junctions[selected_junction]].size()];

        for (int i=0; i<grain_junction_index[grain_junctions[selected_junction]].size(); i++)
        {
            boundary_pos[i]=std::min(pixels_average,(int)grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1].size());
        }

        int diff_x=junction_ranges[selected_junction].x_high-junction_ranges[selected_junction].x_low;
        int diff_y=junction_ranges[selected_junction].y_high-junction_ranges[selected_junction].y_low;

        float scaling=1.0f;
        while (diff_x+50>scaling*display_x || diff_y+50>scaling*display_y) scaling+=0.1f;

        // canvas for our gui to draw on
        // CImg<type> name(dimx, dimy, dimz, colors)
        cimg_library::CImg<unsigned char> canvas(display_x, display_y, 1, 3);

        posx=std::max(0.0f,junction_ranges[selected_junction].x_low - (scaling*display_x-diff_x)/2.0f);
        posy=std::max(0.0f,junction_ranges[selected_junction].y_low - (scaling*display_y-diff_y)/2.0f);

        if (posx+scaling*display_x>image.dimx() && posy+scaling*display_y>image.dimy()) canvas.assign(image.dimx()-posx,image.dimy()-posy, 1, 3);
        else if (posy+scaling*display_y>image.dimy()) canvas.assign(display_x,image.dimy()-posy, 1, 3);
        else if (posx+scaling*display_x>image.dimx()) canvas.assign(image.dimx()-posx,display_y, 1, 3);

        //compute the "center" of pixels
        for (int i=0; i<grain_junction_index[grain_junctions[selected_junction]].size(); i++)
        {
            int number=0;
            int sum_x=0;
            int sum_y=0;

            for(int p=0; p<boundary_pos[i]; p++)
            {
                int index;
                if(grain_junction_index[grain_junctions[selected_junction]][i]>0) index=p;
                else index=grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1].size()-1-p;

                sum_x+=grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1][index].x;
                sum_y+=grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1][index].y;
                number++;
            }

            average_pos[i].x=sum_x/number;
            average_pos[i].y=sum_y/number;
        }

        show_orientation(image, grain_junctions, junctions, dim_x, dim_y, grain_junction_index, grain_boundary_pixels, selected_junction, posx, posy, pixels_average,
            average_pos, canvas.dimx(), canvas.dimy());

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
        {
            for (int i=0; i<grain_junction_index[grain_junctions[selected_junction]].size(); i++)
            {
                for (int p=0; p<grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1].size(); p++)
                {
                    int x=std::max(0,grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1][p].x-posx)/scaling;
                    int y=std::max(0,grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1][p].y-posy)/scaling;
                    canvas(x,y,0,0)=255;
                    canvas(x,y,0,1)=0;
                    canvas(x,y,0,2)=0;
                }

                float dx=average_pos[i].x-junctions[grain_junctions[selected_junction]].x;
                float dy=average_pos[i].y-junctions[grain_junctions[selected_junction]].y;

                if (dx!=0)
                {
                    float m=dy/dx;
                    float t=average_pos[i].y-m*average_pos[i].x;

                    if (dx>0)
                    {
                        for (int xx=junctions[grain_junctions[selected_junction]].x; xx<=posx+canvas.dimx(); xx++)
                        {            
                            int yy=m*(float)xx+t;
                            if (yy>=0 && yy<dim_y && square(junctions[grain_junctions[selected_junction]].x-xx)+
                                square(junctions[grain_junctions[selected_junction]].y-yy)<square(pixels_average)) 
                            {
                                canvas((xx-posx)/scaling,(yy-posy)/scaling,0,0)=0;
                                canvas((xx-posx)/scaling,(yy-posy)/scaling,0,1)=255;
                                canvas((xx-posx)/scaling,(yy-posy)/scaling,0,2)=0;
                            }
                        }
                    }
                    else
                    {
                        for (int xx=posx; xx<=junctions[grain_junctions[selected_junction]].x; xx++)
                        {            
                            int yy=m*(float)xx+t;
                            if (yy>=0 && yy<dim_y && square(junctions[grain_junctions[selected_junction]].x-xx)+
                                square(junctions[grain_junctions[selected_junction]].y-yy)<square(pixels_average)) 
                            {
                                canvas((xx-posx)/scaling,(yy-posy)/scaling,0,0)=0;
                                canvas((xx-posx)/scaling,(yy-posy)/scaling,0,1)=255;
                                canvas((xx-posx)/scaling,(yy-posy)/scaling,0,2)=0;
                            }
                        }
                    }
                }            

                if (dy!=0)
                {
                    float m=dx/dy;
                    float t=average_pos[i].x-m*average_pos[i].y;

                    if (dy>0)
                    {
                        for (int yy=junctions[grain_junctions[selected_junction]].y; yy<=posy+canvas.dimy(); yy++)
                        {            
                            int xx=m*(float)yy+t;
                            if (xx>=0 && xx<dim_x && square(junctions[grain_junctions[selected_junction]].x-xx)+
                                square(junctions[grain_junctions[selected_junction]].y-yy)<square(pixels_average)) 
                            {
                                canvas((xx-posx)/scaling,(yy-posy)/scaling,0,0)=0;
                                canvas((xx-posx)/scaling,(yy-posy)/scaling,0,1)=255;
                                canvas((xx-posx)/scaling,(yy-posy)/scaling,0,2)=0;
                            }
                        }
                    }
                    else
                    {
                        for (int yy=posy; yy<=junctions[grain_junctions[selected_junction]].y; yy++)
                        {            
                            int xx=m*(float)yy+t;
                            if (xx>=0 && xx<dim_x && square(junctions[grain_junctions[selected_junction]].x-xx)+
                                square(junctions[grain_junctions[selected_junction]].y-yy)<square(pixels_average)) 
                            {
                                canvas((xx-posx)/scaling,(yy-posy)/scaling,0,0)=0;
                                canvas((xx-posx)/scaling,(yy-posy)/scaling,0,1)=255;
                                canvas((xx-posx)/scaling,(yy-posy)/scaling,0,2)=0;
                            }
                        }
                    }
                }
            }
        }

        // Initialization of menu
        cimg_library::CImg<unsigned char> menu(500,42,1,3,0);
        cimg_library::CImg<unsigned char> men=menu;
        cimg_library::CImg<unsigned char> output_menu(500,22,1,3,0);
        cimg_library::CImg<unsigned char> output_men=output_menu;

        if (auto_print)
        {
            std::ostringstream Str;
            Str << grain_junctions[selected_junction];
            std::string junction_string(Str.str());

            junction_menu(men, output_men, selected_junction, junction_string, mark, color, scaling, grain_junctions, junctions, grain_junction_index,
                pixels_average, average_pos);

            //print unmarked
            print_junction(scaling, canvas.dimx(), canvas.dimy(), output_men, unmarked_image, posx, posy, filepath_plots, true, junction_string, suffix);

            //print marked
            if (!mark) print_junction(scaling, canvas.dimx(), canvas.dimy(), output_men, image, posx, posy, filepath_plots, false, junction_string, suffix);

            return;
        }

        cimg_library::CImgDisplay main_menu(menu,"Menu");

        // show the current canvas
        cimg_library::CImgDisplay selectionDisplay(canvas,"View single grain junction");

        main_menu.move(0,0);
        selectionDisplay.move(0,main_menu.dimy()+57);

        bool end = false;
        bool unmarked = false;

        while(!end)
        {
               //this disables to change the menu size
                if (main_menu.is_resized)
            {
                    main_menu.resize(main_menu);
            }
    
            //MENU STUFF
            men=menu;
            output_men=output_menu;

            std::ostringstream Str;
            Str << grain_junctions[selected_junction];
            std::string junction_string(Str.str());

            junction_menu(men, output_men, selected_junction, junction_string, unmarked, color, scaling, grain_junctions, junctions, grain_junction_index,
                pixels_average, average_pos);

            main_menu.display(men);
            canvas.display(selectionDisplay);

            selectionDisplay.wait();

            // change pixels_average
            if (selectionDisplay.is_keyARROWUP)
            {
                pixels_average++;

                for (int i=0; i<grain_junction_index[grain_junctions[selected_junction]].size(); i++)
                {
                    boundary_pos[i]=std::min(pixels_average,
                        (int)grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1].size());
                }
            }
            if (selectionDisplay.is_keyARROWDOWN)
            {
                pixels_average--;
                if (pixels_average<2) pixels_average=2;

                for (int i=0; i<grain_junction_index[grain_junctions[selected_junction]].size(); i++)
                {
                    boundary_pos[i]=std::min(pixels_average,
                        (int)grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1].size());
                }
            }

            // change selected grain junction
            if (selectionDisplay.is_keyARROWRIGHT)
            {
                selected_junction++;
                if (selected_junction>grain_junctions.size()-1) selected_junction=grain_junctions.size()-1;

                delete boundary_pos;
                delete average_pos;
                boundary_pos = new int[grain_junction_index[grain_junctions[selected_junction]].size()];
                average_pos = new point[grain_junction_index[grain_junctions[selected_junction]].size()];

                for (int i=0; i<grain_junction_index[grain_junctions[selected_junction]].size(); i++)
                {
                    boundary_pos[i]=std::min(pixels_average,
                        (int)grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1].size());
                }
            }

            if (selectionDisplay.is_keyARROWLEFT)
            {
                selected_junction--;
                if (selected_junction<0) selected_junction=0;

                delete boundary_pos;
                delete average_pos;
                boundary_pos = new int[grain_junction_index[grain_junctions[selected_junction]].size()];
                average_pos = new point[grain_junction_index[grain_junctions[selected_junction]].size()];

                for (int i=0; i<grain_junction_index[grain_junctions[selected_junction]].size(); i++)
                {
                    boundary_pos[i]=std::min(pixels_average,
                        (int)grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1].size());
                }
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
                end = true;
            }      

            // press "ENTER" to print
            if (selectionDisplay.is_keyENTER) print_junction(scaling, canvas.dimx(), canvas.dimy(), output_men, image, posx, posy, filepath_plots,
                unmarked, junction_string, suffix);

            //SELECTION DISPLAY
            diff_x=junction_ranges[selected_junction].x_high-junction_ranges[selected_junction].x_low;
            diff_y=junction_ranges[selected_junction].y_high-junction_ranges[selected_junction].y_low;

            scaling=1.0f;
            while (diff_x+50>scaling*display_x || diff_y+50>scaling*display_y) scaling+=0.1f;

            if (!unmarked) image=original_image;
            else image=unmarked_image;

            posx=std::max(0.0f,junction_ranges[selected_junction].x_low - (scaling*display_x-diff_x)/2.0f);
            posy=std::max(0.0f,junction_ranges[selected_junction].y_low - (scaling*display_y-diff_y)/2.0f);

            if (posx+scaling*display_x>image.dimx() && posy+scaling*display_y>image.dimy()) canvas.assign(image.dimx()-posx,image.dimy()-posy,1,3);
            else if (posy+scaling*display_y>image.dimy()) canvas.assign(display_x,image.dimy()-posy,1,3);
            else if (posx+scaling*display_x>image.dimx()) canvas.assign(image.dimx()-posx,display_y,1,3);
            else canvas.assign(display_x,display_y,1,3);

            //compute the "center" of pixels
            for (int i=0; i<grain_junction_index[grain_junctions[selected_junction]].size(); i++)
            {
                int number=0;
                int sum_x=0;
                int sum_y=0;

                for(int p=0; p<boundary_pos[i]; p++)
                {
                    int index;
                    if(grain_junction_index[grain_junctions[selected_junction]][i]>0) index=p;
                    else index=grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1].size()-1-p;

                    sum_x+=grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1][index].x;
                    sum_y+=grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1][index].y;
                    number++;
                }

                average_pos[i].x=sum_x/number;
                average_pos[i].y=sum_y/number;
            }

            if (!unmarked) show_orientation(image, grain_junctions, junctions, dim_x, dim_y, grain_junction_index, grain_boundary_pixels, selected_junction, posx, posy,
                pixels_average, average_pos, selectionDisplay.dimx(), selectionDisplay.dimy());

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
                for (int i=0; i<grain_junction_index[grain_junctions[selected_junction]].size(); i++)
                {
                    for (int p=0; p<grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1].size(); p++)
                    {
                        int x=std::max(0,grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1][p].x-posx)/scaling;
                        int y=std::max(0,grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1][p].y-posy)/scaling;
                        canvas(x,y,0,0)=255;
                        canvas(x,y,0,1)=0;
                        canvas(x,y,0,2)=0;
                    }

                    float dx=average_pos[i].x-junctions[grain_junctions[selected_junction]].x;
                    float dy=average_pos[i].y-junctions[grain_junctions[selected_junction]].y;

                    if (dx!=0)
                    {
                        float m=dy/dx;
                        float t=average_pos[i].y-m*average_pos[i].x;

                        if (dx>0)
                        {
                            for (int xx=junctions[grain_junctions[selected_junction]].x; xx<=posx+selectionDisplay.dimx(); xx++)
                            {            
                                int yy=m*(float)xx+t;
                                if (yy>=0 && yy<dim_y && square(junctions[grain_junctions[selected_junction]].x-xx)+
                                    square(junctions[grain_junctions[selected_junction]].y-yy)<square(pixels_average)) 
                                {
                                    canvas((xx-posx)/scaling,(yy-posy)/scaling,0,0)=0;
                                    canvas((xx-posx)/scaling,(yy-posy)/scaling,0,1)=255;
                                    canvas((xx-posx)/scaling,(yy-posy)/scaling,0,2)=0;
                                }
                            }
                        }
                        else
                        {
                            for (int xx=posx; xx<=junctions[grain_junctions[selected_junction]].x; xx++)
                            {            
                                int yy=m*(float)xx+t;
                                if (yy>=0 && yy<dim_y && square(junctions[grain_junctions[selected_junction]].x-xx)+
                                    square(junctions[grain_junctions[selected_junction]].y-yy)<square(pixels_average)) 
                                {
                                    canvas((xx-posx)/scaling,(yy-posy)/scaling,0,0)=0;
                                    canvas((xx-posx)/scaling,(yy-posy)/scaling,0,1)=255;
                                    canvas((xx-posx)/scaling,(yy-posy)/scaling,0,2)=0;
                                }
                            }
                        }
                    }            

                    if (dy!=0)
                    {
                        float m=dx/dy;
                        float t=average_pos[i].x-m*average_pos[i].y;

                        if (dy>0)
                        {
                            for (int yy=junctions[grain_junctions[selected_junction]].y; yy<=posy+selectionDisplay.dimy(); yy++)
                            {            
                                int xx=m*(float)yy+t;
                                if (xx>=0 && xx<dim_x && square(junctions[grain_junctions[selected_junction]].x-xx)+
                                    square(junctions[grain_junctions[selected_junction]].y-yy)<square(pixels_average)) 
                                {
                                    canvas((xx-posx)/scaling,(yy-posy)/scaling,0,0)=0;
                                    canvas((xx-posx)/scaling,(yy-posy)/scaling,0,1)=255;
                                    canvas((xx-posx)/scaling,(yy-posy)/scaling,0,2)=0;
                                }
                            }
                        }
                        else
                        {
                            for (int yy=posy; yy<=junctions[grain_junctions[selected_junction]].y; yy++)
                            {            
                                int xx=m*(float)yy+t;
                                if (xx>=0 && xx<dim_x && square(junctions[grain_junctions[selected_junction]].x-xx)+
                                    square(junctions[grain_junctions[selected_junction]].y-yy)<square(pixels_average)) 
                                {
                                    canvas((xx-posx)/scaling,(yy-posy)/scaling,0,0)=0;
                                    canvas((xx-posx)/scaling,(yy-posy)/scaling,0,1)=255;
                                    canvas((xx-posx)/scaling,(yy-posy)/scaling,0,2)=0;
                                }
                            }
                        }
                    }
                }
            }

            selectionDisplay.resize(canvas.dimx(),canvas.dimy());
        }

        delete boundary_pos;
        delete average_pos;
    }

    else if (mode==3) //grain triple junctions analysis
    {
        //initialise plplot class
        plplot plot = plplot();

        int nr_grain_boundaries=grain_boundary_index.size();
        std::vector< std::vector<point> > grain_boundary_pixels(nr_grain_boundaries);

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

                        //fill grain boundary pixels
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

        for(int j=0; j<grain_junctions.size(); j++)
        {
            for (int i=0; i<grain_junction_index[grain_junctions[j]].size(); i++)
            {
                //check which side of grain boundary fits to junction
                point front=grain_boundary_pixels[grain_junction_index[grain_junctions[j]][i]][0];
                point back=grain_boundary_pixels[grain_junction_index[grain_junctions[j]][i]].back();

                if((fabs(front.x-junctions[grain_junctions[j]].x)+fabs(front.y-junctions[grain_junctions[j]].y))<2)
                {
                    grain_junction_index[grain_junctions[j]][i]=grain_junction_index[grain_junctions[j]][i]+1;
                }
                else if((fabs(back.x-junctions[grain_junctions[j]].x)+fabs(back.y-junctions[grain_junctions[j]].y))<2)
                {
                    grain_junction_index[grain_junctions[j]][i]=-grain_junction_index[grain_junctions[j]][i]-1;
                }
                else
                {
                    std::cout<<"Error: Junction do not fit to boundaries!"<<std::endl;
                    exit(-1);
                }
            }
        }

        std::vector<int> x_values;
        std::vector<float> angle_deviations;

        for (int pixels_average=2; pixels_average<100; pixels_average++)
        {
            std::vector<float> angles;
            float angles_mean;
            float angles_standard_deviation;

            for (int selected_junction=0; selected_junction<grain_junctions.size(); selected_junction++)
            {
                int * boundary_pos = new int[grain_junction_index[grain_junctions[selected_junction]].size()];
                point * average_pos = new point [grain_junction_index[grain_junctions[selected_junction]].size()];

                for (int i=0; i<grain_junction_index[grain_junctions[selected_junction]].size(); i++)
                {
                    boundary_pos[i]=std::min(pixels_average,
                        (int)grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1].size());
                }

                //compute the "center" of pixels
                for (int i=0; i<grain_junction_index[grain_junctions[selected_junction]].size(); i++)
                {
                    int number=0;
                    int sum_x=0;
                    int sum_y=0;

                    for(int p=0; p<boundary_pos[i]; p++)
                    {
                        int index;
                        if(grain_junction_index[grain_junctions[selected_junction]][i]>0) index=p;
                        else index=grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1].size()-1-p;

                        sum_x+=grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1][index].x;
                        sum_y+=grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1][index].y;
                        number++;
                    }

                    average_pos[i].x=sum_x/number;
                    average_pos[i].y=sum_y/number;
                }

                for (int i=0; i<3; i++)
                {
                    float angle_front=atan2((average_pos[i].y-junctions[grain_junctions[selected_junction]].y),
                                            (average_pos[i].x-junctions[grain_junctions[selected_junction]].x));
                    float angle_back=atan2((average_pos[(i+1)%3].y-junctions[grain_junctions[selected_junction]].y),
                                           (average_pos[(i+1)%3].x-junctions[grain_junctions[selected_junction]].x));
                    float angle_other=atan2((average_pos[(i+2)%3].y-junctions[grain_junctions[selected_junction]].y),
                                            (average_pos[(i+2)%3].x-junctions[grain_junctions[selected_junction]].x));

                    //if angle of other arc is in between, add 2pi to smaller angle
                    if(angle_front<angle_other && angle_other<angle_back) angle_front+=2.0f*PI;
                    if(angle_back<angle_other && angle_other<angle_front) angle_back+=2.0f*PI;

                    if (angle_front>angle_back) angles.push_back(180.0f*(angle_front-angle_back)/PI);
                    else angles.push_back(180.0f*(angle_back-angle_front)/PI);
                }

                delete boundary_pos;
                delete average_pos;
            }
            calculate_mean_standard_deviation(angles, angles_mean, angles_standard_deviation);
            x_values.push_back(pixels_average);
            angle_deviations.push_back(angles_standard_deviation);
        }

        std::string filepath_output=filepath_plots;
        filepath_output.append("_deviation_angles");
        filepath_output.append(suffix.c_str());
        filepath_output.append(".svg");

        std::cout<<"Export standard deviation to: "<<filepath_output<<std::endl;
        plot.draw_curv_cross("Pixels to average", "Standard deviation [degree]", "Standard deviation in dihedral angles", x_values, angle_deviations,
                             filepath_output.c_str());
    }

    else if (mode==4) //create grain overview
    {
        cimg_library::CImg<unsigned char> grain_image(dim_x, dim_y, 1, 3);
        std::vector< std::vector<point> > areas(nr_areas);

        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
            {
                //fill areas vector
                point p;
                p.x=x;
                p.y=y;
                areas[region_labels[ws_region_image(x,y)-1]-1].push_back(p);

                grain_image(x,y,0,0)=255;
                grain_image(x,y,0,1)=255;
                grain_image(x,y,0,2)=255;
            }
        }

        //bubble areas to output
        for(int area=0; area<nr_areas; area++)
        {
            if (bubble_area_size[area]>0)//bubble area
            {
                for(size_t k=0; k<areas[area].size(); ++k)
                {
                    int x=areas[area][k].x;
                    int y=areas[area][k].y;
                    grain_image(x,y,0,0)=0;
                    grain_image(x,y,0,1)=0;
                    grain_image(x,y,0,2)=0;
                }
            }
        }

        //grain areas and boundaries to output
        for(size_t area=0;area<nr_areas;area++)
        {
            if(grain_arc_index[area].size()>0)
            {
                for(size_t k=0; k<areas[area].size(); ++k)
                {
                    int x=areas[area][k].x;
                    int y=areas[area][k].y;
                    grain_image(x,y,0,0)=85;
                    grain_image(x,y,0,1)=85;
                    grain_image(x,y,0,2)=85;
                }

                for (int a=0;a<(int)grain_arc_index[area].size();a++)//grain boundaries
                {
                    //now we loop over the points in this arc
                    for(int p=0;p<(int)arcs[grain_arc_index[area][a]].size();p++)
                    {
                        int x=arcs[grain_arc_index[area][a]][p].x;
                        int y=arcs[grain_arc_index[area][a]][p].y;
                        grain_image(x,y,0,0)=0;
                        grain_image(x,y,0,1)=0;
                        grain_image(x,y,0,2)=0;
                    }
                }
            }
        }

        //grain center of mass positions to output image        
        for (int area=0; area<nr_areas; area++)
        {
            if (grain_area_size[area]>0)
            {
                std::ostringstream Str;
                Str << area;
                std::string string(Str.str());

                int xx=grain_area_center_mass[area].x;
                int yy=grain_area_center_mass[area].y;
                if (area<10) grain_image.draw_text(std::max(0,xx-3),std::max(0,yy-10),string.c_str(),color[1],0,1,20,1,1);
                else if (area<100) grain_image.draw_text(std::max(0,xx-10),std::max(0,yy-10),string.c_str(),color[1],0,1,20,1,1);
                else if (area<1000) grain_image.draw_text(std::max(0,xx-17),std::max(0,yy-10),string.c_str(),color[1],0,1,20,1,1);
                else grain_image.draw_text(std::max(0,xx-24),std::max(0,yy-10),string.c_str(),color[1],0,1,20,1,1);
            }
        }

        //export overview
        vigra::FRGBImage output_image(dim_x,dim_y);

        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
            {
                output_image(x,y)[0]=grain_image(x,y,0,0);
                output_image(x,y)[1]=grain_image(x,y,0,1);
                output_image(x,y)[2]=grain_image(x,y,0,2);
            }
        }

        std::string filepath_output=filepath_plots;
        filepath_output.append("_overview");
        filepath_output.append(suffix.c_str());
        filepath_output.append(".bmp");

        std::cout<<"Export grain overview to: "<<filepath_output<<std::endl;
        exportImage(srcImageRange(output_image), vigra::ImageExportInfo(filepath_output.c_str()));
    }

    else if (mode==5) //foldings interactive 
    {
        vigra::FRGBImage orientation_image(dim_x,dim_y);
        std::vector< std::vector<point> > areas(nr_areas);

        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
            {
                //fill areas vector
                point p;
                p.x=x;
                p.y=y;
                areas[region_labels[ws_region_image(x,y)-1]-1].push_back(p);

                orientation_image(x,y)[0]=255;
                orientation_image(x,y)[1]=255;
                orientation_image(x,y)[2]=255;
            }
        }

        float scaling=(float)dim_x/(float)display_x;

        std::vector<int> selected_grains(areas.size(),false);
        int nr_selected_grains=0;

        std::string filepath_selected_grains=filepath_plots;
        filepath_selected_grains.append("_selected_grains");
        filepath_selected_grains.append(suffix.c_str());
        filepath_selected_grains.append(".dat");

        std::ifstream selection_file(filepath_selected_grains.c_str());
        std::ifstream temp_selection_file(filepath_selected_grains.c_str());

        //string is just for testing stuff
        std::string teststring;
        temp_selection_file>>teststring;
        if(selection_file)
        {
            if(teststring.size()!=0)
            {
                std::cout<<"Load selected grains from file: "<<filepath_selected_grains.c_str()<<std::endl;

                int temp;
                int old=-1;
                while(!selection_file.eof())
                {
                    selection_file>>temp;
                    selected_grains[temp]=true;
                    if (old==-1 || old!=temp)
                    {
                        nr_selected_grains++;

                        for(int p=0;p<areas[temp].size();p++)
                        {
                            int xx=areas[temp][p].x;
                            int yy=areas[temp][p].y;

                            image(xx,yy,0,0)=255;
                            image(xx,yy,0,1)=255;
                            image(xx,yy,0,2)=255;
                        }
                    }

                    old=temp;
                }
            }

            selection_file.close();
            temp_selection_file.close();
        }

        for (int arc=0; arc<arcs.size(); arc++)
            if (grain_arc[arc])
                for (int p=0; p<arcs[arc].size(); p++)
                {
                    int x=arcs[arc][p].x;
                    int y=arcs[arc][p].y;

                    for(int a=-scaling/2; a<=scaling/2; a++)
                    {
                        if(x+1-a>0 && x+1-a<dim_x && y>0) image(x-a,y-1,0,0)=255;
                        if(x+1-a>0 && x+1-a<dim_x && y>0) image(x-a,y-1,0,1)=0;
                        if(x+1-a>0 && x+1-a<dim_x && y>0) image(x-a,y-1,0,2)=0;
                        if(x+1-a>0 && x+1-a<dim_x ) image(x-a,y,0,0)=255;
                        if(x+1-a>0 && x+1-a<dim_x ) image(x-a,y,0,1)=0;
                        if(x+1-a>0 && x+1-a<dim_x ) image(x-a,y,0,2)=0;
                        if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) image(x-a,y+1,0,0)=255;
                        if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) image(x-a,y+1,0,1)=0;
                        if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) image(x-a,y+1,0,2)=0;
                    }

                    for(int a=-scaling/2; a<=scaling/2; a++)
                    {
                        if(y+1-a>0 && y+1-a<dim_y && x>0) image(x-1,y-a,0,0)=255;
                        if(y+1-a>0 && y+1-a<dim_y && x>0) image(x-1,y-a,0,1)=0;
                        if(y+1-a>0 && y+1-a<dim_y && x>0) image(x-1,y-a,0,2)=0;
                        if(y+1-a>0 && y+1-a<dim_y ) image(x,y-a,0,0)=255;
                        if(y+1-a>0 && y+1-a<dim_y ) image(x,y-a,0,1)=0;
                        if(y+1-a>0 && y+1-a<dim_y ) image(x,y-a,0,2)=0;
                        if(y+1-a>0 && y+1-a<dim_y && x<dim_x-1) image(x+1,y-a,0,0)=255;
                        if(y+1-a>0 && y+1-a<dim_y && x<dim_x-1) image(x+1,y-a,0,1)=0;
                        if(y+1-a>0 && y+1-a<dim_y && x<dim_x-1) image(x+1,y-a,0,2)=0;
                    }

                }

        for(int area=0; area<nr_areas; area++)
            for (int arc=0; arc<bubble_arc_index[area].size(); arc++)
                for (int p=0; p<arcs[bubble_arc_index[area][arc]].size(); p++)
                {
                    int x=arcs[bubble_arc_index[area][arc]][p].x;
                    int y=arcs[bubble_arc_index[area][arc]][p].y;

                    for(int a=-scaling/2; a<=scaling/2; a++)
                    {
                        if(x+1-a>0 && x+1-a<dim_x && y>0) image(x-a,y-1,0,0)=0;
                        if(x+1-a>0 && x+1-a<dim_x && y>0) image(x-a,y-1,0,1)=0;
                        if(x+1-a>0 && x+1-a<dim_x && y>0) image(x-a,y-1,0,2)=255;
                        if(x+1-a>0 && x+1-a<dim_x ) image(x-a,y,0,0)=0;
                        if(x+1-a>0 && x+1-a<dim_x ) image(x-a,y,0,1)=0;
                        if(x+1-a>0 && x+1-a<dim_x ) image(x-a,y,0,2)=255;
                        if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) image(x-a,y+1,0,0)=0;
                        if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) image(x-a,y+1,0,1)=0;
                        if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) image(x-a,y+1,0,2)=255;
                    }

                    for(int a=-scaling/2; a<=scaling/2; a++)
                    {
                        if(y+1-a>0 && y+1-a<dim_y && x>0) image(x-1,y-a,0,0)=0;
                        if(y+1-a>0 && y+1-a<dim_y && x>0) image(x-1,y-a,0,1)=0;
                        if(y+1-a>0 && y+1-a<dim_y && x>0) image(x-1,y-a,0,2)=255;
                        if(y+1-a>0 && y+1-a<dim_y ) image(x,y-a,0,0)=0;
                        if(y+1-a>0 && y+1-a<dim_y ) image(x,y-a,0,1)=0;
                        if(y+1-a>0 && y+1-a<dim_y ) image(x,y-a,0,2)=255;
                        if(y+1-a>0 && y+1-a<dim_y && x<dim_x-1) image(x+1,y-a,0,0)=0;
                        if(y+1-a>0 && y+1-a<dim_y && x<dim_x-1) image(x+1,y-a,0,1)=0;
                        if(y+1-a>0 && y+1-a<dim_y && x<dim_x-1) image(x+1,y-a,0,2)=255;
                    }
                }

        // canvas for our gui to draw on
        // CImg<type> name(dimx, dimy, dimz, colors)
        cimg_library::CImg<unsigned char> canvas(display_x, display_y, 1, 3);

        // draw the current selection on canvas
        cimg_forXY(canvas,x,y)
        {
            if( scaling*x + posx >= 0 && scaling*x + posx <  image.dimx() && scaling*y + posy >= 0 && scaling*y + posy < image.dimy() )
            {
                canvas(x,y,0,0) = image(scaling*x + posx,scaling*y + posy,0,0);
                canvas(x,y,0,1) = image(scaling*x + posx,scaling*y + posy,0,1);
                canvas(x,y,0,2) = image(scaling*x + posx,scaling*y + posy,0,2);
            }
        }

        bool saved = true;
        bool selected = false;

        // Initialization of menu
        cimg_library::CImg<unsigned char> menu(250,40,1,3,0);
        cimg_library::CImg<unsigned char> men=menu;

        cimg_library::CImgDisplay  main_menu(menu,"Menu");
        cimg_library::CImgDisplay  selectionDisplay(canvas,"Select grains for orientation plot!");

        main_menu.move(0,0);
        selectionDisplay.move(0,main_menu.dimy()+57);

        while(!selected)
        {
            //MENU STUFF
            men=menu;

            if(saved) men.draw_text(2,20,"Saved",color[0],0,1,11,1,1);
            else men.draw_text(2,20,"Unsaved",color[1],0,1,11,1,1);

            {   
                std::ostringstream Str;
                Str << nr_selected_grains;
                std::string temp_string(Str.str());
                std::string name="Selected grains: ";
                name.append(temp_string);
                men.draw_text(2,10,name.c_str(),color[3],0,1,11,1,1);
            }

            main_menu.display(men);

            selectionDisplay.display(canvas);
            selectionDisplay.wait();
            bool changed=false;

            //this disables to change the window sizes
            if (main_menu.is_resized)
            {
                main_menu.resize(main_menu);
            }

            // move selection
            if( selectionDisplay.is_keyARROWUP )
            {
                posy -= 0.2*scaling*canvas.dimy();
                if (posy<0) posy=0;
                changed=true;
            }
            if( selectionDisplay.is_keyARROWDOWN )
            {
                posy += 0.2*scaling*canvas.dimy();
                if (posy>(image.dimy()-scaling*canvas.dimy())) posy=image.dimy()-scaling*canvas.dimy();
                changed=true;
            }

                //TESTS IF S-KEY (SAVE) IS PRESSED
                if(selectionDisplay.key==115||selectionDisplay.key==83)
                {
                //save selected grains
                std::ofstream selection_file_out(filepath_selected_grains.c_str());

                for (int selected_grain=0; selected_grain<selected_grains.size(); selected_grain++)
                {
                    if(selected_grains[selected_grain]) selection_file_out <<selected_grain<< "\n";
                }

                std::cout<<"Save selected grains to file: "<<filepath_selected_grains.c_str()<<std::endl;
                selection_file_out.close();

                    selected=false;
                saved=true;
            }

                //TESTS IF Q-KEY (QUIT) IS PRESSED OR WINDOWS IS CLOSED
                if(selectionDisplay.key==113||selectionDisplay.key==81 || selectionDisplay.is_closed)
                {
                if(selectionDisplay.key==113||selectionDisplay.key==81) selectionDisplay.close();
                if (!main_menu.is_closed) main_menu.close();
                    selected=true;
                }

            //TEST IF MOUSE IS CLICKED
            if(selectionDisplay.button)
            {
                //Test if the mouse button is clicked on the image area 
                 if (selectionDisplay.mouse_y>=0 && selectionDisplay.mouse_x>=0)
                {
                    int x = scaling*selectionDisplay.mouse_x + posx;
                    int y = scaling*selectionDisplay.mouse_y + posy;

                    if(!selected_grains[region_labels[ws_region_image(x,y)-1]-1] && grain_area_size[region_labels[ws_region_image(x,y)-1]-1]>0)
                    {
                        selected_grains[region_labels[ws_region_image(x,y)-1]-1]=true;
                        nr_selected_grains++;

                        for(int p=0;p<areas[region_labels[ws_region_image(x,y)-1]-1].size();p++)
                        {
                            int xx=areas[region_labels[ws_region_image(x,y)-1]-1][p].x;
                            int yy=areas[region_labels[ws_region_image(x,y)-1]-1][p].y;

                            image(xx,yy,0,0)=255;
                            image(xx,yy,0,1)=255;
                            image(xx,yy,0,2)=255;
                        }

                        changed=true;
                        saved=false;
                    }
                    else if(selected_grains[region_labels[ws_region_image(x,y)-1]-1] && grain_area_size[region_labels[ws_region_image(x,y)-1]-1]>0)
                    {
                        selected_grains[region_labels[ws_region_image(x,y)-1]-1]=false;
                        nr_selected_grains--;

                        for(int p=0;p<areas[region_labels[ws_region_image(x,y)-1]-1].size();p++)
                        {
                            int xx=areas[region_labels[ws_region_image(x,y)-1]-1][p].x;
                            int yy=areas[region_labels[ws_region_image(x,y)-1]-1][p].y;

                            image(xx,yy,0,0)=original_image(xx,yy,0,0);
                            image(xx,yy,0,1)=original_image(xx,yy,0,1);
                            image(xx,yy,0,2)=original_image(xx,yy,0,2);
                        }

                        changed=true;
                        saved=false;
                    }

                    if(changed)
                    {
                        for (int arc=0; arc<grain_arc_index[region_labels[ws_region_image(x,y)-1]-1].size(); arc++)
                            for (int p=0; p<arcs[grain_arc_index[region_labels[ws_region_image(x,y)-1]-1][arc]].size(); p++)
                            {
                                int xx=arcs[grain_arc_index[region_labels[ws_region_image(x,y)-1]-1][arc]][p].x;
                                int yy=arcs[grain_arc_index[region_labels[ws_region_image(x,y)-1]-1][arc]][p].y;

                                if(grain_arc[grain_arc_index[region_labels[ws_region_image(x,y)-1]-1][arc]])
                                {
                                    for(int a=-scaling/2; a<=scaling/2; a++)
                                    {
                                        if(xx+1-a>0 && xx+1-a<dim_x && yy>0) image(xx-a,yy-1,0,0)=255;
                                        if(xx+1-a>0 && xx+1-a<dim_x && yy>0) image(xx-a,yy-1,0,1)=0;
                                        if(xx+1-a>0 && xx+1-a<dim_x && yy>0) image(xx-a,yy-1,0,2)=0;
                                        if(xx+1-a>0 && xx+1-a<dim_x ) image(xx-a,yy,0,0)=255;
                                        if(xx+1-a>0 && xx+1-a<dim_x ) image(xx-a,yy,0,1)=0;
                                        if(xx+1-a>0 && xx+1-a<dim_x ) image(xx-a,yy,0,2)=0;
                                        if(xx+1-a>0 && xx+1-a<dim_x && yy<dim_y-1) image(xx-a,yy+1,0,0)=255;
                                        if(xx+1-a>0 && xx+1-a<dim_x && yy<dim_y-1) image(xx-a,yy+1,0,1)=0;
                                        if(xx+1-a>0 && xx+1-a<dim_x && yy<dim_y-1) image(xx-a,yy+1,0,2)=0;
                                    }

                                    for(int a=-scaling/2; a<=scaling/2; a++)
                                    {
                                        if(yy+1-a>0 && yy+1-a<dim_y && xx>0) image(xx-1,yy-a,0,0)=255;
                                        if(yy+1-a>0 && yy+1-a<dim_y && xx>0) image(xx-1,yy-a,0,1)=0;
                                        if(yy+1-a>0 && yy+1-a<dim_y && xx>0) image(xx-1,yy-a,0,2)=0;
                                        if(yy+1-a>0 && yy+1-a<dim_y ) image(xx,yy-a,0,0)=255;
                                        if(yy+1-a>0 && yy+1-a<dim_y ) image(xx,yy-a,0,1)=0;
                                        if(yy+1-a>0 && yy+1-a<dim_y ) image(xx,yy-a,0,2)=0;
                                        if(yy+1-a>0 && yy+1-a<dim_y && xx<dim_x-1) image(xx+1,yy-a,0,0)=255;
                                        if(yy+1-a>0 && yy+1-a<dim_y && xx<dim_x-1) image(xx+1,yy-a,0,1)=0;
                                        if(yy+1-a>0 && yy+1-a<dim_y && xx<dim_x-1) image(xx+1,yy-a,0,2)=0;
                                    }
                                }
                                else
                                {
                                    for(int a=-scaling/2; a<=scaling/2; a++)
                                    {
                                        if(xx+1-a>0 && xx+1-a<dim_x && yy>0) image(xx-a,yy-1,0,0)=0;
                                        if(xx+1-a>0 && xx+1-a<dim_x && yy>0) image(xx-a,yy-1,0,1)=0;
                                        if(xx+1-a>0 && xx+1-a<dim_x && yy>0) image(xx-a,yy-1,0,2)=255;
                                        if(xx+1-a>0 && xx+1-a<dim_x ) image(xx-a,yy,0,0)=0;
                                        if(xx+1-a>0 && xx+1-a<dim_x ) image(xx-a,yy,0,1)=0;
                                        if(xx+1-a>0 && xx+1-a<dim_x ) image(xx-a,yy,0,2)=255;
                                        if(xx+1-a>0 && xx+1-a<dim_x && yy<dim_y-1) image(xx-a,yy+1,0,0)=0;
                                        if(xx+1-a>0 && xx+1-a<dim_x && yy<dim_y-1) image(xx-a,yy+1,0,1)=0;
                                        if(xx+1-a>0 && xx+1-a<dim_x && yy<dim_y-1) image(xx-a,yy+1,0,2)=255;
                                    }

                                    for(int a=-scaling/2; a<=scaling/2; a++)
                                    {
                                        if(yy+1-a>0 && yy+1-a<dim_y && xx>0) image(xx-1,yy-a,0,0)=0;
                                        if(yy+1-a>0 && yy+1-a<dim_y && xx>0) image(xx-1,yy-a,0,1)=0;
                                        if(yy+1-a>0 && yy+1-a<dim_y && xx>0) image(xx-1,yy-a,0,2)=255;
                                        if(yy+1-a>0 && yy+1-a<dim_y ) image(xx,yy-a,0,0)=0;
                                        if(yy+1-a>0 && yy+1-a<dim_y ) image(xx,yy-a,0,1)=0;
                                        if(yy+1-a>0 && yy+1-a<dim_y ) image(xx,yy-a,0,2)=255;
                                        if(yy+1-a>0 && yy+1-a<dim_y && xx<dim_x-1) image(xx+1,yy-a,0,0)=0;
                                        if(yy+1-a>0 && yy+1-a<dim_y && xx<dim_x-1) image(xx+1,yy-a,0,1)=0;
                                        if(yy+1-a>0 && yy+1-a<dim_y && xx<dim_x-1) image(xx+1,yy-a,0,2)=255;
                                    }
                                }
                            }
                    }
                    }
                }

            // press "ENTER" if you are content with your selection
            if( selectionDisplay.is_keyENTER)
            {
                selectionDisplay.close();
                if (!main_menu.is_closed) main_menu.close();

                if (!saved)
                {
                    //save selected grains
                    std::ofstream selection_file_out(filepath_selected_grains.c_str());

                    for (int selected_grain=0; selected_grain<selected_grains.size(); selected_grain++)
                    {
                        if(selected_grains[selected_grain]) selection_file_out <<selected_grain<< "\n";
                    }

                    std::cout<<"Save selected grains to file: "<<filepath_selected_grains.c_str()<<std::endl;
                    selection_file_out.close();
                }

                for (int selected_grain=0; selected_grain<selected_grains.size(); selected_grain++)
                {
                    //grain boundary pixels
                    std::vector<point> all_boundary_pixels;

                    for (int arc=0; arc<grain_arc_index[selected_grain].size(); arc++)
                        for (int p=0; p<arcs[grain_arc_index[selected_grain][arc]].size(); p++)
                        {
                            int x=arcs[grain_arc_index[selected_grain][arc]][p].x;
                            int y=arcs[grain_arc_index[selected_grain][arc]][p].y;

                            orientation_image(x,y)[0]=0;
                            orientation_image(x,y)[1]=0;
                            orientation_image(x,y)[2]=0;

                            all_boundary_pixels.push_back(arcs[grain_arc_index[selected_grain][arc]][p]);
                        }

                        //if(ellipse_flattening[selected_grain]>2.0f)
                        if(selected_grains[selected_grain])
                        {
                            area_range range;
                            range.x_low=dim_x;
                            range.y_low=dim_y;
                            range.x_high=0;
                            range.y_high=0;

                            //rotate grain boundary pixels
                            for(int p=0; p<all_boundary_pixels.size(); p++)
                            {
                                point pixel;
                                pixel.x=all_boundary_pixels[p].x-grain_area_center_mass[selected_grain].x;
                                pixel.y=all_boundary_pixels[p].y-grain_area_center_mass[selected_grain].y;

                                all_boundary_pixels[p].x=grain_area_center_mass[selected_grain].x+
                                    pixel.x*cos(PI*ellipse_long_axis_angle[selected_grain]/180.0f)-
                                    pixel.y*sin(PI*ellipse_long_axis_angle[selected_grain]/180.0f);
                                all_boundary_pixels[p].y=grain_area_center_mass[selected_grain].y+
                                    pixel.x*sin(PI*ellipse_long_axis_angle[selected_grain]/180.0f)+
                                    pixel.y*cos(PI*ellipse_long_axis_angle[selected_grain]/180.0f);

                                //find ranges
                                if(all_boundary_pixels[p].x<range.x_low) range.x_low=all_boundary_pixels[p].x;
                                if(all_boundary_pixels[p].y<range.y_low) range.y_low=all_boundary_pixels[p].y;
                                if(all_boundary_pixels[p].x>range.x_high) range.x_high=all_boundary_pixels[p].x;
                                if(all_boundary_pixels[p].y>range.y_high) range.y_high=all_boundary_pixels[p].y;
                            }

                            //show lines
                            for(int x=range.x_low; x<=range.x_high; x++)
                            {
                                for(int y=range.y_low; y<=range.y_high; y+=20)
                                {
                                    point pixel;
                                    pixel.x=x-grain_area_center_mass[selected_grain].x;
                                    pixel.y=y-grain_area_center_mass[selected_grain].y;

                                    int xx=grain_area_center_mass[selected_grain].x+
                                        pixel.x*cos(PI*ellipse_long_axis_angle[selected_grain]/180.0f)+
                                        pixel.y*sin(PI*ellipse_long_axis_angle[selected_grain]/180.0f);
                                    int yy=grain_area_center_mass[selected_grain].y-
                                        pixel.x*sin(PI*ellipse_long_axis_angle[selected_grain]/180.0f)+
                                        pixel.y*cos(PI*ellipse_long_axis_angle[selected_grain]/180.0f);

                                    if (region_labels[ws_region_image(xx,yy)-1]==selected_grain+1)
                                    {
                                        orientation_image(xx,yy)[0]=255;
                                        orientation_image(xx,yy)[1]=0;
                                        orientation_image(xx,yy)[2]=0;
                                    }
                                }
                            }
                    }
                }

                //export orientation
                std::string filepath_output=filepath_plots;
                filepath_output.append("_orientation");
                filepath_output.append(suffix.c_str());
                filepath_output.append(".bmp");

                std::cout<<"Export grain orientation to: "<<filepath_output<<std::endl;
                exportImage(srcImageRange(orientation_image), vigra::ImageExportInfo(filepath_output.c_str()));

                selected = true;
            }

            // resize selection
            if( selectionDisplay.is_resized )
            {
                selectionDisplay.resize();
                scaling=(float)dim_x/(float)selectionDisplay.dimx();

                image=original_image;

                for(int area=0; area<nr_areas; area++)
                    if(selected_grains[area])
                        for(int p=0;p<areas[area].size();p++)
                        {
                            int xx=areas[area][p].x;
                            int yy=areas[area][p].y;

                            image(xx,yy,0,0)=255;
                            image(xx,yy,0,1)=255;
                            image(xx,yy,0,2)=255;
                        }


                for (int arc=0; arc<arcs.size(); arc++)
                    if (grain_arc[arc])
                        for (int p=0; p<arcs[arc].size(); p++)
                        {
                            int x=arcs[arc][p].x;
                            int y=arcs[arc][p].y;

                            for(int a=-scaling/2; a<=scaling/2; a++)
                            {
                                if(x+1-a>0 && x+1-a<dim_x && y>0) image(x-a,y-1,0,0)=255;
                                if(x+1-a>0 && x+1-a<dim_x && y>0) image(x-a,y-1,0,1)=0;
                                if(x+1-a>0 && x+1-a<dim_x && y>0) image(x-a,y-1,0,2)=0;
                                if(x+1-a>0 && x+1-a<dim_x ) image(x-a,y,0,0)=255;
                                if(x+1-a>0 && x+1-a<dim_x ) image(x-a,y,0,1)=0;
                                if(x+1-a>0 && x+1-a<dim_x ) image(x-a,y,0,2)=0;
                                if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) image(x-a,y+1,0,0)=255;
                                if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) image(x-a,y+1,0,1)=0;
                                if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) image(x-a,y+1,0,2)=0;
                            }

                            for(int a=-scaling/2; a<=scaling/2; a++)
                            {
                                if(y+1-a>0 && y+1-a<dim_y && x>0) image(x-1,y-a,0,0)=255;
                                if(y+1-a>0 && y+1-a<dim_y && x>0) image(x-1,y-a,0,1)=0;
                                if(y+1-a>0 && y+1-a<dim_y && x>0) image(x-1,y-a,0,2)=0;
                                if(y+1-a>0 && y+1-a<dim_y ) image(x,y-a,0,0)=255;
                                if(y+1-a>0 && y+1-a<dim_y ) image(x,y-a,0,1)=0;
                                if(y+1-a>0 && y+1-a<dim_y ) image(x,y-a,0,2)=0;
                                if(y+1-a>0 && y+1-a<dim_y && x<dim_x-1) image(x+1,y-a,0,0)=255;
                                if(y+1-a>0 && y+1-a<dim_y && x<dim_x-1) image(x+1,y-a,0,1)=0;
                                if(y+1-a>0 && y+1-a<dim_y && x<dim_x-1) image(x+1,y-a,0,2)=0;
                            }

                        }

                for(int area=0; area<nr_areas; area++)
                    for (int arc=0; arc<bubble_arc_index[area].size(); arc++)
                        for (int p=0; p<arcs[bubble_arc_index[area][arc]].size(); p++)
                        {
                            int x=arcs[bubble_arc_index[area][arc]][p].x;
                            int y=arcs[bubble_arc_index[area][arc]][p].y;

                            for(int a=-scaling/2; a<=scaling/2; a++)
                            {
                                if(x+1-a>0 && x+1-a<dim_x && y>0) image(x-a,y-1,0,0)=0;
                                if(x+1-a>0 && x+1-a<dim_x && y>0) image(x-a,y-1,0,1)=0;
                                if(x+1-a>0 && x+1-a<dim_x && y>0) image(x-a,y-1,0,2)=255;
                                if(x+1-a>0 && x+1-a<dim_x ) image(x-a,y,0,0)=0;
                                if(x+1-a>0 && x+1-a<dim_x ) image(x-a,y,0,1)=0;
                                if(x+1-a>0 && x+1-a<dim_x ) image(x-a,y,0,2)=255;
                                if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) image(x-a,y+1,0,0)=0;
                                if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) image(x-a,y+1,0,1)=0;
                                if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) image(x-a,y+1,0,2)=255;
                            }

                            for(int a=-scaling/2; a<=scaling/2; a++)
                            {
                                if(y+1-a>0 && y+1-a<dim_y && x>0) image(x-1,y-a,0,0)=0;
                                if(y+1-a>0 && y+1-a<dim_y && x>0) image(x-1,y-a,0,1)=0;
                                if(y+1-a>0 && y+1-a<dim_y && x>0) image(x-1,y-a,0,2)=255;
                                if(y+1-a>0 && y+1-a<dim_y ) image(x,y-a,0,0)=0;
                                if(y+1-a>0 && y+1-a<dim_y ) image(x,y-a,0,1)=0;
                                if(y+1-a>0 && y+1-a<dim_y ) image(x,y-a,0,2)=255;
                                if(y+1-a>0 && y+1-a<dim_y && x<dim_x-1) image(x+1,y-a,0,0)=0;
                                if(y+1-a>0 && y+1-a<dim_y && x<dim_x-1) image(x+1,y-a,0,1)=0;
                                if(y+1-a>0 && y+1-a<dim_y && x<dim_x-1) image(x+1,y-a,0,2)=255;
                            }
                        }

                changed=true;
            }

            if(changed)
            {
                canvas.assign(selectionDisplay.dimx(), selectionDisplay.dimy(), 1, 3);
                cimg_forXY(canvas,x,y)
                {
                    if( scaling*x + posx >= 0 && scaling*x + posx <  image.dimx() && scaling*y + posy >= 0 && scaling*y + posy < image.dimy() )
                    {
                        canvas(x,y,0,0) = image(scaling*x + posx,scaling*y + posy,0,0);
                        canvas(x,y,0,1) = image(scaling*x + posx,scaling*y + posy,0,1);
                        canvas(x,y,0,2) = image(scaling*x + posx,scaling*y + posy,0,2);
                    }
                }
            }
        }
    }

    else if (mode==6) //create local curvature map
    {
        int nr_grain_boundaries=grain_boundary_index.size();
        std::vector< std::vector<point> > grain_boundary_pixels(nr_grain_boundaries);

        vigra::IRGBImage curvature_image(dim_x,dim_y);
        std::vector< std::vector<point> > areas(nr_areas);

        float curvature_min=0.005f;
        float dislocation_density_bin_width=20.0f;

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

        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
            {
                //fill areas vector
                point po;
                po.x=x;
                po.y=y;
                areas[region_labels[ws_region_image(x,y)-1]-1].push_back(po);

                curvature_image(x,y)[0]=255;
                curvature_image(x,y)[1]=255;
                curvature_image(x,y)[2]=255;
            }
        }

        for (int arc=0; arc<grain_boundary_pixels.size(); arc++)
        {
            for (int p=0; p<grain_boundary_curvs[arc].size(); p++)
            {
                float curv=fabs(grain_boundary_curvs[arc][p]);

                if (curv>curvature_min)
                {
                    int bin=(int)(dislocation_density_bin_width*log10(curv/curvature_min));

                    int red=(255*std::min(20,bin))/20;
                    int green=255-(255*std::min(20,bin))/20;
                    int blue=255-(255*std::min(20,bin))/20;

                    int x = grain_boundary_pixels[arc][p].x;
                    int y = grain_boundary_pixels[arc][p].y;

                    for(int a=-5; a<6; a++)
                    {
                        if(x+1-a>0 && x+1-a<dim_x && y>0) curvature_image(x-a,y-1)[0]=red;
                        if(x+1-a>0 && x+1-a<dim_x && y>0) curvature_image(x-a,y-1)[1]=green;
                        if(x+1-a>0 && x+1-a<dim_x && y>0) curvature_image(x-a,y-1)[2]=blue;
                        if(x+1-a>0 && x+1-a<dim_x ) curvature_image(x-a,y)[0]=red;
                        if(x+1-a>0 && x+1-a<dim_x ) curvature_image(x-a,y)[1]=green;
                        if(x+1-a>0 && x+1-a<dim_x ) curvature_image(x-a,y)[2]=blue;
                        if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) curvature_image(x-a,y+1)[0]=red;
                        if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) curvature_image(x-a,y+1)[1]=green;
                        if(x+1-a>0 && x+1-a<dim_x && y<dim_y-1) curvature_image(x-a,y+1)[2]=blue;
                    }
                }
                else if (p<grain_boundary_pixels[arc].size())//special case with grain_boundary_pixels[].size()=1
                {
                    curvature_image(grain_boundary_pixels[arc][p].x,grain_boundary_pixels[arc][p].y)[0]=0;
                    curvature_image(grain_boundary_pixels[arc][p].x,grain_boundary_pixels[arc][p].y)[1]=0;
                    curvature_image(grain_boundary_pixels[arc][p].x,grain_boundary_pixels[arc][p].y)[2]=0;
                }
            }
        }

        //bubble areas to output
        for(int area=0; area<nr_areas; area++)
        {
            if (bubble_area_size[area]>0)//bubble area
            {
                for(size_t k=0; k<areas[area].size(); ++k)
                {
                    int x=areas[area][k].x;
                    int y=areas[area][k].y;
                    curvature_image(x,y)[0]=0;
                    curvature_image(x,y)[1]=0;
                    curvature_image(x,y)[2]=0;
                }
            }
        }

        //export curvature map
        std::string filepath_output=filepath_plots;
        filepath_output.append(".dislocation_density_map");
        filepath_output.append(suffix.c_str());
        filepath_output.append(".jpg");

        std::cout<<"Export curvature map to: "<<filepath_output<<std::endl;
        exportImage(srcImageRange(curvature_image), vigra::ImageExportInfo(filepath_output.c_str()));
    }

    else if  (mode==9) //create grain boundary overview
    {   
        cimg_library::CImg<unsigned char> boundary_image(dim_x, dim_y, 1, 3);
        std::vector< std::vector<point> > areas(nr_areas);

        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
            {
                //fill areas vector
                point p;
                p.x=x;
                p.y=y;
                areas[region_labels[ws_region_image(x,y)-1]-1].push_back(p);

                boundary_image(x,y,0,0)=255;
                boundary_image(x,y,0,1)=255;
                boundary_image(x,y,0,2)=255;
            }
        }

        //bubble areas to output
        for(int area=0; area<nr_areas; area++)
        {
            if (bubble_area_size[area]>0)//bubble area
            {
                for(size_t k=0; k<areas[area].size(); ++k)
                {
                    int x=areas[area][k].x;
                    int y=areas[area][k].y;
                    boundary_image(x,y,0,0)=0;
                    boundary_image(x,y,0,1)=0;
                    boundary_image(x,y,0,2)=0;
                }
            }
        }

        //grain areas and boundaries to output
        for(size_t area=0;area<nr_areas;area++)
        {
            if(grain_arc_index[area].size()>0)
            {
                for(size_t k=0; k<areas[area].size(); ++k)
                {
                    int x=areas[area][k].x;
                    int y=areas[area][k].y;
                    boundary_image(x,y,0,0)=85;
                    boundary_image(x,y,0,1)=85;
                    boundary_image(x,y,0,2)=85;
                }

                for (int a=0;a<(int)grain_arc_index[area].size();a++)//grain boundaries
                {
                    //now we loop over the points in this arc
                    for(int p=0;p<(int)arcs[grain_arc_index[area][a]].size();p++)
                    {
                        int x=arcs[grain_arc_index[area][a]][p].x;
                        int y=arcs[grain_arc_index[area][a]][p].y;
                        boundary_image(x,y,0,0)=0;
                        boundary_image(x,y,0,1)=0;
                        boundary_image(x,y,0,2)=0;
                    }
                }
            }
        }

        for(int boundary=0; boundary<grain_boundary_index.size(); boundary++)
        {
            std::vector<point> grain_boundary_pixels;

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

                        grain_boundary_pixels.push_back(po);
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

                        grain_boundary_pixels.push_back(po);
                    }
                }
            }

            std::ostringstream Str;
            Str << boundary;
            std::string string(Str.str());

            int xx = grain_boundary_pixels[grain_boundary_pixels.size()/2].x;
            int yy = grain_boundary_pixels[grain_boundary_pixels.size()/2].y;
            if (boundary<10) boundary_image.draw_text(std::max(0,xx-3),std::max(0,yy-10),string.c_str(),color[1],0,1,20,1,1);
            else if (boundary<100) boundary_image.draw_text(std::max(0,xx-10),std::max(0,yy-10),string.c_str(),color[1],0,1,20,1,1);
            else if (boundary<1000) boundary_image.draw_text(std::max(0,xx-17),std::max(0,yy-10),string.c_str(),color[1],0,1,20,1,1);
            else boundary_image.draw_text(std::max(0,xx-24),std::max(0,yy-10),string.c_str(),color[1],0,1,20,1,1);
        }

        //export overview
        vigra::FRGBImage output_image(dim_x,dim_y);

        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
            {
                output_image(x,y)[0]=boundary_image(x,y,0,0);
                output_image(x,y)[1]=boundary_image(x,y,0,1);
                output_image(x,y)[2]=boundary_image(x,y,0,2);
            }
        }

        std::string filepath_output=filepath_plots;
        filepath_output.append("_overview2");
        filepath_output.append(suffix.c_str());
        filepath_output.append(".bmp");

        std::cout<<"Export grain overview to: "<<filepath_output<<std::endl;
        exportImage(srcImageRange(output_image), vigra::ImageExportInfo(filepath_output.c_str()));
    }
}

void grain_menu(cimg_library::CImg<unsigned char> & men, cimg_library::CImg<unsigned char> & output_men, int selected_grain, bool only_inner_grains,
    std::vector<bool> inner_grain_areas, std::string grain_string, int unmarked, color_type * color, int diff_x, int diff_y, float scaling,
    std::vector< std::vector<point> > arcs, int dim_x, int dim_y, long * grain_area_size, std::vector< std::vector<int> >  grain_arc_index,
    std::vector<point> grain_area_center_mass, std::vector<int> grain_areas, std::vector<float> grain_roundness, std::vector<float> grain_box_flattening,
    std::vector<float> grain_box_width, std::vector<float> grain_box_height, std::vector< std::vector<float> > ellipse_params, std::vector<float> ellipse_long_axis,
    std::vector<float> ellipse_flattening, std::vector<float> ellipse_long_axis_angle, std::vector<float> grain_arc_number, std::vector<float> grain_neighbors,
    std::vector<float> grain_longest_arc_length)
{
    if (only_inner_grains) men.draw_text(2,0,"Only inner grains - Change selection or print!",color[1],0,1,11,1,1);
    else men.draw_text(2,0,"Change selection or print!",color[1],0,1,11,1,1);
    if (unmarked==0) men.draw_text(330,10, "Boundaries are shown",color[1],0,1,11,1,1);
    else if (unmarked==1) men.draw_text(330,10, "Boundaries and boxes are shown",color[1],0,1,11,1,1);
    else if (unmarked==2) men.draw_text(330,10, "Boundaries and ellipses are shown",color[1],0,1,11,1,1);
    else if (unmarked==3) men.draw_text(330,10, "Boundaries, boxes and ellipses are shown",color[1],0,1,11,1,1);
    else if (unmarked==4) men.draw_text(330,10, "Boundaries and fixed boxes are shown",color[1],0,1,11,1,1);
    else if (unmarked==5) men.draw_text(330,10, "Boundaries are not shown",color[1],0,1,11,1,1);

    std::string name="Selected grain: ";
    name.append(grain_string);
    men.draw_text(2,10,name.c_str(),color[3],0,1,11,1,1);

    {   
        std::ostringstream Str;
        Str << scaling;
        std::string temp_string(Str.str());
        std::string name="Scaling: ";
        name.append(temp_string);
        men.draw_text(164,10,name.c_str(),color[3],0,1,11,1,1);
    }

    {   
        std::ostringstream Str;
        Str << std::fixed << std::setprecision(2) << grain_area_size[grain_areas[selected_grain]]/area_scaling << " mm^2";
        std::string temp_string(Str.str());
        std::string name="Size: ";
        name.append(temp_string);
        men.draw_text(2,20,name.c_str(),color[0],0,1,11,1,1);
        output_men.draw_text(2,0,name.c_str(),color[3],0,1,11,1,1);
    }

    {   
        std::ostringstream Str;
        Str << std::fixed << std::setprecision(2) << diff_x/length_scaling << " mm";
        std::string temp_string(Str.str());
        std::string name="Width: ";
        name.append(temp_string);
        men.draw_text(2,30,name.c_str(),color[3],0,1,11,1,1);
        output_men.draw_text(2,10,name.c_str(),color[3],0,1,11,1,1);
    }

    {   
        std::ostringstream Str;
        Str << std::fixed << std::setprecision(2) << diff_y/length_scaling << " mm";
        std::string temp_string(Str.str());
        std::string name="Height: ";
        name.append(temp_string);
        men.draw_text(164,30,name.c_str(),color[3],0,1,11,1,1);
        output_men.draw_text(164,10,name.c_str(),color[3],0,1,11,1,1);
    }

    {   
        std::ostringstream Str;
        Str << std::fixed << std::setprecision(2) << grain_roundness[grain_areas[selected_grain]];
        std::string temp_string(Str.str());
        std::string name="Roundness factor: ";
        name.append(temp_string);
        men.draw_text(330,30,name.c_str(),color[3],0,1,11,1,1);
        output_men.draw_text(330,10,name.c_str(),color[3],0,1,11,1,1);
    }

    if(fabs(1.0f-ellipse_params[grain_areas[selected_grain]][5])<0.1f)//ellipse fitting gave result
    {
        if (ellipse_long_axis[grain_areas[selected_grain]]>0.0f)//axes found correctly
        {   
            std::ostringstream Str;
            Str << std::fixed << std::setprecision(2) << ellipse_flattening[grain_areas[selected_grain]];
            std::string temp_string(Str.str());
            std::string name="Ellipse flattening: ";
            name.append(temp_string);
            men.draw_text(500,40,name.c_str(),color[0],0,1,11,1,1);
            output_men.draw_text(500,20,name.c_str(),color[3],0,1,11,1,1);
        }

        if(ellipse_params[grain_areas[selected_grain]][0]!=ellipse_params[grain_areas[selected_grain]][2])//ellipse is not degenerated
        {   
            std::ostringstream Str;
            Str << std::fixed << std::setprecision(1) << ellipse_long_axis_angle[grain_areas[selected_grain]];
            std::string temp_string(Str.str());
            std::string name="Long axis angle: ";
            name.append(temp_string);
            men.draw_text(500,30,name.c_str(),color[3],0,1,11,1,1);
            output_men.draw_text(500,10,name.c_str(),color[3],0,1,11,1,1);

            area_range range2;
            range2.x_low=dim_x;
            range2.y_low=dim_y;
            range2.x_high=0;
            range2.y_high=0;

            for (int arc=0; arc<grain_arc_index[grain_areas[selected_grain]].size(); arc++)
                for (int p=0; p<arcs[grain_arc_index[grain_areas[selected_grain]][arc]].size(); p++)
                {
                    point pixel;
                    pixel.x=arcs[grain_arc_index[grain_areas[selected_grain]][arc]][p].x-grain_area_center_mass[grain_areas[selected_grain]].x;
                    pixel.y=arcs[grain_arc_index[grain_areas[selected_grain]][arc]][p].y-grain_area_center_mass[grain_areas[selected_grain]].y;

                    int x=grain_area_center_mass[grain_areas[selected_grain]].x+
                        pixel.x*cos(PI*ellipse_long_axis_angle[grain_areas[selected_grain]]/180.0f)-
                        pixel.y*sin(PI*ellipse_long_axis_angle[grain_areas[selected_grain]]/180.0f);
                    int y=grain_area_center_mass[grain_areas[selected_grain]].y+
                        pixel.x*sin(PI*ellipse_long_axis_angle[grain_areas[selected_grain]]/180.0f)+
                        pixel.y*cos(PI*ellipse_long_axis_angle[grain_areas[selected_grain]]/180.0f);

                    //find ranges
                    if(x<range2.x_low) range2.x_low=x;
                    if(y<range2.y_low) range2.y_low=y;
                    if(x>range2.x_high) range2.x_high=x;
                    if(y>range2.y_high) range2.y_high=y;
                }

            {   
                std::ostringstream Str;
                Str << std::fixed << std::setprecision(2) << grain_box_width[grain_areas[selected_grain]] << " mm";
                std::string temp_string(Str.str());
                std::string name="Box width: ";
                name.append(temp_string);
                men.draw_text(2,40,name.c_str(),color[0],0,1,11,1,1);
                output_men.draw_text(2,20,name.c_str(),color[3],0,1,11,1,1);
            }

            {   
                std::ostringstream Str;
                Str << std::fixed << std::setprecision(2) << grain_box_height[grain_areas[selected_grain]] << " mm";
                std::string temp_string(Str.str());
                std::string name="Box height: ";
                name.append(temp_string);
                men.draw_text(164,40,name.c_str(),color[0],0,1,11,1,1);
                output_men.draw_text(164,20,name.c_str(),color[3],0,1,11,1,1);
            }

            {   
                std::ostringstream Str;
                Str << std::fixed << std::setprecision(2) << grain_box_flattening[grain_areas[selected_grain]];
                std::string temp_string(Str.str());
                std::string name="Box flattening: ";
                name.append(temp_string);
                men.draw_text(330,40,name.c_str(),color[0],0,1,11,1,1);
                output_men.draw_text(330,20,name.c_str(),color[3],0,1,11,1,1);
            }
        }
    }

    if (inner_grain_areas[selected_grain])
    {   
        {   
            std::ostringstream Str;
            Str << std::fixed << std::setprecision(0) << grain_arc_number[grain_areas[selected_grain]];
            std::string temp_string(Str.str());
            std::string name="Number boundaries: ";
            name.append(temp_string);
            men.draw_text(330,20,name.c_str(),color[0],0,1,11,1,1);
            output_men.draw_text(330,0,name.c_str(),color[3],0,1,11,1,1);
        }

        {   
            std::ostringstream Str;
            Str << std::fixed << std::setprecision(0) << grain_neighbors[grain_areas[selected_grain]];
            std::string temp_string(Str.str());
            std::string name="Number neighbors: ";
            name.append(temp_string);
            men.draw_text(164,20,name.c_str(),color[0],0,1,11,1,1);
            output_men.draw_text(164,0,name.c_str(),color[3],0,1,11,1,1);
        }

        {   
            std::ostringstream Str;
            Str << std::fixed << std::setprecision(2) << grain_longest_arc_length[grain_areas[selected_grain]] << " mm";
            std::string temp_string(Str.str());
            std::string name="Longest boundary: ";
            name.append(temp_string);
            men.draw_text(500,20,name.c_str(),color[0],0,1,11,1,1);
            output_men.draw_text(500,0,name.c_str(),color[3],0,1,11,1,1);
        }
    }
}

void show_ellipse_box(cimg_library::CImg<unsigned char> & image, int selected_grain, int unmarked, float angle, area_range & range, FitEllipse fitEllipse,
    std::vector<point> & ellipse_points, std::vector< std::vector<point> > arcs, int dim_x, int dim_y, std::vector< std::vector<int> > grain_arc_index,
    std::vector<point> grain_area_center_mass, std::vector<bool> grain_arc, std::vector<int> grain_areas, std::vector< std::vector<float> > ellipse_params)
{
    //grain boundary pixels
    std::vector<point> all_boundary_pixels;

    for (int arc=0; arc<grain_arc_index[grain_areas[selected_grain]].size(); arc++)
        for (int p=0; p<arcs[grain_arc_index[grain_areas[selected_grain]][arc]].size(); p++)
        {
            int x=arcs[grain_arc_index[grain_areas[selected_grain]][arc]][p].x;
            int y=arcs[grain_arc_index[grain_areas[selected_grain]][arc]][p].y;
            if (grain_arc[grain_arc_index[grain_areas[selected_grain]][arc]])
            {
                image(x,y,0,0)=255;
                image(x,y,0,1)=0;
                image(x,y,0,2)=0;
            }
            else
            {
                image(x,y,0,0)=0;
                image(x,y,0,1)=0;
                image(x,y,0,2)=255;
            }

            all_boundary_pixels.push_back(arcs[grain_arc_index[grain_areas[selected_grain]][arc]][p]);
        }

    if (unmarked==1 || unmarked==3 || unmarked==4)
    {
        //rotate grain boundary pixels
        for(int p=0; p<all_boundary_pixels.size(); p++)
        {
            point pixel;
            pixel.x=all_boundary_pixels[p].x-grain_area_center_mass[grain_areas[selected_grain]].x;
            pixel.y=all_boundary_pixels[p].y-grain_area_center_mass[grain_areas[selected_grain]].y;

            all_boundary_pixels[p].x=grain_area_center_mass[grain_areas[selected_grain]].x+
                pixel.x*cos(PI*angle/180.0f)-
                pixel.y*sin(PI*angle/180.0f);
            all_boundary_pixels[p].y=grain_area_center_mass[grain_areas[selected_grain]].y+
                pixel.x*sin(PI*angle/180.0f)+
                pixel.y*cos(PI*angle/180.0f);

            //find ranges
            if(all_boundary_pixels[p].x<range.x_low) range.x_low=all_boundary_pixels[p].x;
            if(all_boundary_pixels[p].y<range.y_low) range.y_low=all_boundary_pixels[p].y;
            if(all_boundary_pixels[p].x>range.x_high) range.x_high=all_boundary_pixels[p].x;
            if(all_boundary_pixels[p].y>range.y_high) range.y_high=all_boundary_pixels[p].y;
        }

        //show box
        for(int x=range.x_low; x<=range.x_high; x++)
        {
            point pixel;
            pixel.x=x-grain_area_center_mass[grain_areas[selected_grain]].x;
            pixel.y=range.y_low-grain_area_center_mass[grain_areas[selected_grain]].y;

            int xx=grain_area_center_mass[grain_areas[selected_grain]].x+
                pixel.x*cos(PI*angle/180.0f)+
                pixel.y*sin(PI*angle/180.0f);
            int yy=grain_area_center_mass[grain_areas[selected_grain]].y-
                pixel.x*sin(PI*angle/180.0f)+
                pixel.y*cos(PI*angle/180.0f);

            if (xx>=0 && xx<dim_x && yy>=0 && yy<dim_y)
            {
                image(xx,yy,0,0)=0;
                image(xx,yy,0,1)=255;
                image(xx,yy,0,2)=0;
            }

            pixel.x=x-grain_area_center_mass[grain_areas[selected_grain]].x;
            pixel.y=range.y_high-grain_area_center_mass[grain_areas[selected_grain]].y;

            xx=grain_area_center_mass[grain_areas[selected_grain]].x+
                pixel.x*cos(PI*angle/180.0f)+
                pixel.y*sin(PI*angle/180.0f);
            yy=grain_area_center_mass[grain_areas[selected_grain]].y-
                pixel.x*sin(PI*angle/180.0f)+
                pixel.y*cos(PI*angle/180.0f);

            if (xx>=0 && xx<dim_x && yy>=0 && yy<dim_y)
            {
                image(xx,yy,0,0)=0;
                image(xx,yy,0,1)=255;
                image(xx,yy,0,2)=0;
            }
        }

        for(int y=range.y_low; y<=range.y_high; y++)
        {
            point pixel;
            pixel.x=range.x_low-grain_area_center_mass[grain_areas[selected_grain]].x;
            pixel.y=y-grain_area_center_mass[grain_areas[selected_grain]].y;

            int xx=grain_area_center_mass[grain_areas[selected_grain]].x+
                pixel.x*cos(PI*angle/180.0f)+
                pixel.y*sin(PI*angle/180.0f);
            int yy=grain_area_center_mass[grain_areas[selected_grain]].y-
                pixel.x*sin(PI*angle/180.0f)+
                pixel.y*cos(PI*angle/180.0f);

            if (xx>=0 && xx<dim_x && yy>=0 && yy<dim_y)
            {
                image(xx,yy,0,0)=0;
                image(xx,yy,0,1)=255;
                image(xx,yy,0,2)=0;
            }

            pixel.x=range.x_high-grain_area_center_mass[grain_areas[selected_grain]].x;
            pixel.y=y-grain_area_center_mass[grain_areas[selected_grain]].y;

            xx=grain_area_center_mass[grain_areas[selected_grain]].x+
                pixel.x*cos(PI*angle/180.0f)+
                pixel.y*sin(PI*angle/180.0f);
            yy=grain_area_center_mass[grain_areas[selected_grain]].y-
                pixel.x*sin(PI*angle/180.0f)+
                pixel.y*cos(PI*angle/180.0f);

            if (xx>=0 && xx<dim_x && yy>=0 && yy<dim_y)
            {
                image(xx,yy,0,0)=0;
                image(xx,yy,0,1)=255;
                image(xx,yy,0,2)=0;
            }
        }
    }

    if ((unmarked==2 || unmarked==3) && fabs(1.0f-ellipse_params[grain_areas[selected_grain]][5])<0.1f)
    {
        fitEllipse.drawConic(ellipse_params[grain_areas[selected_grain]],10*all_boundary_pixels.size(),ellipse_points);

        for (int i=0; i<ellipse_points.size(); i++)
        {
            if (ellipse_points[i].x>=0 && ellipse_points[i].x<dim_x && ellipse_points[i].y>=0 && ellipse_points[i].y<dim_y)
            {
                image(ellipse_points[i].x,ellipse_points[i].y,0,0)=0;
                image(ellipse_points[i].x,ellipse_points[i].y,0,1)=255;
                image(ellipse_points[i].x,ellipse_points[i].y,0,2)=0;
            }
        }
    }
}

void print_grain(float scaling, int dim_x, int dim_y, cimg_library::CImg<unsigned char> output_men, cimg_library::CImg<unsigned char> image, int posx, int posy,
    std::string filepath_plots, int unmarked, std::string grain_string, std::string suffix)
{
    vigra::BRGBImage output(std::max((float)output_men.dimx(),scaling*dim_x),scaling*dim_y+output_men.dimy());
    vigra::BRGBImage output2(std::max((float)output_men.dimx(),scaling*dim_x),scaling*dim_y+output_men.dimy());

    for (int x=scaling*dim_x; x<output.width(); x++)
        for (int y=0; y<output.height()-output_men.dimy(); y++)
        {  
            output(x,y)[0]=0;
            output(x,y)[1]=0;
            output(x,y)[2]=0;
        }

    for (int x=0; x<scaling*dim_x; x++)
        for (int y=0; y<output.height()-output_men.dimy(); y++)
        {  
            output(x,y)[0]=image(x+posx,y+posy,0,0);
            output(x,y)[1]=image(x+posx,y+posy,0,1);
            output(x,y)[2]=image(x+posx,y+posy,0,2);

            output2(x,y)[0]=image(x+posx,y+posy,0,0);
            output2(x,y)[1]=image(x+posx,y+posy,0,1);
            output2(x,y)[2]=image(x+posx,y+posy,0,2);
        }

    for (int x=0; x<output_men.dimx(); x++)
        for (int y=0; y<output_men.dimy(); y++)
        {  
            output(x,y+output.height()-output_men.dimy())[0]=output_men(x,y,0,0);
            output(x,y+output.height()-output_men.dimy())[1]=output_men(x,y,0,1);
            output(x,y+output.height()-output_men.dimy())[2]=output_men(x,y,0,2);

            output2(x,y+output.height()-output_men.dimy())[0]=output_men(x,y,0,0);
            output2(x,y+output.height()-output_men.dimy())[1]=output_men(x,y,0,1);
            output2(x,y+output.height()-output_men.dimy())[2]=output_men(x,y,0,2);
        }

    for (int x=0; x<scaling*dim_x; x++)
        for (int y=0; y<output.height()-output_men.dimy(); y++)
            for (int i=0; i<3; i++)
                if (output2(x,y)[i]==255 && unmarked!=5)
                {  
                    if (x>0) output(x-1,y)[i]=255;
                    if (x>1) output(x-2,y)[i]=255;
                    if (x+1<scaling*dim_x) output(x+1,y)[i]=255;
                    if (x+2<scaling*dim_x) output(x+2,y)[i]=255;

                    if (y>0)
                    {  
                        output(x,y-1)[i]=255;
                        if (x>0) output(x-1,y-1)[i]=255;
                        if (x>1) output(x-2,y-1)[i]=255;
                        if (x+1<scaling*dim_x) output(x+1,y-1)[i]=255;
                        if (x+2<scaling*dim_x) output(x+2,y-1)[i]=255;
                    }

                    if (y>1)
                    {  
                        output(x,y-2)[i]=255;
                        if (x>0) output(x-1,y-2)[i]=255;
                        if (x>1) output(x-2,y-2)[i]=255;
                        if (x+1<scaling*dim_x) output(x+1,y-2)[i]=255;
                        if (x+2<scaling*dim_x) output(x+2,y-2)[i]=255;
                    }

                    if (y+1<output.height()-output_men.dimy())
                    {  
                        output(x,y+1)[i]=255;
                        if (x>0) output(x-1,y+1)[i]=255;
                        if (x>1) output(x-2,y+1)[i]=255;
                        if (x+1<scaling*dim_x) output(x+1,y+1)[i]=255;
                        if (x+2<scaling*dim_x) output(x+2,y+1)[i]=255;
                    }

                    if (y+2<output.height()-output_men.dimy())
                    {  
                        output(x,y+2)[i]=255;
                        if (x>0) output(x-1,y+2)[i]=255;
                        if (x>1) output(x-2,y+2)[i]=255;
                        if (x+1<scaling*dim_x) output(x+1,y+2)[i]=255;
                        if (x+2<scaling*dim_x) output(x+2,y+2)[i]=255;
                    }
                }

    for (int x=0; x<scaling*dim_x; x++)
        for (int y=0; y<output.height()-output_men.dimy(); y++)
            if (output(x,y)[0]==255 && output(x,y)[1]==255 && output(x,y)[2]<255)
            {  
                output(x,y)[0]=0;
                output(x,y)[2]=0;
            }

    std::string filepath_output=filepath_plots;

    filepath_output.append("_grain");
    filepath_output.append(grain_string);
    if (unmarked==1) filepath_output.append("_box");
    if (unmarked==2) filepath_output.append("_ellipse");
    if (unmarked==3) filepath_output.append("_box_ellipse");
    if (unmarked==4) filepath_output.append("_fixed_box");
    if (unmarked==5) filepath_output.append("_unmarked");

    filepath_output.append(suffix.c_str());
    filepath_output.append(".bmp");

    std::cout<<"Export grain to: "<<filepath_output<<std::endl;
    exportImage(srcImageRange(output), vigra::ImageExportInfo(filepath_output.c_str()));
}

void boundary_menu(cimg_library::CImg<unsigned char> & men, cimg_library::CImg<unsigned char> & output_men, int selected_boundary, bool only_inner_boundaries,
    std::vector<bool> inner_boundary, std::string boundary_string, int unmarked, color_type * color, float scaling,
    std::vector< std::vector<int> > grain_boundary_index, std::vector< std::vector<point> > grain_boundary_pixels, std::vector<float> turning_point)
{
    if (only_inner_boundaries) men.draw_text(2,0,"Only inner boundaries - Change selection or print!",color[1],0,1,11,1,1);
    else men.draw_text(2,0,"Change selection or print!",color[1],0,1,11,1,1);
    if (unmarked) men.draw_text(250,10, "Boundaries are not shown",color[1],0,1,11,1,1);
    else men.draw_text(250,10, "Boundaries are shown",color[1],0,1,11,1,1);

    std::string name="Selected boundary: ";
    name.append(boundary_string);
    men.draw_text(2,10,name.c_str(),color[3],0,1,11,1,1);

    {
        std::ostringstream Str;
        Str << scaling;
        std::string temp_string(Str.str());
        std::string name="Scaling: ";
        name.append(temp_string);
        men.draw_text(150,10,name.c_str(),color[0],0,1,11,1,1);
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
        output_men.draw_text(2,0,name.c_str(),color[3],0,1,11,1,1);
    }

    if (inner_boundary[selected_boundary])
    {   
        std::ostringstream Str;
        Str << std::fixed << std::setprecision(2) << grain_boundary_pixels[selected_boundary].size()/length_scaling << " mm";
        std::string temp_string(Str.str());
        std::string name="Length: ";
        name.append(temp_string);
        men.draw_text(2,30,name.c_str(),color[3],0,1,11,1,1);
        output_men.draw_text(2,10,name.c_str(),color[3],0,1,11,1,1);
    }

    {   
        std::ostringstream Str;
        Str << std::fixed << std::setprecision(0) << turning_point[selected_boundary];
        std::string temp_string(Str.str());
        std::string name="Number turning points: ";
        name.append(temp_string);
        men.draw_text(150,30,name.c_str(),color[0],0,1,11,1,1);
        output_men.draw_text(150,10,name.c_str(),color[3],0,1,11,1,1);
    }
}

void print_boundary(float scaling, int dim_x, int dim_y, cimg_library::CImg<unsigned char> output_men, cimg_library::CImg<unsigned char> image, int posx, int posy,
    std::string filepath_plots, int unmarked, std::string boundary_string, int selected_boundary, std::vector< std::vector<float> > grain_boundary_curvs, plplot plot,
    std::string suffix)
{
    vigra::BRGBImage output(std::max((float)output_men.dimx(),scaling*dim_x),scaling*dim_y+output_men.dimy());
    vigra::BRGBImage output2(std::max((float)output_men.dimx(),scaling*dim_x),scaling*dim_y+output_men.dimy());

    for (int x=scaling*dim_x; x<output.width(); x++)
        for (int y=0; y<output.height()-output_men.dimy(); y++)
        {  
            output(x,y)[0]=0;
            output(x,y)[1]=0;
            output(x,y)[2]=0;
        }

    for (int x=0; x<scaling*dim_x; x++)
        for (int y=0; y<output.height()-output_men.dimy(); y++)
        {  
            output(x,y)[0]=image(x+posx,y+posy,0,0);
            output(x,y)[1]=image(x+posx,y+posy,0,1);
            output(x,y)[2]=image(x+posx,y+posy,0,2);

            output2(x,y)[0]=image(x+posx,y+posy,0,0);
            output2(x,y)[1]=image(x+posx,y+posy,0,1);
            output2(x,y)[2]=image(x+posx,y+posy,0,2);
        }

    for (int x=0; x<output_men.dimx(); x++)
        for (int y=0; y<output_men.dimy(); y++)
        {  
            output(x,y+output.height()-output_men.dimy())[0]=output_men(x,y,0,0);
            output(x,y+output.height()-output_men.dimy())[1]=output_men(x,y,0,1);
            output(x,y+output.height()-output_men.dimy())[2]=output_men(x,y,0,2);

            output2(x,y+output.height()-output_men.dimy())[0]=output_men(x,y,0,0);
            output2(x,y+output.height()-output_men.dimy())[1]=output_men(x,y,0,1);
            output2(x,y+output.height()-output_men.dimy())[2]=output_men(x,y,0,2);
        }

    for (int x=0; x<scaling*dim_x; x++)
        for (int y=0; y<output.height()-output_men.dimy(); y++)
            for (int i=0; i<3; i++)
                if (output2(x,y)[i]==255)
                {  
                    if (x>0) output(x-1,y)[i]=255;
                    if (x>1) output(x-2,y)[i]=255;
                    if (x+1<scaling*dim_x) output(x+1,y)[i]=255;
                    if (x+2<scaling*dim_x) output(x+2,y)[i]=255;

                    if (y>0)
                    {  
                        output(x,y-1)[i]=255;
                        if (x>0) output(x-1,y-1)[i]=255;
                        if (x>1) output(x-2,y-1)[i]=255;
                        if (x+1<scaling*dim_x) output(x+1,y-1)[i]=255;
                        if (x+2<scaling*dim_x) output(x+2,y-1)[i]=255;
                    }

                    if (y>1)
                    {  
                        output(x,y-2)[i]=255;
                        if (x>0) output(x-1,y-2)[i]=255;
                        if (x>1) output(x-2,y-2)[i]=255;
                        if (x+1<scaling*dim_x) output(x+1,y-2)[i]=255;
                        if (x+2<scaling*dim_x) output(x+2,y-2)[i]=255;
                    }

                    if (y+1<output.height()-output_men.dimy())
                    {  
                        output(x,y+1)[i]=255;
                        if (x>0) output(x-1,y+1)[i]=255;
                        if (x>1) output(x-2,y+1)[i]=255;
                        if (x+1<scaling*dim_x) output(x+1,y+1)[i]=255;
                        if (x+2<scaling*dim_x) output(x+2,y+1)[i]=255;
                    }

                    if (y+2<output.height()-output_men.dimy())
                    {  
                        output(x,y+2)[i]=255;
                        if (x>0) output(x-1,y+2)[i]=255;
                        if (x>1) output(x-2,y+2)[i]=255;
                        if (x+1<scaling*dim_x) output(x+1,y+2)[i]=255;
                        if (x+2<scaling*dim_x) output(x+2,y+2)[i]=255;
                    }
                }

    for (int x=0; x<scaling*dim_x; x++)
        for (int y=0; y<output.height()-output_men.dimy(); y++)
        {  
            if (output(x,y)[0]<255 && output(x,y)[1]<255 && output(x,y)[2]==255)
            {  
                output(x,y)[0]=0;
                output(x,y)[1]=0;
            }

            if (output2(x,y)[0]==255 && output2(x,y)[1]<255 && output2(x,y)[2]<255)// make short boundaries better visible
            {  
                for (int i=1; i<3; i++)
                {  
                    output(x,y)[i]=0;
                    if (x>0) output(x-1,y)[i]=0;
                    if (x>1) output(x-2,y)[i]=0;
                    if (x+1<scaling*dim_x) output(x+1,y)[i]=0;
                    if (x+2<scaling*dim_x) output(x+2,y)[i]=0;

                    if (y>0)
                    {  
                        output(x,y-1)[i]=0;
                        if (x>0) output(x-1,y-1)[i]=0;
                        if (x>1) output(x-2,y-1)[i]=0;
                        if (x+1<scaling*dim_x) output(x+1,y-1)[i]=0;
                        if (x+2<scaling*dim_x) output(x+2,y-1)[i]=0;
                    }

                    if (y>1)
                    {  
                        output(x,y-2)[i]=0;
                        if (x>0) output(x-1,y-2)[i]=0;
                        if (x>1) output(x-2,y-2)[i]=0;
                        if (x+1<scaling*dim_x) output(x+1,y-2)[i]=0;
                        if (x+2<scaling*dim_x) output(x+2,y-2)[i]=0;
                    }

                    if (y+1<output.height()-output_men.dimy())
                    {  
                        output(x,y+1)[i]=0;
                        if (x>0) output(x-1,y+1)[i]=0;
                        if (x>1) output(x-2,y+1)[i]=0;
                        if (x+1<scaling*dim_x) output(x+1,y+1)[i]=0;
                        if (x+2<scaling*dim_x) output(x+2,y+1)[i]=0;
                    }

                    if (y+2<output.height()-output_men.dimy())
                    {  
                        output(x,y+2)[i]=0;
                        if (x>0) output(x-1,y+2)[i]=0;
                        if (x>1) output(x-2,y+2)[i]=0;
                        if (x+1<scaling*dim_x) output(x+1,y+2)[i]=0;
                        if (x+2<scaling*dim_x) output(x+2,y+2)[i]=0;
                    }
                }
            }
        }

    std::string filepath_output=filepath_plots;

    filepath_output.append("_boundary");
    filepath_output.append(boundary_string);
    if (unmarked) filepath_output.append("_unmarked");
    filepath_output.append(suffix.c_str());
    filepath_output.append(".bmp");

    std::cout<<"Export grain boundary to: "<<filepath_output<<std::endl;
    exportImage(srcImageRange(output), vigra::ImageExportInfo(filepath_output.c_str()));

    if (unmarked) filepath_output.resize(filepath_output.size()-9);

    filepath_output.resize(filepath_output.size()-4-suffix.size());
    filepath_output.append("_curv");
    filepath_output.append(suffix.c_str());
    filepath_output.append(".svg");

    std::vector<int> x_values;
    for (int i=0; i<grain_boundary_curvs[selected_boundary].size(); i++) x_values.push_back(i);

    std::cout<<"Export curvature to: "<<filepath_output<<std::endl;
    plot.draw_curv_cross("Position [pixel]", "Curvature [rad/pixel]", "Curvature", x_values, grain_boundary_curvs[selected_boundary],
                         filepath_output.c_str());
}

void junction_menu(cimg_library::CImg<unsigned char> & men, cimg_library::CImg<unsigned char> & output_men, int selected_junction, std::string junction_string,
    int unmarked, color_type * color, float scaling, std::vector<int> grain_junctions, std::vector<point> junctions, std::vector< std::vector<int> > grain_junction_index,
    int pixels_average, point * average_pos)
{
    men.draw_text(2,0,"Change selection or print!",color[1],0,1,11,1,1);
    if (unmarked) men.draw_text(350,10, "Boundaries are not shown",color[1],0,1,11,1,1);
    else men.draw_text(300,10, "Boundaries are shown",color[1],0,1,11,1,1);

    std::string name="Selected grain junction: ";
    name.append(junction_string);
    men.draw_text(2,10,name.c_str(),color[3],0,1,11,1,1);

    {   
        std::ostringstream Str;
        Str << scaling;
        std::string temp_string(Str.str());
        std::string name="Scaling: ";
        name.append(temp_string);
        men.draw_text(200,10,name.c_str(),color[0],0,1,11,1,1);
    }

    {   
        std::ostringstream Str;

        for (int i=0; i<grain_junction_index[grain_junctions[selected_junction]].size(); i++)
        {
            Str << fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1 << " ";
        }

        std::string temp_string(Str.str());
        std::string name="Grain boundaries: ";
        name.append(temp_string);
        men.draw_text(2,20,name.c_str(),color[0],0,1,11,1,1);
        output_men.draw_text(2,0,name.c_str(),color[3],0,1,11,1,1);
    }

    {   
        std::ostringstream Str;
        Str << pixels_average;
        std::string temp_string(Str.str());
        std::string name="Pixels to average: ";
        name.append(temp_string);
        men.draw_text(2,30,name.c_str(),color[3],0,1,11,1,1);
        output_men.draw_text(2,10,name.c_str(),color[3],0,1,11,1,1);
    }

    {   
        std::ostringstream Str;

        for (int i=0; i<3; i++)
        {
            float angle_front=atan2((average_pos[i].y-junctions[grain_junctions[selected_junction]].y),
                                    (average_pos[i].x-junctions[grain_junctions[selected_junction]].x));
            float angle_back=atan2((average_pos[(i+1)%3].y-junctions[grain_junctions[selected_junction]].y),
                                   (average_pos[(i+1)%3].x-junctions[grain_junctions[selected_junction]].x));
            float angle_other=atan2((average_pos[(i+2)%3].y-junctions[grain_junctions[selected_junction]].y),
                                    (average_pos[(i+2)%3].x-junctions[grain_junctions[selected_junction]].x));

            //if angle of other arc is in between, add 2pi to smaller angle
            if(angle_front<angle_other && angle_other<angle_back) angle_front+=2.0f*PI;
            if(angle_back<angle_other && angle_other<angle_front) angle_back+=2.0f*PI;

            if (angle_front>angle_back) Str << std::fixed << std::setprecision(1) << 180.0f*(angle_front-angle_back)/PI << " ";
            else Str << std::fixed << std::setprecision(1) << 180.0f*(angle_back-angle_front)/PI << " ";
        }

        std::string temp_string(Str.str());
        std::string name="Angles (degree): ";
        name.append(temp_string);
        men.draw_text(200,30,name.c_str(),color[0],0,1,11,1,1);
        output_men.draw_text(200,10,name.c_str(),color[3],0,1,11,1,1);
    }
}

void show_orientation(cimg_library::CImg<unsigned char> & image, std::vector<int> grain_junctions, std::vector<point> junctions, int dim_x, int dim_y,
    std::vector< std::vector<int> > grain_junction_index, std::vector< std::vector<point> > grain_boundary_pixels, int selected_junction, int posx, int posy,
    int pixels_average, point * average_pos, int dimx, int dimy)
{
    for (int i=0; i<grain_junction_index[grain_junctions[selected_junction]].size(); i++)
    {
        for (int p=0; p<grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1].size(); p++)
        {
            int x=grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1][p].x;
            int y=grain_boundary_pixels[fabs(grain_junction_index[grain_junctions[selected_junction]][i])-1][p].y;
            image(x,y,0,0)=255;
            image(x,y,0,1)=0;
            image(x,y,0,2)=0;
        }

        float dx=average_pos[i].x-junctions[grain_junctions[selected_junction]].x;
        float dy=average_pos[i].y-junctions[grain_junctions[selected_junction]].y;

        if (dx!=0)
        {
            float m=dy/dx;
            float t=average_pos[i].y-m*average_pos[i].x;

            if (dx>0)
            {
                for (int xx=junctions[grain_junctions[selected_junction]].x; xx<=posx+dimx; xx++)
                {            
                    int yy=m*(float)xx+t;
                    if (yy>=0 && yy<dim_y && square(junctions[grain_junctions[selected_junction]].x-xx)+
                        square(junctions[grain_junctions[selected_junction]].y-yy)<square(pixels_average)) 
                    {
                        image(xx,yy,0,0)=0;
                        image(xx,yy,0,1)=255;
                        image(xx,yy,0,2)=0;
                    }
                }
            }
            else
            {
                for (int xx=posx; xx<=junctions[grain_junctions[selected_junction]].x; xx++)
                {            
                    int yy=m*(float)xx+t;
                    if (yy>=0 && yy<dim_y && square(junctions[grain_junctions[selected_junction]].x-xx)+
                        square(junctions[grain_junctions[selected_junction]].y-yy)<square(pixels_average)) 
                    {
                        image(xx,yy,0,0)=0;
                        image(xx,yy,0,1)=255;
                        image(xx,yy,0,2)=0;
                    }
                }
            }
        }            

        if (dy!=0)
        {
            float m=dx/dy;
            float t=average_pos[i].x-m*average_pos[i].y;

            if (dy>0)
            {
                for (int yy=junctions[grain_junctions[selected_junction]].y; yy<=posy+dimy; yy++)
                {            
                    int xx=m*(float)yy+t;
                    if (xx>=0 && xx<dim_x && square(junctions[grain_junctions[selected_junction]].x-xx)+
                        square(junctions[grain_junctions[selected_junction]].y-yy)<square(pixels_average)) 
                    {
                        image(xx,yy,0,0)=0;
                        image(xx,yy,0,1)=255;
                        image(xx,yy,0,2)=0;
                    }
                }
            }
            else
            {
                for (int yy=posy; yy<=junctions[grain_junctions[selected_junction]].y; yy++)
                {            
                    int xx=m*(float)yy+t;
                    if (xx>=0 && xx<dim_x && square(junctions[grain_junctions[selected_junction]].x-xx)+
                        square(junctions[grain_junctions[selected_junction]].y-yy)<square(pixels_average)) 
                    {
                        image(xx,yy,0,0)=0;
                        image(xx,yy,0,1)=255;
                        image(xx,yy,0,2)=0;
                    }
                }
            }
        }
    }
}

void print_junction(float scaling, int dim_x, int dim_y, cimg_library::CImg<unsigned char> output_men, cimg_library::CImg<unsigned char> image, int posx, int posy,
    std::string filepath_plots, int unmarked, std::string junction_string, std::string suffix)
{
    vigra::BRGBImage output(std::max((float)output_men.dimx(),scaling*dim_x),scaling*dim_y+output_men.dimy());
    vigra::BRGBImage output2(std::max((float)output_men.dimx(),scaling*dim_x),scaling*dim_y+output_men.dimy());

    for (int x=scaling*dim_x; x<output.width(); x++)
        for (int y=0; y<output.height()-output_men.dimy(); y++)
        {  
            output(x,y)[0]=0;
            output(x,y)[1]=0;
            output(x,y)[2]=0;
        }

    for (int x=0; x<scaling*dim_x; x++)
        for (int y=0; y<output.height()-output_men.dimy(); y++)
        {  
            output(x,y)[0]=image(x+posx,y+posy,0,0);
            output(x,y)[1]=image(x+posx,y+posy,0,1);
            output(x,y)[2]=image(x+posx,y+posy,0,2);

            output2(x,y)[0]=image(x+posx,y+posy,0,0);
            output2(x,y)[1]=image(x+posx,y+posy,0,1);
            output2(x,y)[2]=image(x+posx,y+posy,0,2);
        }

    for (int x=0; x<output_men.dimx(); x++)
        for (int y=0; y<output_men.dimy(); y++)
        {  
            output(x,y+output.height()-output_men.dimy())[0]=output_men(x,y,0,0);
            output(x,y+output.height()-output_men.dimy())[1]=output_men(x,y,0,1);
            output(x,y+output.height()-output_men.dimy())[2]=output_men(x,y,0,2);

            output2(x,y+output.height()-output_men.dimy())[0]=output_men(x,y,0,0);
            output2(x,y+output.height()-output_men.dimy())[1]=output_men(x,y,0,1);
            output2(x,y+output.height()-output_men.dimy())[2]=output_men(x,y,0,2);
        }

    for (int x=0; x<scaling*dim_x; x++)
        for (int y=0; y<output.height()-output_men.dimy(); y++)
            for (int i=0; i<3; i++)
                if (output2(x,y)[i]==255)
                {  
                    if (x>0) output(x-1,y)[i]=255;
                    if (x>1) output(x-2,y)[i]=255;
                    if (x+1<scaling*dim_x) output(x+1,y)[i]=255;
                    if (x+2<scaling*dim_x) output(x+2,y)[i]=255;

                    if (y>0)
                    {  
                        output(x,y-1)[i]=255;
                        if (x>0) output(x-1,y-1)[i]=255;
                        if (x>1) output(x-2,y-1)[i]=255;
                        if (x+1<scaling*dim_x) output(x+1,y-1)[i]=255;
                        if (x+2<scaling*dim_x) output(x+2,y-1)[i]=255;
                    }

                    if (y>1)
                    {  
                        output(x,y-2)[i]=255;
                        if (x>0) output(x-1,y-2)[i]=255;
                        if (x>1) output(x-2,y-2)[i]=255;
                        if (x+1<scaling*dim_x) output(x+1,y-2)[i]=255;
                        if (x+2<scaling*dim_x) output(x+2,y-2)[i]=255;
                    }

                    if (y+1<output.height()-output_men.dimy())
                    {  
                        output(x,y+1)[i]=255;
                        if (x>0) output(x-1,y+1)[i]=255;
                        if (x>1) output(x-2,y+1)[i]=255;
                        if (x+1<scaling*dim_x) output(x+1,y+1)[i]=255;
                        if (x+2<scaling*dim_x) output(x+2,y+1)[i]=255;
                    }

                    if (y+2<output.height()-output_men.dimy())
                    {  
                        output(x,y+2)[i]=255;
                        if (x>0) output(x-1,y+2)[i]=255;
                        if (x>1) output(x-2,y+2)[i]=255;
                        if (x+1<scaling*dim_x) output(x+1,y+2)[i]=255;
                        if (x+2<scaling*dim_x) output(x+2,y+2)[i]=255;
                    }
                }

    for (int x=0; x<scaling*dim_x; x++)
        for (int y=0; y<output.height()-output_men.dimy(); y++)
            if (output(x,y)[0]==255 && output(x,y)[1]==255 && output(x,y)[2]<255)
            {  
                output(x,y)[0]=0;
                output(x,y)[2]=0;
            }

    std::string filepath_output=filepath_plots;

    filepath_output.append("_junction");
    filepath_output.append(junction_string);
    if (unmarked) filepath_output.append("_unmarked");
    filepath_output.append(suffix.c_str());
    filepath_output.append(".bmp");

    std::cout<<"Export grain junction to: "<<filepath_output<<std::endl;
    exportImage(srcImageRange(output), vigra::ImageExportInfo(filepath_output.c_str()));
}
