/*! \file correct_prediction.h
 *  \brief Used for network correction.
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
#include "CImg.h"

typedef unsigned char color_type[3];
/*! \fn select_image_section(cimg_library::CImg<unsigned char> image, int display_x, int display_y, int & posx, int & posy)
  \brief Selects an image section.
  \param image the test image
  \param display_x width
  \param display_y height
  \param posx X-Position
  \param posy Y-Position
 */
cimg_library::CImg<unsigned char> select_image_section(cimg_library::CImg<unsigned char> image, int display_x, int display_y, int & posx, int & posy)
{
    // canvas for our gui to draw on
    // CImg<type> name(dimx, dimy, dimz, colors)
    cimg_library::CImg<unsigned char> canvas(display_x, display_y, 1, 3);

    // draw the current selection on canvas
    cimg_forXY(canvas,x,y)
    {
        if( x + posx >= 0 && x + posx <  image.dimx() && y + posy >= 0 && y + posy < image.dimy() )
        {
            canvas(x,y,0,0) = image(x + posx,y + posy,0,0);
            canvas(x,y,0,1) = image(x + posx,y + posy,0,1);
            canvas(x,y,0,2) = image(x + posx,y + posy,0,2);
        }
    }

    // show the current canvas
    cimg_library::CImgDisplay selectionDisplay(canvas,"Select image section for correction!");

    bool selected = false;

    // this loop allows to change the relative position and size of the image part which we are
    // selecting
    while(!selected)
    {
        canvas.assign( selectionDisplay.dimx(), selectionDisplay.dimy(), 1, 3);
        cimg_forXY(canvas,x,y)
        {
            if( x + posx >= 0 && x + posx <  image.dimx() && y + posy >= 0 && y + posy < image.dimy() )
            {
                canvas(x,y,0,0) = image(x + posx,y + posy,0,0);
                canvas(x,y,0,1) = image(x + posx,y + posy,0,1);
                canvas(x,y,0,2) = image(x + posx,y + posy,0,2);
            }
        }

        canvas.display(selectionDisplay);

        selectionDisplay.wait();

        // move selection
        if( selectionDisplay.is_keyARROWUP )
        {
            posy -= 100;
            if (posy<0) posy=0;
        }
        if( selectionDisplay.is_keyARROWDOWN )
        {
            posy += 100;
            if (posy>(image.dimy()-canvas.dimy())) posy=image.dimy()-canvas.dimy();
        }
        if( selectionDisplay.is_keyARROWRIGHT )
        {
            posx += 100;
            if (posx>(image.dimx()-canvas.dimx())) posx=image.dimx()-canvas.dimx();
        }
        if( selectionDisplay.is_keyARROWLEFT )
        {
            posx -= 100;
            if (posx<0) posx=0;
        }

        // press "ENTER" if you are content with your selection
        if( selectionDisplay.is_keyENTER || selectionDisplay.is_closed)
        {
            selected = true;
        }

        // resize selection
        if( selectionDisplay.is_resized )
        {
            selectionDisplay.resize();
            if (selectionDisplay.dimx()>image.dimx()) selectionDisplay.resize(image.dimx(),selectionDisplay.dimy());
            if (selectionDisplay.dimy()>image.dimy()) selectionDisplay.resize(selectionDisplay.dimx(),image.dimy());
        }
    }

    // image part to process
    cimg_library::CImg<unsigned char> selected_image( selectionDisplay.dimx(), selectionDisplay.dimy(), 1, 3);
    cimg_forXY(selected_image, x, y)
    {
        if( x + posx >= 0 && x + posx <  image.dimx() && y + posy >= 0 && y + posy < image.dimy() )
        {
            selected_image(x,y,0,0) = image(x + posx,y + posy,0,0);
            selected_image(x,y,0,1) = image(x + posx,y + posy,0,1);
            selected_image(x,y,0,2) = image(x + posx,y + posy,0,2);
        }
    }

    return selected_image;
}

void border_mode(cimg_library::CImg<unsigned char> original_image, int display_x, int display_y, color_type * color, std::string filepath_to_image_selection,
                 bool & end_of_correction)
{
    int dim_x=original_image.dimx();
    int dim_y=original_image.dimy();

    //import image selection
    std::ifstream selection_file(filepath_to_image_selection.c_str());
    std::cout<<"Importing image selection from file:"<<std::endl;
    std::cout<<filepath_to_image_selection<<std::endl;

    int low_x, low_y, high_x, high_y;

    if (selection_file.is_open())
    {
        selection_file>>low_x;
        selection_file>>high_x;
        selection_file>>low_y;
        selection_file>>high_y;

        selection_file.close();

        std::cout<<"...done"<<std::endl;
    }
    else
    {
        std::cout << "File not found, default values used"<<std::endl;;

        low_x=0;
        high_x=dim_x-1;
        low_y=0;
        high_y=dim_y-1;
    }

    int border_position=0;
    bool end_border_mode=false;
    bool saved=true;

    // Initialization of menu
    cimg_library::CImg<unsigned char> menu(250,40,1,3,0);
    cimg_library::CImg<unsigned char> men=menu;

    cimg_library::CImgDisplay  main_menu(menu,"Menu");
    cimg_library::CImgDisplay  zoom_disp(250,250,"ZOOM",0);
    cimg_library::CImgDisplay  main_disp(display_x, display_y,"Border Mode");

    main_menu.move(0,0);
    main_disp.move(0,main_menu.dimy()+57);
    zoom_disp.move(main_disp.dimx()+10,main_menu.dimy()+57);

    float scaling=dim_x/display_x;
    int y_min=0;
    int y_max=scaling*display_y-1;
    int factor=125/scaling;

    while (end_border_mode==false)
    {
        cimg_library::CImg<unsigned char> selected_image=original_image.get_crop(0, y_min, dim_x-1, y_max);
        cimg_library::CImg<unsigned char> selected_image_detail=original_image.get_crop(0, y_min, dim_x-1, y_max);

        int line_width=(int)(scaling-1.0f)/2;

        if(((scaling-1.0f)/2.0f)-line_width>0.0f) line_width++;

        //low y
        if(y_min<=low_y && low_y<=y_max)
        for(int x=0; x<selected_image.dimx(); x++)
        {
            for(int y=std::max(0,low_y-line_width-y_min); y<=std::min(y_max-y_min,low_y+line_width-y_min); y++)
            {
                selected_image(x,y,0,0)=color[3][0];
                selected_image(x,y,0,1)=color[3][1];
                selected_image(x,y,0,2)=color[3][2];
            }

            selected_image_detail(x,low_y-y_min,0,0)=color[3][0];
            selected_image_detail(x,low_y-y_min,0,1)=color[3][1];
            selected_image_detail(x,low_y-y_min,0,2)=color[3][2];
        }

        //high y
        if(y_min<=high_y && high_y<=y_max)
        for(int x=0; x<selected_image.dimx(); x++)
        {
            for(int y=std::max(0,high_y-line_width-y_min); y<=std::min(y_max-y_min,high_y+line_width-y_min); y++)
            {
                selected_image(x,y,0,0)=color[3][0];
                selected_image(x,y,0,1)=color[3][1];
                selected_image(x,y,0,2)=color[3][2];
            }

            selected_image_detail(x,high_y-y_min,0,0)=color[3][0];
            selected_image_detail(x,high_y-y_min,0,1)=color[3][1];
            selected_image_detail(x,high_y-y_min,0,2)=color[3][2];
        }

        //low x
        for(int y=0; y<selected_image.dimy(); y++)
        {
            for(int x=std::max(0,low_x-line_width); x<=std::min(dim_x-1,low_x+line_width); x++)
            {
                selected_image(x,y,0,0)=color[3][0];
                selected_image(x,y,0,1)=color[3][1];
                selected_image(x,y,0,2)=color[3][2];
            }

            selected_image_detail(low_x,y,0,0)=color[3][0];
            selected_image_detail(low_x,y,0,1)=color[3][1];
            selected_image_detail(low_x,y,0,2)=color[3][2];
        }

        //high x
        for(int y=0; y<selected_image.dimy(); y++)
        {
            for(int x=std::max(0,high_x-line_width); x<=std::min(dim_x-1,high_x+line_width); x++)
            {
                selected_image(x,y,0,0)=color[3][0];
                selected_image(x,y,0,1)=color[3][1];
                selected_image(x,y,0,2)=color[3][2];
            }

            selected_image_detail(high_x,y,0,0)=color[3][0];
            selected_image_detail(high_x,y,0,1)=color[3][1];
            selected_image_detail(high_x,y,0,2)=color[3][2];
        }

        //this disables to change the window sizes
        if (main_disp.is_resized)
        {
            main_disp.resize(main_disp);
        }
        if (main_menu.is_resized)
        {
            main_menu.resize(main_menu);
        }
        if (zoom_disp.is_resized)
        {
            zoom_disp.resize(zoom_disp);
        }

        //ZOOM IMAGE STUFF
        //DEFINES THE ZOOM REGION TO GET A BIGGER ZOOM REGION YOU CAN MAKE factor bigger than 15
        int x_0=scaling*(main_disp.mouse_x-factor);
        int y_0=scaling*(main_disp.mouse_y-factor);
        int x_1=scaling*(main_disp.mouse_x+factor);
        int y_1=scaling*(main_disp.mouse_y+factor);
        cimg_library::CImg<unsigned char> visu;

        //VISU IS THE ZOOMED IMAGE
        if (main_disp.mouse_y>=0 && main_disp.mouse_x>=0)
        {
            visu=selected_image_detail.get_crop(x_0,y_0,x_1,y_1);
            zoom_disp.display(visu);
        }

        //MENU STUFF
        men=menu;
        men.draw_text(2,0,"Move image border, or change selection!",color[3],0,1,11,1,1);

        std::string name;
        if(border_position==0) name="low y: ";
        if(border_position==1) name="high y: ";
        if(border_position==2) name="low x: ";
        if(border_position==3) name="high x: ";

        std::ostringstream Str;
        if(border_position==0) Str << low_y;
        if(border_position==1) Str << high_y;
        if(border_position==2) Str << low_x;
        if(border_position==3) Str << high_x;
        std::string temp_string(Str.str());

        name.append(temp_string);
        men.draw_text(2,10,name.c_str(),color[0],0,1,11,1,1);

        if(saved) men.draw_text(2,20,"Saved",color[0],0,1,11,1,1);
        else men.draw_text(2,20,"Unsaved",color[1],0,1,11,1,1);

        main_menu.display(men);

        main_disp.display(selected_image);
        main_disp.wait();

        //this disables to change the window sizes
        if (main_disp.is_resized)
        {
            main_disp.resize(main_disp);
        }
        if (main_menu.is_resized)
        {
            main_menu.resize(main_menu);
        }
        if (zoom_disp.is_resized)
        {
            zoom_disp.resize(zoom_disp);
        }

        //TESTS IF B-KEY (BORDER CHANGE) IS PRESSED
        if(main_disp.key==98||main_disp.key==66)
        {
            border_position=(border_position+1)%4;
            end_border_mode=false;
        }

        //TESTS IF C-KEY (CORRECTION MODE) IS PRESSED
        if(main_disp.key==99||main_disp.key==67)
        {
            end_of_correction=false;
            end_border_mode=true;
        }

        //TEST IF Q-KEY (QUIT) IS PRESSED
        if(main_disp.key==113||main_disp.key==81 || main_disp.is_closed)
        {
            end_of_correction=true;
            end_border_mode=true;
        }

        //TESTS IF S-KEY (SAVE) IS PRESSED
        if(main_disp.key==115||main_disp.key==83)
        {
            std::ofstream selection_file(filepath_to_image_selection.c_str());
            selection_file <<low_x<<" "<<high_x<<" "<<low_y<<" "<<high_y<< "\n";
            selection_file.close();

            end_border_mode=false;
            saved=true;
        }

        //TESTS IF U-KEY (UP) IS PRESSED
        if(main_disp.key==117||main_disp.key==85)
        {
            y_min=std::max(0,y_min-(int)(scaling*display_y));
            y_max=y_min+scaling*display_y;
        }

        //TESTS IF D-KEY (DOWN) IS PRESSED
        if(main_disp.key==100||main_disp.key==68)
        {
            y_max=std::min(dim_y-1,y_max+(int)(scaling*display_y));
            y_min=y_max-scaling*display_y;
        }

        //lower value by 1
        if(main_disp.is_keyARROWLEFT || main_disp.is_keyARROWUP)
        {
            if(border_position==0) low_y=std::max(0,low_y-1);
            if(border_position==1) high_y=std::max(0,high_y-1);
            if(border_position==2) low_x=std::max(0,low_x-1);
            if(border_position==3) high_x=std::max(0,high_x-1);
            saved=false;
        }

        //higher value by 1
        if(main_disp.is_keyARROWRIGHT || main_disp.is_keyARROWDOWN)
        {
            if(border_position==0) low_y=std::min(dim_y-1,low_y+1);
            if(border_position==1) high_y=std::min(dim_y-1,high_y+1);
            if(border_position==2) low_x=std::min(dim_x-1,low_x+1);
            if(border_position==3) high_x=std::min(dim_x-1,high_x+1);
            saved=false;
        }

        //lower value by 10
        if(main_disp.is_keyINSERT)
        {
            if(border_position==0) low_y=std::max(0,low_y-10);
            if(border_position==1) high_y=std::max(0,high_y-10);
            if(border_position==2) low_x=std::max(0,low_x-10);
            if(border_position==3) high_x=std::max(0,high_x-10);
            saved=false;
        }

        //higher value by 10
        if(main_disp.is_keyDELETE)
        {
            if(border_position==0) low_y=std::min(dim_y-1,low_y+10);
            if(border_position==1) high_y=std::min(dim_y-1,high_y+10);
            if(border_position==2) low_x=std::min(dim_x-1,low_x+10);
            if(border_position==3) high_x=std::min(dim_x-1,high_x+10);
            saved=false;
        }

        //lower value by 100
        if(main_disp.is_keyPAGEUP)
        {
            if(border_position==0) low_y=std::max(0,low_y-100);
            if(border_position==1) high_y=std::max(0,high_y-100);
            if(border_position==2) low_x=std::max(0,low_x-100);
            if(border_position==3) high_x=std::max(0,high_x-100);
            saved=false;
        }

        //higher value by 100
        if(main_disp.is_keyPAGEDOWN)
        {
            if(border_position==0) low_y=std::min(dim_y-1,low_y+100);
            if(border_position==1) high_y=std::min(dim_y-1,high_y+100);
            if(border_position==2) low_x=std::min(dim_x-1,low_x+100);
            if(border_position==3) high_x=std::min(dim_x-1,high_x+100);
            saved=false;
        }
    }
}

void correct_prediction(std::string filepath_to_feature_file,std::string path_to_image,std::string path_to_ws_image,
                        std::string path_to_output_folder,std::string param_file, ParameterFile paramFile, int mode=0)
{
    std::string filepath_to_ws_image=path_to_ws_image;
    filepath_to_ws_image.append(get_filename(filepath_to_feature_file));
    //remove the ".bin"
    filepath_to_ws_image.resize(filepath_to_ws_image.size()-4);
    filepath_to_ws_image.append(".h5");

    std::string filepath_to_image_selection=filepath_to_feature_file;
    //remove the ".bin"
    filepath_to_image_selection.resize(filepath_to_image_selection.size()-4);
    filepath_to_image_selection.append(".selection.dat");

    //IMPORT RESULTS FROM HDF5 file
    //all variables represent segmentation and don't need changes here
    seg segment(true);
    segment.load_cgp_data_structure(filepath_to_ws_image);

    vigra::BasicImage<unsigned int> & ws_region_image = segment.ws_region_image;
    marray::Marray<unsigned int> & one_boundings = segment.one_boundings;
    marray::Marray<unsigned int> & two_boundings = segment.two_boundings;
    std::vector< std::vector<point> > & arcs = segment.arcs;
    std::vector<point> & junctions = segment.junctions;
    int & dim_x = segment.dim_x;
    int & dim_y = segment.dim_y;

    color_type * color = new color_type[6];
    color[0][0] = 0;    color[0][1] = 255;  color[0][2] = 0;
    color[1][0] = 255;  color[1][1] = 0;    color[1][2] = 0;
    color[2][0] = 0;    color[2][1] = 100;  color[2][2] = 255;
    color[3][0] = 255;  color[3][1] = 255;  color[3][2] = 255;//white
    color[4][0] = 0;    color[4][1] = 0;    color[4][2] = 0;//black
    color[5][0] = 255;  color[5][1] = 255;  color[5][2] = 0;//area selection

    //images for GUI, will be updated
    std::string filepath_to_image=path_to_image;
    filepath_to_image.append(get_filename(filepath_to_feature_file));
    //remove the ".bin"
    filepath_to_image.resize(filepath_to_image.size()-4);

    //temp CImg file is used to avoid an error
    cimg_library::CImg<unsigned char> temp_image(filepath_to_image.c_str());
    cimg_library::CImg<unsigned char> image=temp_image;
    cimg_library::CImg<unsigned char> original_image=temp_image;

    Parameter<int> display_x;
    display_x.assign("", "display_x", 900);
    display_x.load(paramFile,"config");

    Parameter<int> display_y;
    display_y.assign("", "display_y", 600);
    display_y.load(paramFile,"config");

    display_x=std::min(display_x(),original_image.dimx());
    display_y=std::min(display_y(),original_image.dimy());

    bool border_modification=false;

    //only border modification selected
    if(mode==1)
    {
        border_mode(original_image, display_x, display_y, color, filepath_to_image_selection, border_modification);
        return;
    }
    
    gbn grainBoundNet;
    size_t & nr_new_areas =                                 grainBoundNet.nr_new_areas;
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
    std::vector<bool> & subgrain =                          grainBoundNet.subgrain_arcs;

    std::string param_file_name=get_filename(param_file);
    if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
    else param_file_name = "";

    std::string filepath_bubbles=path_to_output_folder;
    std::string filepath_grains=path_to_output_folder;
    std::string filepath_new_classification=path_to_output_folder;

    filepath_bubbles.append("bubbles/");
    filepath_grains.append("grains/");

    filepath_bubbles.append(param_file_name.c_str());
    filepath_grains.append(param_file_name.c_str());
    filepath_new_classification.append(param_file_name.c_str());

    filepath_to_ws_image.resize(filepath_to_ws_image.size()-3);
    filepath_bubbles.append(get_filename(filepath_to_ws_image));
    filepath_grains.append(get_filename(filepath_to_ws_image));
    filepath_new_classification.append(get_filename(filepath_to_ws_image));

    filepath_bubbles.append(".jpg");
    filepath_grains.append(".jpg");
    filepath_new_classification.append(".h5");

    //check if file exists
    FILE *new_classification;
    new_classification=fopen(filepath_new_classification.c_str(),"rb");
    if(new_classification==NULL)//file does NOT exist
    {
        std::cout<<"IceGrain prediction "<<filepath_new_classification<<" does not exist, all arcs will be unclassified."<<std::endl;

        nr_new_areas=1;
        bubble_area_size = new long[nr_new_areas];
        bubble_area_size[0]=0;
        grain_area_size = new long[nr_new_areas];
        grain_area_size[0]=0;

        std::vector<int> arc_index;
        grain_arc_index.push_back(arc_index);
        bubble_arc_index.push_back(arc_index);

        int nr_old_labels=1;
        for(int y=0;y<dim_y;y++)
            for(int x=0;x<dim_x;x++)
               if (ws_region_image(x,y)>nr_old_labels) nr_old_labels=ws_region_image(x,y);
        region_labels.resize(nr_old_labels,1);

        point p;
        p.x=0;
        p.y=0;

        grain_area_center_mass.push_back(p);
        grain_area_center_mass.push_back(p);
        bubble_area_center_mass.push_back(p);
        bubble_area_center_mass.push_back(p);

        grain_arc.resize(1,false);

        if (mode==0) border_modification=true;
    }
    else//file exists
    {
        fclose(new_classification);
        //IMPORT RESULTS FROM HDF5 file
        //all variables need to be updated!!
        grainBoundNet.load_final_structure(filepath_new_classification);
    }

    filepath_new_classification.resize(filepath_new_classification.size()-3);

    std::vector< std::vector<int> > region_labels_reverse;

    std::vector<int> arc_class;
    arc_class.resize(arcs.size());

    std::vector<int> found_bubble_areas;//in new areas

    std::vector < std::vector<point> > new_areas;
    new_areas.resize(nr_new_areas);

    //create mapping from new to old labels
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
            if(grain_arc[arc])
            {                
                arc_class[arc]=1;//grain boundary
            }                           
            else arc_class[arc]=0;//no boundary            
        }        
        else
        {
            arc_class[arc]=0;//no boundary
            grain_arc.push_back(false);
        }

        if (arc>=subgrain.size())
        {
            subgrain.push_back(false);
        }
        if(subgrain[arc])
        {            
            arc_class[arc]= 5; //Subgrain boundary
        }    
    }

    //set found bubble areas
    for (int new_area=0; new_area<nr_new_areas; new_area++)
    {
        if (bubble_arc_index[new_area].size()>0)
        {
            found_bubble_areas.push_back(new_area);
        }

        for (int a=0; a<bubble_arc_index[new_area].size(); a++)      
        {
            arc_class[bubble_arc_index[new_area][a]]=2;//bubble boundary
        }
    }

    //create new areas vector
    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
           point p;
           p.x=x;
           p.y=y;

           if (bubble_area_size[region_labels[ws_region_image(x,y)-1]-1]!=grain_area_size[region_labels[ws_region_image(x,y)-1]-1])
           {
               new_areas[region_labels[ws_region_image(x,y)-1]-1].push_back(p);//fill new areas vector
           }

           if (bubble_area_size[region_labels[ws_region_image(x,y)-1]-1]>0)
           {
               image(x,y,0,0)=color[3][0];//display bubble areas
               image(x,y,0,1)=color[3][1];
               image(x,y,0,2)=color[3][2];
           }
        }
    }

    for(int arc=0; arc<arcs.size(); arc++)//display arcs
    {
        for(int i=0; i<arcs[arc].size(); i++)
        {
            int x=arcs[arc][i].x;
            int y=arcs[arc][i].y;
            image(x,y,0,0)=color[arc_class[arc]][0];
            image(x,y,0,1)=color[arc_class[arc]][1];
            image(x,y,0,2)=color[arc_class[arc]][2];
        }
    }

    int x_mouse=0;
    int y_mouse=0;
    int posx = 0;
    int posy = 0;

    //wenn end of correction true ist ist die korrektur zuende
    bool end_of_correction=false;
    int  border=1;
    bool modify_section=true;
    bool saved=true;

    // Initialization of menu
    cimg_library::CImg<unsigned char> menu(250,40,1,3,0);

    cimg_library::CImg<unsigned char> selected_image=image;
    cimg_library::CImg<unsigned char> selected_original_image=original_image;
    cimg_library::CImg<unsigned char> men=menu;
    int factor=10;

    while(end_of_correction==false)
    {
        if(border_modification)//change image border
        {
            border_mode(original_image, display_x, display_y, color, filepath_to_image_selection, end_of_correction);
            border_modification=false;
        }

        if(end_of_correction) exit(0);

        if(modify_section)//modify image selection for arc/area correction
        {
            if(image.dimx()>display_x || image.dimy()>display_y)
            {
                selected_image=select_image_section(image, display_x, display_y, posx, posy);
                selected_original_image.resize(selected_image.dimx(),selected_image.dimy());

                for (int x=0; x<selected_image.dimx(); x++)
                    for (int y=0; y<selected_image.dimy(); y++)
                        if( x + posx >= 0 && x + posx <  image.dimx() && y + posy >= 0 && y + posy < image.dimy() )
                        {
                            selected_original_image(x,y,0,0) = selected_original_image(x,y,0,1) = selected_original_image(x,y,0,2) =
                                original_image(x+posx,y+posy);
                        }
            }

            modify_section=false;
        }

        cimg_library::CImgDisplay  main_menu(menu,"Menu");
        cimg_library::CImgDisplay  zoom_disp(250,250,"ZOOM",0);
        cimg_library::CImgDisplay  zoom_disp_img(250,250,"Original Image",0);
        cimg_library::CImgDisplay  main_disp(selected_image,"Click a boundary or area");

        main_menu.move(0,0);
        main_disp.move(0,main_menu.dimy()+57);
        zoom_disp.move(main_disp.dimx()+10,main_menu.dimy()+57);
        zoom_disp_img.move(main_disp.dimx()+10,main_menu.dimy()+zoom_disp.dimy()+90);

        while (modify_section==false && end_of_correction==false && border_modification==false)
        {
            //this disables to change the window sizes
            if (main_disp.is_resized)
            {
                main_disp.resize(main_disp);
            }
            if (main_menu.is_resized)
            {
                main_menu.resize(main_menu);
            }
            if (zoom_disp.is_resized)
            {
                zoom_disp.resize(zoom_disp);
            }
            if (zoom_disp_img.is_resized)
            {
                zoom_disp_img.resize(zoom_disp_img);
            }

            //ZOOM IMAGE STUFF
            //DEFINES THE ZOOM REGION TO GET A BIGGER ZOOM REGION YOU CAN MAKE factor bigger than 15
            int x_0=x_mouse-factor;
            int y_0=y_mouse-factor;
            int x_1=x_mouse+factor;
            int y_1=y_mouse+factor;
            cimg_library::CImg<unsigned char> visu;
            cimg_library::CImg<unsigned char> visu_img;

            //VISU IS THE ZOOMED IMAGE
            if(border==0) visu =selected_image.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[0],0.5f);
            if(border==1) visu =selected_image.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[1],0.5f);
            if(border==2) visu =selected_image.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[2],0.5f);
            if(border==3) visu =selected_image.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[1],0.5f);
            if(border==4) visu =selected_image.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[5],0.5f);
            if(border==5) visu =selected_image.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[5],0.5f);                   
 
            //VISU_IMG IS ORIGINAL IMAGE
            visu_img=selected_original_image.get_crop(x_0,y_0,x_1,y_1).draw_point(x_mouse-x_0,y_mouse-y_0,color[5],0.5f);

            zoom_disp.display(visu);
            zoom_disp_img.display(visu_img);

            //MENU STUFF
            men=menu;

            if(border==0) men.draw_text(2,0,"Select grain boundary for merging!",color[0],0,1,11,1,1);
            if(border==1) men.draw_text(2,0,"Add new grain area!",color[1],0,1,11,1,1);
            if(border==2) men.draw_text(2,0,"Grain to bubble area!",color[2],0,1,11,1,1);
            if(border==3) men.draw_text(2,0,"Bubble to grain area!",color[1],0,1,11,1,1);
            if(border==4) men.draw_text(2,0,"Boundary to subgrain boundary!", color[5],0,1,11,1,1);
                            //men.draw_text(2,20,"Klick on a boundary to convert it into a subgrain boundary", color[3],0.1,1,11,1,1); }
            if(border==5) men.draw_text(2,0,"Subgrain boundary to no boundary!",color[0],0,1,11,1,1);

            if(saved) men.draw_text(2,10,"Saved",color[0],0,1,11,1,1);
            else men.draw_text(2,10,"Unsaved",color[1],0,1,11,1,1);

            main_menu.display(men);

            main_disp.display(selected_image);
            main_disp.wait();

            //this disables to change the window sizes
            if (main_disp.is_resized)
            {
                main_disp.resize(main_disp);
            }
            if (main_menu.is_resized)
            {
                main_menu.resize(main_menu);
            }
            if (zoom_disp.is_resized)
            {
                zoom_disp.resize(zoom_disp);
            }
            if (zoom_disp_img.is_resized)
            {
                zoom_disp_img.resize(zoom_disp_img);
            }

            //0
            if(main_disp.key==48)
            {
                border=0;
                end_of_correction=false;
            }

            //1
            if(main_disp.key==49)
            {
                border=1;
                end_of_correction=false;
            }

            //2
            if(main_disp.key==50)
            {
                border=2;
                end_of_correction=false;
            }

            //3
            if(main_disp.key==51)
            {
                border=3;
                end_of_correction=false;
            }

            //4
            if(main_disp.key==52)
            {
                border = 4;
                end_of_correction = false;
            }

            //5
            if(main_disp.key==53)
            {
                border = 5;
                end_of_correction = false;
            }
            
            //TESTS IF S-KEY (SAVE) IS PRESSED
            if(main_disp.key==115||main_disp.key==83)
            {
                end_of_correction=false;
                saved=true;

                std::cout<<"Please wait! Calculation of new structure and saving..."<<std::endl;
                
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

                        for(size_t area=0;area<nr_new_areas && !found;area++)
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

                        for(size_t area=0;area<nr_new_areas && !found;area++)
                        {
                            for (int a=0;a<(int)grain_arc_index[area].size() && !found;a++)//grain boundaries
                            {
                                int arc_index=grain_arc_index[area][a]+1;
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
                    if ((nr_bubble_arcs==2 || nr_bubble_arcs==3) && (nr_grain_arcs==1 || nr_grain_arcs==2))
                    {
                        grain_bubble_junctions.push_back(y);
                    }
                }
               
                //create output images
                vigra::FImage bubble_image(dim_x,dim_y);//bubbles areas and bubble boundaries
                vigra::FImage grain_image(dim_x,dim_y);//grain areas gray, bubbles and grain boundaries black
                
                for(int y=0;y<dim_y;y++)
                {
                    for(int x=0;x<dim_x;x++)
                    {
                        grain_image(x,y)=3;
                    }
                }
 
                //loop over new areas, bubble/grain areas and grain center of mass to output image
                for(size_t area=0;area<nr_new_areas;area++)
                {                    
                    if (bubble_area_size[area]>0)//bubble area
                    {                        
                        for(size_t k=0; k<new_areas[area].size(); ++k)
                        {
                            size_t x=new_areas[area][k].x;
                            size_t y=new_areas[area][k].y;
                            bubble_image(x,y)=1;
                            grain_image(x,y)=0;
                        }
                    }
                    else if (grain_area_size[area]>0)//grain area
                    {                        
                        for(size_t k=0; k<new_areas[area].size(); ++k)
                        {
                            size_t x=new_areas[area][k].x;
                            size_t y=new_areas[area][k].y;
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
                /*
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
                */
/*DEBUG OUTPUT_SAVE*/std::cout << "Loop over arc lists..."<< std::endl;                  
                //loop over arc lists, boundaries of bubble/grain boundaries to output image
                for (int area=0;area<nr_new_areas;area++)//loop over areas
                {
                    for (int a=0;a<(int)bubble_arc_index[area].size();a++)//bubble boundaries
                    {
                        //now we loop over the points in this arc
                        for(int p=0;p<(int)arcs[bubble_arc_index[area][a]].size();p++)
                        {
                            int x=arcs[bubble_arc_index[area][a]][p].x;
                            int y=arcs[bubble_arc_index[area][a]][p].y;
                            bubble_image(x,y)=2;
                        }
                    }
  
                    for (int a=0;a<(int)grain_arc_index[area].size();a++)//grain boundaries
                    {
                        //if(grain_arc_index[area][a] == NULL || grain_arc_index[area] == NULL) std::cout << "NOOO!" << std::endl;
                        //now we loop over the points in this arc
                        for(int p=0;p<(int)arcs[grain_arc_index[area][a]].size();p++)
                        {
                            int x=arcs[grain_arc_index[area][a]][p].x;
                            int y=arcs[grain_arc_index[area][a]][p].y;
                            grain_image(x,y)=0;
                        }
                    }
                }
                                 
                //if subgrain boundaries are labeled show them in output grain image
                for(int arc=0; arc<arcs.size(); arc++)
                {
                    if (subgrain[arc])
                    {
                        for (int i=0; i<arcs[arc].size(); i++)
                        {
                            int x=arcs[arc][i].x;
                            int y=arcs[arc][i].y;
                            grain_image(x,y)=2;
                        }
                    }
                }
                                                
                exportImage(srcImageRange(bubble_image), vigra::ImageExportInfo(filepath_bubbles.c_str()));
                exportImage(srcImageRange(grain_image), vigra::ImageExportInfo(filepath_grains.c_str()));
             
                grainBoundNet.save_final_structure(filepath_new_classification);

                std::cout<<"...done"<<std::endl;
            }

            //TESTS IF Q-KEY (QUIT) IS PRESSED OR WINDOWS IS CLOSED
            if(main_disp.key==113||main_disp.key==81 || main_disp.is_closed)
            {
                end_of_correction=true;
            }

            //TESTS IF B-KEY (BORDER MODIFICATION) IS PRESSED
            if(main_disp.key==98||main_disp.key==66)
            {
                main_disp.close();
                main_menu.close();
                border_modification=true;
            }

            //TEST IF M-KEY (MODIFY SECTION) IS PRESSED
            if(main_disp.key==109||main_disp.key==77)
            {
                main_disp.close();
                main_menu.close();
                modify_section=true;
            }

            //Test if the mouse button is clicked on the image area 
            if (main_disp.mouse_y>=0 && main_disp.mouse_x>=0)
            {
                //mouse stored to  x_mouse and y_mouse to use it out of the loop
                x_mouse = main_disp.mouse_x;
                y_mouse = main_disp.mouse_y;
                if(main_disp.button||main_disp.key==97||main_disp.key==65)//mouse clicked or A-KEY (add) pressed
                {
                    if (border==0)//select grain arc to be merged
                    {
                        int second_x=main_disp.mouse_x;
                        int second_y=main_disp.mouse_y;

                        while(main_disp.button)
                        {
                            if(main_disp.mouse_y>=0 && main_disp.mouse_x>=0)
                            {
                                second_x=main_disp.mouse_x;
                                second_y=main_disp.mouse_y;

                                cimg_library::CImg<unsigned char> temp_selected_image=selected_image;

                                for(int x=std::min(x_mouse,second_x); x<=std::max(x_mouse,second_x); x++)
                                {
                                    temp_selected_image(x,y_mouse,0,0)=color[0][0];
                                    temp_selected_image(x,y_mouse,0,1)=color[0][1];
                                    temp_selected_image(x,y_mouse,0,2)=color[0][2];
                                }
                                for(int x=std::min(x_mouse,second_x); x<=std::max(x_mouse,second_x); x++)
                                {
                                    temp_selected_image(x,second_y,0,0)=color[0][0];
                                    temp_selected_image(x,second_y,0,1)=color[0][1];
                                    temp_selected_image(x,second_y,0,2)=color[0][2];
                                }
                                for(int y=std::min(y_mouse,second_y); y<=std::max(y_mouse,second_y); y++)
                                {
                                    temp_selected_image(x_mouse,y,0,0)=color[0][0];
                                    temp_selected_image(x_mouse,y,0,1)=color[0][1];
                                    temp_selected_image(x_mouse,y,0,2)=color[0][2];
                                }
                                for(int y=std::min(y_mouse,second_y); y<=std::max(y_mouse,second_y); y++)
                                {
                                    temp_selected_image(second_x,y,0,0)=color[0][0];
                                    temp_selected_image(second_x,y,0,1)=color[0][1];
                                    temp_selected_image(second_x,y,0,2)=color[0][2];
                                }

                                main_disp.display(temp_selected_image);
                            }
                        }

                        //HERE WE HAVE TO CHECK IF THE POINT (X,Y) IS A POINT OF A GRAIN ARC
                        for(int arc_count=0; arc_count<(int)arcs.size(); arc_count++)
                        {
                            bool found_point=false;

                            for(int point_count=0; point_count<(int)arcs[arc_count].size() && found_point==false; point_count++)
                            {
                                if(arcs[arc_count][point_count].x>=std::min(x_mouse,second_x)+posx &&
                                    arcs[arc_count][point_count].y>=std::min(y_mouse,second_y)+posy &&
                                    arcs[arc_count][point_count].x<=std::max(x_mouse,second_x)+posx &&
                                    arcs[arc_count][point_count].y<=std::max(y_mouse,second_y)+posy)
                                {
                                    found_point=true;

                                    if(arc_class[arc_count]==1)//grain arc
                                    {
                                        bool outside=false;
                                        int outside_label;
                                        bool bubble_merged=false;

                                        //these two labels are selected for merging
                                        int label_merged = std::min(region_labels[two_boundings(arc_count,0)-1]-1,
                                                                    region_labels[two_boundings(arc_count,1)-1]-1);
                                        int label_removed = std::max(region_labels[two_boundings(arc_count,0)-1]-1,
                                                                     region_labels[two_boundings(arc_count,1)-1]-1);

                                        if (bubble_area_size[label_removed]>0 || bubble_area_size[label_merged]>0) bubble_merged=true;

                                        //look for arcs of these two labels -> no longer grain arcs
                                        for(int y=0;y<(int)two_boundings.shape(0);y++)
                                        {
                                            if ((region_labels[two_boundings(y,0)-1]-1==label_merged &&
                                                    region_labels[two_boundings(y,1)-1]-1==label_removed) ||
                                                (region_labels[two_boundings(y,1)-1]-1==label_merged &&
                                                    region_labels[two_boundings(y,0)-1]-1==label_removed))
                                            {
                                                arc_class[y]=0;
                                                grain_arc[y]=false;
                                                std::cout<<"remove grain arc "<<y<<std::endl;
                                            }
                                        }

                                        //grain belongs now to outside -> remove as grain
                                        if(bubble_area_size[label_merged]==grain_area_size[label_merged])
                                        {
                                            outside=true;
                                            grain_arc_index[label_merged].clear();

                                            outside_label=label_merged+1;

                                            std::cout<<"label "<<label_removed+1<<" will be removed and merged with outside label "<<outside_label
                                                <<std::endl;
                                        }
                                        else if(bubble_area_size[label_removed]==grain_area_size[label_removed])
                                        {
                                            outside=true;
                                            grain_arc_index[label_removed].clear();

                                            outside_label=label_removed+1;
                                            label_removed=label_merged;

                                            std::cout<<"label "<<label_removed+1<<" will be removed and merged with outside label "<<outside_label
                                                <<std::endl;
                                            outside_label--;
                                        }
                                        else
                                        {
                                            //combine arc lists
                                            std::list<int> combined_arc_list;

                                            if (grain_arc_index[label_merged].size()>0)
                                                for (int a=0; a<grain_arc_index[label_merged].size(); a++)
                                                    combined_arc_list.push_back(grain_arc_index[label_merged][a]);
                                            else for (int a=0; a<bubble_arc_index[label_merged].size(); a++)
                                                    combined_arc_list.push_back(bubble_arc_index[label_merged][a]);

                                            if (grain_arc_index[label_removed].size()>0)
                                                for (int a=0; a<grain_arc_index[label_removed].size(); a++)
                                                    combined_arc_list.push_back(grain_arc_index[label_removed][a]);
                                            else for (int a=0; a<bubble_arc_index[label_removed].size(); a++)
                                                    combined_arc_list.push_back(bubble_arc_index[label_removed][a]);

                                            //sort combined list und remove double entries
                                            combined_arc_list.sort();
                                            combined_arc_list.unique();

                                            //replace arc list of label_merged
                                            if (!bubble_merged)
                                            {
                                                grain_arc_index[label_merged].clear();
                            
                                                for (std::list<int>::iterator a=combined_arc_list.begin(); a!=combined_arc_list.end();++a) 
                                                {
                                                    if (arc_class[*a]>0) grain_arc_index[label_merged].push_back(*a);//check if bubble or grain arc
                                                }
                                            }
                                            else
                                            {
                                                bubble_arc_index[label_merged].clear();
                            
                                                for (std::list<int>::iterator a=combined_arc_list.begin(); a!=combined_arc_list.end();++a) 
                                                {
                                                    if (arc_class[*a]>0) bubble_arc_index[label_merged].push_back(*a);//check if bubble or grain arc
                                                }
                                            }

                                            std::cout<<"label "<<label_removed+1<<" will be removed and merged with label "<<label_merged+1<<std::endl;
                                        }

                                        //erase label_removed -> decrease all labels higher than label_removed by one
                                        grain_arc_index.erase(grain_arc_index.begin()+label_removed);
                                        bubble_arc_index.erase(bubble_arc_index.begin()+label_removed);

                                        for(int l=0; l<region_labels.size(); l++)
                                            if (region_labels[l]>label_removed) region_labels[l]--;

                                        //update areas and area size
                                        if (!outside)
                                        {
                                            for (int p=0; p<new_areas[label_removed].size(); p++)
                                                new_areas[label_merged].push_back(new_areas[label_removed][p]);
                                            grain_area_size[label_merged]=new_areas[label_merged].size();
                                        }

                                        //new_areas.erase(new_areas.begin()+label_removed);
                                        //"erase" doesn't free memory, so swaping is necessary
                                        std::vector< std::vector<point> > areas2;

                                        for (int area=label_removed+1; area<new_areas.size(); area++)
                                            areas2.push_back(new_areas[area]);

                                        new_areas.resize(label_removed);

                                        for (int area=0; area<areas2.size(); area++)
                                            new_areas.push_back(areas2[area]);

                                        areas2.clear();

                                        grain_area_center_mass.erase(grain_area_center_mass.begin()+label_removed);
                                        bubble_area_center_mass.erase(bubble_area_center_mass.begin()+label_removed);
                                        found_bubble_areas.clear();
                                        nr_new_areas--;

                                        delete grain_area_size;
                                        delete bubble_area_size;
                                        grain_area_size = new long[nr_new_areas];
                                        bubble_area_size = new long[nr_new_areas];

                                        for (int area=0; area<nr_new_areas; area++)
                                        {
                                            if (grain_arc_index[area].size()>0) grain_area_size[area]=new_areas[area].size();
                                            else grain_area_size[area]=0;

                                            if (bubble_arc_index[area].size()>0)
                                            {
                                                bubble_area_size[area]=new_areas[area].size();
                                                found_bubble_areas.push_back(area);
                                            }
                                            else bubble_area_size[area]=0;
                                        }

                                        if (!outside)
                                        {
                                            //calculate new center of mass of merged grain/bubble
                                            long area_x_sum=0;
                                            long area_y_sum=0;

                                            //loop over all pixels with label
                                            for (int pixel=0;pixel<new_areas[label_merged].size();pixel++)
                                            {
                                                point p=new_areas[label_merged][pixel];
                                                area_x_sum+=p.x;
                                                area_y_sum+=p.y;
                                            }

                                            int xx=area_x_sum/new_areas[label_merged].size();
                                            int yy=area_y_sum/new_areas[label_merged].size();

                                            if (!bubble_merged)
                                            {
                                                grain_area_center_mass[label_merged].x=xx;
                                                grain_area_center_mass[label_merged].y=yy;
                                            }
                                            else
                                            {
                                                bubble_area_center_mass[label_merged].x=xx;
                                                bubble_area_center_mass[label_merged].y=yy;
                                            }

                                            //update region_labels and region_labels_reverse
                                            for(int l=0; l<region_labels_reverse[label_removed].size(); l++)
                                            {
                                                region_labels_reverse[label_merged].push_back(region_labels_reverse[label_removed][l]);
                                                region_labels[region_labels_reverse[label_removed][l]-1]=label_merged+1;
                                            }

                                            region_labels_reverse.erase(region_labels_reverse.begin()+label_removed);
                                        }
                                        else
                                        {
                                            //update region_labels and region_labels_reverse
                                            for(int l=0; l<region_labels_reverse[label_removed].size(); l++)
                                            {
                                                region_labels[region_labels_reverse[label_removed][l]-1]=outside_label;
                                            }

                                            region_labels_reverse.erase(region_labels_reverse.begin()+label_removed);
                                        }

                                        saved=false;
                                    }
                                }
                            }
                        }
                    }

                    else if (border==1)//select area to be split of as grain
                    {
                        int old_label=ws_region_image(x_mouse+posx,y_mouse+posy)-1;
                        int label=region_labels[ws_region_image(x_mouse+posx,y_mouse+posy)-1]-1;

                        //HERE WE HAVE TO CHECK IF THE POINT (X,Y) IS A POINT OF A GRAIN AREA OR OUTSIDE
                        if (grain_arc_index[label].size()>0 || bubble_area_size[label]==grain_area_size[label])
                        {
                            bool connection=true;

                            if (region_labels_reverse[label].size()>1)
                            {
                                //when the selected area would disconnect the background or a grain it can not be split off
                                vigra::BasicImage<bool> temp_img(dim_x,dim_y);
                                int min_x=0;
                                int min_y=0;

                                //find a minimal region that completly contains all pixels with label
                                if(bubble_area_size[label]!=grain_area_size[label])
                                {
                                    min_x=dim_x;
                                    min_y=dim_y;
                                    int max_x=0;
                                    int max_y=0;

                                    //loop over all pixels with label
                                    for (int pixel=0;pixel<new_areas[label].size();pixel++)
                                    {
                                        point p=new_areas[label][pixel];
                                        if (p.x<min_x) min_x=p.x;
                                        if (p.y<min_y) min_y=p.y;
                                        if (p.x>max_x) max_x=p.x;
                                        if (p.y>max_y) max_y=p.y;
                                    }

                                    temp_img.resize(max_x-min_x+1,max_y-min_y+1);
                                }

                                for(int y=min_y;y<min_y+temp_img.height();y++)
                                {
                                    for (int x=min_x;x<min_x+temp_img.width();x++)
                                    {
                                        temp_img(x-min_x,y-min_y)=false;
                                    }
                                }

                                if(bubble_area_size[label]!=grain_area_size[label])
                                {
                                    //loop over all pixels with label
                                    for (int pixel=0;pixel<new_areas[label].size();pixel++)
                                    {
                                        point p=new_areas[label][pixel];

                                        if (ws_region_image(p.x,p.y)-1!=old_label)//pixel belongs NOT to area selected to be split off
                                        {
                                            temp_img(p.x-min_x,p.y-min_y)=true;
                                        }
                                    }
                                }
                                else
                                {
                                    std::cout<<"Please wait! Outside connectivity is being checked..."<<std::endl;

                                    //loop over all pixels, since area vector is not filled for outside, loop over all pixels
                                    for(int y=0;y<dim_y;y++)
                                    {
                                        for(int x=0;x<dim_x;x++)
                                        {
                                            //pixel belongs NOT to area selected to be split off
                                            if (region_labels[ws_region_image(x,y)-1]-1==label && ws_region_image(x,y)-1!=old_label)
                                            {
                                                temp_img(x,y)=true;
                                            }
                                        }
                                    }
                                }

                                int nr_labels;

                                // find connected regions
                                nr_labels=vigra::labelImageWithBackground(srcImageRange(temp_img), destImage(temp_img), true, 0);

                                if (nr_labels>1)
                                {
                                    std::cout<<"Error: this area would disconnect background or another grain!"<<std::endl;
                                    connection=false;
                                }
                                else if(bubble_area_size[label]==grain_area_size[label]) std::cout<<"...done"<<std::endl;
                            }
                            else if (grain_area_size[label] > 0) connection = false;
                            else std::cout<<"Single area of this label, connectivity check is skipped"<<std::endl;

                            if (connection)
                            {
                                //create new label
                                int added_label=nr_new_areas;
                                nr_new_areas++;

                                std::cout<<"create new label "<<added_label+1<<std::endl;

                                new_areas.resize(nr_new_areas);
                                grain_arc_index.resize(nr_new_areas);
                                bubble_arc_index.resize(nr_new_areas);

                                long * new_grain_area_size = new long[nr_new_areas];
                                long * new_bubble_area_size = new long[nr_new_areas];

                                for (int area=0; area<nr_new_areas-1; area++)
                                {
                                    new_grain_area_size[area]=grain_area_size[area];
                                    new_bubble_area_size[area]=bubble_area_size[area];
                                }

                                delete grain_area_size;
                                delete bubble_area_size;

                                grain_area_size = new_grain_area_size;
                                bubble_area_size = new_bubble_area_size;
                                bubble_area_size[added_label]=0;
                                point p;
                                p.x=0;
                                p.y=0;
                                bubble_area_center_mass.push_back(p);

                                grain_area_center_mass.resize(nr_new_areas);
                                region_labels_reverse.resize(nr_new_areas);

                                //loop over all arcs and look for arcs of the area selected
                                for(int y=0;y<(int)two_boundings.shape(0);y++)
                                {
                                    for(int x=0;x<(int)two_boundings.shape(1);x++)
                                    {
                                        if (two_boundings(y,x)==old_label+1)
                                        {
                                            //resize grain_arc if necessary            
                                            if (y+1>grain_arc.size())
                                            {
                                                grain_arc.resize(y+1,false);
                                            }

                                            //check if arc is in grain arc list
                                            int found_arc=-1;

                                            for(int a=0; a<(int)grain_arc_index[label].size() && found_arc==-1; a++)
                                            {
                                                if (grain_arc_index[label][a]==y) found_arc=a;
                                            }

                                            if (found_arc>-1)//arc found in grain arc list
                                            {
                                                grain_arc_index[label].erase(grain_arc_index[label].begin()+found_arc);//erase arc from common list
                                                grain_arc_index[added_label].push_back(y);//add arc to seperate list
                                                if (arc_class[y]==1) std::cout<<"moved grain arc "<<y<<std::endl;
                                                else std::cout<<"moved bubble arc "<<y<<std::endl;
                                            }
                                            //no grain arc list for this label (outside)
                                            else if(bubble_area_size[label]==grain_area_size[label] && grain_arc[y])
                                            {
                                                std::cout<<"copy grain arc "<<y<<std::endl;

                                                //add to grain arc lists of added_label
                                                grain_arc_index[added_label].push_back(y);
                                            }
                                            else if (arc_class[y]==0)//arc is of class 0->new grain arc
                                            {
                                                arc_class[y]=1;
                                                grain_arc[y]=true;
                                                std::cout<<"add grain arc "<<y<<std::endl;

                                                //add to both grain arc lists (if label is not outside label)
                                                if(bubble_area_size[label]!=grain_area_size[label]) grain_arc_index[label].push_back(y);
                                                grain_arc_index[added_label].push_back(y);
                                            }
                                            else if (arc_class[y]==2)//arc is of bubble arc->new grain arc
                                            {
                                                std::cout<<"add bubble arc "<<y<<std::endl;

                                                //add to grain arc lists of added_label
                                                grain_arc_index[added_label].push_back(y);
                                            }
                                            else if (arc_class[y]==5)//Arc is of a subgrain arc->new grain arc
                                            {
                                                arc_class[y] = 1;
                                                grain_arc[y] = true;

                                                //add to both grain arc lists (if label is not outside label)
                                                if(bubble_area_size[label]!=grain_area_size[label]) grain_arc_index[label].push_back(y);
                                                grain_arc_index[added_label].push_back(y);
                                            }

                                            subgrain[y]=false;
                                        }
                                    }
                                }

                                //calculate new area size and center of mass of both grains
                                if(bubble_area_size[label]!=grain_area_size[label])
                                {
                                    grain_area_size[label]=0;
                                    grain_area_size[added_label]=0;

                                    long grain_area_x_sum=0;
                                    long grain_area_y_sum=0;
                                    long grain_area_x_sum_added=0;
                                    long grain_area_y_sum_added=0;

                                    std::vector<point> temp;
                                    temp=new_areas[label];
                                    new_areas[label].clear();

                                    //loop over all pixels with label
                                    for (int pixel=0;pixel<temp.size();pixel++)
                                    {
                                        point p=temp[pixel];

                                        if (ws_region_image(p.x,p.y)-1==old_label)//pixel belongs to area selected to be split off
                                        {
                                            new_areas[added_label].push_back(p);//add to new area
                                            grain_area_size[added_label]++;
                                            grain_area_x_sum_added+=p.x;
                                            grain_area_y_sum_added+=p.y;
                                        }
                                        else
                                        {
                                            grain_area_size[label]++;
                                            grain_area_x_sum+=p.x;
                                            grain_area_y_sum+=p.y;

                                            new_areas[label].push_back(p);
                                        }
                                    }

                                    temp.clear();

                                    if (grain_area_size[added_label]==0) grain_area_size[added_label]++;
                                    if (grain_area_size[label]==0) grain_area_size[label]++;

                                    int xx=grain_area_x_sum_added/grain_area_size[added_label];
                                    int yy=grain_area_y_sum_added/grain_area_size[added_label];
                                    grain_area_center_mass[added_label].x=xx;
                                    grain_area_center_mass[added_label].y=yy;

                                    xx=grain_area_x_sum/grain_area_size[label];
                                    yy=grain_area_y_sum/grain_area_size[label];
                                    grain_area_center_mass[label].x=xx;
                                    grain_area_center_mass[label].y=yy;
                                }
                                else//new label is split off from outside
                                {
                                    grain_area_size[added_label]=0;

                                    long grain_area_x_sum_added=0;
                                    long grain_area_y_sum_added=0;

                                    //loop over all pixels, since area vector is not filled for outside, loop over all pixels
                                    for(int y=0;y<dim_y;y++)
                                    {
                                        for(int x=0;x<dim_x;x++)
                                        {
                                            if (ws_region_image(x,y)-1==old_label)//pixel belongs to area selected to be split off
                                            {
                                                point p;
                                                p.x=x;
                                                p.y=y;

                                                new_areas[added_label].push_back(p);//add to new area
                                                grain_area_size[added_label]++;
                                                grain_area_x_sum_added+=x;
                                                grain_area_y_sum_added+=y;
                                            }
                                        }
                                    }
                                                                                         
                                    int xx=grain_area_x_sum_added/grain_area_size[added_label];
                                    int yy=grain_area_y_sum_added/grain_area_size[added_label];
                                    grain_area_center_mass[added_label].x=xx;
                                    grain_area_center_mass[added_label].y=yy;
                                }

                                //update region_labels and region_labels_reverse
                                region_labels[old_label]=added_label+1;

                                for(int l=0; l<region_labels_reverse[label].size(); l++)
                                {
                                    if(region_labels_reverse[label][l]==old_label+1)
                                        region_labels_reverse[label].erase(region_labels_reverse[label].begin()+l);
                                }

                                region_labels_reverse[added_label].push_back(old_label+1);

                                saved=false;
                            }
                        }

                        //HERE WE HAVE TO CHECK IF THE POINT (X,Y) IS A POINT OF A BUBBLE
                        if (bubble_arc_index[label].size()>0 && region_labels_reverse[label].size()>1)
                        {
                            bool bubble_border=false;

                            for (int arc=0; arc<bubble_arc_index[label].size() && !bubble_border; arc++)
                            {
                                int y=bubble_arc_index[label][arc];
                                if (two_boundings(y,0)-1==old_label || two_boundings(y,1)-1==old_label) bubble_border=true;
                            }

                            if (!bubble_border)
                            {
                                //create new label
                                int added_label=nr_new_areas;
                                nr_new_areas++;

                                std::cout<<"create new label "<<added_label+1<<std::endl;

                                new_areas.resize(nr_new_areas);
                                grain_arc_index.resize(nr_new_areas);
                                bubble_arc_index.resize(nr_new_areas);

                                long * new_grain_area_size = new long[nr_new_areas];
                                long * new_bubble_area_size = new long[nr_new_areas];

                                for (int area=0; area<nr_new_areas-1; area++)
                                {
                                    new_grain_area_size[area]=grain_area_size[area];
                                    new_bubble_area_size[area]=bubble_area_size[area];
                                }

                                delete grain_area_size;
                                delete bubble_area_size;

                                grain_area_size = new_grain_area_size;
                                bubble_area_size = new_bubble_area_size;
                                bubble_area_size[added_label]=0;

                                grain_area_center_mass.resize(nr_new_areas);
                                bubble_area_center_mass.resize(nr_new_areas);
                                region_labels_reverse.resize(nr_new_areas);

                                //loop over all arcs and look for arcs of the area selected
                                for(int y=0;y<(int)two_boundings.shape(0);y++)
                                {
                                    for(int x=0;x<(int)two_boundings.shape(1);x++)
                                    {
                                        if (two_boundings(y,x)==old_label+1)
                                        {
                                            //resize grain_arc if necessary            
                                            if (y+1>grain_arc.size())
                                            {
                                                grain_arc.resize(y+1,false);
                                            }

                                            if (arc_class[y]==0)//arc is of class 0->new grain arc
                                            {
                                                arc_class[y]=1;
                                                grain_arc[y]=true;
                                                std::cout<<"add grain arc "<<y<<std::endl;

                                                //add to grain arc lists of added_label
                                                grain_arc_index[added_label].push_back(y);
                                            }
                                            else if (arc_class[y]==1)//arc is of class 1->next to other inside grain arc
                                            {
                                                std::cout<<"copy grain arc "<<y<<std::endl;

                                                //add to grain arc lists of added_label
                                                grain_arc_index[added_label].push_back(y);
                                            }
                                            else
                                            {
                                                std::cout<<"Error: Connection to bubble border!"<<y<<std::endl;
                                                exit(-1);
                                            }

                                            subgrain[y]=false;
                                        }
                                    }
                                }

                                //calculate new area size and center of mass of grain and bubble
                                bubble_area_size[label]=0;
                                grain_area_size[added_label]=0;

                                long bubble_area_x_sum=0;
                                long bubble_area_y_sum=0;
                                long grain_area_x_sum_added=0;
                                long grain_area_y_sum_added=0;

                                std::vector<point> temp;
                                temp=new_areas[label];
                                new_areas[label].clear();

                                //loop over all pixels with label
                                for (int pixel=0;pixel<temp.size();pixel++)
                                {
                                    point p=temp[pixel];

                                    if (ws_region_image(p.x,p.y)-1==old_label)//pixel belongs to area selected to be split off
                                    {
                                        new_areas[added_label].push_back(p);//add to new area
                                        grain_area_size[added_label]++;
                                        grain_area_x_sum_added+=p.x;
                                        grain_area_y_sum_added+=p.y;
                                    }
                                    else
                                    {
                                        bubble_area_size[label]++;
                                        bubble_area_x_sum+=p.x;
                                        bubble_area_y_sum+=p.y;

                                        new_areas[label].push_back(p);
                                    }
                                }

                                temp.clear();

                                //loop over all pixels with added label
                                for (int pixel=0;pixel<new_areas[added_label].size();pixel++)
                                {
                                    int x=new_areas[added_label][pixel].x;
                                    int y=new_areas[added_label][pixel].y;

                                    image(x,y,0,0)=original_image(x,y);
                                    image(x,y,0,1)=original_image(x,y);
                                    image(x,y,0,2)=original_image(x,y);
                                }

                                int xx=grain_area_x_sum_added/grain_area_size[added_label];
                                int yy=grain_area_y_sum_added/grain_area_size[added_label];
                                grain_area_center_mass[added_label].x=xx;
                                grain_area_center_mass[added_label].y=yy;

                                xx=bubble_area_x_sum/bubble_area_size[label];
                                yy=bubble_area_y_sum/bubble_area_size[label];
                                bubble_area_center_mass[label].x=xx;
                                bubble_area_center_mass[label].y=yy;

                                //update region_labels and region_labels_reverse
                                region_labels[old_label]=added_label+1;

                                for(int l=0; l<region_labels_reverse[label].size(); l++)
                                {
                                    if(region_labels_reverse[label][l]==old_label+1)
                                        region_labels_reverse[label].erase(region_labels_reverse[label].begin()+l);
                                }

                                region_labels_reverse[added_label].push_back(old_label+1);

                                saved=false;
                            }
                        }
                    }

                    else if (border==2)//grain to bubble
                    {
                        int label=region_labels[ws_region_image(x_mouse+posx,y_mouse+posy)-1]-1;

                        //HERE WE HAVE TO CHECK IF THE POINT (X,Y) IS A POINT OF A GRAIN AREA
                        if (grain_arc_index[label].size()>0)
                        {
                            //when the grain has inside arcs they have to be removed
                            if (region_labels_reverse[label].size()>1)
                            {
                                //find a minimal region that completly contains all pixels with label
                                int min_x=dim_x;
                                int min_y=dim_y;
                                int max_x=0;
                                int max_y=0;

                                //loop over all pixels with label
                                for (int pixel=0;pixel<new_areas[label].size();pixel++)
                                {
                                    point p=new_areas[label][pixel];
                                    if (p.x<min_x) min_x=p.x;
                                    if (p.y<min_y) min_y=p.y;
                                    if (p.x>max_x) max_x=p.x;
                                    if (p.y>max_y) max_y=p.y;
                                }

                                vigra::BasicImage<bool> temp_img(max_x-min_x+3,max_y-min_y+3);

                                for(int y=0;y<temp_img.height();y++)
                                {
                                    for (int x=0;x<temp_img.width();x++)
                                    {
                                        temp_img(x,y)=true;
                                    }
                                }

                                //loop over all pixels with label
                                for (int pixel=0;pixel<new_areas[label].size();pixel++)
                                {
                                    point p=new_areas[label][pixel];
                                    temp_img(p.x-min_x+1,p.y-min_y+1)=false;
                                }

                                vigra::IImage labels_check(temp_img.width(),temp_img.height());
                                int nr_labels;

                                // find connected regions
                                nr_labels=vigra::labelImageWithBackground(srcImageRange(temp_img), destImage(labels_check), true, 0);

                                if (nr_labels>1)
                                {
                                    std::list<int> inside_labels;

                                    for(int y=0;y<temp_img.height();y++)
                                    {
                                        for (int x=0;x<temp_img.width();x++)
                                        {
                                            if (labels_check(x,y)>0 && labels_check(x,y)!=labels_check(0,0))
                                                inside_labels.push_back(region_labels[ws_region_image(x+min_x-1,y+min_y-1)-1]-1);
                                        }
                                    }

                                    inside_labels.sort();
                                    inside_labels.unique();

                                    for (std::list<int>::iterator it=inside_labels.begin(); it!=inside_labels.end(); ++it)
                                    {
                                        std::cout << "found inside label "<< *it+1 << std::endl;

                                        //if inside label is bubble remove bubble
                                        if (bubble_area_size[*it]>0)
                                        {
                                            //loop over all pixels with inside label
                                            for (int pixel=0;pixel<new_areas[*it].size();pixel++)
                                            {
                                                int x=new_areas[*it][pixel].x;
                                                int y=new_areas[*it][pixel].y;

                                                image(x,y,0,0)=original_image(x,y);
                                                image(x,y,0,1)=original_image(x,y);
                                                image(x,y,0,2)=original_image(x,y);
                                            }

                                            //look for arcs of these two labels -> no longer bubble arcs
                                            for(int y=0;y<(int)two_boundings.shape(0);y++)
                                            {
                                                if ((region_labels[two_boundings(y,0)-1]-1==label &&
                                                        region_labels[two_boundings(y,1)-1]-1==*it) ||
                                                    (region_labels[two_boundings(y,1)-1]-1==label &&
                                                        region_labels[two_boundings(y,0)-1]-1==*it))
                                                {
                                                    arc_class[y]=0;
                                                    std::cout<<"remove bubble arc "<<y<<std::endl;
                                                }
                                            }

                                            //new grain arcs
                                            std::list<int> new_grain_arcs;

                                            for (int a=0; a<grain_arc_index[label].size(); a++)
                                                new_grain_arcs.push_back(grain_arc_index[label][a]);

                                            //replace arc list
                                            grain_arc_index[label].clear();
                                            for (std::list<int>::iterator a=new_grain_arcs.begin(); a!=new_grain_arcs.end();++a) 
                                            {
                                                if (arc_class[*a]>0) grain_arc_index[label].push_back(*a);//check if arc is bubble or grain arc
                                            }
                                        
                                            std::cout<<"label "<<*it+1<<" will be removed and merged with label "<<label+1<<std::endl;

                                            //erase *it -> decrease all labels higher than *it by one
                                            grain_arc_index.erase(grain_arc_index.begin()+*it);
                                            bubble_arc_index.erase(bubble_arc_index.begin()+*it);

                                            for(int l=0; l<region_labels.size(); l++)
                                                if (region_labels[l]>*it) region_labels[l]--;

                                            //update areas and area size
                                            for (int p=0; p<new_areas[*it].size(); p++)
                                                new_areas[label].push_back(new_areas[*it][p]);
                                            grain_area_size[label]=new_areas[label].size();

                                            //"erase" doesn't free memory, so swaping is necessary
                                            std::vector< std::vector<point> > areas2;

                                            for (int area=*it+1; area<new_areas.size(); area++)
                                                areas2.push_back(new_areas[area]);

                                            new_areas.resize(*it);

                                            for (int area=0; area<areas2.size(); area++)
                                                new_areas.push_back(areas2[area]);

                                            areas2.clear();

                                            grain_area_center_mass.erase(grain_area_center_mass.begin()+*it);
                                            bubble_area_center_mass.erase(bubble_area_center_mass.begin()+*it);
                                            found_bubble_areas.clear();
                                            nr_new_areas--;

                                            delete grain_area_size;
                                            delete bubble_area_size;
                                            grain_area_size = new long[nr_new_areas];
                                            bubble_area_size = new long[nr_new_areas];

                                            for (int area=0; area<nr_new_areas; area++)
                                            {
                                                if (grain_arc_index[area].size()>0) grain_area_size[area]=new_areas[area].size();
                                                else grain_area_size[area]=0;

                                                if (bubble_arc_index[area].size()>0)
                                                {
                                                    bubble_area_size[area]=new_areas[area].size();
                                                    found_bubble_areas.push_back(area);
                                                }
                                                else bubble_area_size[area]=0;
                                            }

                                            //calculate new center of mass of merged grain
                                            long grain_area_x_sum=0;
                                            long grain_area_y_sum=0;

                                            //loop over all pixels with label
                                            for (int pixel=0;pixel<new_areas[label].size();pixel++)
                                            {
                                                point p=new_areas[label][pixel];
                                                grain_area_x_sum+=p.x;
                                                grain_area_y_sum+=p.y;
                                            }

                                            int xx=grain_area_x_sum/grain_area_size[label];
                                            int yy=grain_area_y_sum/grain_area_size[label];
                                            grain_area_center_mass[label].x=xx;
                                            grain_area_center_mass[label].y=yy;

                                            //update region_labels and region_labels_reverse
                                            for(int l=0; l<region_labels_reverse[*it].size(); l++)
                                            {
                                                region_labels_reverse[label].push_back(region_labels_reverse[*it][l]);
                                                region_labels[region_labels_reverse[*it][l]-1]=label+1;
                                            }

                                            region_labels_reverse.erase(region_labels_reverse.begin()+*it);
                                        }
                                        else //remove grain
                                        {
                                            //look for arcs of these two labels -> no longer grain arcs
                                            for(int y=0;y<(int)two_boundings.shape(0);y++)
                                            {
                                                if ((region_labels[two_boundings(y,0)-1]-1==label &&
                                                        region_labels[two_boundings(y,1)-1]-1==*it) ||
                                                    (region_labels[two_boundings(y,1)-1]-1==label &&
                                                        region_labels[two_boundings(y,0)-1]-1==*it))
                                                {
                                                    arc_class[y]=0;
                                                    grain_arc[y]=false;
                                                    std::cout<<"remove grain arc "<<y<<std::endl;
                                                }
                                            }

                                            //new grain arcs
                                            std::list<int> new_grain_arcs;

                                            for (int a=0; a<grain_arc_index[label].size(); a++)
                                                new_grain_arcs.push_back(grain_arc_index[label][a]);

                                            //replace arc list
                                            grain_arc_index[label].clear();
                                            for (std::list<int>::iterator a=new_grain_arcs.begin(); a!=new_grain_arcs.end();++a) 
                                            {
                                                if (arc_class[*a]>0) grain_arc_index[label].push_back(*a);//check if arc is bubble or grain arc
                                            }
                                        
                                            std::cout<<"label "<<*it+1<<" will be removed and merged with label "<<label+1<<std::endl;

                                            //erase *it -> decrease all labels higher than *it by one
                                            grain_arc_index.erase(grain_arc_index.begin()+*it);
                                            bubble_arc_index.erase(bubble_arc_index.begin()+*it);

                                            for(int l=0; l<region_labels.size(); l++)
                                                if (region_labels[l]>*it) region_labels[l]--;

                                            //update areas and area size
                                            for (int p=0; p<new_areas[*it].size(); p++)
                                                new_areas[label].push_back(new_areas[*it][p]);
                                            grain_area_size[label]=new_areas[label].size();

                                            //"erase" doesn't free memory, so swaping is necessary
                                            std::vector< std::vector<point> > areas2;

                                            for (int area=*it+1; area<new_areas.size(); area++)
                                                areas2.push_back(new_areas[area]);

                                            new_areas.resize(*it);

                                            for (int area=0; area<areas2.size(); area++)
                                                new_areas.push_back(areas2[area]);

                                            areas2.clear();

                                            grain_area_center_mass.erase(grain_area_center_mass.begin()+*it);
                                            bubble_area_center_mass.erase(bubble_area_center_mass.begin()+*it);
                                            found_bubble_areas.clear();
                                            nr_new_areas--;

                                            delete grain_area_size;
                                            delete bubble_area_size;
                                            grain_area_size = new long[nr_new_areas];
                                            bubble_area_size = new long[nr_new_areas];

                                            for (int area=0; area<nr_new_areas; area++)
                                            {
                                                if (grain_arc_index[area].size()>0) grain_area_size[area]=new_areas[area].size();
                                                else grain_area_size[area]=0;

                                                if (bubble_arc_index[area].size()>0)
                                                {
                                                    bubble_area_size[area]=new_areas[area].size();
                                                    found_bubble_areas.push_back(area);
                                                }
                                                else bubble_area_size[area]=0;
                                            }

                                            //calculate new center of mass of merged grain
                                            long grain_area_x_sum=0;
                                            long grain_area_y_sum=0;

                                            //loop over all pixels with label
                                            for (int pixel=0;pixel<new_areas[label].size();pixel++)
                                            {
                                                point p=new_areas[label][pixel];
                                                grain_area_x_sum+=p.x;
                                                grain_area_y_sum+=p.y;
                                            }

                                            int xx=grain_area_x_sum/grain_area_size[label];
                                            int yy=grain_area_y_sum/grain_area_size[label];
                                            grain_area_center_mass[label].x=xx;
                                            grain_area_center_mass[label].y=yy;

                                            //update region_labels and region_labels_reverse
                                            for(int l=0; l<region_labels_reverse[*it].size(); l++)
                                            {
                                                region_labels_reverse[label].push_back(region_labels_reverse[*it][l]);
                                                region_labels[region_labels_reverse[*it][l]-1]=label+1;
                                            }

                                            region_labels_reverse.erase(region_labels_reverse.begin()+*it);
                                        }

                                        //update inside labels
                                        for (std::list<int>::iterator it2=inside_labels.begin(); it2!=inside_labels.end(); ++it2)
                                        {
                                            *it2 = *it2-1; //correct as entries are sorted
                                        }
                                    }
                                }
                            }

                            //grain arcs become bubble arcs, bubble arcs are unchanged
                            for(int arc=0; arc<grain_arc_index[label].size(); arc++)
                            {
                                int y=grain_arc_index[label][arc];
                                if (grain_arc[y])
                                {
                                    arc_class[y]=2;
                                    grain_arc[y]=false;
                                    std::cout<<"grain arc to bubble arc "<<y<<std::endl;
                                }
                                else std::cout<<"unchange bubble arc "<<y<<std::endl;
                                bubble_arc_index[label].push_back(y);

                            }

                            grain_arc_index[label].clear();
                            found_bubble_areas.push_back(label);

                            bubble_area_size[label]=grain_area_size[label];
                            grain_area_size[label]=0;

                            bubble_area_center_mass[label].x=grain_area_center_mass[label].x;
                            bubble_area_center_mass[label].y=grain_area_center_mass[label].x;

                            grain_area_center_mass[label].x=0;
                            grain_area_center_mass[label].y=0;

                            saved=false;
                        }
                    }

                    else if (border==3)//bubble to grain
                    {
                        int label=region_labels[ws_region_image(x_mouse+posx,y_mouse+posy)-1]-1;

                        //HERE WE HAVE TO CHECK IF THE POINT (X,Y) IS A POINT OF A BUBBLE AREA
                        if (bubble_arc_index[label].size()>0)
                        {
                            //when the bubble has inside arcs they have to be removed
                            if (region_labels_reverse[label].size()>1)
                            {
                                //find a minimal region that completly contains all pixels with label
                                int min_x=dim_x;
                                int min_y=dim_y;
                                int max_x=0;
                                int max_y=0;

                                //loop over all pixels with label
                                for (int pixel=0;pixel<new_areas[label].size();pixel++)
                                {
                                    point p=new_areas[label][pixel];
                                    if (p.x<min_x) min_x=p.x;
                                    if (p.y<min_y) min_y=p.y;
                                    if (p.x>max_x) max_x=p.x;
                                    if (p.y>max_y) max_y=p.y;
                                }

                                vigra::BasicImage<bool> temp_img(max_x-min_x+3,max_y-min_y+3);

                                for(int y=0;y<temp_img.height();y++)
                                {
                                    for (int x=0;x<temp_img.width();x++)
                                    {
                                        temp_img(x,y)=true;
                                    }
                                }

                                //loop over all pixels with label
                                for (int pixel=0;pixel<new_areas[label].size();pixel++)
                                {
                                    point p=new_areas[label][pixel];
                                    temp_img(p.x-min_x+1,p.y-min_y+1)=false;
                                }

                                vigra::IImage labels_check(temp_img.width(),temp_img.height());
                                int nr_labels;

                                // find connected regions
                                nr_labels=vigra::labelImageWithBackground(srcImageRange(temp_img), destImage(labels_check), true, 0);

                                if (nr_labels>1)
                                {
                                    std::list<int> inside_labels;

                                    for(int y=0;y<temp_img.height();y++)
                                    {
                                        for (int x=0;x<temp_img.width();x++)
                                        {
                                            if (labels_check(x,y)>0 && labels_check(x,y)!=labels_check(0,0))
                                                inside_labels.push_back(region_labels[ws_region_image(x+min_x-1,y+min_y-1)-1]-1);
                                        }
                                    }

                                    inside_labels.sort();
                                    inside_labels.unique();

                                    for (std::list<int>::iterator it=inside_labels.begin(); it!=inside_labels.end(); ++it)
                                    {
                                        std::cout << "found inside label "<< *it+1 << std::endl;

                                        //if inside label is bubble remove bubble
                                        if (bubble_area_size[*it]>0)
                                        {
                                            //loop over all pixels with inside label
                                            for (int pixel=0;pixel<new_areas[*it].size();pixel++)
                                            {
                                                int x=new_areas[*it][pixel].x;
                                                int y=new_areas[*it][pixel].y;

                                                image(x,y,0,0)=original_image(x,y);
                                                image(x,y,0,1)=original_image(x,y);
                                                image(x,y,0,2)=original_image(x,y);
                                            }

                                            //look for arcs of these two labels -> no longer bubble arcs
                                            for(int y=0;y<(int)two_boundings.shape(0);y++)
                                            {
                                                if ((region_labels[two_boundings(y,0)-1]-1==label &&
                                                        region_labels[two_boundings(y,1)-1]-1==*it) ||
                                                    (region_labels[two_boundings(y,1)-1]-1==label &&
                                                        region_labels[two_boundings(y,0)-1]-1==*it))
                                                {
                                                    arc_class[y]=0;
                                                    std::cout<<"remove bubble arc "<<y<<std::endl;
                                                }
                                            }

                                            //new bubble arcs
                                            std::list<int> new_bubble_arcs;

                                            for (int a=0; a<bubble_arc_index[label].size(); a++)
                                                new_bubble_arcs.push_back(bubble_arc_index[label][a]);

                                            //replace arc list
                                            bubble_arc_index[label].clear();
                                            for (std::list<int>::iterator a=new_bubble_arcs.begin(); a!=new_bubble_arcs.end();++a) 
                                            {
                                                if (arc_class[*a]>0) bubble_arc_index[label].push_back(*a);//check if arc is bubble arc
                                            }
                                        
                                            std::cout<<"label "<<*it+1<<" will be removed and merged with label "<<label+1<<std::endl;

                                            //erase *it -> decrease all labels higher than *it by one
                                            grain_arc_index.erase(grain_arc_index.begin()+*it);
                                            bubble_arc_index.erase(bubble_arc_index.begin()+*it);

                                            for(int l=0; l<region_labels.size(); l++)
                                                if (region_labels[l]>*it) region_labels[l]--;

                                            //update areas and area size
                                            for (int p=0; p<new_areas[*it].size(); p++)
                                                new_areas[label].push_back(new_areas[*it][p]);
                                            bubble_area_size[label]=new_areas[label].size();

                                            //"erase" doesn't free memory, so swaping is necessary
                                            std::vector< std::vector<point> > areas2;

                                            for (int area=*it+1; area<new_areas.size(); area++)
                                                areas2.push_back(new_areas[area]);

                                            new_areas.resize(*it);

                                            for (int area=0; area<areas2.size(); area++)
                                                new_areas.push_back(areas2[area]);

                                            areas2.clear();

                                            grain_area_center_mass.erase(grain_area_center_mass.begin()+*it);
                                            bubble_area_center_mass.erase(bubble_area_center_mass.begin()+*it);
                                            found_bubble_areas.clear();
                                            nr_new_areas--;

                                            delete grain_area_size;
                                            delete bubble_area_size;
                                            grain_area_size = new long[nr_new_areas];
                                            bubble_area_size = new long[nr_new_areas];

                                            for (int area=0; area<nr_new_areas; area++)
                                            {
                                                if (grain_arc_index[area].size()>0) grain_area_size[area]=new_areas[area].size();
                                                else grain_area_size[area]=0;

                                                if (bubble_arc_index[area].size()>0)
                                                {
                                                    bubble_area_size[area]=new_areas[area].size();
                                                    found_bubble_areas.push_back(area);
                                                }
                                                else bubble_area_size[area]=0;
                                            }

                                            //calculate new center of mass of merged bubble
                                            long bubble_area_x_sum=0;
                                            long bubble_area_y_sum=0;

                                            //loop over all pixels with label
                                            for (int pixel=0;pixel<new_areas[label].size();pixel++)
                                            {
                                                point p=new_areas[label][pixel];
                                                bubble_area_x_sum+=p.x;
                                                bubble_area_y_sum+=p.y;
                                            }

                                            int xx=bubble_area_x_sum/bubble_area_size[label];
                                            int yy=bubble_area_y_sum/bubble_area_size[label];
                                            bubble_area_center_mass[label].x=xx;
                                            bubble_area_center_mass[label].y=yy;

                                            //update region_labels and region_labels_reverse
                                            for(int l=0; l<region_labels_reverse[*it].size(); l++)
                                            {
                                                region_labels_reverse[label].push_back(region_labels_reverse[*it][l]);
                                                region_labels[region_labels_reverse[*it][l]-1]=label+1;
                                            }

                                            region_labels_reverse.erase(region_labels_reverse.begin()+*it);
                                        }
                                        else //remove grain
                                        {
                                            //look for arcs of these two labels -> no longer grain arcs
                                            for(int y=0;y<(int)two_boundings.shape(0);y++)
                                            {
                                                if ((region_labels[two_boundings(y,0)-1]-1==label &&
                                                        region_labels[two_boundings(y,1)-1]-1==*it) ||
                                                    (region_labels[two_boundings(y,1)-1]-1==label &&
                                                        region_labels[two_boundings(y,0)-1]-1==*it))
                                                {
                                                    arc_class[y]=0;
                                                    grain_arc[y]=false;
                                                    std::cout<<"remove grain arc "<<y<<std::endl;
                                                }
                                            }

                                            //new bubble arcs
                                            std::list<int> new_bubble_arcs;

                                            for (int a=0; a<bubble_arc_index[label].size(); a++)
                                                new_bubble_arcs.push_back(bubble_arc_index[label][a]);

                                            //replace arc list
                                            bubble_arc_index[label].clear();
                                            for (std::list<int>::iterator a=new_bubble_arcs.begin(); a!=new_bubble_arcs.end();++a) 
                                            {
                                                if (arc_class[*a]>0) bubble_arc_index[label].push_back(*a);//check if arc is bubble arc
                                            }
                                        
                                            std::cout<<"label "<<*it+1<<" will be removed and merged with label "<<label+1<<std::endl;

                                            //erase *it -> decrease all labels higher than *it by one
                                            grain_arc_index.erase(grain_arc_index.begin()+*it);
                                            bubble_arc_index.erase(bubble_arc_index.begin()+*it);

                                            for(int l=0; l<region_labels.size(); l++)
                                                if (region_labels[l]>*it) region_labels[l]--;

                                            //update areas and area size
                                            for (int p=0; p<new_areas[*it].size(); p++)
                                                new_areas[label].push_back(new_areas[*it][p]);
                                            bubble_area_size[label]=new_areas[label].size();

                                            //"erase" doesn't free memory, so swaping is necessary
                                            std::vector< std::vector<point> > areas2;

                                            for (int area=*it+1; area<new_areas.size(); area++)
                                                areas2.push_back(new_areas[area]);

                                            new_areas.resize(*it);

                                            for (int area=0; area<areas2.size(); area++)
                                                new_areas.push_back(areas2[area]);

                                            areas2.clear();

                                            grain_area_center_mass.erase(grain_area_center_mass.begin()+*it);
                                            bubble_area_center_mass.erase(bubble_area_center_mass.begin()+*it);
                                            found_bubble_areas.clear();
                                            nr_new_areas--;

                                            delete grain_area_size;
                                            delete bubble_area_size;
                                            grain_area_size = new long[nr_new_areas];
                                            bubble_area_size = new long[nr_new_areas];

                                            for (int area=0; area<nr_new_areas; area++)
                                            {
                                                if (grain_arc_index[area].size()>0) grain_area_size[area]=new_areas[area].size();
                                                else grain_area_size[area]=0;

                                                if (bubble_arc_index[area].size()>0)
                                                {
                                                    bubble_area_size[area]=new_areas[area].size();
                                                    found_bubble_areas.push_back(area);
                                                }
                                                else bubble_area_size[area]=0;
                                            }

                                            //calculate new center of mass of merged bubble
                                            long bubble_area_x_sum=0;
                                            long bubble_area_y_sum=0;

                                            //loop over all pixels with label
                                            for (int pixel=0;pixel<new_areas[label].size();pixel++)
                                            {
                                                point p=new_areas[label][pixel];
                                                bubble_area_x_sum+=p.x;
                                                bubble_area_y_sum+=p.y;
                                            }

                                            int xx=bubble_area_x_sum/bubble_area_size[label];
                                            int yy=bubble_area_y_sum/bubble_area_size[label];
                                            bubble_area_center_mass[label].x=xx;
                                            bubble_area_center_mass[label].y=yy;

                                            //update region_labels and region_labels_reverse
                                            for(int l=0; l<region_labels_reverse[*it].size(); l++)
                                            {
                                                region_labels_reverse[label].push_back(region_labels_reverse[*it][l]);
                                                region_labels[region_labels_reverse[*it][l]-1]=label+1;
                                            }

                                            region_labels_reverse.erase(region_labels_reverse.begin()+*it);
                                        }

                                        //update inside labels
                                        for (std::list<int>::iterator it2=inside_labels.begin(); it2!=inside_labels.end(); ++it2)
                                        {
                                            *it2 = *it2-1; //correct as entries are sorted
                                        }
                                    }
                                }
                            }

                            //bubble arcs become grain arcs
                            for(int arc=0; arc<bubble_arc_index[label].size(); arc++)
                            {
                                int y=bubble_arc_index[label][arc];

                                bool other_bubble_arc=false;

                                for (int bubble=0; bubble<found_bubble_areas.size() && !other_bubble_arc; bubble++)
                                {
                                    int other_label=found_bubble_areas[bubble];

                                    if (other_label!=label)
                                    {
                                        for (int arc=0; arc<bubble_arc_index[other_label].size() && !other_bubble_arc; arc++)
                                        {
                                            if (bubble_arc_index[other_label][arc]==y) other_bubble_arc=true;
                                        }
                                    }
                                }

                                if (!other_bubble_arc)
                                {
                                    arc_class[y]=1;
                                    grain_arc[y]=true;
                                    std::cout<<"bubble arc to grain arc "<<y<<std::endl;
                                }
                                else std::cout<<"unchange bubble arc "<<y<<std::endl;

                                grain_arc_index[label].push_back(y);
                            }

                            bubble_arc_index[label].clear();

                            grain_area_size[label]=bubble_area_size[label];
                            bubble_area_size[label]=0;

                            grain_area_center_mass[label].x=bubble_area_center_mass[label].x;
                            grain_area_center_mass[label].y=bubble_area_center_mass[label].x;

                            bubble_area_center_mass[label].x=0;
                            bubble_area_center_mass[label].y=0;

                            found_bubble_areas.clear();

                            for (int area=0; area<nr_new_areas; area++)
                            {
                                if (bubble_arc_index[area].size()>0)
                                {
                                       found_bubble_areas.push_back(area);
                                }
                            }

                            //loop over all pixels with label
                            for (int pixel=0;pixel<new_areas[label].size();pixel++)
                            {
                                int x=new_areas[label][pixel].x;
                                int y=new_areas[label][pixel].y;

                                image(x,y,0,0)=original_image(x,y);
                                image(x,y,0,1)=original_image(x,y);
                                image(x,y,0,2)=original_image(x,y);
                            }

                            saved=false;
                        }
                    }

                    /* Boundary to subgrain boundary. Procedure:
                     * (1) If the clicked arc is part of a grain boundary, the same procedure as for border = 0 applies:
                     *     Merge if there are two grain boundaries to connect, or destroy if the grain boundary would connect to the outside.
                     *     Since the merging in border = 0 works just fine, and my version does not, all of the merging code has been copy-and-pasted into this version.
                     * (2) If the clicked arc is not part of a grain boundary, make the boundary a subgrain boundary.
                     */
                    else if(border == 4) //**********************************************************************************************************************************************+
                    {
                        int second_x = main_disp.mouse_x;
                        int second_y = main_disp.mouse_y;

                        while(main_disp.button)
                        {
                            if(main_disp.mouse_y>=0 && main_disp.mouse_x>=0)
                            {
                                second_x=main_disp.mouse_x;
                                second_y=main_disp.mouse_y;

                                cimg_library::CImg<unsigned char> temp_selected_image=selected_image;

                                for(int x=std::min(x_mouse,second_x); x<=std::max(x_mouse,second_x); x++)
                                {
                                    temp_selected_image(x,y_mouse,0,0)=color[0][0];
                                    temp_selected_image(x,y_mouse,0,1)=color[0][1];
                                    temp_selected_image(x,y_mouse,0,2)=color[0][2];
                                }
                                for(int x=std::min(x_mouse,second_x); x<=std::max(x_mouse,second_x); x++)
                                {
                                    temp_selected_image(x,second_y,0,0)=color[0][0];
                                    temp_selected_image(x,second_y,0,1)=color[0][1];
                                    temp_selected_image(x,second_y,0,2)=color[0][2];
                                }
                                for(int y=std::min(y_mouse,second_y); y<=std::max(y_mouse,second_y); y++)
                                {
                                    temp_selected_image(x_mouse,y,0,0)=color[0][0];
                                    temp_selected_image(x_mouse,y,0,1)=color[0][1];
                                    temp_selected_image(x_mouse,y,0,2)=color[0][2];
                                }
                                for(int y=std::min(y_mouse,second_y); y<=std::max(y_mouse,second_y); y++)
                                {
                                    temp_selected_image(second_x,y,0,0)=color[0][0];
                                    temp_selected_image(second_x,y,0,1)=color[0][1];
                                    temp_selected_image(second_x,y,0,2)=color[0][2];
                                }

                                main_disp.display(temp_selected_image);
                            }
                        }
                        
                        //Check if the point (x,y) is a point of an arc
                        for(int arc_count = 0;  arc_count < (int)arcs.size(); arc_count++)
                        {
                            bool found_point = false;

                            for(int point_count = 0; point_count < (int)arcs[arc_count].size() && found_point == false; point_count++)
                            {
                                 if(arcs[arc_count][point_count].x>=std::min(x_mouse,second_x)+posx &&
                                    arcs[arc_count][point_count].y>=std::min(y_mouse,second_y)+posy &&
                                    arcs[arc_count][point_count].x<=std::max(x_mouse,second_x)+posx &&
                                    arcs[arc_count][point_count].y<=std::max(y_mouse,second_y)+posy)
                                {
                                    found_point = true;

                                    if(arc_class[arc_count] == 1) //Grain arc
                                    {
                                        bool outside=false;
                                        int outside_label;
                                        bool bubble_merged=false;

                                        //these two labels are selected for merging
                                        int label_merged = std::min(region_labels[two_boundings(arc_count,0)-1]-1,
                                                                    region_labels[two_boundings(arc_count,1)-1]-1);
                                        int label_removed = std::max(region_labels[two_boundings(arc_count,0)-1]-1,
                                                                     region_labels[two_boundings(arc_count,1)-1]-1);

                                        if (bubble_area_size[label_removed]>0 || bubble_area_size[label_merged]>0) bubble_merged=true;

                                        //look for arcs of these two labels -> no longer grain arcs
                                        for(int y=0;y<(int)two_boundings.shape(0);y++)
                                        {
                                            if ((region_labels[two_boundings(y,0)-1]-1==label_merged &&
                                                    region_labels[two_boundings(y,1)-1]-1==label_removed) ||
                                                (region_labels[two_boundings(y,1)-1]-1==label_merged &&
                                                    region_labels[two_boundings(y,0)-1]-1==label_removed))
                                            {
                                                arc_class[y]=0;
                                                grain_arc[y]=false;
                                                std::cout<<"remove grain arc "<<y<<std::endl;
                                            }
                                        }

                                        //grain belongs now to outside -> remove as grain
                                        if(bubble_area_size[label_merged]==grain_area_size[label_merged])
                                        {
                                            outside=true;
                                            grain_arc_index[label_merged].clear();

                                            outside_label=label_merged+1;

                                            std::cout<<"label "<<label_removed+1<<" will be removed and merged with outside label "<<outside_label
                                                <<std::endl;
                                        }
                                        else if(bubble_area_size[label_removed]==grain_area_size[label_removed])
                                        {
                                            outside=true;
                                            grain_arc_index[label_removed].clear();

                                            outside_label=label_removed+1;
                                            label_removed=label_merged;

                                            std::cout<<"label "<<label_removed+1<<" will be removed and merged with outside label "<<outside_label
                                                <<std::endl;
                                            outside_label--;
                                        }
                                        else
                                        {
                                            //combine arc lists
                                            std::list<int> combined_arc_list;

                                            if (grain_arc_index[label_merged].size()>0)
                                                for (int a=0; a<grain_arc_index[label_merged].size(); a++)
                                                    combined_arc_list.push_back(grain_arc_index[label_merged][a]);
                                            else for (int a=0; a<bubble_arc_index[label_merged].size(); a++)
                                                    combined_arc_list.push_back(bubble_arc_index[label_merged][a]);

                                            if (grain_arc_index[label_removed].size()>0)
                                                for (int a=0; a<grain_arc_index[label_removed].size(); a++)
                                                    combined_arc_list.push_back(grain_arc_index[label_removed][a]);
                                            else for (int a=0; a<bubble_arc_index[label_removed].size(); a++)
                                                    combined_arc_list.push_back(bubble_arc_index[label_removed][a]);

                                            //sort combined list und remove double entries
                                            combined_arc_list.sort();
                                            combined_arc_list.unique();

                                            //replace arc list of label_merged
                                            if (!bubble_merged)
                                            {
                                                grain_arc_index[label_merged].clear();
                            
                                                for (std::list<int>::iterator a=combined_arc_list.begin(); a!=combined_arc_list.end();++a) 
                                                {
                                                    if (arc_class[*a]>0) grain_arc_index[label_merged].push_back(*a);//check if bubble or grain arc
                                                }
                                            }
                                            else
                                            {
                                                bubble_arc_index[label_merged].clear();
                            
                                                for (std::list<int>::iterator a=combined_arc_list.begin(); a!=combined_arc_list.end();++a) 
                                                {
                                                    if (arc_class[*a]>0) bubble_arc_index[label_merged].push_back(*a);//check if bubble or grain arc
                                                }
                                            }

                                            std::cout<<"label "<<label_removed+1<<" will be removed and merged with label "<<label_merged+1<<std::endl;
                                        }

                                        //erase label_removed -> decrease all labels higher than label_removed by one
                                        grain_arc_index.erase(grain_arc_index.begin()+label_removed);
                                        bubble_arc_index.erase(bubble_arc_index.begin()+label_removed);

                                        for(int l=0; l<region_labels.size(); l++)
                                            if (region_labels[l]>label_removed) region_labels[l]--;

                                        //update areas and area size
                                        if (!outside)
                                        {
                                            for (int p=0; p<new_areas[label_removed].size(); p++)
                                                new_areas[label_merged].push_back(new_areas[label_removed][p]);
                                            grain_area_size[label_merged]=new_areas[label_merged].size();
                                        }

                                        //new_areas.erase(new_areas.begin()+label_removed);
                                        //"erase" doesn't free memory, so swaping is necessary
                                        std::vector< std::vector<point> > areas2;

                                        for (int area=label_removed+1; area<new_areas.size(); area++)
                                            areas2.push_back(new_areas[area]);

                                        new_areas.resize(label_removed);

                                        for (int area=0; area<areas2.size(); area++)
                                            new_areas.push_back(areas2[area]);

                                        areas2.clear();

                                        grain_area_center_mass.erase(grain_area_center_mass.begin()+label_removed);
                                        bubble_area_center_mass.erase(bubble_area_center_mass.begin()+label_removed);
                                        found_bubble_areas.clear();
                                        nr_new_areas--;

                                        delete grain_area_size;
                                        delete bubble_area_size;
                                        grain_area_size = new long[nr_new_areas];
                                        bubble_area_size = new long[nr_new_areas];

                                        for (int area=0; area<nr_new_areas; area++)
                                        {
                                            if (grain_arc_index[area].size()>0) grain_area_size[area]=new_areas[area].size();
                                            else grain_area_size[area]=0;

                                            if (bubble_arc_index[area].size()>0)
                                            {
                                                bubble_area_size[area]=new_areas[area].size();
                                                found_bubble_areas.push_back(area);
                                            }
                                            else bubble_area_size[area]=0;
                                        }

                                        if (!outside)
                                        {
                                            //calculate new center of mass of merged grain/bubble
                                            long area_x_sum=0;
                                            long area_y_sum=0;

                                            //loop over all pixels with label
                                            for (int pixel=0;pixel<new_areas[label_merged].size();pixel++)
                                            {
                                                point p=new_areas[label_merged][pixel];
                                                area_x_sum+=p.x;
                                                area_y_sum+=p.y;
                                            }

                                            int xx=area_x_sum/new_areas[label_merged].size();
                                            int yy=area_y_sum/new_areas[label_merged].size();

                                            if (!bubble_merged)
                                            {
                                                grain_area_center_mass[label_merged].x=xx;
                                                grain_area_center_mass[label_merged].y=yy;
                                            }
                                            else
                                            {
                                                bubble_area_center_mass[label_merged].x=xx;
                                                bubble_area_center_mass[label_merged].y=yy;
                                            }

                                            //update region_labels and region_labels_reverse
                                            for(int l=0; l<region_labels_reverse[label_removed].size(); l++)
                                            {
                                                region_labels_reverse[label_merged].push_back(region_labels_reverse[label_removed][l]);
                                                region_labels[region_labels_reverse[label_removed][l]-1]=label_merged+1;
                                            }

                                            region_labels_reverse.erase(region_labels_reverse.begin()+label_removed);
                                        }
                                        else
                                        {
                                            //update region_labels and region_labels_reverse
                                            for(int l=0; l<region_labels_reverse[label_removed].size(); l++)
                                            {
                                                region_labels[region_labels_reverse[label_removed][l]-1]=outside_label;
                                            }

                                            region_labels_reverse.erase(region_labels_reverse.begin()+label_removed);
                                        }

                                        //Finally, make the grain arc a subgrain arc
                                        subgrain[arc_count] = true;
                                        grain_arc[arc_count] = false;
                                        
                                        //Change the arc's color appropriately
                                        arc_class[arc_count] = 5;

                                        saved=false;                                                                            
                                    }

                                    if(arc_class[arc_count] == 0) //No grain arc
                                    {                                
                                        //Make the grain arc a subgrain arc
                                        subgrain[arc_count] = true;
                                        grain_arc[arc_count] = false;
                                        
                                        //Change the arc's color appropriately
                                        arc_class[arc_count] = 5;
                                    }
                                }
                            }
                        }
                    }

                    /* Subgrain boundary to no boundary
                     */
                    else if (border==5)
                    {
                        int second_x = main_disp.mouse_x;
                        int second_y = main_disp.mouse_y;

                        while(main_disp.button)
                        {
                            if(main_disp.mouse_y>=0 && main_disp.mouse_x>=0)
                            {
                                second_x=main_disp.mouse_x;
                                second_y=main_disp.mouse_y;

                                cimg_library::CImg<unsigned char> temp_selected_image=selected_image;

                                for(int x=std::min(x_mouse,second_x); x<=std::max(x_mouse,second_x); x++)
                                {
                                    temp_selected_image(x,y_mouse,0,0)=color[0][0];
                                    temp_selected_image(x,y_mouse,0,1)=color[0][1];
                                    temp_selected_image(x,y_mouse,0,2)=color[0][2];
                                }
                                for(int x=std::min(x_mouse,second_x); x<=std::max(x_mouse,second_x); x++)
                                {
                                    temp_selected_image(x,second_y,0,0)=color[0][0];
                                    temp_selected_image(x,second_y,0,1)=color[0][1];
                                    temp_selected_image(x,second_y,0,2)=color[0][2];
                                }
                                for(int y=std::min(y_mouse,second_y); y<=std::max(y_mouse,second_y); y++)
                                {
                                    temp_selected_image(x_mouse,y,0,0)=color[0][0];
                                    temp_selected_image(x_mouse,y,0,1)=color[0][1];
                                    temp_selected_image(x_mouse,y,0,2)=color[0][2];
                                }
                                for(int y=std::min(y_mouse,second_y); y<=std::max(y_mouse,second_y); y++)
                                {
                                    temp_selected_image(second_x,y,0,0)=color[0][0];
                                    temp_selected_image(second_x,y,0,1)=color[0][1];
                                    temp_selected_image(second_x,y,0,2)=color[0][2];
                                }

                                main_disp.display(temp_selected_image);
                            }
                        }
                        
                        //Check if the point (x,y) is a point of an arc
                        for(int arc_count = 0;  arc_count < (int)arcs.size(); arc_count++)
                        {
                            bool found_point = false;

                            for(int point_count = 0; point_count < (int)arcs[arc_count].size() && found_point == false; point_count++)
                            {
                                 if(arcs[arc_count][point_count].x>=std::min(x_mouse,second_x)+posx &&
                                    arcs[arc_count][point_count].y>=std::min(y_mouse,second_y)+posy &&
                                    arcs[arc_count][point_count].x<=std::max(x_mouse,second_x)+posx &&
                                    arcs[arc_count][point_count].y<=std::max(y_mouse,second_y)+posy)
                                {
                                    found_point = true;

                                    if(arc_class[arc_count] == 5) //Subgrain arc
                                    {
                                        //Make the subgrain arc an arc that belongs to no class
                                        subgrain[arc_count] = false;
                                        grain_arc[arc_count] = false;
                                        
                                        //Change the arc's color appropriately
                                        arc_class[arc_count] = 0;
                                    }
                                }
                            }
                        }
                        saved = false;
                    }
                                                                  
                    for(int i=0; i<found_bubble_areas.size(); i++)//display bubble areas
                    {
                        int bubble_area=found_bubble_areas[i];
                        for (int j=0; j<new_areas[bubble_area].size(); j++)
                        {
                            int x=new_areas[bubble_area][j].x;
                            int y=new_areas[bubble_area][j].y;
                            image(x,y,0,0)=color[3][0];
                            image(x,y,0,1)=color[3][1];
                            image(x,y,0,2)=color[3][2];
                        }
                    }

                    for(int arc=0; arc<arcs.size(); arc++)//display arcs
                    {
                        for(int i=0; i<arcs[arc].size(); i++)
                        {
                            int x=arcs[arc][i].x;
                            int y=arcs[arc][i].y;
                            image(x,y,0,0)=color[arc_class[arc]][0];
                            image(x,y,0,1)=color[arc_class[arc]][1];
                            image(x,y,0,2)=color[arc_class[arc]][2];
                        }
                    }

                    cimg_forXY(selected_image, x, y)//copy to selection
                    {
                        if( x + posx >= 0 && x + posx <  image.dimx() && y + posy >= 0 && y + posy < image.dimy() )
                        {
                            selected_image(x,y,0,0) = image(x + posx,y + posy,0,0);
                            selected_image(x,y,0,1) = image(x + posx,y + posy,0,1);
                            selected_image(x,y,0,2) = image(x + posx,y + posy,0,2);
                        }
                    }

                    main_disp.display(selected_image);
                    end_of_correction=false;

                }
            }
        }
    }
}
