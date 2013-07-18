/*! \file boundary_probabilities.h
 *  \brief Boundary probability calculation.
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

#include <vector>
#include <cstring>
#include <fstream>

#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/impex.hxx>

#include <vigra/multi_array.hxx>
#include <vigra/random_forest.hxx>

#include <vigra/hdf5impex.hxx>
#include <vigra/random_forest_hdf5_impex.hxx>

#include "gbn.h"
#include "seg.h"
#include "marray.hxx"
#include <math.h>

#include "bubble_grain_structure.h"
#include "mosaic.cpp"

void save_boundary_probabilities(vigra::MultiArray<2,float> const probability,
                                 std::string filepath_to_ws_region_image,
                                 std::string path_to_output_folder,
                                 seg segment,
                                 std::string param_file_name
                                 )
{
    vigra::FImage * probability_image = new vigra::FImage[nr_of_classes];

    for(int y=0;y<segment.dim_y;y++)  
    {
        for(int x=0;x<segment.dim_x;x++)
        {
            for (int i=0; i<nr_of_classes; i++)
            {
                probability_image[i].resize(segment.dim_x,segment.dim_y);
                probability_image[i](x,y)=255;
            }
        }
    }

    //load validation training data
    int last_read_arc=-1;
    std::vector<int> vector_training_labels;
    vector_training_labels.resize((int)segment.arcs.size(),10);

    std::string validation_file="/home/akuehlwe/src/CIS/boundary-labels/";
    validation_file.append(get_filename(filepath_to_ws_region_image));
    validation_file.resize(validation_file.size()-3);
    validation_file.append(".dat");

    //string is read from temp file to check whether file is empty
    std::string teststring;

    std::ifstream training_file(validation_file.c_str());
    std::ifstream temp_training_file(validation_file.c_str());
    temp_training_file>>teststring;

    bool validation=false;
    if(training_file && teststring.size()!=0) validation=true;

    std::vector<int> class_count;
    std::vector<float> class_error;

    if(validation)
    {
        class_count.resize(nr_of_classes,0);
        class_error.resize(nr_of_classes,0.0f);

        //temp file is used to get one line out of the training file
        std::vector<int> temp;
        temp.resize(2);
        while(!training_file.eof())
        {
            //the vector of the training data is
            training_file>>temp[0];
            training_file>>temp[1];

            if(last_read_arc!=temp[0])
            {
                vector_training_labels[temp[0]]=temp[1];
                //HERE WE SET LAST READ TO THE LAST READED /WRITTEN POINTS TO AVOID THE LAST LINE
                //IS READED TWICE
                last_read_arc=temp[0];
                //we want to know how many data is from wich class
                class_count[temp[1]]++;
            }
        }
        training_file.close();
    }

    for(int a=0;a<(int)segment.arcs.size();a++)
    {
        std::vector<point>  this_arc;
        this_arc=segment.arcs[a];

        //now we loop over the points in this arc
        for(int p=0;p<(int)this_arc.size();p++)
        {
            int x=this_arc[p].x;
            int y=this_arc[p].y;
            for (int i=0; i<nr_of_classes; i++)
            {
                probability_image[i](x,y)=255*(1-probability(a,i));
            }
        }

        if(validation) for (int i=0; i<nr_of_classes; i++)
        {
            if (vector_training_labels[a]==i) class_error[i]=class_error[i]+((1-probability(a,i))/class_count[i]);
        }
    }

    if(validation)
    {
        std::string filepath_log_file="/home/akuehlwe/src/CIS/error_log.txt";
        std::ofstream log_file(filepath_log_file.c_str(), std::ios_base::out | std::ios_base::app);

        log_file <<param_file_name.c_str();

        for (int i=0; i<nr_of_classes; i++)
        {
            std::cout<<"Error Class "<<i<<": "<<class_error[i]<<std::endl;
            log_file <<" "<<class_error[i];
        }

        log_file << "\n";
        log_file.close();
    }

    std::string * filepath_prob = new std::string[nr_of_classes];

    for (int i=0; i<nr_of_classes; i++)
    {
        filepath_prob[i]=path_to_output_folder;

        std::stringstream s;
        s << i <<"/";
        filepath_prob[i].append(s.str());
        //filepath_prob[i].append(param_file_name.c_str());
        filepath_prob[i].append(get_filename(filepath_to_ws_region_image));
        filepath_prob[i].resize(filepath_prob[i].size()-6);
        filepath_prob[i].append("jpg");

        exportImage(srcImageRange(probability_image[i]), vigra::ImageExportInfo(filepath_prob[i].c_str()));
    }
}

void find_vertical_arcs(std::vector< std::vector<point> > arcs,std::vector<int> & vertical_arc_index,std::vector<bool> found_bubble_arcs,
                        vigra::BasicImage<unsigned int> ws_region_image, int dim_x, int dim_y)
{
    std::cout<<"start hough transform..."<<std::endl;

    //do a hough transform to find vertical arcs
    int max_d = 2*sqrt((dim_x*dim_x)+(dim_y*dim_y));

    size_t size[] = {5,max_d};
    marray::Marray<int> hough_space(size,size+2);

    //pixels of found lines above threshold are "hough labeled"
    vigra::BasicImage<bool> hough_labeled(dim_x,dim_y);

    //DEBUG image to show found lines
    //vigra::IRGBImage hough_result_image(dim_x,dim_y);

    for(int x=0; x<dim_x; x++)
    {
        for(int y=0; y<dim_y; y++)
        {
            hough_labeled(x,y)=false;
            //hough_result_image(x,y)[0]=2;
            //hough_result_image(x,y)[1]=2;
            //hough_result_image(x,y)[2]=2;
        }
    }

    /*
    for(int arcindex=0; arcindex<arcs.size(); arcindex++)
    {
        for (int i=0; i<arcs[arcindex].size(); i++)
        {
            hough_result_image(arcs[arcindex][i].x,arcs[arcindex][i].y)[0]=0;
            hough_result_image(arcs[arcindex][i].x,arcs[arcindex][i].y)[1]=0;
            hough_result_image(arcs[arcindex][i].x,arcs[arcindex][i].y)[2]=0;
        }
    }

    //DEBUG
    vigra::IRGBImage vertical_arc_image(dim_x,dim_y);
    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
            vertical_arc_image(x,y)[0]=1;
            vertical_arc_image(x,y)[1]=1;
            vertical_arc_image(x,y)[2]=1;
        }
    }

    for(int arc=0; arc<arcs.size(); arc++)
    {
        for (int i=0; i<arcs[arc].size(); i++)
        {
            int x=arcs[arc][i].x;
            int y=arcs[arc][i].y;

            if(x>0 && y>0) vertical_arc_image(x-1,y-1)[0]=0;
            if(x>0) vertical_arc_image(x-1,y)[0]=0;
            if(x>0 && y<dim_y-1) vertical_arc_image(x-1,y+1)[0]=0;

            if(y>0) vertical_arc_image(x,y-1)[0]=0;
            vertical_arc_image(x,y)[0]=0;
            if(y<dim_y-1) vertical_arc_image(x,y+1)[0]=0;

            if(x<dim_x-1 && y>0) vertical_arc_image(x+1,y-1)[0]=0;
            if(x<dim_x-1) vertical_arc_image(x+1,y)[0]=0;
            if(x<dim_x-1 && y<dim_y-1) vertical_arc_image(x+1,y+1)[0]=0;

            if(x>0 && y>0) vertical_arc_image(x-1,y-1)[1]=0;
            if(x>0) vertical_arc_image(x-1,y)[1]=0;
            if(x>0 && y<dim_y-1) vertical_arc_image(x-1,y+1)[1]=0;

            if(y>0) vertical_arc_image(x,y-1)[1]=0;
            vertical_arc_image(x,y)[1]=0;
            if(y<dim_y-1) vertical_arc_image(x,y+1)[1]=0;

            if(x<dim_x-1 && y>0) vertical_arc_image(x+1,y-1)[1]=0;
            if(x<dim_x-1) vertical_arc_image(x+1,y)[1]=0;
            if(x<dim_x-1 && y<dim_y-1) vertical_arc_image(x+1,y+1)[1]=0;

            if(x>0 && y>0) vertical_arc_image(x-1,y-1)[2]=0;
            if(x>0) vertical_arc_image(x-1,y)[2]=0;
            if(x>0 && y<dim_y-1) vertical_arc_image(x-1,y+1)[2]=0;

            if(y>0) vertical_arc_image(x,y-1)[2]=0;
            vertical_arc_image(x,y)[2]=0;
            if(y<dim_y-1) vertical_arc_image(x,y+1)[2]=0;

            if(x<dim_x-1 && y>0) vertical_arc_image(x+1,y-1)[2]=0;
            if(x<dim_x-1) vertical_arc_image(x+1,y)[2]=0;
            if(x<dim_x-1 && y<dim_y-1) vertical_arc_image(x+1,y+1)[2]=0;
        }
    }
    */

    //window for hough transform
    int nr_steps=10;
    for (int step=0; step<nr_steps; step++)
    {
        int low_y=step*dim_y/(nr_steps+1);
        int high_y=std::min(dim_y,(step+2)*dim_y/(nr_steps+1));

        for (int i=0; i<5; i++)
            for (int j=0; j<max_d; j++)
                hough_space(i,j)=0;

        for(int arcindex=0; arcindex<arcs.size(); arcindex++)
        {
            for (int i=0; i<arcs[arcindex].size(); i++)
            {
                int x=arcs[arcindex][i].x;
                int y=arcs[arcindex][i].y;

                if(low_y<=y && y<high_y)
                    for (int alpha=-2; alpha<3; alpha++)  //angle region between -1 degrees and +1 degrees
                    {
                          int d = 0.5 * max_d + x * cos((float)(alpha*PI)/360.0f) + y * sin((float)(alpha*PI)/360.0f);
                          hough_space(alpha+2,d)++;
                    }
            }
        }

        /*
        //to determine the angle with highest number of votes does not improve the result
        int best_alpha=0;
        int highest_alpha=0;

        for (int i=0; i<6; i++)
        {
            int alpha_sum=0;
           
            for (int j=0; j<max_d; j++)
            {
                alpha_sum+=hough_space(i,j);
            }

            if (alpha_sum>highest_alpha)
            {
                std::cout<<"best alpha "<<i-2<<" "<<alpha_sum<<std::endl;
                best_alpha=i-2;
                highest_alpha=alpha_sum;
            }
        }
        */

        for (int alpha=-2; alpha<3; alpha++)
            for (int d=0; d<max_d; d++)
            {
                if (hough_space(alpha+2,d)>2200/(nr_steps+1))//empirical threshold for probable lines
                {
                    for(int y=low_y; y<high_y; y++)
                    {
                        int x=(d-0.5*max_d)/cos((float)(alpha*PI)/360.0f) -y*tan((float)(alpha*PI)/360.0f);

                        if(x>=0 && x<dim_x)
                        {
                            //hough_result_image(x,y)[0]=1;
                            hough_labeled(x,y)=true;
                            if (x>0) hough_labeled(x-1,y)=true;
                            if (x>1) hough_labeled(x-2,y)=true;
                            if (x>2) hough_labeled(x-3,y)=true;
                            if (x>3) hough_labeled(x-4,y)=true;
                            if (x>4) hough_labeled(x-5,y)=true;
                            if (x>5) hough_labeled(x-6,y)=true;
                            if (x<dim_x-1) hough_labeled(x+1,y)=true;
                            if (x<dim_x-2) hough_labeled(x+2,y)=true;
                            if (x<dim_x-3) hough_labeled(x+3,y)=true;
                            if (x<dim_x-4) hough_labeled(x+4,y)=true;
                            if (x<dim_x-5) hough_labeled(x+5,y)=true;
                            if (x<dim_x-6) hough_labeled(x+6,y)=true;
                        }
                    }
                }
            }

        /*
        //find hough labeled grain center of mass points, how to merge them?
        //problematic to identiy arcs and areas at the same time
        int nr_areas=1;
        std::vector<long> areas_sum_x(nr_areas,0);
        std::vector<long> areas_sum_y(nr_areas,0);
        std::vector<long> areas_size(nr_areas,0);

        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
            {
                int area=ws_region_image(x,y);

                if(area>nr_areas)
                {
                    nr_areas=area;
                    areas_sum_x.resize(nr_areas,0);
                    areas_sum_y.resize(nr_areas,0);
                    areas_size.resize(nr_areas,0);
                }

                areas_sum_x[area-1]+=x;
                areas_sum_y[area-1]+=y;
                areas_size[area-1]++;
            }
        }

        std::vector<point> center_mass_points(nr_areas);

        for(int area=0; area<nr_areas; area++)
        {
            center_mass_points[area].x=areas_sum_x[area]/areas_size[area];
            center_mass_points[area].y=areas_sum_y[area]/areas_size[area];

            if(areas_size[area]<10000 && hough_labeled(center_mass_points[area].x,center_mass_points[area].y))
            {
                hough_result_image(center_mass_points[area].x,center_mass_points[area].y)[0]=0;
                hough_result_image(center_mass_points[area].x,center_mass_points[area].y)[1]=0;
                hough_result_image(center_mass_points[area].x,center_mass_points[area].y)[2]=0;
            }
        }
        */

        std::cout<<"window "<<step+1<<"/"<<nr_steps<<std::endl;

        for(int arc=0; arc<arcs.size(); arc++)
        {
            if(((low_y<=arcs[arc][0].y && arcs[arc][0].y<high_y)||(low_y<=arcs[arc].back().y && arcs[arc].back().y<high_y)) && //arc in current window
                fabs(arcs[arc][0].x-arcs[arc].back().x)<15 && //start and end point are close to vertical
                hough_labeled(arcs[arc][0].x,arcs[arc][0].y) && //start point is hough labeled
                hough_labeled(arcs[arc].back().x,arcs[arc].back().y) &&//end point is hough labeled
                !found_bubble_arcs[arc])//no bubble arc
            {
                int min=dim_x-1;
                int max=0;
                bool trend_pos=true;
                bool trend_neg=true;
                for(int pixel=0;pixel<arcs[arc].size();pixel++)
                {
                    if (arcs[arc][pixel].x<min) min=arcs[arc][pixel].x;
                    if (arcs[arc][pixel].x>max) max=arcs[arc][pixel].x;
                    if (pixel>1)
                    {
                        if (arcs[arc][pixel].x<arcs[arc][pixel-1].x) trend_pos=false;
                        if (arcs[arc][pixel].x>arcs[arc][pixel-1].x) trend_neg=false;
                    }
                }
                int y_distance=fabs(arcs[arc][0].y-arcs[arc].back().y);

                if (max-min<0.1*y_distance && max-min<15 && !trend_pos && !trend_neg) //small deviation from straight line
                {
                     //std::cout<<"vertical arc "<<arc<<" found!"<<std::endl;
                     vertical_arc_index.push_back(arc);
                }
            }
        }

        /*
        for(int found=0; found<vertical_arc_index.size(); found++)
        {
            for (int i=0; i<arcs[vertical_arc_index[found]].size(); i++)
            {
                int x=arcs[vertical_arc_index[found]][i].x;
                int y=arcs[vertical_arc_index[found]][i].y;

                if(x>0 && y>0) vertical_arc_image(x-1,y-1)[0]=1;
                if(x>0) vertical_arc_image(x-1,y)[0]=1;
                if(x>0 && y<dim_y-1) vertical_arc_image(x-1,y+1)[0]=1;

                if(y>0) vertical_arc_image(x,y-1)[0]=1;
                vertical_arc_image(x,y)[0]=1;
                if(y<dim_y-1) vertical_arc_image(x,y+1)[0]=1;

                if(x<dim_x-1 && y>0) vertical_arc_image(x+1,y-1)[0]=1;
                if(x<dim_x-1) vertical_arc_image(x+1,y)[0]=1;
                if(x<dim_x-1 && y<dim_y-1) vertical_arc_image(x+1,y+1)[0]=1;
            }
        }
        */
    }
    std::cout<<"...done"<<std::endl;
    //exportImage(srcImageRange(hough_result_image), vigra::ImageExportInfo("hough_result.jpg"));
    //exportImage(srcImageRange(vertical_arc_image), vigra::ImageExportInfo("vertical_test.png"));
}

//fn is the filepath to the image file
void extract_boundary_probabilities(std::string filepath_to_feature_file,std::string path_to_ws_image,std::string filepath_to_random_forest_file,
                                    std::string path_to_output_folder,std::string path_to_gm_output_folder,std::string param_file,ParameterFile paramFile,
                                    std::string filepath_thumbs="no", std::string filepath_image="", bool originalImageExists = true)
{
    std::string param_file_name=get_filename(param_file);
    if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
    else param_file_name = "";

    if(originalImageExists == false)
    {
        std::cout << "Not using the original image for calculations" << std::endl;
    }

    std::string filepath_to_ws_image=path_to_ws_image;
    filepath_to_ws_image.append(get_filename(filepath_to_feature_file));
    if (filepath_image!="") filepath_image.append(get_filename(filepath_to_feature_file));
    //remove the ".bin"
    filepath_to_ws_image.resize(filepath_to_ws_image.size()-4);
    if (filepath_image!="") filepath_image.resize(filepath_image.size()-4);

    seg segment(true);
    
    //IMPORT RESULTS FROM HDF5 file
    std::string filepath_to_ws_region_image=filepath_to_ws_image;
    filepath_to_ws_region_image.append(".h5");

    segment.load_cgp_data_structure(filepath_to_ws_region_image);
    
    vigra::BasicImage<unsigned int> & ws_region_image = segment.ws_region_image;
    marray::Marray<unsigned int> & one_boundings =      segment.one_boundings;
    marray::Marray<unsigned int> & two_boundings =      segment.two_boundings;
    std::vector< std::vector<point> > & arcs =          segment.arcs;
    std::vector<point> & junctions =                    segment.junctions;
    int & dim_x =                                       segment.dim_x;
    int & dim_y =                                       segment.dim_y;
    size_t & numberOfBinningGroups =                    segment.numberOfBinningGroups;

    bool thumb=false;
    vigra::BasicImage<bool> selection_image(dim_x,dim_y);
    for (int x=0; x<dim_x; x++)
    {
        for (int y=0; y<dim_y; y++)
        {
            selection_image(x,y)=true;
        }
    }

    if (filepath_thumbs!="no")
    {
        filepath_thumbs.append(get_filename(filepath_to_feature_file));
        //remove the ".bin"
        filepath_thumbs.resize(filepath_thumbs.size()-4);
        filepath_thumbs.append(".reduced.bmp");

        std::ifstream thumb_file(filepath_thumbs.c_str());

        if(thumb_file)
        {
            thumb_file.close();
            std::cout<<"Reduced thumb found!"<<std::endl;

            //open thumb
            vigra::ImageImportInfo info(filepath_thumbs.c_str());
            int thumb_dim_x=0.2*dim_x;
            int thumb_dim_y=0.2*dim_y;

            if (info.width()!=thumb_dim_x || info.height()!=thumb_dim_y)
                std::cout << "Error: thumb has not expected size of "<< thumb_dim_x <<" x "<<thumb_dim_y<<"!"<<std::endl;
            else
            {
                vigra::IImage thumb_grayvalues(thumb_dim_x,thumb_dim_y);
                vigra::IImage label_image(thumb_dim_x,thumb_dim_y);
                vigra::BasicImage<bool> thumb_bool(thumb_dim_x,thumb_dim_y);

                importImage(info, destImage(thumb_grayvalues));

                for (int x=0;x<thumb_dim_x;x++)
                    for (int y=0;y<thumb_dim_y;y++)
                    {
                        if (thumb_grayvalues(x,y)<255) thumb_grayvalues(x,y)=0;
                        else thumb_grayvalues(x,y)=1;
                    }

                int nr_regions = vigra::labelImageWithBackground(vigra::srcImageRange(thumb_grayvalues), vigra::destImage(label_image), 0, 0);
                std::vector<bool> regions(nr_regions,false);

                //check for regions connected to border
                for (int y=0; y<thumb_dim_y; y++)
                {
                    if (label_image(0,y)>0) regions[label_image(0,y)-1]=true;
                    if (label_image(thumb_dim_x-1,y)>0) regions[label_image(thumb_dim_x-1,y)-1]=true;
                }

                for (int x=0; x<thumb_dim_x; x++)
                {
                    if (label_image(x,0)>0) regions[label_image(x,0)-1]=true;
                    if (label_image(x,thumb_dim_y-1)>0) regions[label_image(x,thumb_dim_y-1)-1]=true;
                }

                for (int x=0;x<thumb_dim_x;x++)
                    for (int y=0;y<thumb_dim_y;y++)
                        if (label_image(x,y)>0)
                        {
                            if (regions[label_image(x,y)-1]) thumb_bool(x,y)=false;
                            else thumb_bool(x,y)=true;
                        }
                        else thumb_bool(x,y)=true;

                for (int region=0; region<nr_regions; region++)
                    if (!regions[region])
                    {
                        std::cout<<"Inside white region ignored"<<std::endl;
                    }

                resizeImageNoInterpolation(srcImageRange(thumb_bool),destImageRange(selection_image));

                for (int x=0; x<dim_x; x++)
                {
                    selection_image(x,0)=false;
                    selection_image(x,dim_y-1)=false;
                }

                for (int y=0; y<dim_y; y++)
                {
                    selection_image(0,y)=false;
                    selection_image(dim_x-1,y)=false;
                }

                thumb=true;
            }
        }
        else std::cout<<"Thumb "<<filepath_thumbs<<" not found"<<std::endl;
    }

    if (!thumb)
    {
        std::string filepath_to_image_selection=filepath_to_feature_file;
        filepath_to_image_selection.resize(filepath_to_image_selection.size()-3);
        filepath_to_image_selection.append("selection.dat");

        std::ifstream selection_file(filepath_to_image_selection.c_str());
        std::cout<<"Importing image selection from file:"<<std::endl;
        std::cout<<filepath_to_image_selection<<std::endl;

        if (selection_file.is_open())
        {
            int low_x, high_x, low_y, high_y;

            selection_file>>low_x;
            selection_file>>high_x;
            selection_file>>low_y;
            selection_file>>high_y;

            selection_file.close();

            for (int x=0; x<=low_x; x++)
            {
                for (int y=0; y<dim_y; y++)
                {
                    selection_image(x,y)=false;
                }
            }

            for (int x=low_x+1; x<high_x; x++)
            {
                for (int y=0; y<=low_y; y++)
                {
                    selection_image(x,y)=false;
                }

                for (int y=high_y; y<dim_y; y++)
                {
                    selection_image(x,y)=false;
                }
            }

            for (int x=high_x; x<dim_x; x++)
            {
                for (int y=0; y<dim_y; y++)
                {
                    selection_image(x,y)=false;
                }
            }

            std::cout<<"...done"<<std::endl;
        }
        else std::cout << "Selection file not found, default values used"<<std::endl;
    }

    int nr_of_features=0;
    int nr_of_arcs=0;

    //First we open the random_forest file
    filepath_to_random_forest_file.append("_rf.hdf5");
    std::ifstream check_file(filepath_to_random_forest_file.c_str());
    if(!check_file)
    {
        std::cout<<"extract_boundary_probabilities(..)  Error: Random Forest File "<<filepath_to_random_forest_file<<" is NOT existend"<<std::endl;
        exit(-1);
    }
    check_file.close();

    vigra::RandomForest<> rf;
    rf_import_HDF5(rf,filepath_to_random_forest_file.c_str());

    //NOW FOR UNKNOWN IMAGES/FEATURES
    FILE *fp_unknown_features;
    fp_unknown_features =fopen(filepath_to_feature_file.c_str(),"rb");

    if(fp_unknown_features==NULL)
    {
        std::cout<<"extract_boundary_probabilities(..)  Error: Feature File "<<filepath_to_feature_file<<" is NOT existend"<<std::endl;
        exit(-1);
    }

    float float_nr_of_arcs=0;
    float float_nr_of_features=0;

    int read_in_arc_size=fread(&float_nr_of_arcs,sizeof(float),1,fp_unknown_features);
    if(read_in_arc_size!=1)
    {
        std::cout<<" extract_boundary_probabilities(..) error, could not read in arc size"<<std::endl;
        exit(-1);
    }

    int read_in_nr_of_features=fread(&float_nr_of_features,sizeof(float),1,fp_unknown_features);
    if(read_in_nr_of_features!=1)
    {
        std::cout<<"extract_boundary_probabilities(..) error, could not read in boundary feature size"<<std::endl;
        exit(-1);
    }

    nr_of_arcs=(int)float_nr_of_arcs;
    if(nr_of_arcs!=arcs.size())
    {
        std::cout<<"extract_boundary_probabilities(..) error, number of arcs in segmentation does not fit to boundary features"<<std::endl;
        exit(-1);
    }

    //NOW WE KNOW HOW BIG THE unknow_feature_array has to be
    nr_of_features=(int)float_nr_of_features;
    nr_of_arcs=(int)float_nr_of_arcs;
    std::cout<<"nr of boundary features: "<<nr_of_features<<std::endl;
    float * unknown_features_array = new float[nr_of_features*nr_of_arcs];

    int read_in=fread(unknown_features_array,sizeof(float),nr_of_features*nr_of_arcs,fp_unknown_features);
    if(read_in!=nr_of_features*nr_of_arcs)
    {
        std::cout<<"extract_boundary_probabilities() error, in unknown features array, could not load all the floats"<<std::endl;
        exit(-1);
    }

    for(size_t jj=0;jj<nr_of_arcs*nr_of_features;jj++)
    {
        if (std::isnan(unknown_features_array[jj]))
        {
            unknown_features_array[jj] = 0;
        }
    }

    fclose(fp_unknown_features);

    vigra::MultiArray<2, float> unknown_features(vigra::MultiArrayShape<2>::type(nr_of_arcs,nr_of_features),unknown_features_array);
    for(int x=0;x<nr_of_arcs;x++)
    {
        for(int y=0;y<nr_of_features;y++)
        {
            unknown_features(x,y)=unknown_features_array[nr_of_features*x+y];
        }
    }

    delete unknown_features_array;

    vigra::MultiArray<2, double> unknown_probability(vigra::MultiArrayShape<2>::type(nr_of_arcs,nr_of_classes));  //nr_of_classes is a global variable

    std::cout<<"Probability prediction...."<<std::endl;
    rf.predictProbabilities(unknown_features,unknown_probability);
    std::cout<<"....done"<<std::endl;
/*
    save_boundary_probabilities(unknown_probability,
                                filepath_to_ws_region_image,
                                path_to_output_folder,
                                segment,
                                param_file_name);
*/
    //Create a new instance of class gbn
    gbn GrainBoundNet; 
    GrainBoundNet.found_bubble_arcs.resize((int)arcs.size(), false);
    std::vector<bool> & found_bubble_arcs = GrainBoundNet.found_bubble_arcs;
    std::vector<int> & found_bubble_areas = GrainBoundNet.found_bubble_areas;
    std::vector< std::vector<int> > & bubble_arc_index = GrainBoundNet.bubble_arc_index;
    std::vector<int> remove_bubble_arcs((int)arcs.size(), false);
    GrainBoundNet.paramFileGBN = &paramFile;

    //Do not perform bubble recognition if the original image does not exist
    if(originalImageExists == true)
    {
        find_bubble_arcs(unknown_probability,
                           filepath_to_ws_image,
                           path_to_output_folder,
                           segment,
                           GrainBoundNet,
                           param_file,
                           paramFile,
                           selection_image,
                           filepath_image,
                           remove_bubble_arcs);
   }

    std::list<int> found_border_areas;

    //combined grain/bubble boundary probability
    std::vector<double> boundary_probability(nr_of_arcs);
    for (int arcindex=0; arcindex<nr_of_arcs; arcindex++)
    {
        if (found_bubble_arcs[arcindex]==true)
        {
            boundary_probability[arcindex]=1.0f;
        }
        else if (remove_bubble_arcs[arcindex]==true)
        {
            boundary_probability[arcindex]=0.0f;
        }
        else
        {   //both arc endpoint are within selected region
            if (selection_image(arcs[arcindex][0].x,arcs[arcindex][0].y) &&
                selection_image(arcs[arcindex][arcs[arcindex].size()-1].x,arcs[arcindex][arcs[arcindex].size()-1].y))
            {
                //probability 0.00 means forced no-boundary, 1.00 means bubble arc
                boundary_probability[arcindex]=std::max(0.001,std::min(0.999,unknown_probability(arcindex,1)));
            }
            else
            {
                boundary_probability[arcindex]=0.0f;
                //set as border area that have to merge with outside label
                found_border_areas.push_back(two_boundings(arcindex,0));
                found_border_areas.push_back(two_boundings(arcindex,1));
            }
        }
    }

    found_border_areas.sort();
    found_border_areas.unique();

    //REDUCE BOUNDARY PROBABILITY FOR VERTICAL HOUGH LABELED ARCS
    std::vector<int> vertical_arc_index;
    find_vertical_arcs(arcs,vertical_arc_index,found_bubble_arcs,ws_region_image,dim_x,dim_y);

    //REDUCE BOUNDARY PROBABILITY FOR MOSAIC BORDERS
    find_mosaic_borders(arcs,vertical_arc_index,found_bubble_arcs,two_boundings,dim_x,dim_y,get_filename(filepath_to_feature_file));

    for(int found=0; found<vertical_arc_index.size(); found++)
    {
         boundary_probability[vertical_arc_index[found]]=0.0f;
    }

    std::vector<size_t> & region_labels = GrainBoundNet.region_labels;
    std::vector< std::vector<int> > & grain_arc_index = GrainBoundNet.grain_arc_index;

    //Set no_boundary_threshold to 1.0 and minimal_grain_size to 1 if the original image does not exist
    if(originalImageExists == false)
    {
        Parameter<float> no_boundary_threshold;
        no_boundary_threshold.assign("", "no_boundary_threshold", 1.0);
        no_boundary_threshold.save(paramFile, "config");

        Parameter<int> minimal_grain_size;
        minimal_grain_size.assign("", "minimal_grain_size", 1);
        minimal_grain_size.save(paramFile, "config");
    }
    
    find_grain_arcs(segment,
                    GrainBoundNet,
                    found_border_areas,
                    unknown_probability,
                    boundary_probability,
                    paramFile,
                    selection_image,
                    filepath_image,
                    unknown_features,
                    rf);

    combine_bubbles_grains(filepath_to_ws_image,
                           path_to_output_folder,
                           segment,
                           GrainBoundNet,
                           param_file,
                           paramFile,
                           selection_image);
}
