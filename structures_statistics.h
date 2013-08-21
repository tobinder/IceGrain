/*! \file structures_statistics.h
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
#include "param.h"

void junctions_and_outer_circle(size_t nr_areas,
                                long * grain_area_size,
                                std::vector< std::vector< std::vector<int> > > & grain_arc_index,
                                std::vector<int> grain_junctions,
                                std::vector<int> grain_bubble_junctions,
                                std::vector<point> grain_area_center_mass,
                                std::vector<point> bubble_area_center_mass,
                                std::vector<bool> grain_arc,
                                int minimal_bubble_distance,
                                int close_bubble_grain_size,
                                marray::Marray<unsigned int> one_boundings,
                                marray::Marray<unsigned int> two_boundings,
                                std::vector< std::vector<point> > arcs,
                                std::vector<point> junctions,
                                int pixels_average,
                                std::vector< std::vector<point> > areas,
                                std::vector<area_range> & area_ranges,
                                std::vector<int> found_bubble_areas,
                                std::list<int> found_border_areas,
                                std::vector<int> & close_bubble_areas,
                                std::vector< std::vector<int> > & grain_area_junctions,
                                std::vector<bool> & grain_junction,
                                std::vector<std::vector<int> > arc_junctions,
                                std::vector<int> & grain_perimeter,
                                std::vector<int> & min_bubble_distance,
                                std::vector<int> & grain_longest_arc_length,
                                std::vector< std::vector<float> > & grain_junction_angles,
                                std::vector< std::vector<float> > & grain_junction_angles2,
                                std::vector< std::vector<point> > & grain_boundary_pixels,
                                std::vector< std::vector<int> > & grain_boundary_index,
                                std::vector< std::vector<float> > & grain_boundary_phis,
                                std::vector< std::vector<float> > & grain_boundary_curvs,
                                std::vector< std::vector<unsigned int> > & grain_area_boundaries,
                                std::vector< std::vector<int> > & grain_junction_index,
                                std::vector<bool> & arc_to_segment,
                                vigra::BasicImage<unsigned int> ws_region_image,
                                std::vector<size_t> region_labels)
{
    //fill/update grain_area_junctions and grain_junctions
    for (int j=0; j<grain_junctions.size(); j++)
    {
        int y=grain_junctions[j];//junction index

        for (int area=0; area<nr_areas; area++)
        {
            if (grain_arc_index[area].size()>0)//if area is grain, loop over grain arcs
                for (int a=0; a<grain_arc_index[area][0].size(); a++)
                {
                    int arc_index=grain_arc_index[area][0][a]+1;

                    if (one_boundings(y,0)==arc_index || one_boundings(y,1)==arc_index || one_boundings(y,2)==arc_index || one_boundings(y,3)==arc_index)
                    {
                        grain_area_junctions[area].push_back(y);
                        grain_junction[y]=true;
                        break;
                    }
                }
        }
    }

    for (int j=0; j<grain_bubble_junctions.size(); j++)
    {
        int y=grain_bubble_junctions[j];//junction index

        for (int area=0; area<nr_areas; area++)
        {
            if (grain_arc_index[area].size()>0)//if area is grain, loop over grain arcs
                for (int a=0; a<grain_arc_index[area][0].size(); a++)
                {
                    int arc_index=grain_arc_index[area][0][a]+1;

                    if (one_boundings(y,0)==arc_index || one_boundings(y,1)==arc_index || one_boundings(y,2)==arc_index || one_boundings(y,3)==arc_index)
                    {
                        grain_area_junctions[area].push_back(y);
                        break;
                    }
                }
        }
    }

    //find/update an outer circle of grain arcs, calculate min bubble distance, grain segment length, grain perimeter and grain junction angles
    //fill close bubble areas vector if selected
    for (int area=0; area<nr_areas; area++)
    {
        grain_longest_arc_length[area]=0;
        grain_perimeter[area]=0;
        grain_junction_angles[area].clear();
        grain_junction_angles2[area].clear();

        //SPECIAL CASE: grain with no junction -> only one segment / no bubble arcs in arc list
        if (grain_arc_index[area].size()>0 && grain_area_junctions[area].size()==0)
        {
            //store pixels/indeces along all grain arcs
            std::vector<point> this_grain;
            grain_boundary_pixels.push_back(this_grain);
            std::vector<int> segment_index;
            grain_boundary_index.push_back(segment_index);
            grain_area_boundaries[area].push_back(grain_boundary_index.size()-1);

            //calculate min bubble distance
            min_bubble_distance[area]=-1;
            for (int bubble=0; bubble<bubble_area_center_mass.size(); bubble++)
            {
                int distance=sqrt(square(grain_area_center_mass[area].x-bubble_area_center_mass[bubble].x)
                                  + square(grain_area_center_mass[area].y-bubble_area_center_mass[bubble].y));
                if(distance<min_bubble_distance[area] || min_bubble_distance[area]==-1) min_bubble_distance[area]=distance;
            }

            //if there are no bubbles
            if (bubble_area_center_mass.size()==0) min_bubble_distance[area]=minimal_bubble_distance;

            //if grain is closer to bubble than threshold and smaller than close_bubble_grain_size
            if(min_bubble_distance[area]<minimal_bubble_distance && areas[area].size()<close_bubble_grain_size)
            {
                close_bubble_areas.push_back(area);
            }

            //calculate grain perimeter
            for(int a=0; a<grain_arc_index[area][0].size(); a++)
            {
                grain_perimeter[area]+=get_length(arcs[grain_arc_index[area][0][a]]);
            }

            grain_longest_arc_length[area]=grain_perimeter[area];//length of single segment is perimeter

            //add first arc to grain boundary pixels
            for(int p=0; p<arcs[grain_arc_index[area][0][0]].size(); p++)
            {
                grain_boundary_pixels.back().push_back(arcs[grain_arc_index[area][0][0]][p]);
            }

            grain_boundary_index.back().push_back(grain_arc_index[area][0][0]+1);

            //if first arc is not isolated find a closed circle
            if(arc_junctions[grain_arc_index[area][0][0]].size()==2)
            {
                int start_junction, junction_back;

                //check which junction fits to end of start arc
                point back=grain_boundary_pixels.back().back();

                if((fabs(back.x-junctions[arc_junctions[grain_arc_index[area][0][0]][0]].x)+
                    fabs(back.y-junctions[arc_junctions[grain_arc_index[area][0][0]][0]].y))<2)
                {
                    start_junction=arc_junctions[grain_arc_index[area][0][0]][1];
                    junction_back=arc_junctions[grain_arc_index[area][0][0]][0];
                }
                else if((fabs(back.x-junctions[arc_junctions[grain_arc_index[area][0][0]][1]].x)+
                         fabs(back.y-junctions[arc_junctions[grain_arc_index[area][0][0]][1]].y))<2)
                {
                    start_junction=arc_junctions[grain_arc_index[area][0][0]][0];
                    junction_back=arc_junctions[grain_arc_index[area][0][0]][1];
                }
                else
                {
                    std::cout<<"Error: No junction for start arc found!"<<std::endl;
                    exit(-1);
                }

                std::vector<bool> arc_used;
                arc_used.resize(grain_arc_index[area][0].size(),false);
                arc_used[0]=true;

                while(junction_back!=start_junction)//runs until start junction is reached again
                {
                    bool next_arc_found=false;

                    for (int arc=0; arc<grain_arc_index[area][0].size() && next_arc_found==false; arc++)//loop over grain arcs of this area
                    {
                        int arc_index=grain_arc_index[area][0][arc];

                        for (int j=0; j<arc_junctions[arc_index].size() && next_arc_found==false; j++)//loop over junctions of this arc
                        {
                            if (arc_junctions[arc_index][j]==junction_back && arc_used[arc]==false)
                            {
                                int junction_front=arc_junctions[arc_index][j];//junction we looked for
                                junction_back=arc_junctions[arc_index][(j+1)%2];//junction on the other side of the arc
                                next_arc_found=true;
                                arc_used[arc]=true;

                                //check which side of found arc fits to junction
                                point front=arcs[arc_index][0];
                                point back=arcs[arc_index].back();

                                if((fabs(front.x-junctions[junction_front].x)+fabs(front.y-junctions[junction_front].y))<2)
                                {
                                    for(int p=0; p<arcs[arc_index].size(); p++)
                                    {
                                        grain_boundary_pixels.back().push_back(arcs[arc_index][p]);
                                    }
                                    grain_boundary_index.back().push_back(arc_index+1);
                                }
                                else if((fabs(back.x-junctions[junction_front].x)+fabs(back.y-junctions[junction_front].y))<2)
                                {
                                    for(int p=arcs[arc_index].size(); p>0; p--)
                                    {
                                        grain_boundary_pixels.back().push_back(arcs[arc_index][p-1]);
                                    }
                                    grain_boundary_index.back().push_back(-arc_index-1);
                                }
                                else
                                {
                                    std::cout<<"Error: Junction do not fit to arcs!"<<std::endl;
                                    exit(-1);
                                }
                            }
                        }
                    }
                }
            }

        }
        //NORMAL CASE: at least one junction, combine arcs to segments
        else if (grain_arc_index[area].size()>0 && grain_area_junctions[area].size()>0)//find start arc with junction
        {
            std::vector< std::vector< std::vector<point> > > junction_pixels;
            std::vector< std::vector< std::vector<point> > > junction_pixels2;
            junction_pixels.resize(grain_area_junctions[area].size());
            junction_pixels2.resize(grain_area_junctions[area].size());
            std::vector< std::vector<int> > junction_arcs;
            junction_arcs.resize(grain_area_junctions[area].size());
            std::vector<int> junction_index;
            junction_index.resize(grain_area_junctions[area].size());
            std::vector<bool> junction_used;
            junction_used.resize(junctions.size(),false);
            std::vector<int> junction_back;
            std::vector<point> last_segment;
            std::vector<bool> arc_used;
            arc_used.resize(two_boundings.shape(0),false);
            int last_arc;

            //store the pixels of all segments of this grain
            std::vector< std::vector<point> > this_grain;
            std::vector<point> one_segment;
            this_grain.push_back(one_segment);
            std::vector< std::vector<int> > this_grain_arcs;
            std::vector<int> segment_index;
            this_grain_arcs.push_back(segment_index);
            grain_longest_arc_length[area]=0;

            //to make the detection of an outer circle 100% safe junction_used is set true for all inner junctions
            vigra::BasicImage<bool> temp_img(area_ranges[area].x_high-area_ranges[area].x_low+3, area_ranges[area].y_high-area_ranges[area].y_low+3);

            for(int y=0;y<temp_img.height();y++)
            {
                for (int x=0;x<temp_img.width();x++)
                {
                    temp_img(x,y)=true;
                }
            }

            //loop over all area pixels
            for (int pixel=0;pixel<areas[area].size();pixel++)
            {
                point p=areas[area][pixel];
                temp_img(p.x-area_ranges[area].x_low+1,p.y-area_ranges[area].y_low+1)=false;
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
                            inside_labels.push_back(region_labels[ws_region_image(x+area_ranges[area].x_low-1,y+area_ranges[area].y_low-1)-1]-1);
                    }
                }

                inside_labels.sort();
                inside_labels.unique();

                for (std::list<int>::iterator it=inside_labels.begin(); it!=inside_labels.end(); ++it)
                {
                    //std::cout << "exclude inside junctions of label "<< *it+1 << std::endl;

                    for (int j=0;j<grain_area_junctions[*it].size();j++)
                    {
                        junction_used[grain_area_junctions[*it][j]]=true;
                    }
                }

                bool junction_left=false;

                for (int j=0;j<grain_area_junctions[area].size() && !junction_left; j++)
                {
                    if (!junction_used[grain_area_junctions[area][j]]) junction_left=true;
                }

                if(!junction_left)//this is special case, but it must be checked that only outer arcs are considered
                {
                    std::vector<int> reduced_grain_arc_index;

                    for(int a=0; a<grain_arc_index[area][0].size(); a++)
                    {
                        if (labels_check(arcs[grain_arc_index[area][0][a]][0].x-area_ranges[area].x_low+1,
                            arcs[grain_arc_index[area][0][a]][0].y-area_ranges[area].y_low+1)==0 ||
                            labels_check(arcs[grain_arc_index[area][0][a]][0].x-area_ranges[area].x_low+1,
                            arcs[grain_arc_index[area][0][a]][0].y-area_ranges[area].y_low+1)==labels_check(0,0))
                        reduced_grain_arc_index.push_back(grain_arc_index[area][0][a]);
                    }

                    //store pixels/indeces along all grain arcs
                    std::vector<point> this_grain;
                    grain_boundary_pixels.push_back(this_grain);
                    std::vector<int> segment_index;
                    grain_boundary_index.push_back(segment_index);
                    grain_area_boundaries[area].push_back(grain_boundary_index.size()-1);

                    //calculate min bubble distance
                    min_bubble_distance[area]=-1;
                    for (int bubble=0; bubble<bubble_area_center_mass.size(); bubble++)
                    {
                        int distance=sqrt(square(grain_area_center_mass[area].x-bubble_area_center_mass[bubble].x)
                                          + square(grain_area_center_mass[area].y-bubble_area_center_mass[bubble].y));
                        if(distance<min_bubble_distance[area] || min_bubble_distance[area]==-1) min_bubble_distance[area]=distance;
                    }

                    //if there are no bubbles
                    if (bubble_area_center_mass.size()==0) min_bubble_distance[area]=minimal_bubble_distance;

                    //if grain is closer to bubble than threshold and smaller than close_bubble_grain_size
                    if(min_bubble_distance[area]<minimal_bubble_distance && areas[area].size()<close_bubble_grain_size)
                    {
                        close_bubble_areas.push_back(area);
                    }

                    //calculate grain perimeter
                    for(int a=0; a<reduced_grain_arc_index.size(); a++)
                    {
                        grain_perimeter[area]+=get_length(arcs[reduced_grain_arc_index[a]]);
                    }

                    grain_longest_arc_length[area]=grain_perimeter[area];//length of single segment is perimeter

                    //add first arc to grain boundary pixels
                    for(int p=0; p<arcs[reduced_grain_arc_index[0]].size(); p++)
                    {
                        grain_boundary_pixels.back().push_back(arcs[reduced_grain_arc_index[0]][p]);
                    }

                    grain_boundary_index.back().push_back(reduced_grain_arc_index[0]+1);

                    //if first arc is not isolated find a closed circle
                    if(arc_junctions[reduced_grain_arc_index[0]].size()==2)
                    {
                        int start_junction, junction_back;

                        //check which junction fits to end of start arc
                        point back=grain_boundary_pixels.back().back();

                        if((fabs(back.x-junctions[arc_junctions[reduced_grain_arc_index[0]][0]].x)+
                            fabs(back.y-junctions[arc_junctions[reduced_grain_arc_index[0]][0]].y))<2)
                        {
                            start_junction=arc_junctions[reduced_grain_arc_index[0]][1];
                            junction_back=arc_junctions[reduced_grain_arc_index[0]][0];
                        }
                        else if((fabs(back.x-junctions[arc_junctions[reduced_grain_arc_index[0]][1]].x)+
                                 fabs(back.y-junctions[arc_junctions[reduced_grain_arc_index[0]][1]].y))<2)
                        {
                            start_junction=arc_junctions[reduced_grain_arc_index[0]][0];
                            junction_back=arc_junctions[reduced_grain_arc_index[0]][1];
                        }
                        else
                        {
                            std::cout<<"Error: No junction for start arc found!"<<std::endl;
                            exit(-1);
                        }

                        std::vector<bool> arc_used;
                        arc_used.resize(reduced_grain_arc_index.size(),false);
                        arc_used[0]=true;

                        while(junction_back!=start_junction)//runs until start junction is reached again
                        {
                            bool next_arc_found=false;

                            for (int arc=0; arc<reduced_grain_arc_index.size() && next_arc_found==false; arc++)//loop over grain arcs of this area
                            {
                                int arc_index=reduced_grain_arc_index[arc];

                                for (int j=0; j<arc_junctions[arc_index].size() && next_arc_found==false; j++)//loop over junctions of this arc
                                {
                                    if (arc_junctions[arc_index][j]==junction_back && arc_used[arc]==false)
                                    {
                                        int junction_front=arc_junctions[arc_index][j];//junction we looked for
                                        junction_back=arc_junctions[arc_index][(j+1)%2];//junction on the other side of the arc
                                        next_arc_found=true;
                                        arc_used[arc]=true;

                                        //check which side of found arc fits to junction
                                        point front=arcs[arc_index][0];
                                        point back=arcs[arc_index].back();

                                        if((fabs(front.x-junctions[junction_front].x)+fabs(front.y-junctions[junction_front].y))<2)
                                        {
                                            for(int p=0; p<arcs[arc_index].size(); p++)
                                            {
                                                grain_boundary_pixels.back().push_back(arcs[arc_index][p]);
                                            }
                                            grain_boundary_index.back().push_back(arc_index+1);
                                        }
                                        else if((fabs(back.x-junctions[junction_front].x)+fabs(back.y-junctions[junction_front].y))<2)
                                        {
                                            for(int p=arcs[arc_index].size(); p>0; p--)
                                            {
                                                grain_boundary_pixels.back().push_back(arcs[arc_index][p-1]);
                                            }
                                            grain_boundary_index.back().push_back(-arc_index-1);
                                        }
                                        else
                                        {
                                            std::cout<<"Error: Junction do not fit to arcs!"<<std::endl;
                                            exit(-1);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    continue;
                }
            }

            int start_arc=-1;
            int start_junction;
            int junction_front;
            int last_junction;
            int junction_pos=grain_area_junctions[area].size()-1;

            for (int arc=0; arc<grain_arc_index[area][0].size() && start_arc==-1; arc++)//loop over grain arcs of this area
            {
                int arc_index=grain_arc_index[area][0][arc];

                for (int j=0; j<arc_junctions[arc_index].size() && start_arc==-1; j++)//loop over junctions of this arc
                {
                    for (int k=0; k<grain_area_junctions[area].size() && start_arc==-1; k++)//loop over junctions of this area
                    {
                        if (arc_junctions[arc_index][j]==grain_area_junctions[area][k] && !junction_used[grain_area_junctions[area][k]])
                        {
                            start_arc=arc_index;//start arc found
                            junction_front=arc_junctions[arc_index][j];//junction we looked for
                            start_junction=junction_front;
                            junction_back.push_back(arc_junctions[arc_index][(j+1)%2]);//junction on the other side of the arc
                            arc_used[arc_index]=true;
                            junction_used[grain_area_junctions[area][k]]=true;
                            last_arc=start_arc;
                        }
                    }
                }
            }

            if(start_arc==-1)
            {
                std::cout<<"Error: No start junction found!"<<std::endl;
                exit(-1);
            }

            //fill junction_pixels for second side of end junction
            std::vector<point> second_side;
            std::vector<point> second_side2;

            //check which side of found arc fits to junction
            point front=arcs[start_arc][0];
            point back=arcs[start_arc].back();

            if((fabs(front.x-junctions[start_junction].x)+fabs(front.y-junctions[start_junction].y))<2)
            {
                for(int p=0; p<arcs[start_arc].size(); p++)
                {
                    last_segment.push_back(arcs[start_arc][p]);
                    if(grain_arc[start_arc])
                    {
                        this_grain.back().push_back(arcs[start_arc][p]);
                    }
                    if (second_side.size()<pixels_average) second_side.push_back(arcs[start_arc][p]);
                    second_side2.push_back(arcs[start_arc][p]);
                }
                //add start arc to first segment
                this_grain_arcs.back().push_back(start_arc+1);
            }
            else if((fabs(back.x-junctions[start_junction].x)+fabs(back.y-junctions[start_junction].y))<2)
            {
                for(int p=arcs[start_arc].size(); p>0; p--)
                {
                    last_segment.push_back(arcs[start_arc][p-1]);
                    if(grain_arc[start_arc])
                    {
                        this_grain.back().push_back(arcs[start_arc][p-1]);
                    }
                    if (second_side.size()<pixels_average) second_side.push_back(arcs[start_arc][p-1]);
                    second_side2.push_back(arcs[start_arc][p-1]);
                }
                //add start arc to first segment
                this_grain_arcs.back().push_back(-start_arc-1);
            }
            else
            {
                std::cout<<"Error: Start junction do not fit to arcs!"<<std::endl;
                exit(-1);
            }

            junction_pixels[junction_pos].resize(2);
            junction_pixels2[junction_pos].resize(2);
            junction_pixels[junction_pos][1]=second_side;
            junction_pixels2[junction_pos][1]=second_side2;
            junction_arcs[junction_pos].resize(2);
            junction_arcs[junction_pos][1]=start_arc;
            junction_index[junction_pos]=start_junction;

            while(junction_back.back()!=start_junction)//runs until start junction is reached again
            {
                bool next_arc_found=false;

                for (int k=0; k<grain_area_junctions[area].size() && next_arc_found==false; k++)//loop over junctions of this area
                {
                    if (junction_back.back()==grain_area_junctions[area][k])//real junction reached
                    {
                        if(junction_pos<grain_area_junctions[area].size()-1)
                        {
                            //find other arc in this junction
                            point other_arc_midpoint;
                            other_arc_midpoint.x=-1;

                            for(int x=0;x<(int)one_boundings.shape(1);x++)
                            {
                                if (one_boundings(last_junction,x)-1!=junction_arcs[junction_pos][0] &&
                                    one_boundings(last_junction,x)-1!=junction_arcs[junction_pos][1])
                                {
                                    other_arc_midpoint=arcs[one_boundings(last_junction,x)-1][arcs[one_boundings(last_junction,x)-1].size()/2];
                                    break;
                                }
                            }

                            //calculate angle in this junction
                            float angle=calculate_angle(junction_pixels[junction_pos],junctions[last_junction],other_arc_midpoint);
                            float angle2=calculate_angle2(junction_pixels2[junction_pos],junctions[last_junction],other_arc_midpoint);
                            if(grain_junction[last_junction]) grain_junction_angles[area].push_back(angle);
                            if(grain_junction[last_junction]) grain_junction_angles2[area].push_back(angle2);
                        }

                        junction_pos=(junction_pos+1)%grain_area_junctions[area].size();
                        last_junction=junction_back.back();

                        //fill first side in backwards direction
                        std::vector<point> first_side;
                        std::vector<point> first_side2;

                        for(int p=last_segment.size(); p>0 && first_side.size()<pixels_average; p--)
                        {
                            first_side.push_back(last_segment[p-1]);
                        }

                        for(int p=last_segment.size(); p>0; p--)
                        {
                            first_side2.push_back(last_segment[p-1]);
                        }

                        junction_pixels[junction_pos].resize(2);
                        junction_pixels2[junction_pos].resize(2);
                        junction_pixels[junction_pos][0]=first_side;
                        junction_pixels2[junction_pos][0]=first_side2;
                        junction_arcs[junction_pos].resize(2);
                        junction_arcs[junction_pos][0]=last_arc;
                        junction_arcs[junction_pos][1]=-1;
                        junction_index[junction_pos]=last_junction;

                        //save last segment length and clear last segment pixels
                        grain_perimeter[area]+=get_length(last_segment);
                        if(get_length(last_segment)>grain_longest_arc_length[area]) grain_longest_arc_length[area]=get_length(last_segment);
                        last_segment.clear();

                        std::vector<point> one_segment;
                        this_grain.push_back(one_segment);
                        std::vector<int> segment_index;
                        this_grain_arcs.push_back(segment_index);

                        junction_used[grain_area_junctions[area][k]]=true;
                    }
                }

                for (int arc=0; arc<grain_arc_index[area][0].size() && next_arc_found==false; arc++)//loop over grain arcs of this area
                {
                    int arc_index=grain_arc_index[area][0][arc];

                    for (int j=0; j<arc_junctions[arc_index].size() && next_arc_found==false; j++)//loop over junctions of this arc
                    {
                        if (arc_junctions[arc_index][j]==junction_back.back() && arc_used[arc_index]==false)
                        {
                            junction_front=arc_junctions[arc_index][j];//junction we looked for
                            junction_back.push_back(arc_junctions[arc_index][(j+1)%2]);//junction on the other side of the arc
                            next_arc_found=true;
                            arc_used[arc_index]=true;
                            last_arc=arc_index;
                            if(junction_arcs[junction_pos][1]==-1) junction_arcs[junction_pos][1]=arc_index;

                            //check which side of found arc fits to junction
                            point front=arcs[arc_index][0];
                            point back=arcs[arc_index].back();

                            if((fabs(front.x-junctions[junction_front].x)+fabs(front.y-junctions[junction_front].y))<2)
                            {
                                for(int p=0; p<arcs[arc_index].size(); p++)
                                {
                                    last_segment.push_back(arcs[arc_index][p]);
                                    if(grain_arc[arc_index])
                                    {
                                        this_grain.back().push_back(arcs[arc_index][p]);
                                    }
                                    if(junction_pixels[junction_pos][1].size()<pixels_average)
                                        junction_pixels[junction_pos][1].push_back(arcs[arc_index][p]);
                                    junction_pixels2[junction_pos][1].push_back(arcs[arc_index][p]);
                                }
                                this_grain_arcs.back().push_back(arc_index+1);
                            }
                            else if((fabs(back.x-junctions[junction_front].x)+fabs(back.y-junctions[junction_front].y))<2)
                            {
                                for(int p=arcs[arc_index].size(); p>0; p--)
                                {
                                    last_segment.push_back(arcs[arc_index][p-1]);
                                    if(grain_arc[arc_index])
                                    {
                                        this_grain.back().push_back(arcs[arc_index][p-1]);
                                    }
                                    if(junction_pixels[junction_pos][1].size()<pixels_average)
                                        junction_pixels[junction_pos][1].push_back(arcs[arc_index][p-1]);
                                    junction_pixels2[junction_pos][1].push_back(arcs[arc_index][p-1]);
                                }
                                this_grain_arcs.back().push_back(-arc_index-1);
                            }
                            else
                            {
                                std::cout<<"Error: Junction do not fit to arcs!"<<std::endl;
                                exit(-1);
                            }
                        }
                    }
                }

                if (next_arc_found==false)
                {
                    if(junction_back.size()==1)
                    {
                        std::cout<<"Error: No connection in grain arcs found!"<<std::endl;
                        exit(-1);
                    }
                    junction_back.pop_back();
                }
            }

            if(grain_area_junctions[area].size()>1)
            {
                //find other arc in last but one junction
                point other_arc_midpoint;
                other_arc_midpoint.x=-1;

                for(int x=0;x<(int)one_boundings.shape(1);x++)
                {
                    if (one_boundings(last_junction,x)-1!=junction_arcs[junction_pos][0] &&
                        one_boundings(last_junction,x)-1!=junction_arcs[junction_pos][1])
                    {
                        other_arc_midpoint=arcs[one_boundings(last_junction,x)-1][arcs[one_boundings(last_junction,x)-1].size()/2];
                        break;
                    }
                }

                //calculate angle in last but one junction
                float angle=calculate_angle(junction_pixels[junction_pos],junctions[last_junction],other_arc_midpoint);
                float angle2=calculate_angle2(junction_pixels2[junction_pos],junctions[last_junction],other_arc_midpoint);
                if(grain_junction[last_junction]) grain_junction_angles[area].push_back(angle);
                if(grain_junction[last_junction]) grain_junction_angles2[area].push_back(angle2);
            }

            //fill first side of last junction in backwards direction
            std::vector<point> first_side;
            std::vector<point> first_side2;

            for(int p=last_segment.size(); p>0 && first_side.size()<pixels_average; p--)
            {
                first_side.push_back(last_segment[p-1]);
            }

            for(int p=last_segment.size(); p>0; p--)
            {
                first_side2.push_back(last_segment[p-1]);
            }

            junction_pixels[junction_pixels.size()-1][0]=first_side;
            junction_pixels2[junction_pixels.size()-1][0]=first_side2;
            junction_arcs[junction_pixels.size()-1][0]=last_arc;

            //find other arc in last junction
            point other_arc_midpoint;
            other_arc_midpoint.x=-1;

            for(int x=0;x<(int)one_boundings.shape(1);x++)
            {
                if (one_boundings(junction_back.back(),x)-1!=junction_arcs[junction_pixels.size()-1][0] &&
                    one_boundings(junction_back.back(),x)-1!=junction_arcs[junction_pixels.size()-1][1])
                {
                    other_arc_midpoint=arcs[one_boundings(junction_back.back(),x)-1][arcs[one_boundings(junction_back.back(),x)-1].size()/2];
                    break;
                }
            }

            //calculate angle in last junction
            float angle=calculate_angle(junction_pixels[junction_pixels.size()-1],junctions[junction_back.back()],other_arc_midpoint);
            float angle2=calculate_angle2(junction_pixels2[junction_pixels2.size()-1],junctions[junction_back.back()],other_arc_midpoint);
            if(grain_junction[junction_back.back()]) grain_junction_angles[area].push_back(angle);
            if(grain_junction[junction_back.back()]) grain_junction_angles2[area].push_back(angle2);

            grain_perimeter[area]+=get_length(last_segment);
            if(get_length(last_segment)>grain_longest_arc_length[area]) grain_longest_arc_length[area]=get_length(last_segment);

            //if there are unused arcs an inside structure exists
            bool inside_structure=false;

            for (int arc=0; arc<grain_arc_index[area][0].size() && !inside_structure; arc++)//loop over grain arcs of this area
            {
                int arc_index=grain_arc_index[area][0][arc];
                if (!arc_used[arc_index]) inside_structure=true;
            }

            if(inside_structure)
            {
                //std::cout<<"handle inside structure of grain area "<<area<<std::endl;

                //find all areas and arcs inside of this area by minima growing
                vigra::BasicImage<bool> circle_arcs(area_ranges[area].x_high+5-area_ranges[area].x_low, area_ranges[area].y_high+5-
                    area_ranges[area].y_low);
                circle_arcs = true;
                vigra::BasicImage<short int> circle_result(area_ranges[area].x_high+5-area_ranges[area].x_low,
                    area_ranges[area].y_high+5-area_ranges[area].y_low);

                for (int arc=0; arc<grain_arc_index[area][0].size(); arc++)//loop over grain arcs of this area
                {
                    int arc_index=grain_arc_index[area][0][arc];
                    if (arc_used[arc_index])
                    {
                        for(int p=0; p<arcs[arc_index].size(); p++)
                        {
                            circle_arcs(arcs[arc_index][p].x+2-area_ranges[area].x_low,
                                        arcs[arc_index][p].y+2-area_ranges[area].y_low)=false;
                        }
                    }
                }

                // label the areas inside and outside different, assume that inside get maximal label
                int max_region_label = vigra::labelImageWithBackground(vigra::srcImageRange(circle_arcs), vigra::destImage(circle_result), false, 0);

                //calculate min bubble distance
                min_bubble_distance[area]=-1;
                for (int bubble=0; bubble<bubble_area_center_mass.size(); bubble++)
                {
                    bool inside_bubble=false;

                    if(area_ranges[area].x_low<=bubble_area_center_mass[bubble].x+2 && bubble_area_center_mass[bubble].x+2<=area_ranges[area].x_high &&
                       area_ranges[area].y_low<=bubble_area_center_mass[bubble].y+2 && bubble_area_center_mass[bubble].y+2<=area_ranges[area].y_high)
                    {
                        //check that bubble is not inside this grain
                        if(circle_result(bubble_area_center_mass[bubble].x+2-area_ranges[area].x_low,
                                         bubble_area_center_mass[bubble].y+2-area_ranges[area].y_low)==max_region_label)
                            inside_bubble=true;
                    }

                    int distance=sqrt(square(grain_area_center_mass[area].x-bubble_area_center_mass[bubble].x)
                                      + square(grain_area_center_mass[area].y-bubble_area_center_mass[bubble].y));
                    if((distance<min_bubble_distance[area] || min_bubble_distance[area]==-1) && !inside_bubble) min_bubble_distance[area]=distance;
                }

            }
            else//no inside structure
            {
                //calculate min bubble distance
                min_bubble_distance[area]=-1;
                for (int bubble=0; bubble<bubble_area_center_mass.size(); bubble++)
                {
                    int distance=sqrt(square(grain_area_center_mass[area].x-bubble_area_center_mass[bubble].x)
                                      + square(grain_area_center_mass[area].y-bubble_area_center_mass[bubble].y));
                    if(distance<min_bubble_distance[area] || min_bubble_distance[area]==-1) min_bubble_distance[area]=distance;
                }
            }

            //if there are no bubbles
            if (bubble_area_center_mass.size()==0) min_bubble_distance[area]=minimal_bubble_distance;

            //if grain is closer to bubble than threshold and smaller than close_bubble_grain_size
            if(min_bubble_distance[area]<minimal_bubble_distance && areas[area].size()<close_bubble_grain_size)
            {
                close_bubble_areas.push_back(area);
            }

            //copy segments of this grain to global vector
            for(int segment=0; segment<this_grain_arcs.size(); segment++)
            {
                //add segments only once to global vector
                if(!arc_to_segment[fabs(this_grain_arcs[segment][0])-1])
                {
                    if(this_grain[segment].size()>0)//there can be empty segments = bubble arcs
                    {
                        grain_boundary_pixels.push_back(this_grain[segment]);
    
                        //save indeces to global vector
                        grain_boundary_index.push_back(this_grain_arcs[segment]);
                        grain_area_boundaries[area].push_back(grain_boundary_index.size()-1);

                        //correlate index of grain boundary to junctions
                        int arc_front = fabs(this_grain_arcs[segment][0])-1;
                        int arc_back = fabs(this_grain_arcs[segment].back())-1;
          
                        //loop over junctions of this grain
                        for (int j=0; j<junction_arcs.size(); j++)
                        {
                            if(junction_arcs[j].size()==2)//if junction is reached
                            {
                                if(junction_arcs[j][0]==arc_front) grain_junction_index[junction_index[j]].push_back(grain_boundary_index.size()-1);
                                else if(junction_arcs[j][1]==arc_front) grain_junction_index[junction_index[j]].push_back(grain_boundary_index.size()-1);

                                if (grain_junction_index[junction_index[j]].size()>1)
                                    for(int k=0; k<grain_junction_index[junction_index[j]].size()-1; k++)//check if index is already existing
                                    {
                                        if(grain_junction_index[junction_index[j]][k]==grain_junction_index[junction_index[j]].back())
                                        {
                                            grain_junction_index[junction_index[j]].pop_back();
                                            break;
                                        }
                                    }

                                if(junction_arcs[j][0]==arc_back) grain_junction_index[junction_index[j]].push_back(grain_boundary_index.size()-1);
                                else if(junction_arcs[j][1]==arc_back) grain_junction_index[junction_index[j]].push_back(grain_boundary_index.size()-1);

                                if (grain_junction_index[junction_index[j]].size()>1)
                                    for(int k=0; k<grain_junction_index[junction_index[j]].size()-1; k++)
                                    {
                                        if(grain_junction_index[junction_index[j]][k]==grain_junction_index[junction_index[j]].back())
                                        {
                                            grain_junction_index[junction_index[j]].pop_back();
                                            break;
                                        }
                                    }
                            }
                        }
                    }
                    for(int arc=0; arc<this_grain_arcs[segment].size(); arc++)
                    {
                        //save which arcs has been used
                        arc_to_segment[fabs(this_grain_arcs[segment][arc])-1]=true;
                    }
                }
                else for(int boundary=0; boundary<grain_boundary_index.size(); boundary++) //find existing boundary index
                {
                    if (fabs(grain_boundary_index[boundary][0]) == fabs(this_grain_arcs[segment][0]) ||
                        fabs(grain_boundary_index[boundary].back()) == fabs(this_grain_arcs[segment][0]))
                    {
                        grain_area_boundaries[area].push_back(boundary);
                        break;
                    }
                }
            }

            //up to here all arc indeces have been associated with the first segment, here we have more than one segment
            grain_arc_index[area].clear();

            //update grain arc index
            for(int segment=0; segment<this_grain_arcs.size(); segment++)
            {
                std::vector<int> segment_index;
                grain_arc_index[area].push_back(segment_index);

                for(int arc=0; arc<this_grain_arcs[segment].size(); arc++)
                {
                    grain_arc_index[area].back().push_back(fabs(this_grain_arcs[segment][arc])-1);                
                }
            }
        }
    }

    // calculate angle phi of local normal vectors
    for(int i=0; i<grain_boundary_pixels.size(); i++)
    {
        std::vector<float> phi;
        grain_boundary_phis.push_back(phi);

        if(grain_boundary_pixels[i].size() < 5 )
        {
            // if the size of an arc is smaller than 5, use a simpler approximation of the derivative

            //std::cout << "Short arc... " << grain_boundary_pixels[i].size() << " ";
            double phi0 = atan2(grain_boundary_pixels[i][1].y - grain_boundary_pixels[i][0].y,
                                grain_boundary_pixels[i][1].x - grain_boundary_pixels[i][0].x );
            grain_boundary_phis.back().push_back( phi0 + PI/2.0f );

            for(int j=1; j<grain_boundary_pixels[i].size()-1; j++)
            {
                double phi_temp = atan2(grain_boundary_pixels[i][j+1].y - grain_boundary_pixels[i][j-1].y,
                                        grain_boundary_pixels[i][j+1].x - grain_boundary_pixels[i][j-1].x );
                grain_boundary_phis.back().push_back( phi_temp + PI/2.0f );
            }

            double phi_end = atan2(grain_boundary_pixels[i].back().y - grain_boundary_pixels[i][grain_boundary_pixels[i].size()-2].y,
                                   grain_boundary_pixels[i].back().x - grain_boundary_pixels[i][grain_boundary_pixels[i].size()-2].x );
            grain_boundary_phis.back().push_back( phi_end + PI/2.0f );
            //std::cout << "done!" << std::endl;
        }
        else
        {
            //std::cout << "Long arc... " << grain_boundary_pixels[i].size() << " ";
            calculatePhi( &(grain_boundary_pixels[i]), &(grain_boundary_phis[i]));
            //std::cout << "done!"<<std::endl;
        }

        if(grain_boundary_inverted(grain_boundary_pixels[i]))
            for(int j=0; j<grain_boundary_phis.back().size(); j++)
            {                
                if (grain_boundary_phis.back()[0]>2.0f*PI || grain_boundary_phis.back().back()>2.0f*PI) grain_boundary_phis.back()[j]-=PI;
                else grain_boundary_phis.back()[j]+=PI;
            }
    }

    //calculate the curvature
    for(int i=0; i<grain_boundary_pixels.size(); i++)
    {
        std::vector<float> curv;
        grain_boundary_curvs.push_back(curv);
        calculateCurv( &(grain_boundary_phis[i]), &(grain_boundary_curvs[i]) );
    }
}

void initialize_structures(seg & segment,
                           gbn & grainBoundNet,
                           param & para,
                           std::vector< std::vector< std::vector<int> > > & grain_arc_index,
                           int minimal_bubble_distance,
                           int close_bubble_grain_size,
                           int pixels_average,
                           std::vector< std::vector<point> > & areas,
                           std::vector<area_range> & area_ranges,
                           std::list<int> & found_border_areas,
                           std::vector<int> & arc_state,
                           std::vector< std::vector<int> > & region_arc_index,
                           marray::Marray<unsigned int> & twoCellNeighbors,
                           std::vector<int> & close_bubble_areas,
                           std::vector< std::vector<int> > & grain_area_junctions,
                           std::vector<std::vector<int> > & arc_junctions,
                           std::vector<int> & grain_perimeter,
                           std::vector<int> & min_bubble_distance,
                           std::vector<int> & grain_longest_arc_length,
                           std::vector< std::vector<float> > & grain_junction_angles,
                           std::vector< std::vector<float> > & grain_junction_angles2,
                           std::vector< std::vector<point> > & grain_boundary_pixels,
                           vigra::BasicImage<unsigned int> & region_image,
                           std::vector<bool> & arc_to_segment)   
{
    vigra::BasicImage<unsigned int> & ws_region_image = segment.ws_region_image;
    marray::Marray<unsigned int> & one_boundings = segment.one_boundings;
    marray::Marray<unsigned int> & two_boundings = segment.two_boundings;
    std::vector< std::vector<point> > & arcs = segment.arcs;
    std::vector<point> & junctions = segment.junctions;
    int & dim_x = segment.dim_x;
    int & dim_y = segment.dim_y;

    size_t & nr_areas = grainBoundNet.nr_new_areas;        
    long * & bubble_area_size = grainBoundNet.bubble_area_size;
    long * & grain_area_size = grainBoundNet.grain_area_size;
    std::vector< std::vector<int> > bubble_arc_index = grainBoundNet.bubble_arc_index;
    std::vector<size_t> & region_labels = grainBoundNet.region_labels;
    std::vector<int> & grain_junctions = grainBoundNet.grain_junctions;
    std::vector<int> & grain_bubble_junctions = grainBoundNet.grain_bubble_junctions;
    std::vector<point> & grain_area_center_mass = grainBoundNet.grain_area_center_mass;
    std::vector<point> & bubble_area_center_mass = grainBoundNet.bubble_area_center_mass;
    std::vector<bool> & grain_arc = grainBoundNet.grain_arc;
    std::vector<int> & found_bubble_areas = grainBoundNet.found_bubble_areas;
    
    std::vector< std::vector<int> > & grain_boundary_index = para.grain_boundary_index;
    std::vector< std::vector<float> > & grain_boundary_phis = para.grain_boundary_phis;
    std::vector< std::vector<float> > & grain_boundary_curvs = para.grain_boundary_curvs;
    std::vector< std::vector<unsigned int> > & grain_area_boundaries = para.grain_area_boundaries;
    std::vector< std::vector<int> > & grain_junction_index = para.grain_junction_index;
    std::vector<bool> & grain_junction = para.grain_junction;

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
            region_image(x,y)=area;

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

    //fill vector of found bubble areas/vector of found border areas and region arc index
    for(int area=0; area<nr_areas; area++)
    {
        if (bubble_area_size[area]>0)//bubble area
        {
            region_arc_index[area]=bubble_arc_index[area];
        }
        else if (grain_area_size[area]>0)//grain area
        {
            region_arc_index[area]=grain_arc_index[area][0];
        }
        else found_border_areas.push_back(area+1);//outside
    }

    //adjust twoCellNeighbors
    size_t size[] = {arcs.size(),2};
    twoCellNeighbors.resize(size,size+2);

    for(size_t arc=0;arc<arcs.size();arc++)
    {
        //set arc state: grain arc can be merged, all others not
        if (grain_arc[arc]) arc_state[arc]=1;
        else arc_state[arc]=2;

        //FIND THE INDEX OF THE NEIGHBOR REGIONS (2)
        size_t n0=two_boundings(arc,0);
        size_t n1=two_boundings(arc,1);

        twoCellNeighbors(arc,0)=region_labels[n0-1];
        twoCellNeighbors(arc,1)=region_labels[n1-1];
    }

    //erase wrong bubble center of mass positions
    int i=0;
    while (i<bubble_area_center_mass.size())
    {
        if(bubble_area_center_mass[i].x==0 && bubble_area_center_mass[i].y==0)
        {
            bubble_area_center_mass.erase(bubble_area_center_mass.begin()+i);
        }
        else i++;
    }

   //all junctions for all arcs
    for(int y=0;y<(int)one_boundings.shape(0);y++)
    {
        for(int x=0;x<(int)one_boundings.shape(1);x++)
        {
            if (one_boundings(y,x)!=0) arc_junctions[one_boundings(y,x)-1].push_back(y);//start and end junction of every arc
        }
    }

    grain_area_boundaries.resize(nr_areas);

    junctions_and_outer_circle(nr_areas,
                               grain_area_size,
                               grain_arc_index,
                               grain_junctions,
                               grain_bubble_junctions,
                               grain_area_center_mass,
                               bubble_area_center_mass,
                               grain_arc,
                               minimal_bubble_distance,
                               close_bubble_grain_size,
                               one_boundings,
                               two_boundings,
                               arcs,
                               junctions,
                               pixels_average,
                               areas,
                               area_ranges,
                               found_bubble_areas,
                               found_border_areas,
                               close_bubble_areas,
                               grain_area_junctions,
                               grain_junction,
                               arc_junctions,
                               grain_perimeter,
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
}

void update_structures(seg & segment,
                       gbn & grainBoundNet,
                       param & para,
                       std::vector< std::vector< std::vector<int> > > & grain_arc_index,
                       int minimal_bubble_distance,
                       int close_bubble_grain_size,
                       int pixels_average,
                       std::vector< std::vector<point> > & areas,
                       std::vector<area_range> & area_ranges,
                       std::list<int> & found_border_areas,
                       std::vector<int> & arc_state,
                       std::vector< std::vector<int> > & region_arc_index,
                       marray::Marray<unsigned int> & twoCellNeighbors,
                       std::vector<int> & close_bubble_areas,
                       std::vector< std::vector<int> > & grain_area_junctions,
                       std::vector<std::vector<int> > & arc_junctions,
                       std::vector<int> & grain_perimeter,
                       std::vector<int> & min_bubble_distance,
                       std::vector<int> & grain_longest_arc_length,
                       std::vector< std::vector<float> > & grain_junction_angles,
                       std::vector< std::vector<float> > & grain_junction_angles2,
                       std::vector< std::vector<point> > & grain_boundary_pixels,
                       vigra::BasicImage<unsigned int> & region_image,
                       std::vector<bool> & arc_to_segment)   
{
    vigra::BasicImage<unsigned int> & ws_region_image = segment.ws_region_image;
    marray::Marray<unsigned int> & one_boundings = segment.one_boundings;
    marray::Marray<unsigned int> & two_boundings = segment.two_boundings;
    std::vector< std::vector<point> > & arcs = segment.arcs;
    std::vector<point> & junctions = segment.junctions;
    int & dim_x = segment.dim_x;
    int & dim_y = segment.dim_y;

    size_t & nr_areas = grainBoundNet.nr_new_areas;        
    long * & bubble_area_size = grainBoundNet.bubble_area_size;
    long * & grain_area_size = grainBoundNet.grain_area_size;
    std::vector< std::vector<int> > bubble_arc_index = grainBoundNet.bubble_arc_index;
    std::vector<size_t> & region_labels = grainBoundNet.region_labels;
    std::vector<int> & grain_junctions = grainBoundNet.grain_junctions;
    std::vector<int> & grain_bubble_junctions = grainBoundNet.grain_bubble_junctions;
    std::vector<point> & grain_area_center_mass = grainBoundNet.grain_area_center_mass;
    std::vector<point> & bubble_area_center_mass = grainBoundNet.bubble_area_center_mass;
    std::vector<bool> & grain_arc = grainBoundNet.grain_arc;
    std::vector<int> & found_bubble_areas = grainBoundNet.found_bubble_areas;
    
    std::vector< std::vector<int> > & grain_boundary_index = para.grain_boundary_index;
    std::vector< std::vector<float> > & grain_boundary_phis = para.grain_boundary_phis;
    std::vector< std::vector<float> > & grain_boundary_curvs = para.grain_boundary_curvs;
    std::vector< std::vector<unsigned int> > & grain_area_boundaries = para.grain_area_boundaries;
    std::vector< std::vector<int> > & grain_junction_index = para.grain_junction_index;
    std::vector<bool> & grain_junction = para.grain_junction;
    
    //update grain_arc
    for(size_t arc=0;arc<arcs.size();arc++)
    {
        if (arc_state[arc]==1) grain_arc[arc]=true;
        else grain_arc[arc]=false;
    }

    //update regions labels from region image
    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
           region_labels[ws_region_image(x,y)-1]=region_image(x,y);
        }
    }

    //update grain and bubble arc index, first copy from region arc index
    grain_arc_index.clear();
    grain_arc_index.resize(nr_areas);
    bubble_arc_index.resize(nr_areas);

    for (int area=0; area<nr_areas; area++)
    {
        std::vector<int> segment_index;
        grain_arc_index[area].push_back(segment_index);//all arcs are added to first segment

        grain_arc_index[area][0]=region_arc_index[area];
        bubble_arc_index[area]=region_arc_index[area];
    }


    //clear grain arc list for bubbles
    for (int a=0; a<found_bubble_areas.size(); a++)
    {
        grain_arc_index[found_bubble_areas[a]-1].clear();
    }

    //clear arc lists for outside label
    for (std::list<int>::iterator it=found_border_areas.begin(); it!=found_border_areas.end(); ++it)
    {
        int area=*it;
        grain_arc_index[area-1].clear();
        bubble_arc_index[area-1].clear();
    }

    //clear bubble arc list for grains
    for (int area=0; area<nr_areas; area++)
    {
        if (grain_arc_index[area].size()>0) bubble_arc_index[area].clear();
    }

    //update grain and bubble area size, grain center of mass, area ranges
    delete grain_area_size;
    delete bubble_area_size;
    grain_area_size = new long[nr_areas];
    bubble_area_size = new long[nr_areas];
    grain_area_center_mass.resize(nr_areas);
    area_ranges.resize(nr_areas);

    for (int area=0; area<nr_areas; area++)
    {
        long grain_area_x_sum=0;
        long grain_area_y_sum=0;
        area_ranges[area].x_low=dim_x;
        area_ranges[area].y_low=dim_y;
        area_ranges[area].x_high=0;
        area_ranges[area].y_high=0;

        if (grain_arc_index[area].size()>0)//grain area
        {
            grain_area_size[area]=areas[area].size();

            for(size_t k=0; k<areas[area].size(); ++k)
            {
                int x=areas[area][k].x;
                int y=areas[area][k].y;
                grain_area_x_sum+=x;
                grain_area_y_sum+=y;

                if(x<area_ranges[area].x_low) area_ranges[area].x_low=x;
                if(y<area_ranges[area].y_low) area_ranges[area].y_low=y;
                if(x>area_ranges[area].x_high) area_ranges[area].x_high=x;
                if(y>area_ranges[area].y_high) area_ranges[area].y_high=y;
            }

            int xx=grain_area_x_sum/grain_area_size[area];
            int yy=grain_area_y_sum/grain_area_size[area];
            grain_area_center_mass[area].x=xx;
            grain_area_center_mass[area].y=yy;
        }
        else grain_area_size[area]=0;

        if (bubble_arc_index[area].size()>0)//bubble area
        {
            bubble_area_size[area]=areas[area].size();

        }
        else bubble_area_size[area]=0;
    }

    //update junctions
    grain_junctions.clear();
    grain_bubble_junctions.clear();
    
    for(int y=0;y<(int)one_boundings.shape(0);y++)//loop over junctions
    {
        int nr_grain_arcs=0;
        int nr_bubble_arcs=0;

        for(int x=0;x<(int)one_boundings.shape(1);x++)
        {
            bool found=false;

            for(size_t area=0;area<nr_areas && !found; area++)
            {
                for (int a=0;a<(int)bubble_arc_index[area].size() && !found; a++)//bubble boundaries
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
                if (grain_arc_index[area].size()>0)
                    for (int a=0;a<(int)grain_arc_index[area][0].size() && !found; a++)//grain boundaries
                    {
                        int arc_index=grain_arc_index[area][0][a]+1;
                        if (one_boundings(y,x)==arc_index)
                        {
                            nr_grain_arcs++;
                            found=true;
                        }
                    }
            }

        }

        //grain junctions to vector
        if (nr_grain_arcs==3 || nr_grain_arcs==4)
        {
            grain_junctions.push_back(y);
        }

        //grain bubble junction to vector
        if (nr_bubble_arcs==2 && (nr_grain_arcs==1 || nr_grain_arcs==2))
        {
            grain_bubble_junctions.push_back(y);
        }

    }

    //resizing and clearing structures describing junctions and outer circle of grain boundaries
    grain_area_junctions.clear();
    grain_area_junctions.resize(nr_areas);
    for(int j=0; j<one_boundings.shape(0); j++) grain_junction[j]=false;
    grain_perimeter.resize(nr_areas);
    min_bubble_distance.resize(nr_areas);
    grain_longest_arc_length.resize(nr_areas);
    grain_junction_angles.resize(nr_areas);
    grain_junction_angles2.resize(nr_areas);
    grain_boundary_pixels.clear();
    grain_boundary_index.clear();
    grain_area_boundaries.clear();
    grain_area_boundaries.resize(nr_areas);
    grain_boundary_phis.clear();
    grain_boundary_curvs.clear();
    arc_to_segment.clear();
    arc_to_segment.resize(arcs.size(),false);
    grain_junction_index.clear();
    grain_junction_index.resize(junctions.size());

    junctions_and_outer_circle(nr_areas,
                               grain_area_size,
                               grain_arc_index,
                               grain_junctions,
                               grain_bubble_junctions,
                               grain_area_center_mass,
                               bubble_area_center_mass,
                               grain_arc,
                               0,
                               0,
                               one_boundings,
                               two_boundings,
                               arcs,
                               junctions,
                               pixels_average,
                               areas,
                               area_ranges,
                               found_bubble_areas,
                               found_border_areas,
                               close_bubble_areas,
                               grain_area_junctions,
                               grain_junction,
                               arc_junctions,
                               grain_perimeter,
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
}
