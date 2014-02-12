/*! \file bubble_grain_structure.h
 * \brief Bubble-grain structure.
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
#include <list>

#include <vigra/localminmax.hxx>

#include "gbn.h"

/*! \struct two_int
 * \brief A structure consisting of two integers.
 */
struct two_int
{
    int index;   /*!< Index*/
    int length;  /*!< Length*/
};

//mean
float get_abs_mean(std::vector<float> values)
{
    float mean=0;
    for(int i=0;i<(int)values.size();i++)
    {
        mean=mean+fabs(values[i]);
    }

    return (mean/values.size());
}

//median
float get_median(std::vector<float> values)
{

    int size=values.size();
    if(size!=1)
    {
        //sort the "values"-vector
        sort( values.begin(), values.end() );

        //f_index is the index of the median
        float f_index=((float)size-1)*0.5f;
        //if we make f_index-floor(f_index) we know how to round
        //the median
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

void find_bubble_arcs(vigra::MultiArray<2,float> const & probability,
                      std::string filepath_to_ws_region_image,
                      std::string path_to_output_folder,
                      seg & segment,
                      gbn & GrainBoundNet,
                      std::string param_file,
                      ParameterFile paramFile,
                      vigra::BasicImage<bool> selection_image,
                      std::string fp_image,
                      std::vector<int> & remove_bubble_arcs)
{
    //Define references to attributes from seg class
    vigra::BasicImage<unsigned int> & ws_region_image =     segment.ws_region_image;
    marray::Marray<unsigned int> & one_boundings =          segment.one_boundings;
    marray::Marray<unsigned int> & two_boundings =          segment.two_boundings;
    std::vector< std::vector<point> > & arcs =              segment.arcs;
    std::vector<point> & junctions =                        segment.junctions;
    int & dim_x =                                           segment.dim_x;
    int & dim_y =                                           segment.dim_y;

    //Define references to attributes from gbn class
    std::vector<bool> & found_bubble_arcs =                 GrainBoundNet.found_bubble_arcs;
    std::vector<int> & found_bubble_areas =                 GrainBoundNet.found_bubble_areas;
    std::vector< std::vector<int> > & bubble_arc_index =    GrainBoundNet.bubble_arc_index;

    std::string filepath_b_image = fp_image;

    //Check if there exists a bubble image
    filepath_b_image.resize(filepath_b_image.size()-3);
    filepath_b_image.append("b.bmp");  
    
    FILE *info_bubbles;
    info_bubbles = fopen(filepath_b_image.c_str(), "rb");

    //Load bubble image
    vigra::BasicImage<bool> bubbles(dim_x, dim_y);

    bool bubble_img = false;

    if(info_bubbles == NULL)
    {
        std::cout << "No bubble image found" << std::endl;
        //Set the bubble image to black then
        for(int y = 0; y < dim_y; y++)
        {
            for(int x = 0; x < dim_x; x++)
            {
                bubbles(x,y) = false;
            }
        }
    }
    else
    {
        std::cout << "Bubble image found" << std::endl;
        vigra::ImageImportInfo info_b_image(filepath_b_image.c_str());
        importImage(info_b_image, destImage(bubbles));
        fclose(info_bubbles);
        bubble_img = true;
    }   

    std::cout<<"Identify bubbles..."<<std::endl;

    Parameter<float> threshold;
    threshold.assign("", "bubble_boundary_threshold", 0.6);
    threshold.load(paramFile,"config");

    std::cout << "Parameter bubble_boundary_threshold: " << threshold << std::endl;

    Parameter<float> max_bubble_arc_length;
    max_bubble_arc_length.assign("", "max_bubble_arc_length", 500);
    max_bubble_arc_length.load(paramFile,"config");

    std::cout << "Parameter max_bubble_arc_length: " << max_bubble_arc_length << std::endl;

    std::string param_file_name=get_filename(param_file);
    if (param_file_name != "parameters.txt") param_file_name.resize(param_file_name.size()-4);
    else param_file_name = "";

    std::vector<float> bubble_prob;//of areas
    std::vector<int> bubble_arcs;//arc size to calculate weighted average probablity
    std::vector<float> grain_prob;//of areas
    std::vector<int> grain_arcs;//arc size to calculate weighted average probablity

    int nr_areas=0;

    //determine nr of areas
    //loop over all arcs
    for(int y=0;y<(int)two_boundings.shape(0);y++)
    {
        for(int x=0;x<(int)two_boundings.shape(1);x++)
        {
            if (two_boundings(y,x)>nr_areas) nr_areas=two_boundings(y,x);
        }
    }

    bubble_prob.resize(nr_areas+1);
    bubble_arcs.resize(nr_areas+1);
    grain_prob.resize(nr_areas+1);
    grain_arcs.resize(nr_areas+1);
    bubble_arc_index.resize(nr_areas+1);
    
    //loop over all arcs
    for(int a=0;a<(int)arcs.size();a++)
    {
        //sum bubble/grain propability for the areas
        bubble_prob[two_boundings(a,0)]+=arcs[a].size()*probability(a,2);
        bubble_prob[two_boundings(a,1)]+=arcs[a].size()*probability(a,2);
        bubble_arcs[two_boundings(a,0)]+=arcs[a].size();
        bubble_arcs[two_boundings(a,1)]+=arcs[a].size();

        grain_prob[two_boundings(a,0)]+=arcs[a].size()*std::max(probability(a,1),probability(a,2));
        grain_prob[two_boundings(a,1)]+=arcs[a].size()*std::max(probability(a,1),probability(a,2));
        grain_arcs[two_boundings(a,0)]+=arcs[a].size();
        grain_arcs[two_boundings(a,1)]+=arcs[a].size();
    }

    for (int i=1;i<=nr_areas;i++)
    {
        bubble_prob[i]=bubble_prob[i]/(float)bubble_arcs[i];
        grain_prob[i]=grain_prob[i]/(float)grain_arcs[i];
        //std::cout<<bubble_prob[i]<<" "<<grain_prob[i]<<std::endl;
    }

    if(bubble_img)
    {
        std::vector<long> marked_pixels(nr_areas+1,0);
        std::vector<long> all_pixels(nr_areas+1,0);

        for(int y=0; y<dim_y; y++)
        {
            for(int x=0; x<dim_x; x++)
            {
                if(bubbles(x,y) == true) marked_pixels[ws_region_image(x,y)]++;
                all_pixels[ws_region_image(x,y)]++;
            }
        }
        
        FILE *info_main;
        info_main = fopen(fp_image.c_str(), "rb");

        if(info_main == NULL)
        {
            std::cout << "Error: Original image not found!" << std::endl;

            for (int area=1; area<=nr_areas; area++)
            {
                if((double)marked_pixels[area]/(double)all_pixels[area]>0.9)
                {
                    bubble_prob[area]=1.0f;
                }
                else bubble_prob[area]=0.0f;
            }
        }
        else
        {
            fclose(info_main);

            //Load image
            vigra::BImage main_image(dim_x, dim_y);
            vigra::ImageImportInfo info_org_image(fp_image.c_str());
            importImage(info_org_image, destImage(main_image));

            for (int area=1; area<=nr_areas; area++)
            {
                if((double)marked_pixels[area]/(double)all_pixels[area]>0.9)
                {
                    //calculate cross-section difference
                    std::vector<float> cross_section_diff;
                    int found_arcs=0;
                    int arc_index=-1;

                    for(int i=0; i<arcs.size(); i++)//loop over arcs, find bubble boundaries
                    {
                        if (two_boundings(i,0)==area || two_boundings(i,1)==area)//arc belongs to area
                        {
                            found_arcs++;
                            arc_index=i;

                            // calculate angle phi of local normal vectors
                            std::vector<float> phi;
                            if(arcs[i].size() < 5 )
                            {
                                // if the size of an arc is smaller than 5, use a simpler approximation of the derivative
                                double phi0 = atan2(arcs[i][1].y - arcs[i][0].y, arcs[i][1].x - arcs[i][0].x );
                                phi.push_back( phi0 + PI/2.0f );
                 
                                for(int j=1; j<arcs[i].size()-1; j++)
                                {
                                    double phi_temp = atan2(arcs[i][j+1].y - arcs[i][j-1].y, arcs[i][j+1].x - arcs[i][j-1].x );
                                    phi.push_back( phi_temp + PI/2.0f );
                                }
                 
                                double phi_end = atan2(arcs[i].back().y - arcs[i][arcs[i].size()-2].y, arcs[i].back().x - arcs[i][arcs[i].size()-2].x );
                                phi.push_back( phi_end + PI/2.0f );
                            }
                            else
                            {
                                calculatePhi( &(arcs[i]), &(phi));
                            }

                            // loop through the pixels of the arc
                            for(int j=0; j<arcs[i].size(); j++)
                            {
                                // sample the cross section using information on the local normal
                                std::vector<float> crossSection;
                                for(int k=-10; k<11; k++)
                                {
                                    int x =  (int)(arcs[i][j].x + cos(phi[j])*k);
                                    int y =  (int)(arcs[i][j].y + sin(phi[j])*k);

                                    if( 0 <= x && main_image.width() > x + 1 && 0 <= y && main_image.height() > y + 1 )
                                        crossSection.push_back((float)main_image(x,y));
                                }	
                                if(crossSection.size()>0) cross_section_diff.push_back(fabs(crossSection[0]-crossSection.back()));
                            }
                        }
                    }

                    if(get_median(cross_section_diff)>10 || found_arcs>1)//do not remove connected bubbles
                    {
                        bubble_prob[area]=1.0f;
                    }
                    else
                    {
                        bubble_prob[area]=0.0f;
                        remove_bubble_arcs[arc_index]=true;
                    }
                }
                else bubble_prob[area]=0.0f;
            }
        }

        threshold=0.40f;
        std::cout<<"lower bubble_boundary_threshold to: "<<threshold<<std::endl;
    }

    std::vector<std::vector<int> > arc_junctions;//all junctions for all arcs within selected region
    arc_junctions.resize(two_boundings.shape(0));
    
    for(int y=0;y<(int)one_boundings.shape(0);y++)
    {
        for(int x=0;x<(int)one_boundings.shape(1);x++)
        {
            if (one_boundings(y,x)!=0 && selection_image(junctions[y].x,junctions[y].y))
                arc_junctions[one_boundings(y,x)-1].push_back(y);//start and end junction of every arc
        }
    }

    //loop over areas, store arcs of classified area to arc lists
    for (int area=1;area<=nr_areas;area++)//loop over areas
    {
        if ( (bubble_prob[area]>threshold) && ((bubble_prob[area]+0.2)>grain_prob[area]) && (bubble_arcs[area]<max_bubble_arc_length || bubble_img))//area is bubble
        {
            for(int a=0;a<(int)arcs.size();a++)//loop over arcs, find bubble boundaries
            {
                if (two_boundings(a,0)==area || two_boundings(a,1)==area)//arc belongs to area
                {
                    if (((bubble_prob[two_boundings(a,0)]>threshold) && (grain_prob[two_boundings(a,1)]>threshold) &&
                        ((bubble_prob[two_boundings(a,0)]+0.2)>grain_prob[two_boundings(a,0)]) && (probability(a,2)>threshold) ) ||
                        ((bubble_prob[two_boundings(a,1)]>threshold) && (grain_prob[two_boundings(a,0)]>threshold) &&
                        ((bubble_prob[two_boundings(a,1)]+0.2)>grain_prob[two_boundings(a,1)]) && (probability(a,2)>threshold) )) //arc is bubble arc
                    {
                        if (selection_image(arcs[a][0].x,arcs[a][0].y) || selection_image(arcs[a].back().x,arcs[a].back().y))
                        {//at least one arc endpoint is within selected region
                            bubble_arc_index[area].push_back(a);
                        }
                    }

                    //check for bubbles with "border junctions"
                    if (arc_junctions[a].size()==1)//delete all bubbles with connection to the image border
                    {
                        bubble_arc_index[area].clear();
                        //std::cout<<"delete bubble area "<<area<<std::endl;
                        break;
                    }
                }
            }
        }
    }

    for(int y=0;y<(int)two_boundings.shape(0);y++)
    {
        if (arc_junctions[y].size()==1) arc_junctions[y].push_back(one_boundings.shape(0));//"border junction" has label max_label+1, isolated bubbles have no junction
    }

    std::cout<<"start dijkstra bubble closing"<<std::endl;

    //check that all bubble arc lists are closed, otherwise close them
    for (int area=1;area<=nr_areas;area++)//loop over areas
    {
        std::vector<int> bubble_arc_junctions;//junctions for this area
        //std::cout<<"area: "<<area<<" nr arcs: "<<bubble_arc_index[area].size()<<std::endl;
        //loop over arcs and write arc junctions to vector
        for (int a=0;a<(int)bubble_arc_index[area].size();a++)
        if (arc_junctions[bubble_arc_index[area][a]].size()==2)
        {
            bubble_arc_junctions.push_back(arc_junctions[bubble_arc_index[area][a]][0]);
            bubble_arc_junctions.push_back(arc_junctions[bubble_arc_index[area][a]][1]);
        }

        if ((bubble_arc_index[area].size()>0 && bubble_arc_index[area].size()<3) ||
            bubble_arc_junctions.empty())//isolated bubble/area with only one or two arcs classified -> junctions can be close -> dijkstra closing would fail
        {
            if (!bubble_arc_junctions.empty())//not isolated bubble
            {
                for(int a=0;a<(int)arcs.size();a++)//loop over arcs, select all arcs bordering to this area as bubble boundaries
                {
                    if ((two_boundings(a,0)==area || two_boundings(a,1)==area) && arc_junctions[a].size()==2)//arc belongs to area and is not isolated
                    {
                        if (arc_junctions[a][0]<one_boundings.shape(0) && arc_junctions[a][1]<one_boundings.shape(0))
                        {
                            bubble_arc_index[area].push_back(a);
                        }
                    }
                }
            }

            for (int a=0;a<(int)bubble_arc_index[area].size();a++)//bubble boundaries
            {
                int arc_index=bubble_arc_index[area][a];
                found_bubble_arcs[arc_index]=true;
            }

            if (bubble_arc_index[area].size()>0) found_bubble_areas.push_back(area);
        }
        else if (bubble_arc_index[area].size()>2)//do dijkstra bubble closing
        {
            int max_distance=0;
            int end_junction=bubble_arc_junctions.size()-1;

            point start_junction;
            start_junction.x=junctions[bubble_arc_junctions.front()].x;
            start_junction.y=junctions[bubble_arc_junctions.front()].y;

            for (int j=1; j<bubble_arc_junctions.size(); j++)
            {
                int new_distance=sqrt(fabs(junctions[bubble_arc_junctions[j]].x-start_junction.x)*fabs(junctions[bubble_arc_junctions[j]].x-
                    start_junction.x)+ fabs(junctions[bubble_arc_junctions[j]].y-start_junction.y)*fabs(junctions[bubble_arc_junctions[j]].y-
                    start_junction.y));
                if(new_distance>max_distance)
                {
                    max_distance=new_distance;
                    end_junction=j;
                }
            }

            //exchange end_junction(max distance) with last one;
            int temp=bubble_arc_junctions.back();
            bubble_arc_junctions.back()=bubble_arc_junctions[end_junction];
            bubble_arc_junctions[end_junction]=temp;

            int start=bubble_arc_junctions.front();//start junction
            int end=bubble_arc_junctions.back();//end junction
            bubble_arc_junctions.clear();

            std::vector<two_int> arc_length;
            arc_length.resize(two_boundings.shape(0));

            //loop over all arcs
            for(int y=0;y<(int)two_boundings.shape(0);y++)
            {
                two_int arc_temp;
                arc_temp.index=y;//index
                arc_temp.length=arcs[y].size();//length

                //If there exists a bubble image, check if the arc's midpoint is within a certain neighborhood
                //to the actual bubble area. Set the arc's length to zero if that is the case.
                if(bubble_img == true)
                {
                    int xx = arcs[y][arcs[y].size()/2].x;
                    int yy = arcs[y][arcs[y].size()/2].y;

                    int a_loop = 0;
                    int b_loop = 0;

                    if(arcs[y].size()/2 > 2)
                    {
                        a_loop = 2;
                        b_loop = 2;
                    }

                    for(int a = -a_loop; a <= a_loop && (xx+a >= 0 && xx+a < dim_x); a++)
                    {
                        for(int b = -b_loop; b <= b_loop && (yy + b >= 0 && yy+b < dim_y); b++)
                        {
                            if(bubbles(xx+a, yy+b) == true && ws_region_image(xx+a, yy+b)==area)
                            {
                                arc_temp.length = 0;
                            }
                        }
                    }              
                }

                for (int a=0;a<(int)bubble_arc_index[area].size();a++)
                {
                    if (bubble_arc_index[area][a]==y) arc_temp.length=0;//classified arcs have length 0
                }

                arc_length[y]=arc_temp;
            }

            //find dijkstra from start to end 
            std::vector<int> distance;
            std::vector<int> last_junction;
            std::vector<int> last_arc;
            std::list<int> not_visited_junctions;

            distance.resize(one_boundings.shape(0)+1,10000);//distance from start to all junctions, 10000 is infinity
            last_junction.resize(one_boundings.shape(0)+1);//opposite junction of last arc that updated the distance
            last_arc.resize(one_boundings.shape(0)+1);//last arc that updated the distance

            distance[start]=0;
            for (int i=0;i<one_boundings.shape(0)+1;i++) not_visited_junctions.push_back(i);
            bool visited=false;
            bool no_connection=false;

            while(!visited && !no_connection)//stop when end junction is reached or no connection between start and end junction is found
            {
                //find not visited junction with minimal distance
                int min=10000;
                int visit_junction=-1;
                for (std::list<int>::iterator min_junc=not_visited_junctions.begin();min_junc!=not_visited_junctions.end();++min_junc)
                {
                    int min_junction=*min_junc;
                    if (distance[min_junction]<min)
                    {
                        min=distance[min_junction];
                        visit_junction=min_junction;
                    }
                }
                if(visit_junction>-1)//if another unreached junction is found
                {
                    not_visited_junctions.remove(visit_junction);
                    if (visit_junction==end) visited=true;

                    if (visit_junction==one_boundings.shape(0))//junction is border
                    {
                        for(int y=0;y<(int)two_boundings.shape(0);y++)
                        {
                            if (arc_junctions[y].size()==2)
                            {
                                if (arc_junctions[y][1]==visit_junction &&
                                (distance[arc_junctions[y][0]]>distance[arc_junctions[y][1]]+arc_length[y].length))//junction 0 can be updated
                                {
                                    distance[arc_junctions[y][0]]=distance[arc_junctions[y][1]]+arc_length[y].length;
                                    last_junction[arc_junctions[y][0]]=arc_junctions[y][1];
                                    last_arc[arc_junctions[y][0]]=arc_length[y].index;
                                }
                                if (arc_junctions[y][0]==visit_junction &&
                                (distance[arc_junctions[y][1]]>distance[arc_junctions[y][0]]+arc_length[y].length))//junction 1 can be updated
                                {
                                    distance[arc_junctions[y][1]]=distance[arc_junctions[y][0]]+arc_length[y].length;
                                    last_junction[arc_junctions[y][1]]=arc_junctions[y][0];
                                    last_arc[arc_junctions[y][1]]=arc_length[y].index;
                                }
                            }
                        }
                    }
                    else for(int x=0;x<(int)one_boundings.shape(1);x++)
                    {
                        if (one_boundings(visit_junction,x)>0)
                        {
                            if (arc_junctions[one_boundings(visit_junction,x)-1][1]==visit_junction &&
                            (distance[arc_junctions[one_boundings(visit_junction,x)-1][0]]>
                             distance[arc_junctions[one_boundings(visit_junction,x)-1][1]]+arc_length[one_boundings(visit_junction,x)-1].length))//junction 0 can be updated
                            {
                                distance[arc_junctions[one_boundings(visit_junction,x)-1][0]]=distance[arc_junctions[one_boundings(visit_junction,x)-1][1]]
                                    +arc_length[one_boundings(visit_junction,x)-1].length;
                                last_junction[arc_junctions[one_boundings(visit_junction,x)-1][0]]=arc_junctions[one_boundings(visit_junction,x)-1][1];
                                last_arc[arc_junctions[one_boundings(visit_junction,x)-1][0]]=arc_length[one_boundings(visit_junction,x)-1].index;
                            }

                            if (arc_junctions[one_boundings(visit_junction,x)-1][0]==visit_junction &&
                            (distance[arc_junctions[one_boundings(visit_junction,x)-1][1]]>
                             distance[arc_junctions[one_boundings(visit_junction,x)-1][0]]+arc_length[one_boundings(visit_junction,x)-1].length))//junction 1 can be updated
                            {
                                distance[arc_junctions[one_boundings(visit_junction,x)-1][1]]=distance[arc_junctions[one_boundings(visit_junction,x)-1][0]]
                                    +arc_length[one_boundings(visit_junction,x)-1].length;
                                last_junction[arc_junctions[one_boundings(visit_junction,x)-1][1]]=arc_junctions[one_boundings(visit_junction,x)-1][0];
                                last_arc[arc_junctions[one_boundings(visit_junction,x)-1][1]]=arc_length[one_boundings(visit_junction,x)-1].index;
                            }
                        }
                    }
                }
                else no_connection=true;
            }
            //std::cout<<"start "<<start<<", end "<<end<<", distance "<<distance[end]<<std::endl;

            if(!no_connection)
            {
                std::vector<int> all_arcs;

                int check_junction=last_junction[end];
                int check_arc=last_arc[end];
                all_arcs.push_back(check_arc);
                while (check_junction!=start)
                {
                    check_arc=last_arc[check_junction];
                    check_junction=last_junction[check_junction];
                    all_arcs.push_back(check_arc);
                }

                //loop over all arcs
                for(int y=0;y<(int)two_boundings.shape(0);y++)
                {
                    two_int arc_temp;
                    arc_temp.index=y;//index
                    arc_temp.length=arcs[y].size();//length
            
                    //If there exists a bubble image, check if the arc's midpoint is within a certain neighborhood
                    //to the actual bubble area. Set the arc's length to zero if that is the case.
                    if(bubble_img == true)
                    {
                        int xx = arcs[y][arcs[y].size()/2].x;
                        int yy = arcs[y][arcs[y].size()/2].y;

                        int a_loop = 0;
                        int b_loop = 0;

                        if(arcs[y].size()/2 > 2)
                        {
                            a_loop = 2;
                            b_loop = 2;
                        }

                        for(int a = -a_loop; a <= a_loop && (xx+a >= 0 && xx+a < dim_x); a++)
                        {
                            for(int b = -b_loop; b <= b_loop && (yy + b >= 0 && yy+b < dim_y); b++)
                            {
                                if(bubbles(xx+a, yy+b) == true && ws_region_image(xx+a, yy+b)==area)
                                {
                                    arc_temp.length = 0;
                                }
                            }
                        }              
                    }

                    for (int a=0;a<(int)bubble_arc_index[area].size();a++)
                    {
                        if (bubble_arc_index[area][a]==y) arc_temp.length=0;//classified arcs have length 0
                    }

                    for (int b=0;b<(int)all_arcs.size();b++)
                    {
                        if (all_arcs[b]==y) arc_temp.length=10000;//used arcs have length infinity
                    }

                    arc_length[y]=arc_temp;
                }

                //find dijkstra from end to start and set distance of used arcs to infinity
                distance.clear();
                last_junction.clear();
                last_arc.clear();
                distance.resize(one_boundings.shape(0)+1,10000);//distance from start to all junctions, 10000 is infinity
                last_junction.resize(one_boundings.shape(0)+1);//opposite junction of last arc that updated the distance
                last_arc.resize(one_boundings.shape(0)+1);//last arc that updated the distance

                distance[end]=0;

                not_visited_junctions.clear();
                for (int i=0;i<one_boundings.shape(0)+1;i++) not_visited_junctions.push_back(i);
                visited=false;

                while(!visited && !no_connection)//stop when start junction is reached or no connection between end and start junction is found
                {
                    //find not visited junction with minimal distance
                    int min=10000;
                    int visit_junction=-1;
                    for (std::list<int>::iterator min_junc=not_visited_junctions.begin();min_junc!=not_visited_junctions.end();++min_junc)
                    {
                        int min_junction=*min_junc;
                        if (distance[min_junction]<min)
                        {
                            min=distance[min_junction];
                            visit_junction=min_junction;
                        }
                    }
                    if(visit_junction>-1)//if another unreached junction is found
                    {
                        not_visited_junctions.remove(visit_junction);
                        if (visit_junction==start) visited=true;

                        if (visit_junction==one_boundings.shape(0))//junction is border
                        {
                            for(int y=0;y<(int)two_boundings.shape(0);y++)
                            {
                                if (arc_junctions[y].size()==2)
                                {
                                    if (arc_junctions[y][1]==visit_junction &&
                                    (distance[arc_junctions[y][0]]>distance[arc_junctions[y][1]]+arc_length[y].length))//junction 0 can be updated
                                    {
                                        distance[arc_junctions[y][0]]=distance[arc_junctions[y][1]]+arc_length[y].length;
                                        last_junction[arc_junctions[y][0]]=arc_junctions[y][1];
                                        last_arc[arc_junctions[y][0]]=arc_length[y].index;
                                    }
                                    if (arc_junctions[y][0]==visit_junction &&
                                    (distance[arc_junctions[y][1]]>distance[arc_junctions[y][0]]+arc_length[y].length))//junction 1 can be updated
                                    {
                                        distance[arc_junctions[y][1]]=distance[arc_junctions[y][0]]+arc_length[y].length;
                                        last_junction[arc_junctions[y][1]]=arc_junctions[y][0];
                                        last_arc[arc_junctions[y][1]]=arc_length[y].index;
                                    }
                                }
                            }
                        }
                        else for(int x=0;x<(int)one_boundings.shape(1);x++)
                        {
                            if (one_boundings(visit_junction,x)>0)
                            {
                                if (arc_junctions[one_boundings(visit_junction,x)-1][1]==visit_junction &&
                                (distance[arc_junctions[one_boundings(visit_junction,x)-1][0]]>
                                 distance[arc_junctions[one_boundings(visit_junction,x)-1][1]]+arc_length[one_boundings(visit_junction,x)-1].length))//junction 0 can be updated
                                {
                                    distance[arc_junctions[one_boundings(visit_junction,x)-1][0]]=distance[arc_junctions[one_boundings(visit_junction,x)-1][1]]
                                        +arc_length[one_boundings(visit_junction,x)-1].length;
                                    last_junction[arc_junctions[one_boundings(visit_junction,x)-1][0]]=arc_junctions[one_boundings(visit_junction,x)-1][1];
                                    last_arc[arc_junctions[one_boundings(visit_junction,x)-1][0]]=arc_length[one_boundings(visit_junction,x)-1].index;
                                }

                                if (arc_junctions[one_boundings(visit_junction,x)-1][0]==visit_junction &&
                                (distance[arc_junctions[one_boundings(visit_junction,x)-1][1]]>
                                 distance[arc_junctions[one_boundings(visit_junction,x)-1][0]]+arc_length[one_boundings(visit_junction,x)-1].length))//junction 1 can be updated
                                {
                                    distance[arc_junctions[one_boundings(visit_junction,x)-1][1]]=distance[arc_junctions[one_boundings(visit_junction,x)-1][0]]
                                        +arc_length[one_boundings(visit_junction,x)-1].length;
                                    last_junction[arc_junctions[one_boundings(visit_junction,x)-1][1]]=arc_junctions[one_boundings(visit_junction,x)-1][0];
                                    last_arc[arc_junctions[one_boundings(visit_junction,x)-1][1]]=arc_length[one_boundings(visit_junction,x)-1].index;
                                }
                            }
                        }
                    }
                    else no_connection=true;
                }
                if(!no_connection)
                {
                    //std::cout<<"end "<<end<<", start "<<start<<", distance "<<distance[start]<<std::endl;

                    check_junction=last_junction[start];
                    check_arc=last_arc[start];
                    all_arcs.push_back(check_arc);
                    while (check_junction!=end)
                    {
                        check_arc=last_arc[check_junction];
                        check_junction=last_junction[check_junction];
                        all_arcs.push_back(check_arc);
                    }

                    //std::cout<<"all arcs: "<<all_arcs.size()<<std::endl;
                    bubble_arc_index[area]=all_arcs;
                }
            }
            if(no_connection)//no connection found!
            {
                std::cout<<"Error: No connection between bubble arcs found!"<<std::endl;
                bubble_arc_index[area].clear();

                for(int a=0;a<(int)arcs.size();a++)//loop over arcs, select all arcs bordering to this area as bubble boundaries
                {
                    if ((two_boundings(a,0)==area || two_boundings(a,1)==area) && arc_junctions[a].size()==2)//arc belongs to area and is not isolated
                    {
                        if (arc_junctions[a][0]<one_boundings.shape(0) && arc_junctions[a][1]<one_boundings.shape(0))
                        {
                            bubble_arc_index[area].push_back(a);
                        }
                        else//check for bubbles arcs with "border junctions"
                        {
                            for (int a=0;a<(int)bubble_arc_index[area].size();a++)
                            {
                                int arc_index=bubble_arc_index[area][a];
                                found_bubble_arcs[arc_index]=false;
                            }
                            bubble_arc_index[area].clear();

                            break;
                        }
                    }
                }
            }

            for (int a=0;a<(int)bubble_arc_index[area].size();a++)//bubble boundaries
            {
                int arc_index=bubble_arc_index[area][a];
                found_bubble_arcs[arc_index]=true;
            }

            if (bubble_arc_index[area].size()>0) found_bubble_areas.push_back(area);
        }
    }

    /*
    //DEBUG
    vigra::FImage bubble_test(dim_x,dim_y);
    for(int a=0;a<(int)arcs.size();a++)
    {
        if (found_bubble_arcs[a])
        {
            for (int i=0; i<arcs[a].size(); i++)
            {
                int x=arcs[a][i].x;
                int y=arcs[a][i].y;
                bubble_test(x,y)=2;
            }
        }
    }

    std::vector< std::vector<point> > areas(nr_areas);

    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
            int area=ws_region_image(x,y);

            //fill areas vector
            point p;
            p.x=x;
            p.y=y;
            areas[area-1].push_back(p);
        }
    }

    for(int i=0;i<found_bubble_areas.size();i++)
    {
        for(int p=0; p<areas[found_bubble_areas[i]-1].size(); p++)
        {
            int x=areas[found_bubble_areas[i]-1][p].x;
            int y=areas[found_bubble_areas[i]-1][p].y;
            bubble_test(x,y)=1;
        }
    }

    std::string filepath_bubbles=path_to_output_folder;
    filepath_bubbles.append("bubbles/");
    filepath_bubbles.append(param_file_name.c_str());
    filepath_bubbles.append(get_filename(filepath_to_ws_region_image));
    filepath_bubbles.resize(filepath_bubbles.size()-3);
    filepath_bubbles.append("ini.jpg");

    exportImage(srcImageRange(bubble_test), vigra::ImageExportInfo(filepath_bubbles.c_str()));
    */

    std::cout<<"...done"<<std::endl;
}

//merge areas separated by arc of state 0
void merge_areas(vigra::BasicImage<unsigned int> & region_image,
                 std::vector<std::vector<point> > arcs,
                 int dim_x,
                 int dim_y,
                 std::vector<int> & state,
                 marray::Marray<unsigned int> & twoCellNeighbors,
                 size_t & nr_of_regions,
                 std::vector< std::vector<point> > & areas,
                 std::vector< std::vector<int> > & region_arc_index,
                 int minimal_grain_size,
                 std::vector<int> & found_bubble_areas,
                 std::list<int> & found_border_areas,
                 std::vector<int> & areas_to_merge,
                 vigra::MultiArray<2, double> unknown_probability,
                 int version,
                 std::string fp_image = "")
{    
    // used_areas is used to avoid one area being used in merging several times which could lead to merging of a large area
    std::vector<bool> used_areas(nr_of_regions+1,false);

    //do not use bubble areas for merging
    for (int a=0; a<found_bubble_areas.size(); a++)
    {
        used_areas[found_bubble_areas[a]]=true;
    }

    //do not use border areas for merging (with inside areas)
    for (std::list<int>::iterator it=found_border_areas.begin(); it!=found_border_areas.end(); ++it)
    {
        int area=*it;
        used_areas[area]=true;
    }

    //mark certain areas for merging
    std::vector<bool> merge(nr_of_regions,false);

    for (int a=0; a<areas_to_merge.size(); a++)
    {
        merge[areas_to_merge[a]]=true;
    }

    //images required for version 5 is loaded once
    vigra::FImage original_image;

    if(version == 5)
    {
        if(fp_image.compare("") == 0)
        {
            std::cout << "The filepath to the orginal image must be specified for version 5!" << std::endl;
            std::cout << "Using version 2 instead"<< std::endl;
            version=2;
        }
        else
        {
            FILE *info_orig_fp;
            info_orig_fp = fopen(fp_image.c_str(), "rb");

            if(info_orig_fp == NULL)
            {
                std::cout << "Could not find the original image" << std::endl;
                std::cout << "Using version 2 instead"<< std::endl;
                version=2;
            }
            else
            {
                //Load the original image
                vigra::ImageImportInfo info_orig(fp_image.c_str());
                original_image.resize(dim_x, dim_y);
                importImage(info_orig, destImage(original_image));
            }
        }
    }

    for(int area=1; area<nr_of_regions+1; area++)
    {
        //select grain areas with size smaller than minimal grain size, areas with size 1 are marked for merging
        if ((areas[area-1].size()<minimal_grain_size || merge[area-1])  && !used_areas[area])
        {
            merge[area-1]=false;

            //VERSION 1:LONGEST ARC
            if (version==1)
            {
                int max_arc_length=0;
                int longest_arc=-1;

                for(int a=0; a<region_arc_index[area-1].size(); a++)
                    if (arcs[region_arc_index[area-1][a]].size()>max_arc_length && state[region_arc_index[area-1][a]]==1)
                    {
                        longest_arc=region_arc_index[area-1][a];
                        max_arc_length=arcs[region_arc_index[area-1][a]].size();
                    }

                //merge with the area of the longest arc
                if (longest_arc>-1)
                {
                    state[longest_arc]=0;//state 0 will be merged now
                    used_areas[twoCellNeighbors(longest_arc,0)]=true;
                    used_areas[twoCellNeighbors(longest_arc,1)]=true;
                }
            }

            //VERSION 2:LONGEST SEGMENT
            if (version==2)
            {
                std::vector<int> area_arc_length;
                std::vector<int> area_arc_index;
                area_arc_length.resize(nr_of_regions,0);
                area_arc_index.resize(nr_of_regions,0);

                int max_area=-1;
                int max_length=0;

                for(int a=0; a<region_arc_index[area-1].size(); a++)
                    if (state[region_arc_index[area-1][a]]==1)
                    {
                        if (twoCellNeighbors(region_arc_index[area-1][a],0)==area)
                        {
                            area_arc_length[twoCellNeighbors(region_arc_index[area-1][a],1)-1]+=arcs[region_arc_index[area-1][a]].size();
                            area_arc_index[twoCellNeighbors(region_arc_index[area-1][a],1)-1]=region_arc_index[area-1][a];
                        }
                        else
                        {
                            area_arc_length[twoCellNeighbors(region_arc_index[area-1][a],0)-1]+=arcs[region_arc_index[area-1][a]].size();
                            area_arc_index[twoCellNeighbors(region_arc_index[area-1][a],0)-1]=region_arc_index[area-1][a];
                        }
                    }

                //find the area with the longest arc length to current area
                for(int aa=0; aa<nr_of_regions; aa++)
                    if (area_arc_length[aa]>max_length)
                    {
                        max_length=area_arc_length[aa];
                        max_area=aa;
                    }

                //merge with the area of max arc length, take last arc of this area stored in area_arc_index
                if (max_area>-1)
                {
                    state[area_arc_index[max_area]]=0;//state 0 will be merged now
                    used_areas[twoCellNeighbors(area_arc_index[max_area],0)]=true;
                    used_areas[twoCellNeighbors(area_arc_index[max_area],1)]=true;
                }
            }

            //VERSION 3:HIGHEST CURVATURE
            if (version==3)
            {
                float highest_curv;
                int found_arc=-1;

                for(int a=0; a<region_arc_index[area-1].size(); a++)
                    if(state[region_arc_index[area-1][a]]==1)
                    {
                        std::vector<float> this_arc_phi;
                        std::vector<float> curv;

                        //calculate orientation angle phi
                        if(arcs[region_arc_index[area-1][a]].size() < 5 )
                        {
                            // if the size of an arc is smaller than 5, use a simpler approximation of the derivative

                            float phi0 = atan2(arcs[region_arc_index[area-1][a]][1].y - arcs[region_arc_index[area-1][a]][0].y,
                                                arcs[region_arc_index[area-1][a]][1].x - arcs[region_arc_index[area-1][a]][0].x );
                            this_arc_phi.push_back( phi0 + PI/2.0f );

                            for(int j=1; j<arcs[region_arc_index[area-1][a]].size()-1; j++)
                            {
                                float phi_temp = atan2(arcs[region_arc_index[area-1][a]][j+1].y - arcs[region_arc_index[area-1][a]][j-1].y,
                                                        arcs[region_arc_index[area-1][a]][j+1].x - arcs[region_arc_index[area-1][a]][j-1].x );
                                this_arc_phi.push_back( phi_temp + PI/2.0f );
                            }

                            float phi_end = atan2(arcs[region_arc_index[area-1][a]].back().y
                                        - arcs[region_arc_index[area-1][a]][arcs[region_arc_index[area-1][a]].size()-2].y,
                                          arcs[region_arc_index[area-1][a]].back().x
                                        - arcs[region_arc_index[area-1][a]][arcs[region_arc_index[area-1][a]].size()-2].x );
                            this_arc_phi.push_back( phi_end + PI/2.0f );
                        }
                        else
                        {
                            calculatePhi( &(arcs[region_arc_index[area-1][a]]), &(this_arc_phi));
                        }

                        //calculate the curvature
                        calculateCurv( &(this_arc_phi), &(curv) );
                        float mean_curv=get_abs_mean(curv);

                        if (mean_curv>highest_curv || found_arc==-1)
                        {
                            found_arc=region_arc_index[area-1][a];
                            highest_curv=mean_curv;
                        }
                    }

                //merge with the area of the arc of lowest curvature
                if (found_arc>-1)
                {
                    state[found_arc]=0;//state 0 will be merged now
                    used_areas[twoCellNeighbors(found_arc,0)]=true;
                    used_areas[twoCellNeighbors(found_arc,1)]=true;
                }
            }

            //VERSION 4:LOWEST GRAIN PROBABILITY
            if (version==4)
            {
                int lowest_prob=1.0f;
                int most_prob_arc=-1;

                for(int a=0; a<region_arc_index[area-1].size(); a++)
                    if (unknown_probability(region_arc_index[area-1][a],1)<lowest_prob && state[region_arc_index[area-1][a]]==1)
                    {
                        most_prob_arc=region_arc_index[area-1][a];
                        lowest_prob=unknown_probability(region_arc_index[area-1][a],1);
                    }

                //merge with the area of the most probable arc
                if (most_prob_arc>-1)
                {
                    state[most_prob_arc]=0;//state 0 will be merged now
                    used_areas[twoCellNeighbors(most_prob_arc,0)]=true;
                    used_areas[twoCellNeighbors(most_prob_arc,1)]=true;
                }
            }

            //VERSION 5: HIGHEST GRAY VALUE
            if(version == 5)
            {
                int highestGrayVal = 0;
                int highestGrayValArc = -1;

                //Search for the area of the arc with the highest gray value
                for(int a=0; a<region_arc_index[area-1].size(); a++)
                    if (state[region_arc_index[area-1][a]]==1)
                    {
                        for(int j = 0; j < arcs[region_arc_index[area-1][a]].size(); j++)
                        if(original_image(arcs[region_arc_index[area-1][a]][j].x, arcs[region_arc_index[area-1][a]][j].y) > highestGrayVal)
                        {
                            highestGrayVal = original_image(arcs[region_arc_index[area-1][a]][j].x, arcs[region_arc_index[area-1][a]][j].y);
                            highestGrayValArc = region_arc_index[area-1][a];
                        }
                    }
                
                //Merge with the area of the arc with the highest gray value
                if (highestGrayValArc > -1)
                {
                    state[highestGrayValArc]=0;//state 0 will be merged now
                    used_areas[twoCellNeighbors(highestGrayValArc,0)]=true;
                    used_areas[twoCellNeighbors(highestGrayValArc,1)]=true;
                }
            }
        }
    }

    //UFD
    ufd::Partition<size_t> ufd_partition(nr_of_regions+1);

    //NOW WE LOOP OVER ALL STATES (STATES OF THE BOUNDARY´S)
    for(size_t arc=0;arc<state.size();arc++)
    {
        if(state[arc]==0)
        {
           //FIND THE INDEX OF THE NEIGBOUR REGIONS (2)
           size_t n0=twoCellNeighbors(arc,0);
           size_t n1=twoCellNeighbors(arc,1);
           ufd_partition.merge(n0,n1);
        }
    }

    std::map<size_t, size_t> ufd_out;

    //create new labels
    ufd_partition.representativeLabeling(ufd_out);

    std::vector<size_t> region_labels(nr_of_regions);
    for(size_t area=1;area<nr_of_regions+1;area++)
    {
        size_t region_index=ufd_partition.find(area);
        size_t label=ufd_out[ufd_partition.find(area)];
        region_labels[area-1]=label;
    }

    std::vector<size_t>::iterator pos;
    pos = max_element (region_labels.begin(), region_labels.end());
    nr_of_regions=*pos;

    areas.clear();
    areas.resize(nr_of_regions);

    //update bubble areas (nr of bubble areas stays the same)
    for (int a=0; a<found_bubble_areas.size(); a++)
    {
        found_bubble_areas[a]=region_labels[found_bubble_areas[a]-1];
    }

    //update border areas (nr of border areas stays NOT the same)
    for (std::list<int>::iterator it=found_border_areas.begin(); it!=found_border_areas.end(); ++it)
    {
        int area=*it;
        *it=region_labels[area-1];
    }

    found_border_areas.sort();
    found_border_areas.unique();

    //update remaining areas to merge
    areas_to_merge.clear();
    for (int area=0; area<merge.size(); area++)
    {
        if(merge[area]) areas_to_merge.push_back(region_labels[area]-1);
    }

    //update the areas vector and region image
    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
            point p;
            p.x=x;
            p.y=y;
            region_image(x,y)=region_labels[region_image(x,y)-1];
            areas[region_image(x,y)-1].push_back(p);
        }
    }

    //list of arc indeces for new region labels
    region_arc_index.clear();
    region_arc_index.resize(nr_of_regions);

    //update twoCellNeighbors
    size_t size[] = {state.size(),2};
    twoCellNeighbors.resize(size,size+2);

    for(size_t arc=0;arc<state.size();arc++)
    {
        //FIND THE INDEX OF THE NEIGBOUR REGIONS (2)
        size_t n0=twoCellNeighbors(arc,0);
        size_t n1=twoCellNeighbors(arc,1);
        if (region_labels[n0-1]!=region_labels[n1-1])
        {
            region_arc_index[region_labels[n0-1]-1].push_back(arc);
            region_arc_index[region_labels[n1-1]-1].push_back(arc);
        }
        else
        {
            //if subgrains are merged make arcs of same neigbors to subgrain boundaries
            //in any case arc is not added to arc list, state 2 means that arc cannot be selected for merging
            if(version==-1 && state[arc]==1) state[arc]=3; else state[arc]=2;
        }

        twoCellNeighbors(arc,0)=region_labels[n0-1];
        twoCellNeighbors(arc,1)=region_labels[n1-1];
    }
}

//mean
float get_mean(std::vector<float> values)
{
    float mean=0;
    for(int i=0;i<(int)values.size();i++)
    {
        mean=mean+values[i];
    }

    return (mean/values.size());
}

void find_grain_arcs(seg & segment,
                     gbn & GrainBoundNet,
                     std::list<int> & found_border_areas,
                     vigra::MultiArray<2, double> & unknown_probability,
                     std::vector<double> boundary_probability,
                     ParameterFile paramFile,
                     vigra::BasicImage<bool> selection_image,
                     std::string filepath_image,
                     vigra::MultiArray<2, float> & unknown_features,
                     vigra::RandomForest<> rf)
{
    vigra::BasicImage<unsigned int> & ws_region_image = segment.ws_region_image;
    std::vector<size_t> & region_labels = GrainBoundNet.region_labels;
    std::vector<int> & found_bubble_areas = GrainBoundNet.found_bubble_areas;
    std::vector< std::vector<int> > & region_arc_index = GrainBoundNet.grain_arc_index;
    marray::Marray<unsigned int> & twoCellNeighbors = segment.two_boundings;
    std::vector<std::vector<point> > & arcs = segment.arcs;
    int & dim_x = segment.dim_x;
    int & dim_y = segment.dim_y;
/*
    std::string filepath_b_image = filepath_image;

    //Check if there exists a bubble image
    filepath_b_image.resize(filepath_b_image.size()-3);
    filepath_b_image.append("b.bmp");    

    FILE *info_bubbles;
    info_bubbles = fopen(filepath_b_image.c_str(), "rb");

    //Set colours
    //true
    vigra::RGBValue<unsigned int> true;
    true.settrue(255);  

    //yellow
    vigra::RGBValue<unsigned int> yellow;
    yellow.setRed(255);
    yellow.settrue(255);

    //red
    vigra::RGBValue<unsigned int> red;
    red.setRed(255);
    
    //black
    vigra::RGBValue<unsigned int> black;

    vigra::BRGBImage bubbles(dim_x, dim_y);

    bool bubble_img = false;

    if(info_bubbles == NULL)
    {
        std::cout << "No bubble image found" << std::endl;
        //Set the bubble image to black then
        for(int y = 0; y < dim_y; y++)
        {
            for(int x = 0; x < dim_x; x++)
            {
                bubbles(x,y) = 0;
            }
        }
    }
    else
    {
        std::cout << "Bubble image found" << std::endl;
        vigra::ImageImportInfo info_b_image(filepath_b_image.c_str());
        importImage(info_b_image, destImage(bubbles));
        fclose(info_bubbles);
        bubble_img = true;
    }
*/
    std::cout<<"Identify grains..."<<std::endl;
    std::cout<<"find inside of bubbles"<<std::endl;
   
    //find all arcs inside of bubbles by minima growing
    vigra::BasicImage<bool> bubble_arc(dim_x,dim_y);
    bubble_arc = false;
    vigra::BasicImage<bool> bubble_result(dim_x,dim_y);

    for(size_t arc=0; arc<arcs.size(); arc++)
    {
        if (boundary_probability[arc]==1.0f)//arcs classified by bubble search
        {
            for(int p=0; p<arcs[arc].size(); p++)
            {
                bubble_arc(arcs[arc][p].x,arcs[arc][p].y)=true;
            }
        }
    }

    //all bubble pixels become true
    vigra::extendedLocalMinima(vigra::srcImageRange(bubble_arc),vigra::destImage(bubble_result));

    std::vector < std::vector<point> > areas;
    size_t nr_of_regions=0;

    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
           point p;
           p.x=x;
           p.y=y;
           if (ws_region_image(x,y)>nr_of_regions)
           {
               nr_of_regions=ws_region_image(x,y);
               areas.resize(nr_of_regions);
           } 
           areas[ws_region_image(x,y)-1].push_back(p);
        }
    }
/*
    //Set the binary image that resulted from the minima growing to false if there exists a bubble image
    if(bubble_img)
    {
        for(int x = 0; x < dim_x; x++)
        {
            for(int y = 0; y < dim_y; y++)
            {
                bubble_result(x,y) = false;
            }
        }
    }
*/
    region_arc_index.resize(nr_of_regions);
    region_labels.resize(nr_of_regions);

    std::vector<int> arc_state(arcs.size());

    for(size_t arc=0; arc<arcs.size(); arc++)
    {
        point midpoint=arcs[arc][arcs[arc].size()/2];//midpoint of arc

        //arcs classified by image selection, hough labeling or bubble growing
        if (boundary_probability[arc]==0.0f || bubble_result(midpoint.x,midpoint.y))
        {
            arc_state[arc]=0;//will be merged immediately
        }
        else if (boundary_probability[arc]==1.0f)//arcs classified by bubble search
        {
            region_arc_index[twoCellNeighbors(arc,0)-1].push_back(arc);
            region_arc_index[twoCellNeighbors(arc,1)-1].push_back(arc);
            arc_state[arc]=2;//will never be merged
        }
        else
        {
            region_arc_index[twoCellNeighbors(arc,0)-1].push_back(arc);
            region_arc_index[twoCellNeighbors(arc,1)-1].push_back(arc);
            arc_state[arc]=1;//can be selected for merging
        }
    }

    //region image contains current regions and will be updated during merging
    vigra::BasicImage<unsigned int> region_image(dim_x,dim_y);
    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
           region_image(x,y)=ws_region_image(x,y);
        }
    }

    std::vector<int> not_used;
    not_used.clear();

    std::cout<<"merge border, vertical and bubble inside arcs"<<std::endl;

    //we merge initially all arcs of arc state 0, so minimal grain size is zero
    merge_areas(region_image, arcs, dim_x, dim_y, arc_state, twoCellNeighbors, nr_of_regions, areas, region_arc_index, 0,
                found_bubble_areas, found_border_areas, not_used, unknown_probability, 0);


    //Skip if there are no bubbles
    if(!found_bubble_areas.empty())
    {
        std::sort(found_bubble_areas.begin(),found_bubble_areas.end());

        //check if there are still grains inside bubbles
        int bubble_index=0;
        for(int region=0; region<nr_of_regions; region++)
        {
            while (found_bubble_areas[bubble_index]-1<region && bubble_index<found_bubble_areas.size()-1) bubble_index++;

            if(found_bubble_areas[bubble_index]-1>region)//no bubble area with only bubble arcs
            {
                bool only_bubble_arcs=true;

                for(size_t arc=0; arc<region_arc_index[region].size(); arc++)
                {
                    if(boundary_probability[region_arc_index[region][arc]]!=1.0f) only_bubble_arcs=false;
                }

                if(only_bubble_arcs) arc_state[region_arc_index[region][0]]=0;
            }

            if (bubble_index>0) bubble_index--;
        }
    }
    else
    {
        std::cout << "No bubbles found!" << std::endl;
    }

    merge_areas(region_image, arcs, dim_x, dim_y, arc_state, twoCellNeighbors, nr_of_regions, areas, region_arc_index, 0,
                found_bubble_areas, found_border_areas, not_used, unknown_probability, 0);
    std::cout<<"nr of regions after merging: "<<nr_of_regions<<std::endl;

    Parameter<int> minimal_grain_size;
    minimal_grain_size.assign("", "minimal_grain_size", 500);
    minimal_grain_size.load(paramFile,"config");

    std::cout << "Parameter minimal_grain_size: " << minimal_grain_size << std::endl;

    //we merge regions until nr of regions stays the same
    size_t old_nr_of_regions=nr_of_regions;
    size_t diff=1;

    while (diff>0)
    {
        //this is the new version, consistently applied to EDML data set
        merge_areas(region_image, arcs, dim_x, dim_y, arc_state, twoCellNeighbors, nr_of_regions, areas, region_arc_index, minimal_grain_size,
                    found_bubble_areas, found_border_areas, not_used, unknown_probability, 5, filepath_image);

        //this is the old version, consistently applied to NEEM data set
//        merge_areas(region_image, arcs, dim_x, dim_y, arc_state, twoCellNeighbors, nr_of_regions, areas, region_arc_index, minimal_grain_size,
//                    found_bubble_areas, found_border_areas, not_used, unknown_probability, 2);

        std::cout<<"nr of regions after merging: "<<nr_of_regions<<std::endl;
        diff=nr_of_regions-old_nr_of_regions;
        old_nr_of_regions=nr_of_regions;
    }

    Parameter<float> no_boundary_threshold;
    no_boundary_threshold.assign("", "no_boundary_threshold", 1.0f);
    no_boundary_threshold.load(paramFile,"config");

    std::cout << "Parameter no_boundary_threshold: " << no_boundary_threshold << std::endl;
/*
    //boundary region features have to be updated as neighbors might have changed
    if (no_boundary_threshold<=1.0f && filepath_image!="")
    {
        //counts the number of pixels outside image selection
        std::vector<long> outside_pixels(nr_of_regions,0);

        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
            {
               if (!selection_image(x,y)) outside_pixels[region_image(x,y)-1]++;
            }
        }

        //in order to save memory all area mean gray values are computed now instead of doing this for the neighbors of each arc
        std::vector <float> area_mean;
        area_mean.resize(nr_of_regions);
        std::vector <long> area_size;
        area_size.resize(nr_of_regions);

        //original image needs to be loaded
        vigra::ImageImportInfo info(filepath_image.c_str());
        vigra::BImage main_image;
        main_image.resize(dim_x,dim_y);
        importImage(info, destImage(main_image));

        std::cout << "calculating mean gray-values and area size..." << std::endl;
     
        //visit all areas and compute the mean gray-value
        for(int area=0; area<nr_of_regions; area++)
        {
            //vector to store all the pixel values of the area
            std::vector<float> region_pixels;

            //if outside ratio is very low or very high, take all area pixels
            if(outside_pixels[area]<0.1*areas[area].size() || outside_pixels[area]>0.9*areas[area].size() || areas[area].size()-outside_pixels[area]<500)
            { 
                //fill the region_pixels with values
                for(size_t i=0;i<areas[area].size();i++)
                {
                    region_pixels.push_back(   (float)(main_image(areas[area][i].x, areas[area][i].y))     );
                }

                //mean of this area
                area_mean[area]=get_mean(region_pixels);

                //size of this area
                area_size[area]=(long)areas[area].size();
            }
            else//take only inside pixels (which are at least 500)
            {
                //fill the region_pixels with values
                for(size_t i=0;i<areas[area].size();i++)
                {
                    int x=areas[area][i].x;
                    int y=areas[area][i].y;

                    if (selection_image(x,y)) region_pixels.push_back((float)(main_image(x,y)));
                }

                //mean of inside pixels of this area
                area_mean[area]=get_mean(region_pixels);

                //size of inside pixels of this area
                area_size[area]=(long)(areas[area].size()-outside_pixels[area]);
            }
        }

        for(size_t arc=0; arc<arcs.size(); arc++)
        {
            //look in the two boundings which regions are the neigbours
            size_t n0=twoCellNeighbors(arc,0);
            size_t n1=twoCellNeighbors(arc,1);
            if (arc_state[arc]!=1) continue;

            //unknown_features(arc,2)=fabs(area_mean[n0-1]-area_mean[n1-1]);
            unknown_features(arc,3)=area_size[n0-1]+area_size[n1-1];
        }

        std::cout<<"Probability prediction with updated features...."<<std::endl;
        rf.predictProbabilities(unknown_features,unknown_probability);
        std::cout<<"....done"<<std::endl;
*/
    {
        std::vector< std::vector<double> > neighborhood_class0_prob;
        std::vector< std::vector<size_t> > neighborhood_class0_size;
        std::vector< std::vector<int> > neighborhood_class0_arc;

        //we merge regions until nr of regions stays the same
        old_nr_of_regions=nr_of_regions;
        diff=1;
        while (diff>0)
        {
            neighborhood_class0_prob.clear();
            neighborhood_class0_size.clear();
            neighborhood_class0_arc.clear();

            neighborhood_class0_prob.resize(nr_of_regions);
            neighborhood_class0_size.resize(nr_of_regions);
            neighborhood_class0_arc.resize(nr_of_regions);

            for(size_t arc=0; arc<arcs.size(); arc++)
            {
                //average the class0 probability for all segments (combined arcs)
                if (arc_state[arc]==1)
                {
                    //index of neighbor regions
                    size_t n0=std::min(twoCellNeighbors(arc,0),twoCellNeighbors(arc,1))-1;
                    size_t n1=std::max(twoCellNeighbors(arc,0),twoCellNeighbors(arc,1))-1;

                    if(neighborhood_class0_prob[n1].size()<n0+1)
                    {
                        neighborhood_class0_prob[n1].resize(n0+1,0.0f);
                        neighborhood_class0_size[n1].resize(n0+1,0);
                        neighborhood_class0_arc[n1].resize(n0+1,0);
                    }

                    neighborhood_class0_prob[n1][n0]+=arcs[arc].size()*unknown_probability(arc,0);
                    neighborhood_class0_size[n1][n0]+=arcs[arc].size();
                    neighborhood_class0_arc[n1][n0]=arc;
                }
            }

            for(int region=0; region<nr_of_regions; region++)
            {
                for(int n=0; n<(int)neighborhood_class0_prob[region].size(); n++)
                {
                    if(neighborhood_class0_size[region][n]>0)
                        neighborhood_class0_prob[region][n]=neighborhood_class0_prob[region][n]/neighborhood_class0_size[region][n];

                    if(neighborhood_class0_prob[region][n]>no_boundary_threshold)
                    {
                        arc_state[neighborhood_class0_arc[region][n]]=0;//select one arc with these neighbors for merging
                    }
                }
            }

            //merge all arcs of arc state 0, so minimal grain size is zero
            merge_areas(region_image, arcs, dim_x, dim_y, arc_state, twoCellNeighbors, nr_of_regions, areas, region_arc_index, 0,
                        found_bubble_areas, found_border_areas, not_used, unknown_probability, 0);
            std::cout<<"nr of regions after merging: "<<nr_of_regions<<std::endl;
            diff=nr_of_regions-old_nr_of_regions;
            old_nr_of_regions=nr_of_regions;
        }
    }

    //create regions labels from region image
    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
           region_labels[ws_region_image(x,y)-1]=region_image(x,y);
        }
    }

    std::cout<<"...done"<<std::endl;
}

void combine_bubbles_grains(std::string filepath_to_ws_region_image,
                            std::string path_to_output_folder,
                            seg & segment,
                            gbn & GrainBoundNet,
                            std::string param_file,
                            ParameterFile paramFile,
                            vigra::BasicImage<bool> selection_image)
{
    //References to attributes from instance 'GrainBoundNet' of 'gbn'
    size_t & nr_new_areas =                                 GrainBoundNet.nr_new_areas;
    long * bubble_area_size;
    long * grain_area_size;    
    std::vector<int> & grain_junctions =                    GrainBoundNet.grain_junctions;
    std::vector<int> & grain_bubble_junctions =             GrainBoundNet.grain_bubble_junctions;
    std::vector<bool> & grain_arc =                         GrainBoundNet.grain_arc;
    std::vector<point> & grain_area_center_mass =           GrainBoundNet.grain_area_center_mass;
    std::vector<point> & bubble_area_center_mass =          GrainBoundNet.bubble_area_center_mass;
    std::vector<bool> & subgrain_arcs =                     GrainBoundNet.subgrain_arcs;
    std::vector<float> & values =                           GrainBoundNet.values;
    std::vector<size_t> & region_labels =                   GrainBoundNet.region_labels;
    std::vector<int> & found_bubble_areas =                 GrainBoundNet.found_bubble_areas;
    std::vector< std::vector<int> > & bubble_arc_index =    GrainBoundNet.bubble_arc_index;
    std::vector< std::vector<int> > & new_grain_arc_index = GrainBoundNet.grain_arc_index;

    //References to attributes from instance 'segment' of 'seg'
    vigra::BasicImage<unsigned int> & ws_region_image =     segment.ws_region_image;
    marray::Marray<unsigned int> & one_boundings =          segment.one_boundings;
    marray::Marray<unsigned int> & two_boundings =          segment.two_boundings;
    std::vector< std::vector<point> > & arcs =              segment.arcs;
    std::vector<point> & junctions =                        segment.junctions;
    int & dim_x =                                           segment.dim_x;
    int & dim_y =                                           segment.dim_y;

    std::cout<<"Combine bubbles and grains..."<<std::endl;

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

    filepath_bubbles.append(get_filename(filepath_to_ws_region_image));
    filepath_grains.append(get_filename(filepath_to_ws_region_image));
    filepath_new_classification.append(get_filename(filepath_to_ws_region_image));

    //nr of new areas
    std::vector<size_t>::iterator pos;
    pos = max_element (region_labels.begin(), region_labels.end());
    nr_new_areas=*pos;

    //nr of old areas
    size_t nr_old_areas=0;
    std::vector < std::vector<point> > new_areas;
    new_areas.resize(nr_new_areas);

    //check for area outside the selected area
    int outside_label=-1;
    if (!selection_image(0,0)) outside_label=ws_region_image(0,0);
    else if (!selection_image(dim_x-1,0)) outside_label=ws_region_image(dim_x-1,0);
    else if (!selection_image(0,dim_y-1)) outside_label=ws_region_image(0,dim_y-1);
    else if (!selection_image(dim_x-1,dim_y-1)) outside_label=ws_region_image(dim_x-1,dim_y-1);
    if (outside_label>-1)
    {
        outside_label=region_labels[outside_label-1];
    }
    else
    {
        outside_label=region_labels[ws_region_image(0,0)-1];
    }

    //check how detected regions of the network are connected
    vigra::IImage network_region(dim_x,dim_y);
    vigra::IImage region_label_image(dim_x,dim_y);

    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
            network_region(x,y)=1;

            if (region_labels[ws_region_image(x,y)-1]==outside_label)
            {
                network_region(x,y)=0;
            }
        }
    }

    int nr_network_regions = vigra::labelImageWithBackground(vigra::srcImageRange(network_region), vigra::destImage(region_label_image), 0, 0);
    std::vector<std::vector<point> > network_regions(nr_network_regions);

    for (int x=0; x<dim_x; x++)
        for (int y=0; y<dim_y; y++)
        {
            if (region_label_image(x,y)>0)
            {
                point p;
                p.x=x;
                p.y=y;

                network_regions[region_label_image(x,y)-1].push_back(p);
            }
        }

    //find grains of small network regions to exclude in the following
    std::vector<bool> isolated_areas(nr_new_areas,false);

    Parameter<float> minimal_region_size;
    minimal_region_size.assign("", "minimal_region_size", 50000);
    minimal_region_size.load(paramFile,"config");

    std::cout << "Parameter minimal_region_size: " << minimal_region_size << std::endl;
    
    if(!network_regions.empty())
    {
        for(int y=0;y<dim_y;y++)
        {
            for(int x=0;x<dim_x;x++)
            {
                if (network_regions[region_label_image(x,y)-1].size()<minimal_region_size)
                { 
                    isolated_areas[region_labels[ws_region_image(x,y)-1]-1]=true;
                }
            }
        }
    }
    else
    {
        std::cout << "No network regions found!" << std::endl;
    }    
    network_regions.clear();

    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
            point p;
            p.x=x;
            p.y=y;

            if (ws_region_image(x,y)>nr_old_areas)
            {
                nr_old_areas=ws_region_image(x,y);//determine nr of old areas
            } 
            if (region_labels[ws_region_image(x,y)-1]!=outside_label && !isolated_areas[region_labels[ws_region_image(x,y)-1]-1])
            {
                new_areas[region_labels[ws_region_image(x,y)-1]-1].push_back(p);//fill areas vector
            }
        }
    }

    std::vector<bool> new_bubble_labels;
    new_bubble_labels.resize(nr_new_areas,false);

    for(size_t i=0;i<found_bubble_areas.size();i++)
    {
        if (!isolated_areas[found_bubble_areas[i]-1]) new_bubble_labels[found_bubble_areas[i]-1]=true;
    }

    bubble_area_size = new long[nr_new_areas];
    grain_area_size = new long[nr_new_areas];

    std::cout << "calculating area size and center of mass" << std::endl;

    //calculate the center of mass for all grains and bubbles
    //std::vector<point> grain_area_center_mass;
    grain_area_center_mass.resize(nr_new_areas);

    //std::vector<point> bubble_area_center_mass;
    bubble_area_center_mass.resize(nr_new_areas);

    //calculate area size
    for(size_t area=0;area<nr_new_areas;area++)
    {
        long area_x_sum=0;
        long area_y_sum=0;

        grain_area_size[area]=0;
        grain_area_center_mass[area].x=0;
        grain_area_center_mass[area].y=0;

        bubble_area_size[area]=0;
        bubble_area_center_mass[area].x=0;
        bubble_area_center_mass[area].y=0;

        if (area==outside_label-1 || isolated_areas[area])
        {
            new_grain_arc_index[area].clear();
        }
        else if (new_bubble_labels[area]==true)//bubble area
        for(size_t k=0; k<new_areas[area].size(); ++k)
        {
            size_t x=new_areas[area][k].x;
            size_t y=new_areas[area][k].y;
            new_grain_arc_index[area].clear();
            bubble_area_size[area]++;
            area_x_sum+=x;
            area_y_sum+=y;
        }
        else//grain area
        for(size_t k=0; k<new_areas[area].size(); ++k)
        {
            size_t x=new_areas[area][k].x;
            size_t y=new_areas[area][k].y;
            grain_area_size[area]++;
            area_x_sum+=x;
            area_y_sum+=y;
        }

        if (grain_area_size[area]>0)
        {
            int xx=area_x_sum/grain_area_size[area];
            int yy=area_y_sum/grain_area_size[area];
            grain_area_center_mass[area].x=xx;
            grain_area_center_mass[area].y=yy;
        }
        else if (bubble_area_size[area]>0)
        {
            int xx=area_x_sum/bubble_area_size[area];
            int yy=area_y_sum/bubble_area_size[area];
            bubble_area_center_mass[area].x=xx;
            bubble_area_center_mass[area].y=yy;
        }
    }

    std::cout << "creating new arc lists" << std::endl;

    std::vector< std::vector< std::vector<point> > > bubble_arc_list;//all arc points of all arcs of all new bubble areas
    std::vector< std::vector< std::vector<point> > > grain_arc_list;//all arc points of all arcs of all new grain areas

    bubble_arc_list.resize(nr_new_areas);
    grain_arc_list.resize(nr_new_areas);

    std::vector< std::vector<int> > new_bubble_arc_index;
    new_bubble_arc_index.resize(nr_new_areas);

    //Skip if there are no bubbles. 
    if(!bubble_arc_index.empty())
    {
        for(size_t area=1;area<=nr_old_areas;area++)
        {
            for (size_t a=0;a<bubble_arc_index[area].size();a++)
            {
                new_bubble_arc_index[region_labels[area-1]-1].push_back(bubble_arc_index[area][a]);
            }
        }
    }
    else
    {
        std::cout << "No bubbles found!" << std::endl;
    }

    //std::vector<bool> grain_arc;
    grain_arc.resize(arcs.size(),false);

    for(size_t area=0;area<nr_new_areas;area++)
    {
        if (area==outside_label-1 || isolated_areas[area])
        {
            new_bubble_arc_index[area].clear();
        }
        else if (new_bubble_arc_index[area].size()>1)//make entries in bubble arc index unique
        {
            std::sort(new_bubble_arc_index[area].begin(), new_bubble_arc_index[area].end());

            int old_arc_index=new_bubble_arc_index[area][0];
            int i=1;

            while(i<new_bubble_arc_index[area].size())
            {
                if(new_bubble_arc_index[area][i]==old_arc_index)
                {
                    new_bubble_arc_index[area].erase(new_bubble_arc_index[area].begin()+i);
                }
                else
                {
                    old_arc_index=new_bubble_arc_index[area][i];
                    i++;
                }
            }
        }
        for (int a=0;a<(int)new_grain_arc_index[area].size();a++)//grain boundaries
        {
            int arc_index=new_grain_arc_index[area][a];
            grain_arc_list[area].push_back(arcs[arc_index]);
            grain_arc[arc_index]=true;//label arcs as grain arc
        }

        for (int a=0;a<(int)new_bubble_arc_index[area].size();a++)//bubble boundaries
        {
            int arc_index=new_bubble_arc_index[area][a];
            bubble_arc_list[area].push_back(arcs[arc_index]);
            grain_arc[arc_index]=false;//bubble arcs are included in grain arc list but grain arc is false
        }
    }

    std::cout << "find junctions" << std::endl;

    //std::vector<int> grain_junctions;//junctions between 3 or 4 grain arcs
    //std::vector<int> grain_bubble_junctions;//junctions between 1 or 2 grain arc and 2 bubble arcs

    for(int y=0;y<(int)one_boundings.shape(0);y++)//loop over junctions
    {
        int nr_grain_arcs=0;
        int nr_bubble_arcs=0;

        for(int x=0;x<(int)one_boundings.shape(1);x++)
        {
            bool found=false;

            for(size_t area=0;area<nr_new_areas && !found;area++)
            {
                for (int a=0;a<(int)new_bubble_arc_index[area].size() && !found;a++)//bubble boundaries
                {
                    int arc_index=new_bubble_arc_index[area][a]+1;
                    if (one_boundings(y,x)==arc_index)
                    {
                        nr_bubble_arcs++;
                        found=true;
                    }
                }
            }

            for(size_t area=0;area<nr_new_areas && !found;area++)
            {
                for (int a=0;a<(int)new_grain_arc_index[area].size() && !found;a++)//grain boundaries
                {
                    int arc_index=new_grain_arc_index[area][a]+1;
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
        if ((nr_bubble_arcs==2 || nr_bubble_arcs==3) && (nr_grain_arcs==1 || nr_grain_arcs==2))
        {
            grain_bubble_junctions.push_back(y);
        }
    }

    //grain_image: grain areas gray, bubbles and grain boundaries black
    //bubble_image: bubbles areas and bubble boundaries
    vigra::FImage bubble_image(dim_x,dim_y);
    vigra::FImage grain_image(dim_x,dim_y);

    filepath_bubbles.append(".jpg");
    filepath_grains.append(".jpg");

    for(int y=0;y<dim_y;y++)
    {
        for(int x=0;x<dim_x;x++)
        {
            grain_image(x,y)=3;
        }
    }

    for(size_t area=0;area<nr_new_areas;area++)
    {
        if (area==outside_label-1 || isolated_areas[area]) continue;
        else if (new_bubble_labels[area]==true)//bubble area
        for(size_t k=0; k<new_areas[area].size(); ++k)
        {
            size_t x=new_areas[area][k].x;
            size_t y=new_areas[area][k].y;
            bubble_image(x,y)=1;
            grain_image(x,y)=0;
        }
        else//grain area
        for(size_t k=0; k<new_areas[area].size(); ++k)
        {
            size_t x=new_areas[area][k].x;
            size_t y=new_areas[area][k].y;
            grain_image(x,y)=1;
        }

        if (grain_area_size[area]>0)
        {
            int xx=grain_area_center_mass[area].x;
            int yy=grain_area_center_mass[area].y;
            //center of mass positions to grain output image
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
    for (int area=0;area<nr_new_areas;area++)//loop over areas
    {
        for (int a=0;a<(int)bubble_arc_list[area].size();a++)//bubble boundaries
        {
            //now we loop over the points in this arc
            for(int p=0;p<(int)bubble_arc_list[area][a].size();p++)
            {
                int x=bubble_arc_list[area][a][p].x;
                int y=bubble_arc_list[area][a][p].y;
                bubble_image(x,y)=2;
            }
        }

        for (int a=0;a<(int)grain_arc_list[area].size();a++)//grain boundaries
        {
            //now we loop over the points in this arc
            for(int p=0;p<(int)grain_arc_list[area][a].size();p++)
            {
                int x=grain_arc_list[area][a][p].x;
                int y=grain_arc_list[area][a][p].y;
                grain_image(x,y)=0;
            }
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
    GrainBoundNet.bubble_area_size = bubble_area_size;
    GrainBoundNet.grain_area_size = grain_area_size;
    GrainBoundNet.bubble_arc_index = new_bubble_arc_index;

    std::cout<<"export result images"<<std::endl;
    exportImage(srcImageRange(bubble_image), vigra::ImageExportInfo(filepath_bubbles.c_str()));
    exportImage(srcImageRange(grain_image), vigra::ImageExportInfo(filepath_grains.c_str()));
    std::cout<<"...done"<<std::endl;
    GrainBoundNet.save_final_structure(filepath_new_classification);
}
