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

using namespace cimg_library;

struct edge
{
    int dx;
    int dy;
 
    int active;
 
    float corr;
    float mass;
	float coherence;
};
 
struct graphNode
{
    edge up;
    edge down;
    edge right;
    edge left;
 
    int x; //Absolute Position der linken oberen Ecke des Bildes
    int y; //x-Achse l�uft nach rechts, y-Achse nach unten.
 
    int imageWidth;
    int imageHeight;
 
    int loaded;
    int index;
};

struct conEdge
{
    int posr;
    int posc;
 
    int dr;
    int dc;
 
    int dir;
 
    float corr;
    float mass;
};

void connectCoords(graphNode** gData, int i, int j)
{
 
	if(gData[i][j].right.active==1)
        if(gData[i][j+1].index != gData[i][j].index)
        {
            gData[i][j+1].index=gData[i][j].index;
            gData[i][j+1].x=gData[i][j].x + gData[i][j].right.dx;
            gData[i][j+1].y=gData[i][j].y + gData[i][j].right.dy;
            connectCoords(gData,i,j+1);
        }
    if(gData[i][j].down.active==1)
        if(gData[i+1][j].index != gData[i][j].index)
        {
            gData[i+1][j].index=gData[i][j].index;
            gData[i+1][j].x=gData[i][j].x + gData[i][j].down.dx;
            gData[i+1][j].y=gData[i][j].y + gData[i][j].down.dy;
            connectCoords(gData,i+1,j);
        }
    if(gData[i][j].left.active==1)
        if(gData[i][j-1].index != gData[i][j].index)
        {
            gData[i][j-1].index=gData[i][j].index;
            gData[i][j-1].x=gData[i][j].x - gData[i][j].left.dx;
            gData[i][j-1].y=gData[i][j].y - gData[i][j].left.dy;
            connectCoords(gData,i,j-1);
        }
    if(gData[i][j].up.active==1)
        if(gData[i-1][j].index != gData[i][j].index)
        {
            gData[i-1][j].index=gData[i][j].index;
            gData[i-1][j].x=gData[i][j].x - gData[i][j].up.dx;
            gData[i-1][j].y=gData[i][j].y - gData[i][j].up.dy;
            connectCoords(gData,i-1,j);
        }
 
	return;
}

int sortEdges(struct conEdge* edges, int rows, int columns)
{
    int moved=1;
    struct conEdge temp;
 
    while(moved)
    {
        moved=0;
        for(int i=0; i<2*(rows-1)*(columns-1)-1; i++)
        {
			if( edges[i+1].mass*edges[i+1].corr > edges[i].mass*edges[i].corr )
            {
                moved=1;
                temp=edges[i];
                edges[i]=edges[i+1];
                edges[i+1]=temp;
            }
        }
    }
 
    return 0;
 
}
 
int addEdge(graphNode** coords, conEdge edge, int rows, int columns)
{
    if(edge.dir == 0 && coords[edge.posr][edge.posc].right.active != 1)
    {
        if(edge.dc > 1000 || edge.dr > 1000) printf("EDGE KORRUPT!!!!\n");
 
        coords[edge.posr][edge.posc].right.active = 1;
        coords[edge.posr][edge.posc].right.corr = edge.corr;
        coords[edge.posr][edge.posc].right.mass = edge.mass;
 
        coords[edge.posr][edge.posc].right.dx = edge.dc;
        coords[edge.posr][edge.posc].right.dy = edge.dr;
 
        coords[edge.posr][edge.posc+1].left.active = 1;
        coords[edge.posr][edge.posc+1].left.corr = edge.corr;
        coords[edge.posr][edge.posc+1].left.mass = edge.mass;
 
        coords[edge.posr][edge.posc+1].left.dx = edge.dc;
        coords[edge.posr][edge.posc+1].left.dy = edge.dr;
 
         if(coords[edge.posr][edge.posc].right.dx > 1000 || coords[edge.posr][edge.posc].right.dy > 1000)
        {
                printf("ERROR: Zuweisung von korrupter Kante [%d][%d].right: dx = %d, dy = %d !!\n", edge.posr, edge.posc, coords[edge.posr][edge.posc].right.dx, coords[edge.posr][edge.posc].right.dy);
                printf("Korrelation und Masse: %f %f\n", coords[edge.posr][edge.posc].right.corr, coords[edge.posr][edge.posc].right.mass);
                printf("Zugewiesen von: %d %d %f %f\n", edge.dc, edge.dr, edge.corr, edge.mass);
        }
    }
 
    if(edge.dir == 1 && coords[edge.posr][edge.posc].down.active != 1)
    {
        if(edge.dc > 1000 || edge.dr > 1000) printf("EDGE KORRUPT!!!!\n");
 
        coords[edge.posr][edge.posc].down.active = 1;
        coords[edge.posr][edge.posc].down.corr = edge.corr;
        coords[edge.posr][edge.posc].down.mass = edge.mass;
 
        coords[edge.posr][edge.posc].down.dx = edge.dc;
        coords[edge.posr][edge.posc].down.dy = edge.dr;
 
        coords[edge.posr+1][edge.posc].up.active = 1;
        coords[edge.posr+1][edge.posc].up.corr = edge.corr;
        coords[edge.posr+1][edge.posc].up.mass = edge.mass;
 
        coords[edge.posr+1][edge.posc].up.dx = edge.dc;
        coords[edge.posr+1][edge.posc].up.dy = edge.dr;
 
        if(coords[edge.posr][edge.posc].down.dx > 1000 || coords[edge.posr][edge.posc].down.dy > 1000)
        {
                printf("ERROR: Zuweisung von korrupter Kante [%d][%d].down: dx = %d, dy = %d !!\n", edge.posr, edge.posc,  coords[edge.posr][edge.posc].down.dx, coords[edge.posr][edge.posc].down.dy);
                printf("Korrelation und Masse: %f %f\n", coords[edge.posr][edge.posc].down.corr, coords[edge.posr][edge.posc].down.mass);
                printf("Zugewiesen von: %d %d %f %f\n", edge.dc, edge.dr, edge.corr, edge.mass);
        }
    }
 
    for(int i=0; i<rows; i++)
    {
        for(int j=0; j<columns; j++)
        {
            coords[i][j].index=0;
        }
    }
 
    int currIndex=1;
 
    for(int i=0; i<rows; i++)
    {
        for(int j=0; j<columns; j++)
        {
            if(coords[i][j].index==0)
            {
                coords[i][j].index=currIndex;
                connectCoords(coords, i, j);
                currIndex++;
            }
        }
    }
 
    for(int i=0; i<rows-1; i++)
        for(int j=0; j<columns-1; j++)
        {
            if( coords[i][j].index == coords[i][j+1].index && coords[i][j].right.active!=1 )
            {
                coords[i][j].right.active=1;
                coords[i][j].right.dx=coords[i][j+1].x-coords[i][j].x;
                coords[i][j].right.dy=coords[i][j+1].y-coords[i][j].y;
 
                coords[i][j+1].left.active=1;
                coords[i][j+1].left.dx=coords[i][j+1].x-coords[i][j].x;
                coords[i][j+1].left.dy=coords[i][j+1].y-coords[i][j].y;
 
                if(coords[i][j].right.dx > 1000 || coords[i][j].right.dy > 1000)
                    printf("ERROR ERROR [%d][%d].right: dx = %d dy= %d \n\n\n",i,j,coords[i][j].right.dx, coords[i][j].right.dy);
            }
 
            if( coords[i][j].index == coords[i+1][j].index && coords[i][j].down.active!=1 )
            {
                coords[i][j].down.active=1;
                coords[i][j].down.dx=coords[i+1][j].x-coords[i][j].x;
                coords[i][j].down.dy=coords[i+1][j].y-coords[i][j].y;
 
                coords[i+1][j].up.active=1;
                coords[i+1][j].up.dx=coords[i+1][j].x-coords[i][j].x;
                coords[i+1][j].up.dy=coords[i+1][j].y-coords[i][j].y;
 
                if(coords[i][j].down.dx > 1000 || coords[i][j].down.dy > 1000)
                    printf("ERROR ERROR [%d][%d].down: dx = %d dy= %d \n\n\n",i,j,coords[i][j].down.dx, coords[i][j].down.dy);
            }
        }
 
        return 0;
}

void find_mosaic_borders(std::vector< std::vector<point> > arcs, std::vector<int> & vertical_arc_index,std::vector<bool> found_bubble_arcs,
                        marray::Marray<unsigned int> two_boundings, int dim_x, int dim_y, std::string fileName, std::string relativeCoords)
{
    vigra::BasicImage<bool> border_image(dim_x, dim_y);
        for (int x=0;x<dim_x;x++)
            for (int y=0;y<dim_y;y++)
                border_image(x,y)=false;

    std::string filePath=relativeCoords;
    filePath.append(fileName.c_str());
    //remove ".bmp.bin"
    filePath.resize(filePath.size()-8);
    filePath.append(".txt");
    std::ifstream info_file(filePath.c_str());
    std::ifstream temp_info_file(filePath.c_str());

    //string is just for testing stuff
    std::string teststring;
    temp_info_file>>teststring;
    if(info_file && teststring.size()!=0)//file exists and is not empty
    {
        info_file.close();
        temp_info_file.close();
        std::cout<<"Relative coordinates found"<<std::endl;

        int columns, rows;

        FILE* file=fopen(filePath.c_str(), "r");
        fscanf(file, "%d %d",&columns, &rows);

        struct conEdge* edges;
        edges=new struct conEdge[2*(rows-1)*(columns-1)];

        graphNode** coords;
        coords=new graphNode*[rows];
        for(int i=0; i<rows; i++)
        {
            coords[i]=new graphNode[columns];
            for(int j=0; j<columns; j++)
            {
                coords[i][j].up.active=0;
                coords[i][j].down.active=0;
                coords[i][j].right.active=0;
                coords[i][j].left.active=0;
                coords[i][j].imageWidth = 768;
                coords[i][j].imageHeight = 512;
                coords[i][j].loaded=0;
                coords[i][j].index=0;
            }
        }

        coords[0][0].x=0;
        coords[0][0].y=0;

        for(int i=0; i<rows-1; i++)
        {
            for(int j=0; j<columns-1; j++)
            {
                edges[i*(columns-1) + j].dir=0;
                edges[i*(columns-1) + j].posr=i;
                edges[i*(columns-1) + j].posc=j;
                fscanf(file,"%d %d %f %f", &edges[i*(columns-1) + j].dc, &edges[i*(columns-1) + j].dr,
                    &edges[i*(columns-1) + j].corr, &edges[i*(columns-1) + j].mass);

                if( i < 10 && (edges[i*(columns-1) + j].dc < 600 || edges[i*(columns-1) + j].dc > 700) )
                    printf("Error: [%d][%d].hor.dc== %d\n", i, j, edges[i*(columns-1) + j].dc);
                if( i < 10 && (edges[i*(columns-1) + j].dr < -20 || edges[i*(columns-1) + j].dr > 20) )
                    printf("Error: [%d][%d].hor.dr== %d\n", i, j, edges[i*(columns-1) + j].dr);

                edges[i*(columns-1) + j + (rows-1)*(columns-1)].dir=1;
                edges[i*(columns-1) + j + (rows-1)*(columns-1)].posr=i;
                edges[i*(columns-1) + j + (rows-1)*(columns-1)].posc=j;
                fscanf(file,"%d %d %f %f", &edges[i*(columns-1) + j + (rows-1)*(columns-1)].dc,
                    &edges[i*(columns-1) + j + (rows-1)*(columns-1)].dr, &edges[i*(columns-1) + j + (rows-1)*(columns-1)].corr,
                    &edges[i*(columns-1) + j + (rows-1)*(columns-1)].mass);

                if( i < 10 && (edges[i*(columns-1) + j + (rows-1)*(columns-1)].dc < -20 ||
                    edges[i*(columns-1) + j + (rows-1)*(columns-1)].dc > 20) )
                    printf("Error: [%d][%d].vert.dc== %d\n", i, j,
                    edges[i*(columns-1) + j + (rows-1)*(columns-1)].dc);
                if( i < 10 && (edges[i*(columns-1) + j + (rows-1)*(columns-1)].dr < 400 ||
                    edges[i*(columns-1) + j + (rows-1)*(columns-1)].dr > 500) )
                    printf("Error: [%d][%d].vert.dr== %d\n", i, j,
                    edges[i*(columns-1) + j + (rows-1)*(columns-1)].dr);
            }
        }

        sortEdges(edges, rows, columns);

        for(int i=0; i<2*(rows-1)*(columns-1); i++)
        {
            addEdge(coords,edges[i],rows, columns);
        }

        fscanf(file,"%d %d %f %f", &edges[rows*columns-rows-1].dc, &edges[rows*columns-rows-1].dr,
            &edges[rows*columns-rows-1].corr, &edges[rows*columns-rows-1].mass);
        fclose(file);

        coords[rows-1][columns-1].x=coords[rows-1][columns-2].x+edges[rows*columns-rows-1].dc;
        coords[rows-1][columns-1].y=coords[rows-1][columns-2].y+edges[rows*columns-rows-1].dr;

        for(int i=0; i<rows-1; i++)
        {
            for(int j=0; j<columns-1; j++)
            {
                for (int x=std::max(0,coords[i][j].x); x<std::min(coords[i][j+1].x,coords[i][j].x+768); x++)
                {
                    border_image(x,std::max(0,coords[i][j].y))=true;
                }
                for (int y=std::max(0,coords[i][j].y); y<std::min(coords[i+1][j].y,coords[i][j].y+512); y++)
                {
                    border_image(std::max(0,coords[i][j].x),y)=true;
                }
            }

            //last colomn
            for (int x=std::max(0,coords[i][columns-1].x); x<coords[i][columns-1].x+768; x++)
            {
                border_image(x,std::max(0,coords[i][columns-1].y))=true;
            }
            for (int y=std::max(0,coords[i][columns-1].y); y<std::min(coords[i+1][columns-1].y,coords[i][columns-1].y+512); y++)
            {
                border_image(std::max(0,coords[i][columns-1].x),y)=true;
            }
        }

        //last row
        for(int j=0; j<columns-1; j++)
        {
            for (int x=std::max(0,coords[rows-1][j].x); x<std::min(coords[rows-1][j+1].x,coords[rows-1][j].x+768); x++)
            {
                border_image(x,std::max(0,coords[rows-1][j].y))=true;
            }
            for (int y=std::max(0,coords[rows-1][j].y); y<coords[rows-1][j].y+512; y++)
            {
                border_image(std::max(0,coords[rows-1][j].x),y)=true;
            }
        }

        //bottom right
        for (int x=std::max(0,coords[rows-1][columns-1].x); x<coords[rows-1][columns-1].x+768; x++)
        {
            border_image(x,std::max(0,coords[rows-1][columns-1].y))=true;
        }
        for (int y=std::max(0,coords[rows-1][columns-1].y); y<coords[rows-1][columns-1].y+512; y++)
        {
            border_image(std::max(0,coords[rows-1][columns-1].x),y)=true;
        }

        delete [] edges;
 
        for(int i=0; i<rows; i++) delete [] coords[i];
            delete [] coords;

        //Create temp image
        vigra::BasicImage<bool> temp_image(dim_x, dim_y);
        for (int x=0;x<dim_x;x++)
            for (int y=0;y<dim_y;y++)
                temp_image(x,y)=border_image(x,y);

        // Disc Dilation with radius r
        int r=3;
        for(int i=0; i<dim_x; i++)
            for(int j=0; j<dim_y; j++)
                for(int k = -r; k<=r; k++)
                    for(int l= -r; l<=r; l++)
                    {
                        if( k*k + l*l <= r*r )
                        {
                            if( k+i >= 0 && k+i < dim_x && j+l >=0 && j+l < dim_y )
                            {
                                if( temp_image( k+i, j+l ) && !border_image( i, j ) )
                                    border_image( i, j) = temp_image( k+i, j+l );
                            }
                        }
                    }
        /*
        vigra::IRGBImage debug_img(dim_x,dim_y);
        for (int x=0;x<dim_x;x++)
            for (int y=0;y<dim_y;y++)
            {
                //if (border_image(x,y))
                {
                    //debug_img(x,y)[0]=0;
                    //debug_img(x,y)[1]=0;
                    //debug_img(x,y)[2]=0;
                }
                //else
                {
                    debug_img(x,y)[0]=1;
                    debug_img(x,y)[1]=1;
                    debug_img(x,y)[2]=1;
                }
            }
        */

        std::vector<bool> used;

        for(int arc=0; arc<arcs.size(); arc++)
        {
            int max_label=std::max(two_boundings(arc,0),two_boundings(arc,1));
            if (max_label>used.size()) used.resize(max_label+1,false);

            int pixels_on_line=0;

            //calculate fraction at mosaic border
            for(int p=0; p<arcs[arc].size(); p++)
            {
                if (border_image(arcs[arc][p].x,arcs[arc][p].y)) pixels_on_line++;
            }

            if (pixels_on_line>0.7*arcs[arc].size() && !found_bubble_arcs[arc] &&
                //!used[two_boundings(arc,0)] && !used[two_boundings(arc,1)])
                arcs[arc].size()>100)
                //true)
            {
                //std::cout<<"mosaic border arc "<<arc<<" found!"<<std::endl;
                //std::cout<<pixels_on_line<<"/"<<arcs[arc].size()<<std::endl;
                //for(int p=0; p<arcs[arc].size(); p++)
                {
                    //debug_img(arcs[arc][p].x,arcs[arc][p].y)[0]=1;
                    //debug_img(arcs[arc][p].x,arcs[arc][p].y)[1]=0;
                    //debug_img(arcs[arc][p].x,arcs[arc][p].y)[2]=0;
                }

                used[two_boundings(arc,0)]=true;
                used[two_boundings(arc,1)]=true;
                vertical_arc_index.push_back(arc);
            }
        }
    }
    else if (info_file)
    {
        std::cout<<"Relative coordinates file is empty!"<<std::endl;
        info_file.close();
        temp_info_file.close();
    }
    else std::cout<<"Relative coordinates not found!"<<std::endl;
}
