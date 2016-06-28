// ======================================================================== //
// Copyright 2009-2016 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //




# include "dcStream.h"

// ospray
#include "DisplayWall.h"
#include "FrameBuffer.h"
#include <mpi.h>

#include <string>
#include <iostream>
#include <sstream>
#include <ostream>

namespace ospray {

DisplayWallPO::Instance::Instance(DisplayWallPO *po, FrameBuffer *fb)
    : fb(fb), frameIndex(0), dcSocket(NULL), streamName("")
{
#ifdef OSPRAY_DISPLAYCLUSTER
    const char *hostname = po->getParamString("hostname", "localhost");
    streamName = po->getParamString("streamName", "ospray");

    int rank=-1;
    //MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    std::cerr << "connecting to hostname " << hostname << " for stream name " << streamName << " from rank " << rank << std::endl;

    // connect to DisplayCluster at given hostname.
    dcSocket = dcStreamConnect(hostname);

    if(!dcSocket)
        std::cerr << "could not connect to DisplayCluster at " << hostname << std::endl;
#else
    std::cout << "#osp:dw: display cluster support not compiled in" << std::endl;
#endif
}

/*! helper function to convert float-color into rgba-uint format */
// The definition for this is in FrameBuffer.ih.
inline uint32 cvt_uint32(const float f)
{
    return (uint32)(255.f * clamp(f, 0.f, 1.f));
}

/*! helper function to convert float-color into rgba-uint format */
inline uint32 cvt_uint32(const vec4f &v)
{
    return  (cvt_uint32(v.x) << 0)  |
            (cvt_uint32(v.y) << 8)  |
            (cvt_uint32(v.z) << 16) |
            (cvt_uint32(v.w) << 24);
}

/*! helper function to convert float-color into rgba-uint format */
inline uint32 cvt_uint32(const vec3f &v)
{
  return
    (cvt_uint32(v.x) << 0)  |
    (cvt_uint32(v.y) << 8)  |
    (cvt_uint32(v.z) << 16);
}

void DisplayWallPO::Instance::beginFrame()
{
    frameIndex++;

#ifdef OSPRAY_DISPLAYCLUSTER
    dcStreamSetFrameIndex(frameIndex);
#endif
}

std::string IntegerToString (int i)
{
    std::ostringstream convert;
    convert << i;
    return convert.str();
}

void DisplayWallPO::Instance::postAccum(Tile &tile)
{
#ifdef OSPRAY_DISPLAYCLUSTER
    if(!dcSocket)
        return;

    //dcStreamDisconnect(dcSocket);
    //const char *hostname = "tcg-vis-ivb-00";//getParamString("hostname", "localhost");
    //const char *streamName = getParamString("streamName", "ospray");
    //dcSocket = dcStreamConnect(hostname);

    int rank=-1;
    //MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    //std::cout << "postAccum rank = " << rank << "\n";

    //dcStreamReset(dcSocket);

    // Crashing
    // DisplayWallPO::Instance::postAccum --> dcStreamSend --> dcStreamSendJpeg --> vectorIiSaIiEE19_M_emplace_back_aux
    //return;

    int ID = tile.accumID;
    int imageX = tile.region.lower.x;
    int imageY = tile.fbSize.y-tile.region.upper.y;
    int imageWidth = tile.fbSize.x;
    int imageHeight = tile.fbSize.y;
    int imagePitch = 4*TILE_SIZE;
    PIXEL_FORMAT pixelFormat = RGBA;
    int width = TILE_SIZE;
    int height = TILE_SIZE;
    int sourceIndex = tile.region.lower.y*tile.fbSize.x + tile.region.lower.x;
    //int sourceIndex = (int)((long int)(tile.region.lower.y)* (long int)(tile.fbSize.x) + (long int)(tile.region.lower.x));


    if (0 && ID % 1 == 0) {

        int nx = imageWidth / TILE_SIZE;
        int i_tile = imageX / TILE_SIZE;
        int j_tile = imageY / TILE_SIZE;
        int dig = 1 + i_tile + j_tile * nx;

        //unsigned char ppm[3*TILE_SIZE*TILE_SIZE];

        FILE * i_file;
        std::string file_base = "/nfshome/lmlediae/Desktop/tiles/", file_type = ".ppm", file_name;
        file_name = file_base + "ID" + IntegerToString(ID) + "_tile" + IntegerToString(dig) + "_rank" + IntegerToString(rank) + file_type;
        i_file = fopen(file_name.c_str(), "wb");

        bool flipped = true;
        //bool flipped = false;

        if (i_file) {
            fprintf(i_file,"P6\n%d %d\n255\n",TILE_SIZE,TILE_SIZE);
            size_t pixelID = 0;


            if (flipped) {
                int j_high = 63;
                for (int iy=j_high; iy>=0; iy--) {
                    for (int ix=0; ix<TILE_SIZE; ix++) {
                        pixelID = ix + iy*TILE_SIZE;

                        if (pixelID < 0 || pixelID >= TILE_SIZE*TILE_SIZE) {
                            printf("pixel index out of bounds i=%d, j=%d, id=%d, tile size = %d\n",ix,iy, pixelID, TILE_SIZE);
                            exit(1);
                        }

                        float r = tile.r[pixelID];
                        float g = tile.g[pixelID];
                        float b = tile.b[pixelID];
                        float a = tile.a[pixelID];

                        float gamma = 1.0/2.2;
                        r = powf(r,gamma);
                        g = powf(g,gamma);
                        b = powf(b,gamma);

                        int c_r = r*255;
                        int c_g = g*255;
                        int c_b = b*255;
                        if (c_r < 0) c_r = 0;
                        if (c_r > 255) c_r = 255;
                        if (c_g < 0) c_g = 0;
                        if (c_g > 255) c_g = 255;
                        if (c_b < 0) c_b = 0;
                        if (c_b > 255) c_b = 255;

                        unsigned char d[3] = {c_r, c_g, c_b};

                        fwrite(d,3,1,i_file);
                    }
                }
            } else {
                for (size_t iy=0; iy<TILE_SIZE; iy++) {
                    for (size_t ix=0; ix<TILE_SIZE; ix++) {
                        float r = tile.r[pixelID];
                        float g = tile.g[pixelID];
                        float b = tile.b[pixelID];
                        float a = tile.a[pixelID];

                        float gamma = 1.0/2.2;
                        r = powf(r,gamma);
                        g = powf(g,gamma);
                        b = powf(b,gamma);

                        int c_r = r*255;
                        int c_g = g*255;
                        int c_b = b*255;
                        if (c_r < 0) c_r = 0;
                        if (c_r > 255) c_r = 255;
                        if (c_g < 0) c_g = 0;
                        if (c_g > 255) c_g = 255;
                        if (c_b < 0) c_b = 0;
                        if (c_b > 255) c_b = 255;

                        unsigned char d[3] = {c_r, c_g, c_b};

                        fwrite(d,3,1,i_file);

                        ++pixelID;
                    }
                }
            }

            fclose(i_file);
        } else {
            std::cout << "Error opening tile iamge file " << file_name << "\n";
        }
    }


    uint32 colorBuffer[TILE_SIZE*TILE_SIZE];

    size_t pixelID = 0;
    for (size_t iy=0; iy<TILE_SIZE; iy++) {
        for (size_t ix=0; ix<TILE_SIZE; ix++) {
            vec4f col = vec4f(tile.r[pixelID],
                              tile.g[pixelID],
                              tile.b[pixelID],
                              tile.a[pixelID]);

            // Float colors, normal range [0.0 1.0].
            float c_r = tile.r[pixelID];
            float c_g = tile.g[pixelID];
            float c_b = tile.b[pixelID];
            float c_a = tile.a[pixelID];

            float gamma = 1.0/2.2;// Correction for RGBA8 display.;
            c_r = powf(c_r,gamma);
            c_g = powf(c_g,gamma);
            c_b = powf(c_b,gamma);

            // For converting to int color.
            c_r *= 255.0;
            c_g *= 255.0;
            c_b *= 255.0;
            c_a *= 255.0;

            if (1) {// Clamp to range [0 255].
                if (c_r < 0.0) c_r = 0.0;
                if (c_r > 255.0) c_r = 255.0;

                if (c_g < 0.0) c_g = 0.0;
                if (c_g > 255.0) c_g = 255.0;

                if (c_b < 0.0) c_b = 0.0;
                if (c_b > 255.0) c_b = 255.0;

                if (c_a < 0.0) c_a = 0.0;
                if (c_a > 255.0) c_a = 255.0;
            }

            /*
            uint32 i_r = (uint32)c_r;
            uint32 i_g = (((uint32)c_g) << 8);
            uint32 i_b = (((uint32)c_b) << 16);
            uint32 i_a = (((uint32)c_a) << 24);
            */

            uint32 i_r =   (uint32)c_r & 0x000000ff;
            uint32 i_g = (((uint32)c_g & 0x000000ff) << 8);
            uint32 i_b = (((uint32)c_b & 0x000000ff) << 16);
            uint32 i_a = (((uint32)c_a & 0x000000ff) << 24);


            // Convert.
            //colorBuffer[pixelID] = cvt_uint32(col);

            uint32 c_packed = i_r | i_g | i_b | i_a;// Pack them together, 1 byte (8 bits) per color component.

            colorBuffer[pixelID] = c_packed;

            // Unpack, for comparison.
            if (0) {
                float c1,c2,c3,c4;
                c1 = (c_packed & 0x000000ff);
                c2 = (c_packed & 0x0000ff00) >> 8;
                c3 = (c_packed & 0x00ff0000) >> 16;
                c4 = (c_packed & 0xff000000) >> 24;
                printf("r[%f, %f] g[%f, %f] b[%f, %f] a[%f, %f]\n",c_r, c1, c_g, c2, c_b, c3, c_a, c4);
            }


            ++pixelID;
        }
    }

    /*
    float colorBufferFloats[4*TILE_SIZE*TILE_SIZE];
    size_t pixelID = 0;
    for (size_t iy=0;iy<TILE_SIZE;iy++)
        for (size_t ix=0;ix<TILE_SIZE;ix++, pixelID++) {
            colorBufferFloats[4*pixelID + 0] = tile.r[pixelID];
            colorBufferFloats[4*pixelID + 1] = tile.g[pixelID];
            colorBufferFloats[4*pixelID + 2] = tile.b[pixelID];
            colorBufferFloats[4*pixelID + 3] = tile.a[pixelID];
        }*/




    //imageX = 128;
    //imageY = 128;

    if (rank == 1) {
        //imageY -= 128;
    }

    if (0) {
    printf("\nimageX = %d, imageY = %d, imageWidth = %d, imageHeight = %d, imagePitch = %d, width=height = %d, sourceIndex = %d\n\n",
              imageX,      imageY,      imageWidth,      imageHeight,      imagePitch,      width,             sourceIndex);
    }
    //DcStreamParameters dcStreamGenerateParameters(             std::string name, int sourceIndex,  int x,  int y, int width, int height, int totalWidth, int totalHeight)
    DcStreamParameters dcStreamParameters = dcStreamGenerateParameters(streamName,     sourceIndex, imageX, imageY,     width,     height,     imageWidth,     imageHeight);




    if (0) {

        int nx = imageWidth / TILE_SIZE;
        int i_tile = imageX / TILE_SIZE;
        int j_tile = imageY / TILE_SIZE;
        int dig = 1 + i_tile + j_tile * nx;
        unsigned long int uc;

        if (1 <= dig && dig <= 16) {
            // Read in buffer from digit image.
            FILE * i_file;
            std::string file_base = "/nfshome/lmlediae/Desktop/digits_out/n", file_type = ".txt", file_name;
            file_name = file_base + IntegerToString(dig) + file_type;
            i_file = fopen(file_name.c_str(), "r");

            if (i_file) {
                for (int k=0; k<TILE_SIZE*TILE_SIZE; k++) {
                    fscanf(i_file,"%lu", &uc);
                    colorBuffer[k] = (uint32)uc;
                }
                fclose(i_file);
            } else {
                std::cout << "Error opening digit file " << file_name << "\n";
            }
        }
    }

    //bool dcStreamSend(DcSocket * socket,  unsigned char * imageBuffer, int imageX, int imageY, int imageWidth, int imagePitch, int imageHeight, PIXEL_FORMAT pixelFormat, DcStreamParameters parameters)
    bool success = dcStreamSend( dcSocket, (unsigned char *)colorBuffer,      imageX,    imageY,     imageWidth,     imagePitch,     imageHeight,              pixelFormat,            dcStreamParameters);

    if(!success) {
        std::cerr << "error sending tile to DisplayCluster, disconnecting." << std::endl;
        if (dcSocket != NULL) {
            dcStreamDisconnect(dcSocket);
            dcSocket = NULL;
        }
    }
#else 
    PRINT(tile.region);
#endif
}


OSP_REGISTER_PIXEL_OP(DisplayWallPO,display_wall);

}
