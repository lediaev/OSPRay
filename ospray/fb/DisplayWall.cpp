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

void DisplayWallPO::Instance::postAccum(Tile &tile)
{
#ifdef OSPRAY_DISPLAYCLUSTER
    if(!dcSocket)
        return;

    int rank=-1;
    //MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    //std::cout << "postAccum rank = " << rank << "\n";

    //dcStreamReset(dcSocket);

    // Crashing
    // DisplayWallPO::Instance::postAccum --> dcStreamSend --> dcStreamSendJpeg --> vectorIiSaIiEE19_M_emplace_back_aux
    //return;

    uint32 colorBuffer[TILE_SIZE*TILE_SIZE];

    size_t pixelID = 0;
    for (size_t iy=0; iy<TILE_SIZE; iy++)
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


    //return;

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
