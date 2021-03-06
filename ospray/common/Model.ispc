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

#include "Model.ih"
#include "embree2/rtcore_scene.isph"

typedef uniform Geometry *uniform uniGeomPtr;
typedef uniform Material *uniform uniMaterialPtr;
typedef uniform Volume   *uniform uniVolumePtr;

export void *uniform Model_create(void *uniform cppE)
{
  uniform Model *uniform model = uniform new uniform Model;
  model->cppEquivalent     = cppE;
  model->embreeSceneHandle = NULL;
  model->geometry          = NULL;
  model->volumes           = NULL;
  return (void *uniform)model;
}

export void Model_init(void *uniform _model, 
                       void *uniform embreeDevice,
                       uniform int32 numGeometries, 
                       uniform int32 numVolumes)
{
  uniform Model *uniform model = (uniform Model *uniform)_model;
  if (model->embreeSceneHandle)
    rtcDeleteScene(model->embreeSceneHandle);

  model->embreeSceneHandle = rtcDeviceNewScene((RTCDevice)embreeDevice,
                                               //RTC_SCENE_STATIC|RTC_SCENE_HIGH_QUALITY,
                                               RTC_SCENE_STATIC,//|RTC_SCENE_COMPACT,
                                               //RTC_SCENE_DYNAMIC,
                                               //RTC_SCENE_DYNAMIC|RTC_SCENE_COMPACT,
                                               RTC_INTERSECT_UNIFORM|RTC_INTERSECT_VARYING);
  
  if (model->geometry) delete[] model->geometry;
  model->geometryCount = numGeometries;
  if (numGeometries > 0)
    model->geometry = uniform new uniform uniGeomPtr[numGeometries];
  else 
    model->geometry = NULL;

  if (model->volumes) delete[] model->volumes;
  model->volumeCount = numVolumes;
  if (numVolumes > 0) 
    model->volumes = uniform new uniform uniVolumePtr[numVolumes];
  else
    model->volumes = NULL;
}

export void *uniform Model_getEmbreeSceneHandle(void *uniform _model)
{
  uniform Model *uniform model = (uniform Model *uniform)_model;
  return (void *uniform)model->embreeSceneHandle;
}

export void Model_setGeometry(void *uniform _model,
                              uniform int32 geomID,
                              void *uniform _geom)
{
  uniform Model *uniform model = (uniform Model *uniform)_model;
  uniform Geometry *uniform geom = (uniform Geometry *uniform)_geom;
  model->geometry[geomID] = geom;
}

export void Model_setVolume(void *uniform pointer,
                            uniform int32 index,
                            void *uniform volume)
{
  Model *uniform model = (Model *uniform) pointer;
  model->volumes[index] = (Volume *uniform) volume;
}

