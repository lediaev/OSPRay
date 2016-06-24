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

#include "Material.ih"
#include "../bsdfs/Dielectric.ih"
#include "ospray/texture/TextureParam.ih"

struct Glass
{
  PathTraceMaterial super;

  Medium mediumInside;
  Medium mediumOutside;

  float absorptionDistance;
  vec3f absorptionColor;
  vec3f absorption_k;

  TextureParam colorMap;
};

///////////////////////////////////////////////////////////////////////////////
// Implementation

const varying BSDF* uniform Glass_getBSDF(const uniform PathTraceMaterial* uniform super,
                                          uniform ShadingContext* uniform ctx,
                                          const DifferentialGeometry& dg,
                                          const Ray& ray,
                                          const Medium& currentMedium)
{
  uniform const Glass* uniform self = (uniform const Glass* uniform)super;

  float eta = eq(currentMedium, self->mediumOutside) ? self->mediumOutside.ior*rcp(self->mediumInside.ior)
                                                     : self->mediumInside.ior*rcp(self->mediumOutside.ior);

  varying BSDF* uniform bsdf = Dielectric_create(ctx, dg, eta);
  return bsdf;
}

vec3f Glass_getTransparency(const uniform PathTraceMaterial* uniform material,
                            const DifferentialGeometry& dg,
                            const Ray& ray,
                            const Medium& currentMedium)
{
  const uniform Glass* uniform self = (const uniform Glass* uniform)material;

  float eta = eq(currentMedium, self->mediumOutside) ? self->mediumOutside.ior*rcp(self->mediumInside.ior)
                                                     : self->mediumInside.ior*rcp(self->mediumOutside.ior);

  float cosThetaO = max(-dot(ray.dir, dg.Ns), 0.0f);
  return make_vec3f(1.0f-fresnelDielectric(cosThetaO, eta));
}

void Glass_selectNextMedium(const uniform PathTraceMaterial* uniform super,
                            Medium& currentMedium)
{
  const uniform Glass* uniform self = (const uniform Glass* uniform)super;

  if (eq(currentMedium, self->mediumOutside))
    currentMedium = self->mediumInside;
  else
    currentMedium = self->mediumOutside;
}

vec3f Glass_getAbsorption(const uniform PathTraceMaterial* uniform material,
                          const DifferentialGeometry& current_dg, 
                          const DifferentialGeometry& previous_dg, 
                          const float & distance)
{
    const uniform Glass* uniform self = (const uniform Glass* uniform)material;

    vec3f r = make_vec3f(1.0f);

    // Beer's law for glass absorption.

    vec3f k = make_vec3f(0.0f);

    if (self->colorMap.map != NULL) {
        // Using a 2D texture for absorption color will cause the glass color
        // to be view-dependent. As you move the camera the colors 
        // will change dynamically. This is fun, but entirely not physically correct.
        // The correct thing is to have a 3D texture with ray-marching through
        // the solid to integrate the absorption, similar (or identical?) to volume rendering.

        vec3f color1 = get3f(self->colorMap, previous_dg.st);
        vec3f color2 = get3f(self->colorMap, current_dg.st);

        vec3f color;// = 0.5*(color1 + color2);
        //color = color1;
        //color = color2;

        color = make_vec3f(0.5*(color1.x + color2.x),0.5*(color1.y + color2.y),0.5*(color1.z + color2.z));

        //print("color = % % % , color1 = % % % , color2 = % % %\n\n\n", color.x, color.y, color.z, color1.x, color1.y, color1.z, color2.x, color2.y, color2.z);

        if (self->absorptionDistance > 0) {
            if (color.x > 0)
                k.x = log(color.x)/self->absorptionDistance;

            if (color.y > 0)
                k.y = log(color.y)/self->absorptionDistance;

            if (color.z > 0)
                k.z = log(color.z)/self->absorptionDistance;
        }
    } else {
        k = self->absorption_k;
    }

    if (k.x != 0)
        r.x = expf(k.x * distance);
    if (k.y != 0)
        r.y = expf(k.y * distance);
    if (k.z != 0)
        r.z = expf(k.z * distance);

    return r;
}

///////////////////////////////////////////////////////////////////////////////
// External API

export void PathTracer_Glass_set(
  void* uniform _self,
  const uniform float iorInside,
  const uniform vec3f & transmissionInside,
  const uniform float iorOutside,
  const uniform vec3f & transmissionOutside,
  const uniform float absDistance,
  const uniform vec3f & absColor,
  void* uniform map_Kd, 
  const uniform affine2f & xform_Kd)
{
  uniform Glass* uniform self = (uniform Glass* uniform)_self;

  self->mediumInside.ior           = iorInside;
  self->mediumInside.transmission  = transmissionInside;
  self->mediumOutside.ior          = iorOutside;
  self->mediumOutside.transmission = transmissionOutside;
  self->absorptionDistance         = absDistance;
  self->absorptionColor            = absColor;
  self->colorMap = make_TextureParam((uniform Texture2D*)map_Kd, xform_Kd);

  if (absDistance > 0) {
    if (absColor.x > 0)
      self->absorption_k.x = log(absColor.x)/absDistance;
    else 
      self->absorption_k.x = 0;

    if (absColor.y > 0)
      self->absorption_k.y = log(absColor.y)/absDistance;
    else 
      self->absorption_k.y = 0;

    if (absColor.z > 0)
      self->absorption_k.z = log(absColor.z)/absDistance;
    else 
      self->absorption_k.z = 0;
  }
}

void Glass_Constructor(uniform Glass* uniform self)
{
  PathTraceMaterial_Constructor(&self->super, Glass_getBSDF, 
    Glass_getTransparency, Glass_selectNextMedium, Glass_getAbsorption);
}

export void* uniform PathTracer_Glass_create()
{
  uniform Glass* uniform self = uniform new uniform Glass;
  Glass_Constructor(self);
  return self;
}
