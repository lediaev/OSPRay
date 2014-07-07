#include "pointlight.h"

namespace ospray {
  //!Construct a new PointLight object
  PointLight::PointLight()
    : position(0.f, 0.f, 0.f)
    , color(1.f, 1.f, 1.f)
    , range(-1.f)
  {
  }

  //!Commit parameters understood by the base PointLight class
  void PointLight::commit() {
    position = getParam3f("position", vec3f(0.f, 0.f, 0.f));
    color    = getParam3f("color", vec3f(1.f, 1.f, 1.f));
    range    = getParam1f("range", -1.f);
  }
}