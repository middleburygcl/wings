//
//  wings: web interface for graphics applications
//
//  Copyright 2023 Philip Claude Caplan
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
#pragma once

#include <memory>
#include <vector>

#include "glm.h"

namespace wings {

class Vertices;
template <typename T> class Topology;
class ShaderProgram;

struct AABB {
  AABB() {}
  vec3f min{1e20f, 1e20f, 1e20f};
  vec3f max{-1e20f, -1e20f, -1e20f};
};

struct GLClipPlane;

struct PickableObject {
  template <typename T>
  PickableObject(const Vertices& vertices, const Topology<T>& topology,
                 uint64_t k, const std::string& name);

  template <typename T>
  void save_points(const Vertices& vertices, const Topology<T>& topology,
                   uint64_t k);

  double intersection(const vec3f& point, const vec3f& ray,
                      const mat4f& model_matrix) const;
  double intersection(int k, const vec3f& point, const vec3f& ray,
                      const mat4f& model_matrix) const;

  int n_triangles() const { return triangles.size() / 3; }

  bool visible(const GLClipPlane& plane) const;

  std::string name;
  std::vector<vec4f> points;
  std::vector<uint64_t> triangles;
  std::vector<uint64_t> nodes;
  uint64_t index;
};

struct GLClipPlane {
  GLClipPlane();
  ~GLClipPlane();

  void initialize();

  void define(const AABB& aabb);

  void get(vec3f& point, vec3f& normal) const;

  void update();

  void draw(const mat4f& model_matrix, const mat4f& view_matrix,
            const mat4f& perspective_matrix);

  vec3f length;
  bool visible;
  float distance;
  bool active;

  vec3f center;
  vec3f coordinates[4];
  float direction;
  int dimension;
  mat4f transformation;

  int vertex_array;
  int buffer;
  std::shared_ptr<ShaderProgram> shader;
};

}  // namespace wings