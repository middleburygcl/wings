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

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "glm.h"
#include "log.h"
#include "opengl.h"

#ifndef WINGS_SOURCE_DIR
#define WINGS_SOURCE_DIR "undefined"  // to silence vs-code warning
#endif

namespace wings {

class ShaderProgram {
 public:
  ShaderProgram() {
    handle_ = glCreateProgram();
    ASSERT(handle_ >= 0);
  }
  void set_source(const std::string& base, const std::string& name,
                  bool with_geometry, bool with_tessellation,
                  const std::vector<std::string>& macros);
  void compile(const std::string& src, GLenum type);
  void compile(const std::string& vs, const std::string& fs,
               const std::string& gs, const std::string& tcs,
               const std::string& tes);

  template <typename T>
  void set_uniform(const char* name, const T& value) const;
  void use() const;
  int handle() const { return handle_; }

 private:
  int handle_{-1};
  std::vector<std::string> macros_;
};

namespace {
std::string get_shader(const std::string& filename) {
  std::string content;
  std::ifstream file_stream(filename.c_str(), std::ios::in);

  ASSERT(file_stream.is_open()) << "could not read file " << filename;

  std::string line = "";
  while (!file_stream.eof()) {
    std::getline(file_stream, line);
    content.append(line + "\n");
  }

  file_stream.close();
  return content;
}
}  // namespace

void ShaderProgram::set_source(const std::string& dir, const std::string& name,
                               bool with_geometry, bool with_tessellation,
                               const std::vector<std::string>& macros) {
  macros_ = macros;

  std::string base = dir + name;

  // compile the program using each stage, and then link
  std::string vs = get_shader(base + "-vs.glsl");
  std::string fs = get_shader(base + "-fs.glsl");
  std::string gs, tcs, tes;
  if (with_geometry) gs = get_shader(base + "-gs.glsl");
  if (with_tessellation) {
    tcs = get_shader(base + "-tcs.glsl");
    tes = get_shader(base + "-tes.glsl");
  }
  compile(vs, fs, gs, tcs, tes);
}

void ShaderProgram::compile(const std::string& src, GLenum type) {
  ASSERT(handle_ >= 0);

  GLuint shader = glCreateShader(type);

  std::string source;
  if (macros_.size() > 0) {
    std::size_t idx = src.find("#version");
    idx = src.find("\n");

    std::string source1 = "//" + src.substr(0, idx);  // comment out version
    std::string source2 = "\n";
    for (auto& s : macros_) source2 += s + "\n";
    std::string source3 = src.substr(idx + 1, src.size());
    source = source1 + source2 + source3;
  } else
    source = src;

  const auto* c = source.data();
  GL_CALL(glShaderSource(shader, 1, &c, NULL));
  GL_CALL(glCompileShader(shader));

  // check for errors
  int result;
  GL_CALL(glGetShaderiv(shader, GL_COMPILE_STATUS, &result));
  if (result == GL_FALSE) {
    int length = 0;
    GL_CALL(glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length));
    if (length > 0) {
      std::string log(length, ' ');
      int written = 0;
      GL_CALL(glGetShaderInfoLog(shader, length, &written, &log[0]));
      std::cout << log << std::endl;
    }
    // std::cout << "shader source:\n" << source << std::endl;
  }
  ASSERT(result == GL_TRUE) << "failed to compile shader stage " << type;

  // compile succeeded, attach shader to program
  GL_CALL(glAttachShader(handle_, shader));
}

void ShaderProgram::compile(const std::string& vs, const std::string& fs,
                            const std::string& gs, const std::string& tcs,
                            const std::string& tes) {
  // compile each shader stage
  compile(vs, GL_VERTEX_SHADER);
  compile(fs, GL_FRAGMENT_SHADER);
  if (!gs.empty()) compile(gs, GL_GEOMETRY_SHADER);
  if (!tcs.empty()) {
    ASSERT(!tes.empty());
    compile(tcs, GL_TESS_CONTROL_SHADER);
    compile(tes, GL_TESS_EVALUATION_SHADER);
  }

  // link the program
  GL_CALL(glLinkProgram(handle_));

  int status = 0;
  GL_CALL(glGetProgramiv(handle_, GL_LINK_STATUS, &status));
  if (status == GL_FALSE) {
    int length = 0;
    GL_CALL(glGetProgramiv(handle_, GL_INFO_LOG_LENGTH, &length));
    if (length > 0) {
      std::string log(length, '\0');
      int written = 0;
      GL_CALL(glGetProgramInfoLog(handle_, length, &written, &log[0]));
      LOG << "GLSL compiler log:\n" << log;
    }
  }
  ASSERT(status == GL_TRUE) << "failed to link program";
}

void ShaderProgram::use() const {
  ASSERT(handle_ >= 0);
  GL_CALL(glUseProgram(handle_));
}

template <>
inline void ShaderProgram::set_uniform(const char* name, const vec3f& v) const {
  GLint location = glGetUniformLocation(handle_, name);
  if (location >= 0) GL_CALL(glUniform3f(location, v[0], v[1], v[2]));
}

template <>
inline void ShaderProgram::set_uniform(const char* name, const vec2f& v) const {
  GLint location = glGetUniformLocation(handle_, name);
  if (location >= 0) GL_CALL(glUniform2f(location, v[0], v[1]));
}

template <>
inline void ShaderProgram::set_uniform(const char* name, const int& x) const {
  GLint location = glGetUniformLocation(handle_, name);
  if (location >= 0) GL_CALL(glUniform1i(location, x));
}

template <>
inline void ShaderProgram::set_uniform(const char* name, const float& x) const {
  GLint location = glGetUniformLocation(handle_, name);
  if (location >= 0) GL_CALL(glUniform1f(location, x));
}

template <>
inline void ShaderProgram::set_uniform(const char* name, const mat4f& m) const {
  GLint location = glGetUniformLocation(handle_, name);
  if (location >= 0)
    GL_CALL(glUniformMatrix4fv(location, 1, GL_FALSE, &m(0, 0)));
}

class ShaderLibrary {
 public:
  ShaderLibrary(const std::string& base) : base_(base) {}
  void create();

  void add(const std::string& name, const std::string& prefix,
           bool with_geometry, bool with_tessellation,
           const std::vector<std::string>& macros = {}) {
    shaders_.insert({name, ShaderProgram()});
    shaders_[name].set_source(base_, prefix, with_geometry, with_tessellation,
                              macros);
  }

  const ShaderProgram& operator[](const std::string& name) const {
    ASSERT(shaders_.find(name) != shaders_.end())
        << "could not find shader " << name;
    return shaders_.at(name);
  }

 private:
  std::string base_;
  std::map<std::string, ShaderProgram> shaders_;
};

void ShaderLibrary::create() {
  std::string version = "#version " +
                        std::to_string(WINGS360_GL_VERSION_MAJOR) +
                        std::to_string(WINGS360_GL_VERSION_MINOR) + "0";
  add("points", "points", false, false, {version, "#define WITH_GS 0"});
  add("nodes", "points", true, false, {version, "#define WITH_GS 1"});
  add("edges-q1-p0", "edges", true, false, {version, "#define ORDER 0"});
  add("triangles-q1-p0", "triangles", true, false,
      {version, "#define ORDER 0"});
  add("triangles-q1-p1", "triangles", true, false,
      {version, "#define ORDER 1"});
  add("triangles-q1-pt", "triangles", true, false,
      {version, "#define ORDER -1"});
  add("quads-q1-p0", "quads", true, false, {version, "#define ORDER 0"});
  add("quads-q1-p1", "quads", true, false, {version, "#define ORDER 1"});
  add("quads-q1-pt", "quads", true, false, {version, "#define ORDER -1"});
  add("polygons-q1-p0", "polygons", true, false, {version, "#define ORDER 0"});
  add("prisms-q1-p0", "prisms", true, false, {version, "#define ORDER 0"});
  add("prisms-q1-p1", "prisms", true, false, {version, "#define ORDER 1"});
  add("pyramids-q1-p0", "pyramids", true, false, {version, "#define ORDER 0"});
  add("pyramids-q1-p1", "pyramids", true, false, {version, "#define ORDER 1"});
  add("tetrahedra-q1-p0", "tet", true, false, {version, "#define ORDER 0"});
  add("tetrahedra-q1-p1", "tet", true, false, {version, "#define ORDER 1"});
  add("polyhedra-q1-p0", "polygons", true, false,
      {version, "#define POLYHEDRA"});
}

}  // namespace wings