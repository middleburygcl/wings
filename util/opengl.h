#pragma once

#define GL_GLEXT_PROTOTYPES
#ifdef __APPLE__
#define __gl_h_
#define GL_DO_NOT_WARN_IF_MULTI_GL_VERSION_HEADERS_INCLUDED
#define GL_SILENCE_DEPRECATION
#include <OpenGL/OpenGL.h>
#include <OpenGL/gl3.h>
#else
#include <GL/gl.h>
#endif

#define WINGS360_GL_VERSION_MAJOR 3
#define WINGS360_GL_VERSION_MINOR 3

#define GL_CALL(X)                                                          \
  {                                                                         \
    (X);                                                                    \
    [[maybe_unused]] bool error = false;                                    \
    for (GLenum glerr = glGetError(); glerr != GL_NO_ERROR;                 \
         glerr = glGetError()) {                                            \
      const char* message = "";                                             \
      switch (glerr) {                                                      \
        case GL_INVALID_ENUM:                                               \
          message = "invalid enum";                                         \
          break;                                                            \
        case GL_INVALID_VALUE:                                              \
          message = "invalid value";                                        \
          break;                                                            \
        case GL_INVALID_OPERATION:                                          \
          message = "invalid operation";                                    \
          break;                                                            \
        case GL_INVALID_FRAMEBUFFER_OPERATION:                              \
          message = "invalid framebuffer operation";                        \
          break;                                                            \
        case GL_OUT_OF_MEMORY:                                              \
          message = "out of memory";                                        \
          break;                                                            \
        default:                                                            \
          message = "unknown error";                                        \
      }                                                                     \
      std::cout << "OpenGL error at " << __FILE__ << "(" << __LINE__ << ")" \
                << message << std::endl;                                    \
      error = true;                                                         \
    }                                                                       \
    ASSERT(!error);                                                         \
  }
  