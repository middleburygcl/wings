#include <algorithm>
#include <string>

#include "mesh.h"

namespace wings {

class Mesh;

namespace meshb {
void read(const std::string& filename, Mesh& mesh);
void write(const Mesh& mesh, const std::string& filename, bool twod = false);
}  // namespace meshb

namespace obj {
void read(const std::string& filename, Mesh& mesh);
void write(const Mesh& mesh, const std::string& filename);
}  // namespace obj

inline std::string get_file_ext(const std::string& filename) {
  std::string::size_type idx;
  idx = filename.rfind('.');  // find the '.' in reverse order
  if (idx != std::string::npos) return filename.substr(idx + 1);
  return "";
}

inline void read_mesh(const std::string& filename, Mesh& mesh) {
  std::string ext = get_file_ext(filename);
  std::transform(ext.begin(), ext.end(), ext.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  if (ext == "mesh" || ext == "meshb")
    meshb::read(filename, mesh);
  else if (ext == "obj")
    obj::read(filename, mesh);
  else
    NOT_IMPLEMENTED;
}

}  // namespace wings