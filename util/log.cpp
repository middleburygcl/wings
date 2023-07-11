#include "log.h"

#include <cxxabi.h>
#include <dlfcn.h>
#include <execinfo.h>
#include <sys/types.h>
#include <unistd.h>

#include <csignal>
#include <cstdarg>
#include <cstdlib>

namespace wings {

std::string Exception::get_backtrace(int start_frame) {
  int i;
  enum { MAX_DEPTH = 50 };
  void *trace[MAX_DEPTH];
  char *demangled;
  int trace_size, status = 0;
  Dl_info dlinfo;
  const char *symname;

  std::string message = "\n";

  trace_size = backtrace(trace, MAX_DEPTH);

  for (i = start_frame; i < trace_size; i++) {
    if (!dladdr(trace[i], &dlinfo)) continue;
    symname = dlinfo.dli_sname;

    demangled = abi::__cxa_demangle(symname, NULL, 0, &status);
    if (status == 0 && demangled) symname = demangled;

    if (symname) {
      message += symname;
      message += '\n';
    }

    if (demangled) {
      free(demangled);
    }
  }
  return message;
}

}  // namespace wings
