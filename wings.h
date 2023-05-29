#pragma once

#include <cstdint>
#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace wings {

enum class StatusCode {
  kCreateSocketError,
  kBindError,
  kListenError,
  kConnectError,
  kHandshakeError,
  kFrameError,
  kSuccess
};

static const std::map<StatusCode, std::string> RenderStatusString = {
    {StatusCode::kCreateSocketError, "cannot create socket"},
    {StatusCode::kBindError, "cannot bind socket on port"},
    {StatusCode::kListenError, "cannot listen on port"},
    {StatusCode::kConnectError, "cannot connect to client"},
    {StatusCode::kSuccess, "success"}};

enum class InputType : uint32_t {
  MouseMotion,
  SingleClick,
  Scroll,
  DoubleClick,
  KeyValueBool,
  KeyValueInt,
  KeyValueFloat,
  KeyValueStr,
  KeyPressed,
  KeyReleased
};

enum class KeyModifier { Shift, Control, Alt, Option };

struct ClientInput {
  InputType type;
  bool dragging{false};
  bool modifier{false};
  float x{0.0f};
  float y{0.0f};
  char key{'\0'};
  bool bvalue{false};
  float fvalue{0.0f};
  int ivalue{0};
  const char* svalue{nullptr};
};

enum RenderingContextType : uint8_t { kOpenGL };

struct RenderingContext {
  RenderingContext(RenderingContextType _type) : type(_type) {}
  void print(bool with_extensions = false);
  virtual ~RenderingContext(){};
  static std::unique_ptr<RenderingContext> create(const RenderingContext& ctx);
  static std::unique_ptr<RenderingContext> create(RenderingContextType type);

  void enter_render_section() {
    lock_.lock();
    make_context_current();
  }

  void leave_render_section() {
    swap_buffers();
    release_context();
    lock_.unlock();
  }

  virtual void make_context_current() = 0;
  virtual void release_context() = 0;
  virtual void swap_buffers() = 0;

  std::mutex lock_;
  RenderingContextType type;
};

class Scene {
 public:
  virtual bool render(const ClientInput& input, int client_idx) = 0;
  virtual void onconnect() = 0;
  virtual ~Scene() {}

  const RenderingContext& context() const { return *context_; }

  int quality() const { return quality_; }
  const std::vector<unsigned char>& pixels() const { return pixels_; }

  RenderingContext& context() { return *context_; }

 protected:
  std::unique_ptr<RenderingContext> context_;
  std::vector<unsigned char> pixels_;
  int quality_{80};
};

class WebsocketRenderer;
class RenderingServer {
 public:
  RenderingServer(Scene& scene, int port);
  ~RenderingServer();
  StatusCode start(const std::string& html_file, int port);
  StatusCode stop();

 private:
  std::unique_ptr<WebsocketRenderer> renderer_;
  std::atomic<bool> listen_;
  std::future<StatusCode> server_;
  Scene& scene_;
};

}  // namespace wings