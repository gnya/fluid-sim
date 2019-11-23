/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 09
 */

#include <iostream>
#include <string>
#include <memory>
#include <regex>
#include <thread>
#include <filesystem>
#include <algorithm>
#include <GL/glut.h>
#include <boost/program_options.hpp>

#define __STDC_LIB_EXT1__
#include "lib/ink_fluid.h"
#include "lib/map_fluid.h"

#define MAP_FLUID_MODE 0
#define INK_FLUID_MODE 1

// setting
int width, height;
float scale_x, scale_y;
int simulation_mode = -1;
bool save_enable = false;

// simulator
std::shared_ptr<fluid::Fluid> sim;
auto map_sim = []{return std::dynamic_pointer_cast<fluid::MapFluid>(sim);};
auto ink_sim = []{return std::dynamic_pointer_cast<fluid::InkFluid>(sim);};

// image
namespace png {
  template<typename pixel>
  using solid_image = image<pixel, solid_pixel_buffer<pixel>>;
}
int src_frame_number = 0;
std::vector<png::pixel_buffer<png::rgb_pixel>> src_buf_seq;
std::string dst_img_dir;
std::shared_ptr<png::solid_image<png::rgb_pixel>> dst_img_ptr;

// mouse
float mouse_pos[2] = {};
float mouse_pos_old[2] = {};
float mouse_drag[2] = {};
bool mouse_pressed = false;
int mouse_button = -1;

// save frame
int frame = 0;

void mouse(int button , int state , int x , int y) {
  if (state == GLUT_DOWN) {
    mouse_pressed = true;
    mouse_button = button;
  } else if (state == GLUT_UP) {
    mouse_pressed = false;
    mouse_button = -1;
  }
}

void motion(int x, int y) {
  mouse_pos[0] = (float) x;
  mouse_pos[1] = (float) (height - y);
}

void passive_motion(int x, int y) {
  mouse_pos[0] = (float) x;
  mouse_pos[1] = (float) (height - y);
}

void update_drag() {
  if (mouse_pressed) {
    mouse_drag[0] = mouse_pos[0] - mouse_pos_old[0];
    mouse_drag[1] = mouse_pos[1] - mouse_pos_old[1];
  }
  mouse_pos_old[0] = mouse_pos[0];
  mouse_pos_old[1] = mouse_pos[1];
}

void key(unsigned char key, int x, int y) {
  switch (key) {
    case 'q':
    case 'Q':
    case '\033':  // escape key
      std::exit(0);
    default:
      break;
  }
}

void update_fluid() {
  using namespace fluid::util;

  timer([]{sim->advect_velocity();}, "advect_velocity", 1, false);
  //sim->advect_velocity();
  timer([]{sim->diffuse(20);}, "diffuse", 1, false);
  //sim->diffuse(20);
  if (mouse_pressed) {
    // add force
    float p[] = {scale_x * mouse_pos[0] , scale_y * mouse_pos[1] };
    float f[] = {500     * mouse_drag[0], 500     * mouse_drag[1]};
    timer([&]{sim->add_force(p, f, 50);}, "add_force");
    //sim->add_force(p, f, 100);

    if (simulation_mode == MAP_FLUID_MODE) {
      // nothing to do
    } else if (simulation_mode == INK_FLUID_MODE) {
      // draw point
      float color[] = {
          mouse_button == GLUT_LEFT_BUTTON ? 255.0f * 10 : 0,
          mouse_button == GLUT_RIGHT_BUTTON ? 255.0f * 10 : 0
      };
      float pos[] = {mouse_pos[0], mouse_pos[1]};
      timer([&]{ink_sim()->add_ink(pos, color, 100);}, "add_ink");
      //ink_sim()->add_ink(pos, color, 100);
    }
  }
  timer([]{sim->projection(40);}, "projection", 1, false);
  //sim->projection(40);
  timer([]{sim->boundary();}, "boundary", 1, false);
  //sim->boundary();

  if (simulation_mode == MAP_FLUID_MODE) {
    timer([]{map_sim()->advect_map();}, "advect_map", 1, false);
    //map_sim()->advect_map();
  } else if (simulation_mode == INK_FLUID_MODE) {
    timer([]{ink_sim()->advect_ink();}, "advect_ink", 1, false);
    //ink_sim()->advect_ink();
  }
}

template<typename pixel>
void get_dst_buf(png::solid_pixel_buffer<pixel>& dst_buf) {
  if (simulation_mode == MAP_FLUID_MODE) {
    auto src_buf = src_buf_seq[src_frame_number++];
    map_sim()->mapping(src_buf, dst_buf);

    if (src_frame_number >= src_buf_seq.size()) {
      src_frame_number = 0;
    }
  } else if (simulation_mode == INK_FLUID_MODE) {
    ink_sim()->get_ink(dst_buf);
  }
}

template<typename pixel>
void save_dst_img(const png::solid_image<pixel>& dst_img) {
  std::cout << "[info] saving frame:" << frame << std::endl;
  std::string filename = dst_img_dir + "/frame-" + std::to_string(frame) + ".png";
  std::thread th([](png::solid_image<pixel> dst_img, const std::string& filename){
    dst_img.write(filename);
    std::cout << "[info] saved:" << filename << std::endl;
  }, dst_img, filename);

  th.detach();
  frame++;
}

void rendering(const GLvoid *dst_data) {
  glClearColor(1.0, 1.0, 1.0, 10.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, dst_data);
  glFlush();
}

void display() {
  using namespace fluid::util;

  std::cout << "[info] calc fluid" << std::endl;

  // update phase
  update_drag();
  fluid::util::timer([]{update_fluid();}, "AMOUNT OF FLUID SIM");
  //update_fluid();

  std::cout << "[info] draw phase" << std::endl;

  // init dst img
  timer([]{get_dst_buf(dst_img_ptr->get_pixbuf());}, "AMOUNT OF GET DST IMG", 1, false);
  //get_dst_buf(dst_img.get_pixbuf());
  timer([]{if (save_enable) save_dst_img(*dst_img_ptr);}, "AMOUNT OF SAVING", 1, false);
  //if (save_enable) save_dst_img(dst_img);

  // rendering
  timer([]{rendering(dst_img_ptr->get_pixbuf().get_bytes().data());}, "AMOUNT OF RENDERING", 1, false);
  //rendering(dst_img.get_pixbuf().get_bytes().data());
  glutPostRedisplay();
}

void load_png_seq(std::string prefix) {
  namespace fs = std::filesystem;

  // unify separator
  char separator = (char) fs::path::preferred_separator;
  std::replace(prefix.begin(), prefix.end(), '/', separator);
  std::replace(prefix.begin(), prefix.end(), '\\', separator);

  fs::path p(prefix);
  fs::path dir = p.parent_path();

  // escape
  std::regex esc(R"([\\*\+\.\?\|\{\}\(\)\[\]\^\$\-])");
  auto prefix_esc = std::regex_replace(p.string(), esc, "\\$&");

  // find file
  std::regex seq_file(prefix_esc + "[0-9]*\\.png");
  fs::directory_iterator iter(dir), end;
  std::error_code err;
  for (; iter != end && !err; iter.increment(err)) {
    std::string filename = iter->path().string();

    if (std::regex_match(filename, seq_file)) {
      png::image<png::rgb_pixel> src_img(filename);
      src_buf_seq.push_back(src_img.get_pixbuf());

      std::cout << "[info] load: " << filename << std::endl;
    }
  }
}

int main(int argc, char *argv[]) {
  using namespace boost::program_options;

  options_description opt("Option");
  opt.add_options()
    ("help,H", "Print help message")
    ("save,s", "Save destination image")
    ("sequence", "Use png sequence file for source image")
    ("src_img", value<std::string>()->default_value("image"), "Source image")
    ("dst_img", value<std::string>()->default_value("image"), "Destination image")
    ("scale_x", value<float>()->default_value(0.5f), "Field x scale")
    ("scale_y", value<float>()->default_value(0.5f), "Field y scale")
    ("mode,m", value<std::string>()->default_value("Map"), "Simulation mode")
    ("init_from_src", "Init ink from source image (Only Ink Fluid)")
    ("noise_amp", value<float>()->default_value(10), "Perlin noise amplitude")
    ("vel_x", value<float>()->default_value(0), "Fluid x velocity")
    ("vel_y", value<float>()->default_value(0), "Fluid y velocity")
    ("delta_t", value<float>()->default_value(0.03f), "Simulation delta t")
    ("delta_x", value<float>()->default_value(1), "Simulation delta x")
    ("viscosity", value<float>()->default_value(0.5f), "Fluid viscosity");

  variables_map map;
  try {
    store(parse_command_line(argc, argv, opt), map);
  } catch (const error_with_option_name& e) {
    std::cout << "[error] " << e.what() << std::endl;

    return -1;
  }
  notify(map);

  if (map.count("help")) {
    // print help message
    std::cout << "[help] " << opt << std::endl;

    return 0;
  }

  try {
    // select save enable or disable
    save_enable = map.count("save");

    std::string src_img_prefix = map["src_img"].as<std::string>();
    if (map.count("sequence")) {
      // load source image sequence
      load_png_seq(src_img_prefix);
      if (src_buf_seq.empty()) {
        std::cout << "[error] source image sequence was not exist." << std::endl;

        return -1;
      }
    } else {
      // load source image
      std::string filename = src_img_prefix + ".png";
      png::image<png::rgb_pixel> src_img(filename);
      src_buf_seq.push_back(src_img.get_pixbuf());

      std::cout << "[info] load: " << filename << std::endl;
    }

    // set width and height
    width = src_buf_seq[0].get_width();
    height = src_buf_seq[0].get_height();

    // set scale_x and scale_y
    scale_x = map["scale_x"].as<float>();
    scale_y = map["scale_y"].as<float>();

    // fluid property
    float dt = map["delta_t"].as<float>();
    float dx = map["delta_x"].as<float>();
    float viscosity = map["viscosity"].as<float>();

    // select simulation mode
    auto mode_name = map["mode"].as<std::string>();
    if (mode_name == "Map" || mode_name == "map") {
      simulation_mode = MAP_FLUID_MODE;

      // init simulator
      sim = std::make_shared<fluid::MapFluid>(width, height, scale_x, scale_y, dt, dx, viscosity);
    } else if (mode_name == "Ink" || mode_name == "ink") {
      simulation_mode = INK_FLUID_MODE;

      // init simulator
      sim = std::make_shared<fluid::InkFluid>(width, height, scale_x, scale_y, dt, dx, viscosity);

      // init ink from image
      if (map.count("init_from_src")) {
        auto src_buf = src_buf_seq[0];
        ink_sim()->set_ink(src_buf);
      }
    } else {
      std::cout << "[error] mode \"" << mode_name << "\" was not exist." << std::endl;

      return -1;
    }

    // set perlin noise amplitude
    fluid::Fluid::accelerate_by_perlin_noise(*sim, 0, 1, map["noise_amp"].as<float>());

    // set fluid velocity
    float vel[] = {map["vel_x"].as<float>(), map["vel_y"].as<float>()};
    fluid::Fluid::accelerate_by_single_vector(*sim, vel);

    // init destination image
    dst_img_ptr = std::make_shared<png::solid_image<png::rgb_pixel>>(width, height);
    if (save_enable) {
      namespace fs = std::filesystem;

      // set destination image filename
      dst_img_dir = map["dst_img"].as<std::string>();

      fs::path dir(dst_img_dir);
      if (!fs::exists(dir)) {
        // create directory
        fs::create_directories(dir);

        std::cout << "[info] create directories: " << dir << std::endl;
      }
    }
  } catch (const boost::bad_any_cast& e) {
    std::cout << "[error] " << e.what() << std::endl;

    return -1;
  }

  // initialize OpenGL
  glutInit(&argc, argv);

  glutInitWindowSize(width, height);
  glutCreateWindow("fluid sim");
  glutInitDisplayMode(GLUT_RGBA);

  glutPassiveMotionFunc(passive_motion);
  glutMotionFunc(motion);
  glutMouseFunc(mouse);
  glutKeyboardFunc(key);

  glutDisplayFunc(display);
  glutMainLoop();

  return 0;
}