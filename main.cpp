/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 09
 */

#include <iostream>
#include <string>
#include <memory>
#include <GL/glut.h>
#include <boost/program_options.hpp>

#define __STDC_LIB_EXT1__
#include "lib/ink_fluid.h"
#include "lib/map_fluid.h"

class Timer {
public:
  clock_t _start{};

  Timer() = default;

  void start() {
    _start = clock();
  }

  void stop(const char* prefix) {
    clock_t stop = clock();
    double time = (double) (stop - _start) / CLOCKS_PER_SEC * 1000;

    std::cout << "[" << prefix << "] time: " << time << "ms" << std::endl;
  }

  template<class Process>
  void operator()(Process p, const char* prefix, int loop = 1) {
    start();
    for (int i = 0; i < loop; i++) p();
    stop(prefix);
  }
};

#define MAP_FLUID_MODE 0
#define INK_FLUID_MODE 1

// setting
int width, height;
float scale_x, scale_y;
int simulation_mode = -1;
bool save_enable = false;

// simulator
std::shared_ptr<fluid::Fluid> sim;
Timer timer, sim_timer;

// image
namespace png {
  template< typename pixel>
  using solid_image = image<pixel, solid_pixel_buffer<pixel>>;
}
png::image<png::rgb_pixel> src_img;
std::shared_ptr<png::solid_image<png::rgb_pixel>> dst_img;

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
  std::cout << "calc fluid" << std::endl;
  sim_timer.start();

  timer([]{sim->advect_velocity();}, "advect_velocity");
  //sim->advect_velocity();
  timer([]{sim->diffuse(20);}, "diffuse");
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
      /*float color[] = {255, 255, 255};
      float c_pos[] = {mouse_pos[0], mouse_pos[1]};
      sim.add_color(c_pos, color, 100);*/
    }
  }
  timer([]{sim->projection(40);}, "projection");
  //sim->projection(40);
  timer([]{sim->boundary();}, "boundary");
  //sim->boundary();

  if (simulation_mode == MAP_FLUID_MODE) {
    timer([]{std::dynamic_pointer_cast<fluid::MapFluid>(sim)->advect_map();}, "advect_map");
    //timer([]{std::dynamic_pointer_cast<fluid::MapFluid>(sim)->advect_map_new();}, "advect_map_new", 1);
    std::dynamic_pointer_cast<fluid::MapFluid>(sim)->advect_map();
  } else if (simulation_mode == INK_FLUID_MODE) {
    timer([]{std::dynamic_pointer_cast<fluid::InkFluid>(sim)->advect_ink();}, "advect_ink");
    //std::dynamic_pointer_cast<fluid::InkFluid>(sim)->advect_ink();
  }
}

void display() {
  // update phase
  update_drag();
  sim_timer([]{update_fluid();}, "AMOUNT OF FLUID SIM");
  //update_fluid();

  // update dst_img
  if (simulation_mode == MAP_FLUID_MODE) {
    std::dynamic_pointer_cast<fluid::MapFluid>(sim)->mapping(src_img.get_pixbuf(), dst_img->get_pixbuf());
  } else if (simulation_mode == INK_FLUID_MODE) {
    std::dynamic_pointer_cast<fluid::InkFluid>(sim)->get_ink(dst_img->get_pixbuf());
  }

  // draw image
  std::cout << "draw phase" << std::endl;
  glClearColor(1.0, 1.0, 1.0, 10.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, dst_img->get_pixbuf().get_bytes().data());
  glFlush();

  // save frame
  if (save_enable) {
    dst_img->write("img_output/frame-" + std::to_string(frame++) + ".png");
  }

  glutPostRedisplay();
}

int main(int argc, char *argv[]) {
  using namespace boost::program_options;

  options_description opt("Option");
  opt.add_options()
    ("help,H", "Print help message")
    ("save,s", "Save destination image")
    ("width,w", value<int>()->default_value(1440), "Window width")
    ("height,h", value<int>()->default_value(1080), "Window height")
    ("scale_x", value<float>()->default_value(0.5f), "Field x scale")
    ("scale_y", value<float>()->default_value(0.5f), "Field y scale")
    ("src_img", value<std::string>()->default_value("20191028.png"), "Source image")
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

    // set width and height
    width = map["width"].as<int>();
    height = map["height"].as<int>();

    // set scale_x and scale_y
    scale_x = map["scale_x"].as<float>();
    scale_y = map["scale_y"].as<float>();

    // load source image
    src_img.read(map["src_img"].as<std::string>());

    // init destination image
    dst_img = std::make_shared<png::solid_image<png::rgb_pixel>>(width, height);

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
        std::dynamic_pointer_cast<fluid::InkFluid>(sim)->set_ink(src_img.get_pixbuf());
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