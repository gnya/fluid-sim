/**
 * author: gnya / https://twitter.com/gnya_a
 * date:   2019 11 09
 */

#include <iostream>
#include <string>
#include <GL/glut.h>
#include <memory>

#include "lib/ink_fluid.h"
#include "lib/map_fluid.h"

class Timer {
public:
  clock_t _start{};

  Timer() = default;

  void start() {
    _start = clock();
  }

  void stop(const char* str) {
    clock_t stop = clock();
    double time = (double) (stop - _start) / CLOCKS_PER_SEC * 1000;

    std::cout << "* " << str << " time:" << time << "ms" << std::endl;
  }
};

#define MAP_FLUID_MODE 0
#define INK_FLUID_MODE 1

// setting
int width, height;
float scale_x, scale_y;
float noise_amplitude = 0;
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

void display() {
  update_drag();

  std::cout << "calc fluid" << std::endl;
  sim_timer.start();

  sim->advect_velocity();
  sim->diffuse(20);
  if (mouse_pressed) {
    // add force
    float p[] = {scale_x * mouse_pos[0] , scale_y * mouse_pos[1] };
    float f[] = {1500    * mouse_drag[0], 1500    * mouse_drag[1]};
    sim->add_force(p, f, 100);

    if (simulation_mode == MAP_FLUID_MODE) {
      // nothing to do
    } else if (simulation_mode == INK_FLUID_MODE) {
      // draw point
      /*float color[] = {255, 255, 255};
      float c_pos[] = {mouse_pos[0], mouse_pos[1]};
      sim.add_color(c_pos, color, 100);*/
    }
  }
  sim->projection(40);
  sim->boundary();

  if (simulation_mode == MAP_FLUID_MODE) {
    std::dynamic_pointer_cast<fluid::MapFluid>(sim)->advect_map();
    std::dynamic_pointer_cast<fluid::MapFluid>(sim)->mapping(src_img.get_pixbuf(), dst_img->get_pixbuf());
  } else if (simulation_mode == INK_FLUID_MODE) {
    std::dynamic_pointer_cast<fluid::InkFluid>(sim)->advect_ink();
    std::dynamic_pointer_cast<fluid::InkFluid>(sim)->get_ink(dst_img->get_pixbuf());
  }

  sim_timer.stop("SIM AMOUNT OF TIME:");

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

int main(int argc, char **argv) {
  // select save enable or disable
  std::string save_enable_input;
  std::cout << "If you want enable frame save, press 'Y' key. (Default: disable)" << std::endl << ">> ";
  getline(std::cin, save_enable_input);
  save_enable = save_enable_input[0] == 'Y' || save_enable_input[0] == 'y';

  // set width and height
  std::string width_input, height_input;
  std::cout << "Please enter width. (Default: 1000)" << std::endl << ">> ";
  getline(std::cin, width_input);
  width = width_input.empty() ? 1000 : std::stoi(width_input);
  std::cout << "Please enter height. (Default: 1000)" << std::endl << ">> ";
  getline(std::cin, height_input);
  height = height_input.empty() ? 1000 : std::stoi(height_input);
  std::cout << "Size: width x height = " << width << " x " << height << std::endl;

  // set scale_x and scale_y
  std::string scale_x_input, scale_y_input;
  std::cout << "Please enter scale_x. (Default: 0.5f)" << std::endl << ">> ";
  getline(std::cin, scale_x_input);
  scale_x = scale_x_input.empty() ? 0.5f : std::stof(scale_x_input);
  std::cout << "Please enter scale_y. (Default: 0.5f)" << std::endl << ">> ";
  getline(std::cin, scale_y_input);
  scale_y = scale_y_input.empty() ? 0.5f : std::stof(scale_y_input);
  std::cout << "Field scale: scale_x x scale_y = " << scale_x << " x " << scale_y << std::endl;

  // load source image
  std::string filename_input;
  std::cout << "Please enter source image name. (Default: text)" << std::endl << ">> ";
  getline(std::cin, filename_input);
  if (filename_input.empty()) filename_input = "text";
  src_img.read(filename_input + ".png");

  // init destination image
  dst_img = std::make_shared<png::solid_image<png::rgb_pixel>>(width, height);

  // select simulation mode
  std::string sim_mode_input;
  std::cout << "Please select simulation mode. Map Fluid Mode: 0 / Ink Fluid Mode: 1" << std::endl << ">> ";
  getline(std::cin, sim_mode_input);
  if (sim_mode_input.empty() || sim_mode_input[0] == '0') {
    std::cout << "Map Fluid Mode was selected" << std::endl;
    simulation_mode = MAP_FLUID_MODE;

    // init simulator
    sim = std::make_shared<fluid::MapFluid>(width, height, scale_x, scale_y, 0.03f, 1.0f, 0.0001f);
  } else if (sim_mode_input[0] == '1') {
    std::cout << "Ink Fluid Mode was selected" << std::endl;
    simulation_mode = INK_FLUID_MODE;

    // init simulator
    sim = std::make_shared<fluid::InkFluid>(width, height, scale_x, scale_y, 0.03f, 1.0f, 0.0001f);

    // init ink from image
    std::string ink_from_image_input;
    std::cout << "If you want init ink from source image, press 'Y' key. (Default: disable)" << std::endl << ">> ";
    getline(std::cin, ink_from_image_input);
    if (ink_from_image_input[0] == 'Y' || ink_from_image_input[0] == 'y') {
      std::dynamic_pointer_cast<fluid::InkFluid>(sim)->set_ink(src_img.get_pixbuf());
    }
  } else {
    return -1;
  }

  // set perlin noise amplitude
  std::string noise_amplitude_input;
  std::cout << "Please enter noise_amplitude. (Default: 10)" << std::endl << ">> ";
  getline(std::cin, noise_amplitude_input);
  noise_amplitude = noise_amplitude_input.empty() ? 10 : std::stof(noise_amplitude_input);
  std::cout << "Noise Amplitude is " << noise_amplitude << std::endl;
  fluid::Fluid::accelerate_by_perlin_noise(*sim, 0, 1, noise_amplitude);

  // set fluid velocity
  std::string vel_x_input, vel_y_input;
  float vel[2];
  std::cout << "Please enter vel_x. (Default: 0)" << std::endl << ">> ";
  getline(std::cin, vel_x_input);
  vel[0] = vel_x_input.empty() ? 0 : std::stoi(vel_x_input);
  std::cout << "Please enter vel_y. (Default: 0)" << std::endl << ">> ";
  getline(std::cin, vel_y_input);
  vel[1] = vel_y_input.empty() ? 0 : std::stoi(vel_y_input);
  std::cout << "Velocity: vel_x x vel_y = " << vel[0] << " x " << vel[1] << std::endl;
  fluid::Fluid::accelerate_by_single_vector(*sim, vel);

  // initialize OpenGL
  glutInit(&argc, argv);

  glutInitWindowSize(width, height);
  glutCreateWindow("fluid sim");
  glutInitDisplayMode(GLUT_RGBA);

  glutPassiveMotionFunc(passive_motion);
  glutMotionFunc(motion);
  glutMouseFunc(mouse);

  glutDisplayFunc(display);
  glutMainLoop();

  return 0;
}